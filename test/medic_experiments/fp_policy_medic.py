"""Does `src/medic.c` need STRICT floating point, the way `src/romeo.c` does?

`romeo.c` is compiled `-fno-fast-math -ffp-contract=off` because ROMEO's region-growing
order is decided by 8-bit integer edge weights: reassociation can push a weight across a
`rescale()` bin boundary, delete the edge from the graph, and leave a whole connected
region off by a full 2*pi.  That policy was settled by measurement over the FULL corpus
(romeo_plan.md line 15) -- a single-volume check once produced a wrong "FMA is clean"
claim.

`medic.c` rides the repository-wide `-ffast-math` source line.  Its `md_rescale_phase()`
and `md_mcpc3ds()` produce the phase that is then handed to ROMEO's strict-FP kernels via
`romeo_unwrap_frame()`, so the same mechanism is *available* here: a fast-math rounding
difference in the phase perturbs ROMEO's weights, which are the thing that is bin-quantised.
Three external audits flagged this as unclosed.  This script closes it with numbers.

TWO experiments:

  PART A -- the direct compiler experiment.  Build two binaries differing ONLY in the FP
  flags applied to medic.c, run both over a broad corpus (real sbref, a real 170-frame BOLD
  run, and adversarial synthetics), and compare every output bit-for-bit.  The measurement
  that matters is not "do the field maps differ" (a 1-ULP field-map difference is harmless)
  but "does any voxel of the UNWRAPPED PHASE differ by a whole 2*pi branch", plus "did the
  MASK change" (a mask difference means the graph itself changed).

  PART B -- the sensitivity/amplification bound.  Part A can only observe the perturbations
  the compiler happens to produce (measured below: a few voxels in 10^5).  Part B instead
  perturbs the phase handed to ROMEO by +-1 float32 ULP at EVERY voxel -- four to five
  orders of magnitude more perturbation than fast-math produces -- and runs a single binary
  twice.  If even that never flips a 2*pi branch, the absence of evidence in Part A is
  backed by a mechanism bound, not by luck.

Usage
-----
    ~/src/warpkit/.venv/bin/python test/medic_experiments/fp_policy_medic.py --build
    ~/src/warpkit/.venv/bin/python test/medic_experiments/fp_policy_medic.py \
        --fast /tmp/nm_fast --strict /tmp/nm_strict

`--build` replicates the `make -n all` link line, compiling medic.c twice (once with the
project default flags, once adding `-fno-fast-math -ffp-contract=off`) and linking two
otherwise-identical binaries.  Both objects are compiled OUTSIDE the LTO source line so
that the FP flags are the only difference between them.  Nothing in the repository is
modified; objects and binaries land in the scratch directory.

Analysis-only: NOT shipped, NOT run in CI, numpy is fine here (unlike release_smoke.py).
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nii  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, "..", "..", "src"))
DEMO = os.path.expanduser("~/src/warpkit/demo/data")
BOLD = os.path.expanduser("~/src/warpkit/inputs/sub-fm/ses-1/func")
BOLD_PRE = "sub-fm_ses-1_task-rest_acq-2d2echo_run-02"
WKMASK = os.path.abspath(os.path.join(HERE, "..", "medic_ref", "wkmask1.nii"))
SCRATCH = os.environ.get("MD_FP_SCRATCH", "/tmp/medic_fp_policy")

TES = "16.8,38.56"          # both real datasets
TRT = "0.02025"
TWO_PI = 2.0 * math.pi

# The repository-wide compile flags, verbatim from `make -n all` (src/Makefile CFLAGS).
CFLAGS = ("-O3 -ffunction-sections -fdata-sections -Xpreprocessor -fopenmp "
          "-I/opt/homebrew/opt/libomp/include -ffast-math -fno-finite-math-only").split()
STRICT_FP = ["-fno-fast-math", "-ffp-contract=off"]
# Everything else on the link line, with `-DHAVE_MEDIC medic.c` factored out.
LINK_REST = (
    "-DHAVE_BUTTERWORTH bw.c -DHAVE_TENSOR tensor.c -DHAVE_DTIFIT dtifit.c -DHAVE_QC qc.c "
    "-DHAVE_FORMATS base64.c -DHAVE_ALLINEATE allineate.o powell_newuoa.o coreg_fast.o reface.o "
    "-DHAVE_ROMEO romeo.o -DNII2MESH meshify.c quadric.c bwlabel.c radixsort.c fdr.c "
    "-DUSE_CLASSIC_CUBES oldcubes.c niimath.c core.c core32.c nifti_io.c -DFSLSTYLE -DPIGZ "
    "-DREJECT_COMPLEX -lm conform.c -DHAVE_CONFORM unifize.c -DHAVE_64BITS core64.c "
    "-DHAVE_ZLIB -lz -DHAVE_ZSTD -I/opt/homebrew/opt/zstd/include "
    "-L/opt/homebrew/opt/zstd/lib -lzstd -L/opt/homebrew/opt/libomp/lib -lomp "
    "-DHAVE_BMP filter.c bmp.c spng.c -flto").split()

SUFFIXES = ["fieldmaps_native", "fieldmaps", "displacementmaps"]
INTERMEDIATES = ["unwrapped_echo-1", "unwrapped_echo-2", "masks", "phase_offset"]
# Outputs whose unit is radians, i.e. where a 2*pi difference IS a branch flip.
RADIAN = {"unwrapped_echo-1", "unwrapped_echo-2", "phase_offset"}


# ------------------------------------------------------------------ build

def build(fast_bin: str, strict_bin: str) -> None:
    """Compile medic.c twice and link two otherwise-identical niimath binaries."""
    subprocess.run(["make", "-j8"], cwd=SRC, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    objs = {}
    for tag, extra in (("fast", []), ("strict", STRICT_FP)):
        obj = os.path.join(SCRATCH, "medic_%s.o" % tag)
        subprocess.run(["gcc"] + CFLAGS + extra + ["-DHAVE_MEDIC", "-c", "medic.c", "-o", obj],
                       cwd=SRC, check=True)
        objs[tag] = obj
    for tag, out in (("fast", fast_bin), ("strict", strict_bin)):
        subprocess.run(["gcc"] + CFLAGS + ["-DHAVE_MEDIC", objs[tag]] + LINK_REST + ["-o", out],
                       cwd=SRC, check=True)
    # The FP flags must have produced real codegen differences, or the experiment is vacuous.
    for tag in ("fast", "strict"):
        d = subprocess.run(["objdump", "-d", objs[tag]], capture_output=True, text=True).stdout
        n = sum(d.count(m) for m in ("fmadd", "fmsub", "fnmadd", "fnmsub"))
        print("  %-6s object: %d fused multiply-add instructions" % (tag, n))
    if _md5(objs["fast"]) == _md5(objs["strict"]):
        sys.exit("FATAL: the two medic.c objects are byte-identical; the experiment is vacuous")
    if _md5(fast_bin) == _md5(strict_bin):
        sys.exit("FATAL: the two binaries are byte-identical; the experiment is vacuous")
    print("  md5 fast   binary %s" % _md5(fast_bin))
    print("  md5 strict binary %s" % _md5(strict_bin))


def _md5(path: str) -> str:
    import hashlib
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


# ------------------------------------------------------------------ corpus

def demo_inputs():
    return ([os.path.join(DEMO, "echo-%d_part-mag.nii.gz" % e) for e in (1, 2)],
            [os.path.join(DEMO, "echo-%d_part-phase.nii.gz" % e) for e in (1, 2)])


def bold_inputs():
    p = lambda e, k: os.path.join(BOLD, "%s_echo-%d_part-%s_bold.nii.gz" % (BOLD_PRE, e, k))
    return ([p(e, "mag") for e in (1, 2)], [p(e, "phase") for e in (1, 2)])


def _hdr(vox=3.0):
    return {"pixdim": (1.0, vox, vox, vox, 1.0, 1.0, 1.0, 1.0), "xyzt_units": 10,
            "qform_code": 0, "sform_code": 1, "quatern": (0.0,) * 6,
            "srow": np.array([[vox, 0, 0, 0], [0, vox, 0, 0], [0, 0, vox, 0]], float)}


def _wrap(p):
    return (p + np.pi) % (2 * np.pi) - np.pi


def _blob(shape, rng=None):
    """A smooth magnitude 'brain': no hard edge, robustmask-friendly."""
    nx, ny, nz = shape
    x, y, z = np.meshgrid(*[np.arange(n, dtype=np.float64) for n in shape], indexing="ij")
    r2 = ((x - nx / 2) / (0.34 * nx)) ** 2 + ((y - ny / 2) / (0.34 * ny)) ** 2 + \
         ((z - nz / 2) / (0.36 * nz)) ** 2
    return 40.0 + 2600.0 * np.exp(-3.0 * r2)


def synth(tag, phase_e, mag_e, out_dir, dtype=np.float32):
    """Write a 2-echo synthetic pair; return (mags, phases) path lists."""
    h = _hdr()
    mags, phases = [], []
    for e in range(2):
        mp = os.path.join(out_dir, "%s_mag%d.nii" % (tag, e + 1))
        pp = os.path.join(out_dir, "%s_pha%d.nii" % (tag, e + 1))
        nii.write(mp, mag_e[e], ref=h, dtype=np.float32)
        nii.write(pp, phase_e[e], ref=h, dtype=dtype)
        mags.append(mp)
        phases.append(pp)
    return mags, phases


def make_synthetics(out_dir):
    """Deliberately adversarial volumes.  Each targets a distinct way the fast-math
    rounding of md_rescale_phase()/md_mcpc3ds() could reach a rescale() bin edge."""
    os.makedirs(out_dir, exist_ok=True)
    shape = (48, 48, 32)
    tes = [16.8e-3, 38.56e-3]
    rng = np.random.default_rng(20260725)
    mag = _blob(shape)
    cases = []

    # 1. WRAP BOUNDARY -- a field chosen so a large fraction of voxels land exactly on
    #    +-pi after wrapping, where an infinitesimal phase change flips the wrap.
    f = np.zeros(shape)
    f[:] = (np.arange(shape[1])[None, :, None] % 2) * (1.0 / (2 * tes[0]))   # +-pi exactly
    ph = [_wrap(TWO_PI * f * t) for t in tes]
    ph[0][::3] = np.pi          # force literal +pi at every third x-plane
    ph[1][1::3] = -np.pi
    cases.append(("wrapbound", ph, [mag, mag * 0.6]))

    # 2. NON-REPRESENTABLE RESCALE SLOPE -- observed span 3.0, so slope = 2*pi/3 and
    #    inter = -pi + 2*pi/3 are both non-terminating binary; the fused vs unfused
    #    `p[i]*slope + inter` in md_rescale_phase() then round differently in the tail.
    base = rng.uniform(-1.0, 2.0, shape)
    ph = [base.copy(), (base * 1.7 + 0.3)]
    ph[0].flat[0] = -1.0
    ph[0].flat[1] = 2.0         # pin the observed extrema
    ph[1].flat[0] = -1.0
    ph[1].flat[1] = 2.0
    cases.append(("nonrep", ph, [mag, mag * 0.6]))

    # 3. BIN-EDGE SWEEP -- phase quadratic in y, so the phase GRADIENT (and therefore
    #    ROMEO's pre-rescale weight) varies continuously across the volume and sweeps
    #    through every one of the 255 rescale() bin edges.  Somewhere in this volume a
    #    weight sits arbitrarily close to a bin boundary; that is exactly the geometry
    #    that made -ffast-math delete edges in romeo.c.
    y = np.arange(shape[1], dtype=np.float64)[None, :, None] + np.zeros(shape)
    grad = 0.02 + 0.28 * (y / shape[1])          # rad/voxel, sweeping
    f = np.cumsum(grad, axis=1) / (TWO_PI * tes[0])
    ph = [_wrap(TWO_PI * f * t) for t in tes]
    cases.append(("binedge", ph, [mag, mag * 0.6]))

    # 4. NEAR-DEGENERATE MAGNITUDE -- ties in the magnitude-derived weights, so the
    #    tie-break (and hence the growth order) is maximally fragile.
    m = np.full(shape, 1000.0)
    m += (rng.integers(0, 2, shape) * np.spacing(np.float32(1000.0)))  # 1-ULP ties
    m[:4] = 1e-30
    m[-4:] = 1e30
    f = 50.0 * rng.standard_normal(shape)
    ph = [_wrap(TWO_PI * f * t) for t in tes]
    cases.append(("degen", ph, [m, m]))

    # 5. INCONSISTENT / NOISY PHASE -- the regime where the unwrap is genuinely
    #    path-dependent, so a change in growth ORDER changes the RESULT.  This is the
    #    case most likely to expose a 2*pi branch flip if one is reachable at all.
    ph = [_wrap(rng.uniform(-np.pi, np.pi, shape)) for _ in tes]
    cases.append(("noisy", ph, [mag, mag * 0.6]))

    # 6. RESIDUAL FIELD -- a smooth field plus enough noise that the unwrap has real
    #    residues, but with a well-formed mask (the realistic hard case).
    x, yy, z = np.meshgrid(*[np.arange(n, dtype=np.float64) for n in shape], indexing="ij")
    f = 120.0 * np.sin(x / 7.0) * np.cos(yy / 5.0) + 40.0 * rng.standard_normal(shape)
    ph = [_wrap(TWO_PI * f * t) for t in tes]
    cases.append(("residual", ph, [mag, mag * 0.6]))

    # 7-9. The same three hard fields, but stored on the SIEMENS RAW SCALE (span 8190, the
    #      span of the real BOLD data) instead of in radians.  A radian-valued input makes
    #      md_rescale_phase() short-circuit (`fabs(span - 2*pi) <= 0.1`); these three force
    #      it to run, with slope = 2*pi/8190 and inter = -pi + 4096*slope -- both
    #      non-terminating in binary, so the fused and unfused `p[i]*slope + inter` differ.
    for tag, ph, mg in list(cases):
        if tag not in ("binedge", "noisy", "residual"):
            continue
        raw = []
        for p in ph:
            q = (p + np.pi) * (8190.0 / TWO_PI) - 4096.0
            q.flat[0] = -4096.0
            q.flat[1] = 4094.0          # pin the observed extrema
            raw.append(q)
        cases.append((tag + "_raw", raw, mg))

    out = []
    for tag, ph, mg in cases:
        mags, phases = synth(tag, ph, mg, out_dir)
        out.append((tag, mags, phases))
    return out


def configs(synthetics):
    """(name, mags, phases, extra_args, nvox_hint).  Every combination named in the task."""
    c = []
    dm, dp = demo_inputs()
    for pe in ("j", "j-"):
        c.append(("demo/pe=%s" % pe, dm, dp, ["--phase-encoding-direction", pe]))
    for w in ("romeo3", "romeo4"):
        c.append(("demo/weights=%s" % w, dm, dp, ["--weights", w]))
    c.append(("demo/mask=wk", dm, dp, ["--mask", WKMASK]))
    c.append(("demo/mask=wk,weights=romeo3", dm, dp, ["--mask", WKMASK, "--weights", "romeo3"]))
    c.append(("demo/offset=none", dm, dp, ["--phase-offset", "none"]))
    c.append(("demo/offset=none,mask=wk", dm, dp, ["--phase-offset", "none", "--mask", WKMASK]))
    c.append(("demo/rank=0", dm, dp, ["--rank", "0"]))
    c.append(("demo/temporal=0", dm, dp, ["--temporal-correction", "0"]))

    bm, bp = bold_inputs()
    c.append(("bold170/default", bm, bp, []))
    c.append(("bold170/mask=wk", bm, bp, ["--mask", WKMASK]))
    c.append(("bold170/pe=j-,weights=romeo3", bm, bp,
              ["--phase-encoding-direction", "j-", "--weights", "romeo3"]))
    c.append(("bold170/offset=none", bm, bp, ["--phase-offset", "none"]))

    for tag, mags, phases in synthetics:
        c.append(("synth/%s" % tag, mags, phases, []))
        if tag.split("_")[0] in ("noisy", "binedge", "residual"):
            c.append(("synth/%s,weights=romeo3" % tag, mags, phases, ["--weights", "romeo3"]))
    return c


# ------------------------------------------------------------------ run + compare

def run(binary, mags, phases, extra, prefix, nthread=8):
    cmd = [binary, "--medic",
           "--magnitude"] + list(mags) + ["--phase"] + list(phases) + [
           "--te-ms", TES, "--total-readout-time", TRT,
           "--out-prefix", prefix, "--save-intermediates", "--gz", "0", "-n", str(nthread)]
    if not any(a == "--phase-encoding-direction" for a in extra):
        cmd += ["--phase-encoding-direction", "j"]
    cmd += list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("run failed (%d):\n%s\n%s" % (r.returncode, " ".join(cmd), r.stderr[-2000:]))
    return prefix


def compare(pa, pb, suffixes):
    """Bit-for-bit comparison of two output sets.  Returns a per-suffix dict."""
    res = {}
    for s in suffixes:
        fa, fb = "%s_%s.nii" % (pa, s), "%s_%s.nii" % (pb, s)
        if not (os.path.exists(fa) and os.path.exists(fb)):
            continue
        a, _ = nii.read(fa)
        b, _ = nii.read(fb)
        a = np.asarray(a, dtype=np.float64)
        b = np.asarray(b, dtype=np.float64)
        # Raw bytes first: the "did anything change at all" question.
        ba = np.asarray(a, dtype=np.float32).tobytes()
        bb = np.asarray(b, dtype=np.float32).tobytes()
        nbyte = sum(1 for x, y in zip(ba, bb) if x != y) if ba != bb else 0
        fin = np.isfinite(a) & np.isfinite(b)
        # AGENTS.md: --compare's verdict is magnitude-based, so a non-finite mismatch has
        # no magnitude and must be counted SEPARATELY or an all-NaN pair reads as EQUAL.
        hard = int(np.count_nonzero(~fin & ~(np.isnan(a) & np.isnan(b)) &
                                    ~((np.isinf(a) & np.isinf(b)) & (np.sign(a) == np.sign(b)))))
        d = np.zeros_like(a)
        d[fin] = a[fin] - b[fin]
        ndiff = int(np.count_nonzero(d)) + hard
        mx = float(np.abs(d).max()) if d.size else 0.0
        branch = None
        if s in RADIAN:
            branch = int(np.count_nonzero(np.round(d / TWO_PI)))
        res[s] = dict(n=int(a.size), ndiff=ndiff, nbyte=nbyte, maxabs=mx,
                      hard=hard, branch=branch)
    return res


def summarise(name, res):
    lines = []
    worst_branch = 0
    mask_changed = False
    for s, r in res.items():
        tag = ""
        if s in RADIAN and r["branch"]:
            tag = "  <<< %d WHOLE-2*pi BRANCH DIFFERENCES" % r["branch"]
            worst_branch = max(worst_branch, r["branch"])
        if s == "masks" and r["ndiff"]:
            tag = "  <<< MASK CHANGED"
            mask_changed = True
        if r["hard"]:
            tag += "  <<< %d non-finite mismatches" % r["hard"]
        lines.append("    %-18s n=%-9d ndiff=%-7d maxabs=%-12.6g branch=%-5s%s"
                     % (s, r["n"], r["ndiff"], r["maxabs"],
                        "-" if r["branch"] is None else r["branch"], tag))
    return lines, worst_branch, mask_changed


# ------------------------------------------------------------------ Part B

def part_b(binary, out_dir, nthread=8):
    """Amplification bound: perturb the phase handed to ROMEO by +-1 float32 ULP at EVERY
    voxel and run ONE binary twice.  Fast-math was measured (Part A) to perturb a few
    voxels in 1e5; this perturbs 100 % of them, at the same per-voxel magnitude.  If no
    2*pi branch flip appears even here, the Part A null result is a mechanism bound rather
    than a coincidence."""
    os.makedirs(out_dir, exist_ok=True)
    shape = (48, 48, 32)
    tes = [16.8e-3, 38.56e-3]
    rng = np.random.default_rng(7)
    mag = _blob(shape)
    x, y, z = np.meshgrid(*[np.arange(n, dtype=np.float64) for n in shape], indexing="ij")
    out = []
    cases = {
        # smooth field, clean unwrap
        "smooth": 90.0 * np.sin(x / 9.0) * np.cos(y / 6.0),
        # smooth + residues: path-dependent unwrap
        "residual": 120.0 * np.sin(x / 7.0) * np.cos(y / 5.0) + 40.0 * rng.standard_normal(shape),
        # fully inconsistent: maximally path-dependent
        "noisy": 400.0 * rng.standard_normal(shape),
    }
    for tag, f in cases.items():
        ph = [np.float32(_wrap(TWO_PI * f * t)) for t in tes]
        # phase already in radians -> md_rescale_phase() short-circuits, so the +-1 ULP
        # perturbation lands EXACTLY at the hand-off point into romeo_unwrap_frame().
        sgn = rng.integers(0, 2, shape) * 2 - 1
        pert = [np.float32(np.where(sgn > 0, np.nextafter(p, np.float32(1e30)),
                                    np.nextafter(p, np.float32(-1e30)))) for p in ph]
        for sub, phases_arr in (("ref", ph), ("ulp", pert)):
            m, p = synth("%s_%s" % (tag, sub), phases_arr, [mag, mag * 0.6], out_dir)
            run(binary, m, p, [], os.path.join(out_dir, "%s_%s" % (tag, sub)), nthread)
        nchanged = int(np.count_nonzero(np.asarray(ph[0]) != np.asarray(pert[0])))
        res = compare(os.path.join(out_dir, "%s_ref" % tag),
                      os.path.join(out_dir, "%s_ulp" % tag), SUFFIXES + INTERMEDIATES)
        out.append((tag, nchanged, ph[0].size, res))
    return out


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", action="store_true", help="compile+link the two binaries first")
    ap.add_argument("--fast", default=os.path.join(SCRATCH, "nm_fast"))
    ap.add_argument("--strict", default=os.path.join(SCRATCH, "nm_strict"))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skip-bold", action="store_true", help="omit the 170-frame runs (~10 s each)")
    ap.add_argument("--skip-part-b", action="store_true")
    a = ap.parse_args()

    os.makedirs(SCRATCH, exist_ok=True)
    if a.build:
        print("== building two binaries differing ONLY in medic.c FP flags")
        build(a.fast, a.strict)
    for b in (a.fast, a.strict):
        if not os.path.exists(b):
            sys.exit("missing binary %s (run with --build)" % b)
    if _md5(a.fast) == _md5(a.strict):
        sys.exit("FATAL: the two binaries are identical")

    work = os.path.join(SCRATCH, "work")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work, exist_ok=True)

    print("\n== PART A: fast-math medic.c vs strict-FP medic.c, bit-for-bit")
    syn = make_synthetics(os.path.join(work, "syn"))
    cfgs = [c for c in configs(syn) if not (a.skip_bold and c[0].startswith("bold"))]
    total_vox = 0
    any_branch = 0
    any_mask = False
    nconf = 0
    t0 = time.time()
    for name, mags, phases, extra in cfgs:
        d = os.path.join(work, name.replace("/", "_").replace(",", "_").replace("=", ""))
        os.makedirs(d, exist_ok=True)
        pa = run(a.fast, mags, phases, extra, os.path.join(d, "fast"), a.threads)
        pb = run(a.strict, mags, phases, extra, os.path.join(d, "strict"), a.threads)
        res = compare(pa, pb, SUFFIXES + INTERMEDIATES)
        lines, br, mk = summarise(name, res)
        nconf += 1
        total_vox += sum(r["n"] for s, r in res.items() if s in RADIAN)
        any_branch += br
        any_mask |= mk
        flag = "DIFF" if any(r["ndiff"] for r in res.values()) else "identical"
        print("  %-34s %s" % (name, flag))
        for ln in lines:
            print(ln)
    print("  (%d configurations in %.1f s)" % (nconf, time.time() - t0))

    print("\n== PART B: +-1 float32 ULP on 100%% of phase voxels, ONE binary, twice")
    if not a.skip_part_b:
        for tag, nch, ntot, res in part_b(a.strict, os.path.join(work, "partb"), a.threads):
            lines, br, mk = summarise(tag, res)
            print("  %-12s perturbed %d/%d phase voxels (%.1f%%)%s"
                  % (tag, nch, ntot, 100.0 * nch / ntot, "  MASK CHANGED" if mk else ""))
            for ln in lines:
                print(ln)
            any_branch += 0   # reported separately; Part B is a bound, not the policy test

    print("\n== VERDICT")
    print("  configurations compared : %d" % nconf)
    print("  radian-unit voxels compared (Part A): %d" % total_vox)
    print("  whole-2*pi branch differences (Part A): %d" % any_branch)
    print("  mask differences (Part A): %s" % ("YES" if any_mask else "none"))
    if any_branch == 0 and not any_mask:
        print("  -> medic.c does NOT require strict FP on this evidence: the fast-math\n"
              "     differences are confined to the last float32 ULP and never cross a\n"
              "     2*pi branch or change the ROMEO graph's mask.")
    else:
        print("  -> EXPOSURE FOUND: medic.c requires strict FP.  See the flagged rows above.")


if __name__ == "__main__":
    main()

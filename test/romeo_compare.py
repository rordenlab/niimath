#!/usr/bin/env python3
"""Compare niimath's -romeo parity dumps against the pinned Julia oracle (romeo_plan.md §6).

    test/romeo_compare.py <niimath-binary> [--ref test/romeo_ref] [--case e0 ...] [--verbose]

The oracle (test/romeo_oracle.sh) writes raw little-endian arrays named "<case>_<what>"; niimath
writes the same arrays as "c_<what>" into the directory given to the hidden -romeo-dump option.
Raw binaries are used on both sides so the comparison is value-exact and never depends on NIfTI
container bytes.

Three tolerance tiers, used deliberately:
  * exact bytes  - integer-valued intermediates: UInt8 weights, mask stages, region labels
  * 1e-6         - the quality map and other pure float32 reductions
  * 1e-4         - the unwrapped phase (far below the 6.28 a single wrap error would produce);
                   the max |diff| is always reported, and a wrap-count difference always fails.

stdlib only, no numpy: the arrays are a few hundred thousand elements and `array`/`struct` are
fast enough, and this must run anywhere the repository's other test scripts do.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from array import array

HERE = os.path.dirname(os.path.abspath(__file__))
TWO_PI = 6.283185307179586


def read_raw(path: str, typecode: str) -> array:
    a = array(typecode)
    with open(path, "rb") as f:
        data = f.read()
    a.frombytes(data)
    if sys.byteorder != "little":
        a.byteswap()
    return a


class Result:
    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str]] = []
        self.failed = 0

    def add(self, name: str, ok: bool, detail: str) -> None:
        self.rows.append(("PASS" if ok else "FAIL", name, detail))
        if not ok:
            self.failed += 1

    def report(self) -> int:
        width = max((len(r[1]) for r in self.rows), default=10)
        for status, name, detail in self.rows:
            print(f"  {status}  {name:<{width}}  {detail}")
        return self.failed


def cmp_exact(res: Result, label: str, ref: str, got: str) -> None:
    if not os.path.exists(ref):
        res.add(label, False, f"missing oracle file {os.path.basename(ref)}")
        return
    if not os.path.exists(got):
        res.add(label, False, f"missing niimath dump {os.path.basename(got)}")
        return
    with open(ref, "rb") as f:
        a = f.read()
    with open(got, "rb") as f:
        b = f.read()
    if a == b:
        res.add(label, True, f"byte-identical ({len(a)} bytes)")
        return
    if len(a) != len(b):
        res.add(label, False, f"length {len(b)} != oracle {len(a)}")
        return
    ndiff = sum(1 for x, y in zip(a, b) if x != y)
    first = next(i for i, (x, y) in enumerate(zip(a, b)) if x != y)
    res.add(label, False, f"{ndiff}/{len(a)} bytes differ, first at {first} (oracle {a[first]} vs {b[first]})")


def cmp_float(res: Result, label: str, ref: str, got: str, tol: float, wrap_check: bool = False,
              relative: bool = False) -> None:
    if not os.path.exists(ref) or not os.path.exists(got):
        missing = os.path.basename(ref if not os.path.exists(ref) else got)
        res.add(label, False, f"missing {missing}")
        return
    a = read_raw(ref, "f")
    b = read_raw(got, "f")
    if len(a) == 1 and len(b) > 1:
        # Upstream's B0 SNR collapses to ONE value when no magnitude is supplied: the app
        # substitutes a voxel-independent exp(-TE/20) decay, so `sum(mag.*w)/sum(w)` has no
        # spatial extent and ROMEO writes a 1x1x1 image. niimath writes the same constant across
        # the working grid instead (every other side output lives on that grid). Values identical,
        # shape deliberately different - compare the constant against every voxel.
        a = array("f", [a[0]] * len(b))
    if len(a) != len(b):
        res.add(label, False, f"length {len(b)} != oracle {len(a)}")
        return
    maxdiff = 0.0
    argmax = -1
    nonfinite = 0
    maxresid = 0.0
    maxwrap = 0
    for i, (x, y) in enumerate(zip(a, b)):
        fx, fy = math.isfinite(x), math.isfinite(y)
        if not fx or not fy:
            # a non-finite difference has no magnitude; count it separately so an all-NaN pair
            # can never read as "equal" (the --compare gotcha, see AGENTS.md)
            if fx != fy:
                nonfinite += 1
            elif not (math.isnan(x) and math.isnan(y)) and x != y:
                nonfinite += 1
            continue
        d = abs(x - y)
        if d > maxdiff:
            maxdiff, argmax = d, i
        if wrap_check:
            k = round((x - y) / TWO_PI)
            if k != 0:
                maxwrap = max(maxwrap, abs(int(k)))
            r = abs((x - y) - k * TWO_PI)
            maxresid = max(maxresid, r)
    scale = 1.0
    if relative:
        scale = max(1.0, max((abs(x) for x in a if math.isfinite(x)), default=1.0))
    ok = maxdiff <= tol * scale and nonfinite == 0 and maxwrap == 0
    detail = (f"max|diff|={maxdiff:.3e} at {argmax}, tol={tol:g}"
              + (f"*{scale:.3g}" if relative else "") + f", nonfinite-mismatch={nonfinite}")
    if wrap_check:
        detail += f", wraps={maxwrap}, residual={maxresid:.3e}"
    res.add(label, ok, detail)


def cmp_seed_scalars(res: Result, tag: str, ref: str, cmanifest: str) -> None:
    """Seed voxel index and new_seed_thresh (plan M5 gate) — compared directly, not merely
    implied by the byte-exact weights they are derived from."""
    label = f"{tag}: seed scalars"
    orapath = os.path.join(ref, "manifest.txt")
    if not os.path.exists(orapath) or not os.path.exists(cmanifest):
        res.add(label, False, "missing manifest")
        return
    want = {}
    for line in open(orapath):
        parts = line.split()
        if len(parts) >= 3 and parts[0] == "scalar" and parts[1].startswith(tag + "_"):
            key = parts[1][len(tag) + 1:]
            if key in ("seed_index", "seed_w1", "seed_w2", "seed_w3", "new_seed_thresh"):
                want[key] = float(parts[2])
    got = {}
    for line in open(cmanifest):
        parts = line.split()
        if len(parts) == 2 and parts[0] in ("seed_index", "seed_w1", "seed_w2", "seed_w3", "new_seed_thresh"):
            got[parts[0]] = float(parts[1])
    if not want:
        res.add(label, False, "oracle has no seed scalars")
        return
    bad = [k for k in want if k not in got or got[k] != want[k]]
    if bad:
        res.add(label, False, "mismatch: " + ", ".join("%s oracle=%g c=%s" % (k, want[k], got.get(k)) for k in bad))
    else:
        res.add(label, True, "seed=%d new_seed_thresh=%g exact" % (int(want["seed_index"]), want.get("new_seed_thresh", 0)))


def cmp_ulp(res: Result, label: str, ref: str, got: str, max_ulp: int) -> None:
    """Compare two Float64 arrays in ULPs (plan M2: the pre-rescale weights must agree to <=4 ULP).
    Bit patterns are mapped to a monotone signed ordering so the distance is exact."""
    if not os.path.exists(ref) or not os.path.exists(got):
        res.add(label, False, "missing " + os.path.basename(ref if not os.path.exists(ref) else got))
        return
    a = read_raw(ref, "d")
    b = read_raw(got, "d")
    if len(a) != len(b):
        res.add(label, False, f"length {len(b)} != oracle {len(a)}")
        return

    def ordered(x: float) -> int:
        (bits,) = struct.unpack("<q", struct.pack("<d", x))
        return bits if bits >= 0 else -(bits & 0x7FFFFFFFFFFFFFFF) - (1 << 63)

    worst = 0
    argworst = -1
    nonfinite = 0
    for i, (x, y) in enumerate(zip(a, b)):
        if not (math.isfinite(x) and math.isfinite(y)):
            if not (math.isnan(x) and math.isnan(y)) and x != y:
                nonfinite += 1
            continue
        d = abs(ordered(x) - ordered(y))
        if d > worst:
            worst, argworst = d, i
    ok = worst <= max_ulp and nonfinite == 0
    res.add(label, ok, f"max {worst} ULP at {argworst} (limit {max_ulp}), nonfinite-mismatch={nonfinite}")


def cmp_b0_snr(res: Result, label: str, ref: str, got: str) -> None:
    """Exact, but tolerating the deliberate SHAPE difference: without a magnitude upstream's SNR
    collapses to one value (the substituted exp(-TE/20) decay is voxel-independent), while
    niimath writes that constant across the working grid."""
    if not os.path.exists(ref) or not os.path.exists(got):
        res.add(label, False, "missing " + os.path.basename(ref if not os.path.exists(ref) else got))
        return
    a = read_raw(ref, "f")
    b = read_raw(got, "f")
    if len(a) == 1 and len(b) > 1:
        bad = sum(1 for x in b if x != a[0])
        res.add(label, bad == 0, "constant %g broadcast over %d voxels, %d mismatches" % (a[0], len(b), bad))
        return
    cmp_exact(res, label, ref, got)


def run_niimath(binary: str, args: list[str], cwd: str) -> tuple[int, str]:
    p = subprocess.run([binary] + args, cwd=cwd, capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr)


CASES = {
    # case tag -> (phase, magnitude|None, extra -romeo args)
    "e0": ("phase0.nii.gz", "mag0.nii.gz", ["-t", "16.8"]),
    "e1": ("phase1.nii.gz", "mag1.nii.gz", ["-t", "38.56"]),
    "me": ("phase.nii", "mag.nii", ["-t", "[16.8,38.56]"]),
    "e0n": ("phase0.nii.gz", None, ["-t", "16.8"]),
}

SYNTH = ["small_a", "small_b", "line_x", "plane_xy", "degen"]

WEIGHT_SELECTIONS = ["romeo3", "romeo2", "romeo4", "romeo6",
                     "100000", "110000", "101000", "000100", "111111", "010000"]
WEIGHT_TAGS = {"romeo3": "romeo3", "romeo2": "romeo2", "romeo4": "romeo4", "romeo6": "romeo6",
               "100000": "f100000", "110000": "f110000", "101000": "f101000",
               "000100": "f000100", "111111": "f111111", "010000": "f010000"}


def check_case(binary: str, ref: str, tag: str, datadir: str, phase: str, mag: str | None,
               extra: list[str], res: Result, weights_all: bool) -> None:
    tmp = tempfile.mkdtemp(prefix="romeo_cmp_")
    try:
        out = os.path.join(tmp, "out")
        args = [phase, "-romeo", mag if mag else "none"] + extra + ["-romeo-dump", tmp, out]
        rc, log = run_niimath(binary, args, datadir)
        if rc != 0:
            res.add(f"{tag}: run", False, f"exit {rc}: {log.strip()[:200]}")
            return
        res.add(f"{tag}: run", True, "ok")

        def R(name: str) -> str:
            return os.path.join(ref, f"{tag}_{name}")

        def C(name: str) -> str:
            return os.path.join(tmp, f"c_{name}")

        cmp_exact(res, f"{tag}: phase rescale", R("phase_rescaled.f32"), C("phase_rescaled.f32"))
        cmp_exact(res, f"{tag}: weights (u8)", R("weights_romeo3.u8") if mag else R("weights_romeo4.u8"),
                  C("weights.u8"))
        if mag:
            for st in ("mask_s1_thresh.u8", "mask_s2_smooth1.u8", "mask_s3_fill.u8", "mask_s4_final.u8"):
                cmp_exact(res, f"{tag}: {st.split('.')[0]}", R(st), C(st))
        cmp_exact(res, f"{tag}: visited", R("visited.u8"), C("visited.u8"))
        cmp_seed_scalars(res, tag, ref, os.path.join(tmp, "c_manifest.txt"))
        cmp_float(res, f"{tag}: qmap", R("qmap.f32"), C("qmap.f32"), 1e-6)
        for qi in range(1, 7):
            cmp_float(res, f"{tag}: qmap_{qi}", R(f"qmap_{qi}.f32"), C(f"qmap_{qi}.f32"), 1e-6)
        cmp_ulp(res, f"{tag}: weights pre-rescale", R("weights_prerescale.f64"),
                C("weights_prerescale.f64"), 4)
        cmp_float(res, f"{tag}: unwrapped", R("unwrapped.f32"), C("unwrapped.f32"), 1e-4, wrap_check=True)

        # B0 (plan M8): every weighting mode, with and without a magnitude.
        for wmode in ("phase_snr", "phase_var", "average", "TEs", "mag", "simulated_mag"):
            if not os.path.exists(R(f"b0_{wmode}.f32")):
                continue
            tmpb = tempfile.mkdtemp(prefix="romeo_b0_")
            try:
                ab = [phase, "-romeo", mag if mag else "none"] + extra + \
                     ["-B", "-B0-phase-weighting", wmode, "-romeo-dump", tmpb, os.path.join(tmpb, "o")]
                rcb, logb = run_niimath(binary, ab, datadir)
                if rcb != 0:
                    res.add(f"{tag}: -B {wmode}", False, f"exit {rcb}: {logb.strip()[:150]}")
                    continue
                cmp_exact(res, f"{tag}: -B {wmode}", R(f"b0_{wmode}.f32"),
                          os.path.join(tmpb, "c_b0.f32"))
                cmp_b0_snr(res, f"{tag}: -B {wmode} snr", R(f"b0_snr_{wmode}.f32"),
                           os.path.join(tmpb, "c_b0_snr.f32"))
            finally:
                shutil.rmtree(tmpb, ignore_errors=True)

        # Variant runs. The oracle dumps these arrays (M5/M6 exit criteria); each needs its own
        # niimath invocation because they change the unwrapping, not just an output.
        variants = [("correct-global", ["-g"], "unwrapped_correctglobal.f32")]
        if tag in ("me", "me_a"):
            variants += [
                ("individual", ["-i"], "unwrapped_individual.f32"),
                ("tuu 0.5", ["-temporal-uncertain-unwrapping", "0.5"], "unwrapped_tuu05.f32"),
                ("template 2", ["-template", "2"], "unwrapped_template2.f32"),
            ]
        for label, extra_args, refname in variants:
            if not os.path.exists(R(refname)):
                continue
            tmpv = tempfile.mkdtemp(prefix="romeo_v_")
            try:
                av = [phase, "-romeo", mag if mag else "none"] + extra + extra_args + \
                     ["-romeo-dump", tmpv, os.path.join(tmpv, "o")]
                rcv, logv = run_niimath(binary, av, datadir)
                if rcv != 0:
                    res.add(f"{tag}: {label}", False, f"exit {rcv}: {logv.strip()[:150]}")
                    continue
                cmp_float(res, f"{tag}: {label}", R(refname),
                          os.path.join(tmpv, "c_unwrapped.f32"), 1e-4, wrap_check=True)
            finally:
                shutil.rmtree(tmpv, ignore_errors=True)

        if weights_all:
            for sel in WEIGHT_SELECTIONS:
                tmp2 = tempfile.mkdtemp(prefix="romeo_w_")
                try:
                    a2 = [phase, "-romeo", mag if mag else "none"] + extra + \
                         ["-w", sel, "-romeo-dump", tmp2, os.path.join(tmp2, "o")]
                    rc2, log2 = run_niimath(binary, a2, datadir)
                    if rc2 != 0:
                        res.add(f"{tag}: -w {sel}", False, f"exit {rc2}: {log2.strip()[:150]}")
                        continue
                    cmp_exact(res, f"{tag}: -w {sel}",
                              os.path.join(ref, f"{tag}_weights_{WEIGHT_TAGS[sel]}.u8"),
                              os.path.join(tmp2, "c_weights.u8"))
                finally:
                    shutil.rmtree(tmp2, ignore_errors=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def write_f32_nifti(path: str, dims: tuple[int, int, int], data: list[float]) -> None:
    """Minimal NIfTI-1 (n+1) float32 writer for the synthetic property fixtures."""
    hdr = bytearray(348)
    struct.pack_into("<i", hdr, 0, 348)
    struct.pack_into("<8h", hdr, 40, 3, dims[0], dims[1], dims[2], 1, 1, 1, 1)
    struct.pack_into("<h", hdr, 70, 16)   # DT_FLOAT32
    struct.pack_into("<h", hdr, 72, 32)   # bitpix
    struct.pack_into("<8f", hdr, 76, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    struct.pack_into("<f", hdr, 108, 352.0)
    struct.pack_into("<f", hdr, 112, 1.0)
    hdr[123] = 2                          # xyz_units = mm
    struct.pack_into("<h", hdr, 252, 1)
    struct.pack_into("<h", hdr, 254, 1)
    struct.pack_into("<4f", hdr, 280, 1.0, 0.0, 0.0, 0.0)
    struct.pack_into("<4f", hdr, 296, 0.0, 1.0, 0.0, 0.0)
    struct.pack_into("<4f", hdr, 312, 0.0, 0.0, 1.0, 0.0)
    hdr[344:348] = b"n+1\0"
    with open(path, "wb") as f:
        f.write(bytes(hdr) + b"\0\0\0\0" + struct.pack("<%df" % len(data), *data))


def check_primitives(binary: str, ref: str, res: Result) -> None:
    """Byte-compare the low-level numeric primitives (rem2pi Float64/Float32, gamma, rescale,
    both unwrapvoxel widths) against the oracle on a fixed input table shared by both sides.
    The table deliberately includes |x| >= 2^20*pi/2 so the Payne-Hanek branch is covered — no
    real phase image reaches it, but -no-phase-rescale lets a caller supply values that do."""
    if not os.path.exists(os.path.join(ref, "prim_rem2pi64.f64")):
        res.add("primitives", False, "oracle has no primitive tables (regenerate with romeo_oracle.sh)")
        return
    tmp = tempfile.mkdtemp(prefix="romeo_prim_")
    try:
        fix = os.path.join(ref, "fixtures")
        phase = os.path.join(fix, "small_a_phase.nii")
        if not os.path.exists(phase):
            res.add("primitives", False, "missing synthetic fixture")
            return
        rc, log = run_niimath(binary, [phase, "-romeo", os.path.join(fix, "small_a_mag.nii"),
                                       "-t", "16.8", "-romeo-dump", tmp,
                                       os.path.join(tmp, "o")], tmp)
        if rc != 0:
            res.add("primitives", False, f"exit {rc}: {log.strip()[:150]}")
            return
        for name in ("rem2pi64.f64", "rem2pi32_gamma.f32", "rescale.u8", "unwrapvoxel.f32"):
            cmp_exact(res, "primitive %s" % name.split(".")[0],
                      os.path.join(ref, "prim_" + name), os.path.join(tmp, "c_prim_" + name))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def check_features_property(binary: str, res: Result) -> None:
    """ROMEO.jl test/features.jl, reimplemented against the C port (plan M5 exit criterion).

        for l in 7:5:20, offset in 2pi*[0,2,-1,10,-50], mag in [nothing, 10:1, 1:10]
            unwrap(rem2pi(phase_uw) + offset; mag, correctglobal=true) ~= phase_uw

    The library call takes the array directly, so the app-level behaviours must be switched off
    to match it: -no-phase-rescale (readphase would rescale this input), -k nomask (the library
    default mask is `trues`; robustmask also needs >=20 voxels), and -w romeo4, which is the
    LIBRARY's :romeo flag set (1-4) — the CLI would otherwise resolve bare `romeo` to romeo3.
    """
    tmp = tempfile.mkdtemp(prefix="romeo_feat_")
    npass = 0
    try:
        for length in range(7, 21, 5):
            uw = [-TWO_PI + 2.0 * TWO_PI * i / (length - 1) for i in range(length)]
            for offset in [TWO_PI * k for k in (0, 2, -1, 10, -50)]:
                wrapped = [x - TWO_PI * round(x / TWO_PI) + offset for x in uw]
                phase_path = os.path.join(tmp, "p.nii")
                write_f32_nifti(phase_path, (length, 1, 1), wrapped)
                mags = [None,
                        [10.0 - 9.0 * i / (length - 1) for i in range(length)],
                        [1.0 + 9.0 * i / (length - 1) for i in range(length)]]
                for mi, magvals in enumerate(mags):
                    if magvals is None:
                        magarg = "none"
                    else:
                        magarg = os.path.join(tmp, "m.nii")
                        write_f32_nifti(magarg, (length, 1, 1), magvals)
                    out = os.path.join(tmp, "o.nii")
                    label = "features l=%d off=%dpi mag=%d" % (length, round(offset / math.pi), mi)
                    rc, log = run_niimath(binary, [phase_path, "-gz", "0", "-romeo", magarg,
                                                   "-w", "romeo4", "-k", "nomask",
                                                   "-no-phase-rescale", "-g", out], tmp)
                    if rc != 0:
                        res.add(label, False, "exit %d: %s" % (rc, log.strip()[:120]))
                        continue
                    with open(out, "rb") as f:
                        blob = f.read()
                    off = int(struct.unpack_from("<f", blob, 108)[0])
                    got = list(struct.unpack_from("<%df" % length, blob, off))
                    worst = max(abs(a - b) for a, b in zip(got, uw))
                    if worst > 1e-4:
                        res.add(label, False, "max|diff| from ground truth = %.3e" % worst)
                    else:
                        npass += 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    expected = 3 * 5 * 3  # l in 7:5:20 x 5 offsets x 3 magnitude variants
    res.add("features.jl property (%d/%d cases)" % (npass, expected), npass == expected,
            "unwrap(...; correctglobal) == ground truth")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("--ref", default=os.path.join(HERE, "romeo_ref"))
    ap.add_argument("--data", default=os.environ.get("ROMEO_JL_DATA", "/Users/chris/src/ROMEO.jl/romeo"))
    ap.add_argument("--case", action="append", default=None)
    ap.add_argument("--weights-all", action="store_true", help="also check every -w selection")
    args = ap.parse_args()

    ref = os.path.abspath(args.ref)
    if not os.path.isdir(ref):
        print(f"SKIP: oracle directory {ref} not found — run test/romeo_oracle.sh first")
        return 77
    binary = os.path.abspath(args.binary)
    res = Result()

    wanted = args.case
    for tag, (phase, mag, extra) in CASES.items():
        if wanted and tag not in wanted:
            continue
        if not os.path.isdir(args.data):
            print(f"SKIP: validation data {args.data} not found")
            break
        check_case(binary, ref, tag, args.data, phase, mag, extra, res, args.weights_all)

    fixdir = os.path.join(ref, "fixtures")
    if os.path.isdir(fixdir):
        for name in SYNTH:
            if wanted and name not in wanted:
                continue
            check_case(binary, ref, name, fixdir, f"{name}_phase.nii", f"{name}_mag.nii",
                       ["-t", "16.8"], res, args.weights_all)
            check_case(binary, ref, f"{name}_nomag", fixdir, f"{name}_phase.nii", None,
                       ["-t", "16.8"], res, args.weights_all)
        if not wanted or "me_a" in (wanted or []):
            check_case(binary, ref, "me_a", fixdir, "me_a_phase.nii", "me_a_mag.nii",
                       ["-t", "[16.8,38.56]"], res, args.weights_all)

    check_primitives(binary, ref, res)
    check_features_property(binary, res)

    failed = res.report()
    total = len(res.rows)
    print(f"\n{total - failed}/{total} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

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


def cmp_float(res: Result, label: str, ref: str, got: str, tol: float, wrap_check: bool = False) -> None:
    if not os.path.exists(ref) or not os.path.exists(got):
        missing = os.path.basename(ref if not os.path.exists(ref) else got)
        res.add(label, False, f"missing {missing}")
        return
    a = read_raw(ref, "f")
    b = read_raw(got, "f")
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
            if fx != fy or (fx and fy and x != y):
                nonfinite += 1
            elif not fx and not fy and not (math.isnan(x) and math.isnan(y)) and x != y:
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
    ok = maxdiff <= tol and nonfinite == 0 and maxwrap == 0
    detail = f"max|diff|={maxdiff:.3e} at {argmax}, tol={tol:g}, nonfinite-mismatch={nonfinite}"
    if wrap_check:
        detail += f", wraps={maxwrap}, residual={maxresid:.3e}"
    res.add(label, ok, detail)


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
        cmp_float(res, f"{tag}: qmap", R("qmap.f32"), C("qmap.f32"), 1e-6)
        for qi in range(1, 7):
            cmp_float(res, f"{tag}: qmap_{qi}", R(f"qmap_{qi}.f32"), C(f"qmap_{qi}.f32"), 1e-6)
        cmp_float(res, f"{tag}: unwrapped", R("unwrapped.f32"), C("unwrapped.f32"), 1e-4, wrap_check=True)

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


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("--ref", default=os.path.join(HERE, "romeo_ref"))
    ap.add_argument("--data", default="/Users/chris/src/ROMEO.jl/romeo")
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

    failed = res.report()
    total = len(res.rows)
    print(f"\n{total - failed}/{total} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

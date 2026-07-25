"""Head-to-head wall time / CPU time / peak RAM: `niimath --medic` vs `wk-medic`.

Workload is demo/run170.sh: the 170-frame two-echo BOLD run, estimate stage then apply stage.
Writes to a scratch directory -- demo/out170/ is the validation oracle and must not be clobbered.

    ~/src/warpkit/.venv/bin/python bench170.py [--ncpus 8] [--gz 1]

CAVEAT that cost a full bogus benchmark once: verify the binary really has OpenMP before
trusting any niimath number here --

    ./src/niimath <img> -p 8 -s 1 /tmp/x.nii    # must print "Using 8 threads"

A stale object left behind by a `make OMP=0 / ROMEO=0 / MEDIC=0` build silently produces a
serial binary and a meaningless comparison.

Build niimath against zlib-ng for the gzip paths (AGENTS.md): plain `make` links system zlib.

    make -C src ZLIBNG_ROOT=/abs/path/to/zlib-ng-build
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HOME = Path.home()
IN = HOME / "src/warpkit/inputs/sub-fm/ses-1/func"
BASE = "sub-fm_ses-1_task-rest_acq-2d2echo_run-02"
WK = HOME / "src/warpkit/.venv/bin"
DIMS = (76, 76, 46, 170)


def timed(cmd: list[str], tag: str, scratch: Path) -> dict:
    """Run under /usr/bin/time -l and return {wall, cpu, rss, footprint}."""
    tf = scratch / f"{tag}.time"
    with open(scratch / f"{tag}.log", "w") as out, open(tf, "w") as err:
        r = subprocess.run(["/usr/bin/time", "-l"] + cmd, stdout=out, stderr=err)
    # A failed command still produces a plausible-looking time; refuse to report it as a result.
    if r.returncode != 0:
        raise SystemExit(f"{tag}: command failed with exit {r.returncode}\n"
                         f"  {' '.join(cmd)}\n  see {scratch}/{tag}.log and {tf}")
    txt = tf.read_text()
    # macOS puts the number BEFORE the label: "  15.50 real  65.50 user  4.06 sys"
    def num(label: str, pat: str = r"([\d.]+)\s+%s") -> float:
        m = re.search(pat % label, txt)
        return float(m.group(1)) if m else float("nan")
    def bytes_of(label: str) -> float:
        m = re.search(r"(\d+)\s+%s" % label, txt)
        return int(m.group(1)) / 2**30 if m else float("nan")
    return {
        "wall": num("real"),
        "cpu": num("user"),
        "rss": bytes_of("maximum resident set size"),
        "peak": bytes_of("peak memory footprint"),
    }


def show(label: str, r: dict) -> None:
    par = r["cpu"] / r["wall"] if r["wall"] else 0.0
    print(f"  {label:<20s} wall {r['wall']:7.2f}s  cpu {r['cpu']:8.2f}s  par {par:5.1f}x  "
          f"maxRSS {r['rss']:6.2f} GB  peak {r['peak']:6.2f} GB")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ncpus", type=int, default=8)
    ap.add_argument("--gz", default="1", choices=["0", "1"], help="niimath output compression")
    ap.add_argument("--niimath", default=os.environ.get("NIIMATH", "./src/niimath"))
    ap.add_argument("--scratch", type=Path, default=Path("/tmp/medic_bench"))
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[2]
    os.chdir(root)
    if args.scratch.exists():
        shutil.rmtree(args.scratch)
    args.scratch.mkdir(parents=True)

    nvox = DIMS[0] * DIMS[1] * DIMS[2] * DIMS[3]
    print("=== workload ===")
    print(f"  {DIMS[0]}x{DIMS[1]}x{DIMS[2]} x {DIMS[3]} frames, 2 echoes, magnitude + phase")
    print(f"  one series float32: {nvox * 4 / 2**20:.0f} MiB     "
          f"all four float32: {4 * nvox * 4 / 2**30:.2f} GiB     float64: {4 * nvox * 8 / 2**30:.2f} GiB")

    mags = [str(IN / f"{BASE}_echo-{e}_part-mag_bold.nii.gz") for e in (1, 2)]
    phas = [str(IN / f"{BASE}_echo-{e}_part-phase_bold.nii.gz") for e in (1, 2)]
    meta = [str(IN / f"{BASE}_echo-{e}_part-phase_bold.json") for e in (1, 2)]

    print(f"\n=== ESTIMATE stage, {args.ncpus} threads ===")
    show("wk-medic", timed(
        [str(WK / "wk-medic"), "--magnitude", *mags, "--phase", *phas, "--metadata", *meta,
         "--out-prefix", str(args.scratch / "wk"), "-n", str(args.ncpus)], "wk", args.scratch))
    show(f"niimath --gz {args.gz}", timed(
        [args.niimath, "--medic", "--magnitude", *mags, "--phase", *phas,
         "--te-ms", "16.8,38.56", "--total-readout-time", "0.02025",
         "--phase-encoding-direction", "j", "--n-cpus", str(args.ncpus),
         "--gz", args.gz, "--out-prefix", str(args.scratch / "nm")], "nm", args.scratch))

    dmap = args.scratch / "wk_displacementmaps.nii"
    if not dmap.exists():
        print("  (wk-medic produced no displacement map; skipping the apply stage)")
        return 0
    print(f"\n=== APPLY stage (one echo, 170 frames), {args.ncpus} threads ===")
    show("wk-apply-warp", timed(
        [str(WK / "wk-apply-warp"), "--input", mags[0], "--transform", str(dmap),
         "--transform-type", "map", "--phase-encoding-axis", "j",
         "--output", str(args.scratch / "wk_e1_undistorted.nii.gz")], "wka", args.scratch))
    show("niimath -unwarp", timed(
        [args.niimath, mags[0], "-p", str(args.ncpus), "-gz", args.gz,
         "-unwarp", str(dmap), "j", str(args.scratch / "nm_e1_undistorted.nii")], "nma", args.scratch))

    print("\n=== outputs ===")
    for f in sorted(args.scratch.glob("*.nii*")):
        print(f"  {f.name:<44s} {f.stat().st_size / 2**20:8.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())

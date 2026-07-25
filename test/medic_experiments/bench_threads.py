"""Like-for-like: niimath --medic vs wk-medic, UNCOMPRESSED output, 1 thread vs N threads.

Both tools write uncompressed .nii, so this compares compute rather than gzip. Estimate stage
only (the 170-frame two-echo run from demo/run170.sh).

    ~/src/warpkit/.venv/bin/python bench_threads.py [--threads 1 8]

Verify OpenMP is really linked before trusting any niimath number:
    ./src/niimath <img> -p 8 -s 1 /tmp/x.nii     # must print "Using 8 threads"
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


def timed(cmd: list[str], tag: str, scratch: Path) -> dict:
    tf = scratch / f"{tag}.time"
    with open(scratch / f"{tag}.log", "w") as out, open(tf, "w") as err:
        r = subprocess.run(["/usr/bin/time", "-l"] + cmd, stdout=out, stderr=err)
    # A failed command still produces a plausible-looking time; refuse to report it as a result.
    if r.returncode != 0:
        raise SystemExit(f"{tag}: command failed with exit {r.returncode}\n"
                         f"  {' '.join(cmd)}\n  see {scratch}/{tag}.log and {tf}")
    txt = tf.read_text()

    def num(label: str) -> float:
        m = re.search(r"([\d.]+)\s+%s" % label, txt)
        return float(m.group(1)) if m else float("nan")

    def gb(label: str) -> float:
        m = re.search(r"(\d+)\s+%s" % label, txt)
        return int(m.group(1)) / 2**30 if m else float("nan")

    return {"wall": num("real"), "cpu": num("user"),
            "rss": gb("maximum resident set size"), "peak": gb("peak memory footprint")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, nargs="+", default=[1, 8])
    ap.add_argument("--niimath", default=os.environ.get("NIIMATH", "./src/niimath"))
    ap.add_argument("--scratch", type=Path, default=Path("/tmp/medic_thr"))
    args = ap.parse_args()

    os.chdir(Path(__file__).resolve().parents[2])
    if args.scratch.exists():
        shutil.rmtree(args.scratch)
    args.scratch.mkdir(parents=True)

    mags = [str(IN / f"{BASE}_echo-{e}_part-mag_bold.nii.gz") for e in (1, 2)]
    phas = [str(IN / f"{BASE}_echo-{e}_part-phase_bold.nii.gz") for e in (1, 2)]
    meta = [str(IN / f"{BASE}_echo-{e}_part-phase_bold.json") for e in (1, 2)]

    print("ESTIMATE stage, 170 frames x 2 echoes, UNCOMPRESSED output (like for like)\n")
    print(f"  {'tool':<18s} {'thr':>4s} {'wall':>9s} {'cpu':>10s} {'par':>6s} {'peak RAM':>10s}")
    rows = {}
    for n in args.threads:
        r = timed([str(WK / "wk-medic"), "--magnitude", *mags, "--phase", *phas,
                   "--metadata", *meta, "--out-prefix", str(args.scratch / f"wk{n}"),
                   "-n", str(n)], f"wk{n}", args.scratch)
        rows[("wk-medic", n)] = r
        print(f"  {'wk-medic':<18s} {n:>4d} {r['wall']:>8.2f}s {r['cpu']:>9.2f}s "
              f"{r['cpu'] / r['wall']:>5.1f}x {r['peak']:>9.2f} GB")
        r = timed([args.niimath, "--medic", "--magnitude", *mags, "--phase", *phas,
                   "--te-ms", "16.8,38.56", "--total-readout-time", "0.02025",
                   "--phase-encoding-direction", "j", "--n-cpus", str(n), "--gz", "0",
                   "--out-prefix", str(args.scratch / f"nm{n}")], f"nm{n}", args.scratch)
        rows[("niimath", n)] = r
        print(f"  {'niimath --medic':<18s} {n:>4d} {r['wall']:>8.2f}s {r['cpu']:>9.2f}s "
              f"{r['cpu'] / r['wall']:>5.1f}x {r['peak']:>9.2f} GB")
        print()

    if len(args.threads) > 1:
        lo, hi = min(args.threads), max(args.threads)
        print(f"  scaling {lo} -> {hi} threads:")
        for tool in ("wk-medic", "niimath"):
            a, b = rows[(tool, lo)]["wall"], rows[(tool, hi)]["wall"]
            print(f"    {tool:<18s} {a:6.2f}s -> {b:6.2f}s   speedup {a / b:.2f}x")
    return 0


if __name__ == "__main__":
    sys.exit(main())

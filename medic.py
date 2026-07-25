#!/usr/bin/env python3
"""Minimal BIDS front end for `niimath --medic`.

Discovers multi-echo magnitude/phase pairs in a BIDS `func/` directory (or a dataset root),
reads the acquisition parameters from the JSON sidecars, and calls niimath once per run to
estimate the field maps, then once per magnitude echo to undistort it.

Deliberately dependency-free: stdlib only (argparse, json, os, pathlib, re, shutil, subprocess,
sys). niimath performs all NIfTI I/O -- this file never opens an image.

    python3 medic.py /path/to/sub-XX/ses-Y/func --out-dir derivatives/medic

See test/medic_reference_manifest.md for the conventions `--medic` implements, and prior_art.md
for the MCPC-3D-S patent analysis.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# sub-XX[_ses-Y]..._echo-<N>_part-{mag,phase}_bold.nii[.gz|.zst]
ENTITY_RE = re.compile(
    r"^(?P<stem>.+?)_echo-(?P<echo>\d+)_part-(?P<part>mag|phase)_bold"
    r"(?P<ext>\.nii(?:\.gz|\.zst)?)$"
)


def find_runs(root: Path) -> dict[str, dict]:
    """Group NIfTIs into runs keyed by the stem with `echo` and `part` removed.

    Only those two entities are stripped, so different tasks/runs/acquisitions stay separate.
    `sbref` and non-`bold` suffixes never match ENTITY_RE and are ignored.
    """
    runs: dict[str, dict] = {}
    files = sorted(root.rglob("*_bold.nii*")) if root.is_dir() else []
    for f in files:
        m = ENTITY_RE.match(f.name)
        if not m:
            continue
        stem = m.group("stem")
        echo = int(m.group("echo"))
        part = m.group("part")
        run = runs.setdefault(stem, {"dir": f.parent, "mag": {}, "phase": {}})
        run[part][echo] = f
    return runs


def sidecar(img: Path) -> Path:
    """The exact JSON sidecar for an image. v1 does not implement BIDS inheritance."""
    name = img.name
    for ext in (".nii.gz", ".nii.zst", ".nii"):
        if name.endswith(ext):
            return img.with_name(name[: -len(ext)] + ".json")
    return img.with_suffix(".json")


def read_meta(phase_files: dict[int, Path]) -> dict:
    """Echo times (ms, ordered by echo number) plus readout time and PE direction.

    EchoTime is read per echo; TotalReadoutTime and PhaseEncodingDirection must agree across
    echoes -- a disagreement means the files are not one acquisition and is an error, not a
    thing to average.
    """
    echoes = sorted(phase_files)
    tes, trt, ped = [], None, None
    for e in echoes:
        js = sidecar(phase_files[e])
        if not js.is_file():
            raise SystemExit(f"missing sidecar: {js}")
        with open(js) as fh:
            meta = json.load(fh)
        if "EchoTime" not in meta:
            raise SystemExit(f"{js}: no EchoTime")
        tes.append(float(meta["EchoTime"]) * 1000.0)  # BIDS seconds -> ms
        t = meta.get("TotalReadoutTime")
        if t is None:
            ees, npe = meta.get("EffectiveEchoSpacing"), meta.get("ReconMatrixPE")
            if ees is not None and npe is not None:
                t = float(ees) * (int(npe) - 1)
        p = meta.get("PhaseEncodingDirection")
        if t is not None:
            if trt is not None and abs(trt - float(t)) > 1e-9:
                raise SystemExit(f"{js}: TotalReadoutTime {t} disagrees with {trt} from an earlier echo")
            trt = float(t)
        if p is not None:
            if ped is not None and p != ped:
                raise SystemExit(f"{js}: PhaseEncodingDirection {p} disagrees with {ped} from an earlier echo")
            ped = p
    if trt is None:
        raise SystemExit("no TotalReadoutTime, and EffectiveEchoSpacing + ReconMatrixPE are not both present")
    if ped is None:
        raise SystemExit("no PhaseEncodingDirection in any phase sidecar")
    return {"tes": tes, "trt": trt, "ped": ped, "echoes": echoes}


def run_cmd(cmd: list[str], dry: bool) -> None:
    printable = " ".join(shlex_quote(c) for c in cmd)
    print(printable, flush=True)
    if dry:
        return
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise SystemExit(f"command failed with exit code {r.returncode}")


def shlex_quote(s: str) -> str:
    # paths with spaces must survive being printed and pasted back into a shell
    if s and all(c.isalnum() or c in "-_./=,:+@" for c in s):
        return s
    return "'" + s.replace("'", "'\\''") + "'"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path, help="a BIDS func/ directory or a dataset root")
    ap.add_argument("--out-dir", type=Path, required=True, help="output directory")
    ap.add_argument("--niimath", default=os.environ.get("NIIMATH", "niimath"), help="niimath executable")
    ap.add_argument("--n-cpus", type=int, default=None)
    ap.add_argument("--noise-frames", type=int, default=None)
    ap.add_argument("--rank", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true", help="print the commands without running them")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--no-apply", action="store_true", help="estimate only; skip the -unwarp step")
    args = ap.parse_args(argv)

    exe = shutil.which(args.niimath) or args.niimath
    if not args.dry_run and not (Path(exe).is_file() or shutil.which(args.niimath)):
        raise SystemExit(f"niimath not found: {args.niimath} (set --niimath or $NIIMATH)")

    runs = find_runs(args.input)
    if not runs:
        raise SystemExit(f"no *_echo-<N>_part-{{mag,phase}}_bold.nii* under {args.input}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for stem in sorted(runs):
        run = runs[stem]
        mag_e, pha_e = sorted(run["mag"]), sorted(run["phase"])
        if mag_e != pha_e:
            print(f"skipping {stem}: magnitude echoes {mag_e} != phase echoes {pha_e}", file=sys.stderr)
            continue
        if len(mag_e) < 2:
            print(f"skipping {stem}: needs at least two echoes, found {len(mag_e)}", file=sys.stderr)
            continue
        meta = read_meta(run["phase"])
        prefix = args.out_dir / stem
        if not args.overwrite and Path(str(prefix) + "_displacementmaps.nii.gz").exists():
            print(f"skipping {stem}: outputs exist (use --overwrite)", file=sys.stderr)
            continue

        cmd = [exe, "--medic",
               "--magnitude", *[str(run["mag"][e]) for e in mag_e],
               "--phase", *[str(run["phase"][e]) for e in pha_e],
               "--te-ms", ",".join(f"{t:g}" for t in meta["tes"]),
               "--total-readout-time", f"{meta['trt']:g}",
               "--phase-encoding-direction", meta["ped"],
               "--out-prefix", str(prefix)]
        if args.n_cpus:
            cmd += ["--n-cpus", str(args.n_cpus)]
        if args.noise_frames is not None:
            cmd += ["--noise-frames", str(args.noise_frames)]
        if args.rank is not None:
            cmd += ["--rank", str(args.rank)]
        run_cmd(cmd, args.dry_run)

        if args.no_apply:
            continue
        # One displacement map series corrects every echo: all echoes share one EPI readout.
        dmap = str(prefix) + "_displacementmaps.nii.gz"
        axis = meta["ped"][0]  # the sign already lives in the map; -unwarp ignores a '-' suffix
        for e in mag_e:
            src = run["mag"][e]
            out = args.out_dir / (src.name.split(".nii")[0] + "_undistorted.nii.gz")
            run_cmd([exe, str(src), "-unwarp", dmap, axis, str(out)], args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())

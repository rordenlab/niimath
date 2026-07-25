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
        run = runs.setdefault(stem, {"mag": {}, "phase": {}, "dupes": []})
        if echo in run[part]:
            # Two files claiming the same echo/part (e.g. a .nii and a .nii.gz side by side).
            # Silently keeping the last one would quietly change which data was processed.
            run["dupes"].append((part, echo, run[part][echo], f))
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


OUTPUT_SUFFIXES = ("_fieldmaps_native", "_fieldmaps", "_displacementmaps")
EXTS = (".nii.gz", ".nii", ".nii.zst")
STALE = "\n  delete the stale file(s) or point --out-dir at a clean directory"


def candidates(base: Path) -> list[Path]:
    """Every existing file `base` could have been written as: niimath re-derives the extension
    from FSLOUTPUTTYPE, so it is not knowable in advance."""
    return [p for p in (Path(str(base) + e) for e in EXTS) if p.is_file()]


def resolve_outputs(prefix: Path) -> tuple[dict, str | None]:
    """Map each --medic output suffix to the one file that IS it, or return an error.

    Never picks among extensions by priority: a leftover `.nii.gz` beside a fresh `.nii` is two
    different runs, and choosing one silently hands -unwarp a stale displacement map.
    """
    got = {}
    for sfx in OUTPUT_SUFFIXES:
        c = candidates(Path(str(prefix) + sfx))
        if len(c) > 1:
            return got, "more than one file is %s%s:\n  %s%s" % (
                prefix.name, sfx, "\n  ".join(str(p) for p in c), STALE)
        if c:
            got[sfx] = c[0]
    if len(set(str(p)[len(str(prefix)) + len(s):] for s, p in got.items())) > 1:
        return got, "outputs left over from different runs (mixed formats):\n  %s%s" % (
            "\n  ".join(str(p) for p in got.values()), STALE)
    return got, None


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
    processed = skipped = failed = 0
    for stem in sorted(runs):
        run = runs[stem]
        mag_e, pha_e = sorted(run["mag"]), sorted(run["phase"])
        if mag_e != pha_e:
            print(f"error: {stem}: magnitude echoes {mag_e} != phase echoes {pha_e}", file=sys.stderr)
            failed += 1
            continue
        if len(mag_e) < 2:
            print(f"error: {stem}: needs at least two echoes, found {len(mag_e)}", file=sys.stderr)
            failed += 1
            continue
        if run["dupes"]:
            for part, echo, first, second in run["dupes"]:
                print(f"error: {stem}: two files claim echo {echo} part {part}:\n  {first}\n  {second}",
                      file=sys.stderr)
            failed += 1
            continue
        meta = read_meta(run["phase"])
        prefix = args.out_dir / stem
        outs, err = resolve_outputs(prefix)
        if err:
            print(f"error: {stem}: {err}", file=sys.stderr)
            failed += 1
            continue
        if args.noise_frames and not args.no_apply:
            # --medic drops N trailing frames, so the map is shorter than the magnitude series and
            # -unwarp correctly refuses the pairing. Say so BEFORE spending an estimate on it.
            print(f"error: {stem}: --noise-frames {args.noise_frames} makes the displacement map "
                  f"shorter than the magnitude series, which -unwarp will reject. Trim the "
                  f"magnitudes first (niimath <mag> -crop 0 N ...), pass --no-apply, or drop "
                  f"--noise-frames.", file=sys.stderr)
            failed += 1
            continue

        # Estimate and apply are independent, separately resumable stages: a complete estimate is
        # reused rather than recomputed, and reusing it does NOT skip the apply.
        worked = False
        if args.overwrite or len(outs) != len(OUTPUT_SUFFIXES):
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
            before = {p: p.stat().st_mtime_ns
                      for sfx in OUTPUT_SUFFIXES for p in candidates(Path(str(prefix) + sfx))}
            run_cmd(cmd, args.dry_run)
            worked = True
            outs = {}  # --dry-run: nothing was written, so nothing can be observed
            if not args.dry_run:
                outs, err = resolve_outputs(prefix)
                # Keep only what THIS estimate wrote: a file left over from an earlier output
                # format is not an output of this run, whatever its name says.
                if not err:
                    outs = {s: p for s, p in outs.items() if before.get(p) != p.stat().st_mtime_ns}
                    if len(outs) != len(OUTPUT_SUFFIXES):
                        err = "--medic did not write " + ", ".join(
                            s for s in OUTPUT_SUFFIXES if s not in outs)
                if err:
                    print(f"error: {stem}: {err}", file=sys.stderr)
                    failed += 1
                    continue
        if not args.no_apply:
            base = {e: args.out_dir / (run["mag"][e].name.split(".nii")[0] + "_undistorted")
                    for e in mag_e}
            todo = [e for e in mag_e if args.overwrite or not candidates(base[e])]
            if todo and not worked:
                print(f"{stem}: reusing the existing estimate (use --overwrite to redo it)",
                      file=sys.stderr)
            # One displacement map series corrects every echo: all echoes share one EPI readout.
            dmap = outs.get("_displacementmaps")
            if dmap is None:  # --dry-run only: nothing was written, so nothing can be observed
                dmap = Path(str(prefix) + "_displacementmaps.nii.gz")
                if todo:
                    print(f"# dry-run: PREDICTED map path (the real extension follows "
                          f"FSLOUTPUTTYPE): {dmap}", file=sys.stderr)
            axis = meta["ped"][0]  # the sign already lives in the map; -unwarp ignores a '-' suffix
            for e in todo:
                cmd = [exe, str(run["mag"][e])]
                if args.n_cpus:
                    cmd += ["-p", str(args.n_cpus)]   # was not propagated to the apply step
                cmd += ["-unwarp", str(dmap), axis, str(base[e]) + ".nii.gz"]
                run_cmd(cmd, args.dry_run)
                worked = True
        if worked:
            processed += 1
        else:
            print(f"skipping {stem}: all outputs exist (use --overwrite)", file=sys.stderr)
            skipped += 1
    if failed:
        print(f"{processed} run(s) processed, {skipped} skipped, {failed} FAILED", file=sys.stderr)
        return 1
    if processed == 0 and skipped == 0:
        print("no runs processed", file=sys.stderr)
        return 1
    print(f"{processed} run(s) processed, {skipped} skipped", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

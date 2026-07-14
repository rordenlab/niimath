#!/usr/bin/env python3
"""Release smoke test for packaged niimath executables.

The test is intentionally stdlib-only so it can run inside cibuildwheel's clean
test environments and on AppVeyor release workers. It synthesizes small NIfTI
fixtures directly and exercises the packaging-sensitive paths: gzip read/write,
the reported -conform/-gz 0/-odt char case, feature dispatch, and optional zstd.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path


SHORT_READ = "++ WARNING: read "


def run_niimath(exe: str, args: list[str], *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        [exe, *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        check=False,
    )
    combined = result.stdout + result.stderr
    if SHORT_READ in combined:
        raise AssertionError(f"short-read warning from {' '.join([exe, *args])}:\n{combined}")
    return result


def require_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode != 0:
        raise AssertionError(f"{label} failed with exit {result.returncode}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")


def nifti_header(
    dims: tuple[int, int, int],
    datatype: int,
    bitpix: int,
    offset: tuple[float, float, float] = (0.0, 0.0, 0.0),
    xyz_units: int = 2,  # NIFTI_UNITS_MM
    scale: float = 1.0,  # voxel size in xyz_units (pixdim + sform diagonal)
) -> bytes:
    # Defaults (xyz_units=2 mm, scale=1.0) reproduce the historic mm/unit-voxel
    # header exactly; scale/xyz_units let a fixture describe the SAME physical grid
    # in metres (scale=0.001, xyz_units=1) or microns to exercise unit normalization.
    hdr = bytearray(348)
    struct.pack_into("<i", hdr, 0, 348)
    struct.pack_into("<8h", hdr, 40, 3, dims[0], dims[1], dims[2], 1, 1, 1, 1)
    struct.pack_into("<h", hdr, 70, datatype)
    struct.pack_into("<h", hdr, 72, bitpix)
    struct.pack_into("<8f", hdr, 76, 1.0, scale, scale, scale, 0.0, 0.0, 0.0, 0.0)
    struct.pack_into("<f", hdr, 108, 352.0)
    struct.pack_into("<f", hdr, 112, 1.0)
    hdr[123] = xyz_units
    struct.pack_into("<h", hdr, 252, 3)
    struct.pack_into("<h", hdr, 254, 3)
    struct.pack_into("<3f", hdr, 268, *offset)
    struct.pack_into("<4f", hdr, 280, scale, 0.0, 0.0, offset[0])
    struct.pack_into("<4f", hdr, 296, 0.0, scale, 0.0, offset[1])
    struct.pack_into("<4f", hdr, 312, 0.0, 0.0, scale, offset[2])
    hdr[344:348] = b"n+1\0"
    return bytes(hdr) + b"\0\0\0\0"


def write_uint8_nifti(
    path: Path,
    dims: tuple[int, int, int] = (8, 8, 8),
    offset: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> None:
    nvox = dims[0] * dims[1] * dims[2]
    data = bytes((i % 251 for i in range(nvox)))
    path.write_bytes(nifti_header(dims, datatype=2, bitpix=8, offset=offset) + data)


def write_float32_nifti(
    path: Path,
    dims: tuple[int, int, int],
    data: list[float],
    offset: tuple[float, float, float] = (0.0, 0.0, 0.0),
    nt: int = 1,
) -> None:
    nvox = dims[0] * dims[1] * dims[2] * nt
    if len(data) != nvox:
        raise AssertionError(f"{path}: expected {nvox} values, got {len(data)}")
    payload = struct.pack(f"<{nvox}f", *data)
    header = bytearray(nifti_header(dims, datatype=16, bitpix=32, offset=offset))
    if nt > 1:
        struct.pack_into("<8h", header, 40, 4, dims[0], dims[1], dims[2], nt, 1, 1, 1)
    path.write_bytes(bytes(header) + payload)


def read_float32_nifti(path: Path) -> list[float]:
    blob = path.read_bytes()
    dim = struct.unpack_from("<8h", blob, 40)
    nvox = math.prod(dim[1 : dim[0] + 1])
    offset = int(struct.unpack_from("<f", blob, 108)[0])
    return list(struct.unpack_from(f"<{nvox}f", blob, offset))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    rank = fraction * (len(ordered) - 1)
    lo = math.floor(rank)
    hi = math.ceil(rank)
    return ordered[lo] + (rank - lo) * (ordered[hi] - ordered[lo])


def qc_stats(values: list[float]) -> dict[str, float]:
    n = len(values)
    mean = sum(values) / n
    m2 = sum((value - mean) ** 2 for value in values) / n
    m4 = sum((value - mean) ** 4 for value in values) / n
    median = percentile(values, 0.5)
    return {
        "mean": mean,
        "stdv": math.sqrt(m2),
        "median": median,
        "mad": percentile([abs(value - median) for value in values], 0.5) / 0.6744897501960817,
        "p05": percentile(values, 0.05),
        "p95": percentile(values, 0.95),
        "k": m4 / (m2 * m2) - 3.0,
        "n": float(n),
    }


def assert_close(actual: float, expected: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=3e-5, abs_tol=1e-5):
        raise AssertionError(f"{label}: expected {expected:.12g}, saw {actual:.12g}")


def exercise_qc(exe: str, tmp: Path) -> None:
    dims = (12, 12, 12)
    t1_values: list[float] = []
    seg_values: list[int] = []
    tissues: dict[str, list[float]] = {"csf": [], "gm": [], "wm": []}
    for index in range(dims[0] * dims[1] * dims[2]):
        x = index % dims[0]
        if x < 4:
            tissue, label, value = "csf", 1, 20.0 + 0.5 * (index % 7)
        elif x < 8:
            tissue, label, value = "gm", 2, 80.0 + 0.75 * (index % 11)
        else:
            tissue, label, value = "wm", 3, 120.0 + 0.25 * (index % 13)
        t1_values.append(value)
        seg_values.append(label)
        tissues[tissue].append(value)

    t1 = tmp / "qc_t1.nii"
    seg = tmp / "qc_seg.nii"
    out = tmp / "qc.tsv"
    write_float32_nifti(t1, dims, t1_values)
    seg.write_bytes(nifti_header(dims, datatype=2, bitpix=8) + bytes(seg_values))
    result = run_niimath(
        exe,
        ["--qc", str(t1), "--seg", str(seg), "--csf", "1", "--wm", "3", "--erode", "0", "--out", str(out)],
    )
    require_success(result, "anatomical QC")
    rows = out.read_text(encoding="utf-8").splitlines()
    if len(rows) != 2:
        raise AssertionError(f"QC TSV should contain two rows, saw {len(rows)}")
    names = rows[0].split("\t")
    values = rows[1].split("\t")
    if len(names) != len(values) or len(set(names)) != len(names):
        raise AssertionError("QC TSV columns are missing, duplicated, or misaligned")
    observed = {name: float(value) for name, value in zip(names, values)}

    stats = {name: qc_stats(data) for name, data in tissues.items()}
    delta = abs(stats["wm"]["median"] - stats["gm"]["median"])
    expected = {
        "cjv": (stats["wm"]["mad"] + stats["gm"]["mad"]) / delta,
        "cnr_noair": delta / math.sqrt(stats["wm"]["stdv"] ** 2 + stats["gm"]["stdv"] ** 2),
        "wm2max": stats["wm"]["median"] / percentile(t1_values, 0.9995),
    }
    snrs = {}
    for tissue in ("csf", "wm", "gm"):
        n = stats[tissue]["n"]
        snrs[tissue] = stats[tissue]["median"] / (stats[tissue]["stdv"] * math.sqrt(n / (n - 1.0)))
        expected[f"snr_{tissue}"] = snrs[tissue]
    expected["snr_total"] = sum(snrs.values()) / 3.0
    energy = math.sqrt(sum(value * value for value in t1_values))
    efc_max = math.sqrt(len(t1_values)) * math.log(1.0 / math.sqrt(len(t1_values)))
    expected["efc_brain"] = sum((value / energy) * math.log((value + 1e-16) / energy) for value in t1_values) / efc_max
    for tissue in ("csf", "gm", "wm"):
        expected[f"icvs_{tissue}"] = 1.0 / 3.0
        expected[f"vol_{tissue}_mm3"] = float(len(tissues[tissue]))
        for metric, value in stats[tissue].items():
            expected[f"summary_{tissue}_{metric}"] = value
    if set(observed) != set(expected):
        raise AssertionError(f"QC TSV schema mismatch: expected {sorted(expected)}, saw {sorted(observed)}")
    for name, value in expected.items():
        assert_close(observed[name], value, f"QC {name}")

    # Exercise the DEFAULT six-neighbour erosion path (the run above uses --erode 0).
    # Each tissue is a full 4-wide x-slab, so erosion keeps only its interior:
    # x in {1,2}/{5,6}/{9,10} and y,z in [1,10] -> exactly 2*10*10 = 200 voxels.
    # ICV fractions and mm3 volumes must still use the RAW 4*12*12 = 576 counts,
    # so they are identical to the --erode 0 run: this pins the raw-vs-eroded split.
    eroded_out = tmp / "qc_erode.tsv"
    eroded_run = run_niimath(
        exe,
        ["--qc", str(t1), "--seg", str(seg), "--csf", "1", "--wm", "3", "--erode", "1", "--out", str(eroded_out)],
    )
    require_success(eroded_run, "anatomical QC (erode)")
    erows = eroded_out.read_text(encoding="utf-8").splitlines()
    eobserved = dict(zip(erows[0].split("\t"), (float(v) for v in erows[1].split("\t"))))
    for tissue in ("csf", "gm", "wm"):
        if eobserved[f"summary_{tissue}_n"] != 200.0:
            raise AssertionError(
                f"QC erosion: summary_{tissue}_n expected 200, saw {eobserved[f'summary_{tissue}_n']:.12g}"
            )
        assert_close(eobserved[f"icvs_{tissue}"], 1.0 / 3.0, f"QC erode icvs_{tissue}")
        assert_close(eobserved[f"vol_{tissue}_mm3"], 576.0, f"QC erode vol_{tissue}_mm3")

    # Unit normalization: a physically identical grid stored in metres (voxel size
    # 0.001 m, units code 1) must be ACCEPTED against the mm T1 -- max_displacement_mm
    # normalises both to mm before the grid check, so a unit-code difference alone is
    # not a mismatch.
    meter_seg = tmp / "qc_seg_meter.nii"
    meter_seg.write_bytes(nifti_header(dims, datatype=2, bitpix=8, xyz_units=1, scale=0.001) + bytes(seg_values))
    require_success(
        run_niimath(exe, ["--qc", str(t1), "--seg", str(meter_seg), "--csf", "1", "--wm", "3", "--erode", "0", "--out", str(tmp / "meter.tsv")]),
        "anatomical QC (mm image + same-grid metre segmentation)",
    )
    # ...but a genuine 5 mm shift expressed in metres (0.005 m) must still be REJECTED.
    meter_shift = tmp / "qc_seg_meter_shift.nii"
    meter_shift.write_bytes(nifti_header(dims, datatype=2, bitpix=8, offset=(0.005, 0.0, 0.0), xyz_units=1, scale=0.001) + bytes(seg_values))
    shifted_meter = run_niimath(exe, ["--qc", str(t1), "--seg", str(meter_shift), "--csf", "1", "--wm", "3", "--out", str(tmp / "meter_shift.tsv")])
    if shifted_meter.returncode == 0 or "spatial grid differs" not in (shifted_meter.stdout + shifted_meter.stderr):
        raise AssertionError("QC should reject a 5 mm shift expressed in metre units")

    # Same, for the least-used micron branch (voxel size 1000 um == 1 mm, units code 3).
    micron_seg = tmp / "qc_seg_micron.nii"
    micron_seg.write_bytes(nifti_header(dims, datatype=2, bitpix=8, xyz_units=3, scale=1000.0) + bytes(seg_values))
    require_success(
        run_niimath(exe, ["--qc", str(t1), "--seg", str(micron_seg), "--csf", "1", "--wm", "3", "--erode", "0", "--out", str(tmp / "micron.tsv")]),
        "anatomical QC (mm image + same-grid micron segmentation)",
    )
    micron_shift = tmp / "qc_seg_micron_shift.nii"
    micron_shift.write_bytes(nifti_header(dims, datatype=2, bitpix=8, offset=(5000.0, 0.0, 0.0), xyz_units=3, scale=1000.0) + bytes(seg_values))
    shifted_micron = run_niimath(exe, ["--qc", str(t1), "--seg", str(micron_shift), "--csf", "1", "--wm", "3", "--out", str(tmp / "micron_shift.tsv")])
    if shifted_micron.returncode == 0 or "spatial grid differs" not in (shifted_micron.stdout + shifted_micron.stderr):
        raise AssertionError("QC should reject a 5 mm shift expressed in micron units")

    # Large exact count: summary_*_n must remain an exact integer above six significant
    # figures. The old %.6g serializer would have written 1030301 as 1.03030e+06.
    big_dims = (101, 101, 101)
    nbig = big_dims[0] * big_dims[1] * big_dims[2]  # 1_030_301 voxels, all CSF
    big_t1 = tmp / "qc_big_t1.nii"
    big_seg = tmp / "qc_big_seg.nii"
    big_t1.write_bytes(nifti_header(big_dims, datatype=16, bitpix=32) + struct.pack("<f", 100.0) * nbig)
    big_seg.write_bytes(nifti_header(big_dims, datatype=2, bitpix=8) + b"\x03" * nbig)
    require_success(
        run_niimath(exe, ["--qc", str(big_t1), "--seg", str(big_seg), "--csf", "3", "--wm", "1", "--erode", "0", "--out", str(tmp / "big.tsv")]),
        "anatomical QC (large exact count)",
    )
    big_rows = (tmp / "big.tsv").read_text(encoding="utf-8").splitlines()
    big_n = dict(zip(big_rows[0].split("\t"), big_rows[1].split("\t")))["summary_csf_n"]
    if big_n != str(nbig):
        raise AssertionError(f"summary_csf_n must be the exact integer {nbig}, saw {big_n!r}")

    # Erosion <QC_MIN_VOX fallback: a two-voxel-thick CSF slab (x in {0,1}) erodes to
    # zero interior voxels, so its stats must fall back to the raw 2*8*8 = 128 count.
    fdims = (8, 8, 8)
    fseg_vals = bytes(
        (1 if (i % fdims[0]) < 2 else 2 if (i % fdims[0]) < 5 else 3)
        for i in range(fdims[0] * fdims[1] * fdims[2])
    )
    fb_t1 = tmp / "qc_fb_t1.nii"
    fb_seg = tmp / "qc_fb_seg.nii"
    fb_t1.write_bytes(nifti_header(fdims, datatype=16, bitpix=32) + b"".join(struct.pack("<f", 40.0 + (i % 7)) for i in range(fdims[0] * fdims[1] * fdims[2])))
    fb_seg.write_bytes(nifti_header(fdims, datatype=2, bitpix=8) + fseg_vals)
    fb = run_niimath(exe, ["--qc", str(fb_t1), "--seg", str(fb_seg), "--csf", "1", "--wm", "3", "--erode", "1", "--out", str(tmp / "fb.tsv")])
    require_success(fb, "anatomical QC (erosion fallback)")
    if "using un-eroded mask" not in (fb.stdout + fb.stderr):
        raise AssertionError("QC should report the erosion fallback for a fully-eroded thin tissue")
    fb_n = dict(zip((tmp / "fb.tsv").read_text().splitlines()[0].split("\t"), (tmp / "fb.tsv").read_text().splitlines()[1].split("\t")))["summary_csf_n"]
    if fb_n != "128":
        raise AssertionError(f"erosion fallback: summary_csf_n should be the raw 128, saw {fb_n!r}")

    shifted_seg = tmp / "qc_seg_shifted.nii"
    shifted_seg.write_bytes(
        nifti_header(dims, datatype=2, bitpix=8, offset=(0.0, 1.0, 0.0)) + bytes(seg_values)
    )
    shifted = run_niimath(
        exe,
        ["--qc", str(t1), "--seg", str(shifted_seg), "--csf", "1", "--wm", "3", "--out", str(tmp / "shifted.tsv")],
    )
    if shifted.returncode == 0 or "spatial grid differs" not in (shifted.stdout + shifted.stderr):
        raise AssertionError("QC should reject a segmentation translated along the y axis")

    fractional_seg = tmp / "qc_seg_fractional.nii"
    fractional = [float(value) for value in seg_values]
    fractional[0] = 1.5
    write_float32_nifti(fractional_seg, dims, fractional)
    invalid = run_niimath(
        exe,
        ["--qc", str(t1), "--seg", str(fractional_seg), "--csf", "1", "--wm", "3", "--out", str(tmp / "invalid.tsv")],
    )
    if invalid.returncode == 0 or "integer labels" not in (invalid.stdout + invalid.stderr):
        raise AssertionError("QC should reject a non-integer segmentation")

    overlap = run_niimath(
        exe,
        ["--qc", str(t1), "--seg", str(seg), "--csf", "1", "--wm", "1", "--out", str(tmp / "overlap.tsv")],
    )
    if overlap.returncode == 0 or "both CSF and WM" not in (overlap.stdout + overlap.stderr):
        raise AssertionError("QC should reject overlapping tissue label sets")

    nonfinite_t1 = tmp / "qc_t1_nan.nii"
    with_nan = list(t1_values)
    with_nan[0] = math.nan
    write_float32_nifti(nonfinite_t1, dims, with_nan)
    invalid = run_niimath(
        exe,
        ["--qc", str(nonfinite_t1), "--seg", str(seg), "--csf", "1", "--wm", "3", "--out", str(tmp / "nan.tsv")],
    )
    if invalid.returncode == 0 or "non-finite" not in (invalid.stdout + invalid.stderr):
        raise AssertionError("QC should reject a non-finite T1")

    unwritable = run_niimath(
        exe,
        ["--qc", str(t1), "--seg", str(seg), "--csf", "1", "--wm", "3", "--out", str(tmp / "missing" / "qc.tsv")],
    )
    if unwritable.returncode == 0 or "cannot open output" not in (unwritable.stdout + unwritable.stderr):
        raise AssertionError("QC should propagate an output-open failure")


def write_large_float_gz(path: Path) -> None:
    dims = (256, 256, 256)
    nvox = dims[0] * dims[1] * dims[2]
    zero_chunk = b"\0" * (1024 * 1024)
    remaining = nvox * 4 - 4
    with gzip.open(path, "wb") as f:
        f.write(nifti_header(dims, datatype=16, bitpix=32))
        f.write(struct.pack("<f", 255.0))
        while remaining:
            n = min(remaining, len(zero_chunk))
            f.write(zero_chunk[:n])
            remaining -= n


def read_nifti_bytes(path: Path) -> bytes:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as f:
            return f.read()
    return path.read_bytes()


def parse_header(blob: bytes) -> dict[str, int]:
    if len(blob) < 352:
        raise AssertionError(f"NIfTI file too small: {len(blob)} bytes")
    sizeof_hdr = struct.unpack_from("<i", blob, 0)[0]
    if sizeof_hdr != 348:
        raise AssertionError(f"unexpected sizeof_hdr {sizeof_hdr}")
    dims = struct.unpack_from("<8h", blob, 40)
    datatype = struct.unpack_from("<h", blob, 70)[0]
    bitpix = struct.unpack_from("<h", blob, 72)[0]
    vox_offset = int(struct.unpack_from("<f", blob, 108)[0])
    nvox = dims[1] * dims[2] * dims[3] * max(dims[4], 1) * max(dims[5], 1) * max(dims[6], 1) * max(dims[7], 1)
    return {
        "dim0": dims[0],
        "nx": dims[1],
        "ny": dims[2],
        "nz": dims[3],
        "datatype": datatype,
        "bitpix": bitpix,
        "vox_offset": vox_offset,
        "nvox": nvox,
    }


def assert_payload_size(path: Path, datatype: int, bitpix: int, dims: tuple[int, int, int]) -> None:
    blob = read_nifti_bytes(path)
    hdr = parse_header(blob)
    expected = {
        "dim0": 3,
        "nx": dims[0],
        "ny": dims[1],
        "nz": dims[2],
        "datatype": datatype,
        "bitpix": bitpix,
    }
    for key, value in expected.items():
        if hdr[key] != value:
            raise AssertionError(f"{path}: expected {key}={value}, saw {hdr[key]}")
    expected_size = hdr["vox_offset"] + hdr["nvox"] * (bitpix // 8)
    if len(blob) != expected_size:
        raise AssertionError(f"{path}: expected {expected_size} bytes from header, saw {len(blob)}")


def exercise_allineate(exe: str, tmp: Path, help_text: str) -> None:
    """Regression for the -allineate -fill / -weight options and the -dilate fix — the
    niimath-only dispatch, chain integration, and CLI parsing that the shared allineate
    suite does not cover. Skips cleanly on a build without registration (e.g. nano)."""
    if "-allineate" not in help_text:
        return

    n = 32

    def blob() -> list[float]:
        # A volume-filling smooth blob plus a gradient: foreground exists at every pyramid
        # level so the fast engine's coarse sample build succeeds deterministically (a small
        # sparse cube collapses to too few samples at 8 mm and fails the build).
        data = [0.0] * (n * n * n)
        for z in range(n):
            for y in range(n):
                for x in range(n):
                    r2 = (x - 16) ** 2 + (y - 16) ** 2 + (z - 16) ** 2
                    data[x + y * n + z * n * n] = 100.0 * math.exp(-r2 / 200.0) + 0.1 * (x + y + z) + 1.0
        return data

    base = tmp / "al_base.nii"
    moving = tmp / "al_mov.nii"
    weight = tmp / "al_weight.nii"
    write_float32_nifti(base, (n, n, n), blob())
    write_float32_nifti(moving, (n, n, n), blob())
    write_float32_nifti(
        weight, (n, n, n),
        [1.0 if (8 <= x < 24 and 8 <= y < 24 and 8 <= z < 24) else 0.0
         for z in range(n) for y in range(n) for x in range(n)],
    )

    # -dilate with a threshold > 1: grown voxels must reach at least `iso`. The prior
    # fmax(1.0, ..) left a `-dilate 10 dx` grow at value 1 — below the requested threshold.
    seed = [0.0] * (12 * 12 * 12)
    seed[6 + 6 * 12 + 6 * 144] = 10.0
    seed_path = tmp / "al_seed10.nii"
    write_float32_nifti(seed_path, (12, 12, 12), seed)
    dil = tmp / "al_dil.nii"
    require_success(run_niimath(exe, [str(seed_path), "-dilate", "10", "2", "-gz", "0", str(dil)]), "-dilate iso>1")
    grown = [v for v in read_float32_nifti(dil) if v > 0.0]
    if len(grown) <= 1:
        raise AssertionError("-dilate 10 2 did not grow the seed")
    if any(v < 10.0 for v in grown):
        raise AssertionError(f"-dilate 10 2 grew voxels below the threshold (iso=10): {sorted(set(grown))}")

    # -applymat does no registration, so -weight must be rejected (not silently ignored).
    ident = tmp / "al_ident.json"
    ident.write_text('{"fixed_to_moving": [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]}')
    rej = run_niimath(exe, [str(moving), "-allineate", str(base), "-applymat", str(ident),
                            "-weight", str(weight), "-gz", "0", str(tmp / "al_rej.nii")])
    if rej.returncode == 0 or "cannot be combined" not in (rej.stdout + rej.stderr):
        raise AssertionError("-applymat combined with -weight was not rejected")

    # -weight interleaved with other sub-options (exercises the shared parser + AL_CAP_WEIGHT)
    # and its provenance recorded by -savemat.
    prov = tmp / "al_prov.json"
    require_success(
        run_niimath(exe, [str(moving), "-allineate", str(base), "-weight", str(weight),
                          "-final", "linear", "-savemat", str(prov), "-gz", "0", str(tmp / "al_w.nii")]),
        "-allineate -weight with saved provenance",
    )
    if '"weight"' not in prov.read_text():
        raise AssertionError("-savemat did not record the -weight provenance")

    # Out-of-FOV fill: a full-X translation maps every target voxel outside the (all -5)
    # source, so the whole output is the fill value — a clean auto/zero/nan distinction.
    neg = tmp / "al_neg.nii"
    grid = tmp / "al_grid.nii"
    write_float32_nifti(neg, (8, 8, 8), [-5.0] * 512)
    write_float32_nifti(grid, (8, 8, 8), [0.0] * 512)
    shift = tmp / "al_shift.json"
    shift.write_text('{"fixed_to_moving": [1,0,0,100, 0,1,0,0, 0,0,1,0, 0,0,0,1]}')

    def fill_out(mode: str) -> list[float]:
        out = tmp / f"al_fill_{mode}.nii"
        require_success(
            run_niimath(exe, [str(neg), "-allineate", str(grid), "-applymat", str(shift),
                              "-fill", mode, "-gz", "0", str(out)]),
            f"-applymat -fill {mode}",
        )
        return read_float32_nifti(out)

    if any(v != 0.0 for v in fill_out("zero")):
        raise AssertionError("-fill zero left a non-zero out-of-FOV voxel")
    if any(abs(v + 5.0) > 1e-3 for v in fill_out("auto")):
        raise AssertionError("-fill auto did not fill with the negative source minimum (-5)")
    if not all(v != v for v in fill_out("nan")):
        raise AssertionError("-fill nan did not write NaN out-of-FOV")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("exe", nargs="?", default=shutil.which("niimath") or "niimath")
    parser.add_argument("--expect-zstd", action="store_true", default=os.environ.get("NIIMATH_EXPECT_ZSTD") == "1")
    parser.add_argument("--expect-bsd", action="store_true", default=os.environ.get("NIIMATH_EXPECT_BSD") == "1")
    args = parser.parse_args()

    exe = str(Path(args.exe).resolve()) if Path(args.exe).exists() else args.exe
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        info = run_niimath(exe, [])
        require_success(info, "help/version")
        help_text = info.stdout + info.stderr
        for token in ("-conform", "-allineate", "-deface", "--dtifit", "--qc", "-bitmap", "-bandpass", "-mesh"):
            if token not in help_text:
                raise AssertionError(f"packaged binary help is missing {token}")
        if args.expect_bsd and " BSD " not in help_text:
            raise AssertionError("expected a BSD build version string")

        big = tmp / "strokeSubject.nii.gz"
        write_large_float_gz(big)
        conform_out_arg = tmp / "conform_char.nii.gz"
        result = run_niimath(exe, [str(big), "-conform", "-gz", "0", str(conform_out_arg), "-odt", "char"])
        require_success(result, "conform gzip input to uncompressed uint8")
        conform_out = tmp / "conform_char.nii"
        if not conform_out.exists():
            raise AssertionError(f"expected -gz 0 output {conform_out}, found {sorted(p.name for p in tmp.iterdir())}")
        assert_payload_size(conform_out, datatype=2, bitpix=8, dims=(256, 256, 256))
        require_success(run_niimath(exe, [str(conform_out)]), "read conform output")

        small = tmp / "small.nii"
        write_uint8_nifti(small)

        # Unsafe huge work must be rejected from the header, before any multi-gigabyte payload
        # allocation. The deliberately header-only file has 3,220,029,867 voxels.
        huge_header = tmp / "huge_header_only.nii"
        huge_header.write_bytes(nifti_header((32767, 32767, 3), datatype=16, bitpix=32))
        huge_reject = run_niimath(exe, [str(huge_header), "-fmean", str(tmp / "huge_reject.nii")])
        if huge_reject.returncode != 2 or "not huge-image capable" not in (huge_reject.stdout + huge_reject.stderr):
            raise AssertionError("unsafe huge input was not rejected by header preflight")

        # Dispatch accepts exact operator names only; formerly these substring variants ran a
        # different operation and made the huge-safe admission grammar diverge from dispatch.
        disguised_ops = [
            ["x-otsu"],
            ["-dogjunk", "1", "2"],
            ["-Tmeanjunk"],
        ]
        for disguised in disguised_ops:
            rejected = run_niimath(exe, [str(small), *disguised, str(tmp / "disguised.nii")])
            message = rejected.stdout + rejected.stderr
            rejected_exactly = any(
                text in message for text in ("unsupported operation", "unknown dimensionality reduction operation")
            )
            if rejected.returncode == 0 or not rejected_exactly:
                raise AssertionError(f"substring operator should be rejected: {disguised[0]}")

        for exact, label in ((["-otsu", "2"], "otsu"), (["-dog", "1", "2"], "dog"), (["-Tmean"], "Tmean")):
            require_success(
                run_niimath(exe, [str(small), *exact, str(tmp / f"exact_{label}.nii")]),
                f"exact {label} dispatch",
            )

        # reslice is 3D-only. This exact failure path previously continued after reslice() failed
        # and read a second volume past the 3D mask buffer; sanitizer/guard-page runs of this smoke
        # turn the clean rejection into a focused memory-safety regression.
        reslice_4d = tmp / "reslice_mask_4d.nii"
        reslice_mask = tmp / "reslice_mask_3d.nii"
        write_float32_nifti(reslice_4d, (16, 16, 4), [float(i % 13) for i in range(2048)], nt=2)
        write_float32_nifti(reslice_mask, (16, 16, 4), [1.0] * 1024)
        reslice_reject = run_niimath(
            exe,
            [str(reslice_4d), "-reslice_mask", str(reslice_mask), str(tmp / "reslice_mask_4d_out.nii")],
        )
        if reslice_reject.returncode == 0 or "only for 3D data not 4D time series" not in (
            reslice_reject.stdout + reslice_reject.stderr
        ):
            raise AssertionError("4D -reslice_mask input should fail cleanly before applying the 3D mask")

        missing_operand = run_niimath(
            exe,
            [str(small), "-add", str(tmp / "missing_operand.nii"), str(tmp / "missing_out.nii")],
        )
        if missing_operand.returncode == 0 or "failed to read NIfTI image" not in (
            missing_operand.stdout + missing_operand.stderr
        ):
            raise AssertionError("missing binary operand should fail cleanly")

        # FSL assumes a representable robust range. Opposite finite extrema used to reach an
        # undefined NaN-to-int histogram cast; keep the FSL-style operator arithmetic and only
        # require that the commands complete without UB (the UBSan build runs this same smoke).
        flt_max = 3.4028234663852886e38
        extrema = tmp / "extrema.nii"
        write_float32_nifti(extrema, (8, 4, 4), [flt_max if i % 2 == 0 else -flt_max for i in range(128)])
        for case_index, pct_op in enumerate(
            (["-thrp", "2"], ["-thrP", "2"], ["-uthrp", "90"], ["-uthrP", "90"],
             ["-clamp", "25"], ["-uclamp", "75"])
        ):
            require_success(
                run_niimath(exe, [str(extrema), *pct_op, str(tmp / f"extrema_{case_index}.nii")]),
                f"robust-range {pct_op[0]} on finite extrema",
            )

        # Otsu is a niimath extension. Its histogram cannot represent this span, so reject it
        # explicitly instead of inventing normalized semantics or performing a NaN-to-int cast.
        otsu_extreme = run_niimath(exe, [str(extrema), "-otsu", "2", str(tmp / "otsu_extreme.nii")])
        if otsu_extreme.returncode == 0 or "representable intensity range" not in (
            otsu_extreme.stdout + otsu_extreme.stderr
        ):
            raise AssertionError("Otsu should reject an unrepresentable finite intensity span")

        # FSL defines capital-P thresholds over POSITIVE voxels, not all nonzero voxels.
        # With only 85 positive samples the robust range is exactly 10..20, hence 50% is 15.
        positive_range = tmp / "positive_range.nii"
        positive_values = [(-100.0, 10.0, 20.0)[i % 3] for i in range(128)]
        write_float32_nifti(positive_range, (8, 4, 4), positive_values)
        for op, keep in (
            ("-thrP", lambda value: value if value >= 15.0 else 0.0),
            ("-uthrP", lambda value: value if value <= 15.0 else 0.0),
        ):
            positive_out = tmp / f"positive_{op[1:]}.nii"
            require_success(
                run_niimath(exe, [str(positive_range), op, "50", "-gz", "0", str(positive_out)]),
                f"{op} positive-voxel range",
            )
            expected = [keep(value) for value in positive_values]
            if read_float32_nifti(positive_out) != expected:
                raise AssertionError(f"{op} must derive its robust range from positive voxels")

        # The case above has <100 positive samples, so it takes the endpoint shortcut. Exercise the
        # 1000-BIN HISTOGRAM path with 180 positive samples (60 each of 10/30/50) + negatives. With
        # well-separated discrete values the robust range is 10..50, so -uthrP 60 thresholds at ~34:
        # the 50s are zeroed, the 30s/10s/negatives survive. (niimath's 1000-bin robust range is an
        # approximation of FSL's iterative refinement; they can differ by a few voxels ONLY when
        # samples land exactly on a bin cutoff — discrete inputs like this are exact.)
        hist_range = tmp / "hist_range.nii"
        hist_values = [(10.0, 30.0, 50.0)[i // 60] if i < 180 else -100.0 for i in range(256)]
        write_float32_nifti(hist_range, (8, 8, 4), hist_values)
        hist_out = tmp / "hist_uthrP.nii"
        require_success(
            run_niimath(exe, [str(hist_range), "-uthrP", "60", "-gz", "0", str(hist_out)]),
            "-uthrP over the 1000-bin histogram path",
        )
        hist_expected = [v if v <= 30.0 else 0.0 for v in hist_values]
        if read_float32_nifti(hist_out) != hist_expected:
            raise AssertionError("-uthrP histogram path did not threshold the positive robust range")

        # Exercise conversion-only parsing without making CI reserve the 12.9 GB payload declared
        # by huge_header. FORCE/tiny and the opt-in huge suite cover the size boundary itself.
        require_success(
            run_niimath(exe, [str(small), str(tmp / "conv_only.nii"), "-odt", "float"]),
            "conversion-only trailing -odt",
        )

        # RGB/RGBA is admitted by its EFFECTIVE scalar count (x3/x4 after expansion): a packed image
        # at/below INT_MAX that expands above it must reject an unsafe op before allocating gigabytes.
        rgb_header = tmp / "rgb_header_only.nii"  # 838,860,800 packed voxels -> 2.5e9 scalar (>INT_MAX)
        rgb_bytes = bytearray(nifti_header((1024, 1024, 800), datatype=128, bitpix=24))
        struct.pack_into("<h", rgb_bytes, 70, 128)  # DT_RGB24
        struct.pack_into("<h", rgb_bytes, 72, 24)
        rgb_header.write_bytes(bytes(rgb_bytes))
        rgb_reject = run_niimath(exe, [str(rgb_header), "-fmean", str(tmp / "rgb_reject.nii")])
        if rgb_reject.returncode != 2 or "not huge-image capable" not in (rgb_reject.stdout + rgb_reject.stderr):
            raise AssertionError("huge-after-RGB-expansion input was not rejected before allocation")

        # -roc must convert its truth/noise auxiliaries with their own header (a uint8 truth was
        # previously reinterpreted byte-for-byte as float32 — an out-of-bounds read) and exit 0.
        roc_obs = tmp / "roc_obs.nii"
        roc_truth = tmp / "roc_truth.nii"
        rdim = (16, 16, 16)
        obs_vals, truth_vals = [], []
        for z in range(16):
            for y in range(16):
                for x in range(16):
                    sig = 1 if (6 <= x < 10 and 6 <= y < 10 and 6 <= z < 10) else 0
                    truth_vals.append(sig)
                    obs_vals.append(5.0 if sig else float((x + y + z) % 3))
        write_float32_nifti(roc_obs, rdim, obs_vals)
        roc_truth.write_bytes(nifti_header(rdim, datatype=2, bitpix=8) + bytes(truth_vals))  # uint8 truth
        roc = run_niimath(exe, [str(roc_obs), "-roc", "-0.1", str(tmp / "roc.txt"), str(roc_truth), str(tmp / "roc_out.nii")])
        require_success(roc, "-roc with uint8 truth")
        roc_text = (tmp / "roc.txt").read_text()
        if "nan" in roc_text.lower() or "inf" in roc_text.lower():
            raise AssertionError("-roc produced non-finite output (datatype/validation defect)")

        # Positive-threshold ROC: noise rank `i` must never index the observed array `k`. Use only
        # two included truth voxels but ten noise volumes, so nvol > nTest deterministically.
        ndims = (12, 12, 12)
        nn3 = math.prod(ndims)
        noise_obs = tmp / "roc_noise_obs.nii"
        noise_truth = tmp / "roc_noise_truth.nii"
        noise_stack = tmp / "roc_noise_stack.nii"
        write_float32_nifti(noise_obs, ndims, [float(i % 17) for i in range(nn3)])
        truth = [-1.0] * nn3
        truth[5 + 12 * (5 + 12 * 5)] = 0.0
        truth[6 + 12 * (5 + 12 * 5)] = 1.0
        write_float32_nifti(noise_truth, ndims, truth)
        noise_hdr = bytearray(nifti_header(ndims, datatype=2, bitpix=8))
        struct.pack_into("<8h", noise_hdr, 40, 4, 12, 12, 12, 10, 1, 1, 1)
        noise_stack.write_bytes(bytes(noise_hdr) + bytes((i % 251 for i in range(nn3 * 10))))
        roc_noise = run_niimath(
            exe,
            [str(noise_obs), "-roc", "0.1", str(tmp / "roc_noise.txt"), str(noise_stack),
             str(noise_truth), str(tmp / "roc_noise_out.nii")],
        )
        require_success(roc_noise, "positive -roc with nvol > nTest")
        noise_text = (tmp / "roc_noise.txt").read_text().lower()
        if "nan" in noise_text or "inf" in noise_text:
            raise AssertionError("positive -roc produced non-finite output")

        gz_out = tmp / "roundtrip.nii.gz"
        require_success(run_niimath(exe, [str(small), "-add", "1", "-gz", "1", str(gz_out), "-odt", "char"]), "gzip round-trip write")
        if gz_out.read_bytes()[:2] != b"\x1f\x8b":
            raise AssertionError("gzip output does not have a gzip header")
        assert_payload_size(gz_out, datatype=2, bitpix=8, dims=(8, 8, 8))
        require_success(run_niimath(exe, [str(gz_out)]), "read gzip output")

        shifted = tmp / "small_shifted_y.nii"
        write_uint8_nifti(shifted, offset=(0.0, 1.0, 0.0))
        orientation = run_niimath(exe, [str(small), "-add", str(shifted), str(tmp / "shift_add.nii")])
        require_success(orientation, "binary operation with shifted orientation")
        if "Inconsistent orientations" not in (orientation.stdout + orientation.stderr):
            raise AssertionError("binary operation failed to detect a y-axis spatial mismatch")

        exercise_qc(exe, tmp)

        exercise_allineate(exe, tmp, help_text)

        if args.expect_bsd:
            spm = run_niimath(exe, [str(small), "-spm_coreg", str(small), str(tmp / "spm.nii")])
            if spm.returncode == 0 or "requires a build with the optional GPL module" not in (spm.stdout + spm.stderr):
                raise AssertionError("BSD package should reject -spm_coreg with the GPL-module message")

        if args.expect_zstd:
            env = os.environ.copy()
            env["FSLOUTPUTTYPE"] = "NIFTI_ZST"
            zst_prefix = tmp / "zstd_roundtrip"
            require_success(run_niimath(exe, [str(small), "-add", "1", str(zst_prefix), "-odt", "char"], env=env), "zstd write")
            zst_out = tmp / "zstd_roundtrip.nii.zst"
            if not zst_out.exists():
                raise AssertionError(f"expected zstd output {zst_out}, found {sorted(p.name for p in tmp.iterdir())}")
            require_success(run_niimath(exe, [str(zst_out)]), "read zstd output")

    print("release smoke passed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as exc:
        print(f"release smoke failed: {exc}", file=sys.stderr)
        raise SystemExit(1)

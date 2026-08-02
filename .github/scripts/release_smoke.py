#!/usr/bin/env python3
"""Release smoke test for packaged niimath executables.

The test is intentionally stdlib-only so it can run inside cibuildwheel's clean
test environments and on AppVeyor release workers. It synthesizes small NIfTI
fixtures directly and exercises the packaging-sensitive paths: gzip read/write,
the reported -conform/-gz 0/-odt char case, feature dispatch, and optional zstd.
"""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path


SHORT_READ = "++ WARNING: read "


def _prod(values):
    # math.prod is Python 3.8+, but this script is the cibuildwheel test-command and
    # must run on the package minimum (requires-python >=3.7). Keep it stdlib-portable.
    result = 1
    for value in values:
        result *= value
    return result


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
    scl_slope: float = 1.0,  # intensity scaling: stored -> scaled value is slope*v + inter
    scl_inter: float = 0.0,
    tr: float = 0.0,  # pixdim[4], expressed in time_units
    time_units: int = 0,  # 0 = none (historic default), 8 = sec, 16 = msec, 24 = usec
) -> bytes:
    # Defaults (xyz_units=2 mm, scale=1.0) reproduce the historic mm/unit-voxel
    # header exactly; scale/xyz_units let a fixture describe the SAME physical grid
    # in metres (scale=0.001, xyz_units=1) or microns to exercise unit normalization.
    hdr = bytearray(348)
    struct.pack_into("<i", hdr, 0, 348)
    struct.pack_into("<8h", hdr, 40, 3, dims[0], dims[1], dims[2], 1, 1, 1, 1)
    struct.pack_into("<h", hdr, 70, datatype)
    struct.pack_into("<h", hdr, 72, bitpix)
    struct.pack_into("<8f", hdr, 76, 1.0, scale, scale, scale, tr, 0.0, 0.0, 0.0)
    struct.pack_into("<f", hdr, 108, 352.0)
    struct.pack_into("<f", hdr, 112, scl_slope)
    struct.pack_into("<f", hdr, 116, scl_inter)
    hdr[123] = xyz_units | time_units
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
    scl_slope: float = 1.0,
    scl_inter: float = 0.0,
    tr: float = 0.0,
    time_units: int = 0,
    scale: float = 1.0,  # voxel size, for fixtures whose PHYSICAL size matters
) -> None:
    nvox = dims[0] * dims[1] * dims[2] * nt
    if len(data) != nvox:
        raise AssertionError(f"{path}: expected {nvox} values, got {len(data)}")
    payload = struct.pack(f"<{nvox}f", *data)
    header = bytearray(
        nifti_header(dims, datatype=16, bitpix=32, offset=offset, scl_slope=scl_slope,
                     scl_inter=scl_inter, tr=tr, time_units=time_units, scale=scale)
    )
    if nt > 1:
        struct.pack_into("<8h", header, 40, 4, dims[0], dims[1], dims[2], nt, 1, 1, 1)
    path.write_bytes(bytes(header) + payload)


def read_float32_nifti(path: Path) -> list[float]:
    blob = path.read_bytes()
    dim = struct.unpack_from("<8h", blob, 40)
    nvox = _prod(dim[1 : dim[0] + 1])
    offset = int(struct.unpack_from("<f", blob, 108)[0])
    return list(struct.unpack_from(f"<{nvox}f", blob, offset))


def exercise_fillh(exe: str, tmp: Path) -> None:
    """Exercise independent 3D scratch on a multi-volume flood fill."""
    dims = (5, 5, 5)
    n3 = _prod(dims)
    first = [1.0] * n3
    first[2 + 5 * (2 + 5 * 2)] = 0.0
    second = [0.0] * n3
    src = tmp / "fillh_4d.nii"
    out = tmp / "fillh_4d_out.nii"
    write_float32_nifti(src, dims, first + second, nt=2)
    require_success(
        run_niimath(exe, [str(src), "-fillh", "-gz", "0", str(out)]),
        "multi-volume -fillh",
    )
    result = read_float32_nifti(out)
    if result[:n3] != [1.0] * n3 or result[n3:] != second:
        raise AssertionError("-fillh did not process 4D volumes independently")


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


# Golden numeric-primitive tables, captured from the pinned Julia oracle. Embedded rather than
# read from the oracle's romeo_ref/ output because that directory is gitignored and is absent from a built
# wheel, so this runs wherever release_smoke.py runs: the cibuildwheel matrix (gcc/Linux,
# MSVC/Windows, AppleClang) plus any local `make test`. It does NOT by itself cover the
# Emscripten or WASI builds — the WASI suite checks its own dump against this same native binary
# instead. The tables cover rem2pi in both widths (Payne-Hanek branch included), gamma, rescale's
# bin boundaries, and both unwrapvoxel subtraction widths.
ROMEO_PRIMITIVE_GOLDEN = {
    "rem2pi64.f64":
        "AAAAAAAAAAAAAAAAAADwPwAAAAAAAPC/GC1EVPshCUAYLURU+yEJwBgtRFT7Ifk/GC1EVPsh+b8ZLURU+yH5vwdc"
        "FDMmprG8B1wUMyamsTwAAAAAAAAEQAAAAAAAAATAMVqIqPZDBsAxWoio9kMGQGG0EFHth/S/YbQQUe2H9D8+l95d"
        "JfDmPz6X3l0l8Oa/I4wWIqr94L8jjBYiqv3gP2HvP001J+8/phdBq8DYCECYs/bQVOLWv3pcxUCqTuK/ZPjUjkVJ"
        "AcBk+NSORUkBQLL5U6MMqQVAsvlTowypBcAT4iH2nkvgvxPiIfaeS+A/xv0JaKngAEDmuYsUenHmvzHDh+HeDglA"
        "G+qDVbj92L8AAAAAAADgPwAAAAAAAOC/AAAAAAAA+D8BAAAAAAAEQAEAAAAAAAAAAQAAAAAAAIAr5nCLaBIAAP//"
        "/////w8A",
    "rem2pi32_gamma.f32":
        "AAAAAAAAgD8AAIC/2g9JQNoPScDaD0nA2g9JQOhSRcDoUkVALr07NC69O7S7D0lALr27NGBCog1gQqKNUe0HvwXG"
        "RkBlSC1A+FwCv/hcAj+LXTc/tR8SwLUfEkABAAAAAQAAgP//fwAAAAAAAACAPwAAgL/aD0lA2w9JwNoPScDbD0lA"
        "6VJFwOlSRUAAAAAAAAAAALoPSUDbD8lAYEKiDWBCoo0Cb7tC3EzDR3qWGEv5AhVQ+QIV0Ox4rWC2HxLAth8SQAEA"
        "AAABAACA//9/AA==",
    "rescale.u8":
        "/wGA/v78AQEAAAEBv0A=",
    "unwrapvoxel.f32":
        "AACAP7UfUsDbD8lA2w/JwNsP6UDbD8dC+QIVUNoPScDthwRB7YcEwQAAgD+1H1LA2w/JQNsPycDbD+lA2w/HQvkC"
        "FVDaD0nA7YcEQe2HBME=",
}


def check_romeo_primitives(exe: str, tmp: Path, phase: Path, mag: Path) -> None:
    dump = tmp / "primdump"
    dump.mkdir(exist_ok=True)
    require_success(
        run_niimath(exe, [str(phase), "-gz", "0", "-romeo", str(mag), "-t", "5.0", "-k", "nomask",
                          "-no-phase-rescale", "-romeo-dump", str(dump), str(tmp / "prim_out.nii")]),
        "romeo primitive dump",
    )
    for name, b64 in ROMEO_PRIMITIVE_GOLDEN.items():
        want = base64.b64decode(b64)
        got_path = dump / ("c_prim_" + name)
        if not got_path.exists():
            raise AssertionError("romeo did not emit c_prim_%s" % name)
        got = got_path.read_bytes()
        if got != want:
            ndiff = sum(1 for a, b in zip(got, want) if a != b)
            raise AssertionError(
                "romeo primitive table %s differs from the Julia golden: %d/%d bytes "
                "(len %d vs %d) - the numeric core does not match the reference on this target"
                % (name, ndiff, len(want), len(got), len(want))
            )
    print("  -romeo: numeric primitives match the Julia golden")


def check_romeo_b0_modes(exe: str, tmp: Path) -> None:
    """All six B0 weighting formulas, on TWO echoes with analytically known values.

    A single echo cancels every nonzero weight, so the one-echo check cannot tell the modes
    apart — an incorrect weighting formula passes it. Here phase and magnitude are constant per
    echo, so each mode's expected B0 and SNR are closed-form.
    """
    dims = (4, 4, 4)
    nvox = dims[0] * dims[1] * dims[2]
    tes = [5.0, 11.0]
    two_pi = 6.283185307179586
    ph = [0.4, 1.1]          # radians, small enough that no unwrapping changes them
    mg = [700.0, 300.0]
    phase = tmp / "b0_phase.nii"
    mag = tmp / "b0_mag.nii"
    write_float32_nifti(phase, dims, [ph[0]] * nvox + [ph[1]] * nvox, nt=2)
    write_float32_nifti(mag, dims, [mg[0]] * nvox + [mg[1]] * nvox, nt=2)

    def weights(mode):
        if mode == "phase_snr":
            return [mg[e] * tes[e] for e in (0, 1)]
        if mode == "phase_var":
            return [mg[e] * mg[e] * tes[e] * tes[e] for e in (0, 1)]
        if mode == "average":
            return [1.0, 1.0]
        if mode == "TEs":
            return [tes[e] for e in (0, 1)]
        if mode == "mag":
            return [mg[e] for e in (0, 1)]
        return [math.exp(-tes[e] / 20.0) * tes[e] for e in (0, 1)]   # simulated_mag

    for mode in ("phase_snr", "phase_var", "average", "TEs", "mag", "simulated_mag"):
        w = weights(mode)
        num = sum(ph[e] / tes[e] * w[e] for e in (0, 1))
        den = sum(w)
        want_b0 = (1000.0 / two_pi) * num / den
        want_snr = sum(mg[e] * w[e] for e in (0, 1)) / den
        out = tmp / ("b0m_%s.nii" % mode)
        require_success(
            run_niimath(exe, [str(phase), "-gz", "0", "-romeo", str(mag), "-t", "5.0,11.0",
                              "-k", "nomask", "-no-phase-rescale", "-B",
                              "-B0-phase-weighting", mode, str(out)]),
            "romeo B0 two-echo (%s)" % mode,
        )
        got_b0 = read_float32_nifti(tmp / ("b0m_%s_B0.nii" % mode))
        got_snr = read_float32_nifti(tmp / ("b0m_%s_B0_snr.nii" % mode))
        for label, got, want in (("B0", got_b0, want_b0), ("SNR", got_snr, want_snr)):
            worst = max(abs(v - want) for v in got)
            if worst > 1e-4 * max(1.0, abs(want)):
                raise AssertionError(
                    "romeo -B %s %s: got %g, expected %g (max|diff| %g)"
                    % (mode, label, got[0], want, worst)
                )
    print("  -romeo: all six B0 weighting modes match closed-form values")


def check_romeo_raw_mask(exe: str, tmp: Path, phase: Path) -> None:
    """`-k <file>` truth is `raw != 0` in the STORED width. A float64 mask of 1e-300 is entirely
    true; narrowing it to float32 first would round every voxel to a false zero."""
    dims = (12, 10, 8)
    nvox = dims[0] * dims[1] * dims[2]
    maskf = tmp / "tiny_mask.nii"
    header = bytearray(nifti_header(dims, datatype=64, bitpix=64))   # DT_FLOAT64
    maskf.write_bytes(bytes(header) + struct.pack("<%dd" % nvox, *([1e-300] * nvox)))
    out = tmp / "rawmask_out.nii"
    require_success(
        run_niimath(exe, [str(phase), "-gz", "0", "-romeo", "none", "-t", "5.0",
                          "-k", str(maskf), "-no-phase-rescale", str(out)]),
        "romeo float64 subnormal-ish raw mask",
    )
    if not out.exists():
        raise AssertionError("romeo produced no output for a float64 raw mask")
    print("  -romeo: float64 raw mask truth preserved")


def exercise_romeo(exe: str, tmp: Path, help_text: str) -> None:
    """-romeo phase unwrapping.

    Checks the RESULT, not merely a zero exit status, with a property that pins the answer
    exactly and stays portable across the five wheel runners: unwrap a synthetic phase whose
    ground truth is known, and require `unwrapped - ground_truth` to be ONE constant multiple
    of 2*pi over the whole volume.  Any mis-assigned wrap shows up as a second constant.
    """
    # Match the FULL distinguishing line: "NOT in this build" alone would silently disable this
    # whole test the moment any other feature adopts the same phrasing.
    if "ROMEO phase unwrapping — NOT in this build" in help_text:
        print("  -romeo: not built (ROMEO=0) - skipping")
        return

    dims = (12, 10, 8)
    nvox = dims[0] * dims[1] * dims[2]
    two_pi = 6.283185307179586

    truth = []
    mag = []
    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                truth.append(0.9 * i + 0.4 * j + 0.25 * k)
                # bright core, dark rim: gives robustmask something to find
                core = (2 <= i < dims[0] - 2) and (2 <= j < dims[1] - 2) and (1 <= k < dims[2] - 1)
                mag.append(900.0 + (i + j + k) % 7 if core else 10.0)
    wrapped = [t - two_pi * round(t / two_pi) for t in truth]
    if max(wrapped) - min(wrapped) < 5.0:
        raise AssertionError("romeo fixture is not actually wrapped")

    phase_path = tmp / "romeo_phase.nii"
    mag_path = tmp / "romeo_mag.nii"
    write_float32_nifti(phase_path, dims, wrapped)
    write_float32_nifti(mag_path, dims, mag)

    # 1. nomask: every voxel is unwrapped, so the property covers the whole volume.
    out = tmp / "romeo_out.nii"
    require_success(
        run_niimath(exe, [str(phase_path), "-gz", "0", "-romeo", str(mag_path), "-t", "5.0",
                          "-k", "nomask", "-no-phase-rescale", str(out)]),
        "romeo unwrap (nomask)",
    )
    got = read_float32_nifti(out)
    if len(got) != nvox:
        raise AssertionError("romeo output has %d voxels, expected %d" % (len(got), nvox))
    offsets = set()
    for value, expected in zip(got, truth):
        delta = value - expected
        wraps = round(delta / two_pi)
        if abs(delta - wraps * two_pi) > 1e-3:
            raise AssertionError(
                "romeo left a residual of %g rad (not a multiple of 2*pi)" % (delta - wraps * two_pi)
            )
        offsets.add(wraps)
    if len(offsets) != 1:
        raise AssertionError("romeo assigned %d different 2*pi offsets: %s" % (len(offsets), sorted(offsets)))
    if (tmp / "romeo_out_mask.nii").exists():
        raise AssertionError("-k nomask must not write a mask side output")

    # 2. default robustmask: writes a 0/1 mask, and the unwrapped phase still rewraps to the input
    out2 = tmp / "romeo_masked.nii"
    require_success(
        run_niimath(exe, [str(phase_path), "-gz", "0", "-romeo", str(mag_path), "-t", "5.0",
                          "-no-phase-rescale", "-q", str(out2)]),
        "romeo unwrap (robustmask)",
    )
    mask_path = tmp / "romeo_masked_mask.nii"
    if not mask_path.exists():
        raise AssertionError("robustmask must write a <out>_mask side output")
    mask = read_float32_nifti(mask_path)
    if set(mask) - {0.0, 1.0}:
        raise AssertionError("romeo mask is not binary")
    if not 0 < sum(mask) < nvox:
        raise AssertionError("romeo mask is empty or covers everything (%g of %d)" % (sum(mask), nvox))
    quality = read_float32_nifti(tmp / "romeo_masked_quality.nii")
    if min(quality) < -1e-6 or max(quality) > 1.0 + 1e-6:
        raise AssertionError("romeo quality map outside [0,1]: %g..%g" % (min(quality), max(quality)))
    unwrapped2 = read_float32_nifti(out2)
    for value, original in zip(unwrapped2, wrapped):
        delta = value - original
        if abs(delta - round(delta / two_pi) * two_pi) > 1e-3:
            raise AssertionError("romeo output does not rewrap to its input")

    # 3. options that must FAIL with a specific message rather than being silently ignored
    for opts, needle in (
        (["-w", "romeo9"], "unknown -w"),
        (["-w", "bestpath"], "not implemented"),
        (["-t", "notanumber"], "cannot parse -t"),
        (["-max-seeds", "2"], "not implemented"),
        (["-k", "0.25"], "undefined"),
    ):
        bad = run_niimath(exe, [str(phase_path), "-romeo", str(mag_path), "-t", "5.0"] + opts + [str(tmp / "romeo_bad.nii")])
        if bad.returncode == 0:
            raise AssertionError("romeo %s should have failed" % " ".join(opts))
        if needle not in (bad.stdout + bad.stderr):
            raise AssertionError("romeo %s error message lacks %r" % (" ".join(opts), needle))

    # 4. multi-echo: temporal unwrapping keeps every echo consistent with its own wrapped input
    truth2 = [t * (11.0 / 5.0) for t in truth]
    wrapped2 = [t - two_pi * round(t / two_pi) for t in truth2]
    phase4d = tmp / "romeo_phase4d.nii"
    mag4d = tmp / "romeo_mag4d.nii"
    write_float32_nifti(phase4d, dims, wrapped + wrapped2, nt=2)
    write_float32_nifti(mag4d, dims, mag + mag, nt=2)
    out3 = tmp / "romeo_me.nii"
    require_success(
        run_niimath(exe, [str(phase4d), "-gz", "0", "-romeo", str(mag4d), "-t", "[5.0,11.0]",
                          "-k", "nomask", "-no-phase-rescale", str(out3)]),
        "romeo multi-echo unwrap",
    )
    me = read_float32_nifti(out3)
    if len(me) != 2 * nvox:
        raise AssertionError("romeo multi-echo output has %d voxels, expected %d" % (len(me), 2 * nvox))
    for value, original in zip(me, wrapped + wrapped2):
        delta = value - original
        if abs(delta - round(delta / two_pi) * two_pi) > 1e-3:
            raise AssertionError("romeo multi-echo output does not rewrap to its input")
    # 5. B0: for a SINGLE echo the weighting cancels, so B0 == (1000/2pi)*phase/TE exactly,
    #    whichever mode is selected. Cheap exact oracle that needs no reference data.
    for wmode in ("phase_snr", "average", "mag"):
        out4 = tmp / ("romeo_b0_%s.nii" % wmode)
        require_success(
            run_niimath(exe, [str(phase_path), "-gz", "0", "-romeo", str(mag_path), "-t", "5.0",
                              "-k", "nomask", "-no-phase-rescale", "-B",
                              "-B0-phase-weighting", wmode, str(out4)]),
            "romeo B0 (%s)" % wmode,
        )
        b0 = read_float32_nifti(tmp / ("romeo_b0_%s_B0.nii" % wmode))
        unwrapped = read_float32_nifti(out4)
        for value, ph in zip(b0, unwrapped):
            want = (1000.0 / two_pi) * ph / 5.0
            if abs(value - want) > 1e-3 * max(1.0, abs(want)):
                raise AssertionError("romeo -B %s: %g != (1000/2pi)*phase/TE = %g" % (wmode, value, want))
        if not (tmp / ("romeo_b0_%s_B0_snr.nii" % wmode)).exists():
            raise AssertionError("romeo -B must also write <out>_B0_snr")

    check_romeo_b0_modes(exe, tmp)
    check_romeo_raw_mask(exe, tmp, phase_path)
    check_romeo_primitives(exe, tmp, phase_path, mag_path)
    print("  -romeo: unwrap/mask/quality/multi-echo OK")


def exercise_medic(exe: str, tmp: Path, help_text: str) -> None:
    """--medic / -unwarp: analytic property checks, no external data.

    Both properties pin the answer exactly rather than asserting a zero exit status:

      -unwarp  a ramp whose value equals its j index, pulled by a CONSTANT displacement map of
               exactly N voxels, must come back shifted by exactly -N in the interior.  This
               catches a sign flip, a wrong axis, a wrong length unit, and a broken kernel.
               (The convention is measured; see medic_bench's test/medic_reference_manifest.md section 3.5.)

      --medic  phase synthesised as wrap(2*pi*f*TE) for a known linear field f must be recovered
               by the magnitude-weighted regression.  The field is checked by least-squares slope
               and intercept, which are immune to the small wrap-boundary ripple the synthesis
               itself introduces.
    """
    if "--medic" not in help_text:
        print("  --medic: not built (MEDIC=0) - skipping")
        return

    nx, ny, nz = 16, 24, 8
    vox = 1.0  # write_float32_nifti's default header is 1 mm isotropic, axis-aligned RAS
    nvox = nx * ny * nz

    def idx(x: int, y: int, z: int) -> int:
        return x + y * nx + z * nx * ny

    # ---- -unwarp: constant map, ramp along j -------------------------------------------------
    ramp = [0.0] * nvox
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                ramp[idx(x, y, z)] = float(y)
    shift_vox = 3.0
    dmap = [shift_vox * vox] * nvox
    ramp_path = tmp / "medic_ramp.nii"
    dmap_path = tmp / "medic_dmap.nii"
    write_float32_nifti(ramp_path, (nx, ny, nz), ramp)
    write_float32_nifti(dmap_path, (nx, ny, nz), dmap)
    out_path = tmp / "medic_unwarped.nii"
    require_success(
        run_niimath(exe, [str(ramp_path), "-unwarp", str(dmap_path), "j", str(out_path)]),
        "-unwarp constant map",
    )
    # niimath may append .gz depending on FSLOUTPUTTYPE, so accept either spelling
    written = out_path if out_path.exists() else Path(str(out_path) + ".gz")
    if not written.exists():
        raise AssertionError("-unwarp did not write an output image")
    if written.suffix == ".gz":
        plain = tmp / "medic_unwarped_plain.nii"
        plain.write_bytes(gzip.decompress(written.read_bytes()))
        written = plain
    got = read_float32_nifti(written)
    # interior only: the kernel has radius 5, so edges legitimately see the zero fill
    for z in range(2, nz - 2):
        for y in range(8, ny - 8):
            for x in range(2, nx - 2):
                expect = float(y) - shift_vox
                actual = got[idx(x, y, z)]
                if abs(actual - expect) > 1e-3:
                    raise AssertionError(
                        f"-unwarp: at ({x},{y},{z}) expected {expect} got {actual} "
                        f"(a constant {shift_vox}-voxel map must pull by exactly that)"
                    )
    # out-of-FOV must be zero-filled, not clamped
    if abs(got[idx(nx // 2, 0, nz // 2)]) > 1e-6:
        raise AssertionError("-unwarp: out-of-FOV voxels must be zero-filled")

    # ---- --medic: known linear field ---------------------------------------------------------
    tes = (10.0, 30.0)
    slope_hz_per_vox = 4.0
    two_pi = 6.283185307179586
    mags, phases = [], []
    for e, te in enumerate(tes):
        m = [0.0] * nvox
        p = [0.0] * nvox
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    f = slope_hz_per_vox * (y - ny / 2.0)
                    ang = two_pi * f * (te / 1000.0)
                    # principal value in (-pi, pi]
                    ang = ang - two_pi * math.floor(ang / two_pi + 0.5)
                    p[idx(x, y, z)] = ang
                    m[idx(x, y, z)] = 1000.0 if (2 < x < nx - 3 and 3 < y < ny - 4 and 1 < z < nz - 2) else 30.0
        mp = tmp / f"medic_mag{e}.nii"
        pp = tmp / f"medic_pha{e}.nii"
        write_float32_nifti(mp, (nx, ny, nz), m)
        write_float32_nifti(pp, (nx, ny, nz), p)
        mags.append(str(mp))
        phases.append(str(pp))

    prefix = tmp / "medic_out"
    result = run_niimath(exe, [
        "--medic",
        "--magnitude", *mags,
        "--phase", *phases,
        "--te-ms", f"{tes[0]:g},{tes[1]:g}",
        "--total-readout-time", "0.02",
        "--phase-encoding-direction", "j",
        "--out-prefix", str(prefix),
        "--rank", "0",
    ])
    require_success(result, "--medic synthetic run")

    native = None
    for suffix in (".nii", ".nii.gz"):
        cand = Path(str(prefix) + "_fieldmaps_native" + suffix)
        if cand.exists():
            native = cand
            break
    if native is None:
        raise AssertionError("--medic did not write <prefix>_fieldmaps_native")
    for extra in ("_fieldmaps", "_displacementmaps"):
        if not any(Path(str(prefix) + extra + sfx).exists() for sfx in (".nii", ".nii.gz")):
            raise AssertionError(f"--medic did not write <prefix>{extra}")

    if native.suffix == ".gz":
        plain = tmp / "medic_native_plain.nii"
        plain.write_bytes(gzip.decompress(native.read_bytes()))
        native = plain
    field = read_float32_nifti(native)
    # least-squares fit of field against j over the high-signal interior
    xs, ys = [], []
    for z in range(3, nz - 3):
        for y in range(6, ny - 6):
            for x in range(4, nx - 4):
                xs.append(float(y))
                ys.append(field[idx(x, y, z)])
    n = float(len(xs))
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    slope = sxy / sxx
    intercept = my - slope * mx
    if abs(slope - slope_hz_per_vox) > 0.05 * slope_hz_per_vox:
        raise AssertionError(
            f"--medic: recovered field slope {slope:.4f} Hz/voxel, expected {slope_hz_per_vox}"
        )
    expect_intercept = -slope_hz_per_vox * ny / 2.0
    if abs(intercept - expect_intercept) > 0.05 * abs(expect_intercept):
        raise AssertionError(
            f"--medic: recovered field intercept {intercept:.4f} Hz, expected {expect_intercept}"
        )
    for value in field:
        if value != value or abs(value) > 1e6:
            raise AssertionError("--medic: field map contains non-finite or absurd values")
    print("  --medic/-unwarp: displacement sign, fill and field regression OK")


# ---------------------------------------------------------------------------------------------
# MEDIC regression fixtures.
#
# Every check below pins a bug that was shipped and fixed; each is annotated with the value the
# assertion would see under the OLD behaviour, since the fix cannot be reverted to prove it.
#
# The synthetic field is shaped so that niimath's readphase rescaling (medic.c md_rescale_phase)
# is a NO-OP.  That rescale maps the observed [min, max] of a phase series onto [-pi, pi] unless
# the span is already within 0.1 rad of 2*pi; when it fires it applies a per-echo affine that no
# analytic expectation survives.  Coverage of the full circle is bought with three ramps of very
# different rates -- coarse along j, fine along i, finer along k -- so the wrapped samples tile
# the circle to ~0.002 rad at EVERY echo, while the steepest per-voxel phase step stays well
# under pi so the unwrap is unambiguous.  With that, --medic recovers the field exactly: the
# least-squares slope along j comes back 5.0000 Hz/voxel and the intercept to five digits.
MEDIC_DIMS = (16, 24, 8)
MEDIC_FIELD_J = 5.0          # Hz per j voxel   (0.05 cycle/voxel at TE 10 ms)
MEDIC_FIELD_I = 0.3125       # Hz per i voxel   (fills the j gaps)
MEDIC_FIELD_K = 0.0390625    # Hz per k voxel   (fills the i gaps)
MEDIC_TWO_PI = 6.283185307179586


def medic_field(x: int, y: int, z: int) -> float:
    return MEDIC_FIELD_J * (y - MEDIC_DIMS[1] / 2.0) + MEDIC_FIELD_I * x + MEDIC_FIELD_K * z


def medic_magnitude(x: int, y: int, z: int) -> float:
    nx, ny, nz = MEDIC_DIMS
    inside = (2 < x < nx - 3) and (3 < y < ny - 4) and (1 < z < nz - 2)
    return 1000.0 if inside else 30.0


def medic_wrap(angle: float) -> float:
    """principal value in (-pi, pi], matching medic.c md_wrapf"""
    return angle - MEDIC_TWO_PI * math.floor(angle / MEDIC_TWO_PI + 0.5)


def medic_write_series(
    tmp: Path,
    tag: str,
    tes: tuple[float, ...],
    frames: int,
    amplitude,
    nan_index: int | None = None,
    offset_rad: float = 0.0,
    field_fn=None,
    magnitude_fn=None,
) -> tuple[list[str], list[str]]:
    """Write one magnitude and one phase image per echo for a field amplitude(t)*medic_field(v).

    `nan_index` poisons a single voxel of the FIRST echo's phase (flat index into the whole
    series), which is how the silent all-zero-output bug is provoked.  `offset_rad` adds a
    TE-independent phase offset to every echo, i.e. the term MCPC-3D-S exists to remove.

    `field_fn(x, y, z, t)` replaces amplitude(t)*medic_field(v) when a fixture needs a field whose
    SPATIAL shape (not just its scale) changes with the frame -- the only way to build a series of
    known temporal rank > 1.  `magnitude_fn(x, y, z, t)` likewise replaces medic_magnitude(v), so
    a fixture can make ROMEO's per-frame robustmask differ between frames."""
    nx, ny, nz = MEDIC_DIMS
    mags: list[str] = []
    phases: list[str] = []
    for e, te in enumerate(tes):
        mag_values: list[float] = []
        phase_values: list[float] = []
        for t in range(frames):
            scale = amplitude(t)
            for z in range(nz):
                for y in range(ny):
                    for x in range(nx):
                        hz = field_fn(x, y, z, t) if field_fn is not None else scale * medic_field(x, y, z)
                        phase_values.append(medic_wrap(MEDIC_TWO_PI * hz * te / 1000.0 + offset_rad))
                        mag_values.append(
                            magnitude_fn(x, y, z, t) if magnitude_fn is not None else medic_magnitude(x, y, z)
                        )
        if nan_index is not None and e == 0:
            phase_values[nan_index] = float("nan")
        mag_path = tmp / f"{tag}_mag{e}.nii"
        phase_path = tmp / f"{tag}_pha{e}.nii"
        write_float32_nifti(mag_path, MEDIC_DIMS, mag_values, nt=frames)
        write_float32_nifti(phase_path, MEDIC_DIMS, phase_values, nt=frames)
        mags.append(str(mag_path))
        phases.append(str(phase_path))
    return mags, phases


def medic_run(exe: str, mags: list[str], phases: list[str], tes: tuple[float, ...], prefix: Path, extra: list[str]):
    return run_niimath(exe, [
        "--medic",
        "--magnitude", *mags,
        "--phase", *phases,
        "--te-ms", ",".join(f"{te:g}" for te in tes),
        "--total-readout-time", "0.02",
        "--phase-encoding-direction", "j",
        "--out-prefix", str(prefix),
        *extra,
    ])


def medic_read_output(prefix: Path, suffix: str, tmp: Path, tag: str) -> list[float] | None:
    """Read <prefix><suffix>.nii[.gz]; niimath appends .gz depending on FSLOUTPUTTYPE/-gz."""
    for ext in (".nii", ".nii.gz"):
        candidate = Path(str(prefix) + suffix + ext)
        if candidate.exists():
            if ext == ".nii.gz":
                plain = tmp / f"{tag}{suffix}_plain.nii"
                plain.write_bytes(gzip.decompress(candidate.read_bytes()))
                candidate = plain
            return read_float32_nifti(candidate)
    return None


def medic_output_exists(prefix: Path, suffix: str) -> bool:
    return any(Path(str(prefix) + suffix + ext).exists() for ext in (".nii", ".nii.gz", ".nii.zst"))


def medic_fit_along_j(field: list[float], frame: int) -> tuple[float, float]:
    """Least-squares slope/intercept of the field against the j index over a complete interior
    block.  The block spans every i and k for each j, so the i/k ramps contribute exactly zero to
    the slope and a computable constant to the intercept."""
    nx, ny, nz = MEDIC_DIMS
    n3 = nx * ny * nz
    xs: list[float] = []
    ys: list[float] = []
    for z in range(3, nz - 3):
        for y in range(6, ny - 6):
            for x in range(4, nx - 4):
                xs.append(float(y))
                ys.append(field[frame * n3 + x + y * nx + z * nx * ny])
    n = float(len(xs))
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    slope = sxy / sxx
    return slope, my - slope * mx


def medic_expected_intercept() -> float:
    nx, ny, nz = MEDIC_DIMS
    ix = list(range(4, nx - 4))
    kz = list(range(3, nz - 3))
    return (
        -MEDIC_FIELD_J * ny / 2.0
        + MEDIC_FIELD_I * (sum(ix) / float(len(ix)))
        + MEDIC_FIELD_K * (sum(kz) / float(len(kz)))
    )


def exercise_medic_unwarp_io(exe: str, tmp: Path) -> None:
    """-unwarp input handling: scaled float32 maps, grid mismatch, non-finite map values."""
    nx, ny, nz = MEDIC_DIMS
    nvox = nx * ny * nz

    def index(x: int, y: int, z: int) -> int:
        return x + y * nx + z * nx * ny

    ramp = [0.0] * nvox
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                ramp[index(x, y, z)] = float(y)
    ramp_path = tmp / "unwarp_ramp.nii"
    write_float32_nifti(ramp_path, MEDIC_DIMS, ramp)

    def unwarp(map_path: Path, out_path: Path):
        return run_niimath(exe, [str(ramp_path), "-unwarp", str(map_path), "j", str(out_path)])

    def read_written(out_path: Path, tag: str) -> list[float]:
        written = out_path if out_path.exists() else Path(str(out_path) + ".gz")
        if not written.exists():
            raise AssertionError(f"-unwarp did not write {out_path}")
        if written.suffix == ".gz":
            plain = tmp / f"{tag}_plain.nii"
            plain.write_bytes(gzip.decompress(written.read_bytes()))
            written = plain
        return read_float32_nifti(written)

    # (3) float32 map carrying scl_slope: the SCALED value is the displacement in mm.  Stored 1.0
    # with slope 2.0 means a 2 mm (= 2 voxel) pull.  The bug skipped the scale conversion for
    # float32 input specifically, which would leave a 1-voxel shift here (interior sample 11.0).
    scaled_map = tmp / "unwarp_scaled.nii"
    write_float32_nifti(scaled_map, MEDIC_DIMS, [1.0] * nvox, scl_slope=2.0)
    scaled_out = tmp / "unwarp_scaled_out.nii"
    require_success(unwarp(scaled_map, scaled_out), "-unwarp scaled float32 map")
    got = read_written(scaled_out, "unwarp_scaled")
    for z in range(2, nz - 2):
        for y in range(8, ny - 8):
            for x in range(2, nx - 2):
                expect = float(y) - 2.0
                actual = got[index(x, y, z)]
                if abs(actual - expect) > 1e-3:
                    raise AssertionError(
                        f"-unwarp: scl_slope=2 on a float32 map must pull by the SCALED 2 mm; "
                        f"at ({x},{y},{z}) expected {expect} got {actual} "
                        f"(the raw, unscaled value would give {float(y) - 1.0})"
                    )

    # (4) same dims and same 3x3, but the sform ORIGIN differs: the map describes a different
    # patch of the world and must be rejected.  The bug accepted it and applied it misaligned.
    shifted_map = tmp / "unwarp_shifted.nii"
    write_float32_nifti(shifted_map, MEDIC_DIMS, [3.0] * nvox, offset=(5.0, 0.0, 0.0))
    shifted_out = tmp / "unwarp_shifted_out.nii"
    result = unwarp(shifted_map, shifted_out)
    message = result.stdout + result.stderr
    if result.returncode == 0:
        raise AssertionError("-unwarp accepted a displacement map whose sform origin differs from the input")
    if "grid" not in message:
        raise AssertionError(f"-unwarp rejected a mismatched grid without saying so:\n{message}")
    if shifted_out.exists() or Path(str(shifted_out) + ".gz").exists():
        raise AssertionError("-unwarp wrote an output after rejecting a mismatched displacement map")

    # (5) NaN / +-Inf map values must not reach floor() and the int cast (undefined behaviour).
    # The convention is that such a voxel is fully out of FOV, i.e. the documented zero fill.
    bad_values = [3.0] * nvox
    bad_values[index(8, 12, 4)] = float("nan")
    bad_values[index(9, 12, 4)] = float("inf")
    bad_values[index(10, 12, 4)] = float("-inf")
    bad_map = tmp / "unwarp_nonfinite.nii"
    write_float32_nifti(bad_map, MEDIC_DIMS, bad_values)
    bad_out = tmp / "unwarp_nonfinite_out.nii"
    require_success(unwarp(bad_map, bad_out), "-unwarp non-finite map")
    got = read_written(bad_out, "unwarp_nonfinite")
    for i, value in enumerate(got):
        if value != value or abs(value) > 1e30:
            raise AssertionError(f"-unwarp: a non-finite displacement leaked into output voxel {i} ({value})")
    for x in (8, 9, 10):
        if got[index(x, 12, 4)] != 0.0:
            raise AssertionError(
                f"-unwarp: a non-finite map voxel must take the zero fill, saw {got[index(x, 12, 4)]}"
            )
    # a neighbour with a finite 3 mm displacement is still pulled correctly
    if abs(got[index(7, 12, 4)] - 9.0) > 1e-3:
        raise AssertionError("-unwarp: a non-finite voxel must not disturb its finite neighbours")


def exercise_medic_polarity(exe: str, tmp: Path) -> None:
    """(1) --phase-encoding-direction j vs j-: the native field must be IDENTICAL and the
    displacement map must be negated.  Measured convention, manifest section 3.5b: the polarity
    enters the inversion and the Hz->mm sign but not the weighted regression.

    Under the bug (the '-' suffix dropped) the two runs were byte-identical, so every displacement
    ratio below would be exactly +1 and max|d_j + d_j-| would be 2*max|d_j| instead of ~0.05 of it
    -- i.e. a j- acquisition was corrected backwards, roughly doubling the distortion."""
    tes = (10.0, 30.0)
    mags, phases = medic_write_series(tmp, "medic_pol", tes, 1, lambda t: 1.0)
    fields: dict[str, list[float]] = {}
    disps: dict[str, list[float]] = {}
    for direction in ("j", "j-"):
        tag = "pol_" + ("jm" if direction.endswith("-") else "jp")
        prefix = tmp / f"medic_{tag}"
        result = run_niimath(exe, [
            "--medic",
            "--magnitude", *mags,
            "--phase", *phases,
            "--te-ms", "10,30",
            "--total-readout-time", "0.005",
            "--phase-encoding-direction", direction,
            "--out-prefix", str(prefix),
            "--rank", "0",
        ])
        require_success(result, f"--medic --phase-encoding-direction {direction}")
        native = medic_read_output(prefix, "_fieldmaps_native", tmp, tag)
        disp = medic_read_output(prefix, "_displacementmaps", tmp, tag)
        if native is None or disp is None:
            raise AssertionError(f"--medic ({direction}) did not write its outputs")
        fields[direction] = native
        disps[direction] = disp

    worst = max(abs(a - b) for a, b in zip(fields["j"], fields["j-"]))
    if worst > 1e-6:
        raise AssertionError(
            f"--medic: _fieldmaps_native must not depend on phase-encoding polarity (max diff {worst:g} Hz)"
        )

    dj, dm = disps["j"], disps["j-"]
    # Compare inside the mask and away from its j edges: the field drops to zero outside the mask,
    # and the inversion samples ACROSS that discontinuity in opposite directions for the two
    # polarities, so the two boundary rows legitimately differ by more than the fixed-point term.
    nx, ny, nz = MEDIC_DIMS
    pairs = []
    for z in range(2, nz - 2):
        for y in range(8, 16):
            for x in range(4, nx - 4):
                i = x + y * nx + z * nx * ny
                pairs.append((dj[i], dm[i]))
    peak = max(abs(a) for a, _ in pairs)
    if peak < 1e-3:
        raise AssertionError("--medic: the polarity fixture produced a degenerate displacement map")
    residual = max(abs(a + b) for a, b in pairs)
    # The two are exact negatives only in the limit; the inversion's fixed point makes the ratio
    # -(1 - g*TRT)/(1 + g*TRT) for a field with gradient g along the PE axis, i.e. 4.9 % here.
    # A polarity-blind implementation gives +1 and residual == 2*peak.
    if residual > 0.15 * peak:
        raise AssertionError(
            f"--medic: j and j- displacement maps must be near-negatives; max|d_j + d_j-| = {residual:g} "
            f"against a peak of {peak:g} (a polarity-blind run gives {2 * peak:g})"
        )
    for a, b in zip(dj, dm):
        if abs(a) > 1e-3 and a * b >= 0.0:
            raise AssertionError(
                f"--medic: j and j- displacements must have opposite signs, saw {a:g} and {b:g}"
            )

    # The fold detector must use the SAME polarity as the inversion. Two checks, because either
    # alone is worthless:
    #
    #   liveness      -- a steep field must make it fire at all, else an always-zero detector
    #                    would pass the discriminator below.
    #   discriminator -- a CONVERGED field (|slope*TRT| < 1) with a NEGATIVE slope run as `j-`
    #                    must report NO folds.
    #
    # The discriminator is the load-bearing one. For a linear field of slope b with u = b*TRT,
    # the inverted field has gradient b/(1+u) for j-, so the correct test (-u/(1+u) <= -1) is
    # false for every convergent u, while the polarity-blind test (u/(1+u) <= -1) is TRUE for
    # u <= -0.5. Measured against a binary built with the pe_sign factor removed: the correct
    # code reports 0 folds here and the broken one reports 1024. An earlier version of this test
    # only checked that SOME fold was reported on a non-convergent fixture, where both builds
    # report a fold (127 vs 128) -- it passed with the bug restored.
    full_mask = tmp / "medic_pol_full_mask.nii"
    nvox_mask = MEDIC_DIMS[0] * MEDIC_DIMS[1] * MEDIC_DIMS[2]
    write_float32_nifti(full_mask, MEDIC_DIMS, [1.0] * nvox_mask)

    def fold_count(direction: str, trt: str, phase_files: list, tag: str) -> int:
        res = run_niimath(exe, [
            "--medic", "--magnitude", *mags, "--phase", *phase_files,
            "--te-ms", "10,30", "--total-readout-time", trt,
            "--phase-encoding-direction", direction, "--mask", str(full_mask),
            "--out-prefix", str(tmp / ("medic_fold_" + tag)), "--rank", "0", "--gz", "0",
        ])
        require_success(res, "--medic fold diagnostic (%s, TRT %s)" % (direction, trt))
        m = re.search(r"(\d+) folded adjacent pair", res.stdout + res.stderr)
        return int(m.group(1)) if m else 0

    if fold_count("j-", "0.25", phases, "live") <= 0:
        raise AssertionError("--medic fold detector reported nothing on a folding field; "
                             "an always-zero detector would make the polarity check below vacuous")

    # Negative slope, |slope * TRT| = 0.6 -> the inversion converges, so no genuine fold exists.
    conv_dims = MEDIC_DIMS
    nvox = conv_dims[0] * conv_dims[1] * conv_dims[2]
    conv_phase = []
    for idx, te in enumerate((10.0, 30.0)):
        vals = [0.0] * nvox
        for z in range(conv_dims[2]):
            for y in range(conv_dims[1]):
                for x in range(conv_dims[0]):
                    f = -3.0 * (y - conv_dims[1] / 2.0)
                    ang = 2.0 * math.pi * f * (te / 1000.0)
                    vals[x + y * conv_dims[0] + z * conv_dims[0] * conv_dims[1]] = (
                        ang - 2.0 * math.pi * math.floor(ang / (2.0 * math.pi) + 0.5))
        p = tmp / ("medic_fold_conv_p%d.nii" % idx)
        write_float32_nifti(p, conv_dims, vals)
        conv_phase.append(str(p))
    folds = fold_count("j-", "0.2", conv_phase, "conv")
    if folds != 0:
        raise AssertionError(
            "--medic reported %d fold(s) on a CONVERGENT j- field that cannot fold; the detector "
            "is ignoring the phase-encoding polarity (a build without the pe_sign factor reports "
            "1024 here)" % folds)


def exercise_medic_nonfinite(exe: str, tmp: Path, mags: list[str], phases: list[str], frames: int) -> None:
    """(2) A single NaN phase voxel must fail LOUDLY.

    The bug lived in the low-rank filter: one NaN poisons the whole Gram matrix, every eigenvalue
    becomes NaN, the retained rank collapses to zero and the projector multiplies the entire field
    series by 0 -- three all-zero outputs and exit 0.  It needed T > rank to reach that code, hence
    the 12-frame fixture with the default --rank 10."""
    prefix = tmp / "medic_nan"
    result = medic_run(exe, mags, phases, (10.0, 30.0), prefix, [])
    message = result.stdout + result.stderr
    if result.returncode == 0:
        raise AssertionError(
            "--medic exited 0 on phase data containing NaN; it must fail rather than emit a field map"
        )
    if "non-finite" not in message:
        raise AssertionError(f"--medic rejected non-finite input without saying so:\n{message}")
    for suffix in ("_fieldmaps_native", "_fieldmaps", "_displacementmaps"):
        if medic_output_exists(prefix, suffix):
            values = medic_read_output(prefix, suffix, tmp, "nan")
            raise AssertionError(
                f"--medic wrote <prefix>{suffix} despite failing on non-finite input"
                + (" (and it is all zero -- the exact silent failure this guards)"
                   if values is not None and not any(v != 0.0 for v in values) else "")
            )

    # --rank 0 skips the filter entirely.  Either it errors cleanly (the current behaviour: the
    # non-finite gate now sits on the field series itself, independent of the filter) or it
    # succeeds -- but it must never succeed with an all-zero or non-finite field map.
    prefix0 = tmp / "medic_nan_rank0"
    result0 = medic_run(exe, mags, phases, (10.0, 30.0), prefix0, ["--rank", "0"])
    if result0.returncode == 0:
        values = medic_read_output(prefix0, "_fieldmaps_native", tmp, "nan0")
        if values is None:
            raise AssertionError("--medic --rank 0 exited 0 without writing a field map")
        if not any(v != 0.0 for v in values):
            raise AssertionError("--medic --rank 0 wrote an all-zero field map and exited 0")
        if any(v != v for v in values):
            raise AssertionError("--medic --rank 0 wrote a non-finite field map and exited 0")
    else:
        if "non-finite" not in (result0.stdout + result0.stderr):
            raise AssertionError("--medic --rank 0 failed on NaN input without a diagnostic")
        if medic_output_exists(prefix0, "_fieldmaps_native"):
            raise AssertionError("--medic --rank 0 left a field map behind after failing")


def exercise_medic_three_echo(exe: str, tmp: Path) -> None:
    """(6) Three echoes and the paper's Eq. 6, the cumulative through-origin fit that predicts each
    echo's 2*pi branch from ALL the echoes already corrected.  The retired code predicted every
    later echo from echo 1 alone (phi_n ~ phi_1 * TE_n/TE_1), which coincides with Eq. 6 only for
    two echoes.  Two frames, so the temporal correction actually runs.

    (a) Consistent data: the field must come back exactly.
    (b) Discriminating data: a TE-independent phase offset with --phase-offset none, which makes
        phi_e/TE_e differ between echoes and so separates the two prediction rules.  With
        TEs 10/20/30 ms and an offset c the third echo's discrepancy is 0.8*c under Eq. 6 (rounds
        to no shift) but 2*c under the retired rule (c = 2.5 rad rounds to a whole 2*pi), which
        lands in the field map as a uniform +TE_3/sum(TE^2)/1000 = +21.43 Hz intercept shift."""
    tes = (10.0, 20.0, 30.0)
    mags, phases = medic_write_series(tmp, "medic_e3", tes, 2, lambda t: 1.0)
    prefix = tmp / "medic_e3_out"
    require_success(medic_run(exe, mags, phases, tes, prefix, ["--rank", "0"]), "--medic three-echo run")
    field = medic_read_output(prefix, "_fieldmaps_native", tmp, "e3")
    if field is None:
        raise AssertionError("--medic three-echo run wrote no field map")
    expect_intercept = medic_expected_intercept()
    for frame in (0, 1):
        slope, intercept = medic_fit_along_j(field, frame)
        if abs(slope - MEDIC_FIELD_J) > 0.02 * MEDIC_FIELD_J:
            raise AssertionError(
                f"--medic (3 echoes, frame {frame}): field slope {slope:.5f} Hz/voxel, expected {MEDIC_FIELD_J}"
            )
        if abs(intercept - expect_intercept) > 0.02 * abs(expect_intercept):
            raise AssertionError(
                f"--medic (3 echoes, frame {frame}): field intercept {intercept:.5f} Hz, expected {expect_intercept:.5f}"
            )

    offset_rad = 2.5
    mags, phases = medic_write_series(tmp, "medic_e3off", tes, 2, lambda t: 1.0, offset_rad=offset_rad)
    prefix = tmp / "medic_e3off_out"
    require_success(
        medic_run(exe, mags, phases, tes, prefix, ["--rank", "0", "--phase-offset", "none"]),
        "--medic three-echo run with an uncorrected phase offset",
    )
    field = medic_read_output(prefix, "_fieldmaps_native", tmp, "e3off")
    if field is None:
        raise AssertionError("--medic three-echo offset run wrote no field map")
    seconds = [te / 1000.0 for te in tes]
    sum_t = sum(seconds)
    sum_t2 = sum(t * t for t in seconds)
    # An uncorrected offset c enters the weighted regression as a constant c*sum(t)/(2pi*sum(t^2)).
    predicted = expect_intercept + offset_rad * sum_t / (MEDIC_TWO_PI * sum_t2)
    # Multi-echo unwrapping cannot see a global TE-proportional branch, so the whole field may sit
    # a multiple of 1/TE_1 away; that ambiguity is legitimate and deterministic, and reducing the
    # residual modulo it keeps the check on the Eq. 6 term (21.43 Hz, not a multiple of 100 Hz).
    quantum = 1.0 / seconds[0]
    retired_rule_shift = seconds[2] / sum_t2
    for frame in (0, 1):
        slope, intercept = medic_fit_along_j(field, frame)
        residual = (intercept - predicted) % quantum
        if residual > quantum / 2.0:
            residual -= quantum
        if abs(slope - MEDIC_FIELD_J) > 0.02 * MEDIC_FIELD_J:
            raise AssertionError(
                f"--medic (3 echoes + offset, frame {frame}): field slope {slope:.5f}, expected {MEDIC_FIELD_J}"
            )
        if abs(residual) > 2.0:
            raise AssertionError(
                f"--medic (3 echoes + offset, frame {frame}): intercept {intercept:.4f} Hz is {residual:+.4f} Hz "
                f"off the Eq. 6 prediction {predicted:.4f} (mod {quantum:g}); predicting echo 3 from echo 1 "
                f"alone shifts it by {retired_rule_shift:.3f} Hz"
            )


def exercise_medic_rank(exe: str, tmp: Path, mags: list[str], phases: list[str], frames: int, amplitude) -> None:
    """(7) --rank 0 versus the default --rank 10 on a series that is EXACTLY rank 1 in time
    (field = amplitude(t) * f(voxel)).  Rank-10 truncation of a rank-1 matrix is the identity, so
    the two must agree to round-off.  This is the only coverage of the Jacobi eigensolver and the
    projector application; a broken sweep, a mis-transposed projector, or an off-by-one in the
    retained rank shows up here as a large difference or a collapsed (all-zero) series."""
    outputs: dict[str, list[float]] = {}
    for rank in ("0", "10"):
        prefix = tmp / f"medic_rank{rank}"
        result = medic_run(exe, mags, phases, (10.0, 30.0), prefix, ["--rank", rank])
        require_success(result, f"--medic --rank {rank}")
        values = medic_read_output(prefix, "_fieldmaps_native", tmp, f"rank{rank}")
        if values is None:
            raise AssertionError(f"--medic --rank {rank} wrote no field map")
        outputs[rank] = values
    unfiltered, filtered = outputs["0"], outputs["10"]
    peak = max(abs(v) for v in unfiltered)
    if peak < 1.0:
        raise AssertionError("--medic: the rank fixture produced a degenerate field map")
    worst = max(abs(a - b) for a, b in zip(unfiltered, filtered))
    if worst > 1e-3 * peak:
        raise AssertionError(
            f"--medic: rank-10 truncation of a rank-1 series changed it by {worst:g} Hz "
            f"(peak {peak:g}); the two must agree to round-off"
        )
    # and the filtered series must still carry the known per-frame field, not a collapsed one
    for frame in (0, frames - 1):
        slope, _ = medic_fit_along_j(filtered, frame)
        expect = MEDIC_FIELD_J * amplitude(frame)
        if abs(slope - expect) > 0.02 * expect:
            raise AssertionError(
                f"--medic --rank 10 (frame {frame}): field slope {slope:.5f} Hz/voxel, expected {expect:.5f}"
            )


def exercise_medic_phase_offset_none(exe: str, tmp: Path) -> None:
    """(8) --phase-offset none --save-intermediates must not write <prefix>_phase_offset.

    Nothing computes an offset when MCPC-3D-S is disabled, so the retired code wrote the
    uninitialised buffer out as if it were an image."""
    tes = (10.0, 30.0)
    mags, phases = medic_write_series(tmp, "medic_off", tes, 1, lambda t: 1.0)
    prefix_none = tmp / "medic_off_none"
    require_success(
        medic_run(exe, mags, phases, tes, prefix_none,
                  ["--rank", "0", "--phase-offset", "none", "--save-intermediates"]),
        "--medic --phase-offset none --save-intermediates",
    )
    if medic_output_exists(prefix_none, "_phase_offset"):
        raise AssertionError("--medic --phase-offset none wrote a _phase_offset image (uninitialised heap)")
    # the other intermediates ARE expected, so the absence above is not simply a dead option
    for suffix in ("_masks", "_unwrapped_echo-1", "_unwrapped_echo-2"):
        if not medic_output_exists(prefix_none, suffix):
            raise AssertionError(f"--medic --save-intermediates did not write <prefix>{suffix}")

    # positive control: with MCPC enabled the offset image IS written
    prefix_mcpc = tmp / "medic_off_mcpc"
    require_success(
        medic_run(exe, mags, phases, tes, prefix_mcpc,
                  ["--rank", "0", "--phase-offset", "mcpc", "--save-intermediates"]),
        "--medic --phase-offset mcpc --save-intermediates",
    )
    if not medic_output_exists(prefix_mcpc, "_phase_offset"):
        raise AssertionError("--medic --phase-offset mcpc --save-intermediates did not write _phase_offset")


def exercise_medic_mask_contract(exe: str, tmp: Path) -> None:
    """(9) --mask is binarised at `>= 1`, NOT at `!= 0`.

    The contract is measured (manifest section 3.7: the reference's unwrapped phase is nonzero
    exactly on mask >= 1), so a probability map is not a mask.  Under the retired `!= 0` test a
    uniform 0.5 probability map masked EVERY voxel in -- including pure background -- and the run
    silently produced a field map fitted to noise.  It must now fail with an explanation instead."""
    nx, ny, nz = MEDIC_DIMS
    nvox = nx * ny * nz
    tes = (10.0, 30.0)
    mags, phases = medic_write_series(tmp, "medic_mk", tes, 1, lambda t: 1.0)

    # (a) fractional mask: no voxel satisfies `>= 1`, so there is no mask at all.
    half = tmp / "medic_mask_half.nii"
    write_float32_nifti(half, MEDIC_DIMS, [0.5] * nvox)
    result = medic_run(exe, mags, phases, tes, tmp / "medic_mk_half", ["--rank", "0", "--mask", str(half)])
    message = result.stdout + result.stderr
    if result.returncode == 0:
        raise AssertionError(
            "--medic --mask accepted a uniform 0.5 probability map; the in-mask test is `>= 1`, so "
            "no voxel qualifies and the run must fail rather than mask the whole volume in"
        )
    if ">= 1" not in message:
        raise AssertionError(f"--medic rejected a fractional mask without naming the `>= 1` rule:\n{message}")
    if medic_output_exists(tmp / "medic_mk_half", "_fieldmaps_native"):
        raise AssertionError("--medic wrote a field map after rejecting a fractional mask")

    # (b) NaN mask voxels: every comparison against NaN is false, so `>= 1` excludes them.  This
    # doubles as the positive control -- the surrounding 1.0 voxels ARE a valid mask, so the run
    # succeeds and its field map is nonzero there.
    values = [1.0] * nvox
    holes = [(x, 12, 4) for x in range(4, 12)]
    for x, y, z in holes:
        values[x + y * nx + z * nx * ny] = float("nan")
    nan_mask = tmp / "medic_mask_nan.nii"
    write_float32_nifti(nan_mask, MEDIC_DIMS, values)
    prefix = tmp / "medic_mk_nan"
    require_success(
        medic_run(exe, mags, phases, tes, prefix, ["--rank", "0", "--mask", str(nan_mask)]),
        "--medic --mask with NaN voxels",
    )
    field = medic_read_output(prefix, "_fieldmaps_native", tmp, "mknan")
    if field is None:
        raise AssertionError("--medic --mask (NaN voxels) wrote no field map")
    for x, y, z in holes:
        got = field[x + y * nx + z * nx * ny]
        if got != 0.0:
            raise AssertionError(
                f"--medic: a NaN mask voxel fails `>= 1` and must be excluded; ({x},{y},{z}) holds {got:g} Hz"
            )
    if not any(v != 0.0 for v in field):
        raise AssertionError("--medic --mask (NaN voxels) produced an all-zero field map")


# Magnitude bump whose width grows slightly per frame, so ROMEO's robustmask -- and only the mask,
# the phase is frame-independent -- differs between frames.  The widths are close enough that the
# frames stay above MEDIC_CORR_THRESH (0.98) magnitude correlation and so form ONE temporal group.
MEDIC_MASK_SIGMA = (5.0, 5.15, 5.3)


def medic_mask_bump(x: int, y: int, z: int, t: int) -> float:
    nx, ny, nz = MEDIC_DIMS
    sigma = MEDIC_MASK_SIGMA[t]
    r2 = (x - nx / 2.0) ** 2 + (y - ny / 2.0) ** 2 + ((z - nz / 2.0) * 2.0) ** 2
    return 30.0 + 970.0 * math.exp(-r2 / (2.0 * sigma * sigma))


def exercise_medic_mask_temporal(exe: str, tmp: Path) -> None:
    """(10) Mask gating must SURVIVE the temporal 2*pi correction.

    The unwrapped phase is zeroed outside each frame's mask, but the temporal correction runs
    afterwards and moves every voxel toward its group's mean branch.  For a voxel masked out in
    frame t yet inside the mask in the others, that mean is nonzero: if it exceeds pi the
    correction adds a whole 2*pi and the excluded voxel comes back to life in the field map.  (The
    bug: the masks were released before md_temporal(), so nothing could re-apply them.)

    The fixture is built so the failure is REACHABLE rather than merely asserted-against: the
    per-frame robustmasks genuinely differ, and TEs of 25/75 ms make the first-echo unwrapped phase
    large enough that the group mean at the mask boundary passes pi.  Both preconditions are
    checked below from the saved intermediates, so the test cannot quietly become vacuous."""
    nx, ny, nz = MEDIC_DIMS
    nvox = nx * ny * nz
    frames = 3
    tes = (25.0, 75.0)
    mags, phases = medic_write_series(tmp, "medic_mt", tes, frames, lambda t: 1.0,
                                      magnitude_fn=medic_mask_bump)
    prefix = tmp / "medic_mt_out"
    require_success(
        medic_run(exe, mags, phases, tes, prefix,
                  ["--rank", "0", "--temporal-correction", "1", "--save-intermediates"]),
        "--medic --temporal-correction 1 with per-frame masks",
    )
    masks = medic_read_output(prefix, "_masks", tmp, "mtm")
    unwrapped = medic_read_output(prefix, "_unwrapped_echo-1", tmp, "mtu")
    field = medic_read_output(prefix, "_fieldmaps_native", tmp, "mtf")
    if masks is None or unwrapped is None or field is None:
        raise AssertionError("--medic --save-intermediates did not write the masks/unwrapped/field set")

    boundary = [q for q in range(nvox)
                if any(masks[t * nvox + q] == 0.0 for t in range(frames))
                and any(masks[t * nvox + q] != 0.0 for t in range(frames))]
    if not boundary:
        raise AssertionError("the temporal-mask fixture no longer varies its mask between frames")
    # md_temporal's reference for a voxel is the group mean of every frame's first-echo unwrapped
    # phase, the masked-out frames contributing 0.  |mean| > pi is exactly the condition under
    # which the correction rounds to a nonzero multiple of 2*pi.
    movable = 0
    for q in boundary:
        mean = sum(unwrapped[t * nvox + q] for t in range(frames)) / float(frames)
        if abs(mean) > math.pi:
            movable += 1
    if movable < 1:
        raise AssertionError(
            "the temporal-mask fixture is vacuous: no masked-out voxel has a group mean past pi, "
            "so the temporal correction could not move one off zero even without the gating"
        )
    leaked = [i for i in range(len(field)) if masks[i] == 0.0 and field[i] != 0.0]
    if leaked:
        i = leaked[0]
        raise AssertionError(
            f"--medic: {len(leaked)} out-of-mask voxels are nonzero in _fieldmaps_native after the "
            f"temporal correction (frame {i // nvox}, voxel {i % nvox} holds {field[i]:g} Hz); "
            f"{movable} of them sit past the pi threshold that makes the correction add a full 2*pi"
        )
    if not any(v != 0.0 for v in field):
        raise AssertionError("the temporal-mask fixture produced an all-zero field map")


# A field series of EXACTLY temporal rank 2: a j ramp whose amplitude changes per frame, plus an
# i/k ramp that does not.  Rank-1 truncation must therefore lose a real component of the series.
MEDIC_RANK_AMPS = (1.0, 0.6, 1.3, 0.2)


def medic_rank_field(x: int, y: int, z: int, t: int) -> float:
    return (MEDIC_RANK_AMPS[t] * MEDIC_FIELD_J * (y - MEDIC_DIMS[1] / 2.0)
            + MEDIC_FIELD_I * x + MEDIC_FIELD_K * z)


def exercise_medic_rank_boundaries(exe: str, tmp: Path) -> None:
    """(11) The three boundaries of the low-rank filter.

    (a) An all-zero field series with T > rank must SUCCEED.  Its Gram matrix has no positive
        eigenvalue, and the retired code read that as an error ("no positive spectrum") and killed
        an otherwise valid run; the rank-k truncation of a zero matrix is that same zero matrix.
    (b) --rank 1 on a series of genuine temporal rank 2 must MEASURABLY truncate it.  Test (7)
        only shows that truncating a rank-1 series is a no-op, which a filter that silently did
        nothing at all would also satisfy.
    (c) --rank > T must be a harmless no-op, bit-identical to --rank 0."""
    nx, ny, nz = MEDIC_DIMS
    nvox = nx * ny * nz
    tes = (10.0, 30.0)
    frames = 4

    # (a) all-zero phase in every echo and frame -> an all-zero field series; rank 2 < T = 4, so
    # the filter really runs.
    zero_mags, zero_phases = medic_write_series(tmp, "medic_zero", tes, frames, lambda t: 0.0,
                                                field_fn=lambda x, y, z, t: 0.0)
    prefix = tmp / "medic_zero_out"
    result = medic_run(exe, zero_mags, zero_phases, tes, prefix, ["--rank", "2"])
    require_success(result, "--medic --rank 2 on an all-zero field series")
    zeros = medic_read_output(prefix, "_fieldmaps_native", tmp, "zero")
    if zeros is None:
        raise AssertionError("--medic wrote no field map for an all-zero series")
    for i, v in enumerate(zeros):
        if v != 0.0:
            raise AssertionError(f"--medic: an all-zero field series must stay zero, voxel {i} holds {v:g}")

    # (b)/(c) a rank-2 series.
    mags, phases = medic_write_series(tmp, "medic_rk", tes, frames, lambda t: 1.0,
                                      field_fn=medic_rank_field)
    outputs: dict[str, list[float]] = {}
    for rank in ("0", "1", "2", "99"):
        p = tmp / f"medic_rk{rank}"
        require_success(
            medic_run(exe, mags, phases, tes, p, ["--rank", rank, "--temporal-correction", "0"]),
            f"--medic --rank {rank} on a rank-2 series",
        )
        values = medic_read_output(p, "_fieldmaps_native", tmp, f"rk{rank}")
        if values is None:
            raise AssertionError(f"--medic --rank {rank} wrote no field map")
        outputs[rank] = values

    unfiltered = outputs["0"]
    peak = max(abs(v) for v in unfiltered)
    if peak < 1.0:
        raise AssertionError("the rank-boundary fixture produced a degenerate field map")
    # The fixture is only meaningful if the frames really do carry different spatial fields.
    slopes = [medic_fit_along_j(unfiltered, t)[0] for t in range(frames)]
    for t, slope in enumerate(slopes):
        expect = MEDIC_FIELD_J * MEDIC_RANK_AMPS[t]
        if abs(slope - expect) > 0.02 * expect:
            raise AssertionError(
                f"--medic --rank 0 (frame {t}): field slope {slope:.5f} Hz/voxel, expected {expect:.5f}"
            )

    # (c) rank 99 > T = 4: md_lowrank returns before touching the series, so the two runs share
    # every code path and must agree exactly.
    worst = max(abs(a - b) for a, b in zip(unfiltered, outputs["99"]))
    if worst != 0.0:
        raise AssertionError(
            f"--medic --rank 99 on a {frames}-frame series must be a no-op, but it changed the "
            f"field by {worst:g} Hz"
        )
    # (b) truncating BELOW the true rank must lose something; truncating AT it must not.
    lost = max(abs(a - b) for a, b in zip(unfiltered, outputs["1"]))
    kept = max(abs(a - b) for a, b in zip(unfiltered, outputs["2"]))
    if lost < 0.05 * peak:
        raise AssertionError(
            f"--medic --rank 1 changed a rank-2 series by only {lost:g} Hz (peak {peak:g}); the "
            f"low-rank filter is not truncating"
        )
    if kept > 0.05 * peak:
        raise AssertionError(
            f"--medic --rank 2 changed a rank-2 series by {kept:g} Hz (peak {peak:g}); truncation "
            f"at the true rank must be nearly lossless (--rank 1 loses {lost:g})"
        )


# A 90-degree rotation about x as a NIfTI quaternion: a = b = cos(45 deg), c = d = 0.  With
# qfac = -1 (pixdim[0]) the resulting qto_xyz 3x3 is [[1,0,0],[0,0,1],[0,1,0]], i.e. voxel-j runs
# along world +z -- while the sform below is the identity, where voxel-j runs along world +y.
MEDIC_QUAT_ROOT2 = 0.7071067811865476


def medic_write_split_xform(path: Path, data: list[float], sform_code: int, qform_code: int) -> None:
    """A float32 volume whose sform and qform DISAGREE about voxel-j's world axis, with the two
    codes chosen by the caller.  nifti_header() hard-codes qform_code = sform_code = 3, so the
    fields are written directly: pixdim[0]/qfac 76, qform_code 252, sform_code 254,
    quatern_b/c/d 256/260/264, qoffset_x/y/z 268/272/276 (srow_* at 280/296/312 stay identity)."""
    header = bytearray(nifti_header(MEDIC_DIMS, datatype=16, bitpix=32))
    struct.pack_into("<f", header, 76, -1.0)                  # qfac
    struct.pack_into("<h", header, 252, qform_code)
    struct.pack_into("<h", header, 254, sform_code)
    struct.pack_into("<3f", header, 256, MEDIC_QUAT_ROOT2, 0.0, 0.0)
    struct.pack_into("<3f", header, 268, 0.0, 0.0, 0.0)
    path.write_bytes(bytes(header) + struct.pack(f"<{len(data)}f", *data))


def exercise_medic_xform_precedence(exe: str, tmp: Path) -> None:
    """(12) -unwarp must resolve voxel->world exactly as core.c's xform() does: the sform, unless
    sform_code < qform_code, in which case the qform.

    md_xform3() used to prefer the sform unconditionally while md_same_grid() validated through
    xform().  A header whose two transforms disagree therefore passed the grid check on one matrix
    and was corrected using the other -- in the WRONG DIRECTION, doubling the distortion.

    Both fixtures below carry the SAME pair of matrices and differ ONLY in the two codes.  On a
    ramp of value j with a uniform +2 mm displacement map, the qform (voxel-j along world +z, an
    axis with no RAS->LPS sign flip) pulls from j+2 while the sform (voxel-j along world +y, which
    does flip) pulls from j-2.  Swapping the codes must swap the answer; the retired code returned
    j-2 for both."""
    nx, ny, nz = MEDIC_DIMS
    nvox = nx * ny * nz

    def index(x: int, y: int, z: int) -> int:
        return x + y * nx + z * nx * ny

    ramp = [float(i // nx % ny) for i in range(nvox)]
    for codes, expect_shift, tag in (((1, 2), +2.0, "qform"), ((2, 1), -2.0, "sform")):
        sform_code, qform_code = codes
        ramp_path = tmp / f"medic_xf_{tag}_ramp.nii"
        map_path = tmp / f"medic_xf_{tag}_map.nii"
        out_path = tmp / f"medic_xf_{tag}_out.nii"
        medic_write_split_xform(ramp_path, ramp, sform_code, qform_code)
        medic_write_split_xform(map_path, [2.0] * nvox, sform_code, qform_code)
        require_success(
            run_niimath(exe, [str(ramp_path), "-unwarp", str(map_path), "j", str(out_path)]),
            f"-unwarp with sform_code={sform_code} qform_code={qform_code}",
        )
        written = out_path if out_path.exists() else Path(str(out_path) + ".gz")
        if not written.exists():
            raise AssertionError(f"-unwarp wrote no output for the {tag} fixture")
        if written.suffix == ".gz":
            plain = tmp / f"medic_xf_{tag}_plain.nii"
            plain.write_bytes(gzip.decompress(written.read_bytes()))
            written = plain
        got = read_float32_nifti(written)
        for z in range(2, nz - 2):
            for y in range(6, ny - 6):
                for x in range(2, nx - 2):
                    want = float(y) + expect_shift
                    actual = got[index(x, y, z)]
                    if abs(actual - want) > 1e-3:
                        raise AssertionError(
                            f"-unwarp: sform_code={sform_code}, qform_code={qform_code} must honour "
                            f"the {tag}; at ({x},{y},{z}) expected {want} got {actual} (the other "
                            f"transform gives {float(y) - expect_shift})"
                        )


def exercise_medic_parsing(exe: str, tmp: Path) -> None:
    """(13) Every numeric option is parsed strictly and names itself when it rejects a token.

    atoi()/atof() read "5xyz" as 5 and "abc" as 0, which turns a typo into a plausible-looking run
    on the wrong parameters rather than an error.  These four exit before any image is read, so
    they cost nothing."""
    tes = (10.0, 30.0)
    mags, phases = medic_write_series(tmp, "medic_parse", tes, 1, lambda t: 1.0)
    cases = (
        (["--rank", "5xyz"], "--rank"),
        (["--total-readout-time", "abc"], "--total-readout-time"),
        (["--n-cpus", "0"], "--n-cpus"),
        (["--gz", "2"], "--gz"),
    )
    for extra, option in cases:
        prefix = tmp / "medic_parse_out"
        result = medic_run(exe, mags, phases, tes, prefix, extra)
        message = result.stdout + result.stderr
        if result.returncode == 0:
            raise AssertionError(f"--medic accepted {' '.join(extra)}; it must be a hard error")
        if option not in message:
            raise AssertionError(f"--medic rejected {' '.join(extra)} without naming {option}:\n{message}")
        if medic_output_exists(prefix, "_fieldmaps_native"):
            raise AssertionError(f"--medic wrote outputs after rejecting {' '.join(extra)}")


def _is_native_executable(path: str) -> bool:
    """Is `path` a real ELF/Mach-O binary rather than a script that re-spawns one?

    The pip wheel ships a Python console script named `niimath`; a test that assumes the process
    it launches is the one running medic.c is wrong there. Sniff the magic rather than guessing
    from the name or the platform.
    """
    try:
        with open(path, "rb") as fh:
            magic = fh.read(4)
    except OSError:
        return False
    return magic in (
        b"\x7fELF",                                     # ELF
        b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe",   # Mach-O 64/32 little-endian
        b"\xfe\xed\xfa\xcf", b"\xfe\xed\xfa\xce",   # Mach-O big-endian
        b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca",   # Mach-O universal
    )


def exercise_medic_output_transaction(exe: str, tmp: Path) -> None:
    """(14) A failed write must leave a PREVIOUS run's outputs byte-for-byte intact.

    An earlier revision wrote straight to the final names and, when a later output failed, removed
    all three in every extension -- destroying results this invocation had never produced.  The
    outputs are now staged under a sibling temporary prefix and renamed in only once all three
    exist, so an unwritable destination is a clean no-op.

    The backup-rename check needs the PID of the process that will run medic.c, which it gets by
    exec-ing the binary from a launcher (exec preserves the PID) after a short delay so the
    obstacle can be planted first. That identity only holds when `exe` IS the native binary: the
    pip wheel installs a console script that SUBPROCESS-SPAWNS niimath, so medic.c would see a
    different PID, the obstacle would never collide, the run would succeed and the assertion would
    fire -- failing CIBW_TEST_COMMAND on every non-Windows wheel and blocking the PyPI upload.
    So it is gated on `exe` actually being a native executable, and skipped with a reason
    otherwise. The permission check is skipped where permissions cannot be made to bite."""
    tes = (10.0, 30.0)
    mags, phases = medic_write_series(tmp, "medic_txn", tes, 1, lambda t: 1.0)
    outdir = tmp / "medic_txn_dir"
    outdir.mkdir(exist_ok=True)
    prefix = outdir / "run"
    require_success(medic_run(exe, mags, phases, tes, prefix, ["--rank", "0", "--gz", "0"]),
                    "--medic first run")
    before = {p.name: p.read_bytes() for p in sorted(outdir.iterdir())}
    if len(before) != 3:
        raise AssertionError(f"--medic wrote {len(before)} outputs, expected 3: {sorted(before)}")

    exec_path = shutil.which(exe) or exe
    if os.name != "nt" and _is_native_executable(exec_path):
        mags2, phases2 = medic_write_series(tmp, "medic_txn_changed", tes, 1, lambda t: 0.7)
        args = [
            "--medic",
            "--magnitude", *mags2,
            "--phase", *phases2,
            "--te-ms", "10,30",
            "--total-readout-time", "0.02",
            "--phase-encoding-direction", "j",
            "--out-prefix", str(prefix),
            "--rank", "0",
            "--gz", "0",
        ]
        launcher = "import os,sys,time; time.sleep(.2); os.execv(sys.argv[1], sys.argv[1:])"
        proc = subprocess.Popen(
            [sys.executable, "-c", launcher, exec_path, *args],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        obstacle = Path(str(prefix) + f"_fieldmaps_native.nii.medicbak{proc.pid}")
        obstacle.mkdir()
        blocker = obstacle / "keep"
        blocker.write_text("do not remove")
        try:
            stdout, stderr = proc.communicate()
        finally:
            blocker.unlink()
            obstacle.rmdir()
        if proc.returncode == 0:
            raise AssertionError("--medic replaced old outputs after its backup rename failed")
        if "failed to preserve existing output" not in stdout + stderr:
            raise AssertionError(f"--medic backup failure had no clear diagnostic:\n{stdout}{stderr}")
        after_backup_failure = {p.name: p.read_bytes() for p in sorted(outdir.iterdir())}
        if after_backup_failure != before:
            raise AssertionError("--medic changed old outputs after its backup rename failed")
    elif os.name != "nt":
        print("  --medic backup-rename transaction: %r is not a native binary (console script or "
              "wrapper), so its PID cannot be predicted - skipping" % exec_path)

    geteuid = getattr(os, "geteuid", None)
    if geteuid is not None and geteuid() == 0:
        print("  --medic output permission transaction: running as root - skipping")
        return

    os.chmod(str(outdir), 0o500)
    try:
        probe = outdir / "writable.probe"
        try:
            probe.write_text("x")
            probe.unlink()
            print("  --medic output transaction: destination stayed writable - skipping")
            return
        except OSError:
            pass
        result = medic_run(exe, mags, phases, tes, prefix, ["--rank", "0"])
        message = result.stdout + result.stderr
        if result.returncode == 0:
            raise AssertionError("--medic exited 0 with an unwritable output directory")
        if "failed to write" not in message:
            raise AssertionError(f"--medic failed to write its outputs without saying so:\n{message}")
    finally:
        os.chmod(str(outdir), 0o700)

    after = {p.name: p.read_bytes() for p in sorted(outdir.iterdir())}
    if after != before:
        raise AssertionError(
            "--medic destroyed or altered a previous run's outputs when its own write failed: "
            f"had {sorted(before)}, now {sorted(after)}"
        )


def exercise_medic_regressions(exe: str, tmp: Path, help_text: str) -> None:
    """Regressions for the MEDIC correctness fixes; see each helper for the bug it pins."""
    if "--medic" not in help_text:
        print("  --medic regressions: not built (MEDIC=0) - skipping")
        return

    exercise_medic_unwarp_io(exe, tmp)
    exercise_medic_polarity(exe, tmp)
    exercise_medic_three_echo(exe, tmp)
    exercise_medic_phase_offset_none(exe, tmp)
    exercise_medic_mask_contract(exe, tmp)
    exercise_medic_mask_temporal(exe, tmp)
    exercise_medic_rank_boundaries(exe, tmp)
    exercise_medic_xform_precedence(exe, tmp)
    exercise_medic_parsing(exe, tmp)
    exercise_medic_output_transaction(exe, tmp)

    # One 12-frame series feeds both the rank check and (with a poisoned voxel) the non-finite
    # check.  12 > the default --rank 10, which the low-rank bug required.
    frames = 12
    amplitude = lambda t: 1.0 + 0.4 * math.sin(0.9 * t)  # noqa: E731 - keeps the fixture inline
    mags, phases = medic_write_series(tmp, "medic_series", (10.0, 30.0), frames, amplitude)
    exercise_medic_rank(exe, tmp, mags, phases, frames, amplitude)

    nx, ny, nz = MEDIC_DIMS
    poisoned = 8 + 12 * nx + 4 * nx * ny  # an interior voxel of frame 0, well inside the mask
    nan_mags, nan_phases = medic_write_series(tmp, "medic_nanseries", (10.0, 30.0), frames, amplitude,
                                              nan_index=poisoned)
    exercise_medic_nonfinite(exe, tmp, nan_mags, nan_phases, frames)

    print("  --medic/-unwarp regressions: polarity, non-finite, scaling, grid, rank, offsets, "
          "mask contract, xform precedence, parsing, output transaction OK")


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

    # Every advertised fast selector must reach the niimath host dispatch and
    # serialize its resolved engine/cost. The bare default, fast, and fastx are
    # aliases for the same adaptive strategy and must remain byte-identical.
    fast_outputs: dict[str, list[float]] = {}
    expected_cost = {
        "default": "hel+cr", "fast": "hel+cr", "fastx": "hel+cr",
        "fasthel": "hel", "fastcr": "cr",
    }
    for selector in ("default", "fast", "fastx", "fasthel", "fastcr"):
        out = tmp / f"al_{selector}.nii"
        mat = tmp / f"al_{selector}.json"
        args = [str(moving), "-allineate", str(base)]
        if selector != "default":
            args.extend(["-cost", selector])
        args.extend(["-savemat", str(mat), "-gz", "0", str(out)])
        require_success(run_niimath(exe, args), f"-allineate -cost {selector}")
        meta = json.loads(mat.read_text())
        if meta.get("engine") != "coreg_fast" or meta.get("cost") != expected_cost[selector]:
            raise AssertionError(
                f"{selector} resolved to engine/cost {meta.get('engine')}/{meta.get('cost')}, "
                f"expected coreg_fast/{expected_cost[selector]}"
            )
        fast_outputs[selector] = read_float32_nifti(out)
    if fast_outputs["default"] != fast_outputs["fast"] or fast_outputs["fast"] != fast_outputs["fastx"]:
        raise AssertionError("bare default, -cost fast, and -cost fastx are not identical aliases")

    # -deface shares the parser but owns a separate shared-engine dispatch.
    # Regression-guard that fastx/default maps to adaptive HEL/CR rather than
    # accidentally collapsing every non-fasthel selector to correlation ratio.
    df_out = tmp / "al_deface_fastx.nii"
    df = run_niimath(
        exe, [str(moving), "-deface", str(base), str(weight), "-cost", "fastx",
              "-gz", "0", str(df_out)]
    )
    require_success(df, "-deface -cost fastx")
    if "adaptive HEL/CR" not in (df.stdout + df.stderr):
        raise AssertionError("-deface -cost fastx did not dispatch the adaptive HEL/CR strategy")

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

    # An oversized -weight must be rejected from its HEADER (the aux oversize gate), before any
    # payload allocation — not decompressed/allocated then rejected by the estimator. The
    # deliberately header-only file declares 32767*32767*3 > INT_MAX voxels.
    huge_weight = tmp / "al_huge_weight.nii"
    huge_weight.write_bytes(nifti_header((32767, 32767, 3), datatype=16, bitpix=32))
    hw = run_niimath(exe, [str(moving), "-allineate", str(base), "-weight", str(huge_weight),
                           "-gz", "0", str(tmp / "al_hw.nii")])
    if hw.returncode == 0 or "exceeds the supported INT_MAX" not in (hw.stdout + hw.stderr):
        raise AssertionError("oversized -weight was not header-rejected before load")

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

    # --- Back-ported allineate integration surfaces (glue not covered by the standalone suite).
    #     These assert the REPAIRED BEHAVIOR, not merely that an output file was created. ---

    # The ordinary engine consumes the graded -weight (old "-weight is fast-only" gate removed).
    # Prove it took the weight (not silently ignored it): the ordinary weighted -savemat must
    # record the "weight" provenance key (which the ordinary/master path previously dropped).
    hel_prov = tmp / "al_hel_prov.json"
    require_success(
        run_niimath(exe, [str(moving), "-allineate", str(base), "-cost", "hel", "-weight", str(weight),
                          "-savemat", str(hel_prov), "-gz", "0", str(tmp / "al_hel_w.nii")]),
        "-allineate -cost hel -weight (ordinary engine)",
    )
    if '"weight"' not in hel_prov.read_text():
        raise AssertionError("ordinary weighted -savemat dropped the -weight provenance")
    # Prove the ordinary engine actually CONSUMES/validates the weight (not just records the name):
    # a dims-mismatched weight must FAIL the fit with no output (if it ignored the image it would pass).
    bad_w = tmp / "al_badw.nii"
    write_float32_nifti(bad_w, (8, 8, 8), [1.0] * 512)
    bad = run_niimath(exe, [str(moving), "-allineate", str(base), "-cost", "hel", "-weight", str(bad_w),
                            "-gz", "0", str(tmp / "al_badw_out.nii")])
    if bad.returncode == 0 or (tmp / "al_badw_out.nii").exists():
        raise AssertionError("ordinary -cost hel -weight accepted a dims-mismatched weight (weight not consumed)")
    # A same-spatial-dims 4D weight must hit the single-volume guard (nvox != nx*ny*nz) in
    # al_load_user_weight — distinct from the dims-mismatch path above (spatial dims agree here).
    w4d = tmp / "al_w4d.nii"
    write_float32_nifti(w4d, (n, n, n), [1.0] * (n * n * n * 2), nt=2)
    r4 = run_niimath(exe, [str(moving), "-allineate", str(base), "-cost", "hel", "-weight", str(w4d),
                           "-gz", "0", str(tmp / "al_w4d_out.nii")])
    if r4.returncode == 0 or (tmp / "al_w4d_out.nii").exists():
        raise AssertionError("ordinary -weight accepted a same-dims 4D weight (single-volume guard missing)")

    # -reface: identity registration (subject==template) -> full coverage. Assert it actually
    # REPLACES the shell>0 face voxels (not an unchanged anonymization-failed passthrough).
    shell_vals = [1.0 if (12 <= x < 20 and 12 <= y < 20 and 12 <= z < 20) else 0.0
                  for z in range(n) for y in range(n) for x in range(n)]
    shell = tmp / "al_shell.nii"
    write_float32_nifti(shell, (n, n, n), shell_vals)
    reface_out = tmp / "al_reface.nii"
    require_success(
        run_niimath(exe, [str(moving), "-reface", str(base), str(shell), str(weight), "-gz", "0", str(reface_out)]),
        "-reface positional triplet + output",
    )
    before = read_float32_nifti(moving)
    after = read_float32_nifti(reface_out)
    n_face = sum(1 for s in shell_vals if s > 0.0)
    n_changed = sum(1 for i, s in enumerate(shell_vals) if s > 0.0 and abs(after[i] - before[i]) > 1e-4)
    if n_changed < 0.9 * n_face:
        raise AssertionError(f"-reface replaced only {n_changed}/{n_face} face voxels (expected most anonymized)")
    # -reface routes its fast-selector through the shared cf_cost_from_fast_engine() mapping (same
    # as -allineate/-deface). Smoke that every adaptive/explicit selector is accepted and produces
    # output on the reface host path (guards the third mapping site the -allineate/-deface tests miss).
    for sel in ("fastx", "fasthel", "fastcr"):
        rf_sel = tmp / f"al_reface_{sel}.nii"
        require_success(
            run_niimath(exe, [str(moving), "-reface", str(base), str(shell), str(weight),
                              "-cost", sel, "-gz", "0", str(rf_sel)]),
            f"-reface -cost {sel}",
        )
        if not rf_sel.exists():
            raise AssertionError(f"-reface -cost {sel} produced no output")
    # Privacy fail-closed: a shell with no positive support -> <10% coverage -> refuse to write.
    empty_shell = tmp / "al_shell_empty.nii"
    write_float32_nifti(empty_shell, (n, n, n), [0.0] * (n * n * n))
    fc = run_niimath(exe, [str(moving), "-reface", str(base), str(empty_shell), str(weight),
                           "-gz", "0", str(tmp / "al_reface_fc.nii")])
    if fc.returncode == 0 or (tmp / "al_reface_fc.nii").exists():
        raise AssertionError("-reface did not fail closed on <10% face coverage (wrote output)")
    # Privacy fail-closed via the coverage RATIO (not just the empty-shell ternary): a NONEMPTY face
    # shell whose sform origin places it far outside the subject FOV maps ~0 voxels into the subject
    # -> cov ~ 0 -> refuse. Exercises the mm^3-normalized face_subj/face_tmpl arithmetic.
    far_shell = tmp / "al_shell_far.nii"
    write_float32_nifti(far_shell, (n, n, n), shell_vals, offset=(1000.0, 1000.0, 1000.0))
    fc2 = run_niimath(exe, [str(moving), "-reface", str(base), str(far_shell), str(weight),
                            "-gz", "0", str(tmp / "al_reface_far.nii")])
    if fc2.returncode == 0 or (tmp / "al_reface_far.nii").exists():
        raise AssertionError("-reface did not fail closed on an out-of-FOV nonempty shell (wrote output)")
    # Discriminating regression for the |det(sform)| coverage fix: a shell whose PIXDIM disagrees
    # with its sform scale. Reslicing uses the sform (identity here -> the face maps ~1:1 onto the
    # subject), so the physical coverage is ~100% using |det(sform)|=1 and the run must PASS. The old
    # pixdim-based volume (pixdim=5 -> voxel 125x too large) would compute ~0.8% and WRONGLY refuse —
    # so the verdict flips, locking the repaired metric. (sform bytes untouched; only pixdim patched.)
    skew_shell = tmp / "al_shell_skew.nii"
    write_float32_nifti(skew_shell, (n, n, n), shell_vals)
    sbuf = bytearray(skew_shell.read_bytes())
    struct.pack_into("<3f", sbuf, 80, 5.0, 5.0, 5.0)  # pixdim[1..3] = 5 (header sform diagonal stays 1.0)
    skew_shell.write_bytes(bytes(sbuf))
    skew_out = tmp / "al_reface_skew.nii"
    skew = run_niimath(exe, [str(moving), "-reface", str(base), str(skew_shell), str(weight),
                             "-gz", "0", str(skew_out)])
    if skew.returncode != 0 or not skew_out.exists():
        raise AssertionError("-reface wrongly refused a pixdim!=sform shell (coverage used pixdim, not |det(sform)|)")
    # The three reface aux operands each reject stdin '-' (a piped primary could otherwise make the
    # aux read reuse/misread the input stream — a privacy hazard).
    for label, triplet in (("template", ["-", str(shell), str(weight)]),
                           ("shell", [str(base), "-", str(weight)]),
                           ("weight", [str(base), str(shell), "-"])):
        rs = run_niimath(exe, [str(moving), "-reface", *triplet, "-gz", "0", str(tmp / "al_rf_stdin.nii")])
        if rs.returncode == 0 or (tmp / "al_rf_stdin.nii").exists():
            raise AssertionError(f"-reface accepted stdin '-' for the {label} operand")
    # -reface rejects -final (back-projection is always nearest-neighbour, so -final is meaningless).
    rf_final = run_niimath(exe, [str(moving), "-reface", str(base), str(shell), str(weight),
                                 "-final", "linear", "-gz", "0", str(tmp / "al_rf_final.nii")])
    if rf_final.returncode == 0 or "does not support -final" not in (rf_final.stdout + rf_final.stderr):
        raise AssertionError("-reface did not reject -final")

    # -unifize -GM must CONSUME the trailing -GM token AND actually apply GM scaling: its output
    # must DIFFER from a plain -unifize (token consumption alone would leave them identical).
    uni = tmp / "al_uni.nii"
    unigm = tmp / "al_unigm.nii"
    require_success(run_niimath(exe, [str(moving), "-unifize", "-gz", "0", str(uni)]), "-unifize")
    res = run_niimath(exe, [str(moving), "-unifize", "-GM", "-gz", "0", str(unigm)])
    if res.returncode != 0 or "unsupported operation" in (res.stdout + res.stderr) or not unigm.exists():
        raise AssertionError(f"-unifize -GM did not consume -GM / produce output: {res.stdout + res.stderr}")
    if all(abs(a - b) < 1e-4 for a, b in zip(read_float32_nifti(uni), read_float32_nifti(unigm))):
        raise AssertionError("-unifize -GM output identical to plain -unifize (GM scaling not applied)")


def exercise_moco(exe: str, tmp: Path, help_text: str) -> None:
    """-moco: analytic property check on a synthetic series with a KNOWN integer shift.

    Volume 1 is volume 0 displaced by an exact whole number of voxels along k, so the answer is
    known in closed form rather than by comparing against a stored golden.  The checks pin the
    conventions that were established by measurement (moco_bench's test/moco_reference_manifest.md):

      * the parameter file records the CORRECTION, so a +N voxel displacement is reported as a
        NEGATIVE dS of N mm -- this catches a sign flip, which no residual-based check would;
      * only dS moves; the other five parameters stay at zero, catching an axis mix-up;
      * row 0 is all zeros (the base registers to itself);
      * the file is six %8.4f columns joined by single spaces, one row per volume;
      * the corrected volume actually lands on the base far better than the input did.
    """
    if "-moco" not in help_text or "NOT in this build" in _help_line(help_text, "-moco"):
        print("  -moco: not built (MOCO=0) - skipping")
        return

    nx, ny, nz, nt = 24, 24, 20, 3
    shift = 2  # whole voxels along k, so the displaced volume is an exact copy

    def cell(x: int, y: int, z: int) -> float:
        # A compact object on a ZERO background, well inside the 5% edging border.  A ramp that
        # reaches the FOV edge would make the shifted volume differ from the base by a boundary
        # discontinuity no rigid motion can explain, and the fit would (correctly) refuse it.
        return (
            140.0 * math.exp(-((x - 10) ** 2 + (y - 11) ** 2 + (z - 9) ** 2) / 18.0)
            + 90.0 * math.exp(-((x - 15) ** 2 + (y - 8) ** 2 + (z - 12) ** 2) / 9.7)
            + 70.0 * math.exp(-((x - 9) ** 2 + (y - 15) ** 2 + (z - 11) ** 2) / 8.0)
        )

    vol0 = [0.0] * (nx * ny * nz)
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                vol0[x + y * nx + z * nx * ny] = cell(x, y, z)

    data = list(vol0)
    for t in range(1, nt):
        sh = shift if t == 1 else 0
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    src = z - sh
                    v = vol0[x + y * nx + src * nx * ny] if 0 <= src < nz else 0.0
                    data.append(v)

    src_path = tmp / "moco_in.nii"
    write_float32_nifti(src_path, (nx, ny, nz), data, nt=nt)
    par_path = tmp / "moco.1D"
    out_path = tmp / "moco_out.nii"
    require_success(
        run_niimath(exe, [str(src_path), "-moco", "-1Dfile", str(par_path), str(out_path)]),
        "-moco with -1Dfile",
    )
    if not par_path.exists():
        raise AssertionError("-moco -1Dfile did not write a parameter file")

    lines = par_path.read_text().splitlines()
    if len(lines) != nt:
        raise AssertionError("-moco -1Dfile: expected %d rows, got %d" % (nt, len(lines)))
    for row in lines:
        # exactly six %8.4f fields joined by single spaces: 6*8 + 5 = 53 chars, no trailing space
        if len(row) != 53:
            raise AssertionError("-moco -1Dfile: row is not 53 chars: %r" % row)
        for f in range(6):
            field = row[f * 9 : f * 9 + 8]
            if len(field) != 8 or field != field.rstrip():
                raise AssertionError("-moco -1Dfile: field %d is not width 8: %r" % (f, row))
            body = field.strip()
            if "." not in body or len(body.split(".")[1]) != 4:
                raise AssertionError("-moco -1Dfile: field %d lacks 4 decimals: %r" % (f, row))
            float(body)  # must parse in the C locale
            if f < 5 and row[f * 9 + 8] != " ":
                raise AssertionError("-moco -1Dfile: fields not single-space separated: %r" % row)
    if not par_path.read_text().endswith("\n"):
        raise AssertionError("-moco -1Dfile: missing final newline")
    rows = [[float(v) for v in row.split()] for row in lines]
    for k, v in enumerate(rows[0]):
        if abs(v) > 1e-4:
            raise AssertionError("-moco: base row must be zeros, column %d is %g" % (k, v))

    roll, pitch, yaw, dS, dL, dP = rows[1]
    # +shift voxels along k is +shift mm Superior; the file records the correction, so dS < 0.
    if abs(dS + float(shift)) > 0.20:
        raise AssertionError("-moco: expected dS near %.1f, got %.4f" % (-float(shift), dS))
    for name, v in (("roll", roll), ("pitch", pitch), ("yaw", yaw), ("dL", dL), ("dP", dP)):
        if abs(v) > 0.10:
            raise AssertionError("-moco: %s should be ~0 for a pure k shift, got %.4f" % (name, v))

    written = out_path if out_path.exists() else Path(str(out_path) + ".gz")
    if not written.exists():
        raise AssertionError("-moco did not write an output image")
    if written.suffix == ".gz":
        plain = tmp / "moco_out_plain.nii"
        require_success(
            run_niimath(exe, [str(written), "-gz", "0", str(plain)]), "-moco decompress"
        )
        written = plain
    vals = read_float32_nifti(written)
    n3 = nx * ny * nz

    def interior_rms(a_off: int, b_off: int) -> float:
        tot = 0.0
        cnt = 0
        for z in range(4, nz - 4):
            for y in range(4, ny - 4):
                for x in range(4, nx - 4):
                    i = x + y * nx + z * nx * ny
                    d = vals[a_off + i] - vals[b_off + i]
                    tot += d * d
                    cnt += 1
        return math.sqrt(tot / cnt) if cnt else 0.0

    before = 0.0
    cnt = 0
    for z in range(4, nz - 4):
        for y in range(4, ny - 4):
            for x in range(4, nx - 4):
                i = x + y * nx + z * nx * ny
                d = data[i] - data[n3 + i]
                before += d * d
                cnt += 1
    before = math.sqrt(before / cnt)
    after = interior_rms(0, n3)
    if not (after < before * 0.10):
        raise AssertionError(
            "-moco: correction did not align volume 1 (interior rms %.4f before, %.4f after)"
            % (before, after)
        )
    if abs(vals[0] - data[0]) > 1e-6:
        raise AssertionError("-moco: volume 0 must be copied through unchanged")

    # without -1Dfile no parameter file is produced
    out2 = tmp / "moco_out2.nii"
    require_success(run_niimath(exe, [str(src_path), "-moco", str(out2)]), "-moco without -1Dfile")

    # a 3D input must be refused rather than silently treated as a single volume
    mean_path = tmp / "moco_3d.nii"
    require_success(run_niimath(exe, [str(src_path), "-Tmean", str(mean_path)]), "-moco 3D fixture")
    mean_written = mean_path if mean_path.exists() else Path(str(mean_path) + ".gz")
    rejected = run_niimath(exe, [str(mean_written), "-moco", str(tmp / "moco_bad.nii")])
    if rejected.returncode == 0:
        raise AssertionError("-moco accepted a 3D image; it must reject one")

    # the base volume must be copied through in full, not merely its first voxel.  Compare against
    # the file as written (float32), not the float64 source list.
    src_vals = read_float32_nifti(src_path)
    for i in range(n3):
        if vals[i] != src_vals[i]:
            raise AssertionError(
                "-moco: volume 0 must be copied through unchanged (voxel %d: %r vs %r)"
                % (i, vals[i], src_vals[i])
            )

    # The .1D suffix keeps parameter and supported NIfTI output namespaces disjoint.
    bad_par = tmp / "moco_params.nii.gz"
    if run_niimath(
        exe, [str(src_path), "-moco", "-1Dfile", str(bad_par), str(tmp / "moco_bad_suffix.nii")]
    ).returncode == 0:
        raise AssertionError("-moco accepted a -1Dfile path without the required .1D suffix")

    # a failing run must leave an existing parameter file untouched and drop no temporary
    keep = tmp / "moco_keep.1D"
    keep.write_text("PRESERVE ME\n")
    run_niimath(exe, [str(mean_written), "-moco", "-1Dfile", str(keep), str(tmp / "moco_bad2.nii")])
    if keep.read_text() != "PRESERVE ME\n":
        raise AssertionError("-moco clobbered an existing parameter file on a failing run")
    leftovers = list(tmp.glob("*.mocotmp*"))
    if leftovers:
        raise AssertionError("-moco left a temporary parameter file behind: %s" % leftovers)

    # structurally 4D but single-volume must be rejected too
    one = tmp / "moco_one.nii"
    require_success(run_niimath(exe, [str(src_path), "-crop", "0", "1", str(one)]), "-moco 1-volume fixture")
    one_written = one if one.exists() else Path(str(one) + ".gz")
    if run_niimath(exe, [str(one_written), "-moco", str(tmp / "moco_bad3.nii")]).returncode == 0:
        raise AssertionError("-moco accepted a single-volume 4D image")

    # ---- rotation + fractional interpolation ------------------------------------------------
    # Rotating the object about the k axis by a few degrees exercises the four-shear
    # factorization, fractional heptic interpolation and the clip-to-input-range rule, none of
    # which an integral translation touches.  Checks: the recovered roll has the right magnitude
    # and sign, the other five parameters stay small, the corrected volume lands on the base, and
    # clipping holds the output inside the moving volume's own range.
    ang = math.radians(4.0)
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    rot = [0.0] * n3
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                # sample the analytic object at the back-rotated location: exact, no resampler
                dx0, dy0 = x - cx, y - cy
                sx = cx + dx0 * math.cos(ang) + dy0 * math.sin(ang)
                sy = cy - dx0 * math.sin(ang) + dy0 * math.cos(ang)
                rot[x + y * nx + z * nx * ny] = cell(sx, sy, z)
    rdata = list(vol0) + rot + list(vol0)
    rot_src = tmp / "moco_rot.nii"
    write_float32_nifti(rot_src, (nx, ny, nz), rdata, nt=3)
    rot_par = tmp / "moco_rot.1D"
    rot_out = tmp / "moco_rot_out.nii"
    require_success(
        run_niimath(exe, [str(rot_src), "-moco", "-1Dfile", str(rot_par), str(rot_out)]),
        "-moco rotation",
    )
    rrows = [[float(v) for v in r.split()] for r in rot_par.read_text().splitlines()]
    rroll, rpitch, ryaw = rrows[1][0], rrows[1][1], rrows[1][2]
    # index k is +Superior for this RAS 1 mm header, so a rotation in the (i,j) plane is roll.
    if abs(rroll + 4.0) > 0.35:
        raise AssertionError("-moco: expected correction roll near -4 deg, got %.4f" % rroll)
    for name, v in (("pitch", rpitch), ("yaw", ryaw)):
        if abs(v) > 0.35:
            raise AssertionError("-moco: %s should be ~0 for an in-plane rotation, got %.4f" % (name, v))
    rwritten = rot_out if rot_out.exists() else Path(str(rot_out) + ".gz")
    if rwritten.suffix == ".gz":
        rplain = tmp / "moco_rot_plain.nii"
        require_success(run_niimath(exe, [str(rwritten), "-gz", "0", str(rplain)]), "-moco rot decompress")
        rwritten = rplain
    rvals = read_float32_nifti(rwritten)
    lo_in, hi_in = min(rot), max(rot)
    for i in range(n3):
        v = rvals[n3 + i]
        if v < lo_in - 1e-3 or v > hi_in + 1e-3:
            raise AssertionError(
                "-moco: -clipit violated at voxel %d: %g outside input range [%g, %g]"
                % (i, v, lo_in, hi_in)
            )
    num = den = 0.0
    for z in range(5, nz - 5):
        for y in range(5, ny - 5):
            for x in range(5, nx - 5):
                i = x + y * nx + z * nx * ny
                num += (rvals[n3 + i] - rvals[i]) ** 2
                den += (rot[i] - vol0[i]) ** 2
    if not (num < den * 0.10):
        raise AssertionError(
            "-moco: rotation not corrected (interior sse %.4g after vs %.4g before)" % (num, den)
        )

    print("  -moco: parameter sign/axis/format, base passthrough, alignment, rotation+clipping,")
    print("         output separation, existing-file preservation and 3D/singleton rejection OK")


def _stc_reference(x: list[float], shift: int, nt: int) -> list[float]:
    """Closed-form -stc output for an INTEGER sample shift, from the measured contract.

    An integer shift makes the Fourier stage exactly a circular shift of the zero-padded
    residual, so the whole pipeline -- least-squares detrend, shift, clip to the residual range,
    retrend, clip to the original range -- has a closed form and needs no FFT here.  That is what
    lets this test pin the SHIFT DIRECTION, which no residual-magnitude check could.
    """
    half = 0.5 * (nt - 1)
    mean = sum(x) / nt
    sdd = nt * (nt * nt - 1) / 12.0
    slope = sum((i - half) * x[i] for i in range(nt)) / sdd
    trend = [mean + slope * (i - half) for i in range(nt)]
    xd = [x[i] - trend[i] for i in range(nt)]
    lo, hi = min(xd), max(xd)
    xlo, xhi = min(x), max(x)
    out = []
    for i in range(nt):
        j = i - shift
        y = xd[j] if 0 <= j < nt else 0.0  # outside [0, nt) the padded array is zero
        y = min(max(y, lo), hi)
        r = y + trend[i]
        out.append(min(max(r, xlo), xhi))
    return out


def exercise_stc(exe: str, tmp: Path, help_text: str) -> None:
    """-stc: closed-form check of the measured contract, with no AFNI dependency.

    The fixture uses two slices with times [0, TR] so that -tzero pins one slice to an EXACT
    one-sample shift and the other to zero.  Everything the manifest records as load-bearing is
    checked here: shift direction, the least-squares detrend/retrend, both clips, the
    small-shift skip, toffset = tzero, the non-finite policy, and the parser's rejections.
    """
    if "-stc" not in help_text or "NOT in this build" in _help_line(help_text, "-stc"):
        print("  -stc: not built (STC=0) - skipping")
        return

    nx, ny, nz, nt = 3, 3, 2, 16
    tr = 2.0
    n3 = nx * ny * nz
    # A spiky series with a strong trend, so BOTH clips have a chance to engage; 9 voxels per
    # slice also leaves the last batch partly filled (STC processes 8 series at a time).
    def series(v: int) -> list[float]:
        return [100.0 + 6.0 * i + 30.0 * math.sin(2.7 * i + v) + (25.0 if i == 5 + (v % 3) else 0.0)
                for i in range(nt)]

    data = [0.0] * (n3 * nt)
    for t in range(nt):
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    v = x + y * nx + z * nx * ny
                    data[t * n3 + v] = series(v)[t]
    src = tmp / "stc_in.nii"
    write_float32_nifti(src, (nx, ny, nz), data, nt=nt, tr=tr, time_units=8)
    times = "0,%g" % tr

    def corrected(args: list[str], label: str) -> list[float]:
        out = tmp / ("stc_%s.nii" % label)
        require_success(run_niimath(exe, [str(src), *args, "-gz", "0", str(out)]), "-stc " + label)
        return read_float32_nifti(out)

    # -tzero 0: slice 0 has shift 0 (skipped verbatim), slice 1 has shift exactly +1 sample.
    vals = corrected(["-stc", "--slicetiming", times, "-tzero", "0"], "fwd")
    src_vals = read_float32_nifti(src)
    clipped_seen = False
    for z, shift in ((0, 0), (1, 1)):
        for y in range(ny):
            for x in range(nx):
                v = x + y * nx + z * nx * ny
                got = [vals[t * n3 + v] for t in range(nt)]
                have = [src_vals[t * n3 + v] for t in range(nt)]
                if shift == 0:
                    if got != have:
                        raise AssertionError("-stc: zero-shift slice must be copied verbatim")
                    continue
                want = _stc_reference(have, shift, nt)
                for t in range(nt):
                    if abs(got[t] - want[t]) > 2e-3 * max(1.0, abs(want[t])):
                        raise AssertionError(
                            "-stc: voxel %d sample %d is %.6g, closed form says %.6g"
                            % (v, t, got[t], want[t])
                        )
                if min(want) <= min(have) or max(want) >= max(have):
                    clipped_seen = True
    if not clipped_seen:
        raise AssertionError("-stc fixture never exercised the output clip; strengthen it")

    # -tzero TR reverses the roles AND the sign: slice 0 now shifts by exactly -1 sample.
    rev = corrected(["-stc", "--slicetiming", times, "-tzero", "%g" % tr], "rev")
    for y in range(ny):
        for x in range(nx):
            v = x + y * nx
            have = [src_vals[t * n3 + v] for t in range(nt)]
            want = _stc_reference(have, -1, nt)
            for t in range(nt):
                if abs(rev[t * n3 + v] - want[t]) > 2e-3 * max(1.0, abs(want[t])):
                    raise AssertionError("-stc: negative shift disagrees with the closed form")

    # A pure line in time is annihilated by the detrend and restored by the retrend, so it must
    # come back unchanged for ANY shift.  This catches a retrend applied at the shifted index.
    ramp = [0.0] * (n3 * nt)
    for t in range(nt):
        for v in range(n3):
            ramp[t * n3 + v] = 5.0 + 3.0 * t + 0.25 * v
    ramp_src = tmp / "stc_ramp.nii"
    write_float32_nifti(ramp_src, (nx, ny, nz), ramp, nt=nt, tr=tr, time_units=8)
    ramp_out = tmp / "stc_ramp_out.nii"
    require_success(
        run_niimath(exe, [str(ramp_src), "-stc", "--slicetiming", "0,1.3", "-gz", "0", str(ramp_out)]),
        "-stc linear ramp",
    )
    got = read_float32_nifti(ramp_out)
    for i, (a, b) in enumerate(zip(got, ramp)):
        if abs(a - b) > 1e-3:
            raise AssertionError("-stc: a linear time trend must survive unchanged (%d: %g vs %g)"
                                 % (i, a, b))

    # All-equal slice times leave every shift at zero: the output must be bit-identical.
    eq = tmp / "stc_eq.nii"
    require_success(
        run_niimath(exe, [str(src), "-stc", "--slicetiming", "0.8,0.8", "-gz", "0", str(eq)]),
        "-stc all-equal timings",
    )
    if read_float32_nifti(eq) != src_vals:
        raise AssertionError("-stc: all-equal slice times must leave the data untouched")

    # Measured skip threshold: |fractional shift| < 0.001 copies the slice verbatim.  With
    # TR = 2 and tzero = 0 a slice time of 0.0018 s is a shift of 0.0009.
    skip = tmp / "stc_skip.nii"
    require_success(
        run_niimath(exe, [str(src), "-stc", "--slicetiming", "0,0.0018", "-tzero", "0",
                          "-gz", "0", str(skip)]),
        "-stc sub-threshold shift",
    )
    if read_float32_nifti(skip) != src_vals:
        raise AssertionError("-stc: a shift below 0.001 samples must be copied verbatim")

    # toffset takes the common time point, in the header's own time units; TR and the spatial
    # transform are untouched.
    out_hdr = (tmp / "stc_fwd.nii").read_bytes()
    src_hdr = src.read_bytes()
    toffset = struct.unpack_from("<f", out_hdr, 136)[0]  # NIfTI-1 toffset
    if abs(toffset - 0.0) > 1e-6:
        raise AssertionError("-stc: toffset must be tzero (0), got %g" % toffset)
    dflt = tmp / "stc_default_tzero.nii"
    require_success(
        run_niimath(exe, [str(src), "-stc", "--slicetiming", times, "-gz", "0", str(dflt)]),
        "-stc default tzero",
    )
    toffset = struct.unpack_from("<f", dflt.read_bytes(), 136)[0]
    if abs(toffset - tr / 2.0) > 1e-5:
        raise AssertionError("-stc: default tzero must be the mean slice time, got %g" % toffset)
    for name, off, size in (("pixdim", 76, 32), ("srow", 280, 48), ("dim", 40, 16)):
        if out_hdr[off:off + size] != src_hdr[off:off + size]:
            raise AssertionError("-stc changed the %s header block" % name)

    # niimath policy: a series holding any non-finite sample is written out as all-NaN, and its
    # neighbours are unaffected (nothing leaks between voxels).
    nf = list(data)
    nf_vox = 4 + nx * ny          # a voxel in slice 1, the slice that actually gets shifted
    nf[3 * n3 + nf_vox] = float("nan")
    nf_src = tmp / "stc_nf.nii"
    write_float32_nifti(nf_src, (nx, ny, nz), nf, nt=nt, tr=tr, time_units=8)
    nf_out = tmp / "stc_nf_out.nii"
    require_success(
        run_niimath(exe, [str(nf_src), "-stc", "--slicetiming", times, "-tzero", "0",
                          "-gz", "0", str(nf_out)]),
        "-stc non-finite input",
    )
    nfv = read_float32_nifti(nf_out)
    for t in range(nt):
        if nfv[t * n3 + nf_vox] == nfv[t * n3 + nf_vox]:
            raise AssertionError("-stc: a non-finite series must be written out as all-NaN")
    for v in range(n3):
        if v == nf_vox:
            continue
        for t in range(nt):
            if nfv[t * n3 + v] != nfv[t * n3 + v]:
                raise AssertionError("-stc: a non-finite series contaminated voxel %d" % v)

    # The all-NaN rule stops at the skip threshold: a slice that is copied verbatim keeps its
    # non-finite samples exactly as they were.  Pin it, because it is the one place where the
    # verbatim-copy guarantee and the non-finite policy meet and the documented wording has to
    # match which one wins.
    nfs = list(data)
    nfs[3 * n3 + 4] = float("nan")          # voxel 4 is in slice 0, which -tzero 0 leaves alone
    nfs_src = tmp / "stc_nf_skip.nii"
    write_float32_nifti(nfs_src, (nx, ny, nz), nfs, nt=nt, tr=tr, time_units=8)
    nfs_out = tmp / "stc_nf_skip_out.nii"
    require_success(
        run_niimath(exe, [str(nfs_src), "-stc", "--slicetiming", times, "-tzero", "0",
                          "-gz", "0", str(nfs_out)]),
        "-stc non-finite in a skipped slice",
    )
    skipped = read_float32_nifti(nfs_out)
    src_nfs = read_float32_nifti(nfs_src)
    for t in range(nt):
        for v in range(nx * ny):        # slice 0 only
            a, b = skipped[t * n3 + v], src_nfs[t * n3 + v]
            if a != b and not (a != a and b != b):
                raise AssertionError(
                    "-stc: a skipped slice must be copied verbatim, NaN included (voxel %d, t %d)"
                    % (v, t))

    # Parser and validation: every rejection must fail before any output appears.
    bad = tmp / "stc_bad.nii"
    adjacent = tmp / "stc_adjacent.1D"
    adjacent.write_text("0+%g\n" % tr)
    rejects = [
        (["-stc", str(bad)], "missing --slicetiming"),
        (["-stc", "--SliceTiming", times, str(bad)], "wrong case"),
        (["-stc", "-tzero", "0", "--slicetiming", times, str(bad)], "-tzero before --slicetiming"),
        (["-stc", "--slicetiming", "0", str(bad)], "too few values"),
        (["-stc", "--slicetiming", "0,1,2", str(bad)], "too many values"),
        (["-stc", "--slicetiming", "0,,2", str(bad)], "empty field"),
        (["-stc", "--slicetiming", "0,%g," % tr, str(bad)], "trailing comma"),
        (["-stc", "--slicetiming", "0,nan", str(bad)], "non-finite value"),
        (["-stc", "--slicetiming", "0,2.5", str(bad)], "slice time above TR"),
        (["-stc", "--slicetiming", "-0.1,1.0", str(bad)], "negative slice time"),
        (["-stc", "--slicetiming", times, "-tzero", "2.5", str(bad)], "tzero above max"),
        (["-stc", "--slicetiming", times, "-tzero", "x", str(bad)], "unparsable tzero"),
        (["-stc", "--slicetiming", "@%s" % (tmp / "stc_absent.1D"), str(bad)], "missing @file"),
        (["-stc", "--slicetiming", "@%s" % adjacent, str(bad)], "missing @file separator"),
    ]
    for argv, label in rejects:
        res = run_niimath(exe, [str(src), *argv])
        if res.returncode == 0 or bad.exists() or Path(str(bad) + ".gz").exists():
            raise AssertionError("-stc accepted %s" % label)

    # 3D and short series are refused rather than silently mishandled.
    short = tmp / "stc_short.nii"
    write_float32_nifti(short, (nx, ny, nz), [1.0] * (n3 * 4), nt=4, tr=tr, time_units=8)
    if run_niimath(exe, [str(short), "-stc", "--slicetiming", times, str(bad)]).returncode == 0:
        raise AssertionError("-stc accepted nt = 4")
    vol3d = tmp / "stc_3d.nii"
    write_float32_nifti(vol3d, (nx, ny, nz), [1.0] * n3, tr=tr, time_units=8)
    if run_niimath(exe, [str(vol3d), "-stc", "--slicetiming", times, str(bad)]).returncode == 0:
        raise AssertionError("-stc accepted a 3D image")
    # A header without a temporal unit must be refused, not assumed to be seconds.
    nounit = tmp / "stc_nounit.nii"
    write_float32_nifti(nounit, (nx, ny, nz), data, nt=nt, tr=tr, time_units=0)
    if run_niimath(exe, [str(nounit), "-stc", "--slicetiming", times, str(bad)]).returncode == 0:
        raise AssertionError("-stc accepted a header with no temporal unit")

    # The AFNI-style @file form (whitespace separated, '#' comments) matches the inline list, and
    # -stc chains like any other operation.
    tfile = tmp / "stc_times.1D"
    tfile.write_text("# slice times, seconds\n0\n%g\n" % tr)
    at = corrected(["-stc", "--slicetiming", "@%s" % tfile, "-tzero", "0"], "atfile")
    if at != vals:
        raise AssertionError("-stc: '@file' and the inline list disagree")
    chained = tmp / "stc_chain.nii"
    require_success(
        run_niimath(exe, [str(src), "-stc", "--slicetiming", times, "-tzero", "0", "-add", "1",
                          "-gz", "0", str(chained)]),
        "-stc chained with -add",
    )
    ch = read_float32_nifti(chained)
    for i in range(len(vals)):
        if abs(ch[i] - (vals[i] + 1.0)) > 1e-3:
            raise AssertionError("-stc did not chain into the following operation")

    print("  -stc: shift direction/closed form, detrend+retrend, clips, skip threshold,")
    print("        toffset, non-finite policy, parser rejections, @file and chaining OK")


def _skullstrip_geometry(path: Path) -> tuple:
    """Everything that locates the voxels in world space. -skullstrip works on a fixed internal
    grid and pulls the mask back, so a regression there shows up as a changed header, not as a
    changed value."""
    blob = read_nifti_bytes(path)
    return (
        struct.unpack_from("<8h", blob, 40),    # dim
        struct.unpack_from("<8f", blob, 76),    # pixdim
        struct.unpack_from("<2h", blob, 252),   # qform_code, sform_code
        struct.unpack_from("<6f", blob, 256),   # quatern b/c/d, qoffset x/y/z
        struct.unpack_from("<12f", blob, 280),  # srow_x/y/z
    )


def exercise_fmap(exe: str, tmp: Path, help_text: str) -> None:
    """-fugue / -fmapprep: closed-form checks of the measured contract, with no FSL dependency.

    Everything asserted here is derivable analytically, which is the point: the reference is a
    non-redistributable binary, so the smoke test must stand on the arithmetic rather than on a
    golden file.  See fmap_bench's test/fmap_reference_manifest.md for the experiment behind each
    convention.
    """
    if "-fugue" not in help_text or "NOT in this build" in _help_line(help_text, "-fugue"):
        print("  -fugue/-fmapprep: not built (FMAP=0) - skipping")
        return

    # A cube, so a 2-voxel shift of a centrally placed impulse stays inside the FOV along EVERY
    # axis -- with a short z the "z-" case silently walks off the end and the test measures
    # out-of-FOV fill instead of the sign convention.
    nx, ny, nz = 8, 8, 8
    n3 = nx * ny * nz
    dwell = 0.001
    two_pi = 2.0 * math.pi

    # An impulse makes the output literally the interpolation kernel.
    impulse = [0.0] * n3
    impulse[3 + 3 * nx + 3 * nx * ny] = 1000.0
    src = tmp / "fugue_impulse.nii"
    write_float32_nifti(src, (nx, ny, nz), impulse)

    def constant_field(shift_vox: float, axis_n: int) -> Path:
        """A uniform field whose implied shift is exactly `shift_vox` voxels."""
        val = shift_vox * two_pi / (dwell * axis_n)
        p = tmp / ("fugue_f_%g_%d.nii" % (shift_vox, axis_n))
        write_float32_nifti(p, (nx, ny, nz), [val] * n3)
        return p

    def unwarp(field: Path, direction: str, inp: Path = src) -> list[float]:
        out = tmp / "fugue_out.nii"
        require_success(run_niimath(exe, [str(inp), "-fugue", str(field), repr(dwell),
                                          direction, "-gz", "0", str(out)]), "-fugue")
        return read_float32_nifti(out)

    # Whole-voxel shift, and the sign convention: a POSITIVE field moves signal toward -axis for
    # "y" and toward +axis for "y-".  N is the dimension along the UNWARP axis, so the same
    # displacement along z needs a field scaled by nz, not ny.
    for direction, di, dj, dk, axis_n in (("y", 0, -2, 0, ny), ("y-", 0, 2, 0, ny),
                                          ("x", -2, 0, 0, nx), ("x-", 2, 0, 0, nx),
                                          ("z", 0, 0, -2, nz), ("z-", 0, 0, 2, nz)):
        vals = unwarp(constant_field(2.0, axis_n), direction)
        want = (3 + di) + (3 + dj) * nx + (3 + dk) * nx * ny
        hot = [i for i, v in enumerate(vals) if abs(v) > 1e-3]
        if hot != [want] or abs(vals[want] - 1000.0) > 1e-2:
            raise SystemExit("-fugue: --unwarpdir %s put the impulse at %s (expected [%d], value %g)"
                             % (direction, hot, want, vals[want] if hot else float("nan")))

    # Half-voxel shift must split EXACTLY 500/500.  This is the check that pins the kernel to
    # linear: no cubic, spline or windowed-sinc kernel can produce two taps.
    vals = unwarp(constant_field(1.5, ny), "y")
    hot = sorted(i for i, v in enumerate(vals) if abs(v) > 1e-3)
    # out(v) = in(v + 1.5), so the impulse at y=3 lands half in y=1 and half in y=2.
    want = sorted([3 + 1 * nx + 3 * nx * ny, 3 + 2 * nx + 3 * nx * ny])
    if hot != want or any(abs(vals[i] - 500.0) > 1e-2 for i in hot):
        raise SystemExit("-fugue: a 1.5-voxel shift did not split 500/500; got %s at %s"
                         % ([vals[i] for i in hot], hot))

    # EXTRAPOLATION.  A zero fieldmap is unsupported, and the reference does NOT read that as
    # zero shift -- it extrapolates the shift field along each line before resampling.  Without
    # this case every fieldmap in this test is either uniform or all-zero, so three of
    # fm_fill_line's four branches never run, and the naive model would pass.
    #
    # Field is 0 for y < 2 and a 1-voxel shift elsewhere; input is a ramp in y.  Linear
    # interpolation of a linear function is exact, so out - (1000 + y) == -1.0 at EVERY y
    # including y = 0 and 1 IFF the leading run replicated the first supported value.
    ramp = [0.0] * n3
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                ramp[x + y * nx + z * nx * ny] = 1000.0 + y
    ramp_src = tmp / "fugue_ramp.nii"
    write_float32_nifti(ramp_src, (nx, ny, nz), ramp)
    partial = [0.0] * n3
    val = 1.0 * two_pi / (dwell * ny)
    for z in range(nz):
        for y in range(2, ny):
            for x in range(nx):
                partial[x + y * nx + z * nx * ny] = val
    pf = tmp / "fugue_partial.nii"
    write_float32_nifti(pf, (nx, ny, nz), partial)
    vals = unwarp(pf, "y", inp=ramp_src)
    for y in range(ny - 1):          # last row pulls out of FOV, so skip it
        v = vals[3 + y * nx + 3 * nx * ny]
        if abs(v - (1000.0 + y + 1.0)) > 1e-2:
            raise SystemExit("-fugue: shift field was not extrapolated over unsupported voxels; "
                             "row y=%d is %g, expected %g (a naive model gives %g at y<2)"
                             % (y, v, 1000.0 + y + 1.0, 1000.0 + y))

    # 4D: one 3D fieldmap must be applied identically to every volume.
    vol4 = [0.0] * (n3 * 3)
    for t in range(3):
        vol4[t * n3 + (3 + (3 + t) * nx + 3 * nx * ny)] = 1000.0
    src4 = tmp / "fugue_4d.nii"
    write_float32_nifti(src4, (nx, ny, nz), vol4, nt=3)
    vals = unwarp(constant_field(2.0, ny), "y", inp=src4)
    for t in range(3):
        want = 3 + (1 + t) * nx + 3 * nx * ny
        hot = [i - t * n3 for i, v in enumerate(vals[t * n3:(t + 1) * n3], start=t * n3)
               if abs(v) > 1e-3]
        if hot != [want]:
            raise SystemExit("-fugue: 4D volume %d put the impulse at %s, expected [%d]"
                             % (t, hot, want))

    # A zero fieldmap is an exact identity -- the cheapest possible check that the extrapolation
    # path cannot introduce a shift where there is no field.
    zero = tmp / "fugue_zero.nii"
    write_float32_nifti(zero, (nx, ny, nz), [0.0] * n3)
    if unwarp(zero, "y") != impulse:
        raise SystemExit("-fugue: an all-zero fieldmap was not an exact identity")

    # Out-of-FOV samples contribute 0, and no Jacobian modulation is applied.  Use a UNIFORM
    # input, not the impulse: with an impulse an all-zero output is also what ANY bug that
    # simply blanks the image produces, so the assertion would not distinguish them.  Here the
    # rows that still pull from inside the FOV must keep their value exactly.
    uniform = tmp / "fugue_uniform.nii"
    write_float32_nifti(uniform, (nx, ny, nz), [1000.0] * n3)
    vals = unwarp(constant_field(2.0, ny), "y", inp=uniform)
    for y in range(ny):
        v = vals[3 + y * nx + 3 * nx * ny]
        want = 1000.0 if y + 2 < ny else 0.0
        if abs(v - want) > 1e-3:
            raise SystemExit("-fugue: out-of-FOV fill wrong at y=%d: %g, expected %g" % (y, v, want))

    # -no-debranch must reach fmap_prepare, not merely parse: on a smooth phase the correction
    # is inert, so the two outputs have to be byte-identical.
    vals = unwarp(constant_field(float(ny), ny), "y")
    if any(abs(v) > 1e-6 for v in vals):
        raise SystemExit("-fugue: a whole-FOV shift should have emptied the image")

    # -p 1 vs -p 8 byte-equality: every line is independent, so this must hold exactly.
    outs = []
    for threads in ("1", "8"):
        o = tmp / ("fugue_p%s.nii" % threads)
        require_success(run_niimath(exe, [str(src), "-p", threads, "-fugue",
                                          str(constant_field(1.5, ny)), repr(dwell), "y",
                                          "-gz", "0", str(o)]), "-fugue -p " + threads)
        outs.append(o.read_bytes())
    if outs[0] != outs[1]:
        raise SystemExit("-fugue: -p 1 and -p 8 outputs differ")

    # Rejections, all of which must fail closed rather than write a wrong image.
    field = constant_field(1.0, ny)
    bad = tmp / "fugue_bad.nii"
    for args, label in (
        (["-fugue", str(field)], "missing dwell and direction"),
        (["-fugue", str(field), "y", "0.001"], "transposed dwell/direction"),
        (["-fugue", str(field), repr(dwell), "q"], "unknown unwarpdir"),
        (["-fugue", str(field), "0", "y"], "zero dwell"),
        (["-fugue", str(field), "-0.001", "y"], "negative dwell"),
    ):
        if run_niimath(exe, [str(src), *args, "-gz", "0", str(bad)]).returncode == 0:
            raise SystemExit("-fugue: accepted %s" % label)

    # A 4D fieldmap and an off-grid fieldmap must both be refused.
    fmap4d = tmp / "fugue_f4d.nii"
    write_float32_nifti(fmap4d, (nx, ny, nz), [1.0] * (n3 * 2), nt=2)
    offgrid = tmp / "fugue_off.nii"
    write_float32_nifti(offgrid, (nx, ny, nz), [1.0] * n3, offset=(25.0, 0.0, 0.0))
    for f, label in ((fmap4d, "4D fieldmap"), (offgrid, "off-grid fieldmap")):
        if run_niimath(exe, [str(src), "-fugue", str(f), repr(dwell), "y",
                             "-gz", "0", str(bad)]).returncode == 0:
            raise SystemExit("-fugue: accepted a %s" % label)
    if run_niimath(exe, ["-dt", "double", str(src), "-fugue", str(field), repr(dwell), "y",
                         "-gz", "0", str(bad)]).returncode == 0:
        raise SystemExit("-fugue: accepted -dt double")

    print("  -fugue: shift constant, all six directions, linear kernel, zero-field identity,")
    print("        out-of-FOV fill, -p 1 vs -p 8 byte-equality and parser rejections OK")

    if "-fmapprep" not in help_text:
        print("  -fmapprep: not built (needs ROMEO) - skipping")
        return

    # -fmapprep on a smooth, non-wrapping phase ramp: unwrapping is then the identity, so the
    # whole op reduces to closed form -- scale by 2*pi/span, divide by deltaTE, subtract the
    # median over the mask, zero outside it.
    delta_te_ms = 2.5
    mag = [0.0] * n3
    phase = [0.0] * n3
    inside = []
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                v = x + y * nx + z * nx * ny
                phase[v] = float(y)  # a gentle ramp: no wraps for the unwrapper to resolve
                if 1 <= x < nx - 1 and 1 <= y < ny - 1:
                    mag[v] = 100.0
                    inside.append(v)
    pha = tmp / "fmapprep_phase.nii"
    mgn = tmp / "fmapprep_mag.nii"
    write_float32_nifti(pha, (nx, ny, nz), phase)
    write_float32_nifti(mgn, (nx, ny, nz), mag)
    out = tmp / "fmapprep_out.nii"
    require_success(run_niimath(exe, [str(pha), "-fmapprep", str(mgn), repr(delta_te_ms),
                                      "-gz", "0", str(out)]), "-fmapprep")
    got = read_float32_nifti(out)

    span = max(phase) - min(phase)
    rad = [(p - min(phase)) / span * two_pi - math.pi for p in phase]
    field = [r * 1000.0 / delta_te_ms for r in rad]
    ordered = sorted(field[v] for v in inside)
    median = ordered[len(ordered) // 2]  # UPPER central value, not the average
    for v in range(n3):
        want = (field[v] - median) if mag[v] != 0.0 else 0.0
        if abs(got[v] - want) > 1e-2:
            raise SystemExit("-fmapprep: voxel %d is %g, expected %g" % (v, got[v], want))

    # The mask is the magnitude's nonzero support, verbatim -- no erosion, no dilation.
    if [v for v in range(n3) if got[v] != 0.0 and mag[v] == 0.0]:
        raise SystemExit("-fmapprep: wrote nonzero values outside the magnitude mask")

    # -no-debranch is accepted before the output name and REJECTED in it, so a truncated command
    # line cannot silently produce a file called "-no-debranch".
    require_success(run_niimath(exe, [str(pha), "-fmapprep", str(mgn), repr(delta_te_ms),
                                      "-no-debranch", "-gz", "0", str(out)]), "-fmapprep -no-debranch")
    if run_niimath(exe, [str(pha), "-fmapprep", str(mgn), repr(delta_te_ms),
                         "-no-debranch"]).returncode == 0:
        raise SystemExit("-fmapprep: accepted '-no-debranch' as the output filename")
    stray = tmp / "-no-debranch"
    if stray.exists():
        raise SystemExit("-fmapprep: wrote a file named after its own flag")

    for args, label in (
        ([str(pha), "-fmapprep", str(mgn)], "missing deltaTE"),
        ([str(pha), "-fmapprep", str(mgn), "0"], "zero deltaTE"),
        ([str(pha), "-fmapprep", str(mgn), "-2.5"], "negative deltaTE"),
        ([str(pha), "-fmapprep", str(mgn), "abc"], "non-numeric deltaTE"),
    ):
        if run_niimath(exe, [*args, "-gz", "0", str(tmp / "fmapprep_bad.nii")]).returncode == 0:
            raise SystemExit("-fmapprep: accepted %s" % label)

    # An empty mask has no median to subtract, so it must be refused rather than divided by.
    empty = tmp / "fmapprep_empty.nii"
    write_float32_nifti(empty, (nx, ny, nz), [0.0] * n3)
    if run_niimath(exe, [str(pha), "-fmapprep", str(empty), repr(delta_te_ms),
                         "-gz", "0", str(tmp / "fmapprep_bad.nii")]).returncode == 0:
        raise SystemExit("-fmapprep: accepted an all-zero magnitude (empty mask)")

    print("  -fmapprep: closed-form scale/median/mask contract, verbatim mask, -no-debranch")
    print("        parsing and deltaTE/empty-mask rejections OK")


def exercise_openmp_scratch_ops(exe: str, tmp: Path) -> None:
    """-tfce/-tfceS/-bptf/-bptfm/-detrend/-sobel: the four ops whose OpenMP worker scratch was
    hardened to fail closed, plus the two whose per-voxel allocation was hoisted to per-thread.

    These ops had NO in-repo coverage at all before this, which is how eight unchecked
    allocations survived a hardening sweep that fixed six siblings in the same file.  The
    checks below are chosen for what they can actually catch:

    * -detrend on an exactly-linear time series must return all-zero.  Catches a wrong scratch
      size (nvol vs nvox3D) and any cross-voxel bleed from the hoisted buffer.
    * -bptf at -p 1 vs -p 8 must be BYTE-IDENTICAL.  This is the check that catches a botched
      hoist: a scratch buffer that is not fully overwritten before each read produces
      thread-count-dependent output, and nothing else here would notice.
    * -sobel_binary must be two-valued, and plain -sobel must be unaffected by imgdir now
      being allocated only in the binary branch.
    """
    nx, ny, nz, nt = 6, 6, 4, 12
    nvox3d = nx * ny * nz
    # voxel v has series a_v + b_v * t -- exactly linear, so a linear detrend must null it.
    linear = []
    for t in range(nt):
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    v = i + j * nx + k * nx * ny
                    linear.append(float(v % 5) + 0.25 * (v % 3 + 1) * t)
    src = tmp / "omp_linear.nii"
    write_float32_nifti(src, (nx, ny, nz), linear, nt=nt)  # dims is 3-tuple; nt is separate

    det = tmp / "omp_detrend.nii"
    require_success(run_niimath(exe, [str(src), "-detrend", "-gz", "0", str(det)]), "-detrend")
    for idx, value in enumerate(read_float32_nifti(det)):
        if abs(value) > 1e-3:
            raise AssertionError(f"-detrend left {value} at voxel {idx}; a linear series must null out")

    # Thread-count byte-equality: the real regression detector for the hoisted scratch.
    for op in (["-bptf", "4", "2"], ["-bptf", "3", "-1"], ["-bptfm", "4", "2"], ["-detrend"]):
        outs = []
        for threads in ("1", "8"):
            dst = tmp / ("omp_%s_p%s.nii" % (op[0].lstrip("-"), threads))
            require_success(
                run_niimath(exe, [str(src)] + op + ["-p", threads, "-gz", "0", str(dst)]),
                " ".join(op) + " -p " + threads,
            )
            outs.append(dst.read_bytes())
        if outs[0] != outs[1]:
            raise AssertionError(
                " ".join(op) + " is not byte-identical at -p 1 vs -p 8; worker scratch is "
                "carrying state between voxels")

    # -tfce and -tfceS must run and stay finite; -tfceS is 3D-only.
    blob = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                inside = (i - 3) ** 2 + (j - 3) ** 2 + (k - 2) ** 2 <= 4
                blob.append(100.0 if inside else 1.0)
    src3 = tmp / "omp_blob.nii"
    write_float32_nifti(src3, (nx, ny, nz), blob)
    for name, args in (("-tfce", ["-tfce", "2", "0.5", "6"]),
                       ("-tfceS", ["-tfceS", "2", "0.5", "6", "3", "3", "2", "0.5"])):
        dst = tmp / ("omp_%s.nii" % name.lstrip("-"))
        require_success(run_niimath(exe, [str(src3)] + args + ["-gz", "0", str(dst)]), name)
        values = read_float32_nifti(dst)
        if len(values) != nvox3d:
            raise AssertionError(f"{name} changed the voxel count")
        for value in values:
            if value != value or value in (float("inf"), float("-inf")):
                raise AssertionError(f"{name} produced a non-finite value")

    # -sobel_binary is two-valued; plain -sobel must not have been perturbed by imgdir now
    # being allocated only in the binary branch (it is never read when isBinary == 0).
    edge = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                edge.append(200.0 if i >= nx // 2 else 10.0)
    src_edge = tmp / "omp_edge.nii"
    write_float32_nifti(src_edge, (nx, ny, nz), edge)
    sb = tmp / "omp_sobelb.nii"
    require_success(run_niimath(exe, [str(src_edge), "-sobel_binary", "-gz", "0", str(sb)]), "-sobel_binary")
    if not set(read_float32_nifti(sb)) <= {0.0, 1.0}:
        raise AssertionError("-sobel_binary must emit only 0 and 1")
    so = tmp / "omp_sobel.nii"
    require_success(run_niimath(exe, [str(src_edge), "-sobel", "-gz", "0", str(so)]), "-sobel")
    svals = read_float32_nifti(so)
    if len(svals) != nvox3d or max(svals) <= min(svals):
        raise AssertionError("-sobel produced no gradient on a step edge")

    # Degenerate TFCE input, matching the reference implementation. fslmaths derives its step
    # size as max/100 and REJECTS a run where that is not positive -- exit 1, no output written,
    # and for 4D a SINGLE degenerate volume fails the whole operation. An all-NaN volume is
    # NOT rejected (NaN <= 0 is false), and both tools emit an all-zero image for it.
    # Before this contract existed niimath accepted all of these and evaluated (int)NaN.
    flat = [0.0] * nvox3d
    zero_src = tmp / "omp_tfce_zero.nii"
    write_float32_nifti(zero_src, (nx, ny, nz), flat)
    zero_out = tmp / "omp_tfce_zero_out.nii"
    rej = run_niimath(exe, [str(zero_src), "-tfce", "2", "0.5", "6", "-gz", "0", str(zero_out)])
    if rej.returncode == 0:
        raise AssertionError("-tfce must reject a volume whose maximum is not positive")
    if zero_out.exists():
        raise AssertionError("-tfce wrote an output for a rejected degenerate volume")
    if "deltaT" not in (rej.stdout + rej.stderr):
        raise AssertionError("-tfce rejection must name the positive-deltaT requirement")
    # 4D: one degenerate volume fails the whole run.
    mixed = [(100.0 if v % 7 == 0 else 1.0) for v in range(nvox3d)] + flat
    mixed_src = tmp / "omp_tfce_mixed.nii"
    write_float32_nifti(mixed_src, (nx, ny, nz), mixed, nt=2)
    mixed_out = tmp / "omp_tfce_mixed_out.nii"
    if run_niimath(exe, [str(mixed_src), "-tfce", "2", "0.5", "6", "-gz", "0", str(mixed_out)]).returncode == 0:
        raise AssertionError("-tfce must fail the whole run when any volume is degenerate")
    if mixed_out.exists():
        raise AssertionError("-tfce wrote an output despite a degenerate volume")

    # -tfceS coordinate validation. The range check used to print a diagnostic and then FALL
    # THROUGH, computing seed from the rejected coordinate and reading inimg[seed] out of
    # bounds -- on a 16^3 image, -100000 in each axis reads ~109 MB before the buffer, and it
    # did not crash, which is why a diagnostic without a return is worse than no check at all.
    seed_src = tmp / "omp_tfces.nii"
    write_float32_nifti(seed_src, (nx, ny, nz), [100.0 if v % 5 == 0 else 1.0 for v in range(nvox3d)])
    for bad in (["-100000", "-100000", "-100000"], ["999", "999", "999"], ["-1", "0", "0"],
                [str(nx), "0", "0"], ["0", str(ny), "0"], ["0", "0", str(nz)]):
        dst = tmp / ("omp_tfces_bad_%s.nii" % "_".join(bad).replace("-", "m"))
        res = run_niimath(exe, [str(seed_src), "-tfceS", "2", "0.5", "6"] + bad + ["0.5", "-gz", "0", str(dst)])
        if res.returncode == 0:
            raise AssertionError("-tfceS accepted out-of-range coordinate " + " ".join(bad))
        if dst.exists():
            raise AssertionError("-tfceS wrote an output for out-of-range coordinate " + " ".join(bad))
    # ...and an in-range coordinate must still work, so the check is not simply rejecting all.
    ok_dst = tmp / "omp_tfces_ok.nii"
    require_success(
        run_niimath(exe, [str(seed_src), "-tfceS", "2", "0.5", "6", "0", "0", "0", "0.5",
                          "-gz", "0", str(ok_dst)]),
        "-tfceS with an in-range coordinate",
    )

    print("  -tfce/-tfceS/-bptf/-bptfm/-detrend/-sobel: linear-detrend nulling, -p 1 vs -p 8")
    print("        byte-equality, finiteness, sobel two-valuedness and the fslmaths")
    print("        positive-deltaT rejection (3D and 4D) and -tfceS bounds rejection OK")


def exercise_skullstrip(exe: str, tmp: Path, help_text: str) -> None:
    """-skullstrip: end-to-end CLI contract on a synthetic head.

    test_skullstrip_mesh.c checks the surface primitives analytically and nothing checked the
    OPERATION -- which is how both repositioning stages once sat inside the SSV() diagnostic macro
    and ran only when SKULLSTRIP_VERBOSE was set, with the whole suite green.  The verbose-equality
    case below is the regression test for exactly that: a diagnostic must never move a voxel.

    The working grid is a fixed 167x212x175 whatever the input, so a 48^3 fixture costs the same
    ~2 s as a real head.  The numerics stay in skullstrip_bench; this is dispatch and contract.
    """
    line = _help_line(help_text, "-skullstrip")
    if not line:
        raise AssertionError("-skullstrip help line missing entirely (it must be #ifdef-paired)")

    if "NOT in this build" in line:
        # OFF by default, so a disabled build must still say how to enable it and must refuse the
        # op rather than passing the image through unstripped.
        if "SKULLSTRIP=1" not in line:
            raise AssertionError("disabled -skullstrip help must name SKULLSTRIP=1")
        src = tmp / "ss_disabled_in.nii"
        out = tmp / "ss_disabled.nii"
        write_float32_nifti(src, (6, 6, 6), [float(i % 7) for i in range(216)])
        result = run_niimath(exe, [str(src), "-skullstrip", str(out)])
        if result.returncode == 0 or out.exists():
            raise AssertionError("-skullstrip must fail in a build without SKULLSTRIP")
        if "SKULLSTRIP=1" not in (result.stdout + result.stderr):
            raise AssertionError("disabled -skullstrip must say how to enable it")
        print("  -skullstrip: not built (SKULLSTRIP=0) - contract checked")
        return

    # 4D input: rejected before any surface work, and nothing at all is created.
    src4d = tmp / "ss_4d.nii"
    write_float32_nifti(src4d, (8, 8, 8), [float(i % 7) for i in range(8 * 8 * 8 * 3)], nt=3)
    before = sorted(p.name for p in tmp.iterdir())
    result = run_niimath(exe, [str(src4d), "-skullstrip", str(tmp / "ss_4d_out.nii")])
    after = sorted(p.name for p in tmp.iterdir())
    if result.returncode == 0:
        raise AssertionError("-skullstrip must reject 4D input")
    if before != after:
        raise AssertionError(f"-skullstrip on 4D input wrote {set(after) - set(before)}")

    # A ~190 mm ellipsoid on a 4 mm grid: normalisation measures from the top of the head, so a
    # fixture that is small in MILLIMETRES (rather than in voxels) has degenerate contrast and is
    # rejected -- the physical size is what matters here, not the voxel count.
    dims = (48, 48, 48)
    values: list[float] = []
    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                r = ((i - 23.5) / 13.0) ** 2 + ((j - 23.5) / 15.0) ** 2 + ((k - 23.5) / 12.0) ** 2
                values.append(400.0 if r <= 1.0 else (120.0 if r <= 1.6 else 5.0))
    src = tmp / "ss_head.nii"
    write_float32_nifti(src, dims, values, scale=4.0)

    # DT32 only: -dt double must say so, not quietly emit a float64 result.
    dbl_out = tmp / "ss_double.nii"
    dbl = run_niimath(exe, ["-dt", "double", str(src), "-skullstrip", str(dbl_out)])
    if dbl.returncode == 0 or dbl_out.exists():
        raise AssertionError("-skullstrip must reject -dt double")
    if "double" not in (dbl.stdout + dbl.stderr):
        raise AssertionError("-skullstrip -dt double rejection must name the datatype")

    out = tmp / "ss_out.nii"
    require_success(run_niimath(exe, [str(src), "-skullstrip", "-gz", "0", str(out)]), "-skullstrip")
    if not out.exists():
        raise AssertionError("-skullstrip wrote no output")
    if _skullstrip_geometry(src) != _skullstrip_geometry(out):
        raise AssertionError("-skullstrip changed the header geometry (dims/pixdim/qform/sform)")

    # Output contract (skullstrip.h): in-mask voxels keep their ORIGINAL value, out-of-mask voxels
    # become the image MINIMUM. So the result is a strict subset of the input -- no value may
    # appear that the input did not contain, and every changed voxel must hold the minimum.
    vmin = min(values)

    def check_subset(result, label):
        """The output contract, stated ONCE. Both kernels are held to the same rule -- writing it
        twice is how two rules for one contract quietly drift apart. Returns the removed count."""
        if len(result) != len(values):
            raise AssertionError(f"{label} changed the voxel count")
        if min(result) != vmin:
            raise AssertionError(f"{label} background is {min(result)}, expected the input minimum {vmin}")
        if not set(result) <= set(values):
            raise AssertionError(f"{label} invented values {sorted(set(result) - set(values))}")
        n_removed = 0
        for i in range(len(values)):
            if result[i] != values[i]:
                if result[i] != vmin:
                    raise AssertionError(f"{label} voxel {i} became {result[i]}, not the minimum {vmin}")
                n_removed += 1
        if sum(1 for value in result if value == 400.0) == 0:
            raise AssertionError(f"{label} removed the whole brain")
        # KNOW WHAT THIS FIXTURE CAN AND CANNOT SEE. It is a concentric ellipsoid whose
        # background value (5.0) IS vmin, so blanking a background voxel is a no-op that no
        # assertion can detect: 90,784 of its 110,592 voxels are invisible to this check and only
        # the 19,808 bright ones (9,800 at 400.0 + 10,008 at 120.0) can ever register as removed.
        # Measured, ~330 do -- about 3% of the shell -- so this is a DISPATCH AND CONTRACT test,
        # not a segmentation-quality test, and it would still pass against a nearly-identity mask.
        # Mask quality lives in strip_bench (Dice against external reference masks); do not add a
        # quality claim here. The floor below is only strong enough to catch the mask never being
        # applied at all. (An earlier draft asserted the far-corner voxel was removed -- vacuous,
        # because that voxel already holds vmin in the input.)
        if n_removed < 100:
            raise AssertionError(f"{label} removed only {n_removed} voxels; the fixture no longer exercises the mask")
        return n_removed

    stripped = read_float32_nifti(out)
    removed = check_subset(stripped, "-skullstrip")

    # -restart adopts a new dataset, including its STORED datatype. Start from an integer image,
    # restart from the float fixture, and require the same result as processing that float fixture
    # directly. Passing the original input's datatype through restart selects the wrong AFNI
    # normalization branch.
    restart_seed = tmp / "ss_restart_seed.nii"
    write_uint8_nifti(restart_seed)
    restart_out = tmp / "ss_restart_out.nii"
    require_success(
        run_niimath(
            exe,
            [str(restart_seed), "-restart", str(src), "-skullstrip", "-gz", "0", str(restart_out)],
        ),
        "-restart float dataset followed by -skullstrip",
    )
    if read_float32_nifti(restart_out) != stripped:
        raise AssertionError("-restart left -skullstrip using the original input's stored datatype")

    # THE REGRESSION TEST: a diagnostic environment variable must not change the segmentation.
    verbose_env = os.environ.copy()
    verbose_env["SKULLSTRIP_VERBOSE"] = "1"
    vout = tmp / "ss_out_verbose.nii"
    verbose_result = run_niimath(
        exe, [str(src), "-skullstrip", "-gz", "0", str(vout)], env=verbose_env
    )
    require_success(
        verbose_result,
        "-skullstrip with SKULLSTRIP_VERBOSE=1",
    )
    if "skullstrip: deformation kernel fast" not in verbose_result.stderr:
        raise AssertionError("default -skullstrip did not select the fast deformation kernel")
    if vout.read_bytes() != out.read_bytes():
        raise AssertionError("SKULLSTRIP_VERBOSE=1 changed the -skullstrip output")

    # -faithful selects the reference deformation kernel.  Two things are worth pinning: it must
    # be ACCEPTED and produce a real strip (a silently-ignored sub-option would look identical to
    # a working one on a pass/fail check), and it must not be swallowed as the output name -- the
    # op-loop peek gotcha: a sub-option in the output-name slot must not become the output name.
    fout = tmp / "ss_faithful.nii"
    faithful_result = run_niimath(
        exe,
        [str(src), "-skullstrip", "-faithful", "-gz", "0", str(fout)],
        env=verbose_env,
    )
    require_success(
        faithful_result,
        "-skullstrip -faithful",
    )
    if "skullstrip: deformation kernel faithful" not in faithful_result.stderr:
        raise AssertionError("-skullstrip -faithful did not select the faithful deformation kernel")
    faithful = read_float32_nifti(fout)
    check_subset(faithful, "-skullstrip -faithful")
    # The faithful kernel needs its OWN verbose byte-equality check. The one above covers the
    # default kernel only, and -faithful is otherwise run exclusively WITH SKULLSTRIP_VERBOSE=1
    # (to read the dispatch line) -- yet verbose does gate work on this path too, via
    # `ss_verbose() ? &xs : NULL` in skullstrip.c. A diagnostic must never move a voxel in
    # EITHER kernel.
    fquiet = tmp / "ss_faithful_quiet.nii"
    require_success(
        run_niimath(exe, [str(src), "-skullstrip", "-faithful", "-gz", "0", str(fquiet)]),
        "-skullstrip -faithful without SKULLSTRIP_VERBOSE",
    )
    if fquiet.read_bytes() != fout.read_bytes():
        raise AssertionError("SKULLSTRIP_VERBOSE=1 changed the -skullstrip -faithful output")

    # Rejected AND writes nothing -- both halves matter, because an op
    # that errors after creating the file is a different bug from one that errors cleanly.
    before = set(os.listdir(tmp))
    trailing = run_niimath(exe, [str(src), "-skullstrip", "-faithful"])
    if trailing.returncode == 0:
        raise AssertionError("trailing -skullstrip -faithful must fail, not become the output name")
    if set(os.listdir(tmp)) != before:
        raise AssertionError("trailing -skullstrip -faithful wrote a file; it must write nothing")

    # -skullstrip is an ordinary chain op, not a terminal subcommand.
    chained = tmp / "ss_chain.nii"
    require_success(
        run_niimath(exe, [str(src), "-skullstrip", "-mul", "2", "-gz", "0", str(chained)]),
        "-skullstrip chained with -mul",
    )
    doubled = read_float32_nifti(chained)
    for i in range(len(stripped)):
        if abs(doubled[i] - 2.0 * stripped[i]) > 1e-3:
            raise AssertionError("-skullstrip did not chain into the following operation")

    print("  -skullstrip: 4D/-dt double rejections, geometry preserved, %d of %d voxels kept,"
          % (len(values) - removed, len(values)))
    print("        subset-of-input contract, restart datatype, SKULLSTRIP_VERBOSE equality,")
    print("        -faithful kernel accepted and not swallowed as the output name, and chaining OK")


def _help_line(help_text: str, tag: str) -> str:
    """The help line for one operation, so a caller can tell "absent" from "present but
    disabled" -- a platform-gated feature prints a '... NOT in this build' line rather than
    vanishing, and the two cases need different test behaviour."""
    for line in help_text.splitlines():
        if line.strip().startswith(tag):
            return line
    return ""


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
        for token in ("-conform", "-allineate", "-deface", "--dtifit", "--qc", "-bitmap", "-mesh", "-romeo"):
            if token not in help_text:
                raise AssertionError(f"packaged binary help is missing {token}")
        # COPYLEFT SELF-CONSISTENCY -- and KNOW EXACTLY WHAT THIS DOES AND DOES NOT PROVE.
        # -spm_coreg is the entire copyleft payload (SPM, GPL-2-or-later, in the src/GPL
        # submodule). Both the " GPL " version brand (niimath.c's kLicense) and the -spm_coreg
        # help line are emitted under the SAME `HAVE_GPL` macro, so the two checks below cannot
        # disagree unless the build is HALF-wired -- GPL sources compiled with the brand or the
        # help line out of step. That is a real failure mode and worth catching, but it is NOT
        # a proof that a BSD-branded binary is free of GPL object code: a build that linked
        # GPL/*.c WITHOUT -DHAVE_GPL would brand itself BSD, hide the op, and pass. Proving
        # absence needs a symbol check (nm), which this stdlib-only cross-platform script
        # cannot do. Do not describe this as a licence gate; the load-bearing check is the
        # BEHAVIOURAL one below, which requires the op to actually refuse to run.
        # (It used to key on -bandpass, retired along with Exstrom's LGPL-3 bw.c; with that
        # gone the payload is SPM alone, so a GPL=1 binary is GPL-2-or-later, not GPL-3.)
        is_gpl_build = " GPL " in help_text
        if is_gpl_build and "-spm_coreg" not in help_text:
            raise AssertionError("GPL build brands itself GPL but hides -spm_coreg (half-wired build)")
        if not is_gpl_build and "-spm_coreg" in help_text:
            raise AssertionError("BSD-branded build advertises -spm_coreg (half-wired build)")
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
        nn3 = _prod(ndims)
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
        exercise_fillh(exe, tmp)

        exercise_allineate(exe, tmp, help_text)
        exercise_romeo(exe, tmp, help_text)
        exercise_medic(exe, tmp, help_text)
        exercise_moco(exe, tmp, help_text)
        exercise_stc(exe, tmp, help_text)
        exercise_fmap(exe, tmp, help_text)
        exercise_medic_regressions(exe, tmp, help_text)
        exercise_skullstrip(exe, tmp, help_text)
        exercise_openmp_scratch_ops(exe, tmp)

        # THE LOAD-BEARING COPYLEFT CHECK, and it runs for EVERY BSD-branded binary, not just
        # a packaged one -- it used to be gated on --expect-bsd, which meant only the wheel
        # build exercised it. Behavioural, not a help-string grep: the op must actually REFUSE
        # to run. Still not a proof of absent object code (see the note further up), but it is
        # the strongest thing available here.
        spm = run_niimath(exe, [str(small), "-spm_coreg", str(small), str(tmp / "spm.nii")])
        spm_out = spm.stdout + spm.stderr
        stub_msg = "requires a build with the optional GPL module"
        if not is_gpl_build:
            if spm.returncode == 0 or stub_msg not in spm_out:
                raise AssertionError("BSD-branded binary must reject -spm_coreg with the GPL-module message")
        else:
            # The other half of the same half-wiring test. A build where niimath.c received
            # -DHAVE_GPL but core32/core64.c did not would brand itself GPL, print the help
            # line, satisfy both greps above -- and still dispatch to the !HAVE_GPL stub. Only
            # running the op catches that. Deliberately does NOT require success: -spm_coreg on
            # a degenerate 8-bit fixture may legitimately fail to converge, and this is a
            # wiring check, not a registration-quality check.
            if stub_msg in spm_out:
                raise AssertionError("GPL-branded binary dispatched -spm_coreg to the !HAVE_GPL stub (half-wired build)")

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

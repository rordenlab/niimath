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


def nifti_header(dims: tuple[int, int, int], datatype: int, bitpix: int) -> bytes:
    hdr = bytearray(348)
    struct.pack_into("<i", hdr, 0, 348)
    struct.pack_into("<8h", hdr, 40, 3, dims[0], dims[1], dims[2], 1, 1, 1, 1)
    struct.pack_into("<h", hdr, 70, datatype)
    struct.pack_into("<h", hdr, 72, bitpix)
    struct.pack_into("<8f", hdr, 76, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    struct.pack_into("<f", hdr, 108, 352.0)
    struct.pack_into("<f", hdr, 112, 1.0)
    struct.pack_into("<h", hdr, 252, 3)
    struct.pack_into("<h", hdr, 254, 3)
    struct.pack_into("<4f", hdr, 280, 1.0, 0.0, 0.0, 0.0)
    struct.pack_into("<4f", hdr, 296, 0.0, 1.0, 0.0, 0.0)
    struct.pack_into("<4f", hdr, 312, 0.0, 0.0, 1.0, 0.0)
    hdr[344:348] = b"n+1\0"
    return bytes(hdr) + b"\0\0\0\0"


def write_uint8_nifti(path: Path, dims: tuple[int, int, int] = (8, 8, 8)) -> None:
    nvox = dims[0] * dims[1] * dims[2]
    data = bytes((i % 251 for i in range(nvox)))
    path.write_bytes(nifti_header(dims, datatype=2, bitpix=8) + data)


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
        for token in ("-conform", "-allineate", "-deface", "--dtifit", "-bitmap", "-bandpass", "-mesh"):
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
        gz_out = tmp / "roundtrip.nii.gz"
        require_success(run_niimath(exe, [str(small), "-add", "1", "-gz", "1", str(gz_out), "-odt", "char"]), "gzip round-trip write")
        if gz_out.read_bytes()[:2] != b"\x1f\x8b":
            raise AssertionError("gzip output does not have a gzip header")
        assert_payload_size(gz_out, datatype=2, bitpix=8, dims=(8, 8, 8))
        require_success(run_niimath(exe, [str(gz_out)]), "read gzip output")

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

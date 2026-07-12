#!/usr/bin/env python3
"""Deterministically generate the WASM-benchmark fixture: a 256x256x192 int16 T2w-like phantom.

No nibabel dependency: writes a raw NIfTI-1 single-file (.nii) by hand so the output is
byte-reproducible on any machine with the same numpy. The volume has real spatial structure
(overlapping smooth blobs + a background gradient + fine texture) so that -edge/-dog/-s/-fmean
produce meaningful, non-degenerate output. Fixed RNG seed; no wall-clock or environment input.

Usage:  python3 gen_fixture.py [out.nii]     (default: fixtures/t2w_256x256x192_int16.nii)
Prints the SHA-256 of the written file to stdout.
"""
import sys, os, struct, hashlib
import numpy as np

NX, NY, NZ = 256, 256, 192
SEED = 20260711  # frozen; matches the release-baseline date, not a wall clock


def build_volume():
    rng = np.random.default_rng(SEED)
    zz, yy, xx = np.meshgrid(
        np.linspace(-1.0, 1.0, NZ),
        np.linspace(-1.0, 1.0, NY),
        np.linspace(-1.0, 1.0, NX),
        indexing="ij",
    )
    vol = np.zeros((NZ, NY, NX), dtype=np.float64)
    # A handful of smooth Gaussian blobs at frozen centers/scales -> structured intensities.
    blobs = [
        (0.10, 0.00, -0.10, 0.55, 900.0),
        (-0.35, 0.30, 0.20, 0.30, 650.0),
        (0.40, -0.25, 0.05, 0.25, 700.0),
        (0.00, 0.45, -0.30, 0.20, 500.0),
        (-0.20, -0.40, 0.35, 0.18, 480.0),
    ]
    for cz, cy, cx, sig, amp in blobs:
        r2 = (zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2
        vol += amp * np.exp(-r2 / (2.0 * sig * sig))
    # Gentle background gradient + fixed-seed fine texture for high-frequency content.
    vol += 120.0 * (xx + yy + zz)
    vol += rng.normal(0.0, 25.0, size=vol.shape)
    # Clip to a positive int16 range typical of a T2w volume.
    vol = np.clip(vol, 0.0, 4095.0)
    return np.rint(vol).astype("<i2")


def nifti1_header(vox_offset=352):
    h = bytearray(348)
    struct.pack_into("<i", h, 0, 348)                    # sizeof_hdr
    struct.pack_into("<h", h, 40, 3)                     # dim[0] = 3
    struct.pack_into("<h", h, 42, NX)                    # dim[1]
    struct.pack_into("<h", h, 44, NY)                    # dim[2]
    struct.pack_into("<h", h, 46, NZ)                    # dim[3]
    for i in range(4, 8):                                # dim[4..7] = 1
        struct.pack_into("<h", h, 40 + 2 * i, 1)
    struct.pack_into("<h", h, 70, 4)                     # datatype = DT_INT16
    struct.pack_into("<h", h, 72, 16)                    # bitpix
    struct.pack_into("<f", h, 76, 0.0)                   # pixdim[0] (qfac)
    struct.pack_into("<f", h, 80, 1.0)                   # pixdim[1] = 1mm
    struct.pack_into("<f", h, 84, 1.0)                   # pixdim[2]
    struct.pack_into("<f", h, 88, 1.0)                   # pixdim[3]
    struct.pack_into("<f", h, 108, float(vox_offset))    # vox_offset
    struct.pack_into("<f", h, 112, 1.0)                  # scl_slope
    struct.pack_into("<f", h, 116, 0.0)                  # scl_inter
    struct.pack_into("<b", h, 123, 2)                    # xyzt_units = NIFTI_UNITS_MM
    struct.pack_into("<h", h, 252, 1)                    # qform_code = SCANNER_ANAT
    struct.pack_into("<h", h, 254, 1)                    # sform_code = SCANNER_ANAT
    # quaternion identity, offsets center the volume
    struct.pack_into("<f", h, 256, 0.0)                  # quatern_b
    struct.pack_into("<f", h, 260, 0.0)                  # quatern_c
    struct.pack_into("<f", h, 264, 0.0)                  # quatern_d
    struct.pack_into("<f", h, 268, -NX / 2.0)            # qoffset_x
    struct.pack_into("<f", h, 272, -NY / 2.0)            # qoffset_y
    struct.pack_into("<f", h, 276, -NZ / 2.0)            # qoffset_z
    # srow: 1mm iso, matching qform
    struct.pack_into("<4f", h, 280, 1.0, 0.0, 0.0, -NX / 2.0)   # srow_x
    struct.pack_into("<4f", h, 296, 0.0, 1.0, 0.0, -NY / 2.0)   # srow_y
    struct.pack_into("<4f", h, 312, 0.0, 0.0, 1.0, -NZ / 2.0)   # srow_z
    h[344:348] = b"n+1\x00"                              # magic
    return bytes(h)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "fixtures", "t2w_256x256x192_int16.nii")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    vol = build_volume()
    payload = nifti1_header() + b"\x00\x00\x00\x00" + vol.tobytes(order="C")
    with open(out, "wb") as f:
        f.write(payload)
    digest = hashlib.sha256(payload).hexdigest()
    print(f"wrote {out} ({len(payload)} bytes)")
    print(f"sha256 {digest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Generate small deterministic fixtures for the WASI feature-parity suite (dtifit, QC, allineate,
NIfTI-2). All pure-numpy NIfTI writers, fixed seed, no wall clock. Writes into fixtures/.

Outputs:
  fixtures/dwi_6x6x6x7.nii, fixtures/dwi.bvec, fixtures/dwi.bval, fixtures/dwi_mask.nii  (dtifit)
  fixtures/t1_small.nii, fixtures/seg_small.nii                                          (qc)
  fixtures/mov_small.nii, fixtures/tmpl_small.nii, fixtures/mask_small.nii               (allineate/deface)
  fixtures/t2w_small_nifti2.nii                                                          (NIfTI-2)
"""
import os, struct, hashlib
import numpy as np

FX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
os.makedirs(FX, exist_ok=True)
SEED = 20260711


def hdr1(dims, datatype, bitpix, vox_offset=352, pixdim=(1, 1, 1), intent_dim4=1):
    h = bytearray(348)
    struct.pack_into("<i", h, 0, 348)
    nd = len(dims)
    struct.pack_into("<h", h, 40, nd)
    for i, d in enumerate(dims):
        struct.pack_into("<h", h, 42 + 2 * i, d)
    for i in range(nd, 7):
        struct.pack_into("<h", h, 42 + 2 * i, 1)
    struct.pack_into("<h", h, 70, datatype)
    struct.pack_into("<h", h, 72, bitpix)
    struct.pack_into("<f", h, 76, 0.0)
    struct.pack_into("<f", h, 80, pixdim[0])
    struct.pack_into("<f", h, 84, pixdim[1])
    struct.pack_into("<f", h, 88, pixdim[2])
    struct.pack_into("<f", h, 108, float(vox_offset))
    struct.pack_into("<f", h, 112, 1.0)   # scl_slope
    struct.pack_into("<b", h, 123, 2)     # units mm
    struct.pack_into("<h", h, 252, 1)     # qform
    struct.pack_into("<h", h, 254, 1)     # sform
    nx, ny, nz = (list(dims) + [1, 1, 1])[:3]
    struct.pack_into("<f", h, 268, -nx / 2.0)
    struct.pack_into("<f", h, 272, -ny / 2.0)
    struct.pack_into("<f", h, 276, -nz / 2.0)
    struct.pack_into("<4f", h, 280, pixdim[0], 0.0, 0.0, -nx / 2.0)
    struct.pack_into("<4f", h, 296, 0.0, pixdim[1], 0.0, -ny / 2.0)
    struct.pack_into("<4f", h, 312, 0.0, 0.0, pixdim[2], -nz / 2.0)
    h[344:348] = b"n+1\x00"
    return bytes(h)


def write1(path, arr, datatype, bitpix):
    # numpy arrays are indexed [ (t,) z, y, x ]; NIfTI dim[] is [x, y, z, t]. Reverse.
    dims = list(reversed(arr.shape))
    payload = hdr1(dims, datatype, bitpix) + b"\x00\x00\x00\x00" + arr.tobytes(order="C")
    open(path, "wb").write(payload)
    return hashlib.sha256(payload).hexdigest()[:16]


def write2(path, arr, datatype, bitpix):
    """Minimal NIfTI-2 single-file writer (540-byte header)."""
    h = bytearray(540)
    struct.pack_into("<i", h, 0, 540)          # sizeof_hdr
    h[4:12] = b"n+2\x00\r\n\x1a\n"             # magic
    struct.pack_into("<h", h, 12, datatype)    # datatype
    struct.pack_into("<h", h, 14, bitpix)      # bitpix
    shape = list(reversed(arr.shape))          # numpy [z,y,x] -> NIfTI [x,y,z]
    nd = len(shape)
    struct.pack_into("<q", h, 16, nd)          # dim[0]
    for i, d in enumerate(shape):
        struct.pack_into("<q", h, 24 + 8 * i, d)
    for i in range(nd, 7):
        struct.pack_into("<q", h, 24 + 8 * i, 1)
    struct.pack_into("<d", h, 104, 1.0)        # scl_slope
    # pixdim: offset 112, 8 doubles (pixdim[0..7])
    struct.pack_into("<d", h, 112, 0.0)
    for i in range(3):
        struct.pack_into("<d", h, 120 + 8 * i, 1.0)
    struct.pack_into("<q", h, 168, 544)        # vox_offset
    struct.pack_into("<i", h, 344, 1)          # qform_code
    struct.pack_into("<i", h, 348, 1)          # sform_code
    nx, ny, nz = (list(arr.shape) + [1, 1, 1])[:3]
    # srow_x/y/z are doubles at 400/432/464 (4 each)
    struct.pack_into("<4d", h, 400, 1.0, 0.0, 0.0, -nx / 2.0)
    struct.pack_into("<4d", h, 432, 0.0, 1.0, 0.0, -ny / 2.0)
    struct.pack_into("<4d", h, 464, 0.0, 0.0, 1.0, -nz / 2.0)
    struct.pack_into("<i", h, 500, 2)          # xyzt_units mm
    payload = bytes(h) + b"\x00\x00\x00\x00" + arr.tobytes(order="C")
    open(path, "wb").write(payload)
    return hashlib.sha256(payload).hexdigest()[:16]


def gen_dtifit():
    rng = np.random.default_rng(SEED)
    nx = ny = nz = 6
    # frozen 6-direction scheme + one b0
    bvecs = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [0.707, 0.707, 0], [0.707, 0, 0.707], [0, 0.707, 0.707],
    ], dtype=float)
    bvals = np.array([0, 1000, 1000, 1000, 1000, 1000, 1000], dtype=float)
    nvol = len(bvals)
    # a plausible anisotropic tensor per voxel -> monoexponential signal
    D = np.array([[1.7e-3, 0, 0], [0, 0.4e-3, 0], [0, 0, 0.3e-3]])
    S0 = 1000.0 + 200.0 * rng.random((nz, ny, nx))
    dwi = np.zeros((nvol, nz, ny, nx), dtype="<f4")
    for v in range(nvol):
        g = bvecs[v]
        adc = g @ D @ g
        dwi[v] = (S0 * np.exp(-bvals[v] * adc)).astype("<f4")
    # NIfTI stores 4D as x-fastest, volume-slowest -> transpose to (nz,ny,nx,nvol) then C-order?
    # niimath expects planar volume-major: reorder to (t, z, y, x) already; write as t-major.
    write1(os.path.join(FX, "dwi_6x6x6x7.nii"), dwi, 16, 32)
    mask = np.ones((nz, ny, nx), dtype="<i2")
    write1(os.path.join(FX, "dwi_mask.nii"), mask, 4, 16)
    with open(os.path.join(FX, "dwi.bval"), "w") as f:
        f.write(" ".join(f"{b:g}" for b in bvals) + "\n")
    with open(os.path.join(FX, "dwi.bvec"), "w") as f:
        for r in range(3):
            f.write(" ".join(f"{bvecs[v, r]:g}" for v in range(nvol)) + "\n")
    print("dtifit fixtures written")


def gen_qc():
    rng = np.random.default_rng(SEED + 1)
    nx = ny = nz = 16
    zz, yy, xx = np.mgrid[0:nz, 0:ny, 0:nx].astype(float)
    cx = cy = cz = 8.0
    r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2 + (zz - cz) ** 2)
    seg = np.zeros((nz, ny, nx), dtype="<i2")  # 0 = non-brain
    seg[r < 6.5] = 2   # GM (any nonzero not in csf/wm sets)
    seg[r < 4.5] = 3   # WM
    seg[(r >= 6.0) & (r < 6.5)] = 1  # CSF rim
    t1 = np.zeros((nz, ny, nx), dtype="<f4")
    t1[seg == 1] = 300.0
    t1[seg == 2] = 600.0
    t1[seg == 3] = 900.0
    t1 += (10.0 * rng.standard_normal((nz, ny, nx))).astype("<f4")
    t1[seg == 0] = 0.0
    write1(os.path.join(FX, "t1_small.nii"), t1.astype("<f4"), 16, 32)
    write1(os.path.join(FX, "seg_small.nii"), seg, 4, 16)
    print("qc fixtures written (csf=1 wm=3 gm=2)")


def gen_allineate():
    rng = np.random.default_rng(SEED + 2)
    n = 24
    zz, yy, xx = np.mgrid[0:n, 0:n, 0:n].astype(float)
    c = n / 2.0
    def blob(dx, dy, dz):
        r2 = (xx - c - dx) ** 2 + (yy - c - dy) ** 2 + (zz - c - dz) ** 2
        return (1000.0 * np.exp(-r2 / (2 * 5.0 ** 2))).astype("<f4")
    tmpl = blob(0, 0, 0)
    mov = blob(2.0, -1.5, 1.0)  # translated -> registration should recover it
    write1(os.path.join(FX, "tmpl_small.nii"), tmpl, 16, 32)
    write1(os.path.join(FX, "mov_small.nii"), mov, 16, 32)
    mask = (np.sqrt((xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2) < 5.0).astype("<i2")
    write1(os.path.join(FX, "mask_small.nii"), mask, 4, 16)
    print("allineate fixtures written")


def gen_nifti2():
    rng = np.random.default_rng(SEED + 3)
    arr = (rng.integers(0, 2000, size=(8, 8, 8))).astype("<i2")
    write2(os.path.join(FX, "t2w_small_nifti2.nii"), arr, 4, 16)
    print("nifti2 fixture written")


def gen_bold():
    # a small 4D timeseries (>=12 timepoints) so temporal filters (-bandpass) are exercisable
    rng = np.random.default_rng(SEED + 4)
    nx = ny = nz = 4
    nt = 20
    t = np.arange(nt)
    # per-voxel: slow drift + a 0.05 Hz oscillation (TR=2s -> fs=0.5Hz) + noise
    bold = np.zeros((nt, nz, ny, nx), dtype="<f4")
    for tt in range(nt):
        drift = 100.0 + 0.5 * tt
        osc = 20.0 * np.sin(2 * np.pi * 0.05 * tt * 2.0)
        bold[tt] = (drift + osc + 5.0 * rng.standard_normal((nz, ny, nx))).astype("<f4")
    write1(os.path.join(FX, "bold_4x4x4x20.nii"), bold, 16, 32)
    print("bold fixture written")


if __name__ == "__main__":
    gen_dtifit()
    gen_bold()
    gen_qc()
    gen_allineate()
    gen_nifti2()

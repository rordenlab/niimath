#!/usr/bin/env python3
"""simbrain — a dependency-free, parameterized synthetic T1w brain phantom generator.

Stdlib only (no numpy / nibabel): writes NIfTI-1 (.nii or gzipped .nii.gz) with a struct
header. Produces a realistic-enough head for exercising affine registration / defacing:
an egg-shaped head with scalp+skull layers, a folded GM cortex over a WM core, lateral
ventricles, a few subcortical nuclei, a smooth multiplicative bias field, and Gaussian
noise. A rigid pose (translation + rotation) can be baked in so a moving/template PAIR
with a KNOWN transform is trivial to generate for registration tests, and the field of
view can be cropped to a sub-volume to exercise the partial-FOV case (where the fast
engine's intersection-based cost beats an AFNI-style full-overlap assumption).

Examples
  # template (canonical pose) + a brain mask
  python3 simbrain.py --out tmpl.nii.gz --mask brain_mask.nii.gz
  # moving: same anatomy, rotated 8 deg about z and shifted, independent noise
  python3 simbrain.py --out mov.nii.gz --rotate 3 -4 8 --translate 2 -1.5 1 --seed 2
  # a partial-FOV moving (brain-only slab, no neck/skull-base) to stress registration
  python3 simbrain.py --out mov_pfov.nii.gz --fov-z 18 46 --rotate 0 0 6 --seed 3
"""
import argparse
import gzip
import math
import random
import struct


# ---- tissue geometry evaluated in the CANONICAL (template) frame ---------------------
def tissue_intensity(xc, yc, zc, P):
    """Return the noise-free intensity at centered voxel coords (xc,yc,zc), template frame.
    Layers are tested outermost-first so inner structures overwrite outer ones."""
    # Anteroposterior-elongated ("egg") head; slight superior taper.
    ap = 1.30
    taper = 1.0 + 0.010 * yc
    r2 = (xc * xc + (yc / ap) ** 2 + zc * zc) * taper
    r = math.sqrt(r2) if r2 > 0 else 0.0

    # Cortical folding: undulate the GM/WM boundary with a few angular harmonics so the
    # brain carries mid-frequency structure (what real registration keys on), not just a
    # smooth envelope. theta = polar (from +z), phi = azimuth.
    if r > 1e-6:
        theta = math.acos(max(-1.0, min(1.0, zc / r)))
        phi = math.atan2(yc, xc)
    else:
        theta = phi = 0.0
    fold = (math.sin(6.0 * theta) * math.cos(5.0 * phi)
            + 0.6 * math.sin(9.0 * phi) * math.sin(4.0 * theta))
    gm_wm = P["wm_r"] + P["fold"] * fold          # WM core boundary (folded)
    brain = P["brain_r"]                           # GM outer boundary
    skull_o = brain + P["skull"]                   # skull outer
    scalp_o = skull_o + P["scalp_t"]               # scalp outer

    # Neck cylinder (inferior), soft tissue.
    if zc <= -0.34 * P["dim"] and math.sqrt(xc * xc + (yc + 0.03 * P["dim"]) ** 2) <= 0.20 * P["dim"]:
        return P["gm"]

    if r > scalp_o:
        return P["air"]
    if r > skull_o:
        return P["scalp"]
    if r > brain:
        return P["skull"]          # dark skull in T1w
    # inside the brain: sulcal CSF just inside the GM rim
    if r > brain - P["csf_rim"]:
        return P["csf"]
    if r > gm_wm:
        return P["gm"]
    # white matter core, with embedded structures
    # lateral ventricles (paired, curved, CSF)
    for sx in (-1.0, 1.0):
        vx = xc - sx * 0.11 * P["dim"]
        vy = yc / 1.6
        vz = zc + 0.06 * P["dim"]
        if math.sqrt(vx * vx + vy * vy + vz * vz) <= 0.075 * P["dim"]:
            return P["csf"]
    # subcortical nuclei (GM islands in WM: thalamus + basal ganglia)
    for (bx, by, bz, br) in (
        (-0.10, 0.02, 0.0, 0.055), (0.10, 0.02, 0.0, 0.055), (0.0, -0.10, 0.05, 0.06),
    ):
        d = math.sqrt((xc - bx * P["dim"]) ** 2 + (yc - by * P["dim"]) ** 2 + (zc - bz * P["dim"]) ** 2)
        if d <= br * P["dim"]:
            return P["gm"]
    return P["wm"]


def brain_mask_value(xc, yc, zc, P):
    """1 inside the brain (GM+WM envelope), else 0 — template frame (for -deface tests)."""
    ap = 1.30
    taper = 1.0 + 0.010 * yc
    r2 = (xc * xc + (yc / ap) ** 2 + zc * zc) * taper
    return 1 if (r2 <= P["brain_r"] ** 2 and zc > -0.34 * P["dim"]) else 0


def rot_matrix(rx, ry, rz):
    """Rz*Ry*Rx from degrees -> 3x3 (row-major tuples)."""
    ax, ay, az = (math.radians(a) for a in (rx, ry, rz))
    cx, sx = math.cos(ax), math.sin(ax)
    cy, sy = math.cos(ay), math.sin(ay)
    cz, sz = math.cos(az), math.sin(az)
    # Rx
    Rx = ((1, 0, 0), (0, cx, -sx), (0, sx, cx))
    Ry = ((cy, 0, sy), (0, 1, 0), (-sy, 0, cy))
    Rz = ((cz, -sz, 0), (sz, cz, 0), (0, 0, 1))
    def mul(A, B):
        return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    return mul(Rz, mul(Ry, Rx))


def nifti1_bytes(vol, dim, pixdim, dtype):
    """Serialize a flat list `vol` (z-slowest, x-fastest) as NIfTI-1 (+sform/qform code 2)."""
    codes = {"uint8": (2, 8, "<B"), "int16": (4, 16, "<h"), "float32": (16, 32, "<f")}
    dt, bp, pk = codes[dtype]
    h = bytearray(352)
    struct.pack_into("<i", h, 0, 348)
    struct.pack_into("<h", h, 40, 3)
    for i in range(3):
        struct.pack_into("<h", h, 42 + 2 * i, dim)
    struct.pack_into("<h", h, 48, 1)  # dim[4]=1
    struct.pack_into("<h", h, 50, 1)
    struct.pack_into("<h", h, 52, 1)
    struct.pack_into("<h", h, 54, 1)
    struct.pack_into("<h", h, 70, dt)
    struct.pack_into("<h", h, 72, bp)
    struct.pack_into("<f", h, 76, 1.0)  # pixdim[0]
    for i in range(3):
        struct.pack_into("<f", h, 80 + 4 * i, float(pixdim))
    struct.pack_into("<f", h, 108, 352.0)  # vox_offset
    struct.pack_into("<h", h, 252, 2)  # qform_code
    struct.pack_into("<h", h, 254, 2)  # sform_code
    c = (dim - 1) / 2.0 * pixdim
    # qform: identity quaternion, offsets = -center
    struct.pack_into("<f", h, 256, 0.0)  # quatern_b
    struct.pack_into("<f", h, 260, 0.0)
    struct.pack_into("<f", h, 264, 0.0)
    struct.pack_into("<f", h, 268, -c)  # qoffset_x
    struct.pack_into("<f", h, 272, -c)
    struct.pack_into("<f", h, 276, -c)
    for r, off in enumerate((280, 296, 312)):  # srow_{x,y,z}
        row = [0.0, 0.0, 0.0, -c]
        row[r] = float(pixdim)
        for cc in range(4):
            struct.pack_into("<f", h, off + 4 * cc, row[cc])
    h[344:348] = b"n+1\x00"
    body = struct.pack("<%d%s" % (len(vol), pk[1]), *vol)
    return bytes(h) + body


def write_nifti(path, vol, dim, pixdim, dtype):
    data = nifti1_bytes(vol, dim, pixdim, dtype)
    if path.endswith(".gz"):
        data = gzip.compress(data, mtime=0)  # mtime=0 -> reproducible bytes
    with open(path, "wb") as f:
        f.write(data)


def clamp_cast(v, dtype):
    if dtype == "float32":
        return v
    hi = 255 if dtype == "uint8" else 32767
    return int(max(0, min(hi, round(v))))


def generate(P, translate, rotate, fov_z, seed, dtype, want_mask):
    dim = P["dim"]
    rng = random.Random(seed)
    R = rot_matrix(*rotate)  # brain pose R applied about the volume center
    tx, ty, tz = translate
    c = (dim - 1) / 2.0
    zmin, zmax = (0, dim) if fov_z is None else fov_z
    vol = [0.0] * (dim * dim * dim)
    mask = [0] * (dim * dim * dim) if want_mask else None
    for k in range(dim):
        zc0 = k - c
        for j in range(dim):
            yc0 = j - c
            base = (k * dim + j) * dim
            for i in range(dim):
                idx = base + i
                # Outside the requested FOV -> background (partial-FOV moving image).
                if k < zmin or k >= zmax:
                    continue
                xc0 = i - c
                # Map this voxel back into the canonical brain frame: q = R^T (p - t).
                px, py, pz = xc0 - tx, yc0 - ty, zc0 - tz
                xc = R[0][0] * px + R[1][0] * py + R[2][0] * pz
                yc = R[0][1] * px + R[1][1] * py + R[2][1] * pz
                zc = R[0][2] * px + R[1][2] * py + R[2][2] * pz
                val = tissue_intensity(xc, yc, zc, P)
                if val > 0.0:
                    # Smooth multiplicative bias field (realistic MRI inhomogeneity).
                    bias = 1.0 + P["bias"] * (0.6 * (xc / dim) - 0.4 * (yc / dim)
                                              + 0.3 * math.sin(zc / (0.35 * dim)))
                    val *= bias
                    val += rng.gauss(0.0, P["noise"])
                vol[idx] = clamp_cast(val, dtype)
                if mask is not None:
                    mask[idx] = brain_mask_value(xc, yc, zc, P)
    return vol, mask


def build_params(a):
    return {
        "dim": a.dim, "air": a.air, "csf": a.csf, "gm": a.gm, "wm": a.wm,
        "skull": a.skull, "scalp": a.scalp, "noise": a.noise, "bias": a.bias,
        "brain_r": 0.34 * a.dim, "wm_r": 0.24 * a.dim, "fold": 0.03 * a.dim,
        "csf_rim": 0.02 * a.dim, "skull": a.skull, "skull_t": 0.02 * a.dim,
        "scalp_t": 0.03 * a.dim,
    }


def main():
    ap = argparse.ArgumentParser(description="Dependency-free synthetic T1w brain phantom.")
    ap.add_argument("--out", required=True, help="output .nii or .nii.gz")
    ap.add_argument("--mask", help="also write a brain (GM+WM) mask here (template frame)")
    ap.add_argument("--dim", type=int, default=64, help="cubic volume size (default 64)")
    ap.add_argument("--pixdim", type=float, default=3.0, help="isotropic mm (default 3.0)")
    ap.add_argument("--air", type=float, default=0.0)
    ap.add_argument("--csf", type=float, default=10.0)
    ap.add_argument("--gm", type=float, default=128.0)
    ap.add_argument("--wm", type=float, default=250.0)
    ap.add_argument("--skull", type=float, default=40.0, help="skull intensity (dark in T1w)")
    ap.add_argument("--scalp", type=float, default=120.0)
    ap.add_argument("--noise", type=float, default=10.0, help="Gaussian noise std")
    ap.add_argument("--bias", type=float, default=0.25, help="bias-field strength (0=off)")
    ap.add_argument("--translate", type=float, nargs=3, default=(0.0, 0.0, 0.0),
                    metavar=("DX", "DY", "DZ"), help="rigid translation, voxels")
    ap.add_argument("--rotate", type=float, nargs=3, default=(0.0, 0.0, 0.0),
                    metavar=("RX", "RY", "RZ"), help="rigid rotation, degrees")
    ap.add_argument("--fov-z", type=int, nargs=2, default=None, metavar=("ZMIN", "ZMAX"),
                    help="keep only slices [ZMIN,ZMAX) (partial-FOV moving image)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dtype", choices=("uint8", "int16", "float32"), default="float32")
    args = ap.parse_args()

    P = build_params(args)
    vol, mask = generate(P, args.translate, args.rotate, args.fov_z, args.seed,
                         args.dtype, want_mask=bool(args.mask))
    write_nifti(args.out, vol, args.dim, args.pixdim, args.dtype)
    print("wrote %s (%d^3, %s)" % (args.out, args.dim, args.dtype))
    if args.mask:
        write_nifti(args.mask, mask, args.dim, args.pixdim, "uint8")
        print("wrote %s (brain mask)" % args.mask)


if __name__ == "__main__":
    main()

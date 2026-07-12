#!/usr/bin/env python3
"""simbrain — a dependency-free, parameterized synthetic T1w head phantom generator.

Stdlib only (no numpy / nibabel): writes NIfTI-1 (.nii or gzipped .nii.gz) with a struct
header.  The default 64^3 phantom is modelled after ``64c.nii`` with 24 analytic objects
for the cranium, face, neck, cerebral and cerebellar tissues, ventricles, deep nuclei,
eyes, and airways.  A rigid pose (translation + rotation) can be baked in so a
moving/template PAIR
with a KNOWN transform is trivial to generate for registration tests, and the field of
view can be cropped to a sub-volume (``--fov-z``) to generate a partial-FOV moving image.
(No engine "wins" this small synthetic partial-FOV case — real-data behavior differs;
this is only a test-fixture generator, not a benchmark.)

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
# One entry per implicit shape.  Keep this list honest: the model deliberately has a
# small "shape budget" so it stays understandable rather than becoming a voxel atlas.
ANATOMICAL_OBJECTS = (
    "neck", "shoulders", "scalp", "outer_skull", "intracranial_space",
    "left_cerebral_gm", "right_cerebral_gm", "left_cerebral_wm",
    "right_cerebral_wm", "left_cerebellar_gm", "right_cerebellar_gm",
    "left_cerebellar_wm", "right_cerebellar_wm", "brainstem",
    "left_lateral_ventricle", "right_lateral_ventricle", "third_ventricle",
    "left_thalamus", "right_thalamus", "left_eye", "right_eye",
    "nasal_cavity", "oral_airway", "face_and_jaw",
)
assert len(ANATOMICAL_OBJECTS) == 24


def ellipsoid(x, y, z, cx, cy, cz, rx, ry, rz, power=2.0, shear=0.0):
    """Implicit (super-)ellipsoid; ``shear`` tilts y as z changes."""
    dz = (z - cz) / rz
    dx = (x - cx) / rx
    dy = (y - cy - shear * (z - cz)) / ry
    return abs(dx) ** power + abs(dy) ** power + abs(dz) ** power


def anatomy_fields(xc, yc, zc, P):
    """Return the reusable implicit fields for the 24-object head model.

    Geometry is expressed in voxels of the reference 64^3 image and scaled for other
    dimensions.  In the NIfTI array, +x is left/right, +y is anterior, and +z is superior.
    """
    s = P["dim"] / 64.0
    # Global calibration obtained by maximizing NCC against a lightly smoothed 64c.nii.
    # Keeping it here (rather than distorting every object separately) also makes the
    # fitted inferior crop and the slightly narrow left/right field explicit.
    x = 1.0495 * xc / s + 0.5063
    y = 0.8957 * yc / s + 0.5216
    z = 0.7500 * zc / s + 7.8146

    # A superiorly tapered cranium matches the reference better than a sphere.  The
    # taper is folded into x so all three nested cranial objects share the same contour.
    taper = max(0.80, min(1.10, 1.0 - 0.0045 * (z - 1.0)))
    head = ellipsoid(x, y, z, -1.0, -3.0, 1.0, 22.0 * taper, 25.0, 25.5,
                     power=2.15, shear=-0.055)
    skull = ellipsoid(x, y, z, -1.0, -3.2, 1.1, 20.7 * taper, 23.6, 24.0,
                      power=2.12, shear=-0.055)
    intracranial = ellipsoid(x, y, z, -1.0, -3.4, 1.2, 19.5 * taper, 22.4, 22.8,
                             power=2.08, shear=-0.055)

    # Paired lobes make a convincing interhemispheric fissure.  A small periodic term
    # wrinkles the WM interface without spending dozens of extra primitive shapes.
    l_cgm = ellipsoid(x, y, z, -8.0, -4.0, 9.7, 10.7, 19.8, 15.2, shear=-0.035)
    r_cgm = ellipsoid(x, y, z, 6.0, -4.0, 9.7, 10.7, 19.8, 15.2, shear=-0.035)
    folds = (0.050 * math.sin(0.62 * y + 0.18 * z)
             + 0.035 * math.sin(0.78 * z - 0.25 * y)
             + 0.025 * math.cos(0.72 * x + 0.31 * z))
    l_cwm = ellipsoid(x, y, z, -7.4, -4.3, 9.3, 8.2, 15.8, 11.4,
                      power=2.05, shear=-0.025) + folds
    r_cwm = ellipsoid(x, y, z, 5.4, -4.3, 9.3, 8.2, 15.8, 11.4,
                      power=2.05, shear=-0.025) + folds

    fields = {
        "head": head, "skull": skull, "intracranial": intracranial,
        "neck": ellipsoid(x, y, z, -1.0, -1.5, -7.0, 14.0, 19.0, 10.0,
                          power=2.55, shear=0.04),
        "shoulders": ellipsoid(x, y, z, -1.0, 0.0, -16.0, 23.0, 20.0, 6.0,
                               power=2.6),
        "face": ellipsoid(x, y, z, -1.0, 15.0, -11.0, 13.5, 13.5, 16.5,
                          power=2.35, shear=-0.08),
        "l_cgm": l_cgm, "r_cgm": r_cgm, "l_cwm": l_cwm, "r_cwm": r_cwm,
        "l_cbl_gm": ellipsoid(x, y, z, -7.2, -10.0, -3.2, 8.8, 9.5, 6.0),
        "r_cbl_gm": ellipsoid(x, y, z, 5.2, -10.0, -3.2, 8.8, 9.5, 6.0),
        "l_cbl_wm": ellipsoid(x, y, z, -6.8, -9.8, -3.0, 5.3, 6.6, 3.6),
        "r_cbl_wm": ellipsoid(x, y, z, 4.8, -9.8, -3.0, 5.3, 6.6, 3.6),
        "brainstem": ellipsoid(x, y, z, -1.0, -1.5, -3.2, 4.7, 5.5, 9.0,
                               power=2.3, shear=-0.09),
        "l_vent": ellipsoid(x, y, z, -4.8, -3.0, 9.5, 2.2, 6.8, 3.1,
                            power=2.5, shear=0.10),
        "r_vent": ellipsoid(x, y, z, 2.8, -3.0, 9.5, 2.2, 6.8, 3.1,
                            power=2.5, shear=0.10),
        "third_vent": ellipsoid(x, y, z, -1.0, -1.8, 7.2, 1.1, 3.8, 3.7,
                                power=2.2),
        "l_thalamus": ellipsoid(x, y, z, -4.7, -1.5, 6.7, 3.1, 4.4, 2.5),
        "r_thalamus": ellipsoid(x, y, z, 2.7, -1.5, 6.7, 3.1, 4.4, 2.5),
        "l_eye": ellipsoid(x, y, z, -7.3, 17.0, -7.5, 3.8, 4.2, 3.7),
        "r_eye": ellipsoid(x, y, z, 5.3, 17.0, -7.5, 3.8, 4.2, 3.7),
        "nasal": ellipsoid(x, y, z, -1.0, 21.0, -12.0, 3.2, 7.4, 6.6,
                           power=2.4, shear=-0.10),
        "oral": ellipsoid(x, y, z, -1.0, 17.0, -20.0, 5.0, 7.7, 4.7,
                          power=2.5, shear=-0.08),
    }
    return x, y, z, fields


def tissue_intensity(xc, yc, zc, P):
    """Return noise-free T1-like intensity in the canonical template frame."""
    x, y, z, f = anatomy_fields(xc, yc, zc, P)

    # Paint large soft-tissue objects first, then increasingly specific anatomy.
    val = P["air"]
    if f["shoulders"] <= 1.0 or f["neck"] <= 1.0 or f["face"] <= 1.0:
        val = P["scalp"] * 0.92
    if f["head"] <= 1.0:
        # Fat is bright on T1.  This spatial term also avoids a featureless shell.
        fat = 0.22 + 0.18 * (0.5 + 0.5 * math.cos(0.30 * y - 0.18 * z))
        val = P["scalp"] * (1.0 + fat)
    if f["skull"] <= 1.0:
        val = P["skull_i"]
    if f["intracranial"] <= 1.0:
        val = P["csf"] * 0.72

    # The anterior face overlaps the cranial ellipsoids in projection.  Repainting its
    # anterior half prevents the skull shell from turning the lower axial slices into an
    # empty ring, while retaining the bony skull behind it.
    if f["face"] <= 1.0 and y > 7.0:
        val = P["scalp"] * 0.78
    if f["neck"] <= 1.0 and z < 2.5:
        val = P["scalp"] * 0.84

    cerebral = f["l_cgm"] <= 1.0 or f["r_cgm"] <= 1.0
    cerebellar = f["l_cbl_gm"] <= 1.0 or f["r_cbl_gm"] <= 1.0
    if cerebral or cerebellar or f["brainstem"] <= 1.0:
        val = P["gm"]
    if f["l_cwm"] <= 1.0 or f["r_cwm"] <= 1.0:
        val = P["wm"]
    elif cerebral:
        # Sparse analytic sulci break up the cortical ribbon.  They are texture on the
        # two cerebral objects, not extra primitives, so the object budget remains 24.
        cortical_depth = min(f["l_cgm"], f["r_cgm"])
        sulcus = (math.sin(0.82 * x + 0.24 * z)
                  + 0.70 * math.sin(0.57 * y - 0.43 * z)
                  + 0.48 * math.cos(0.91 * x - 0.33 * y))
        if cortical_depth > 0.58 and sulcus > 1.18:
            val = 0.58 * P["gm"] + 0.42 * P["csf"]
    if f["l_cbl_wm"] <= 1.0 or f["r_cbl_wm"] <= 1.0:
        # Alternating folia provide the characteristic striped posterior fossa.
        stripe = 0.84 + 0.16 * abs(math.sin(0.88 * y + 0.42 * z))
        val = P["wm"] * stripe
    if f["brainstem"] <= 1.0:
        val = 0.82 * P["wm"] + 0.18 * P["gm"]

    # Deep structures overwrite white matter; ventricles have final priority there.
    if f["l_thalamus"] <= 1.0 or f["r_thalamus"] <= 1.0:
        val = 0.90 * P["gm"] + 0.10 * P["wm"]
    if f["l_vent"] <= 1.0 or f["r_vent"] <= 1.0 or f["third_vent"] <= 1.0:
        val = P["csf"]

    # Facial structures sit in front of the cranium and therefore paint last.
    if f["l_eye"] <= 1.0 or f["r_eye"] <= 1.0:
        val = 0.82 * P["csf"]
    if f["nasal"] <= 1.0 or f["oral"] <= 1.0:
        val = P["air"]

    # Low-amplitude deterministic tissue texture survives even with --noise 0.
    if val > P["air"] and (cerebral or cerebellar):
        val *= 1.0 + 0.035 * math.sin(0.73 * x + 0.41 * y - 0.29 * z)
    return val


def brain_mask_value(xc, yc, zc, P):
    """1 inside cerebrum/cerebellum/brainstem, else 0 (template frame)."""
    _x, _y, _z, f = anatomy_fields(xc, yc, zc, P)
    return int(f["l_cgm"] <= 1.0 or f["r_cgm"] <= 1.0
               or f["l_cbl_gm"] <= 1.0 or f["r_cbl_gm"] <= 1.0
               or f["brainstem"] <= 1.0)


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
        "skull_i": a.skull, "scalp": a.scalp, "noise": a.noise, "bias": a.bias,
    }


def main():
    ap = argparse.ArgumentParser(description="Dependency-free synthetic T1w brain phantom.")
    ap.add_argument("--out", required=True, help="output .nii or .nii.gz")
    ap.add_argument("--mask", help="also write a brain (GM+WM) mask here (template frame)")
    ap.add_argument("--dim", type=int, default=64, help="cubic volume size (default 64)")
    ap.add_argument("--pixdim", type=float, default=3.75,
                    help="isotropic mm (default 3.75, matching 64c.nii)")
    ap.add_argument("--air", type=float, default=0.0)
    ap.add_argument("--csf", type=float, default=45.0)
    ap.add_argument("--gm", type=float, default=105.0)
    ap.add_argument("--wm", type=float, default=145.0)
    ap.add_argument("--skull", type=float, default=45.0, help="skull intensity (dark in T1w)")
    ap.add_argument("--scalp", type=float, default=95.0)
    ap.add_argument("--noise", type=float, default=7.0, help="Gaussian noise std")
    ap.add_argument("--bias", type=float, default=0.18, help="bias-field strength (0=off)")
    ap.add_argument("--translate", type=float, nargs=3, default=(0.0, 0.0, 0.0),
                    metavar=("DX", "DY", "DZ"), help="rigid translation, voxels")
    ap.add_argument("--rotate", type=float, nargs=3, default=(0.0, 0.0, 0.0),
                    metavar=("RX", "RY", "RZ"), help="rigid rotation, degrees")
    ap.add_argument("--fov-z", type=int, nargs=2, default=None, metavar=("ZMIN", "ZMAX"),
                    help="keep only slices [ZMIN,ZMAX) (partial-FOV moving image)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--dtype", choices=("uint8", "int16", "float32"), default="uint8",
                    help="output voxel type (default uint8, matching 64c.nii)")
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

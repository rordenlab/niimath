"""§7.7 (sweep) -- recover the PHYSICAL displacement wk-apply-warp applies.

For each (grid, letter) we resample three orthogonal ramps (value == voxel
index along i, j, k) through a CONSTANT map of D_MM. The measured per-axis
voxel shift s gives the physical displacement directly:

    out(v) = in(v + s)   and   out(p) = in(p + delta),  delta = A @ s   [RAS mm]

Printing delta for every grid makes the convention readable off the table
instead of guessed.

Run:  ~/src/warpkit/.venv/bin/python exp07b_axis_sweep.py
"""

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nii  # noqa: E402

WK = os.path.expanduser("~/src/warpkit/.venv/bin/wk-apply-warp")
OUT = os.path.expanduser("~/src/niimath/test/medic_ref/exp07b")
N, VOX, D_MM = 24, 2.0, 4.0

E = np.eye(3)


def rot(axis, deg):
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


GRIDS = {
    # name: 3x3 voxel->world RAS. Columns are the i, j, k directions.
    "RAS   j->+y": np.diag([VOX, VOX, VOX]),
    "LAS   j->+y": np.diag([-VOX, VOX, VOX]),
    "flipY j->-y": np.diag([VOX, -VOX, VOX]),
    "swapYZ j->+z": VOX * np.array([E[0], E[2], E[1]]).T,
    "swapXY j->+x": VOX * np.array([E[1], E[0], E[2]]).T,
    "j->-z       ": VOX * np.array([E[0], -E[2], E[1]]).T,
    "j->-x       ": VOX * np.array([E[1], -E[0], E[2]]).T,
    # oblique: separates "displace along the voxel column" from "displace along
    # the canonical world axis the column is closest to" -- identical for every
    # axis-aligned grid above, different here.
    "oblique x20 ": VOX * rot("x", 20.0),
    "oblique z25 ": VOX * rot("z", 25.0),
    "obl x20 LAS ": VOX * rot("x", 20.0) @ np.diag([-1.0, 1.0, 1.0]),
}


def hdr(A):
    return {"pixdim": (1.0, VOX, VOX, VOX, 1.0, 1.0, 1.0, 1.0), "xyzt_units": 10,
            "qform_code": 0, "sform_code": 1, "quatern": (0.0,) * 6,
            "srow": np.hstack([A, np.zeros((3, 1))])}


def measure(A, letter):
    """Return voxel-space shift vector s (out(v) == in(v+s))."""
    h = hdr(A)
    os.makedirs(OUT, exist_ok=True)
    s = np.zeros(3)
    d = np.full((N, N, N), D_MM, np.float32)
    fd = f"{OUT}/d.nii"
    nii.write(fd, d, ref=h)
    for ax in range(3):
        shape = [1, 1, 1]
        shape[ax] = N
        ramp = np.broadcast_to(np.arange(N, dtype=np.float32).reshape(shape), (N, N, N))
        fi, fo = f"{OUT}/in{ax}.nii", f"{OUT}/o{ax}.nii"
        nii.write(fi, ramp, ref=h)
        subprocess.run([WK, "--input", fi, "--transform", fd, "--transform-type", "map",
                        "--phase-encoding-axis", letter, "--output", fo],
                       check=True, capture_output=True)
        o = nii.read(fo)[0]
        idx = np.broadcast_to(np.arange(N, dtype=float).reshape(shape), o.shape)
        sl = (slice(7, N - 7),) * 3
        s[ax] = np.median(o[sl] - idx[sl])
    return s


def main():
    print(f"constant map d = {D_MM} mm; delta = physical displacement (RAS mm), "
          f"out(p) = in(p + delta)\n")
    print(f"{'grid':14s} {'letter':6s} {'voxel shift s':>22s} {'delta_RAS mm':>22s} "
          f"{'delta / d':>16s}")
    for gname, A in GRIDS.items():
        for letter in ("i", "j", "k"):
            s = measure(A, letter)
            delta = A @ s
            print(f"{gname:14s} {letter:6s} {np.array2string(s, precision=2, sign='+'):>22s} "
                  f"{np.array2string(delta, precision=2, sign='+'):>22s} "
                  f"{np.array2string(delta / D_MM, precision=2, sign='+'):>16s}")
        print()


if __name__ == "__main__":
    main()

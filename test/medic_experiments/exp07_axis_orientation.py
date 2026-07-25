"""§7.7 -- what does --phase-encoding-axis <letter> actually select?

Two readings are indistinguishable on the demo data (its voxel j happens to be
the world axis closest to +y):

  (A) letter -> WORLD axis        i/x->+x, j/y->+y, k/z->+z  (RAS)
  (B) letter -> VOXEL index axis  displacement along the image's own column

Permute the header so voxel j points along world +z and voxel k along world +y.
A ramp along voxel j then moves only under reading (B).

Run:  ~/src/warpkit/.venv/bin/python exp07_axis_orientation.py
"""

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nii  # noqa: E402

WK = os.path.expanduser("~/src/warpkit/.venv/bin/wk-apply-warp")
OUT = os.path.expanduser("~/src/niimath/test/medic_ref/exp07")
N, VOX = 24, 2.0

# name -> 3x3 voxel->world (RAS) direction*spacing matrix
GRIDS = {
    "identity  (j -> +y)": np.diag([VOX, VOX, VOX]),
    "permuted  (j -> +z)": np.array([[VOX, 0, 0], [0, 0, VOX], [0, VOX, 0]], float),
    "flipped   (j -> -y)": np.diag([VOX, -VOX, VOX]),
}


def hdr(A):
    h = {"pixdim": (1.0, VOX, VOX, VOX, 1.0, 1.0, 1.0, 1.0), "xyzt_units": 10,
         "qform_code": 0, "sform_code": 1, "quatern": (0.0,) * 6,
         "srow": np.hstack([A, np.zeros((3, 1))])}
    return h


def shift_of(out, ramp_axis):
    """Median (out - index) along ramp_axis over the interior."""
    shape = [1, 1, 1]
    shape[ramp_axis] = N
    idx = np.broadcast_to(np.arange(N, dtype=float).reshape(shape), out.shape)
    sl = (slice(6, N - 6),) * 3
    return float(np.median(out[sl] - idx[sl]))


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"{'grid':22s} {'ramp axis':10s} {'axis arg':8s} {'measured shift':>14s}")
    for gname, A in GRIDS.items():
        h = hdr(A)
        for ramp_axis, aname in ((1, "voxel-j"), (2, "voxel-k")):
            ramp = np.zeros((N, N, N), np.float32)
            shape = [1, 1, 1]
            shape[ramp_axis] = N
            ramp[:] = np.arange(N, dtype=np.float32).reshape(shape)
            d = np.full((N, N, N), 2.0 * VOX, np.float32)  # 2 voxels' worth of mm
            tag = f"{gname[:8].strip().replace(' ', '')}_{aname}"
            fi, fd, fo = f"{OUT}/{tag}_in.nii", f"{OUT}/{tag}_d.nii", f"{OUT}/{tag}_o.nii"
            nii.write(fi, ramp, ref=h)
            nii.write(fd, d, ref=h)
            subprocess.run([WK, "--input", fi, "--transform", fd, "--transform-type", "map",
                            "--phase-encoding-axis", "j", "--output", fo],
                           check=True, capture_output=True)
            o = nii.read(fo)[0]
            print(f"{gname:22s} {aname:10s} {'j':8s} {shift_of(o, ramp_axis):+14.4f}")
    print()
    print("reading (A) world-axis : permuted grid moves the voxel-k ramp, not voxel-j")
    print("reading (B) voxel-axis : permuted grid moves the voxel-j ramp in every grid")


if __name__ == "__main__":
    main()

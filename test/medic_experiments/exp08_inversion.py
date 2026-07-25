"""§7.8 / §7.11 -- how is the undistorted-grid field obtained from the native one?

Working hypothesis, from the exact identity  disp_mm = -f_undist * TRT * voxsize_j
measured on the demo outputs:

    f_undist(y) = f_native(y + f_undist(y) * TRT)      [voxels along the PE axis]

i.e. the undistorted-grid field is the native field sampled at the DISTORTED
location -- a fixed point, which is what "invert the displacement field" reduces
to for a scalar map along one axis. This script solves that fixed point and
compares against warpkit's own sub-fm_fieldmaps.nii.

Run:  ~/src/warpkit/.venv/bin/python exp08_inversion.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nii  # noqa: E402

DEMO = os.path.expanduser("~/src/warpkit/demo/out/")
TRT = 0.02025


def sample_j(vol, pos, mode):
    """Sample vol along axis 1 at fractional index `pos` (same shape as vol)."""
    ny = vol.shape[1]
    if mode == "nearest":
        idx = np.clip(np.rint(pos).astype(np.int64), 0, ny - 1)
        return np.take_along_axis(vol, idx, axis=1)
    lo = np.floor(pos).astype(np.int64)
    f = pos - lo
    a = np.take_along_axis(vol, np.clip(lo, 0, ny - 1), axis=1)
    b = np.take_along_axis(vol, np.clip(lo + 1, 0, ny - 1), axis=1)
    out = a * (1 - f) + b * f
    if mode == "linear_zero":
        out = np.where((pos < 0) | (pos > ny - 1), 0.0, out)
    return out


def main():
    fn = nii.read(DEMO + "sub-fm_fieldmaps_native.nii")[0][..., 0]
    fu_ref, h = nii.read(DEMO + "sub-fm_fieldmaps.nii")[0][..., 0], nii.read(DEMO + "sub-fm_fieldmaps.nii")[1]
    dm_ref = nii.read(DEMO + "sub-fm_displacementmaps.nii")[0][..., 0]
    vox = h["pixdim"][2]
    ny = fn.shape[1]
    j = np.broadcast_to(np.arange(ny, dtype=np.float64)[None, :, None], fn.shape)

    print("fixed point  f(y) <- f_native(y + f(y)*TRT)   [voxels along j]")
    for mode in ("linear_clamp", "linear_zero", "nearest"):
        for sgn in (+1, -1):
            f = np.zeros_like(fn)
            for _ in range(60):
                f = sample_j(fn, j + sgn * f * TRT, "linear_clamp" if mode == "nearest" else mode)
                if mode == "nearest":
                    f = sample_j(fn, j + sgn * f * TRT, "nearest")
            e = np.abs(f - fu_ref)
            print(f"  {mode:13s} sgn{sgn:+d}: p50={np.percentile(e, 50):.5f} "
                  f"p95={np.percentile(e, 95):.5f} p99.9={np.percentile(e, 99.9):.4f} max={e.max():.4f} Hz")

    print("\nconvergence of the winning variant (linear, clamp, +1):")
    f = np.zeros_like(fn)
    for it in (1, 2, 3, 5, 10, 20, 40):
        while True:
            f2 = sample_j(fn, j + f * TRT, "linear_clamp")
            break
        f = f2
        e = np.abs(f - fu_ref)
        print(f"  after {it:3d} iters (cumulative): p95={np.percentile(e, 95):.6f} max={e.max():.4f}")

    print("\ndisplacement identity  disp_mm == -f_undist * TRT * pixdim_j:")
    pred = -fu_ref * TRT * vox
    e = np.abs(pred - dm_ref)
    print(f"  using warpkit's own f_undist: p95={np.percentile(e, 95):.3e} max={e.max():.3e} mm "
          f"(uint16 quantum = {abs(h['scl_slope']):.3e})")


if __name__ == "__main__":
    main()

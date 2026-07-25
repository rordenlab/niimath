"""§7.9 / §7.7 / §7.10 -- wk-apply-warp: sign, interpolation, fill, Jacobian.

Analytic probe: resample a known image through a CONSTANT displacement map and
read off where the content landed. A constant map removes every ambiguity that
real data leaves open (field inversion, per-voxel map interpolation).

Run:  ~/src/warpkit/.venv/bin/python exp09_interp_fill.py
"""

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nii  # noqa: E402

WK = os.path.expanduser("~/src/warpkit/.venv/bin/wk-apply-warp")
OUT = os.path.expanduser("~/src/niimath/test/medic_ref/exp09")
N = 32
VOX = 2.0  # isotropic, axis-aligned: keeps the analytic answer trivial


def run(img, disp, axis, tag, ref=None):
    os.makedirs(OUT, exist_ok=True)
    fi, fd, fo = f"{OUT}/{tag}_in.nii", f"{OUT}/{tag}_disp.nii", f"{OUT}/{tag}_out.nii"
    nii.write(fi, img, ref=ref)
    nii.write(fd, disp, ref=ref)
    cmd = [WK, "--input", fi, "--transform", fd, "--transform-type", "map",
           "--phase-encoding-axis", axis, "--output", fo]
    subprocess.run(cmd, check=True, capture_output=True)
    return nii.read(fo)[0]


def base_header():
    """Axis-aligned RAS header, isotropic VOX mm."""
    _, h = nii.read(nii.write(f"{OUT}/_hdr.nii", np.zeros((N, N, N), np.float32)))
    h["pixdim"] = (1.0, VOX, VOX, VOX, 1.0, 1.0, 1.0, 1.0)
    h["srow"] = np.array([[VOX, 0, 0, 0], [0, VOX, 0, 0], [0, 0, VOX, 0]], float)
    h["qform_code"], h["sform_code"] = 0, 1
    h["quatern"] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    h["xyzt_units"] = 10
    return h


def main():
    os.makedirs(OUT, exist_ok=True)
    h = base_header()
    ramp = np.zeros((N, N, N), np.float32)
    ramp[:] = np.arange(N, dtype=np.float32)[None, :, None]  # value == j index
    impulse = np.zeros((N, N, N), np.float32)
    impulse[16, 16, 16] = 1000.0
    const = np.full((N, N, N), 100.0, np.float32)

    print("== sign / magnitude of the shift (ramp, value == j index) ==")
    for shift_vox in (+2.0, -2.0, +0.5):
        d = np.full((N, N, N), shift_vox * VOX, np.float32)
        o = run(ramp, d, "j", f"ramp{shift_vox:+g}", ref=h)
        # out[j] == ramp[j + s] == j + s  =>  s = out - j
        s = float(np.median(o[8:24, 8:24, 8:24] - np.arange(N)[None, 8:24, None]))
        print(f"  map={shift_vox:+g} vox -> sampled at j{s:+.4f}  (pull = j {'-' if s < 0 else '+'} |d|)")

    print("\n== interpolation kernel (impulse response along j) ==")
    d = np.full((N, N, N), 0.5 * VOX, np.float32)
    o = run(impulse, d, "j", "imp0.5", ref=h)
    prof = o[16, 10:23, 16]
    print("  impulse at j=16, |shift|=0.5 vox; out[j=10..22]:")
    print("   ", np.array2string(prof, precision=5, suppress_small=False))
    nz = np.nonzero(np.abs(prof) > 1e-9)[0]
    print(f"  support width = {len(nz)} taps  (linear=2, cubic B-spline/order3=4, sinc/lanczos>=6)")
    print(f"  negative lobes: {(prof < -1e-9).sum()}  (0 => linear or nearest)")

    print("\n== out-of-FOV fill (pull at j-d, so +d walks off the LOW edge) ==")
    d = np.full((N, N, N), 5.0 * VOX, np.float32)
    o = run(ramp, d, "j", "fill", ref=h)
    print(f"  ramp pulled +5 vox, out[j=0..7] = {np.array2string(o[16, 0:8, 16], precision=3)}")
    print("    (analytic in-FOV value is j-5; j<5 is off the low edge)")
    d = np.full((N, N, N), -5.0 * VOX, np.float32)
    o = run(ramp, d, "j", "fill_neg", ref=h)
    print(f"  ramp pulled -5 vox, out[j=N-8..N-1] = {np.array2string(o[16, N - 8:N, 16], precision=3)}")
    print("    (analytic in-FOV value is j+5; j>N-6 is off the high edge)")

    print("\n== Jacobian modulation (constant image, nonuniform field) ==")
    grad = np.zeros((N, N, N), np.float32)
    grad[:] = (np.arange(N, dtype=np.float32)[None, :, None] - N / 2) * 0.1 * VOX  # 0.1 vox/vox
    o = run(const, grad, "j", "jac", ref=h)
    print(f"  interior values (expect flat 100 if NO modulation): "
          f"{np.array2string(o[16, 12:20, 16], precision=4)}")

    print("\n== axis / sign selector (j vs j-) ==")
    for ax in ("j", "j-", "i", "k"):
        d = np.full((N, N, N), 2.0 * VOX, np.float32)
        o = run(ramp, d, ax, f"ax_{ax.replace('-', 'm')}", ref=h)
        s = float(np.median(o[8:24, 8:24, 8:24] - np.arange(N)[None, 8:24, None]))
        print(f"  axis={ax:2s}: median(out-j) = {s:+.4f}  (nonzero only if the shift lands on j)")


if __name__ == "__main__":
    main()

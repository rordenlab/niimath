"""§7.1 + §7.8 -- synthesise data with an ANALYTICALLY known B0 field and run wk-medic.

Construction: phase_e = wrap(2*pi * f * TE_e), zero phase offset, smooth
magnitude. The weighted regression through the origin must then return exactly
f, so:

  * `_fieldmaps_native` vs f            -> validates phase scaling + unwrap + fit
  * `_fieldmaps` vs the analytic inverse -> validates the inversion
  * `_displacementmaps` vs -f_u*TRT*vox -> validates the Hz->mm convention

For a linear ramp f(y) = a + b*y (y in voxels along the PE axis) the fixed point
f_u(y) = f(y + f_u(y)*TRT) has the closed form  f_u(y) = f(y) / (1 - b*TRT).

Phase is written three ways to settle §7.1 (header scaling vs observed extrema).

Run:  ~/src/warpkit/.venv/bin/python exp01_08_known_field.py
"""

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nii  # noqa: E402

WK = os.path.expanduser("~/src/warpkit/.venv/bin/wk-medic")
OUT = os.path.expanduser("~/src/niimath/test/medic_ref/exp01_08")
NX, NY, NZ = 48, 48, 24
VOX = 3.0
TES_MS = [16.8, 38.56]
TRT = 0.02025
B_HZ_PER_VOX = 6.25  # -> +-150 Hz across the volume, ~3 voxels of displacement


def hdr():
    return {"pixdim": (1.0, VOX, VOX, VOX, 1.0, 1.0, 1.0, 1.0), "xyzt_units": 10,
            "qform_code": 0, "sform_code": 1, "quatern": (0.0,) * 6,
            "srow": np.array([[VOX, 0, 0, 0], [0, VOX, 0, 0], [0, 0, VOX, 0]], float)}


def truth():
    j = np.arange(NY, dtype=np.float64)[None, :, None]
    f = B_HZ_PER_VOX * (j - NY / 2.0) + np.zeros((NX, NY, NZ))
    x, y, z = np.meshgrid(*[np.arange(n, dtype=np.float64) for n in (NX, NY, NZ)], indexing="ij")
    r2 = ((x - NX / 2) / (0.34 * NX)) ** 2 + ((y - NY / 2) / (0.34 * NY)) ** 2 + ((z - NZ / 2) / (0.36 * NZ)) ** 2
    mag = 40.0 + 2600.0 * np.exp(-3.0 * r2)  # smooth "brain", no hard edge
    return f, mag


def wrap(p):
    return (p + np.pi) % (2 * np.pi) - np.pi


def write_inputs(tag, phase_style):
    f, mag = truth()
    h = hdr()
    mags, phases = [], []
    for e, te in enumerate(TES_MS):
        ph = wrap(2 * np.pi * f * (te / 1000.0))
        decay = np.exp(-te / 40.0) / np.exp(-TES_MS[0] / 40.0)
        m = mag * decay
        fm = f"{OUT}/{tag}_e{e + 1}_mag.nii"
        nii.write(fm, m, ref=h, dtype=np.float32)
        fp = f"{OUT}/{tag}_e{e + 1}_phase.nii"
        if phase_style == "float_rad":
            nii.write(fp, ph, ref=h, dtype=np.float32)
        elif phase_style == "siemens_u16":
            # stored 0..4095, scl 2/-4096  -> scaled -4096..4094, like the demo data
            raw = np.rint((ph / np.pi) * 2047.5 + 2047.5).clip(0, 4095)
            nii.write(fp, raw, ref=h, dtype=np.uint16, scl=(2.0, -4096.0))
        elif phase_style == "raw_u16_noscl":
            raw = np.rint((ph / np.pi) * 2047.5 + 2047.5).clip(0, 4095)
            nii.write(fp, raw, ref=h, dtype=np.uint16, scl=(1.0, 0.0))
        mags.append(fm)
        phases.append(fp)
    return mags, phases, f


def run(tag, mags, phases):
    pre = f"{OUT}/{tag}"
    cmd = [WK, "--magnitude", *mags, "--phase", *phases,
           "--TEs", *[str(t) for t in TES_MS],
           "--total-readout-time", str(TRT),
           "--phase-encoding-direction", "j",
           "--out-prefix", pre, "-n", "4"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  !! wk-medic failed: {r.stderr.strip().splitlines()[-1] if r.stderr else r.returncode}")
        return None
    return (nii.read(pre + "_fieldmaps_native.nii")[0][..., 0],
            nii.read(pre + "_fieldmaps.nii")[0][..., 0],
            nii.read(pre + "_displacementmaps.nii")[0][..., 0])


def main():
    os.makedirs(OUT, exist_ok=True)
    core = (slice(8, NX - 8), slice(8, NY - 8), slice(4, NZ - 4))
    for style in ("float_rad", "siemens_u16", "raw_u16_noscl"):
        print(f"\n=== phase written as {style} ===")
        mags, phases, f = write_inputs(style, style)
        got = run(style, mags, phases)
        if got is None:
            continue
        fn, fu, dm = got
        e = np.abs(fn - f)[core]
        print(f"  native vs truth   : p50={np.percentile(e, 50):.4f} p95={np.percentile(e, 95):.4f} "
              f"max={e.max():.4f} Hz   (truth span {f.min():.1f}..{f.max():.1f})")
        fu_analytic = f / (1.0 - B_HZ_PER_VOX * TRT)
        e = np.abs(fu - fu_analytic)[core]
        print(f"  undist vs analytic: p50={np.percentile(e, 50):.4f} p95={np.percentile(e, 95):.4f} "
              f"max={e.max():.4f} Hz   (1/(1-b*TRT) = {1 / (1 - B_HZ_PER_VOX * TRT):.5f})")
        # what scale factor did warpkit actually apply?
        ok = np.abs(f[core]) > 20
        print(f"  measured f_undist/f_native ratio: median={np.median((fu[core] / fn[core])[ok]):.5f}")
        e = np.abs(dm - (-fu * TRT * VOX))[core]
        print(f"  disp vs -f_u*TRT*vox: p95={np.percentile(e, 95):.3e} max={e.max():.3e} mm")


if __name__ == "__main__":
    main()

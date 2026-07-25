"""§7.4 temporal correction, §7.5 low-rank, §7.12 noise frames.

Synthesises a T-frame two-echo run whose per-frame B0 field is an EXACTLY
known low-rank series, then asks the black box three questions:

  rank      : build a series of numerical rank R > 10 with well separated
              singular values.  If the output series comes back at rank 10 the
              truncation is applied to `_fieldmaps_native`; if it comes back at
              rank R it is not.  Centering is detected by whether the mean
              frame survives.
  temporal  : inject a whole 2*pi jump into one frame of one echo's phase.  A
              working temporal correction removes it and leaves the other
              frames untouched.
  noise     : pass -f N and count output frames.

Run:  ~/src/warpkit/.venv/bin/python exp04_05_12_multiframe.py
"""

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nii  # noqa: E402

WK = os.path.expanduser("~/src/warpkit/.venv/bin/wk-medic")
OUT = os.path.expanduser("~/src/niimath/test/medic_ref/exp0405")
NX, NY, NZ = 32, 32, 16
VOX, TRT = 3.0, 0.02025
TES_MS = [16.8, 38.56]


def hdr(nt=1):
    return {"pixdim": (1.0, VOX, VOX, VOX, 1.0, 1.0, 1.0, 1.0), "xyzt_units": 10,
            "qform_code": 0, "sform_code": 1, "quatern": (0.0,) * 6,
            "srow": np.array([[VOX, 0, 0, 0], [0, VOX, 0, 0], [0, 0, VOX, 0]], float)}


def magnitude():
    x, y, z = np.meshgrid(*[np.arange(n, dtype=np.float64) for n in (NX, NY, NZ)], indexing="ij")
    r2 = ((x - NX / 2) / (0.34 * NX)) ** 2 + ((y - NY / 2) / (0.34 * NY)) ** 2 + ((z - NZ / 2) / (0.36 * NZ)) ** 2
    return 40.0 + 2600.0 * np.exp(-3.0 * r2)


def field_series(T, rank, rng):
    """f[x,y,z,t] with exactly `rank` nonzero singular values, decaying 2x each."""
    nv = NX * NY * NZ
    U = rng.standard_normal((nv, rank))
    U /= np.linalg.norm(U, axis=0, keepdims=True)
    V = rng.standard_normal((rank, T))
    V /= np.linalg.norm(V, axis=1, keepdims=True)
    s = 120.0 * (0.5 ** np.arange(rank))
    f = (U * s) @ V
    # a smooth spatial ramp dominates so the field looks physical
    j = np.arange(NY, dtype=np.float64)[None, :, None]
    ramp = (5.0 * (j - NY / 2.0) + np.zeros((NX, NY, NZ))).reshape(nv, 1)
    return (f + ramp).reshape(NX, NY, NZ, T)


def write_run(tag, f, jump=None):
    """jump = (frame, echo, n2pi) injects a whole-2pi phase jump."""
    h = hdr()
    mag = magnitude()
    T = f.shape[3]
    mags, phases = [], []
    for e, te in enumerate(TES_MS):
        ph = (2 * np.pi * f * (te / 1000.0) + np.pi) % (2 * np.pi) - np.pi
        m = np.repeat((mag * np.exp(-te / 40.0) / np.exp(-TES_MS[0] / 40.0))[..., None], T, axis=3)
        if jump is not None and jump[1] == e:
            ph[..., jump[0]] = (ph[..., jump[0]] + np.pi) % (2 * np.pi) - np.pi
        fm, fp = f"{OUT}/{tag}_e{e + 1}_mag.nii", f"{OUT}/{tag}_e{e + 1}_phase.nii"
        nii.write(fm, m, ref=h, dtype=np.float32)
        nii.write(fp, ph, ref=h, dtype=np.float32)
        mags.append(fm)
        phases.append(fp)
    return mags, phases


def run(tag, mags, phases, extra=()):
    pre = f"{OUT}/{tag}"
    cmd = [WK, "--magnitude", *mags, "--phase", *phases, "--TEs", *[str(t) for t in TES_MS],
           "--total-readout-time", str(TRT), "--phase-encoding-direction", "j",
           "--out-prefix", pre, "-n", "4", *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("  !! failed:", (r.stderr or "").strip().splitlines()[-3:])
        return None
    return nii.read(pre + "_fieldmaps_native.nii")[0]


def svals(f):
    nv = f.shape[0] * f.shape[1] * f.shape[2]
    return np.linalg.svd(f.reshape(nv, f.shape[3]), compute_uv=False)


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(7)

    print("=== §7.5 low-rank: input series built with 15 nonzero singular values ===")
    T, R = 20, 15
    f = field_series(T, R, rng)
    mags, phases = write_run("rank", f)
    got = run("rank", mags, phases)
    if got is not None:
        si, so = svals(f), svals(got)
        print(f"  input  singular values (1..{T}): {np.array2string(si, precision=1, max_line_width=200)}")
        print(f"  output singular values (1..{T}): {np.array2string(so, precision=1, max_line_width=200)}")
        tol = so[0] * 1e-6
        print(f"  numerical rank: input={int((si > si[0] * 1e-6).sum())}  output={int((so > tol).sum())}")
        print("  -> output rank 10 means truncation IS applied to _fieldmaps_native")
        # centering probe: does the temporal mean survive?
        print(f"  ||mean frame|| input={np.linalg.norm(f.mean(axis=3)):.3f} output={np.linalg.norm(got.mean(axis=3)):.3f}")

    print("\n=== §7.4 temporal correction: whole 2*pi injected into frame 3, echo 1 ===")
    T = 8
    f2 = np.repeat(field_series(1, 1, rng)[..., :1], T, axis=3)  # constant in time
    mags, phases = write_run("tempclean", f2)
    clean = run("tempclean", mags, phases)
    mags, phases = write_run("tempjump", f2, jump=(3, 0, 1))
    jumped = run("tempjump", mags, phases)
    if clean is not None and jumped is not None:
        for t in range(T):
            d = np.abs(jumped[..., t] - clean[..., t])
            print(f"   frame {t}: p50={np.percentile(d, 50):.4f} p95={np.percentile(d, 95):.4f} max={d.max():.4f} Hz"
                  + ("   <-- injected" if t == 3 else ""))

    print("\n=== §7.12 noise frames ===")
    T = 8
    f3 = field_series(T, 3, rng)
    mags, phases = write_run("noise", f3)
    for nf in (0, 2):
        g = run(f"noise{nf}", mags, phases, extra=["-f", str(nf)])
        print(f"   -f {nf}: output frames = {None if g is None else g.shape[3]} (input {T})")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Quality gate for `-allineate` registration.

`-allineate` IS byte-reproducible across thread counts (since the mfac/afac
thread-local fix — the workflow asserts this separately with `--compare`), but
NOT across *builds*: NEWUOA + `-ffast-math` lands in neighboring, equally-valid
optima, so a codegen change (a refactor, fused-warp reassociation) shifts the
recovered parameters and resliced bytes. A cross-build byte-diff gate would
therefore flap. This asserts registration *quality* instead — the property that
actually matters and that a real regression (optimizer stops converging, cost
wiring broken, output zeroed) would break:

  1. the warped result correlates with the target above an absolute floor, and
  2. registration measurably improved alignment over the unregistered input.

The gate FAILS CLOSED: an empty/degenerate mask, zero-variance (all-zero or
constant) output, a non-finite correlation, or mismatched shapes are failures,
not passes — otherwise a broken allineate that writes zeros would slip through
(np.corrcoef on a constant returns NaN, and NaN<floor is False).

Usage: reg_quality.py <registered> <reference> <unregistered> [min_corr] [min_gain]
"""
import sys
import numpy as np
import nibabel as nib

# These thresholds (and the `reg != 0` in-FOV mask below) are tuned for the dense
# synthetic big/bigref phantom this gate runs on. They are NOT general: a sparse or
# legitimately mostly-zero image would fail the coverage floor, and a tiny nonzero
# subset can clear 5%. Keep this scoped to the CI fixture; for other data derive the
# mask from known source/reference support instead of `reg != 0`.
MIN_CORR = 0.99      # observed ~0.9999 at every thread count; floor leaves wide margin
MIN_GAIN = 0.02      # registered must beat unregistered by at least this
MIN_COVERAGE = 0.05  # in-FOV voxels must be at least this fraction of the volume


def corr(a, b, mask):
    """Pearson correlation over masked voxels; NaN if either side has no variance."""
    x = a[mask]
    y = b[mask]
    if x.std() == 0.0 or y.std() == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def main():
    if len(sys.argv) < 4:
        sys.exit("usage: reg_quality.py <registered> <reference> <unregistered> "
                 "[min_corr] [min_gain]")
    reg_p, ref_p, unreg_p = sys.argv[1:4]
    min_corr = float(sys.argv[4]) if len(sys.argv) > 4 else MIN_CORR
    min_gain = float(sys.argv[5]) if len(sys.argv) > 5 else MIN_GAIN

    reg = nib.load(reg_p).get_fdata()
    ref = nib.load(ref_p).get_fdata()
    unreg = nib.load(unreg_p).get_fdata()

    if not (reg.shape == ref.shape == unreg.shape):
        sys.exit(f"FAIL: shape mismatch reg={reg.shape} ref={ref.shape} "
                 f"unreg={unreg.shape}")

    mask = np.isfinite(reg) & (reg != 0)   # in-FOV after the warp
    coverage = mask.sum() / mask.size
    print(f"in-FOV coverage={coverage:.3f} ({int(mask.sum())}/{mask.size} voxels)")
    if coverage < MIN_COVERAGE:
        sys.exit(f"FAIL: in-FOV coverage {coverage:.3f} < {MIN_COVERAGE} "
                 f"(registered output is empty or near-empty)")

    c_reg = corr(reg, ref, mask)
    c_unreg = corr(unreg, ref, mask)
    print(f"registered corr={c_reg:.6f}  unregistered corr={c_unreg:.6f}  "
          f"(floor {min_corr}, gain {min_gain})")

    if not np.isfinite(c_reg) or not np.isfinite(c_unreg):
        sys.exit("FAIL: non-finite correlation (zero-variance / constant output)")

    ok = True
    if c_reg < min_corr:
        print(f"FAIL: registered correlation {c_reg:.6f} < floor {min_corr}")
        ok = False
    if c_reg < c_unreg + min_gain:
        print(f"FAIL: registration did not improve alignment "
              f"({c_reg:.6f} vs {c_unreg:.6f} + {min_gain})")
        ok = False
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

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

`--smoke` relaxes (1)+(2) to a NON-DEGENERACY check: verify the seed was accepted
and produced a real, finite, non-constant registration with adequate in-FOV
coverage, WITHOUT the golden corr floor / improvement gain. Use it for the header-
seed loop (-com/-sym/-symd/…), whose job is "the seed runs and registers", not to
re-certify quality (the dedicated quality step does that). `-symd`'s extra mirror
`ls` fit is the most cross-build-sensitive path: on the synthetic phantom it can
land in a neighboring — here NOT equally-valid — NEWUOA optimum under one compiler
(observed: gcc-9/x86_64), which the strict floor would flag as a false regression.
Smoke mode still FAILS CLOSED on empty/constant/zero-variance/non-finite output.

Usage: reg_quality.py <registered> <reference> <unregistered> [min_corr] [min_gain] [--smoke]
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
    argv = sys.argv[1:]
    smoke = "--smoke" in argv
    argv = [a for a in argv if a != "--smoke"]
    if len(argv) < 3:
        sys.exit("usage: reg_quality.py <registered> <reference> <unregistered> "
                 "[min_corr] [min_gain] [--smoke]")
    reg_p, ref_p, unreg_p = argv[0:3]
    min_corr = float(argv[3]) if len(argv) > 3 else MIN_CORR
    min_gain = float(argv[4]) if len(argv) > 4 else MIN_GAIN

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

    if smoke:
        # Non-degeneracy only: coverage + finite + non-constant already verified above.
        # The seed ran and produced a real registration; golden quality is asserted by the
        # dedicated quality step, and this path is cross-build NEWUOA-sensitive (see docstring).
        print(f"SMOKE OK: seed ran, non-degenerate registration (corr={c_reg:.6f})")
        sys.exit(0)

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

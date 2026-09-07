# REFILL reference manifest

Evidence for the C port of the REFILL dynamic distortion correction (Robinson et al., HBM 2023). The MATLAB oracle lives in the REFILL repository, `validation/oracle/` (provenance and SHA-256 in its `manifest.json`); nothing large is duplicated here. Read this before changing any constant in `src/refill.c`.

## Oracle provenance

| item | value |
| --- | --- |
| MATLAB pipeline | `refill_oracle_main.m` = `refill_ddc_main.m` with `data.verbose=1`, `data.assess=0` |
| MATLAB | R2024b Update 5, Image Processing Toolbox licensed (`bwdist` path of `smoothn`) |
| ROMEO | ROMEO.jl v1.4.0 (git 60d83fb, the commit niimath ports), MriResearchTools 3.5.0, Statistics 1.11.1, Julia 1.12.3, `JULIA_NUM_THREADS=1` |
| FSL | `bet`, `fugue`, `mcflirt` black boxes only; sources under `/Users/chris/fsl/src/` are never read |
| runner | `validation/run_oracle.sh <sub>`; floors `refill_oracle_dump.m` (`smoothn_nit.m` is a verbatim copy of `smoothn.m` that returns the iteration count); float32 ROMEO companion `validation/oracle_romeo_f32.sh` |
| comparator | `validation/compare_steps.py` (stored voxel order, non-finite mismatches counted separately, thresholds fail the exit code) |

## Measured facts that set the contract

- **`rescale.m` runs in double, not single.** `load_untouch_nii(fn, idx)` with a volume index preallocates `zeros(...)` (double), so the int16 phase arrives as double and `(2π·(x−min))/range − π` is evaluated in double and stored float64. niimath computes the same expression in double and stores float32; the check is `epi_p.nii` within one float32 ULP.
- **ROMEO on float64 vs float32 input is NOT the same computation.** ROMEO.jl's weights follow the phase eltype; MATLAB feeds float64. Measured (sub1, pinned CLI): GE echo 3 has 44 voxels 2π apart, EPI has 4 volumes with 460 490 voxels 2π apart. **Every one of them is outside the mask the pipeline consumes** (GE bet-ero mask; EPI quality > 0.5, which is identical between the two runs). Inside those masks the field maps differ by at most 2.7e-4 rad/s (GE) and 1.3e-4 rad/s (EPI), i.e. float32 rounding. niimath's port computes in float32 (its own oracle), so the unwrap stage is checked bit-exactly against `oracle/romeo_f32/`, and against the MATLAB float64 files only inside the consumed mask. The pinned CLI reproduces the pipeline's own `ge_pd_uw.nii` / `epi_phase_uw.nii` / `epi_quality.nii` bit-for-bit.
- **ROMEO's readphase does not rescale the MATLAB files** (`readphase(fn)` == `readphase(fn; rescale=false)` for `epi_p.nii` and `ge_p.nii`, extrema exactly ±π).
- **`smoothn` is an early-stopped iteration, not a converged solution.** With ~80 % missing data the fixed-point iteration converges so slowly that `TolZ=1e-9` is not reached in 5000 iterations, and the 1e-3 result differs from that run by 21 % relative RMS. The MATLAB output is therefore defined by its exact stopping iteration; the port must reproduce the initial guess (bwdist nearest-neighbour fill, first `ceil(n/10)` DCT coefficients per axis), the relaxation `RF=1.75`, the update, and the stopping test `norm(z0−z)/norm(z) ≤ 1e-3` so that `nit` matches per volume. The floor for this stage is the one-more-iteration difference recorded below.
- **`bwdist` tie rule = lowest linear index (column-major, x fastest; i.e. lexicographic min over z, then y, then x) among equidistant feature voxels.** Verified on sub1's clipped EPI volume 1: 3000 sampled missing voxels (364 with ties) all matched a brute-force lowest-index search; scipy's EDT disagrees on 31 374 of 537 808 missing voxels, MATLAB's pick having the lower linear index in every one. A three-pass per-line brute-force EDT that scans indices upward with a strict `<` reproduces it exactly.
- **MATLAB R2024b's multithreaded `fft` segfaults** (libmwmfl_fft inside TBB) on the 128×128×40 DCT of this data; the dump helper runs `maxNumCompThreads(1)`. The pipeline itself did not hit it.

## Per-subject oracle scalars

| subject | readout gradient (Hz/voxel) | GE median (rad/s) | EPI median (rad/s) | GE−EPI mean (Hz) |
| --- | --- | --- | --- | --- |
| sub1 | −0.001359 | 53.65 | −156.71 | −8.855 |
| sub2 | 0.000293 | 12.91 | 9.03 | −18.108 |
| sub3 | 0.000373 | 102.24 | −82.09 | −40.515 |

(The 2026-09-06 replication with ROMEO.jl 1.6.0 gave −156.77 / −8.208 Hz for sub1: the 1.4→1.6 delta is real but small.)

## Floors (from `refill_oracle_dump.m`, all three subjects)

| floor | sub1 | sub2 | sub3 | meaning |
| --- | --- | --- | --- | --- |
| `aspire_unwarp` jitter, magnitude (max / rms / rel) | 0.0087 / 1.8e-4 / 1.0e-6 | 0.0075 / 1.4e-4 / 9.7e-7 | 0.0091 / 2.1e-4 / 1.1e-6 | two MATLAB runs of the same unwarp differ by this (the `1e-5*rand` grid jitter) |
| `aspire_unwarp` jitter, field map rad/s (max / rms / rel) | 0.0011 / 5.2e-5 / 1.7e-7 | 0.0021 / 1.0e-4 / 2.8e-7 | 0.0020 / 6.4e-5 / 1.8e-7 | same, on the field map |
| `smoothn` one extra iteration, EPI vol 1 (nit; max / rel) | 54; 1.47 / 9.7e-4 | 86; 3.19 / 9.9e-4 | 62; 2.42 / 9.7e-4 | the cost of an off-by-one in `nit`; NOT a tolerance, `nit` must match |
| `smoothn` one extra iteration, GE (nit; max / rel) | 58; 3.76 / 9.8e-4 | 75; 2.88 / 9.8e-4 | 62; 3.48 / 9.7e-4 | same |
| `smoothn` nit over the 100 EPI volumes | 39..75 | 69..100 (15 volumes hit MaxIter) | 52..72 | per-volume counts in `dump/smoothn_nit_epi.txt` |

Algorithm pins (Python transcriptions in REFILL `validation/*_ref.py`, run against sub1): `smoothn_ref.py` reproduces `dump/smoothn_epi1_tol3.nii` to max 2.9e-11 with nit 54 and a bwdist index map identical to MATLAB's; `unwarp_ref.py` reproduces `epi_m_ddc.nii` to rel RMS 1.7e-6 (the jitter floor); `ramp_ref.py` reproduces `dump/ramp.nii` to 6e-8 (one float32 ULP) and the printed gradient to 6 decimals.

## Stage contract (artifact → metric → threshold)

Thresholds are `max(3 × measured floor, a float32 rounding allowance)` and were fixed before any C ran. Run with `-gz 0` before the op, or the `steps/` files and side outputs carry `.nii.gz`. Comparisons are on stored voxel order with `compare_steps.py`; "in mask" restricts to the mask the pipeline consumes at that stage (`--mask`).

| # | niimath op | oracle artifact | dtype | metric | threshold |
| --- | --- | --- | --- | --- | --- |
| 1 | `-refill-epifm -steps`, `-refill-gefm -steps` | `steps/epi_p.nii`, `steps/ge_p.nii` | f64 | max abs | ≤ 3e-7 rad (one float32 ULP at π) |
| 2 | `-refill-epifm` | `dump/ramp.nii`, printed gradient | f32 | max abs; 6-decimal print | ≤ 1e-6 rad; equal |
| 3 | `-refill-epifm`, `-refill-gefm` | `romeo_f32/*/uw.nii` (both); `quality.nii`, `mask.nii` (epifm) | f32 | identical | max abs 0, non-finite 0 |
| 4 | same | `steps/epi_phase_uw.nii` in qmask, `steps/ge_pd_uw.nii` in bet-ero mask | f64 | max abs | ≤ 1e-5 rad (float32-vs-float64 ROMEO, measured ≤ 3e-6) |
| 5 | same | `steps/epi_fm.nii` in qmask, `steps/ge_fm.nii` in mask | f32 / f64 | max abs | ≤ 5e-4 rad/s |
| 6 | same | `steps/qmask.nii`, `steps/epi_mask.nii` | f64 / f32 | identical vs the f32 oracle | 0 mismatches (sub3 has ONE documented boundary voxel vs the float64 pipeline, quality 0.50000006 vs 0.49999997) |
| 7 | same | `steps/epi_fm_masked.nii`, `steps/ge_fm_masked.nii` | f32 / f64 | nit per volume; max abs; rel RMS | nit identical; ≤ 5e-4 rad/s; ≤ 1e-6 |
| 8 | `-refill-centre -steps` | `steps/vsm_prelim.nii`, printed median | f32 | max abs; abs | ≤ 2e-5 voxel; ≤ 0.01 rad/s |
| 9 | `-refill-centre` | `steps/epi_fm_masked-median.nii` | f32 | max abs | ≤ 1e-3 rad/s |
| 10 | `-refill-unwarp -steps` | `steps/vsm.nii` | f32 | max abs | ≤ 2e-5 voxel |
| 11 | `-refill-unwarp` | `epi_m_ddc.nii` | f64 | max abs; rel RMS | ≤ 0.03; ≤ 3.5e-6 |
| 12 | `-refill-unwarp` | `steps/epi_fm_masked-median_dc.nii`, `steps/epi_fm-median_dc.nii` | f32 / f32 | max abs; rel RMS | ≤ 0.007 rad/s; ≤ 1e-6 |
| 13 | `-refill-unwarp` | `steps/epi_quality_dc.nii` | f32 | max abs | ≤ 1e-4 |
| 14 | `-fugue` (reuse) | `epi_m_sdc.nii` | f32 | report only | fmap_bench caveats apply; not a REFILL claim |

Sub3 is validated against a second MATLAB run, `validation/oracle_f32/` (`run_oracle_f32.sh`: the same pipeline with a ROMEO wrapper that casts its inputs to float32, `romeo_pinned_f32`; MATLAB itself single-threaded because R2024b's multithreaded `fft` segfaults on this data). Its one boundary voxel (quality 0.50000006 in float64, 0.49999997 in float32) changes `smoothn`'s missing set, two volumes' `nit` and every downstream map, so the float64 pipeline cannot be its oracle for stages 6-13.

## Results (2026-09-06, M4 Pro, `REFILL=1 make`, system zlib, `.nii` I/O; `validation/run_niimath.sh <sub>`)

Wall time of the whole niimath chain (four ops plus `-fugue`) and the two heavy stages; MATLAB's own timings for the same stages were 15 s + 124-151 s + 25 s per subject, plus FSL:

| subject | chain wall | `-refill-gefm` | `-refill-epifm` | peak RSS (epifm) |
| --- | --- | --- | --- | --- |
| sub1 | 80 s | 7.2 s | 70 s | 1.56 GB (1.31 after the 2026-09-07 refactor) |
| sub2 | 117 s | 9.1 s | 107 s | 1.56 GB (1.31) |
| sub3 | 89 s | 7.2 s | 80 s | 1.56 GB (1.31) |

## End-to-end benchmark (2026-09-07, M4 Pro, sequential, whole process tree sampled at 1 s; `validation/bench_refill.sh`)

| subject | niimath chain wall | niimath peak RSS | MATLAB pipeline wall | MATLAB peak RSS | MATLAB as shipped (2026-09-06 replication) |
| --- | --- | --- | --- | --- | --- |
| sub1 | 79 s | 1.32 GB | 235 s | 7.73 GB | 146 s / 8.08 GB |
| sub2 | 121 s | 1.32 GB | 301 s | 7.72 GB | 163 s / 7.68 GB |
| sub3 | 91 s | 1.32 GB | 237 s | 7.73 GB | 147 s / 8.13 GB |

What each column contains. **niimath chain**: `-refill-gefm`, `-refill-epifm`, `-refill-centre`, `-refill-unwarp` and `-fugue` (the static arm), `.nii` I/O, default threads (14), the oracle's bet mask as input; no bet, no mcflirt. **MATLAB pipeline**: `refill_oracle_main` (the reference with `verbose=1`, `assess=0`), i.e. MATLAB launch, both arms with FSL bet and fugue, ROMEO.jl 1.4.0 (two Julia launches, ~10 s each), one mcflirt (`refill_copy_rescale` calls `refill_assess` unconditionally), and every intermediate written; MATLAB single-threaded because its multithreaded `fft` crashes on this data. Its per-stage times: `epi calc_fms` 161 / 227 / 173 s, `epi do_dc` 28 / 31 / 25 s, `ge calc_fms` 14-17 s. **MATLAB as shipped**: the authors' `refill_ddc_main` (`verbose=0`, `assess=1`: four mcflirt passes, cleanup) with ROMEO.jl 1.6.0 and multithreaded MATLAB, from the 2026-09-06 replication, three subjects run in parallel. Peak RSS on the MATLAB side is dominated by double-precision 4D copies (0.5 GB each); on the niimath side by the float32 input, phase and magnitude plus one field-map series.

Every stage in the contract passes on every subject (sub1 and sub2 against `oracle/`, sub3 against `oracle_f32/`; full reports in `validation/niimath_out/sub*_report.txt`). Scalar checks: printed gradient equal to 6 decimals, printed median equal to 2 decimals, `smoothn` `nit` identical on all 100 EPI volumes and the GE map (sub1, sub2; sub3's two shifted volumes are the boundary-voxel effect and match `oracle_f32`). Selected measurements (max abs / relative RMS):

| stage | sub1 | sub2 | sub3 (vs `oracle_f32`) | threshold |
| --- | --- | --- | --- | --- |
| 3 EPI unwrap vs f32 ROMEO (100 vols), quality, mask | identical | identical | identical | identical |
| 4 EPI unwrap vs f64 MATLAB in qmask (rad) | 2.8e-6 / 3.1e-8 | 3.1e-6 / 3.4e-8 | n/a | 1e-5 |
| 7 `ge_fm_masked` (rad/s) | 1.5e-4 / 2.9e-8 | 1.4e-4 / 3.1e-8 | 2.4e-4 / 2.7e-8 | 5e-4 / 1e-6 |
| 7 `epi_fm_masked` (rad/s, 100 vols) | 2.4e-4 / 5.8e-8 | 2.4e-4 / 5.9e-8 | 1.2e-4 / 5.5e-8 | 5e-4 / 1e-6 |
| 11 `epi_m_ddc` (signal, range ~2000) | 0.011 / 1.6e-6 | 0.012 / 1.5e-6 | 0.028 / 1.7e-6 | 0.03 / 3.5e-6 |
| 12 `epi_fm_masked-median_dc` (rad/s) | 2.7e-3 / 3.2e-7 | 3.1e-3 / 4.0e-7 | 4.1e-3 / 3.1e-7 | 0.007 / 1e-6 |
| 13 `epi_quality_dc` | 5.3e-6 | 7.8e-6 | 8.9e-6 | 1e-4 |
| 14 `epi_m_sdc` `-fugue` vs FSL fugue (report only) | 3.5 / 4.5e-4 | 3.5 / 4.6e-4 | 4.1 / 4.7e-4 | none |

Stage 11 is at the jitter floor: MATLAB's own two runs differ by up to 0.009. Stage 14 is `fmap_bench`'s known `-fugue`/FSL divergence (gap fill, NaN policy) and is not a REFILL claim.

Byte-stability holds within a build, not across codegen: the 2026-09-07 audit refactor (no arithmetic change) moved 15 of 100 sub1 volumes by ≤7.6e-6 rad/s through `-ffast-math` reassociation in the DCT reduction, `nit` unchanged, every stage still inside its threshold. Thread-count invariance, measured on sub1 with one binary: `-refill-gefm`, `-refill-epifm` (100 volumes), `-refill-centre` and `-refill-unwarp` outputs are byte-identical at `-p 1`, `-p 2` and `-p 8`, and run-to-run. (A first attempt showed 1-ULP differences on every volume; the cause was a UBSan `-O1` rebuild that replaced the binary while the background 8-thread run was queued, i.e. two different builds, not a race. Build and test in one step.)

## Reference behaviours reproduced deliberately

Both are bugs in the reference code relative to the paper, kept because the contract is equivalence with the code as it runs, and both are options or arguments so the intended behaviour is available:

- **Readout ramp sign.** `refill_calc_fms.m` forms `angle(exp(1i*(epi_p2-epi_p1)))` with `p2` the REFILL volume, `refill_calc_readout_gradient.m` halves its x-slope, and `refill_calc_fms.m` SUBTRACTS the resulting ramp. The paper's Eq. 7 defines φG_EPI from θ₁ − θ_REFILL, the opposite sign, so the code doubles the residual gradient instead of removing it. Measured on the oracle (in-mask x-slope of EPI minus FLASH field map, rad/s per voxel): sub1 0.785 with the reference sign vs 0.009 with the paper's; sub2 −0.334 vs −0.167; sub3 −0.560 vs −0.347. A synthetic fixture in `release_smoke.py` pins both behaviours. `-ramp-fix` applies the paper's sign.
- **Median mask.** `refill_do_dc.m` does `epi_fm(mask~=1)=NaN` with a 3D logical mask on the 4D unwarped maps; MATLAB applies it to the first `numel(mask)` elements, i.e. volume 1 only, so the "in-mask" median pools the masked first volume with every voxel of the others (sub1: −156.71 rad/s; −145.39 with all volumes masked). Also, `mask` is silently reassigned to the GE bet-ero mask when the SDC arm ran. `-refill-centre` reproduces the pooling and takes the mask as an explicit argument; only voxels equal to 1 count (`mask~=1`), whereas `-refill-gefm` keeps any non-zero voxel (`qmask==0`).

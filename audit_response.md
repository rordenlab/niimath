# Audit Response — Round 5

Response to `audit_temp.md` (Round 5). All Round-4 repairs were re-confirmed. Also audited the maintainer's intentional coreg_fast.c partial-FOV change (see end). Verdicts: **Fixed**, **Resolved** (stale — done between the snapshot and now), **Defer** (pre-existing/own pass).

## High severity

### 1. Output-write failure `stat()` workaround was insufficient → **Fixed properly**
Correct — the stat-after-write could not distinguish a fresh write from a stale file (a read-only pre-existing target read as success; a partial write likewise). Replaced with a real status boundary: **`nifti_image_write` now returns 0/1** (covering header conversion, open, short header/data write, and compressor/close/`pclose` errors — disk-full flushes at close). Its only caller `nifti_save` propagates it (the stat hack is gone), and the main-dispatch / pass-through / intermediate-`-save` sites propagate `nifti_save`. Multi-output `--dtifit` (11 files) and `-tensor_decomp` (9 volumes) now OR every write status (`save_rc |= nifti_save(...)`) and fail if any output failed. Verified: bad directory, **stale read-only file** (the reviewer's specific gap), and normal/gz/zst writes all behave correctly; leak-clean; canonical regression still passes. CI added.

### 2. `--compare` reported non-finite mismatches as equal → **Fixed**
Confirmed (all-NaN vs finite returned 0). Root cause: `fabs(nan/inf)` never beat `maxDiff`, so `differentVox` stayed at its sentinel and both the early "equal" test and the final `maxDiff > thresh` exit read as equal. Fixes: `essentiallyEqual` now handles inf (equal only for identical sign) and one-sided nan (differ); the loop records the FIRST mismatch location regardless of magnitude; the equality test uses `nDifferent == 0`; and a **non-finite mismatch always fails** (`nHardMismatch`). Verified: nan-vs-finite / +inf-vs--inf / inf-vs-finite fail; nan-vs-nan / +inf-vs-+inf equal; finite tolerance compare unchanged. CI added (oracle integrity).

### 3. Whole-program fast-math not implemented → **Resolved (stale)**
The snapshot predated the change. Whole-program `-ffast-math -fno-finite-math-only` now applies in Make, CMake (`CMAKE_C_FLAGS`), notarize.sh (inline), and wasm (main `emcc` line). MSVC stays strict. Validated: Make vs CMake fast affine byte-identical; **both the canonical suite AND `close.sh` (tolerance ops) pass** on the whole-program build (the reviewer's requested Make/CMake/canonical checks). GPL SPM-MATLAB + dtifit/tensor goldens flagged in AGENTS.md for re-validation if touched.

## Medium severity

### 4. Affine JSON boundary accepted malformed input → **Fixed**
Confirmed. Fixes: (a) `nii_apply_affine` now validates the **finiteness** of the bottom row (a `[-nan,0,0,1]` row passed `fabsf(nan)>1e-4` and `al_mat44_usable`'s upper-3-rows check — it wrote an all-NaN image); (b) the JSON reader requires the array to actually **close** (`closed` flag on the outer `]` — a 16-number unterminated array is now rejected); (c) the key search **continues to later occurrences** so an earlier value equal to `"fixed_to_moving"` no longer shadows the real key. Verified all three; valid round-trip unchanged. CI added.

### 5. Documented CMake `-DUSE_OPENMP=OFF` did not exist → **Fixed**
Added one `option(USE_OPENMP ... ON)` gating OpenMP for AppleClang (subordinating the legacy `OPENMP_XCODE`) and GNU. Verified: `-DUSE_OPENMP=OFF` configures with no "unused variable" warning and builds single-threaded; default ON builds with OpenMP. Docs (AGENTS.md/README) now accurate.

### 6. `int` overflow in `nx*ny*nz` → **Defer (pre-existing, Known Issue #4)**
Unchanged; needs a shared checked 3D-count helper as its own pass.

## Low severity
- **7. `al_opts` split** → Defer (the capability mask already prevents the silent-no-op class).
- **8. Serial pyramid resample** → Defer (profile-gated; bit-identical if parallelized).

## Intentional change reviewed — coreg_fast.c partial-FOV cost → **Sound**
For HEL/CR the fast engine now excludes out-of-FOV moving samples from its statistics (was: fill with `m_bg`), so a moving image that doesn't fully cover the template isn't biased toward scale/shear that pull all fixed samples into the moving box. An independent agent verified: **determinism preserved** (the exclusion is inside each fixed reduction chunk; `nin` combined in fixed order), **all divisions guarded** by the `nin < 0.10*ns || nin < 16` floor, bins in bounds, and HEL/CR statistics are internally consistent (no full-set/intersection mixing). Validated on the real motivating case (`T2w.nii.gz` → `avg152T1`, a 60-slice partial-coverage cross-modal pair): the recovered transform is **near-rigid — scales ≈1, inter-axis angles ≈90° (no shear)**, i.e. the pre-change distortion is gone. Added a self-contained partial-FOV CI fixture (a z-slab moving image) asserting both overlap correlation and a near-rigid recovered transform, so a regression that re-distorts partial-FOV fits fails CI.

## Verification
- Make / CMake (incl. `-DUSE_OPENMP=OFF`) / nano / tiny builds clean; canonical + close.sh pass.
- Write-failure (bad dir + stale read-only), `--compare` non-finite, malformed-JSON, seeded replay, and partial-FOV registration all behave as intended; leak-clean.
- Version remains `v1.0.20260711`; nothing committed.

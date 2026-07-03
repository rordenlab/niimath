# Audit Response

## Status

Not ready before this pass. The allocator ownership sweep was directionally correct, but several touched paths still failed open or mutated image state before fallible allocations completed. I corrected the release-blocking issues I found.

## Changes Made

- Reworked `nii_calloc` into `nii_calloc(count, size)` with an overflow check via `nii_mul_size`. Updated current call sites to pass element counts and element sizes explicitly.
- Set the allocation policy for `nii_calloc`: allocation overflow or nonzero OOM prints a clear error and exits with `EXIT_FAILURE`.
- Hardened conform reslicing. `doReslice` now validates dimensions, allocates input/output buffers, and only mutates `nim` header/data after all fallible work succeeds. Allocation failure no longer leaves `nim->data` freed or the header half-updated.
- Hardened high-risk core image operations touched by the allocator sweep: `nifti_crop`, `nifti_dim_reduce`, `nifti_tensor_decomp`, `nifti_subsamp2`, `nifti_resize`, `-pval`, and `-cpval`. These now check key output-size products and allocation results before replacing `nim->data`.
- Fixed a `Tar1` scratch-buffer leak in `nifti_dim_reduce` when a voxel time series is constant.
- Removed per-voxel heap allocations from `nifti_tensor_decomp`; the 6-input and 14-output temporary arrays are now stack locals.
- Hardened allineate setup. Base/source voxel counts are overflow-checked, autoweight/source-mask/base-mask/source-automask scratch allocations fail closed, and source automask backup allocation is mandatory when noise fill is enabled.
- Hardened `create_GA_BLOK_set` realloc paths. Realloc now uses a temporary pointer, failed shrink keeps the original allocation, and all surviving per-block arrays are freed if output struct allocation fails.
- Kept the allineate all-`AL_BIGVAL` coarse-search path fail-closed. It no longer reports a successful registration with only the identity transform.
- Made TypeScript declaration generation fatal in `js/esbuild.config.ts`; a package build should not ship missing or stale `.d.ts` files.

## Rationale

- `nim->data` ownership must be boring: plain `malloc`/`calloc` buffers released by plain `free`. The sweep removed the aligned-allocator mismatch risk, but the helper also needed normal `calloc(count, size)` semantics so large count products are checked consistently. Failed allocation is treated as fatal, not as recoverable application state.
- Image transforms must be transactional at the `nim` level. Allocate first, mutate header/data last.
- Registration and defacing should fail closed on setup/cost-path allocation failure. A degraded identity transform is not an acceptable fallback for privacy-sensitive defacing.
- Release packaging should fail on type-generation errors. Warning-only declaration generation hides broken package artifacts.

## Verification

- `make -C src` passed.
- `git diff --check` passed.
- `bun run build` passed after allowing Emscripten to write its cache outside the workspace. BSD and GPL WASM were rebuilt.
- `bun run esbuild.config.ts` passed. It printed the existing non-fatal `pyenv: cannot rehash` warning.
- `bun test ./tests` passed after the full build: 8 tests, 0 failures. It printed the existing non-fatal `pyenv: cannot rehash` warning.
- Native CLI smoke tests passed on synthetic NIfTI fixtures for `-Tmedian`, `-Tar1`, `-pval`, `-resize`, `-subsamp2`, and `-comply`.

## Remaining Release Risks

- This is not a project-wide allocator-wrapper conversion. MarchingCubes, quadric, and older core paths still have direct allocations; they should eventually route through fail-fast helpers or explicitly `exit(EXIT_FAILURE)` on OOM.
- Several legacy voxel-count products still use `int` in code outside the paths fixed here. The largest-image overflow story is improved, not globally solved.
- `nifti_save` still reports success unconditionally; multi-output operations can still hide write failures.
- Release zstd source download verification remains unpinned. Add a SHA-256 check before relying on that workflow for a release.
- Build-source lists remain duplicated across Makefile, CMake, SuperBuild, and release scripts.

## Supervisor audit + deployment (2026-07-03, v1.0.20260703)

A three-agent audit (security, refactor, docs) plus independent verification reviewed the auditor's hardening ahead of the version bump to `v1.0.20260703` and push to `master`.

**Verdict: no release-blocking regression; the hardening is correct.**
- **Byte-identical output** confirmed between a clean `git HEAD` build and the hardened build across the touched ops on both synthetic and real (`niivue-demo-images/register`) data: `-crop`, `-resize` (shrink/grow/Lanczos), `-subsamp2`, `-conform`, `-Tar1`/dim_reduce, `-pval`/`-cpval`, and `-tensor_decomp` (all 9 outputs). The auditor's transactional/overflow-check rework changes no successful-path result.
- **`nii_calloc(count,size)` + `nii_mul_size`**: overflow check is textbook-correct; all 8 call sites pass `(voxel_count, sizeof)` with no double-counting. `exit(EXIT_FAILURE)` on OOM is **not a regression** (the prior `aligned_calloc` had no NULL check → segfault on OOM; `exit()` is cleaner). **Noted, non-blocking:** the auditor also added `if (dat==NULL) return 1;` handlers at every site, which are unreachable while `nii_calloc` exits; switching `nii_calloc` to *return NULL* would make them live and keep the long-lived WASM worker alive on OOM (the refactor agent's recommendation). Left to the auditor's stated "fatal allocation" design; recommended as a follow-up.
- **`nifti_tensor_decomp`** stack locals are fixed 6/14-float arrays (80 B) — safe, output unchanged. conform `doReslice` reslice math is byte-identical; `create_GA_BLOK_set` realloc-temp and the `Tar1` leak fix are correct; allineate mask/automask fail-closed additions don't alter a successful registration (verified real `-allineate`/`-deface`). No double-free/UAF/leak found in the diff.
- **ASan** could not run a full suite on this Apple Silicon host (platform ASan is impractically slow even via Homebrew LLVM clang; documented in agent memory). Correctness rests on byte-identical diffing + code review.

**CI / Windows confidence (the original `0xC0000374` failure):** the root cause — a `_mm_malloc`/`_aligned_malloc` buffer handed to `nim->data` and freed with plain `free()` on MSVC — is eliminated for the whole tree by the `_mm_malloc`→plain sweep (conform now uses `nii_calloc`). `release_smoke.py --expect-bsd --expect-zstd` **passes locally on both the Makefile and the CMake (Windows build-path) binaries**; no build system references the deleted `arm_malloc.h`; the diff contains no MSVC-incompatible constructs (no VLA/`typeof`/`__attribute__`/statement-expressions; `SIZE_MAX` via `<stdint.h>`). Native, CMake, nano, tiny, and WASM all build clean. Pushing to `master` runs `release.yml` (Windows wheels + `release_smoke`) but does **not** publish (PyPI upload is gated on `refs/tags/v*`), so this is CI validation without a release — a maintainer tags `v1.0.20260703` later to publish.

**Highest-value open follow-ups for next session:** (1) finish Known Issue #4 — the remaining `int nvox3D` **divisor/loop-bound** truncation sites (SIGFPE/OOB on >2³¹-voxel images) not covered by the allocation-sizing fix; (2) optionally switch `nii_calloc` to return-NULL for WASM-worker survival; (3) pin the zstd tarball SHA-256 in `release.yml`.

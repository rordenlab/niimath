This file provides guidance for AI agents when working with code in this repository.

## Project Overview

niimath is an open-source clone of FSL's `fslmaths` — a general-purpose NIfTI image calculator for neuroimaging. It extends fslmaths with mesh generation (nii2mesh), additional filters, defacing (`-deface`), affine registration (`-allineate`), anatomical QC (`--qc`), zstd compression (.nii.zst), and cross-platform support (Linux, macOS, Windows, WebAssembly).

**Repository:** `rordenlab/niimath` (BSD-2-Clause license)

## Building

### Quick build (Makefile, recommended for development)
```bash
cd src && make          # Standard optimized build (OpenMP enabled by default)
make debug              # Debug build (-g, no optimization)
make ubsan              # UndefinedBehaviorSanitizer (lightweight; retains OpenMP)
make sanitize           # AddressSanitizer build (serial on Apple Clang; no leak detection on macOS)
make verbose            # All warnings enabled
make static             # Static binary
```

### Prerequisites
- **macOS**: `brew install libomp zstd` (libomp for OpenMP — enabled by default; zstd for .nii.zst support)
- **Linux**: `apt install libzstd-dev` (or equivalent; OpenMP is built-in with gcc)

### Build variants
```bash
make tiny               # Minimal terminal build (emulates WASM constraints)
make nano               # Without mesh (nii2mesh) functions
MESH=0 make             # Disable mesh support
AL=0 make               # Disable allineate registration
DTIFIT=0 make           # Disable dtifit (diffusion tensor fit)
QC=0 make               # Disable anatomical QC (--qc)
GPL=1 make              # Enable optional GPL -spm_coreg (needs git submodule update --init src/GPL)
ZSTD=0 make             # Disable zstd (.nii.zst) compression support
OMP=0 make              # Disable OpenMP (single-threaded build)
CF=1 make               # CloudFlare accelerated zlib
MC=1 make               # Use new MarchingCubes algorithm (handles ambiguities)
STB=1 make              # Use STB image library instead of libspng for bitmaps
make wasm               # Emscripten/WebAssembly target
```

### CMake build
```bash
mkdir build && cd build && cmake .. && make     # top-level SuperBuild
cmake -DUSE_OPENMP=OFF -DENABLE_QC=OFF ..        # feature toggles (forwarded by SuperBuild)
```

### Key compile-time flags
- `NII2MESH` / `DNII2MESH` — enables mesh conversion features
- `USE_CLASSIC_CUBES` — selects oldcubes.h (classic marching cubes) vs MarchingCubes.h
- `HAVE_FORMATS` — enables GIfTI, OBJ, VTK, STL mesh output
- `HAVE_BMP` — bitmap/PNG creation
- `HAVE_BUTTERWORTH` — bandpass temporal filtering
- `HAVE_TENSOR` — tensor decomposition
- `HAVE_DTIFIT` — linear diffusion tensor fit (`--dtifit`), emulates FSL dtifit; needs `HAVE_TENSOR`
- `HAVE_QC` — anatomical QC metrics (`--qc`); self-contained, no tensor/mesh deps
- `HAVE_CONFORM` — image conforming to standard space
- `HAVE_ALLINEATE` — affine registration and defacing (allineate.c + powell_newuoa.c + coreg_fast.c); release builds use whole-program -ffast-math and OpenMP
- `HAVE_GPL` — optional GPL `-spm_coreg`/`-spm_deface` (SPM rigid-body coregistration/defacing); OFF by default so builds stay BSD-2. `-spm_coreg` is the canonical flag (`-spmcoreg` is a silent backward-compat alias). Sources in the `niimath_gpl` submodule (`src/GPL`); enable with `GPL=1 make` / `cmake -DENABLE_GPL=ON` / `GPL=1 make wasm` (needs `git submodule update --init src/GPL`). A `GPL=1` binary is a GPL-2 combined work (version string ends ` GPL` vs ` BSD`). See the **src/GPL** module note below.
- `AL_LPC_MICHO` — enables lpc+ZZ/lpa+ZZ combined cost variant for allineate (default: pure lpc/lpa)
- `HAVE_ZSTD` — zstd compression support for .nii.zst files (auto-detected; set `FSLOUTPUTTYPE=NIFTI_ZST` to write)
- `FSLSTYLE` — FSL-compatible behavior mode

## Source Architecture

### Core computational pipeline
- **niimath.c** — CLI entry point, dispatches to `main32()` or `main64()` based on `-dt` flag. `--dtifit`/`--qc` are self-contained subcommands dispatched here before the op loop; `main32`/`main64` reject `nx*ny*nz` or total `nvox > INT_MAX` (via the shared `nii_nvox3d_int()`) before operator dispatch.
- **coreFLT.c** (6k lines) — Main computational engine, compiled twice via template pattern:
  - **core32.c** — `#define DT32` + `#include "coreFLT.c"` → float32 (SSE 4-wide SIMD on x86_64; scalar on ARM/WASM)
  - **core64.c** — includes coreFLT.c without DT32 → float64 (SSE 2-wide SIMD on x86_64; scalar on ARM/WASM)
- **core.c** — Shared utilities: datatype conversion, kernel creation, Otsu thresholding, resampling filters, NIfTI I/O helpers. **`max_displacement_mm()` gotcha:** it returns the max corner displacement between two images' world transforms and gates the `-add`/binary-op and `--compare` orientation warnings (0.5 mm) and the `--qc` grid check (0.001 mm). It measures true 3D Euclidean distance AND normalises each transform to mm via `xyz_units` (`xyz_units_to_mm`), so metre/micron headers compare correctly against mm thresholds. mm images are unaffected (factor 1.0); touch this helper only with that unit scaling in mind.
- **unifize.c** — Bias field correction via `-unifize` flag (adapted from AFNI 3dUnifize, public domain)
- **allineate.c** (~3.9k lines) — Affine registration, defacing, skull-stripping. **Shared byte-identically with the standalone `/Users/chris/src/allineate` project — it is the canonical source; sync any change both ways** (see Shared-code invariant below). Cost functions via `-cost` (Hellinger default, lpc/lpa/ls; `-DAL_LPC_MICHO` adds lpc+ZZ/lpa+ZZ); DOF via `-warp`; match/output interpolation via `-interp`/`-final`. Twopass coarse-to-fine with a full-resolution AFNI-faithful random-startup coarse search (see the AFNI-fidelity note below) and OpenMP-parallel candidate refinement using thread-local histogram/warp/workspace buffers (adapted from AFNI 3dAllineate, public domain) — this thread-local-cost-eval pattern is the template for any spm_coreg OpenMP work. `al_register()` is shared by `nii_allineate()`/`nii_deface()`; `nii_reslice_affine` + `nii_apply_deface_mask` are factored out and reused by `nii_deface` and the GPL `-spm_coreg`/`-spm_deface` (interpolation lives once in BSD). **`-deface` is the single mask-based removal command** — the supplied mask, not the command, determines what is removed: a brain mask skull-strips, a face mask defaces. **Registration direction matters:** register the SUBJECT onto the TEMPLATE (`al_register(input, tmpl)`, base=template — the well-posed direction), then INVERT the transform (`GA_setup_affine`→`nifti_mat44_inverse`) to pull the template-space mask onto the subject's native grid via `nii_reslice_affine`. Registering a brain-only template ONTO a full-head subject (base=subject) converges to a mislocated transform (the subject's neck/shoulders/FOV dominate the cost) and masks the wrong region — do not do it. **`nii_deface` rejects 4D input** (`al_dims_ok`): the face mask covers one volume, so 4D would silently leave faces in volumes 2..N — a privacy failure. Options in the `al_opts` struct (`allineate.h`).
- **coreg_fast.c/.h** — **Fast affine coreg engine** (`-cost fast` → Hellinger/MI cross-modal default; `-cost fastcr` → correlation-ratio). Independently implemented (clean-room, NOT AFNI/GPL). **Shared with the standalone allineate project** (sync both ways). `coreg_fast_estimate()` does not mutate inputs — returns a world-mm FIXED→MOVING `mat44`; the caller applies it once. Cost reductions use fixed `CF_CR_NCHUNK=64` disjoint chunks combined in fixed order → deterministic per team size. **Build gotcha:** coreFLT.c references `coreg_fast_estimate` unconditionally under `HAVE_ALLINEATE`, so every allineate-enabled build MUST link it (native `coreg_fast.o`, wasm `cf_wasm.o`, CMake `ADDITIONAL_SRCS`, notarize `AL_SRCS`). Blurs its pyramid via the exported `nifti_smooth_gauss_f32`/`_f64`; multi-volume smoothing passes per-volume dims — **do NOT introduce shared `nx`/`ny`/`nz` here (a 4D OpenMP race).** **Partial-FOV cost:** HEL/CR form statistics from the image INTERSECTION — out-of-FOV moving samples are excluded (`if(!ok) continue`), NOT filled with `m_bg` — else a partial-coverage moving image (e.g. a cropped T2w onto a full template) distorts toward scale/shear that pull all fixed samples into the moving box. A floor (`nin < 0.10*ns || nin < 16 → CF_PENALTY`) rejects a degenerate tiny overlap; the exclusion is inside each fixed reduction chunk so determinism holds; LS keeps the `m_bg` fill. **Initialization:** default/`-cmass` scores the supplied affine and exact-COM starts once at 8 mm, maximizing `(1-cost)×covered_fixed_foreground`, then runs only the winning descent; `-com` forces a recentered header, `-nocmass` forces the supplied affine. Translation guard is ±128 mm. **Contract:** fixed 12-DOF schedule; fail-closed rejects `-warp`/`-interp`/`-source_automask`/`-dark_automask`/`-sym`/`-zoom`; `-cost` last-one-wins; both-form-codes-zero input uses the shared pixdim-centered fallback. **Reslice fill safety:** cubic reslicing clamps only in-FOV voxels (`nii_reslice_affine` inline; `al_scalar_warpone` via `al_clip_fused_cubic_infov`) so a positive-only source's out-of-FOV fill stays 0 — clamping it up to the source min is a defacing mask-safety hazard.
- **powell_newuoa.c** (~2.8k lines) — Powell's NEWUOA derivative-free optimizer (f2c translation, used by allineate). **Shared with the standalone allineate project.** Thread-safe statics via `__thread`; `powell_newuoa_free_threadlocal()` releases the per-thread workspace (`pn_w`), called from `al_register`'s end-of-registration teardown. **KEY gotcha — `mfac`/`afac` (NEWUOA sampling factors, set by `powell_set_mfac`) are `__thread`**, so a main-thread `powell_set_mfac()` does NOT configure OpenMP workers. Every parallel `powell_newuoa()` loop MUST capture the factors with `powell_get_mfac()` before the region and re-apply `powell_set_mfac()` inside the body — otherwise workers use default/stale `npt` and the result becomes thread-count-dependent (a real reproducibility trap). **Fail-closed OOM:** `pn_w` alloc failure returns `-7` rather than calling `newuoa_` with NULL; this propagates — `al_scalar_optim()` returns it (and `-3` on its own `wpar` calloc failure), `al_register()` aborts via `goto al_cleanup` (no "Registration complete", no result emitted), and parallel candidate loops mark a failed optimization `AL_BIGVAL` so it cannot be selected. The cost path is likewise fail-closed: `GA_get_warped_values()` returns 1 on a thread-local buffer OOM and `GA_scalar_fitter()` returns `AL_BIGVAL` rather than scoring a stale `avm`; `al_scalar_setup()` leaves `stup->setup = 0` on an image-buffer OOM, which the three `al_register` call sites verify (`stup.setup != AL_SMAGIC → goto al_cleanup`). An optimizer/setup OOM is a clean error exit, never a silent unrefined registration.

### Mesh code (nii2mesh — unique to niimath, not in FSL)
- **meshify.c** — Main mesh pipeline: smoothing → marching cubes → vertex welding → degenerate removal → export
- **MarchingCubes.c/.h** — Newer implementation with ambiguity resolution (14 cases with subcases)
- **oldcubes.c/.h** — Classic/simpler marching cubes (~500 lines, selected with `USE_CLASSIC_CUBES`)
- **quadric.c** — Mesh simplification via quadric error metrics
- **meshtypes.h** — `vec3d` (double xyz), `vec3i` (int xyz) structs

### Supporting modules
- **dtifit.c** — `niimath --dtifit`: linear diffusion tensor fit emulating FSL `dtifit`. Self-contained TU dispatched early in `main()`, FSL-identical flags (`-k/-m/-r/-b/-o`, `-xflip 0|1|auto`). 7-param (S0 + 6 tensor) log-linear OLS (fit math from AFNI 3dDWItoDT linear path, public domain; nonlinear omitted), reuses `EIG_tsfunc` (tensor.c) for FA/MD/L*/V*. Emulates FSL's determinant-based bvec x-flip. Writes `<base>_{FA,MD,L1,L2,L3,V1,V2,V3,S0,MO,tensor}`; unsupported FSL features (`--wls`/`--kurt`) error clearly. **Gotchas:** tensor buffer must use NIfTI planar (volume-major) layout; validate V* with \|cos angle\| not Pearson r (eigenvector sign is arbitrary).
- **qc.c** — `niimath --qc <t1> --seg <seg> --csf <i[,j..]> --wm <i[,j..]> [--erode 0|1] [--out qc.tsv]`: MRIQC-style anatomical Image Quality Metrics from a T1 + a hard integer segmentation, written as a wide TSV (header row + one value row). Self-contained TU dispatched early in `main()` (guarded by `HAVE_QC`; toggle `QC=0 make` or CMake `-DENABLE_QC=OFF`); no tensor/mesh deps; ships in native and WASM builds. Label convention: `0` = non-brain (excluded); `--csf`/`--wm` are disjoint nonzero label sets; every other non-zero label is GM. Validates exact finite integer labels (seg read as float64 to avoid >2^24 aliasing), rejects non-3D/oversized inputs, and requires dims + spatial transforms to match. **Only air-free metrics are computed** (backgrounds are masked to 0): **CJV, cnr_noair, snr_{csf,wm,gm,total}, wm2max, efc_brain, icvs_{csf,gm,wm} + vol_*_mm3, summary_{tissue}_{mean,stdv,median,mad,p05,p95,k,n}**; omits SNR-Dietrich, FBER, Qi1/Qi2 (air), INU (needs N4), rPVE (needs soft pvms). **Gotchas:** (1) hard-segmentation variant, not a numerical MRIQC clone — unrounded intensities, NumPy-linear quantiles, project 6-connected erosion, vs MRIQC's rounded/weighted soft-PVM path. Formula plumbing follows MRIQC: CJV uses median + **normalised MAD** (`median(|x−med|)/0.6744897501960817`), CNR/SNR use median + population **stdv**, SNR applies `sqrt(n/(n−1))`. (2) `cnr_noair` drops MRIQC's `σ_bg²` → biased high, not comparable to MRIQC norms (renamed to flag it). (3) `efc_brain` excludes zero voxels from its support. (4) Erosion affects intensity statistics only and falls back to the raw mask below `QC_MIN_VOX=100`; ICV and mm³ volumes always use RAW counts. (5) `summary_*_n` counts print as exact integers (`%ld`, not `%.6g`). (6) Comparator-free quickselect keeps the WASM hot path free of `qsort`. The stdlib-only `release_smoke.py` independently checks formulas, TSV schema, label validation, unit normalization (mm/metre/micron), large exact counts, and the erosion fallback.
- **src/GPL/ (niimath_gpl submodule)** — Optional GPL-2 SPM coregistration, compiled only with `GPL=1` (`HAVE_GPL`). A normal clone fetches zero GPL bytes; `git submodule update --init src/GPL` is required (else `$(error)`). Two chain ops via `spmcoreg_niimath.c`: `-spm_coreg <ref>` and `-spm_deface <tmpl> <mask>`. **The license split is the design point:** the GPL module computes ONLY the 6-param rigid transform (`coreg_estimate`; ported SPM `spm_coreg`, in-memory `Vol` adapter, zero file I/O), while BSD code applies it — `nii_reslice_affine` reslices and `nii_apply_deface_mask` zeros faces, both in allineate.c (so reslice/`-spm_deface` need `HAVE_ALLINEATE`). `-estimate` is pure BSD header math. The only GPL-aware line in the BSD tree is the gated `extern` hook in `coreFLT.c` (+ a `!HAVE_GPL` stub). **`powell.c` (GPL, SPM direction-set Powell) ≠ BSD `powell_newuoa.c` (NEWUOA) — not interchangeable; the SPM golden match depends on it.** **OpenMP gotchas:** (1) the per-eval joint-histogram build (`coreg_hist2_cached`, `hist2.c`) has each thread fill its own slice of a per-pass `Hist2Scratch` buffer, reduced in fixed thread order (guarded `bs->n > 100000 && nt > 1`). **Chunk by the ACTUAL team size (`omp_get_num_threads()` inside the region), never the requested `hs->nt`** — the runtime can launch fewer threads (`OMP_THREAD_LIMIT`, dynamic teams), so chunking by the requested count silently drops the unlaunched IDs' samples (garbage result). (2) histogram smoothing (`smooth_hist_into`) is a real `#pragma omp for` map → bit-identical to serial. **Do NOT parallelize the `cost_on` log2 reductions** (NMI denominator etc.): that perturbation hits the final cost directly (no smoothing damping) and shifts the optimum by >1 mm / ~1° with thread count — it breaks the SPM golden match. The one-time pre-smoothing (`conv_axis_parallel`, `smooth.c`) is bit-identical parallel. **WASM** (`GPL=1 make wasm`): needs `-s STACK_SIZE=4194304` (the cost function's 512 KB on-stack joint histogram overflows emcc's 64 KB default); worker OOM is caught by an `SC_TLOCAL` setjmp/longjmp guard in `coreg_run()`. **Packaging:** the published `@niivue/niimath` npm package ships TWO WASM builds — the default `.` export (BSD-2, includes `-allineate`/`-deface`) and a `./gpl` export (GPL-2 combined work adding `-spm_coreg`/`-spm_deface`, emitted to `niimath-gpl.js`); `js/esbuild.config.ts` skips the GPL entry gracefully if it was not built. `package.json` `license` is `BSD-2-Clause AND GPL-2.0-only`; the tarball carries `js/LICENSE`, `js/LICENSE.GPL-2.0.txt`, and `js/GPL-NOTICE.md` (the §3(b) written offer, which **hardcodes the pinned `niimath_gpl` submodule SHA — update it in lockstep whenever the `src/GPL` pin is bumped**). SPM-MATLAB golden (gitignored): nmi 0.060 mm / 0.193°, reslice r=0.99996; WASM↔desktop parity corr 0.99988.
- **tensor.c** — Eigenvalue decomposition (Jacobi method from EISPACK)
- **bw.c** — Butterworth IIR filter design (LGPL, Exstrom Laboratories)
- **bwlabel.c** — Connected component labeling (6/18/26 connectivity)
- **conform.c** — Image conforming to 1mm³ standard space
- **filter.c** — Separable resampling filters (box, triangle, B-spline, Lanczos3, Mitchell)
- **bmp.c** — PNG/BMP slice visualization with color LUTs
- **radixsort.c** — Radix sort for median operations
- **fdr.c** — False Discovery Rate correction
- **base64.c** — Base64 encoding for GIfTI export

### External/vendored libraries
- **nifti_io.c/nifti_io.h** — Consolidated NIfTI 1/2 format I/O with integrated zlib and optional zstd wrapper (public domain, ~2k lines; replaces niftilib and znzlib)
- **spng.c** — PNG encoder library (~7k lines)

## Shared-code invariant (allineate)

`allineate.c`, `allineate.h`, `powell_newuoa.c`, `coreg_fast.c`, and `coreg_fast.h` are shared **byte-identically** with the standalone project at `/Users/chris/src/allineate`, which is the canonical source of truth for the registration/coreg engine. Any change to these five files MUST be mirrored both ways and verified identical (`diff`). Everything else — `core.c`, `coreFLT.c`, the JSON `-savemat`/`-applymat` I/O, and the niimath CLI dispatch — is niimath-only. Consequences to preserve: do NOT delete the `g_last_affine`/`nii_last_affine` global (the standalone's `-savemat` reads it); the no-form (both form codes 0) policy is the single exported `al_image_xform_or_pixdim(nim, out, who)` used by `al_register`, coreg_fast, and the seeds (pixdim-centered fallback); a header-mutating seed (`-com`/`-sym`) composes its saved matrix back to the ORIGINAL frame (`M' = S_orig · inv(S_seed) · M`) so `-applymat` on the un-seeded input is correct; `al_parse_subopts(..., caps)` rejects out-of-scope options AT PARSE TIME via an `AL_CAP_*` mask (`-deface` passes `AL_CAP_TUNING|AL_CAP_FINAL`, so seed/matrix/master/fast options are rejected, not silently ignored — a privacy hazard) — this is the single choke point.

### AFNI-fidelity note (allineate coarse search)
The coarse random+grid search (`al_scalar_ransetup`, a port of AFNI's `mri_genalign_scalar_ransetup`) is deliberately AFNI-faithful, matching values verbatim: `PARAM_MAXTRIAL=29` (AFNI `3ddata.h`), `NKEEP=3*PARAM_MAXTRIAL+1`=88 and the `nrand≥NKEEP+13`=101 floor (AFNI `mri_genalign.c`), `DEFAULT_TBEST`=5 / `DEFAULT_TBEST_LPA`=17 and the `vol_src>1.3*vol_base→tbest=PARAM_MAXTRIAL` largeness branch (AFNI `3dAllineate.c`). AFNI runs this full 88/101 coarse search **regardless of `tbest`** (tbest only controls how many candidates get refined), so the coarse pass is the dominant runtime (~8.5 s of ~9.3 s on a T1→MNI-2mm case) and runs at full resolution (no downsampling). This cost IS AFNI fidelity — do not "optimize" it by shrinking `NKEEP`/`nrand` for the default case, which would silently make the coarse search sparser than AFNI. `PARAM_MAXTRIAL+2` sizes the candidate arrays; max `tfdone` is `tbest+1`=30, so 31 has headroom. Two branches (`al_shift_range_mm` padded-FOV expansion, and the largeness `tbest=29`) fire only for edge geometries (off-center subject; source FOV >1.3× base — the whole-head→brain-template/defacing direction); they lack committed CI fixtures (manual + UBSan verified only).

## Memory Management

- **Use plain `malloc`/`calloc`/`free` everywhere.** Do NOT reintroduce `_mm_malloc`/`_mm_free` or `arm_malloc.h`. The 64-byte alignment they gave was never load-bearing: the SIMD is SSE-only (128-bit) with unaligned loads (`_mm_loadu_*`), `malloc`'s 16-byte alignment already prevents cache-line splits, and these ops are memory-bandwidth-bound.
- **`nim->data` ownership invariant:** any buffer assigned to a NIfTI image's `nim->data` MUST come from plain `malloc`/`calloc`, because `nifti_image_free()` and core.c's datatype conversion release it with plain `free()`. On MSVC, `_mm_malloc` resolves to `_aligned_malloc`, which `free()` cannot release — plain `free()` on it corrupts the heap (STATUS_HEAP_CORRUPTION 0xC0000374 at teardown, after output is written). Use the named helper **`nii_calloc(count, size)`** (core.c/core.h) for `nim->data`: a two-argument, fail-closed allocator that runs `count`×`size` through the checked-multiply `nii_mul_size()` and `exit(EXIT_FAILURE)`s on overflow or NULL, so a caller can never get a short/NULL `nim->data`. Callers still pre-compute overflow-prone products with `nii_mul_size()`. (`release_smoke.py`'s conform/pval paths CI-exercise this on MSVC.)
- Kernel arrays: 4-int-per-voxel layout (offset, x, y, weight); plain `malloc`.
- NIfTI images: managed by `nifti_image_free()` from nifti_io.
- MarchingCubes.c: constructor/destructor pattern (`MarchingCubes()`/`FreeMarchingCubes()`) with `clean_all()`.

## Testing & Benchmarking

### Benchmark suite (`benchmark/`)
```bash
cd benchmark
bash benchmark.sh                    # Run niimath through all operations
bash conformance.sh                  # Compare niimath vs fslmaths output
bash close.sh                        # Test operations with allowed FP differences
bash slow_benchmark.sh               # Performance timing on larger datasets
```

### Leak detection (macOS)
LeakSanitizer is unavailable on Apple Silicon; ASan is libomp-deadlocked (so `make sanitize` omits OpenMP and runs serial). For routine diagnostics prefer `make ubsan` (retains OpenMP). For leaks, use `MallocNanoZone=0 leaks --atExit -- ./niimath ...` on a NORMAL (non-ASan) build.

### Test data
- `benchmark/In/` (69 test images), `benchmark/Ref/` (fslmaths reference), `benchmark/New/` (niimath output); `mesh/` (sample NIfTI + mesh).
- Corner-case validation lives in a separate repo: https://github.com/rordenlab/niimath_tests — compares niimath against a reference (fslmaths or a known-good niimath), so it needs a reference to diff. Use it when a change could affect numerical output beyond what `benchmark/` covers.
- Registration/defacing sample pairs (`-allineate`, `-deface`, `-spm_coreg`) live in the `register/` folder of https://github.com/niivue/niivue-demo-images (locally `/Users/chris/src/niivue-demo-images/register`).
- **`--compare` oracle gotcha:** its "equal" verdict and exit code are MAGNITUDE-based (`maxDiff > thresh`), and a non-finite difference (NaN/±inf) has no magnitude — so non-finite mismatches must be counted separately (`nHardMismatch`) or an all-NaN-vs-finite pair reads as EQUAL and exits 0. `essentiallyEqual` treats two same-sign inf and two NaN as equal but any one-sided/opposite non-finite as different. The canonical/close/thread-parity CI all rely on `--compare` as the oracle.

### CI
- AppVeyor for cross-platform builds; `gpl-build.yml` for GPL `-spm_coreg` + allineate thread parity and registration quality; `js.yml` builds the WASM bundle and runs the `bun` runtime test suite (BSD + GPL, GPL self-skips when its wasm isn't built).

## Known Remaining Issues

1. **NULL checks** — ~26 unchecked malloc/calloc calls remain in MarchingCubes.c (6), oldcubes.c (3), and quadric.c (17); meshify.c is largely fixed.
2. **Build feature-list drift** — the feature source/define inventory is duplicated across `src/Makefile` (all/static/debug/verbose/sanitize/wasm), `src/CMakeLists.txt`, `src/notarize.sh`, and `SuperBuild/SuperBuild.cmake` (which forwards each feature's option toggle to `src/CMakeLists.txt` — a fourth place to touch: a new feature needs its `ENABLE_*` option declared+forwarded here, and `OPENMP_XCODE` must default to `${USE_OPENMP}` so the Apple universal build's explicit `-DOPENMP_XCODE=OFF` still disables OpenMP). Adding a feature requires touching all of them. A shared/generated source-list fragment would prevent release mismatches.
3. **Whole-program `-ffast-math` (all build systems).** Make, CMake (non-MSVC `CMAKE_C_FLAGS`), notarize.sh, and wasm compile the entire program `-ffast-math -fno-finite-math-only`; MSVC stays strict. Deliberate: registration must round like the standalone allineate/fslmaths goldens (the shared `nifti_io` mat44 helpers run during the fit); scoping it to the registration TUs made CMake/notarize/wasm diverge numerically from Make. This also covers core ops, GPL `-spm_coreg`, and dtifit/tensor — re-validate the SPM-MATLAB golden and dtifit/tensor if you touch those. `-fno-finite-math-only` preserves NaN/Inf detection; fast-math-sensitive finiteness checks use magnitude guards (e.g. `v >= -DBL_MAX && v <= DBL_MAX`).
4. **Prebuilt-object flag stamp.** `allineate.o`/`powell_newuoa.o`/`coreg_fast.o` depend on a flag-signature stamp (`.al_obj_flags` = `CNAME|CFLAGS|AFLAGS`) so `make` → `make OMP=0 all` relinks cleanly with no `make -B` — otherwise toggling `OMP=`/`GPL=` leaves a stale object that fails to link on `__kmpc_*`/`omp_*`. Extend `AL_OBJ_SIG` for a new object-affecting flag outside `CFLAGS`.
5. **Output-write failures propagate.** `nifti_image_write_status` returns 0/1 (open, short header/data write, compressor/close/pclose errors incl. disk-full at close); `void nifti_image_write` remains a compatibility wrapper. `nifti_save` uses the status API; the main dispatch, pass-through, and intermediate `-save` sites propagate it; multi-output `--dtifit`/`-tensor_decomp` OR every status (`save_rc |= …`). Preserve this — a failed/partial write must report failure, not exit 0.

### macOS universal release (zstd)
The AppVeyor macOS job builds a universal binary by compiling x86_64 + arm64 slices separately and `lipo`-combining them. Homebrew ships only the runner's native arch of libzstd, so the cross slice cannot link homebrew zstd. The job builds a **universal static `libzstd.a` from source** (`-arch x86_64 -arch arm64`) and points both slices at it via `PKG_CONFIG_PATH` (also making the binary self-contained, no runtime `libzstd.dylib`). `src/CMakeLists.txt` resolves the pkg-config result to a full library path so the correct zstd links for a cross/universal build. The universal build also passes `-DOPENMP_XCODE=OFF` (single-arch libomp cannot link into a universal binary) — hence `OPENMP_XCODE` must remain a consumed, override-able option.

## Optimization Constraints

- Voxel operations are lean and SIMD-optimized; most are memory-bandwidth limited. OpenMP is **not the primary optimization target** for core ops — typical usage runs one subject per thread across many subjects. CloudFlare zlib (`CF=1`) already provides major I/O speedup.
- **Allineate reproducibility:** `-allineate` is byte-reproducible **across thread counts** (nt=1/2/4/8, verified on real T1→template) — preserved by the `mfac` thread-local re-apply rule (see powell_newuoa.c). It is **NOT bit-reproducible across builds**: NEWUOA + `-ffast-math` means a codegen change can flip the optimizer to a neighboring, equally-valid optimum, shifting recovered parameters sub-voxel and resliced output by ~1% in a few edge voxels. Quality is unchanged (cost and base-correlation match to 4+ digits); only exact bytes move — do NOT treat `-allineate` output as a byte-stable golden across niimath versions. The `gpl-build.yml` CI gates registration by QUALITY (correlation floor + improvement-over-baseline via `reg_quality.py`, which fails closed on empty/constant/non-finite output), not a byte diff. `-deface`/`-spm_deface` reslice via `nii_reslice_affine` (not the fused path) and stay byte-stable.

### WASM performance pitfalls (two hard-won rules)
1. **No `qsort`+comparator on large or per-voxel arrays in WASM-targeted code.** emscripten/musl `qsort` calls the comparator via `call_indirect` per comparison — ~100× slower than native, enough to make the WASM build appear to *hang* on large images. Use a comparator-free routine: **quickselect** for a single value (`al_select_rank`, `select_kth_flt`, `uf_select`), **insertion sort** for genuinely small fixed N only (NOT a per-voxel neighborhood — that is O(n²) per voxel), or **heapsort** for a full order on large arrays (`heapsort_sortIdx`). The only remaining `qsort` is `qsort_floatint` on ~46-element candidate lists (negligible).
2. **WASM allineate MUST be compiled `-ffast-math` (matching native)**, or the cost-function reduction rounds differently and NEWUOA/Powell converges to a worse optimum (`-deface` masked 37% vs 60% == native). The `wasm:` target and its prebuilt `al_wasm.o`/`pn_wasm.o`/`cf_wasm.o` use whole-program `-ffast-math -fno-finite-math-only`. (`-msimd128` adds only ~3% — not worth the golden-match risk.)

## Code Conventions
- C99 with extensive `#ifdef` for conditional compilation.
- Template pattern: coreFLT.c compiled as both float32 and float64 via macro inclusion.
- Function naming: `nifti_*` for NIfTI operations, `nii_*` for internal helpers.
- Explicit SIMD (Intel intrinsics via `immintrin.h`) is compiled only on `__x86_64__`; ARM/WASM use scalar fallbacks, which clang `-O3` auto-vectorizes to NEON just as fast (these ops are memory-bandwidth bound). No `sse2neon.h` shim (it gave no benefit on Apple Silicon).
- Error returns: `EXIT_SUCCESS`/`EXIT_FAILURE` from stdlib.

## Documentation style
- Do not add artificial end-of-line characters to Markdown/text prose. Let the renderer word-wrap — one paragraph (or list item) is one physical line; do not hard-wrap at a fixed column.
- Use newlines only for genuine structure: between paragraphs, headings, list items, table rows, and code-fence boundaries, and for new lines of code inside fences.
- When editing an existing `.md` file that has hard-wrapped paragraphs, unwrap them to one line per paragraph.

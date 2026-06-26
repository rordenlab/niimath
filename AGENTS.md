# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

niimath is an open-source clone of FSL's `fslmaths` — a general-purpose NIfTI image calculator for neuroimaging. It extends fslmaths with mesh generation features (nii2mesh), additional filters, defacing (`-deface`), zstd compression (.nii.zst), and cross-platform support (Linux, macOS, Windows, WebAssembly).

**Repository:** `rordenlab/niimath` (BSD-2-Clause license)

## Building

### Quick build (Makefile, recommended for development)
```bash
cd src && make          # Standard optimized build (OpenMP enabled by default)
make debug              # Debug build (-g, no optimization)
make sanitize           # AddressSanitizer build (leak/overflow detection)
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
GPL=1 make              # Enable optional GPL -spmcoreg (needs git submodule update --init src/GPL)
ZSTD=0 make             # Disable zstd (.nii.zst) compression support
OMP=0 make              # Disable OpenMP (single-threaded build)
CF=1 make               # CloudFlare accelerated zlib
MC=1 make               # Use new MarchingCubes algorithm (handles ambiguities)
STB=1 make              # Use STB image library instead of libspng for bitmaps
make wasm               # Emscripten/WebAssembly target
```

### CMake build
```bash
mkdir build && cd build && cmake .. && make
cmake -DUSE_OPENMP=OFF ..    # To disable OpenMP (enabled by default)
```

### Key compile-time flags
- `NII2MESH` / `DNII2MESH` — enables mesh conversion features
- `USE_CLASSIC_CUBES` — selects oldcubes.h (classic marching cubes) vs MarchingCubes.h
- `HAVE_FORMATS` — enables GIfTI, OBJ, VTK, STL mesh output
- `HAVE_BMP` — bitmap/PNG creation
- `HAVE_BUTTERWORTH` — bandpass temporal filtering
- `HAVE_TENSOR` — tensor decomposition
- `HAVE_DTIFIT` — linear diffusion tensor fit (`--dtifit`), emulates FSL dtifit; needs `HAVE_TENSOR`
- `HAVE_CONFORM` — image conforming to standard space
- `HAVE_ALLINEATE` — affine image registration, defacing, and skull-stripping (allineate.c + powell_newuoa.c); compiled separately with -ffast-math and OpenMP
- `HAVE_GPL` — optional GPL `-spmcoreg`/`-spm_deface` (SPM rigid-body coregistration/defacing); OFF by default so builds stay BSD-2. Sources in the `niimath_gpl` submodule (`src/GPL`); `GPL=1 make` / `cmake -DENABLE_GPL=ON` / `GPL=1 make wasm` (needs `git submodule update --init src/GPL`). A `GPL=1` binary is a GPL-2 combined work (version string ends ` GPL` vs ` BSD`). See the **src/GPL architecture note** below for the license split, the WASM stack/OOM gotchas, and the BSD-minimal packaging decision
- `AL_LPC_MICHO` — enables lpc+ZZ/lpa+ZZ combined cost variant for allineate (default: pure lpc/lpa)
- `HAVE_ZSTD` — zstd compression support for .nii.zst files (auto-detected; set `FSLOUTPUTTYPE=NIFTI_ZST` to write)
- `FSLSTYLE` — FSL-compatible behavior mode

## Source Architecture

### Core computational pipeline
- **niimath.c** — CLI entry point, dispatches to `main32()` or `main64()` based on `-dt` flag
- **coreFLT.c** (6k lines) — Main computational engine, compiled twice via template pattern:
  - **core32.c** — `#define DT32` + `#include "coreFLT.c"` → float32 (SSE 4-wide SIMD on x86_64; scalar on ARM/WASM)
  - **core64.c** — includes coreFLT.c without DT32 → float64 (SSE 2-wide SIMD on x86_64; scalar on ARM/WASM)
- **core.c** — Shared utilities: datatype conversion, kernel creation, Otsu thresholding, resampling filters, NIfTI I/O helpers
- **unifize.c** — Bias field correction via `-unifize` flag (adapted from AFNI 3dUnifize, public domain)
- **allineate.c** (~3.9k lines) — Affine registration, defacing, skull-stripping. Shared identically with the standalone `allineate/` project. Cost functions via `-cost` (Hellinger default, lpc/lpa/ls; `-DAL_LPC_MICHO` adds lpc+ZZ/lpa+ZZ); DOF via `-warp`; match/output interpolation via `-interp`/`-final`. Twopass coarse-to-fine with 2x-downsampled coarse grid search and OpenMP-parallel candidate refinement using thread-local histogram/warp/workspace buffers (adapted from AFNI 3dAllineate, public domain) — this thread-local-cost-eval pattern is the template for any spm_coreg OpenMP work. Core `al_register()` is shared by `nii_allineate()`/`nii_deface()`/`-skullstrip`; `nii_reslice_affine` + `nii_apply_deface_mask` are factored out and reused by the GPL `-spmcoreg`/`-spm_deface` (interpolation lives once in BSD). **`nii_deface` rejects 4D input** (`al_dims_ok`): the face mask covers a single volume, so 4D would silently leave faces in volumes 2..N — a privacy failure. Options in the `al_opts` struct (`allineate.h`).
- **powell_newuoa.c** (~2.8k lines) — Powell's NEWUOA derivative-free optimizer (f2c translation, used by allineate). Thread-safe statics via `__thread` for parallel use.

### Mesh code (nii2mesh — unique to niimath, not in FSL)
- **meshify.c** — Main mesh pipeline: smoothing → marching cubes → vertex welding → degenerate removal → export
- **MarchingCubes.c/.h** — Newer implementation with ambiguity resolution (14 cases with subcases)
- **oldcubes.c/.h** — Classic/simpler marching cubes (~500 lines, selected with `USE_CLASSIC_CUBES`)
- **quadric.c** — Mesh simplification via quadric error metrics
- **meshtypes.h** — `vec3d` (double xyz), `vec3i` (int xyz) structs

### Supporting modules
- **dtifit.c** — `niimath --dtifit`: linear diffusion tensor fit emulating FSL `dtifit`. Self-contained TU dispatched early in `main()`, FSL-identical flags (`-k/-m/-r/-b/-o`, `-xflip 0|1|auto`). 7-param (S0 + 6 tensor) log-linear OLS (fit math from AFNI 3dDWItoDT linear path, public domain; nonlinear omitted), reuses `EIG_tsfunc` (tensor.c) for FA/MD/L*/V*. Emulates FSL's determinant-based bvec x-flip. Writes `<base>_{FA,MD,L1,L2,L3,V1,V2,V3,S0,MO,tensor}`; unsupported FSL features (`--wls`/`--kurt`) error clearly. Validated against FSL (FA r=0.998, MD/L* r≈0.999, V* \|cos\|≈0.9998; data gitignored, regenerate to re-run). **Gotchas:** tensor buffer must use NIfTI planar (volume-major) layout; validate V* with \|cos angle\| not Pearson r (eigenvector sign is arbitrary).
- **src/GPL/ (niimath_gpl submodule)** — Optional GPL-2 SPM coregistration, compiled only with `GPL=1` (`HAVE_GPL`). A normal clone fetches zero GPL bytes (`.gitmodules` = URL + pinned SHA); `git submodule update --init src/GPL` is required (else `$(error)`). Two chain ops via `spmcoreg_niimath.c`: `-spmcoreg <ref>` and `-spm_deface <tmpl> <mask>`. **The license split is the design point:** the GPL module computes ONLY the 6-param rigid transform (`coreg_estimate`; ported SPM `spm_coreg`, in-memory `Vol` adapter, zero file I/O), while BSD code applies it — `nii_reslice_affine` reslices and `nii_apply_deface_mask` zeros faces, both in allineate.c (so reslice/`-spm_deface` need `HAVE_ALLINEATE`). `-estimate` is pure BSD header math (no allineate dep). The only GPL-aware line in the BSD tree is the gated `extern` hook in `coreFLT.c` (+ a `!HAVE_GPL` stub). **`powell.c` (GPL, SPM direction-set Powell) ≠ BSD `powell_newuoa.c` (NEWUOA) — not interchangeable; the SPM match depends on it.** **OpenMP** parallelizes the two safe per-eval hot loops; Powell itself is sequential. (1) The joint-histogram build (`coreg_hist2_cached`, `hist2.c`): each thread zeroes+fills its own 512 KB slice of a **per-pass** buffer (`Hist2Scratch`, owned by `CostScratch` — allocated once per pass, no per-eval `malloc`), reduced in fixed thread order. Guarded `bs->n > 100000 && nt > 1`. **Chunk by the ACTUAL team size (`omp_get_num_threads()` inside the region), never the requested `hs->nt`** — `num_threads()` is only an upper bound and the runtime can launch fewer (`OMP_THREAD_LIMIT`, dynamic teams), so manual chunking by the requested count silently drops the unlaunched IDs' samples (an audit caught this: `OMP_THREAD_LIMIT=2 OMP_NUM_THREADS=8` → r=0.70, garbage). The reduction runs over the actual size so no stale slice is summed. (2) The histogram smoothing (`smooth_hist_into`, `cost.c`): a map (each output cell is an independent fixed-order convolution), parallelized with a real `#pragma omp for` (immune to the team-size bug by construction) → **bit-identical** to serial. All `#pragma omp` no-op without `-fopenmp`. **Determinism:** the histogram reduction reorders float adds, so it is deterministic *for a given actual team size* but not bit-identical to serial; the result stays within the SPM golden tolerance (the histogram feeds a smoothed cost) — empirically byte-identical across thread counts on the test data, but only *guaranteed* equivalent within tolerance, not bit-exact (smoothing is linear; it propagates rather than erases ULP differences). **Do NOT parallelize the `cost_on` log2 reductions** (the NMI denominator etc.): unlike the histogram, that perturbation hits the final cost *directly* (no smoothing in between) and shifts the optimum by >1 mm / ~1° depending on thread count — it breaks the SPM golden match. Measured **~4.4× at 8 threads** (memory-bandwidth-bound, saturates ~8). Thread parity (reproducibility at fixed `nt`; equivalence across `nt`; the `OMP_THREAD_LIMIT<nt` regression) is covered by the `gpl-build.yml` "OpenMP thread parity" CI step on an anisotropic multi-lobe 64³ fixture that clears the 100000-sample threshold. **WASM** (`GPL=1 make wasm`): needs `-s STACK_SIZE=4194304` (the cost function's 512 KB on-stack joint histogram overflows emcc's 64 KB default); OOM in the long-lived worker is caught by an `SC_TLOCAL` setjmp/longjmp guard in `coreg_run()` (`xalloc.h` allocators jump there instead of `exit(1)`). **Packaging:** `@niivue/niimath` stays BSD-minimal (allineate gated on GPL via `$(WASM_AL)`); the GPL toolset is vendored into the GPL-2 `niivue/deface` app (live at niivue.github.io/deface), not a separate npm package. **Validation vs SPM-MATLAB (golden, gitignored):** nmi 0.060 mm / 0.193°, reslice r=0.99996; WASM↔desktop parity corr 0.99988 through the full chain. Build wiring duplicated across Makefile/CMake/SuperBuild (drift caveat).
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

## Memory Management Patterns

- Heap allocations use `calloc()`/`malloc()` paired with `free()`
- SIMD-aligned allocations: `_mm_malloc(size, 64)` / `_mm_free()`
- Kernel arrays: 4-int-per-voxel layout (offset, x, y, weight) with 64-byte alignment
- NIfTI images: managed by `nifti_image_free()` from nifti_io
- MarchingCubes.c: constructor/destructor pattern (`MarchingCubes()`/`FreeMarchingCubes()`) with `clean_all()` helper

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
```bash
/usr/bin/time -l ./niimath test4D -add 0 tst    # Peak memory + page faults
```

### AddressSanitizer build
```bash
cd src && make sanitize    # Builds with -fsanitize=address
```

### Test data
- `benchmark/In/` — 69 test images (3D anatomical, 4D fMRI, tensors, orientations)
- `benchmark/Ref/` — fslmaths reference output
- `benchmark/New/` — niimath output for comparison
- `mesh/` — Sample NIfTI (`bet.nii.gz`) and mesh output (`mesh.gii`)

### CI
- AppVeyor CI for cross-platform builds (see badge in README.md)

## Known Remaining Issues

1. **NULL checks** — meshify.c is largely fixed; ~26 unchecked malloc/calloc calls remain in MarchingCubes.c (6), oldcubes.c (3), and quadric.c (17)
2. **Build feature-list drift** — the feature source/define inventory is duplicated across `src/Makefile` (all/static/debug/verbose/sanitize/wasm), `src/CMakeLists.txt`, and `src/notarize.sh`; adding a feature requires touching all of them (dtifit hit this). A shared/generated source-list fragment would prevent release mismatches.
3. **`nifti_save` returns 0 unconditionally** (`core.c`) — callers cannot detect write failures. Commands writing many outputs (e.g. `--dtifit`'s 11 files) can report success on a failed write. Project-wide API improvement.

### macOS universal release (zstd)

The AppVeyor macOS job builds a universal binary by compiling x86_64 and arm64 slices separately and `lipo`-combining them. Homebrew only ships the runner's native arch of libzstd, so the cross-compiled slice cannot link homebrew zstd. The job therefore builds a **universal static `libzstd.a` from source** (`-arch x86_64 -arch arm64`) and points both slices at it via `PKG_CONFIG_PATH`; this also makes the released binary self-contained (no runtime `libzstd.dylib`). `src/CMakeLists.txt` resolves the pkg-config result to a full library path so the correct (cross/universal/non-default-prefix) zstd links — without this, zstd fails to link from `/opt/homebrew` or for a cross build.

## Optimization Constraints

- Voxel operations are already lean and SIMD-optimized; most are memory-bandwidth limited
- OpenMP is **not the primary optimization target** for core ops — typical usage runs one subject per thread across many subjects
- CloudFlare zlib (`CF=1`) already provides major I/O speedup
- Allineate: NEWUOA optimizer is inherently sequential (trust-region). The cost function hot path is 3D cubic interpolation (64 scattered array lookups per voxel), memory-bandwidth-limited. 2x downsampling for coarse grid search and interior-loop cubic optimization (no bounds checking for safe voxels) are the practical workarounds. Candidate refinement is parallelized across OpenMP threads.
- Allineate/deface/skullstrip compiled separately with `-ffast-math` and OpenMP, scoped to avoid affecting other translation units. `isfinite()` replaced with magnitude guard for -ffast-math safety. Core files (`allineate.c`, `allineate.h`, `powell_newuoa.c`) are shared identically with the standalone `allineate/` project.

## Code Conventions
- C99 with extensive use of `#ifdef` for conditional compilation
- Template pattern: coreFLT.c compiled as both float32 and float64 via macro inclusion
- Function naming: `nifti_*` for NIfTI operations, `nii_*` for internal helpers
- Explicit SIMD (Intel intrinsics via `immintrin.h`) is compiled only on `__x86_64__`; ARM/WASM use the scalar fallbacks, which clang `-O3` auto-vectorizes to NEON just as fast (these ops are memory-bandwidth bound). The `sse2neon.h` shim was removed after benchmarks showed it gave no benefit on Apple Silicon (bit-identical output, conformance suite passes).
- Error returns: `EXIT_SUCCESS`/`EXIT_FAILURE` from stdlib

## Documentation style
- Do not add artificial end-of-line characters to Markdown/text prose. Let the editor/renderer word-wrap to the window — one paragraph (or list item) is one physical line; do not hard-wrap paragraphs at a fixed column.
- Use newlines only for genuine structure: between paragraphs, headings, list items, table rows, and code-fence boundaries, and for new lines of code inside fences.
- When editing an existing `.md` file that has hard-wrapped paragraphs, unwrap them to one line per paragraph.

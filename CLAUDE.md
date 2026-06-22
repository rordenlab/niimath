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
- `AL_LPC_MICHO` — enables lpc+ZZ/lpa+ZZ combined cost variant for allineate (default: pure lpc/lpa)
- `HAVE_ZSTD` — zstd compression support for .nii.zst files (auto-detected; set `FSLOUTPUTTYPE=NIFTI_ZST` to write)
- `FSLSTYLE` — FSL-compatible behavior mode

## Source Architecture

### Core computational pipeline
- **niimath.c** — CLI entry point, dispatches to `main32()` or `main64()` based on `-dt` flag
- **coreFLT.c** (6k lines) — Main computational engine, compiled twice via template pattern:
  - **core32.c** — `#define DT32` + `#include "coreFLT.c"` → float32 with SSE 4-wide SIMD
  - **core64.c** — includes coreFLT.c without DT32 → float64 with SSE 2-wide SIMD
- **core.c** — Shared utilities: datatype conversion, kernel creation, Otsu thresholding, resampling filters, NIfTI I/O helpers
- **unifize.c** — Bias field correction via `-unifize` flag (adapted from AFNI 3dUnifize, public domain)
- **allineate.c** (~3.9k lines) — Affine image registration, defacing, and skull-stripping. Shared identically with the standalone `allineate/` project. Supports Hellinger (default), lpc, lpa, and Pearson (ls) cost functions via `-cost`; compile with `-DAL_LPC_MICHO` for lpc+ZZ/lpa+ZZ combined cost variant. Variable DOF via `-warp` (sho/3, shr/6, srs/9, aff/12; default: aff). Matching interpolation via `-interp` (NN, linear, cubic; default: linear). Output interpolation via `-final` or `-nearest`/`-linear`/`-cubic` (default: cubic for allineate, linear for deface/skullstrip via `AL_INTERP_DEFAULT`). `-cmass`/`-nocmass` control center-of-mass initial alignment; `-source_automask` fills outside source brain mask with noise for robust cross-modal registration. CLEQWD edge-bin histogram mode (from AFNI's `clipate`/`THD_cliplevel`). TOHD blok type (truncated octahedron, ~555 voxels/blok) for local Pearson correlation. 2x downsampling for coarse grid search. Twopass coarse-to-fine optimization with parallel candidate refinement (adapted from AFNI 3dAllineate, public domain). Core registration in `al_register()` helper, used by `nii_allineate()`, `nii_deface()`, and `-skullstrip`. Options defined in `al_opts` struct in `allineate.h` (cost, cmass, source_automask, interp, final_interp, warp, skullstrip). Reports wall-clock time and thread count on completion. OpenMP parallelization of coarse search and candidate refinement. Thread-local histogram, warp matrix, and workspace buffers.
- **powell_newuoa.c** (~2.8k lines) — Powell's NEWUOA derivative-free optimizer (f2c translation, used by allineate). Thread-safe statics via `__thread` for parallel use.

### Mesh code (nii2mesh — unique to niimath, not in FSL)
- **meshify.c** — Main mesh pipeline: smoothing → marching cubes → vertex welding → degenerate removal → export
- **MarchingCubes.c/.h** — Newer implementation with ambiguity resolution (14 cases with subcases)
- **oldcubes.c/.h** — Classic/simpler marching cubes (~500 lines, selected with `USE_CLASSIC_CUBES`)
- **quadric.c** — Mesh simplification via quadric error metrics
- **meshtypes.h** — `vec3d` (double xyz), `vec3i` (int xyz) structs

### Supporting modules
- **dtifit.c** — `niimath --dtifit`: linear diffusion tensor fit emulating FSL's `dtifit`. Self-contained TU dispatched early in `niimath.c:main()` (like the `.mz3` mesh path), parsing FSL-identical flags (`-k/-m/-r/-b/-o`, plus `-xflip 0|1|auto`). Fits the 6-volume tensor by 7-parameter (S0 + 6 tensor terms) log-linear OLS (`D = (AᵀA)⁻¹Aᵀ·logS`; fit math from AFNI 3dDWItoDT linear path, public domain — nonlinear path intentionally omitted), then reuses `EIG_tsfunc` (tensor.c) for FA/MD/L*/V*. Emulates FSL's determinant-based bvec x-flip (flips gradient X when voxel→world det>0). Writes `<base>_{FA,MD,L1,L2,L3,V1,V2,V3,S0,MO,tensor}`. Unsupported FSL features (`--wls`, `--kurt`, etc.) error out clearly. **Validation (against FSL `dtifit` reference, in brain mask):** FA r=0.998, MD r=0.9998, L1–L3 r≈0.999, S0 r=1.0, MO r=0.998, V1/V2/V3 \|cos\|≈0.9998, tensor r=1.0. bvec determinant x-flip confirmed FSL-correct on real LPS (det>0, auto-flips) and LAS (det<0, no flip) images: both yield identical world-space fiber orientations (\|cos\|=1.0); a flip-disabled negative control drops to 0.55. Validation data and scripts are kept out of the repo (gitignored), so re-running requires regenerating the FSL reference. Notes on gotchas: tensor buffer must use NIfTI planar (volume-major) layout; eigenvector sign is arbitrary, so validate V* with \|cos angle\|, not Pearson r.
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
- **sse2neon.h** — SSE→NEON SIMD translation for ARM (~228k)

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
3. **`src/windows.bat` is stale/non-functional** — it references the removed `niftilib/` and `znzlib/` (replaced by `nifti_io.c`) and lacks most feature defines; it cannot build. The maintained Windows path is CMake (AppVeyor). Remove or rewrite it.
4. **`nifti_save` returns 0 unconditionally** (`core.c`) — callers cannot detect write failures. Commands writing many outputs (e.g. `--dtifit`'s 11 files) can report success on a failed write. Project-wide API improvement.

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
- SIMD code has scalar fallbacks gated on `__x86_64__`, `__aarch64__`, or `myDisableSSE`
- Error returns: `EXIT_SUCCESS`/`EXIT_FAILURE` from stdlib

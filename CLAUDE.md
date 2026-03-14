# CLAUDE.md - niimath Project Guide

## Project Overview

niimath is an open-source clone of FSL's `fslmaths` — a general-purpose NIfTI image calculator for neuroimaging. It extends fslmaths with mesh generation features (nii2mesh), additional filters, and cross-platform support (Linux, macOS, Windows, WebAssembly).

**Repository:** `rordenlab/niimath` (BSD-2-Clause license)

## Building

### Quick build (Makefile, recommended for development)
```bash
cd src && make          # Standard optimized build (includes allineate with OpenMP)
make debug              # Debug build (-g, no optimization)
make sanitize           # AddressSanitizer build (leak/overflow detection)
make verbose            # All warnings enabled
make static             # Static binary
```

### Prerequisites
- **macOS**: `brew install libomp` (required for allineate OpenMP parallelization)
- **Linux**: OpenMP support is built-in with gcc (`-fopenmp`)

### Build variants
```bash
make tiny               # Minimal terminal build (emulates WASM constraints)
make nano               # Without mesh (nii2mesh) functions
MESH=0 make             # Disable mesh support
AL=0 make               # Disable allineate registration
OMP=1 make              # Enable OpenMP for core ops (requires gcc-9 on macOS)
CF=1 make               # CloudFlare accelerated zlib
MC=1 make               # Use new MarchingCubes algorithm (handles ambiguities)
STB=1 make              # Use STB image library instead of libspng for bitmaps
make wasm               # Emscripten/WebAssembly target
```

### CMake build
```bash
mkdir build && cd build && cmake .. && make
```

### Key compile-time flags
- `NII2MESH` / `DNII2MESH` — enables mesh conversion features
- `USE_CLASSIC_CUBES` — selects oldcubes.h (classic marching cubes) vs MarchingCubes.h
- `HAVE_FORMATS` — enables GIfTI, OBJ, VTK, STL mesh output
- `HAVE_BMP` — bitmap/PNG creation
- `HAVE_BUTTERWORTH` — bandpass temporal filtering
- `HAVE_TENSOR` — tensor decomposition
- `HAVE_CONFORM` — image conforming to standard space
- `HAVE_ALLINEATE` — affine image registration (allineate.c + powell_newuoa.c); compiled separately with -ffast-math and OpenMP
- `FSLSTYLE` — FSL-compatible behavior mode

## Source Architecture

### Core computational pipeline
- **niimath.c** — CLI entry point, dispatches to `main32()` or `main64()` based on `-dt` flag
- **coreFLT.c** (6k lines) — Main computational engine, compiled twice via template pattern:
  - **core32.c** — `#define DT32` + `#include "coreFLT.c"` → float32 with SSE 4-wide SIMD
  - **core64.c** — includes coreFLT.c without DT32 → float64 with SSE 2-wide SIMD
- **core.c** — Shared utilities: datatype conversion, kernel creation, Otsu thresholding, resampling filters, NIfTI I/O helpers
- **unifize.c** — Bias field correction via `-unifize` flag (adapted from AFNI 3dUnifize, public domain)
- **allineate.c** (~2.7k lines) — Affine (12 DOF) image registration with lpc+ZZ cost function, twopass coarse-to-fine optimization (adapted from AFNI 3dAllineate, public domain). OpenMP parallelization of coarse search and candidate refinement. Thread-local histogram, warp matrix, and workspace buffers.
- **powell_newuoa.c** (~2.8k lines) — Powell's NEWUOA derivative-free optimizer (f2c translation, used by allineate). Thread-safe statics via `__thread` for parallel use.

### Mesh code (nii2mesh — unique to niimath, not in FSL)
- **meshify.c** — Main mesh pipeline: smoothing → marching cubes → vertex welding → degenerate removal → export
- **MarchingCubes.c/.h** — Newer implementation with ambiguity resolution (14 cases with subcases)
- **oldcubes.c/.h** — Classic/simpler marching cubes (~500 lines, selected with `USE_CLASSIC_CUBES`)
- **quadric.c** — Mesh simplification via quadric error metrics
- **meshtypes.h** — `vec3d` (double xyz), `vec3i` (int xyz) structs

### Supporting modules
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
- **nifti_io.c/nifti_io.h** — Consolidated NIfTI 1/2 format I/O with integrated zlib wrapper (public domain, ~1.9k lines; replaces niftilib and znzlib)
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

## Current Development Focus

### Phase 1: Bug fixing and major optimizations
**Completed fixes:**
- **Endianness bug** — `&littleEndianPlatform` (address, always true) → `littleEndianPlatform()` (12 call sites in meshify.c + MarchingCubes.c); also fixed `#ifdef` structural bug in save_mz3
- **Double free** in `nifti_image_change_datatype()` DT_INT32 section (core.c)
- **sizeof mismatch** — `sizeof(float)` vs `sizeof(flt)` in memcpy for 64-bit mode (4 locations in coreFLT.c)
- **Memory leaks** — `vxs2` never freed in loop (coreFLT.c), MCB struct never freed in `FreeMarchingCubes()` (MarchingCubes.c)
- **Unsafe realloc** — data loss on malloc failure in `add_triangle()`/`test_vertex_addition()` (MarchingCubes.c) and oldcubes.c
- **Division by zero** — `add_x/y/z/c_vertex()` (MarchingCubes.c), `vx()` in EDT (coreFLT.c)
- **Uninitialized fields** — TTriangle structs in quadric.c (switched to calloc)
- **Dead code** — unreachable `return 1.0` in `calculate_error()` (quadric.c), duplicate `loopi` macro
- **Allineate registration** — added affine image registration via `-allineate` flag (adapted from AFNI 3dAllineate, public domain; uses lpc+ZZ cost, twopass optimization, NEWUOA optimizer)
- **Allineate optimization** — 2.0-2.3x wall-clock speedup (23-28s → 11-12s on test data). See validate/README.md for detailed benchmarks.

**Remaining priority areas:**
1. **NULL checks** — ~22 malloc/calloc calls in mesh code lack NULL checks (pre-existing)
2. **File handle leaks** — save_mz3 doesn't close fp/fgz if malloc fails after fopen
3. **Buffer overflows** — fixed-size buffers in save_mesh() (768-char filename limit)

### Allineate optimization details

**Completed optimizations (2.0-2.3x speedup, quality maintained or improved):**
- **OpenMP parallel coarse search**: Pre-generate all grid+random parameter sets, evaluate 7K+ reflections in parallel. Thread-safe via `__thread` histogram statics, warp matrix, and powell_newuoa.c statics.
- **OpenMP parallel candidate refinement**: 46 candidates refined with independent NEWUOA runs in parallel.
- **-ffast-math** for allineate.c + powell_newuoa.c (compiled as separate .o files, scoped to avoid affecting other translation units). `isfinite()` replaced with -ffast-math-safe magnitude guard.
- **Thread-local buffer reuse**: avm, wpar, imw/jmw/kmw buffers allocated once per thread and reused across cost function evaluations (eliminates ~200K+ malloc/free per registration).
- **Fused incremental warp+interpolation kernel**: For the all-voxels case (`im_ar == NULL`), replaces separate warp→interpolate pipeline with a single (k,j,i) loop. The affine warp becomes separable: per-row base computed once, then 3 adds per voxel instead of 12 multiplies. Eliminates intermediate coordinate arrays (~3MB/eval), chunking overhead, and mod/div coordinate generation. With OpenMP on k-slices for fine-pass parallelism.
- **Two-stage fine pass**: Stage 1 optimizes at 1/4 resolution with linear interpolation (fast exploration), stage 2 refines at full resolution with cubic. Reduces total cubic evaluations by ~50% and improves quality (better exploration before cubic refinement).
- **Thread count cap**: Limited to 8 (configurable via `OMP_NUM_THREADS`). Peak memory ~115MB (vs 50MB baseline, mainly OpenMP thread stacks).

**Evaluated but not implemented (diminishing returns):**
- **Apple Accelerate (vDSP/vImage/BLAS)**: No vDSP/vImage support for 3D volumetric interpolation (the core bottleneck at 45% of cost evaluation time). cblas_sgemm overhead exceeds benefit for 3×4 affine matrices. vDSP_vclip could replace cubic clipping loop but saves <0.1%. The fused incremental approach already reduces warp cost below what BLAS could offer.
- **NEON SIMD for cubic interpolation**: Could vectorize Lagrange polynomial evaluation (4 points at a time) for ~1.3x on cubic = ~5-10% overall. Worthwhile but diminishing returns. Main challenge is NEON's lack of native gather; relies on spatial coherence of neighboring warped coordinates for cache sharing.
- **Metal GPU compute**: The entire warp+interpolate+correlate pipeline (~90M FLOPs/eval) could run on Apple M-series GPU in <1ms per evaluation. This is the largest untapped opportunity (potentially 10x+ for fine pass) but requires significant integration work and introduces platform dependency. Best pursued as a dedicated GPU backend.
- **Progressive resolution within NEWUOA**: Expose trust region radius via callback to dynamically adjust npt_match during optimization. Would avoid full-resolution evaluations early when NEWUOA is taking large steps. Requires modifying the f2c-translated optimizer internals.
- **AMX co-processor**: Only accessible through Accelerate functions; not useful for the scatter-gather memory access patterns in interpolation.

**Architectural constraints on further optimization:**
- NEWUOA optimizer is inherently sequential (trust-region method: each iteration depends on previous). Cannot parallelize its iterations. The two-stage approach is the practical workaround.
- The cost function hot path is 3D cubic interpolation (64 scattered array lookups per voxel). This is memory-bandwidth-limited on modern CPUs. SIMD helps with the arithmetic but not the memory access pattern.
- Histogram scatter-add (for MI/NMI metrics) is inherently sequential per element due to bin collisions. Not a bottleneck at current bin counts (~97² bins).

### General optimization notes
- Voxel operations are already lean and SIMD-optimized; most are memory-bandwidth limited
- OpenMP is available but **not the primary optimization target** for core ops — typical usage runs one subject per thread across many subjects
- CloudFlare zlib (`CF=1`) already provides major I/O speedup

## Code Conventions
- C99 with extensive use of `#ifdef` for conditional compilation
- Template pattern: coreFLT.c compiled as both float32 and float64 via macro inclusion
- Function naming: `nifti_*` for NIfTI operations, `nii_*` for internal helpers
- SIMD code has scalar fallbacks gated on `__x86_64__`, `__aarch64__`, or `myDisableSSE`
- Error returns: `EXIT_SUCCESS`/`EXIT_FAILURE` from stdlib

# CLAUDE.md - niimath Project Guide

## Project Overview

niimath is an open-source clone of FSL's `fslmaths` — a general-purpose NIfTI image calculator for neuroimaging. It extends fslmaths with mesh generation features (nii2mesh), additional filters, and cross-platform support (Linux, macOS, Windows, WebAssembly).

**Repository:** `rordenlab/niimath` (BSD-2-Clause license)

## Building

### Quick build (Makefile, recommended for development)
```bash
cd src && make          # Standard optimized build
make debug              # Debug build (-g, no optimization)
make sanitize           # AddressSanitizer build (leak/overflow detection)
make verbose            # All warnings enabled
make static             # Static binary
```

### Build variants
```bash
make tiny               # Minimal terminal build (emulates WASM constraints)
make nano               # Without mesh (nii2mesh) functions
MESH=0 make             # Disable mesh support
OMP=1 make              # Enable OpenMP (requires gcc-9 on macOS)
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
- `HAVE_ALLINEATE` — affine image registration (allineate.c + powell_newuoa.c)
- `FSLSTYLE` — FSL-compatible behavior mode

## Source Architecture

### Core computational pipeline
- **niimath.c** — CLI entry point, dispatches to `main32()` or `main64()` based on `-dt` flag
- **coreFLT.c** (6k lines) — Main computational engine, compiled twice via template pattern:
  - **core32.c** — `#define DT32` + `#include "coreFLT.c"` → float32 with SSE 4-wide SIMD
  - **core64.c** — includes coreFLT.c without DT32 → float64 with SSE 2-wide SIMD
- **core.c** — Shared utilities: datatype conversion, kernel creation, Otsu thresholding, resampling filters, NIfTI I/O helpers
- **unifize.c** — Bias field correction via `-unifize` flag (adapted from AFNI 3dUnifize, public domain)
- **allineate.c** (~2.5k lines) — Affine (12 DOF) image registration with lpc+ZZ cost function, twopass coarse-to-fine optimization (adapted from AFNI 3dAllineate, public domain)
- **powell_newuoa.c** (~2.8k lines) — Powell's NEWUOA derivative-free optimizer (f2c translation, used by allineate)

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

**Remaining priority areas:**
1. **NULL checks** — ~22 malloc/calloc calls in mesh code lack NULL checks (pre-existing)
2. **File handle leaks** — save_mz3 doesn't close fp/fgz if malloc fails after fopen
3. **Buffer overflows** — fixed-size buffers in save_mesh() (768-char filename limit)

### Optimization notes
- Voxel operations are already lean and SIMD-optimized; most are memory-bandwidth limited
- OpenMP is available but **not the primary optimization target** — typical usage runs one subject per thread across many subjects
- Focus optimization effort on operations with large algorithmic gains, not micro-optimization
- CloudFlare zlib (`CF=1`) already provides major I/O speedup

## Code Conventions
- C99 with extensive use of `#ifdef` for conditional compilation
- Template pattern: coreFLT.c compiled as both float32 and float64 via macro inclusion
- Function naming: `nifti_*` for NIfTI operations, `nii_*` for internal helpers
- SIMD code has scalar fallbacks gated on `__x86_64__`, `__aarch64__`, or `myDisableSSE`
- Error returns: `EXIT_SUCCESS`/`EXIT_FAILURE` from stdlib

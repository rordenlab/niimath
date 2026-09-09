# niimath

[![Build status](https://ci.appveyor.com/api/projects/status/7o0xp2fgbhadkgn1?svg=true)](https://ci.appveyor.com/project/neurolabusc/niimath)

## About

niimath is an open-source clone of [fslmaths](https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro3/index.html), the general image calculator of FSL. It accepts the same commands and gives equivalent results. It also adds operations that fslmaths does not have, such as mesh generation, affine and nonlinear registration, defacing, motion correction, slice-time correction, phase unwrapping and distortion correction.

fslmaths is one of the foundations of the FSL pipelines, such as [FEAT](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT), and many other tools call it. This popularity shows that it fills an important role. Scientists are often encouraged to find new solutions, but replication has value too. niimath exists for these reasons:

1. fslmaths is free of charge, but it is not [open source](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence). This limits its use in other projects, in particular commercial ones.
2. An open-source license lets niimath use libraries that the FSL team cannot. The release builds use [zlib-ng](https://github.com/zlib-ng/zlib-ng), an accelerated zlib that is much faster than the public domain zlib used by fslmaths. zlib-ng is also fast on arm64 and correct on Windows. We previously helped update the [CloudFlare zlib](https://github.com/cloudflare/zlib/pull/19) fork, which lets recent FSL releases use an accelerated library and speeds up every FSL tool.
3. niimath has minimal dependencies. This makes it easy to distribute, compile and develop. It compiles for macOS, Linux and Windows. FSL cannot target Windows.
4. niimath was designed from the start for optional parallel processing, with OpenMP and the CloudFlare-enhanced [pigz](https://github.com/madler/pigz).
5. Most programs grow organically as needs arise. A clone starts with a full specification, which permits optimization. niimath uses explicit single and double precision pipelines, so the compiler can use the SIMD instructions that every x86_64 CPU provides but that high-level code rarely exploits. Modern compilers make these operations limited by memory bandwidth, so no [hand tuning](https://github.com/neurolabusc/simd) is needed.
6. A robust regression test set has found edge cases where fslmaths gives anomalous or unexpected answers. See [Compatibility with fslmaths](#compatibility-with-fslmaths). This feedback benefits fslmaths too.
7. The code is fully reverse engineered, but the FSL team allowed us to copy their error messages and help text. This gives true plug-in compatibility. They also supplied pseudo code for poorly documented routines, so the community can understand the actual algorithms.
8. niimath is an open-source base for features that fill gaps in FSL, such as `-unsharp`, `-sobel` and `-resize`. Bob Cox once gave permission to use code from [AFNI's](https://afni.nimh.nih.gov) 3dTshift and 3dBandpass tools, but niimath ships no code from either. The slice-timing command that niimath ships, `-stc`, is a clean-room BSD-2 implementation with no AFNI code, written when `3dTshift.c` and its FFT were GPL-2. Temporal filtering is `-bptf` and `-bptfm`, which are fslmaths-compatible. MCW relicensed its 1994-2000 AFNI code to CC BY 4.0 on 12 May 2026, so that bar is gone. CC BY still asks for attribution and a statement of changes, which original code does not, so `-stc` stays as it is. See [License](#license).

There is one reason to use fslmaths instead of niimath. niimath is newer and less tested, so unknown corner cases may give poor results. fslmaths has been in use for years. In the few cases where fslmaths differs from its own documentation (described below), you can argue that its result is the `correct` one, because it agrees with itself. Other tools may have been built to expect that behavior, such as the loss of high frequency signal, and may perform worse when given the documented result.

## Installation

Choose one of these methods:

- (Recommended) Download the latest compiled release from the [GitHub releases page](https://github.com/rordenlab/niimath/releases).
- Install from [PyPI](https://pypi.org/project/niimath/) with `pip install niimath`.
- Download from the command line on Linux, macOS or Windows:

```
curl -fLO https://github.com/rordenlab/niimath/releases/latest/download/niimath_lnx.zip
curl -fLO https://github.com/rordenlab/niimath/releases/latest/download/niimath_macos.zip
curl -fLO https://github.com/rordenlab/niimath/releases/latest/download/niimath_win.zip
```

- (Developers) Download the source code from [GitHub](https://github.com/rordenlab/niimath) and build it. See [Compilation](#compilation).

## Compilation

### CMake (recommended)

On Linux and macOS, build with CMake. OpenMP is enabled by default. Affine registration uses it, and the core operations can use it:

```
git clone https://github.com/rordenlab/niimath.git
cd niimath; mkdir build; cd build; cmake ..
make
```

On macOS, OpenMP needs Homebrew's libomp. Install it with `brew install libomp`. To disable OpenMP, use `cmake -DUSE_OPENMP=OFF ..`. zstd compression support is detected automatically. To enable it, install zstd with `brew install zstd` (macOS) or `apt install libzstd-dev` (Linux).

On Windows:

```
git clone https://github.com/rordenlab/niimath.git
cd niimath & mkdir build & cd build & cmake ..
cmake --build .
```

### Makefile (alternative)

Run `make` in the `src` folder. This works with Clang/LLVM and gcc on Linux and macOS:

```
git clone https://github.com/rordenlab/niimath.git
cd niimath/src
make
```

The default build includes OpenMP for all operations. On macOS this requires `brew install libomp`. Plain `make` links the system zlib. CMake release builds default to zlib-ng (`ZLIB_IMPLEMENTATION`). Other Makefile options:

```
OMP=0 make             # Disable OpenMP
make debug             # Debug build (-g, no optimization)
make ubsan             # Lightweight undefined-behavior checks (OpenMP-safe on macOS)
make sanitize          # AddressSanitizer build
AL=0 make              # Disable allineate registration
SKULLSTRIP=1 make      # Enable surface skull stripping (-skullstrip); off by default, 64-bit native only
ZSTD=0 make            # Disable zstd compression support
make wasm              # Emscripten/WebAssembly target
make wasm-wasi         # Experimental zlib-free WASI compute backend (needs Zig)
make wasm-emcc-core    # Feature-matched zlib-free Emscripten build (WASI benchmark baseline)
```

### Optional copyleft module

The `-spm_coreg` (SPM rigid-body coregistration) and `-spm_deface` (SPM-based defacing) commands live in the separate [niimath_gpl](https://github.com/rordenlab/niimath_gpl) submodule at `src/GPL`. Both are **GPL-2 or later**. They are off by default, so ordinary builds stay BSD-2-Clause. To build them, initialize the submodule and pass `GPL=1` (Makefile) or `-DENABLE_GPL=ON` (CMake). Run both commands from the repository root so the paths resolve:

```
git submodule update --init src/GPL    # or clone with --recurse-submodules
make -C src GPL=1
```

A binary built this way is a copyleft combined work, covered by SPM's own terms. Distribute it under the GNU GPL-2 or later. Its version string ends in ` GPL`. Without the module, the version string ends in ` BSD`, and `-spm_coreg` and `-spm_deface` report a clear error. The GPL module computes only the rigid transform. niimath's BSD code applies it, with the same reslicing and mask warping that `-allineate` and `-deface` use.

With OpenMP, `-spm_coreg` and `-spm_deface` parallelize one registration (the per-evaluation histogram and smoothing). Results agree across thread counts within the SPM golden tolerance. They are deterministic at a fixed thread count, but they are not guaranteed bit-identical across different thread counts, because the histogram reduction sums in a thread-count-dependent order. Control threads with `OMP_NUM_THREADS=N` or `-p N` (place `-p` before `-spm_coreg`). For batches of many subjects, run one subject per process with `OMP_NUM_THREADS=1` instead of threading each registration.

### Windows (command line)

On Windows, the CMake method above is recommended. You can also compile directly from the command line. This example omits `-DHAVE_ZLIB`, so the binary cannot read or write `.gz` files:

```
cl /Feniimath niimath.c core.c tensor.c bwlabel.c core32.c core64.c fdr.c meshify.c MarchingCubes.c quadric.c base64.c radixsort.c unifize.c nifti_io.c -DNII2MESH
```

### Linux universal binary

`make` in the `src` folder builds a binary that works only on specific Linux versions. To build a binary that runs on many Linux versions, use [holy-build-box](https://github.com/FooBarWidget/holy-build-box). It uses an old gcc (4.8.5), so the binary may not be fully optimized for your system.

```
git clone https://github.com/rordenlab/niimath
sudo docker run -t -i --rm  -v `pwd`:/io ghcr.io/foobarwidget/holy-build-box-x64 /hbb_exe/activate-exec bash
cd /io/niimath/src
make
exit
sudo chown $(whoami) ./niimath/src/niimath
```

## JavaScript and WebAssembly

niimath compiles to WebAssembly, so you can use it in web pages. The `@niivue/niimath` package runs the processing in a Web Worker and is meant for the browser, not for Node.js. See the [live demo](https://niivue.github.io/niivue-niimath/), which links to source code and instructions, and the [@niivue/niimath README](https://github.com/rordenlab/niimath/blob/master/js/README.md). The rest of this README describes the `niimath` command-line program.

## Usage

niimath accepts the same commands as fslmaths, so you can use it as you would use fslmaths. You can even rename it `fslmaths` and use it as a drop-in replacement. Run `niimath` with no arguments to print the full list of operations. The general form is:

```
niimath [-dt <datatype>] <first_input> [operations and inputs] <output> [-odt <datatype>]
```

Two environment variables control the output:

- `FSLOUTPUTTYPE` sets the output file format, as in [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslEnvironmentVariables). On Unix, run `export FSLOUTPUTTYPE=NIFTI_GZ`, `export FSLOUTPUTTYPE=NIFTI` or `export FSLOUTPUTTYPE=NIFTI_ZST` (zstd compressed; requires zstd support) on the command line or in your profile. On Windows, use `set` instead of `export`.
- `AFNI_COMPRESSOR=PIGZ` compresses `.gz` output with the parallel [pigz](https://github.com/madler/pigz) program. If the variable is not set, or has any other value, compression is single threaded.

To set the number of OpenMP threads, use `-p <threads>` or set `OMP_NUM_THREADS`. `-p` is an operation in the chain, so it goes after the input and before the operations it should affect:

```
niimath in.nii -p 4 -add 1 out.nii
export OMP_NUM_THREADS=4
```

To read from standard input or write to standard output, use the filename `-`. Both streams must be uncompressed single-file NIfTI-1 images. Only one image can be piped at a time. See the [library README](https://github.com/rordenlab/niimath/blob/master/library/README.md) for Python examples.

```
niimath - -add 1 -
```

## Operations not in fslmaths

This section lists the operations that niimath adds. Run `niimath` with no arguments for the fslmaths-compatible operations.

### Simple operations

| Operation | Description |
|---|---|
| `-bptfm <hp> <lp>` | Same as `-bptf` but does not remove the mean. Emulates fslmaths before 5.0.7. |
| `-bwlabel <conn>` | Connected component labeling of non-zero voxels. `conn` sets the neighbors: 6, 18 or 26. |
| `-c2h` | Reverse the `-h2c` transform. |
| `-ceil` | Round voxels up to the nearest integer. |
| `-close <thr> <dx1> <dx2>` | Morphological close. Binarizes at `thr`, dilates by `dx1`, erodes by `dx2`. Fills bubbles with `thr`. |
| `-comply <nx> <ny> <nz> <dx> <dy> <dz> <f_high> <isLinear>` | Conform to axial slices with `nx × ny × nz` voxels of `dx × dy × dz` mm. `f_high` clamps bright voxels (0.98 clamps the top 2%). `isLinear` selects linear (1) or nearest-neighbor (0) interpolation. |
| `-conform` | Reslice to 1 mm voxels in the coronal slice direction with 256³ voxels. |
| `-crop <tmin> <tsize>` | Keep `tsize` volumes starting at volume `tmin`. Volume numbers start at 0, not 1. A `tsize` of -1 means the full range. |
| `-dehaze <mode>` | Set dark voxels to zero. `mode` is 1 to 5. A higher mode keeps more voxels. |
| `-demean` | Remove the average signal across volumes. Requires 4D input. |
| `-detrend` | Remove the linear trend and the mean. |
| `-dilate <thr> <dx>` | Morphological dilate. Binarizes at `thr` and grows up to distance `dx`. |
| `-dog <sPos> <sNeg>` | Difference of Gaussians with zero-crossing edges. Sigmas are in mm. |
| `-dogr <sPos> <sNeg>` | As `-dog`, without zero crossing (raw data, not binarized). |
| `-dogx <sPos> <sNeg>` | As `-dog`, with zero crossing in 2D sagittal slices. |
| `-dogy <sPos> <sNeg>` | As `-dog`, with zero crossing in 2D coronal slices. |
| `-dogz <sPos> <sNeg>` | As `-dog`, with zero crossing in 2D axial slices. |
| `-edginess` | Scalar field of local vector contrast: the Euclidean distance between each voxel and its neighbors. Useful for RGB or multi-channel data. |
| `-edt` | Euclidean distance transform (distance field). Assumes isotropic voxels. |
| `-erode <thr> <dx>` | Morphological erode. Binarizes at `thr` and shrinks within distance `dx`. |
| `-floor` | Round voxels down to the nearest integer. |
| `-gz <mode>` | NIfTI gzip mode: 0 = uncompressed, 1 = compressed, any other value = use the FSL environment. Default -1. |
| `-h2c` | Convert CT scans from Hounsfield to Cormack units, which emphasize soft tissue contrast. |
| `-hollow <threshold> <thickness>` | Hollow out a solid object, keeping a wall of `thickness`. Binarizes at `threshold`. Single-volume input only. Useful before `-mesh`. |
| `-mod` | Fractional modulus remainder. Same as `-rem`, but keeps fractions. |
| `-otsu <mode>` | Binarize with Otsu's method. `mode` is 1 to 5. A higher mode keeps more bright voxels. |
| `-p <threads>` | Set the maximum number of parallel threads. 0 uses all available threads. |
| `-power <exponent>` | Raise each voxel to the given exponent. |
| `-qform <code>` | Set `qform_code`. |
| `-ras` | Reorder and flip dimensions to RAS orientation. |
| `-reslice <target>` | Reslice to match image `target` with linear interpolation. |
| `-reslice_nn <target>` | Reslice to match image `target` with nearest-neighbor interpolation. |
| `-reslice_mask <mask>` | Reslice `mask` onto the current image with nearest neighbor. Voxels where the mask is ≤ 0 are set to the minimum intensity. |
| `-resize <X> <Y> <Z> <m>` | Grow (> 1) or shrink (< 1) the image. Method `m`: 0 = nearest, 1 = linear, 2 = spline, 3 = Lanczos, 4 = Mitchell. |
| `-robustfov [mm]` | Crop to a robust field of view (default 170 mm) from the top of the head down. Removes the lower head and neck. Emulates FSL `robustfov`. Adjusts the dimensions and the sform/qform. |
| `-round` | Round voxels to the nearest integer. |
| `-scale01` | Rescale intensities linearly to the range 0 to 1 with the global minimum and maximum. |
| `-sedt` | Signed Euclidean distance transform (distance field). Assumes isotropic voxels. |
| `-sform <code>` | Set `sform_code`. |
| `-sobel` | Fast edge detection. |
| `-sobel_binary` | Sobel edge detection with a binary result. |
| `-tensor_2lower` | Convert an FSL-style upper-triangle tensor image to the NIfTI-standard lower-triangle order. |
| `-tensor_2upper` | Convert a NIfTI-standard lower-triangle tensor image to the FSL-style upper-triangle order. |
| `-tensor_decomp_lower` | As `-tensor_decomp`, but the input stores the lower triangle (AFNI, ANTS and Camino convention). |
| `-trunc` | Remove the fractional part of each voxel and return the integer value. |
| `-unifize [-GM]` | Bias field correction, adapted from AFNI 3dUnifize. The optional `-GM` also scales gray-matter intensity toward a common target. |
| `-unsharp <sigma> <scl>` | Edge-enhancing unsharp mask. `sigma` is in mm, not voxels (1 is typical). `scl` is the amount (0.5 medium, 1.0 heavy). |

### `-allineate <base> [opts]`

Affine registration of the current image to `base`, from AFNI 3dAllineate. The output is the registered image on the `base` grid.

| Option | Values | Description |
|---|---|---|
| `-cost XX` | `fast`, `fastx` (default), `fasthel`, `fastcr`, `hel`, `nmi`, `lpc`, `lpa`, `ls` | Cost function and engine. See the engines below. The last `-cost` wins. |
| `-cmass`, `-nocmass` | | Initialization. See the engines below. |
| `-source_automask` | | Use with `lpc` or `lpa`. |
| `-warp XX` | `sho`, `shr`, `srs`, `aff` (default) | Transform type. |
| `-interp XX` | `NN`, `linear` (default), `cubic` | Interpolation during matching. |
| `-final XX` | `NN`, `linear`, `cubic` (default) | Output interpolation. `-nearest`, `-linear` and `-cubic` are shortcuts. |
| `-fill XX` | `auto` (default), `zero`, `nan` | Fill for output voxels outside the source field of view. `auto` fills with 0, or with the darkest source voxel if that voxel is negative (CT/Hounsfield air is about -1000, so out-of-FOV voxels read as air, not soft tissue). `zero` is always 0 and is byte-identical to the historical behavior for positive-only MRI. `nan` fills with NaN. |
| `-master <grid>` | image | Estimate at the base resolution, but reslice the result onto `<grid>`. `<grid>` must share the base world frame (for example a higher-resolution template). |
| `-savemat out.json` | | Save the fitted world-space `fixed_to_moving` affine and its inverse as self-describing JSON. |
| `-applymat in.json` | | Reslice the moving image onto `base` with a saved affine. Does no registration. Exclusive with the registration and seed options. Use `-nearest` for label or atlas volumes. |
| `-com` | | Seed: reset the origin to the brightness center of mass. |
| `-sym`, `-symd`, `-symb` | | Seed: fold a midsagittal-plane correction into the header. `-symd` de-obliques the frame first. `-symb` auto-competes both. |
| `-nosagseed` | | Disable the in-MSP rigid seed that `-sym` runs by default. |
| `-zoom` | | Relax the scale range for abnormal-size subjects (for example an infant against an adult template). Needs an ordinary cost, not a fast one. |
| `-weight <img>` | image | Graded base-space weight in the style of AFNI 3dAllineate. See the weights below. |

Engines:

- A bare `-allineate`, `-cost fast` and `-cost fastx` select the adaptive clean-room fast engine. It is inspired by SPM and FLIRT, multiresolution, and typically several times faster than the ordinary engine. On a whole-head base it fits HEL and correlation-ratio coarse candidates independently, selects by HEL dependence × overlap, and continues once through the finer HEL stages. A hard-zeroed (skull-stripped) base activates the deeper rigid-HEL, scale-bracketed-HEL and CR-seeded multi-start, with fine-level arbitration before the final affine polish.
- `-cost fasthel` forces the HEL-only fast trajectory. `-cost fastcr` forces correlation ratio only.
- `-cost hel`, `nmi`, `lpc`, `lpa` and `ls` select the ordinary AFNI-style engine.
- The fast engine runs a fixed schedule with internal sampling, so it rejects `-warp`, `-interp`, `-source_automask`, `-dark_automask` and `-zoom`. Use an ordinary cost for those. The `-com`, `-sym`, `-symd` and `-symb` seeds work with every cost. `-final`, `-master` and `-savemat` are honored.
- With the fast engine, the default and `-cmass` choose between the supplied affine and a center-of-mass recentered start by dependence × overlap. `-com` forces the center of mass. `-nocmass` forces the supplied affine.

Weights (`-weight <img>`):

- `<img>` is a 3D image in base space. Its dimensions and world frame must match `base`.
- Values are normalized to `[0, 1]` (divided by the maximum) and applied per base voxel. A voxel weighted 0 is excluded. A voxel near 1 dominates. The weight is graded. It is not an exclusion mask.
- Keep the head outside the region of interest attenuated but nonzero. If the whole exterior is zero, the global scale is underdetermined and the cross-modal fit collapses into the scalp. An AFNI-style whole-head weight anchors the absolute size.
- Both engines honor the weight. The ordinary engine uses it in place of its manufactured autoweight. The fast engine applies it only at the finest 2 mm stage, so coarse capture and global-scale selection stay whole-head.
- niimath rejects `-weight` with stdin (`-`), with `-applymat`, and when the weight has no positive voxel over the fixed foreground.
- If the default fast engine fails and niimath falls back to the ordinary Hellinger engine, the fallback honors the weight. An explicit fast selector never switches engines.

Constraints:

- A 4D moving image registers its first volume only, equivalent to `-crop 0 1`, with a warning. Use `-Tmean` or `-crop` first to choose a volume.
- Skull stripping here is `-deface` with a brain mask, which is template-driven. [`-skullstrip`](#-skullstrip--faithful) does it with no template at all, but is built only with `SKULLSTRIP=1`. To crop the field of view first, chain `-robustfov` before `-allineate`. Together these make niimath a superset of the standalone `allineate` registration tool.

### `-deface <tmpl> <mask> [opts]`

Removes voxels with a template-space mask. niimath registers the input to `tmpl` (affine), inverts the transform, and warps `mask` onto the input's native grid. Voxels where the warped mask is below 0.5 are set to the image's finite minimum, which is about 0 for typical MRI. The input itself is never resampled.

Arguments:

- `tmpl`: the template image.
- `mask`: a mask in `tmpl` space. A value ≥ 0.5 keeps the voxel. A value < 0.5 removes it. The mask decides what is removed: a brain mask keeps the brain (skull stripping), a face mask removes the face (defacing).

Options:

- `-cost XX`: `fast`, `fastx` (default), `fasthel`, `fastcr`, `hel`, `nmi`, `lpc`, `lpa` or `ls`. The engines are the same as for `-allineate`. The fast engine is the default. `-cost hel` selects the ordinary AFNI-style engine.
- `-cmass`, `-nocmass`: initialization tuning.
- `-final XX`, `-nearest`, `-linear`, `-cubic`: output interpolation. Default `linear`.

Constraints:

- The input must be 3D. A 4D or multi-volume input is rejected.
- The fast engine cannot honor `-warp`, `-interp`, `-source_automask` or `-dark_automask`. Use `-cost hel` for those.
- The `-allineate` workflow options `-savemat`, `-applymat`, `-com`, `-sym`, `-symd`, `-symb`, `-nosagseed`, `-zoom`, `-master`, `-fill` and `-weight` are rejected. The fill stays the image minimum.
- Breaking change, and the name has since been re-used. The former `-skullstrip <tmpl> <mask>` command ran this identical operation and was removed. Replace it with `-deface` and supply your brain mask. The result is unchanged. [`-skullstrip`](#-skullstrip--faithful) now names a different operation, template-free surface skull stripping, which takes no arguments. An old two-argument `-skullstrip` line is detected and rejected before any work is done, with a message naming `-deface`. It does not run the new operation, and it does not misread your template and mask.

### `-reface <tmpl> <shell> <weight> [opts]`

Anonymizes a head image by face replacement. Emulates AFNI `afni_refacer2 -mode_reface`. niimath registers the subject to `tmpl`, back-projects the signed template-space `shell` onto the original subject grid, and composites an artificial face. The output stays on the subject grid.

Arguments:

- `tmpl`: the template image.
- `shell`: the signed face-replacement shell in `tmpl` space. Where the shell is > 0, the voxel is replaced by the shell scaled by a brightness-match factor. Where it is 0, the subject is kept. Where it is < 0, the voxel is set to zero. An edge blend is applied inside the replaced region.
- `weight`: required. A registration weight in `tmpl` space, as for `-allineate -weight`.

Options: `-cost XX`, as for `-deface`. Both engines work. The shell back-projection is always nearest neighbor, so `-final` does not apply.

Constraints:

- The command fails closed for privacy. If less than 10% of the shell mapped into the subject field of view, which indicates a registration failure, niimath refuses to write the image and returns an error.
- Check the output visually.

### `-skullstrip [-faithful]`

AFNI-style surface skull stripping. No template, no mask, no network. niimath expands a surface outward from inside the head until it wraps the brain, the method of AFNI `3dSkullStrip -no_use_edge`, then applies the resulting mask. Further operations may follow.

Options:

- `-faithful`: run the reference deformation kernel instead of the optimized default. Same algorithm, same quality. The two differ only in how the ray walk is compiled and in whether the surface's node loop runs on several threads. `-faithful` reproduces the pre-optimization release exactly, which makes it the thing to diff against when you check for a regression. The default is about 1.8 times quicker, and is what you want otherwise. Across nine varied images on an M4 Pro: 17.4 s with `-faithful`, 9.6 s without. `-faithful` must come before the output filename.

Output: in-mask voxels keep their original intensities. Out-of-mask voxels become the image minimum. This is a deliberate divergence from AFNI's default, which rescales intensities.

Constraints:

- Scalar 3D input only, computed in float32.
- Off by default, and 64-bit native only. Build with `SKULLSTRIP=1 make` or `cmake -DENABLE_SKULLSTRIP=ON`. The wasm, tiny and nano targets reject the request rather than silently dropping it. WebAssembly support is a stated goal of the design, but it is not implemented, so today this is a desktop-only command.
- The name is re-used. It previously aliased the template and mask removal now spelled [`-deface`](#-deface-tmpl-mask-opts). The new `-skullstrip` needs no template and no mask. Read the note under `-deface` before you migrate an old command line.
- The surface converges on a discrete test: an integer count of troubled nodes. A difference far below single-precision rounding can therefore add a whole extra pass. Output is not byte-stable between compilers, between build systems, or between the two kernels. Compare masks with Dice, not `cmp`. Within one binary it is exactly reproducible, including across `-p` thread counts. The default kernel spreads the node loop over threads and still gives bit-identical output at any team size. `-faithful` runs that loop on one thread, so it ignores `-p` for the part that dominates its time.
- This is automated research segmentation, not a diagnostic guarantee. Inspect the output.

### `-qwarp <base>`

Nonlinear (deformable) registration to `base`. An attributed port of AFNI `3dQwarp -blur 0 3`. The output is the warped image on the `base` grid. Further operations may follow.

Constraints:

- `-qwarp` takes one argument and has no sub-options.
- The input must already be unifized, skull-stripped, affine-aligned, and on the `base` grid (same dimensions and world frame).
- Built only with `QWARP=1`. It is off by default because it is memory and CPU heavy and impractically slow in WebAssembly. It requires the allineate engine, so it is unavailable in an `AL=0` build.

### `-romeo <mag|none> [opts]`

ROMEO phase unwrapping. A faithful MIT-licensed C port of [ROMEO.jl](https://github.com/korbinian90/ROMEO.jl) and the [MriResearchTools.jl](https://github.com/korbinian90/MriResearchTools.jl) helpers that its command-line app uses. The current image is the wrapped phase, 3D or 4D with echoes on dimension 4. Further operations may follow.

Arguments:

- `mag`: the magnitude image. This positional argument is required. Pass the literal `none` to unwrap without a magnitude.

| Option | Description |
|---|---|
| `-t <TEs>` | Echo times in ms: `16.8`, `16.8,38.56`, `'[16.8,38.56]'` or `epi [te]`. Quote the bracket form for the shell. Required for multi-echo input, optional for a single echo. |
| `-k <spec>` | Mask: `nomask`, `robustmask` (default), `qualitymask [thr]` (default threshold 0.1), or a mask file. |
| `-w <spec>` | Weights: `romeo` (default), `romeo2`, `romeo3`, `romeo4`, `romeo6`, or up to 6 bits such as `1010`. A bare `romeo` resolves to `romeo3` with a magnitude and `romeo4` without one. |
| `-template <n>` | Echo to unwrap spatially (default 1). |
| `-i` | Individual (not temporal) unwrapping. |
| `-temporal-uncertain-unwrapping [x]` | Re-unwrap low-quality voxels spatially. 0.5 when the flag is bare, off otherwise. |
| `-g` | Correct the global n2π offset. |
| `-q` | Write `<out>_quality`. |
| `-Q` | Write `<out>_quality_1..6`. |
| `-B [name]` | Also write a B0 field map in Hz: `<out>_B0` and `<out>_B0_snr`. Needs `-t`. `[name]` replaces the B0 stem. |
| `-B0-phase-weighting <mode>` | `phase_snr` (default), `phase_var`, `average`, `TEs`, `mag` or `simulated_mag`. |
| `-no-phase-rescale` | Do not rescale the phase. `-no-rescale` is an alias. |
| `-no-mask-out` | Do not write `<out>_mask`. |
| `-v` | Verbose. |

Outputs:

- Side outputs use `nifti_save` postfixes on the output name: `<out>_mask` (only when a mask was computed; `-k nomask` writes none), `<out>_quality` (`-q`), `<out>_quality_1..6` (`-Q`; a map that is uniformly 1.0 in the interior is skipped, as upstream), and `<out>_B0` and `<out>_B0_snr` (`-B`).
- Side outputs honor `FSLOUTPUTTYPE`. They honor `-gz` and `-p` only when those precede `-romeo`, because the side outputs are written during the operation. A later `-gz` reaches only the main output.

Constraints:

- The phase is rescaled to `[-π, π]` as in `readphase`, unless `-no-phase-rescale` is given. The rescale re-reads the unscaled stored values, so `-romeo` must be the first computational operation when rescaling is active.
- `-B` computes B0 without MCPC-3D-S phase-offset correction, which ROMEO's own multi-echo `-B` enables silently. The maps correspond to `romeo --compute-B0 --phase-offset-correction off`. niimath says so on stderr. Without a magnitude, ROMEO's SNR map collapses to one value, because the substituted `exp(-TE/20)` decay does not depend on the voxel. niimath writes that constant across the working grid.
- Not yet ported, and rejected with a specific message: `-u`, `-e`, `-threshold`, `-w bestpath`, `-max-seeds > 1`, `-merge-regions`, `-correct-regions`, `-wrap-addition != 0` and `-fix-ge-phase`. MCPC-3D-S phase-offset correction and multi-channel (5D) input are out of scope.
- Enabled by default. `ROMEO=0 make` or `-DENABLE_ROMEO=OFF` omits it.

Citation: Dymerska, B. et al. 2020, *Magnetic Resonance in Medicine*, [doi:10.1002/mrm.28563](https://doi.org/10.1002/mrm.28563).

### `-unwarp <map> <axis>`

EPI distortion correction. Resamples the current image through a scalar displacement map in millimeters, as written by `--medic`. Further operations may follow.

Arguments:

- `map`: a 3D map (applied to every frame) or a 4D map with the same frame count as the input. It must share the input's dimensions and world transform.
- `axis`: `i`, `j` or `k` (or `x`, `y`, `z`). A trailing `-` is accepted and ignored, because the sign is already stored in the map and a second negation would double-correct. `--medic --phase-encoding-direction` is the opposite: it honors the suffix. Pass the full BIDS value there.

Method: an unnormalized Lanczos-5 windowed sinc, applied separably in 3D, with zero fill outside the field of view and no Jacobian intensity modulation. This matches the measured behavior of the reference implementation. See the `medic_bench` repository.

Enabled by default with `--medic`. `MEDIC=0 make` or `-DENABLE_MEDIC=OFF` omits it.

### `-fugue <fieldmap> <dwell> <unwarpdir>`

B0 fieldmap EPI distortion correction. Resamples the current image along one phase-encoding axis to undo susceptibility distortion. Emulates FSL `fugue` in its `--loadfmap`, `--dwell` and `--unwarpdir` mode. Further operations may follow.

Arguments:

- `fieldmap`: a 3D B0 field map in rad/s, on the input's voxel grid. Those are the units `-fmapprep` and FSL `fsl_prepare_fieldmap` write. A 4D input gets the same map on every volume.
- `dwell`: the effective echo spacing in seconds, the BIDS `EffectiveEchoSpacing`.
- `unwarpdir`: `x`, `y` or `z` (or `i`, `j`, `k`), with an optional trailing `-`.

Method: the shift in voxels is `fieldmap / (2π) × dwell × N`, where `N` is the image dimension along the unwarp axis. Resampling is 1D linear interpolation along that axis, with zero fill outside the field of view and no Jacobian intensity modulation. The shift is sampled at the output voxel, so this is a pure pull and signal is not conserved under compression. That is the reference behavior.

A zero field map means no data, not zero shift. The shift field is extrapolated over unsupported voxels along each line parallel to the unwarp axis before resampling. Skipping that step is not an edge refinement. It changes about a third of the voxels in a real image.

Constraints:

- The image runs in float32 and must share the field map's dimensions and world transform.
- The image is transformed in place, one line at a time, so peak memory is the input plus one shift volume.
- Byte-identical across thread counts.
- Enabled by default. `FMAP=0 make` or `-DENABLE_FMAP=OFF` omits `-fugue` and `-fmapprep`.

Provenance: original BSD-2 code. FSL is under a non-commercial license that is incompatible with BSD-2, so its fugue and prelude sources were never read. The executables served only as a black-box oracle.

Citation: Jezzard, P. & Balaban, R. S. 1995, *Magnetic Resonance in Medicine* 34:65-73, [doi:10.1002/mrm.1910340111](https://doi.org/10.1002/mrm.1910340111).

### `-fmapprep <mag> <dTE_ms> [-no-debranch]`

Builds the rad/s field map that `-fugue` consumes. The current image is the wrapped two-echo phase difference. Emulates `fsl_prepare_fieldmap SIEMENS`. Further operations may follow, so you can chain `-fugue` straight after it.

Arguments:

- `mag`: the brain-extracted magnitude belonging to the phase difference, for example `bet` output. Its non-zero voxels are the mask, and it weights the unwrapping.
- `dTE_ms`: the echo time difference `EchoTime2 - EchoTime1`, in milliseconds.

Options:

- `-no-debranch`: skip the 2π branch-outlier correction described below and keep the raw ROMEO field.

Method: the mask is the non-zero magnitude, used verbatim, with no erosion and no dilation. niimath rescales the phase so its observed range spans exactly 2π. It then unwraps in 3D with ROMEO, divides by the echo time difference, and subtracts the median over the mask. The output is 0 outside the mask. No regularization is applied: no median filter, no despiking, no smoothing and no erosion. That is measured reference behavior, not an omission. Unwrapping is ROMEO, not the reference's PRELUDE, so the two disagree at poorly conditioned voxels.

Branch correction: ROMEO leaves isolated voxels a full 2π away from their neighbors where PRELUDE does not, and each becomes a large local shift error. niimath re-projects them. Every change is an integer multiple of 2π/ΔTE, so it cannot alter a genuine field value. The correction runs after ROMEO returns, so `-romeo` and `--medic` keep a faithful ROMEO. `-no-debranch` recovers the raw field for parity work.

Constraints:

- The input must be a single 3D volume in float32, and the magnitude must share its grid.
- `-fmapprep` needs ROMEO. `ROMEO=0` drops `-fmapprep` and keeps `-fugue`.
- Enabled by default. `FMAP=0 make` or `-DENABLE_FMAP=OFF` omits both commands.

### `-moco [-ref <n|image>] [-1Dfile <path.1D>] [-relative]`

Rigid-body motion correction. Registers every volume of a 4D series onto a reference and replaces the image with the corrected series. A clean-room BSD-2 implementation of the method of Cox & Jesmanowicz (*Magnetic Resonance in Medicine* 42:1014-1018, 1999), the algorithm behind AFNI `3dvolreg`. Further operations may follow.

```
niimath bold -moco out
niimath bold -moco -1Dfile out.1D out
niimath bold -moco -ref 5 out
niimath bold -moco -ref ref.nii out
niimath bold -moco -relative rel.1D
```

Options, which may appear in any order:

- `-ref <n>`: register onto volume `n` of the series. `-ref 0` is the default. That volume is copied through unchanged, and its parameter row is all zeros.
- `-ref <image>`: register onto an external reference image instead. It must sit on the input's voxel grid: the same dimensions, and a voxel-to-world transform agreeing within 0.001 mm. A reference on any other grid is rejected, not resliced. A 4D reference contributes its volume 0. With an external reference no volume is copied through: every volume is registered and gets a non-zero parameter row.
- An all-digit argument to `-ref` is a volume number. To name a file called `7`, write `./7`.
- `-1Dfile <path.1D>`: also write the six motion parameters per volume. The filename must end in `.1D`. The file is compatible with AFNI's `-1Dfile`: six columns `roll pitch yaw dS dL dP`, one row per volume. Rotations are in degrees counter-clockwise about the I-S, R-L and A-P axes. Shifts are in mm toward Superior, Left and Posterior. The values record the correction that was applied, not the estimated motion.
- `-relative`: measure each volume against the original previous volume, as a quality-control statistic. niimath registers nothing and writes no image, so the trailing output name is the parameter file itself and must end in `.1D`. `-ref` and `-1Dfile` are rejected with it, and it must be the last operation. Rows are written at `%12.8f`, so `numpy.loadtxt` reads the file. Row 0 is all zeros, because volume 0 has no predecessor. This costs about ten times the work per volume that ordinary correction does.

Constraints:

- The input must be 4D with more than one volume.
- A singleton spatial dimension is rejected in every mode, because the fit is degenerate there.
- Correction runs in float32, so `-dt double` is rejected.
- Enabled by default on every platform, including WebAssembly. `MOCO=0 make` or `-DENABLE_MOCO=OFF` omits it.

### `-stc --slicetiming <t0,t1,...|@file> [-tzero <sec>]`

Slice-time correction. Shifts every voxel time series of a 4D series so that all slices share one temporal origin. A clean-room BSD-2 implementation of the default Fourier method (detrend, interpolate, retrend) of AFNI `3dTshift`. That source was used only as a black-box oracle. It was GPL-2 when this code was written, and CC BY 4.0 since the 12 May 2026 relicense. Further operations may follow.

```
niimath bold -stc --slicetiming @times.1D out
```

Arguments:

- `--slicetiming`: required, case-sensitive, and must come first. It takes one comma-separated list of slice acquisition times in seconds, or an AFNI-style `@file` (whitespace- or comma-separated; `#` starts a comment). Supply exactly one value per slice along storage axis `k`, in slice-index order. A count that does not match `nz` is an error, not a silently truncated list.
- `-tzero <sec>`: the common time point. Default: the arithmetic mean of the supplied times. It must lie within their `[min, max]`.

Output: only `toffset` changes in the header, to the common time point in the header's own time units. Geometry, TR and the spatial transforms are untouched.

Constraints:

- The input must be a scalar 4D image with at least 5 volumes, a finite positive `pixdim[4]`, and a usable temporal unit in `xyzt_units` (seconds, milliseconds or microseconds). niimath does not assume seconds.
- Correction runs in float32, so `-dt double` is rejected.
- This version corrects along storage axis `k` only.
- Enabled by default on every platform, including WebAssembly. `STC=0 make` or `-DENABLE_STC=OFF` omits it.

BIDS helper: `test/stc_slicetiming.py` (standard library only) reads a BIDS sidecar and prints the `--slicetiming` argument. It honors `SliceEncodingDirection`: `k` passes through, `k-` is reversed into slice-index order, and `i` and `j` are rejected. It stops with an error if the sidecar's `RepetitionTime` disagrees with the unit-normalized header TR:

```
niimath bold.nii.gz -stc --slicetiming "$(python3 test/stc_slicetiming.py bold.nii.gz)" out.nii.gz
```

### `-spm_coreg <ref> [opts]`

SPM rigid-body coregistration of the current image to `ref`. Requires the [optional copyleft module](#optional-copyleft-module).

Options:

- `-cost XX`: `nmi` (default), `mi`, `ecc`, `ncc` or `ls`.
- `-sep`, `-fwhm`, `-dither 0|1`, `-coarse sparse|downsample`, `-verbose 0|1`.
- `-interp trilinear|nearest` (default `trilinear`) and `-fill zero|nan` (default `zero`) control the reslicing onto the `ref` grid.
- `-estimate`: rewrite only the source sform/qform instead of reslicing.

### `-spm_deface <tmpl> <mask> [opts]`

The SPM analog of `-deface`. Registers with `-spm_coreg`. Requires the [optional copyleft module](#optional-copyleft-module). Options: the same estimate options as `-spm_coreg`, plus `-interp`.

### `-mesh [opts] <output>`

Converts the current image to a triangulated mesh. The output filename extension selects the mesh format. See [Creating meshes](#creating-meshes) for examples.

| Option | Description |
|---|---|
| `-i <isovalue>` | Isosurface: `d` (dark), `m` (medium), `b` (bright) or a number. The `d`, `m` and `b` values use Otsu's method. Default `m`. |
| `-a <atlasFile>` | Mesh each region of an atlas. One mesh per label, named by the label file. |
| `-b <0\|1>` | Fill bubbles, meaning interior cavities, before meshing. Default 0. |
| `-l <0\|1>` | Keep only the largest connected object. Default 1. |
| `-n <0\|1>` | Use the half-edge simplifier. Only in a `Q2=1` build. |
| `-o <0\|1>` | 1 = classic marching cubes tables, no ambiguity resolution. Default 0. |
| `-p <0\|1>` | Smooth the volume before marching cubes. Default 1. |
| `-q <0\|1\|2>` | Quality: 0 quick, 1 balanced, 2 precise. Default 2, which adds self-intersection guards and a lossless finish. |
| `-s <iterations>` | Post-smooth the mesh with Humphrey's Classes. Default 0. |
| `-r <fraction>` | Reduce to this fraction of the triangles. 1 means no simplification. Default 0.25. |
| `-v <0\|1>` | Verbose. Prints a mesh check after every stage. Default 1. |

A mesh can also be the input. `niimath in.mz3 -s 10 -r 0.5 out.mz3` accepts `-r`, `-s`, `-q` and `-v`.

### `-bitmap [overlay] [opts] <output.png>`

Creates a PNG image from the current volume, with an optional overlay volume. The arguments are inspired by FSL `slicer`. See [niimath-bitmap](https://github.com/rordenlab/niimath-bitmap) for examples and documentation.

Slice selection. Values from 0.0 to 1.0 are fractional positions. Negative values are absolute slice numbers.

| Option | Description |
|---|---|
| `-a` | Axial, coronal and sagittal slices at the midpoint (0.5). |
| `-m` | Mosaic view with slices at 0.25, 0.5 and 0.75 for each axis. |
| `-o` | Select the largest plane orientation automatically. |
| `-x <val1> [val2...]` | Sagittal slices at the given positions. |
| `-y <val1> [val2...]` | Coronal slices at the given positions. |
| `-z <val1> [val2...]` | Axial slices at the given positions. |
| `-X`, `-Y`, `-Z <vals>` | As `-x`, `-y`, `-z`, but draw crosshairs from the other axes. |
| `-r` | Insert a row separator between slice groups. |

Display options:

| Option | Description |
|---|---|
| `-f [0\|1]` | Flip left-right. 0 = neurological, 1 = radiological (default). |
| `-u [0\|1]` | Show L/R labels (default 1). |
| `-n [0\|1]` | Interpolation. 0 = nearest neighbor, 1 = linear. |
| `-s <scale>` | Scale factor for the output image size. |

Color and intensity:

| Option | Description |
|---|---|
| `-t <min> <max>` | Intensity range for the base image. |
| `-T <min> <max>` | Intensity range for the overlay image. |
| `-c <lut> [alpha]` | Base image color lookup table. |
| `-c <R> <G> <B> <A>` | Base image RGBA tint (values 0.0 to 1.0). |
| `-C <lut> [alpha]` | Overlay color lookup table. |
| `-C <R> <G> <B> <A>` | Overlay RGBA color (values 0.0 to 1.0). |
| `-b <R> <G> <B> <A>` | Background RGBA color (values 0.0 to 1.0). |
| `-N [0\|1]` | Use a negative colormap for the overlay (blue-green for negative values). |
| `-e` | Apply edge detection to the overlay. |

Color lookup tables: `gray`, `red`, `green`, `blue`, `cyan`, `yellow`, `bluegreen`, `redyellow`, `viridis`, `inferno`, `magma`, `plasma`.

```
niimath T1.nii -bitmap -a output.png
niimath T1.nii -bitmap fmri.nii -C red 0.5 -z 0.5 overlay.png
niimath T1.nii -bitmap -m -c viridis mosaic.png
```

### `--dtifit -k <dwi> -r <bvec> -b <bval> -o <base> [-m <mask>] [-xflip 0|1|auto]`

Linear diffusion tensor fit. Emulates FSL `dtifit`. This is a self-contained command with its own arguments. The fit math comes from AFNI 3dDWItoDT (public domain).

Outputs: `<base>_{FA,MD,L1,L2,L3,V1,V2,V3,S0,MO,tensor}`.

Options:

- `-m <mask>`: restrict the fit to a mask.
- `-xflip auto` (default) flips the bvec X component when the spatial transform determinant is positive, which matches FSL. `0` never flips. `1` always flips.

### `--qc <t1> --seg <seg> --csf <i[,j..]> --wm <i[,j..]> [--erode 0|1] [--air <template>] [--out qc.tsv] [--json qc.json]`

MRIQC-style anatomical quality metrics from a T1 image and an integer segmentation. This is a self-contained command with its own arguments.

Output: a wide TSV (`--out`, the default `qc.tsv` when neither output is named) and/or an MRIQC-style JSON (`--json`: the same metrics flat at the top level, plus `size_*`, `spacing_*` and `provenance`; non-finite values are `null`). Metrics: CJV, cnr_noair, per-tissue and total SNR, WM2MAX, efc_brain, ICV fractions with mm³ volumes, and per-tissue summary statistics.

Labels: `0` is non-brain and is excluded. `--csf` and `--wm` give disjoint sets of CSF and WM label values. Every other non-zero label is GM.

`--air <template>` adds the background metrics, which need an air region: `summary_bg_*`, Dietrich SNR (`snrd_*`), `fber`, `qi_1`, and `cnr` with its air term (`cnr_noair` stays for comparison). The method follows MRIQC's `ArtifactMask`: the image is made RAS-canonical, registered to the template (`-allineate`, transform only), and the air is the head-free region (Otsu head mask, filled and closed) superior to MRIQC's landmark plane at template z = −14. Voxels brighter than 10 MADs of that region, outside the 10 % shell nearest the head and surviving a 6-connected opening, are artifacts (`qi_1` is their fraction) and are pruned from the background before its statistics. Any T1 template in the MNI frame works; MRIQC's landmarks were placed on `avg152T1`. Needs a build with allineate.

```
niimath --qc T1.nii.gz --seg seg.nii.gz --csf 3,4,11,12 --wm 1,5 --air avg152T1.nii.gz --json qc.json
```

Constraints: this hard-segmentation variant uses unrounded intensities and NumPy-linear percentiles, so its values are not numerically interchangeable with MRIQC's soft partial-volume summaries; the head mask is Otsu-based where MRIQC's is gradient-based. The names `cnr_noair` and `efc_brain` flag their deviation from the MRIQC norms.

### `--medic --magnitude <e1> [<e2> ...] --phase <e1> [<e2> ...] --te-ms <t1,t2,...> --total-readout-time <sec> --phase-encoding-direction <i|j|k|i-|j-|k-> --out-prefix <path> [options]`

MEDIC multi-echo distortion correction. Estimates a B0 field map per frame from multi-echo phase and converts it to an EPI displacement map that `-unwarp` can apply. This is a self-contained command with its own arguments. Run `niimath --medic --help` for the full option list.

Outputs: `<prefix>_fieldmaps_native` (Hz, distorted grid), `<prefix>_fieldmaps` (Hz, undistorted grid) and `<prefix>_displacementmaps` (mm). All are float32.

Method: at least two echoes are required. For each frame, niimath builds a tiered brain mask: core, border ring and outside. It removes the MCPC-3D-S phase offset with a global 2π branch selection, unwraps with ROMEO, and applies a per-echo intra-frame 2π offset. Across frames it then applies a temporal 2π consistency correction to the unwrapped phase. It fits a magnitude-weighted field map and truncates that series to a low rank. The brain interior is truncated at `--rank`. The border ring is smoothed spatially and rebuilt from far fewer components.

| Option | Description |
|---|---|
| `--rank <N>` | Low-rank truncation of the field-map series. Default 10. `0` disables it. It applies to the brain interior. The border ring is handled by `--border-regularization`. |
| `--mask-mode <tiered\|robustmask\|mindgrab>` | Brain mask. Default `tiered`. `robustmask` is the pre-2026 behavior, kept for bisection. `mindgrab` shells out to the external `brainchop-mindgrab` and fails if it is absent. |
| `--branch-correction <0\|1>` | Global 2π handling. Default 1. Covers ROMEO's global correction at both unwrap stages and the MCPC-3D-S branch selection. |
| `--echo-offset <0\|1>` | The intra-frame per-echo 2π offset. Default 1. |
| `--border-regularization <0\|1>` | Graded regularization at the brain border. Default 1. `0` gives the plain global truncation. |
| `--temporal-correction <0\|1>` | Temporal 2π consistency correction. Default 1. |
| `--phase-offset <mcpc\|none>` | MCPC-3D-S phase-offset correction. Default `mcpc`. |
| `--noise-frames <N>`, `-f` | Drop N trailing frames from the outputs. Default 0. |
| `--weights <sel>` | ROMEO weight preset: `romeo`, `romeo2`, `romeo3`, `romeo4` (default) or `romeo6`. Governs both the MCPC-3D-S and the multi-echo unwrap. |
| `--mask <file>` | Use this mask verbatim for both unwrapping stages, instead of the default tiered mask. Mutually exclusive with `--mask-mode`. See the mask rule below. |
| `--save-intermediates` | Also write the per-echo unwrapped phase, the masks, and the estimated phase offset when MCPC-3D-S runs. |
| `--n-cpus <N>`, `-n` | OpenMP threads. |
| `--gz <0\|1>` | Output compression. Default: the `FSLOUTPUTTYPE` environment. |

Bisection: each of `--mask-mode robustmask`, `--branch-correction 0`, `--echo-offset 0` and `--border-regularization 0` restores the pipeline as it stood immediately before that stage was added. That is how a change in output is attributed to a stage. No single combination returns to the original, because the temporal grouping's switch from magnitude to unwrapped phase has no flag. It corrects what the statistic measures rather than exposing a choice. The grouping also presupposes that the global branch has been pinned, so `--branch-correction 0` leaves every frame alone in its group and the temporal correction inactive. niimath says so when that happens.

Phase-encoding polarity: give `--phase-encoding-direction` the BIDS value with its sign (`j-`, not `j`). `--medic` honors the `-` suffix. It drives the inversion and flips the sign of the displacement map. `-unwarp` ignores the suffix, because the sign is already stored in the map it reads.

Mask rule: a supplied `--mask` is binary, so it becomes the core tier with no border ring, and `--border-regularization` has nothing to do. `--mask` counts a voxel as inside the brain when its value is `>= 1`, not merely non-zero. That is the measured convention of the reference tool. A fractional probability map is not a mask. Threshold it first: `niimath p.nii -thr 0.5 -bin mask.nii`. NaN is treated as outside. A mask with no voxel `>= 1` is an error, not an empty result.

Memory: the whole run is held in RAM by design. A 4D `.nii.gz` cannot be seeked, so streaming gains nothing. The work arrays are `nx·ny·nz × frames × (2·echoes + 3) × 4` bytes. niimath prints that budget at startup. The peak adds one echo pair of input (echoes are loaded and released one at a time) and one echo pair in transit. The measured peak for 170 frames × 2 echoes at 76×76×46 is 1.54 GB single-threaded and 2.04 GB at 8 threads writing uncompressed. The reference tool needs 3.50 GB and 4.21 GB for the same run. The default mask and border stages add allocations the printed budget does not cover: a per-frame magnitude mask, a per-thread morphology workspace, and the border filter's own buffers. Those are most of the gap between the banner and the 8-thread peak. Memory grows linearly with frames × echoes. Where memory is tight, estimate natively and apply `-unwarp` separately. A WebAssembly build has a 4 GiB ceiling.

Fidelity: `--medic` is a clean-room emulation. It was developed from the published method, from black-box measurement of the reference tool, and from two non-source documents supplied by the method's author: a theory note on phase-offset ambiguities, and correspondence. No reference source was read at any point. Every convention it implements is recorded in the `medic_bench` repository, which also holds the benchmarks and the patent analysis. **No equivalence with the reference tool is claimed.** `-unwarp` does reproduce it closely: fed the reference's own displacement map, it matches the reference's corrected images to nrmse 3.5e-5. `--medic` now agrees closely over the brain interior while differing at the border, so read the two separately. Restricted to the core tier, the native field maps agree at r 0.995 and p99 3.5 Hz. The displacement maps agree to 0.02 mm median and 0.3 mm p99, excluding folds. The corrected images correlate 0.99991 end to end. Including the three-voxel border ring, the global field p99 is about 60 Hz, because that is where the reference's own field is broadband. Quoting the global figure alone misrepresents both. The reference's brain-mask construction is now emulated, at valid-tier Dice 0.9998, and so is its border-aware low-rank filtering. Its iteration-limited field inversion is the one part deliberately not reproduced.

BIDS wrapper: `medic.py` in the `medic_bench` repository (standard library only) discovers multi-echo runs, reads the parameters from their JSON sidecars, and runs `--medic` and `-unwarp` for you.

Citation: Van et al. 2026, *Imaging Neuroscience* 4, [doi:10.1162/IMAG.a.1262](https://doi.org/10.1162/IMAG.a.1262), and the ROMEO reference above for the unwrapping.

### `--compare [<thresh>] <ref>`

Reports whether the current image and `ref` are identical, then exits without saving an image. With `<thresh>`, the exit code is success if the largest difference is less than `thresh`. See [Compatibility with fslmaths](#compatibility-with-fslmaths) for a sample report.

### `niimath <filename.nii>`

With an input file and no other arguments, niimath prints the header, as `fslhd` does, and exits without saving an image. To save the report to a text file, redirect stderr: `niimath T1.nii 2> T1.txt`.

## Creating meshes

niimath converts NIfTI images to meshes for Surfice, Blender, SUMA, FreeSurfer and other tools. The features come from [nii2mesh](https://github.com/neurolabusc/nii2mesh) and are almost identical. The argument order differs, to match fslmaths and niimath. The call `nii2mesh -r 1 bet.nii.gz r100.ply` becomes `niimath bet.nii.gz -mesh -r 1 r100.ply`.

With niimath you can apply voxel operations before you create the mesh, for example the morphological operations `-close`, `-ero` and `-dilM`. To apply a 4 mm Gaussian smooth before meshing:

```
niimath mni152.nii.gz -s 4 -mesh -i 122 -l 0 -b 1 b1.ply
```

To create one mesh for each region of an atlas, as described on the [nii2mesh](https://github.com/neurolabusc/nii2mesh) page:

```
niimath D99_atlas_v2.0_right.nii.gz -mesh -p 0 -s 10 -a D99_v2.0_labels_semicolon.txt ./gii/D99s10roi.gii
```

Both programs set the isolevel with `-i`. With `-i 128`, the surface encloses the voxels brighter than 128. niimath also accepts `-i d`, `-i m` and `-i b` for dark, medium and bright. These use Otsu's method and usually find pleasing values. If you give no isolevel, nii2mesh uses the midpoint between the darkest and brightest value, while niimath uses the medium Otsu threshold, which is more robust to outliers.

```
niimath bet.nii.gz -mesh -i 128 Isolevel128.gii
niimath bet.nii.gz -mesh -i d darkIsolevel.gii
niimath bet.nii.gz -mesh -i m medIsolevel.gii
niimath bet.nii.gz -mesh -i b brightIsolevel.gii
```

### Mesh quality

Marching cubes uses the Lewiner tables by default. They resolve the ambiguous cube configurations against the trilinear interpolant. `-o 1` selects the classic tables instead.

Simplification never breaks the topology of the surface. A link condition rejects any collapse that would create a non-manifold edge. `-q` sets the quality level. `-q 2`, the default, adds a self-intersection guard to both the `-s` smoothing and the simplification, plus a lossless finishing pass. `-q 1` drops the guards and the lossless finish. It runs in about 40% of the time: 2.3 s against 5.7 s for a smoothed, simplified MNI surface. `-q 0` is fastest and writes uncompressed mz3. `-r 1` leaves the mesh unsimplified at every quality level.

With `-v 1` every stage prints a `mesh check` line: components, Euler characteristic and genus, boundary edges and holes, non-manifold edges and vertices, and self-intersecting triangles.

A second, half-edge simplifier (`-n 1`) exists as a reference. It stops within one face of the requested count and refuses non-manifold input, at the same geometric fidelity and about 30% more time. It is not compiled by default. Build it with `Q2=1 make` or `-DENABLE_QUADRIC2=ON`.

An existing mesh can be processed the same way:

```
niimath in.mz3 -s 10 -r 0.5 out.mz3
```

## Creating bitmaps

Use `-bitmap` to visualize the result of any chain of operations. See the [`-bitmap` reference](#-bitmap-overlay-opts-outputpng) above and the [niimath-bitmap](https://github.com/rordenlab/niimath-bitmap) repository for examples and documentation.

## Compatibility with fslmaths

### Identical versus equivalent results

niimath is designed to give results equivalent to fslmaths. In most cases the results are identical. In almost all other cases they are equivalent. The results are not always identical because both tools compute in floating point, where the precise order of instructions creates small rounding differences. As [Kernighan and Plauger](https://www.amazon.com/Elements-Programming-Style-Brian-Kernighan/dp/0070341990) wrote: `Floating point numbers are like piles of sand; every time you move one you lose a little sand and pick up a little dirt.` Raw brain imaging data are usually stored as 16-bit integers, and the signal-to-noise ratio is usually a fraction of that dynamic range. niimath computes in single (32-bit) or double (64-bit) precision. So niimath may give results that are not identical, but they are intended to be always comparable. For more on floating point accuracy, see [here](https://introcs.cs.princeton.edu/java/91float/) and [here](http://www.freshsources.com/page1/page7/files/Sand-1.pdf).

The `--compare` operation compares the results of niimath and fslmaths directly. A [validation repository](https://github.com/rordenlab/niimath_tests) runs hundreds of commands to check the output. Its `batch.sh` script tests the functions that give identical results. Its `close.sh` script tests the functions that give equivalent but not identical results. For example, in tensor decomposition the vector [1 0 0] is functionally identical to [-1 0 0], because fiber tracking ignores the polarity of the direction. When `--compare` finds a difference, it prints a report so you can judge whether the results are equivalent:

```
Images Differ: Correlation r = 1, identical voxels 73%
 Most different voxel -69.3133 vs -69.3133 (difference 1.52588e-05)
 Most different voxel location 43x17x49 volume 39
Image 1 Descriptives
 Range: -472.393..491.385 Mean -0.00121971 StDev 6.8898
Image 2 Descriptives
 Range: -472.393..491.385 Mean -0.00121971 StDev 6.8898
    86.29 real    41.08 user    23.41 sys
```

### Known differences

These operations give meaningfully different results. Each item gives the reason:

1. The command `fslmaths inputimg -add 0 outputimg -odt input` can convert a uint8 image to float output, despite the explicit request to keep the input type. This happens when the header has a non-unitary scale slope or a non-zero intercept. niimath keeps both the datatype and the intensity scaling parameters.
2. Versions of fslmaths differ for the pass-through `fslmaths in out`, which is useful for copying files. Old versions save losslessly in the input datatype. fslmaths 6.0 converts the data to float. niimath keeps the datatype.
3. The fslmaths function `-fillh26` sometimes fills unconnected regions. An example was sent to the FSL team. niimath gives the correct solution.
4. The fslmaths function `-dilD` does not do what it claims. It introduces a blur that reduces the edge artifacts of iterative morphology. The blur runs in a fixed order, so it shifts the signal spatially. niimath does the dilation as described. [Better solutions](https://github.com/neurolabusc/niiSmooth) exist for these functions. The niimath `-edt` operation can also dilate.
5. The fslmaths function `-roc` works differently than its help describes. It appears to ignore voxels near the image edge, and it reports "given object has non-finite elements" if any dimension is less than 12 voxels. With an external noise file, it adds undocumented columns to the output file. It does not detect the requested `AROC-thresh` precisely, but samples at stepped intervals. niimath emulates the stepped intervals for reporting, but finds the precise cutoff.
6. The fslmaths help says: `If you apply a Binary operation (one that takes the current image and a new image together), when one is 3D and the other is 4D, the 3D image is cloned temporally to match the temporal dimensions of the 4D image.` This is not the case for `-thr` and `-uthr`. If the second image is 4D, only its first volume is used and the output stays 3D. `-uthr` is odd: `fslmaths 3D -uthr 4D out` fills the 3D input with zeros regardless of the mask values.
7. `fslmaths in1 -rem 0 out` throws an exception, which is understandable. `fslmaths in1 -rem in2 out` also throws an exception if any voxel in `in2` is zero. niimath describes this error.
8. The fslmaths function `-rem` returns the **integer** modulus remainder, like the C `%` operator. This may be unexpected: in Python `2.7 % 2` is 0.7, as in Matlab's `mod(2.7, 2)` and the standard C `fmod`. niimath clones the fslmaths behavior and adds `-mod`, which returns the fractional remainder.
9. fslmaths accounts for a negative determinant by flipping the first dimension. fslstats does not, so fslstats coordinates are often misleading. For an image in RAS orientation, `fslstats tfRAS -x` gives coordinates that are incompatible with the fslmaths `tfceS` function. niimath emulates fslmaths for the relevant functions (`-index`, `-roi`, `-tfceS`).
10. Neither `-subsamp2` nor `-subsamp2offc` applies anti-aliasing. `-subsamp2offc` shows odd edge effects. For slices in the middle of a volume, an output slice is weighted 50% from the center slice and 25% each from the slices below and above, which makes sense. At the edges (the first and last slices, rows and columns) the filter weights 75% on the central slice and 25% on the neighbor, so the neighbor's signal is heavily diluted. A better mixture is 66% edge slice and 33% neighbor. niimath uses the latter.
11. fslmaths 6.0.0 to 6.0.3 cannot process files when the string ".nii" appears in a folder name. For the folder "test.niim", `fslmaths ~/test.niim/RAS -add 0 tst` [throws an exception](https://github.com/FCP-INDI/C-PAC/issues/976). niimath recognizes that this is a folder name, not a file extension, and works. niimath helped detect this anomaly. It is an example of how a clone gives feedback to the original project.
12. The fslmaths function [`-ztop`](https://github.com/rordenlab/niimath/issues/8) does not clamp extreme values.

Some edge cases may remain where niimath does not replicate fslmaths. This is new software, and many fslmaths operations are undocumented. If you find a problem, open a GitHub issue.

## Performance

These speedup factors compare niimath with fslmaths. The T1-weighted and resting-state data use the [HCP 3T Imaging Protocol](http://protocols.humanconnectome.org/HCP/3T/imaging-protocols.html) sequences. The first table is from a macOS laptop with four cores (8 threads, 28 W):

| Command : Seconds (GZ)                                 |  Serial (GZ)  | Parallel (GZ) |
|--------------------------------------------------------|--------------:|--------------:|
| fslmaths rest -s 2.548 out : 270 (424)                 | 5.0x (2.9x)   | 8.6x (6.3x)   |
| fslmaths t1 -kernel boxv 7 -dilM out : 216 (228)       | 245x (41x)    | 225x (72x)    |
| fslmaths rest -Tmean -mul -1 -add rest out : 101 (328) | 2.5x (2.5x)   | 2.8x (4.5x)   |
|  niimath rest -demean out (same output as above)       | 3.5x (3.0x)   | 4.6x (6.2x)   |
| fslmaths rest -bptf 77 8.68 out : 998 (1155)           | 2.0x (2.0x)   | 6.8x (6.7x)   |

The second table is the same tests on a desktop with twelve cores (24 threads, Ryzen 3900X):

| Command : Seconds (GZ)                                 |  Serial (GZ)  | Parallel (GZ) |
|--------------------------------------------------------|--------------:|--------------:|
| fslmaths rest -s 2.548 out : 123 (229)                 | 4.2x (2.4x)   | 9.9x (12.1x)  |
| fslmaths t1 -kernel boxv 7 -dilM out : 156 (159)       | 371x (37x)    | 371x (248x)   |
| fslmaths rest -Tmean -mul -1 -add rest out : 32 (186)  | 1.7x (2.5x)   | 1.8x (7.6x)   |
|  niimath rest -demean out (same output as above)       | 2.6x (2.6x)   | 3.0x (10.8x)  |
| fslmaths rest -bptf 77 8.68 out : 887 (1019)           | 2.6x (2.5x)   | 23x (23.0x)   |

Gaussian smoothing (`-s`, `-dog` and `-unsharp`) uses a contiguous vectorizable kernel in every build: native, WASM, and the shared registration pyramid. The neighborhood mean, minimum, maximum and erosion filters keep their local gathers, but they evaluate adjacent interior outputs in SIMD lanes. This avoids the non-finite propagation errors of separable running-sum and deque filters.

## License

<!-- codespell-ignore-line --> niimath is licensed under the 2-Clause BSD License. Except where noted, Chris Rorden wrote the code in 2020-2022. Daniel Glen of the US National Institutes of Health wrote the code in `tensor.c` (2004). It is not copyrighted, and it is included here with the author's permission. The FSL team allowed the text strings (help, warning and error messages) to be copied verbatim. Taylor Hanayik of the FSL group provided pseudo code for functions with little available documentation. The PolygoniseCube function comes from Cory Bloyd's public domain [Marching Cubes example](http://paulbourke.net/geometry/polygonise/) program. Jesper Andersson wrote the bwlabel.cpp file and explicitly allowed it to be shared under the BSD 2-Clause license. Jouni Malinen wrote the [high performance](https://github.com/gaspardpetit/base64) base64.cpp, distributed under the BSD license. [Sven Forstmann](https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification) wrote the mesh simplification, distributed under the MIT license. Chris Rorden ported it from C++ to C. Cameron Hart wrote [radixsort.c](https://github.com/bitshifter/radixsort) (2014) under the zlib license.

The `-romeo` phase-unwrapping command (`src/romeo.c`) is a C port of [ROMEO.jl](https://github.com/korbinian90/ROMEO.jl) and the [MriResearchTools.jl](https://github.com/korbinian90/MriResearchTools.jl) helpers its command-line app uses, by Korbinian Eckstein, Barbara Dymerska and Simon Robinson, together with the 2π range reduction from the Julia standard library. All are distributed under the MIT license. The upstream copyright and permission notices are preserved verbatim in `src/romeo.LICENSE`. Unlike the GPL module below, `-romeo` is compiled in **by default**, so a standard niimath binary contains this MIT-licensed component. MIT is compatible with the 2-Clause BSD License, so the binary as a whole remains BSD-2-Clause. `ROMEO=0 make` (or `-DENABLE_ROMEO=OFF`) omits it.

The `--medic` and `-unwarp` commands (`src/medic.c`) are original BSD-2-Clause code by the niimath authors. They are a clean-room emulation of the MEDIC method published by Van et al. (*Imaging Neuroscience* 4, 2026, [doi:10.1162/IMAG.a.1262](https://doi.org/10.1162/IMAG.a.1262)), developed from the paper and from black-box measurement of the reference tool's public executables. No reference implementation, test, build product or debug symbol was read, and no code from it is included. The measurements that fix each convention are recorded in the `medic_bench` repository. Phase unwrapping uses the MIT-licensed `-romeo` port described above, so `--medic` requires it (`ROMEO=0` implies `MEDIC=0`). `MEDIC=0 make` or `-DENABLE_MEDIC=OFF` omits MEDIC alone.

The `-moco` and `-stc` commands (`src/moco.c`, `src/stc.c`) are original BSD-2-Clause code by the niimath authors. Both emulate a published AFNI method whose reference implementation is copyrighted by the Medical College of Wisconsin: `3dvolreg`, `mri_3dalign`, `thd_rot3d` and `thd_shear3d` for `-moco`, and `3dTshift` and its FFT for `-stc`. Those files were **GPL-2** when this code was written. On 12 May 2026 MCW relicensed its 1994-2000 AFNI code to **CC BY 4.0**, which removes the copyleft bar but adds attribution and change-notice duties. Those sources were **not** read, translated or paraphrased. They served only as black-box oracles. The clean-room specification is the published method (Cox & Jesmanowicz 1999 for `-moco`; AFNI's published `3dTshift -help` and `-verbose` output for `-stc`) together with measured inputs and outputs, recorded in the `moco_bench` repository (`test/moco_reference_manifest.md` and `test/stc_reference_manifest.md`). The FFT in `stc.c` is original niimath code, a batched Stockham autosort kernel. No FFT implementation was read or adapted. Neither command carries any attribution obligation as a result. A binary containing these commands remains BSD-2-Clause.

The `-fugue` and `-fmapprep` commands (`src/fmap.c`) are original BSD-2-Clause code by the niimath authors. They are a clean-room implementation of the observable behavior of FSL's `fugue` and of `fsl_prepare_fieldmap SIEMENS`. FSL is licensed under the University of Oxford's non-commercial licence, which is incompatible with BSD-2-Clause. FSL's sources for fugue, prelude and `fsl_prepare_fieldmap` were not read, grepped, translated or paraphrased. The executables served only as black-box oracles, driven with synthetic inputs built to isolate one convention at a time. The published basis is Jezzard & Balaban, *Magnetic Resonance in Medicine* 34:65-73 (1995), and Jenkinson, *Magnetic Resonance in Medicine* 49:193-197 (2003) for PRELUDE. Phase unwrapping in `-fmapprep` uses the MIT-licensed `-romeo` port described above, so `ROMEO=0` drops `-fmapprep` while keeping `-fugue`. `FMAP=0 make` or `-DENABLE_FMAP=OFF` omits both.

The optional `-skullstrip` command (`src/skullstrip.c`) adapts public-domain AFNI code by Robert W. Cox and colleagues (NIMH): spatial normalization from `thd_brainormalize.c` and `thd_automask.c`, and the surface deformation and touchup stages from `SUMA_BrainWrap.c`. That file carries no copyright notice, so it falls under AFNI's US-Government-work clause, and a US Government work is not copyrightable (17 U.S.C. §105). AFNI states this affirmatively rather than by implication. Its `LICENSE.txt` declares the tree a "United States Government Work" apart from a listed set of exceptions. It also states that "contributions without explicit licensing will be assumed to be entered into the public domain". Its `README.copyright` dates the rule to work after 15 January 2001, and the adapted files' first commits are 2001-2004, by NIH authors. AFNI's `SUMA_3dedge3`, which wraps Malandain's GPL-3.0 `Extract_Gradient_Maxima_3D`, is deliberately out of scope, which is why niimath implements only the `-no_use_edge` behavior. One item was resolved rather than argued. The nearest-neighbor index conversion was originally written after reading `THD_3dmm_to_3dind_warn` in AFNI's `thd_coords.c`, which carries an MCW copyright header and predates the public-domain cutoff. It has been replaced by a clean-room reimplementation, `ss_world_to_index`. An implementer who had read neither AFNI nor niimath reproduced a 19,139-row table of measured input and output behavior. Working only from that table, they produced masks that are byte-identical. The protocol, the table, the sufficiency check and the implementer's derivation account are in the `skullstrip_bench` repository's `clean_room/` directory. The relicense has not changed this: a clean-room result carries no attribution duty, where adapting the CC BY original would. The surface primitives (icosphere, adjacency, intersection testing, rasterization) are original BSD-2-Clause code. A binary containing this command remains BSD-2-Clause. `-skullstrip` is off by default (`SKULLSTRIP=1 make`, `cmake -DENABLE_SKULLSTRIP=ON`) and is absent from every released binary.

The optional `-spm_coreg` and `-spm_deface` commands are the project's entire copyleft payload. They are carried in the [niimath_gpl](https://github.com/rordenlab/niimath_gpl) submodule at `src/GPL` and are enabled only when built with `make GPL=1` or `cmake -DENABLE_GPL=ON` (`-DHAVE_GPL`). They link SPM's `spm_coreg` module, which is **GPL-2 or later**, so distribute a `GPL=1` binary under those terms. The default build contains none of this code and remains BSD-2-Clause. The version string reported by `niimath` ends in ` GPL` or ` BSD` to show which applies, and `release_smoke.py` asserts both directions: a BSD binary that still carried `-spm_coreg` would fail the release gate.

## Links

- [imbibe](https://github.com/jonclayden/imbibe) is an R wrapper for niimath. It gives the performance of tuned code with the convenience of a scripting language.
- [3dcalc](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html) is AFNI's tool for image arithmetic.
- [c3d](https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md) provides mathematical functions and format conversion for medical images.
- [fslmaths](https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro3/index.html) is the inspiration for niimath.

## Citation

- Rorden C, Webster M, Drake C, Jenkinson M, Clayden JD, Li N, Hanayik T ([2024](https://apertureneuro.org/article/94384-niimath-and-fslmaths-replication-as-a-method-to-enhance-popular-neuroimaging-tools)) niimath and fslmaths: replication as a method to enhance popular neuroimaging tools. Aperture Neuro. 4. doi:10.52294/001c.94384

# niimath

[![Build status](https://ci.appveyor.com/api/projects/status/7o0xp2fgbhadkgn1?svg=true)](https://ci.appveyor.com/project/neurolabusc/niimath)

## About

It is said that `imitation is the sincerest form of flattery`. This project emulates the popular [fslmaths](https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro3/index.html) tool. fslmaths is a `general image calculator` and is not only one of the foundational tools for FSL's brain imaging pipelines (such as [FEAT](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT)), but has also been widely adopted by many tools. This popularity suggests that it fulfills an important niche. While scientists are often encouraged to discover novel solutions, it sometimes seems that replication is undervalued. Here are some specific reasons for creating this tool:

1. While fslmaths is provided without charge, it is not [open source](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence). This limits its inclusion in other projects, in particular for commercial exploitation.
2. Using an open source license allows niimath to build with open source libraries that the FSL team can not use. Specifically, an accelerated zlib ([zlib-ng](https://github.com/zlib-ng/zlib-ng), the release baseline) provides dramatically faster performance than the public domain library used by fslmaths. n.b. We previously helped update the [CloudFlare zlib](https://github.com/cloudflare/zlib/pull/19) fork that allows recent FSL releases to use an accelerated library, improving the speed for all FSL tools; niimath's releases now default to zlib-ng, which is also fast on arm64 and correct on Windows.
3. Minimal dependencies allow easy distribution, compilation and development. For example, it can be compiled for MacOS, Linux and Windows (fsl can not target Windows).
4. Designed from ground up to optionally use parallel processing (OpenMP and CloudFlare-enhanced [pigz](https://github.com/madler/pigz)).
5. Most programs are developed organically, with new features added as need arises. Cloning an existing tool provides a full specification, which can lead to optimization. niimath uses explicit single and double precision pipelines that allow the compiler to better use advanced instructions (every x86_64 CPU provides SSE, but high level code has trouble optimizing these routines). The result is that modern compilers are able to create operations that are limited by memory bandwidth, obviating the need for [hand tuning](https://github.com/neurolabusc/simd) the code.
6. Developing a robust regression testing dataset has allowed us to discover a few edge cases where fslmaths provides anomalous or unexpected answers (see below). Therefore, this can benefit the popular tool that is being cloned.
7. While the code is completely reverse engineered, the FSL team has been gracious to allow us to copy their error messages and help information. This allows true plug in compatibility. They have also provided pseudo code for poorly documented routines. This will allow the community to better understand the actual algorithms.
8. This project provides an open-source foundation to introduce new features that fill gaps with the current FSL tools (e.g. unsharp, sobel, resize functions). For future releases, Bob Cox has graciously provided permission to use code from [AFNI's](https://afni.nimh.nih.gov) 3dTshift and 3dBandpass tools that provide performance unavailable within [FSL](https://neurostars.org/t/bandpass-filtering-different-outputs-from-fsl-and-nipype-custom-function/824). Including them in this project ensures they work in a familiar manner to other FSL tools (and leverage the same environment variables). Note that the slice-timing command niimath actually ships, `-stc`, did **not** draw on that permission: it is a clean-room BSD-2 implementation that contains no AFNI code, written when `3dTshift.c` and its FFT were GPL-2. MCW has since relicensed its 1994-2000 AFNI code to CC BY 4.0 (12 May 2026), so that bar is gone — but CC BY still requires attribution and a statement of changes, which original code does not, so `-stc` stays as it is. See the License section.

The Reason to use fslmaths instead of niimath:

1. niimath is new and largely untested software. There may be unknown corner cases where produces poor results. fslmaths has been used for years and therefore has been battle tested. In the few instances where fslmaths generates results that bear no resemblance to its own documentation (as described below), one could argue it is the `correct` result (with comparison to itself). However, many tools may have been developed to assume this loss of high frequency signal and these tools may not perform well when provided with the result specified in the documentation.

## Installation

You can get niimath using several methods:

 - (Recommended) Download latest compiled release from [Github release web page](https://github.com/rordenlab/niimath/releases).
 - (Recommended) Download latest compiled release from [PyPI](https://pypi.org/project/niimath/):
  * `pip install niimath`
 - (Recommended) You can also download from the command line for Linux, MacOS and Windows:
  * `curl -fLO https://github.com/rordenlab/niimath/releases/latest/download/niimath_lnx.zip`
  * `curl -fLO https://github.com/rordenlab/niimath/releases/latest/download/niimath_macos.zip`
  * `curl -fLO https://github.com/rordenlab/niimath/releases/latest/download/niimath_win.zip`
 - (Developers) Download the source code from [GitHub](https://github.com/rordenlab/niimath), the next section describes how to build the software.

## Compilation

### CMake (recommended)

The easiest way to build niimath on a Unix computer is to use cmake. OpenMP is enabled by default (used by affine registration and optionally by core operations):

```
git clone https://github.com/rordenlab/niimath.git
cd niimath; mkdir build; cd build; cmake ..
make
```

On macOS, OpenMP requires Homebrew's libomp: `brew install libomp`. To disable OpenMP, use `cmake -DUSE_OPENMP=OFF ..`. Optional zstd compression support is auto-detected; install with `brew install zstd` (macOS) or `apt install libzstd-dev` (Linux).

Likewise, if you are compiling on Windows using cmake:

```
git clone https://github.com/rordenlab/niimath.git
cd niimath & mkdir build & cd build & cmake ..
cmake --build .
```

### Makefile (alternative)

You can compile the software by running the terminal command `make` from the project's `src` folder. This works with both Clang/LLVM and gcc on Linux and macOS:

```
git clone https://github.com/rordenlab/niimath.git
cd niimath/src
make
```

The default build includes OpenMP for all operations. On macOS this requires `brew install libomp`. To disable OpenMP, use `OMP=0 make`. Other Makefile options:

```
OMP=0 make             # Disable OpenMP
# zlib: CMake release builds default to zlib-ng (ZLIB_IMPLEMENTATION); plain make uses system -lz
make debug             # Debug build (-g, no optimization)
make ubsan             # Lightweight undefined-behavior checks (OpenMP-safe on macOS)
make sanitize          # AddressSanitizer build
AL=0 make              # Disable allineate registration
SKULLSTRIP=1 make      # Enable AFNI-style surface skull stripping (-skullstrip); OFF by default
ZSTD=0 make            # Disable zstd compression support
make wasm              # Emscripten/WebAssembly target
make wasm-wasi         # Experimental zlib-free WASI compute backend (needs Zig)
make wasm-emcc-core    # Feature-matched zlib-free Emscripten build (WASI benchmark baseline)
```

You can also compile this project to Web Assembly so it can be embedded in a web page, as shown in the [live demo](https://niivue.github.io/niivue-niimath/).

#### Optional copyleft module (`-spm_coreg`, `-spm_deface`)

Two commands live in the separate [niimath_gpl](https://github.com/rordenlab/niimath_gpl) submodule at `src/GPL` and are OFF by default, so ordinary builds stay BSD-2-Clause: `-spm_coreg` (SPM rigid-body coregistration) and `-spm_deface` (SPM-based defacing), both **GPL-2-or-later**. To build them, initialize the submodule and pass `GPL=1` (Makefile) or `-DENABLE_GPL=ON` (CMake). Run both commands from the repository root so the paths resolve:

```
git submodule update --init src/GPL    # or clone with --recurse-submodules
make -C src GPL=1
```

A binary built this way is a copyleft combined work (its version string ends in ` GPL`), covered by SPM's own terms: distribute it under the **GNU GPL-2 or later**. Without the module, `-spm_coreg`/`-spm_deface` report a clear error and the build stays BSD-2 (` BSD`). The GPL module computes only the rigid transform; niimath's BSD code applies it (reslicing and mask warping shared with `-allineate`/`-deface`).

When built with OpenMP, `-spm_coreg`/`-spm_deface` parallelize a single registration (the per-evaluation histogram and smoothing); results agree across thread counts to within the SPM golden tolerance (the histogram reduction sums in a thread-count-dependent order, so it is not guaranteed bit-identical across different team sizes, though it is deterministic at a fixed thread count). Control threads with `OMP_NUM_THREADS=N` or `-p N` (place `-p` before `-spm_coreg`). For batch runs of many subjects, prefer one subject per process with `OMP_NUM_THREADS=1` rather than threading each registration.

### Windows (command line)

For Windows, using the cmake method described above is highly recommended. However, you can also compile the project directly from the command line (here without the `-DHAVE_ZLIB` directive, so gz files will not be supported):

```
cl /Feniimath niimath.c core.c tensor.c bwlabel.c core32.c core64.c fdr.c meshify.c MarchingCubes.c quadric.c base64.c radixsort.c unifize.c nifti_io.c -DNII2MESH
```

### Linux universal binary

Simply running `make` in the `src` folder should compile niimath on Linux. However, the resulting executable will only work with specific versions of Linux. If you want to make a universal Linux release you can use [holy-build-box](https://github.com/FooBarWidget/holy-build-box). Be aware that this uses an old version of the gcc compiler (4.8.5), so the resulting performance may not be optimized for your system.

```
git clone https://github.com/rordenlab/niimath
sudo docker run -t -i --rm  -v `pwd`:/io ghcr.io/foobarwidget/holy-build-box-x64 /hbb_exe/activate-exec bash
cd /io/niimath/src
make
exit
sudo chown $(whoami) ./niimath/src/niimath
```

## JavaScript/WebAssembly

To read the WASM specific README, please click [here](./js/README.md). The rest of this README is for the `niimath` CLI program.

## Usage

niimath provides the same commands as fslmaths, so you can use it just as you would fslmaths. If you are brave, you can even rename it fslmaths and use it as a drop in replacement. You can also modify your environment variables to unleash advanced features:

 - Just like fslmaths, it uses your [`FSLOUTPUTTYPE` Environment Variable ](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslEnvironmentVariables) to determine output file format. Unix users can specify `export FSLOUTPUTTYPE=NIFTI_GZ`, `export FSLOUTPUTTYPE=NIFTI`, or `export FSLOUTPUTTYPE=NIFTI_ZST` (zstd compressed, requires zstd support) from the command line or profile. Windows users can use `set` instead of `export`.
 - To turn on parallel processing and threading, you can either set the environment variable `export AFNI_COMPRESSOR=PIGZ`. If the environment variable `AFNI_COMPRESSOR` does not exist, or is set to any value other than `PIGZ` you will get single threaded compresson.

niimath has a few features not provided by fslmaths:

 - `bptfm <hp> <lp>`        : Same as bptf but does not remove mean (emulates fslmaths < 5.0.7)
 - `bwlabel <conn>`         : Connected component labelling for non-zero voxels (conn sets neighbors: 6, 18, 26)
 - `ceil`                   : round voxels upwards to the nearest integer
 - `crop <tmin> <tsize>`    : remove volumes, starts with 0 not 1! Inputting -1 for a size will set it to the full range
 - `dehaze <mode>`          : set dark voxels to zero (mode 1..5; higher yields more surviving voxels)
 - `detrend`                : remove linear trend (and mean) from input
 - `demean`                 : remove average signal across volumes (requires 4D input)
 - `edt`                    : estimate Euler Distance Transform (distance field). Assumes isotropic input
 - `floor`                  : round voxels downwards to the nearest integer
 - `mod`                    : modulus fractional remainder - same as '-rem' but includes fractions
 - `otsu <mode>`            : binarize image using Otsu''s method (mode 1..5; higher yields more bright voxels))
 - `power <exponent>`       : raise the current image by following exponent
 - `resize <X> <Y> <Z> <m>` : grow (>1) or shrink (<1) image. Method <m> (0=nearest,1=linear,2=spline,3=Lanczos,4=Mitchell)
 - `robustfov [mm]`         : crop to a robust field of view (default 170mm) from the top of the head down, removing lower head/neck (emulates FSL robustfov); adjusts dim and sform/qform
 - `round`                  : round voxels to the nearest integer
 - `sobel`                  : fast edge detection
 - `sobel_binary`           : sobel creating binary edge
 - `tensor_2lower`          : convert FSL style upper triangle image to NIfTI standard lower triangle order
 - `tensor_2upper`          : convert NIfTI standard lower triangle image to FSL style upper triangle order
 - `tensor_decomp_lower`    : as tensor_decomp except input stores lower diagonal (AFNI, ANTS, Camino convention)
 - `trunc`                  : truncates the decimal value from floating point value and returns integer value
 - `unsharp  <sigma> <scl>` : edge enhancing unsharp mask (sigma in mm, not voxels; 1.0 is typical for amount (scl))
 - `dog <sPos> <sNeg>`      : difference of gaussian with zero-crossing edges (positive and negative sigma mm)
 - `dogr <sPos> <sNeg>`     : as dog, without zero-crossing (raw rather than binarized data)
 - `dogx <sPos> <sNeg>`    : as dog, zero-crossing for 2D sagittal slices
 - `dogy <sPos> <sNeg>`    : as dog, zero-crossing for 2D coronal slices
 - `dogz <sPos> <sNeg>`    : as dog, zero-crossing for 2D axial slices
 - `mesh`                  : see separate section below
 - `qform <code>`          : set qform code
 - `sform <code>`          : set sform code
 - `unifize [-GM]`         : bias field correction (adapted from AFNI 3dUnifize); optional `-GM` also scales gray-matter intensity toward a common target
 - `allineate <base> [opts]`: affine registration to match 'base' (from AFNI 3dAllineate)
   - opts: `-cost XX` (`fast`/`fastx` [adaptive default], `fasthel`, `fastcr`, `hel`, `nmi`, `lpc`, `lpa`, `ls`) `-cmass` `-nocmass` `-source_automask`
   - `-warp XX` (sho, shr, srs, aff [default]) `-interp XX` (NN, linear [default], cubic)
   - `-final XX` (NN, linear, cubic [default]) or `-nearest` `-linear` `-cubic`
   - `-fill XX` (auto [default], zero, nan): out-of-FOV output fill. `auto` uses 0, or the source's darkest voxel if it is negative (CT/Hounsfield air ≈ -1000, so out-of-FOV reads as air rather than soft tissue); `zero` is always 0 (byte-identical to the historical behavior for positive-only MRI); `nan`. Not available for `deface`, whose fill stays the image minimum
   - `-master <grid>`: estimate at the base resolution but reslice the result onto `<grid>` (must share the base world frame, e.g. a higher-resolution template)
   - `-savemat out.json`: save the fitted world-space `fixed_to_moving` affine (and its inverse) as self-describing JSON. `-applymat in.json`: reslice the moving image onto `base` using a saved affine, doing **no** registration (exclusive with the registration/seed options; use `-nearest` for label/atlas volumes)
   - Registration seeds (applied to the moving image before the fit): `-com` (reset the origin to the brightness center of mass), `-sym`/`-symd`/`-symb` (fold a midsagittal-plane correction into the header; `-symd` de-obliques the frame first, `-symb` auto-competes both), `-nosagseed` (disable the in-MSP rigid seed `-sym` runs by default), `-zoom` (relax the scale range for abnormal-size subjects, e.g. infant vs adult template)
   - Skull-stripping here is `deface` with a brain mask (template-driven); `-skullstrip` (surface-based, `SKULLSTRIP=1`) does it with no template at all. To crop the field of view first, chain `-robustfov` before `-allineate`. Together these make niimath a superset of the standalone `allineate` registration tool.
   - **Fast engine (the default):** a bare `-allineate`, `-cost fast`, and `-cost fastx` all select the adaptive clean-room estimator. On a whole-head base it independently fits HEL and correlation-ratio coarse candidates, selects by HEL dependence×overlap, and continues once through the finer HEL stages. A hard-zeroed/skull-stripped base activates the deeper rigid-HEL, scale-bracketed-HEL, and CR-seeded multi-start, with fine-level arbitration before final affine polish. Use `-cost fasthel` to force the HEL-only fast trajectory or `-cost fastcr` to force correlation-ratio only; use `-cost hel`/`nmi`/`lpc`/`lpa`/`ls` for the ordinary AFNI-style engine. The fast estimator is SPM/FLIRT-inspired, multiresolution, and typically several times faster. It runs a fixed schedule with internal sampling, so it rejects `-warp`/`-interp`/`-source_automask`/`-dark_automask`/`-zoom` (use an ordinary cost for those) — but the `-com`/`-sym`/`-symd`/`-symb` header seeds work with it. Default/`-cmass` chooses the supplied-affine or COM-recentered initialization from dependence×overlap, `-com` forces COM, and `-nocmass` forces the supplied affine; `-final`/`-master`/`-savemat` are also honored. `-cost` is last-one-wins. `-weight <img>` (a base-space region image, 3D with dims and world frame matching `base`) is an AFNI 3dAllineate-style **graded** weight: it is normalized to `[0, 1]` (divide by max) and applied per base voxel — a voxel weighted 0 is excluded, one near 1 dominates. It is **not an exclusion mask**; keep the out-of-ROI head attenuated (nonzero), because zeroing the whole exterior leaves global scale underdetermined and the cross-modal fit collapses into the scalp (an AFNI-style whole-head weight anchors absolute size). BOTH engines honor it: the ordinary engine loads it internally in place of its manufactured autoweight, while the fast engine applies it only at the finest 2 mm stage so coarse capture and global-scale selection stay whole-head. `-weight` is rejected with stdin (`-`) and `-applymat`. A background-only weight is rejected. If the implicit default fast engine fails and falls back to ordinary Hellinger, the weight remains honored; explicit fast selectors never silently switch engines.
 - `deface <tmpl> <mask> [opts]`: remove voxels using a template-space mask. Registers the input to `tmpl` (affine), inverts the transform, warps `mask` onto the input's native grid, and sets voxels where the warped mask < 0.5 to the image's finite minimum (≈0 for typical MRI; the input itself is never resampled). `mask` is in `tmpl` space: ≥0.5 = keep, <0.5 = remove. The mask determines what is removed — supply a brain mask to skull-strip (keep the brain) or a face mask to deface (remove the face).
   - opts: `-cost XX` (`fast`/`fastx` [adaptive default], `fasthel`, `fastcr`, `hel`, `nmi`, `lpc`, `lpa`, `ls`) `-cmass`/`-nocmass` tuning + `-final`/`-nearest`/`-linear`/`-cubic` (default final: linear). Like `-allineate`, deface defaults to the adaptive fast engine; `-cost fasthel`/`fastcr` force one fast cost and `-cost hel` selects the ordinary AFNI-style engine. Fast cannot honor `-warp`/`-interp`/`-source_automask`/`-dark_automask` (use `-cost hel`). The seed/matrix/master workflow options (`-savemat`/`-applymat`/`-com`/`-sym`/`-symd`/`-symb`/`-nosagseed`/`-zoom`/`-master`) apply to `-allineate` only and are **rejected** here
   - **Breaking change, and the name has since been re-used:** the former `-skullstrip <tmpl> <mask>` command ran this identical operation and was removed — replace it with `-deface`, supplying your brain mask, and the result is unchanged. `-skullstrip` now names a **different** operation, template-free surface skull stripping (see below), which takes **no arguments**. An old two-argument `-skullstrip` line is now **detected and rejected before any work is done**, with a message naming `-deface`; it does not run the new operation and does not misread your template and mask. Migrate to `-deface` explicitly.
 - `reface <tmpl> <shell> <weight> [opts]`: **anonymize** by face replacement (emulates AFNI `afni_refacer2 -mode_reface`). Registers the subject to `tmpl`, back-projects the signed template-space face-replacement `shell` onto the **original subject grid**, and composites an anonymized image (shell > 0 → the shell scaled by a brightness-match factor, shell == 0 → keep the subject, shell < 0 → zero), with an edge blend inside the replaced region. Output stays on the subject grid. Works with both engines and honors `-weight` (the `weight` argument is **required** — it is reused as the registration weight); opts are the `-cost` tuning as for `deface` (the shell back-projection is always nearest-neighbour, so `-final` does not apply). For privacy the coverage diagnostic **fails closed**: if little of the shell mapped into the subject FOV (<10%, a likely registration failure) `-reface` refuses to write the possibly-unanonymized image and returns an error.
 - `qwarp <base>` (built only with `QWARP=1`): **nonlinear (deformable) registration** to `base`, an attributed port of AFNI `3dQwarp -blur 0 3`. It takes a single `base` argument and has **no sub-options** of its own, but it is an ordinary chain operation — it produces a warped image on the `base` grid, so you may chain further ops after it. The input must **already** be unifized, skull-stripped, affine-aligned, and share the `base` grid (same dims + world frame). Off by default (memory/CPU-heavy, impractically slow in WebAssembly) and depends on the allineate engine, so it is unavailable in an `AL=0` build.
 - `romeo <mag|none> [opts]`: **ROMEO phase unwrapping** — a faithful MIT-licensed C port of [ROMEO.jl](https://github.com/korbinian90/ROMEO.jl) plus the [MriResearchTools.jl](https://github.com/korbinian90/MriResearchTools.jl) helpers its command-line app uses. The chain image is the wrapped phase (3D, or 4D with echoes on dim 4); the magnitude is a **required positional argument** — pass the literal `none` to unwrap without one. An ordinary chain operation, so further niimath operations may follow. Enabled by default; `ROMEO=0 make` / `-DENABLE_ROMEO=OFF` omits it.
   - opts: `-t <TEs>` (echo times in ms: `16.8`, `16.8,38.56`, `'[16.8,38.56]'`, `epi [te]` — quote the bracket form for the shell; **required for multi-echo input**, optional for a single echo) `-k nomask|robustmask|qualitymask [thr, default 0.1]|<mask file>` `-w romeo|romeo2|romeo3|romeo4|romeo6|<up to 6 bits e.g. 1010>` (bare `romeo` resolves to `romeo3` when a magnitude is supplied and `romeo4` when it is not) `-template <n>` `-i` (individual, not temporal) `-temporal-uncertain-unwrapping [x]` `-g` (correct global n2π offset) `-q`/`-Q` (quality maps) `-B [name]` (B0 field map in Hz) `-B0-phase-weighting phase_snr|phase_var|average|TEs|mag|simulated_mag` `-no-phase-rescale` (alias `-no-rescale`) `-no-mask-out` `-v`
   - side outputs use `nifti_save` postfixes on the output name: `<out>_mask` (only when a mask was actually computed — `-k nomask` writes none), `<out>_quality` (`-q`), `<out>_quality_1..6` (`-Q`; a map that is uniformly 1.0 in the interior is skipped, as upstream), `<out>_B0` and `<out>_B0_snr` (`-B`); they honor `FSLOUTPUTTYPE`, and `-gz` **when it precedes `-romeo`** — the companions are written during the operation, so a later `-gz` only reaches the main output
   - the phase is rescaled to `[-π, π]` following `readphase` unless `-no-phase-rescale` is given; because that inspects the *unscaled stored* values, `-romeo` must be the **first** computational operation when rescaling is active
   - `-B` computes B0 **without** MCPC-3D-S phase-offset correction, which ROMEO's own multi-echo `-B` silently enables; the maps therefore correspond to `romeo --compute-B0 --phase-offset-correction off`, and niimath says so on stderr rather than implying full equivalence. Without a magnitude, ROMEO's SNR map collapses to a single value (the substituted `exp(-TE/20)` decay is voxel-independent); niimath writes that same constant across the working grid instead.
   - **Not yet ported** (rejected with a specific message rather than silently ignored): `-u`, `-e`, `-threshold`, `-w bestpath`, `-max-seeds > 1`, `-merge-regions`, `-correct-regions`, `-wrap-addition != 0`, `-fix-ge-phase`. MCPC-3D-S phase-offset correction and multi-channel (5D) input are out of scope.
   - please cite Dymerska, B. et al. 2020, *Magnetic Resonance in Medicine*, [doi:10.1002/mrm.28563](https://doi.org/10.1002/mrm.28563)
 - `unwarp <map> <axis>`: **EPI distortion correction** — resample the chain image through a scalar displacement map in millimetres (as written by `--medic`, see below). An ordinary chain operation, so further niimath operations may follow. Enabled by default with `--medic`; `MEDIC=0 make` / `-DENABLE_MEDIC=OFF` omits it.
   - `<map>` is 3D (broadcast over every frame) or 4D matching the input's frame count, and must share the input's dimensions and world transform
   - `<axis>` is `i`, `j` or `k` (equivalently `x`, `y`, `z`). A trailing `-` is accepted and **ignored**: the sign already lives in the stored map, so negating again would double-correct. `--medic --phase-encoding-direction` is deliberately the opposite — it **honors** the suffix — so pass the full BIDS value there and do not worry about it here
   - resampling uses an unnormalized Lanczos-5 windowed sinc applied separably in 3D, with zero fill outside the field of view and no Jacobian intensity modulation — matching the reference implementation's measured behavior (measured; see the `medic_bench` repository)
 - `moco [-ref <n|image>] [-1Dfile <path>]`: **rigid-body motion correction** — registers every volume of a 4D series onto a reference and replaces the image with the corrected series: `niimath bold -moco out`. The reference is volume 0 unless `-ref` names another one: `-ref 7` uses volume 7 of the series (`niimath bold -moco -ref 7 out`), while `-ref sbref.nii.gz` registers onto an **external reference image** such as a single-band reference (`niimath bold -moco -ref sbref.nii.gz out`). An all-digit argument is a volume number and anything else is a filename, so a file whose name is only digits must be written as `./7`; an external reference must lie on the input's voxel grid (identical dimensions and a voxel-to-world transform agreeing within 0.001 mm) and is otherwise rejected rather than silently resliced, and a 4D one contributes its volume 0. With a reference drawn from the series that volume is copied through unchanged with an all-zero parameter row; with an external reference every volume is registered. Add `-1Dfile` to also write the six motion parameters per volume: `niimath bold -moco -1Dfile out.1D out`; the parameter filename must end in `.1D`. For quality control there is `-relative`, which measures each volume against the **previous** one rather than against a reference: `niimath bold -moco -relative rel.1D`. That is a measurement-only pass — nothing is registered or resampled and **no image is written at all**, so the trailing filename that would normally name the output image is the parameter file itself and must end in `.1D` (naming an image there is refused rather than silently producing a `.1D` called `out.nii.gz`). `-ref` and `-1Dfile` are rejected with it, since there is no reference fit to report. Each pair is fit against the *original* predecessor, so the numbers are raw frame-to-frame motion (the usual input to framewise-displacement style QC) rather than residual drift after correction, and row 0 is zero because volume 0 has no predecessor. It writes two files and nothing else: the text `.1D` at `%12.8f`, and a float64 companion at `<name>.1D.bin` holding the same `nt`×6 numbers as raw little-endian doubles, row-major and header-free, so `numpy.fromfile('rel.1D.bin').reshape(-1, 6)` reads it exactly. Rebuilding the registration weight and derivative images for every pair makes this roughly an order of magnitude more work per volume than ordinary correction. The input must be 4D with more than one volume (correction runs in float32, so `-dt double` is rejected); an ordinary chain operation, so further niimath operations may follow. A clean-room BSD-2 implementation of the method of Cox & Jesmanowicz (*Magnetic Resonance in Medicine* 42:1014-1018, 1999), the algorithm behind AFNI `3dvolreg`. Enabled by default on every platform, including the WebAssembly build; `MOCO=0 make` (or `-DENABLE_MOCO=OFF`) omits it.
   - the parameter file is compatible with AFNI's `-1Dfile`: six columns `roll pitch yaw dS dL dP`, one row per volume — rotations in degrees counter-clockwise about the I-S, R-L and A-P axes, shifts in mm toward Superior, Left and Posterior — recording the **correction** that was applied, not the estimated motion. Row 0 is all zeros (the base registers to itself)
 - `stc --slicetiming <t0,t1,...|@file> [-tzero <sec>]`: **slice-time correction** — shifts every voxel time series of a 4D series so that all slices share one temporal origin: `niimath bold -stc --slicetiming @times.1D out`. `--slicetiming` is required, case-sensitive, and must come first; it takes either one comma-separated list of slice acquisition times in **seconds** or an AFNI-style `@file` (whitespace- or comma-separated, `#` starts a comment). Supply exactly one value per slice along storage axis `k`, in slice-index order — a count that does not match `nz` is an error, not a silently truncated list. Optional `-tzero <sec>` sets the common time point (default: the arithmetic mean of the supplied times) and must lie within their `[min, max]`. The input must be a scalar 4D image with at least 5 volumes, a finite positive `pixdim[4]`, and a usable temporal unit in `xyzt_units` — seconds, milliseconds or microseconds; niimath will not assume seconds. Correction runs in float32, so `-dt double` is rejected; an ordinary chain operation, so further niimath operations may follow. Only `toffset` changes in the header (to the common time point, in the header's own time units); geometry, TR and the spatial transforms are untouched. A clean-room BSD-2 implementation of the default Fourier method (detrend → interpolate → retrend) of AFNI `3dTshift`, whose source was used only as a black-box oracle. Enabled by default on every platform, including the WebAssembly build; `STC=0 make` (or `-DENABLE_STC=OFF`) omits it.
   - v1 corrects along storage axis `k` only. `test/stc_slicetiming.py` is a standard-library-only helper that reads a BIDS sidecar and prints the `--slicetiming` argument, honoring `SliceEncodingDirection` (`k` passes through, `k-` is reversed into slice-index order, `i`/`j` are rejected) and hard-erroring if the sidecar's `RepetitionTime` disagrees with the unit-normalized header TR: `niimath bold.nii.gz -stc --slicetiming "$(python3 test/stc_slicetiming.py bold.nii.gz)" out.nii.gz`
 - `spm_coreg <ref> [opts]`: SPM rigid-body coregistration of the chain image to `ref` (optional GPL module, see below)
   - opts: `-cost XX` (nmi [default], mi, ecc, ncc, ls) `-sep` `-fwhm` `-dither 0|1` `-coarse sparse|downsample` `-verbose 0|1`
   - default reslices onto the `ref` grid (`-interp trilinear [default]|nearest`, `-fill zero [default]|nan`); `-estimate` instead rewrites only the source sform/qform
 - `spm_deface <tmpl> <mask> [opts]`: SPM analogue of `deface`, registering with spm_coreg (optional GPL module, see below)
   - opts: same estimate sub-options as `spm_coreg`, plus `-interp`
 - `--dtifit -k <dwi> -r <bvec> -b <bval> -o <base> [-m <mask>] [-xflip 0|1|auto]` : linear diffusion tensor fit (emulates FSL `dtifit`)
   - writes `<base>_{FA,MD,L1,L2,L3,V1,V2,V3,S0,MO,tensor}`; fit math from AFNI 3dDWItoDT (public domain)
   - `-xflip auto` (default) flips the bvec X component when the spatial transform determinant is positive, matching FSL
 - `--qc <t1> --seg <seg> --csf <i[,j..]> --wm <i[,j..]> [--erode 0|1] [--out qc.tsv]` : MRIQC-style anatomical quality metrics from a T1 + integer segmentation
   - writes a wide TSV (default `qc.tsv`) with CJV, cnr_noair, per-tissue/total SNR, WM2MAX, efc_brain, ICV fractions + mm³ volumes, and per-tissue summary stats
   - label convention: `0` = non-brain (excluded); `--csf`/`--wm` give disjoint CSF/WM label values, every other non-zero label is GM. Only air-free metrics are computed; this hard-segmentation variant uses unrounded intensities and NumPy-linear percentiles, so it is not numerically interchangeable with MRIQC's soft-PVM summaries. `cnr_noair`/`efc_brain` flag deviations from MRIQC norms
 - `--medic --magnitude <e1> [<e2> ...] --phase <e1> [<e2> ...] --te-ms <t1,t2,...> --total-readout-time <sec> --phase-encoding-direction <i|j|k|i-|j-|k-> --out-prefix <path>` : **MEDIC** multi-echo distortion correction — estimates a B0 field map per frame from multi-echo phase and converts it to an EPI displacement map you can apply with `-unwarp`
   - writes `<prefix>_fieldmaps_native` (Hz, distorted grid), `<prefix>_fieldmaps` (Hz, undistorted grid) and `<prefix>_displacementmaps` (mm), all float32. At least two echoes are required. Per frame it builds a tiered brain mask (core / border ring / outside), removes the MCPC-3D-S phase offset with a global 2π branch selection, unwraps with ROMEO, and applies a per-echo intra-frame 2π offset. Across frames it then applies a temporal 2π consistency correction **to the unwrapped phase**, fits a magnitude-weighted field map, and truncates that series to a low rank — the brain interior at `--rank`, the border ring smoothed spatially and rebuilt from far fewer components
   - give `--phase-encoding-direction` the BIDS value including its sign (`j-`, not `j`): `--medic` **honors** the `-` suffix — it drives the inversion and flips the sign of the displacement map. Note the contrast with `-unwarp`, which **ignores** the suffix because the sign is already stored in the map it reads
   - options: `--rank <N>` (low-rank truncation, default 10; `0` disables — it applies to the brain **interior**, the border ring being handled by `--border-regularization`) `--mask-mode <tiered|robustmask|mindgrab>` (default `tiered`; `robustmask` is the pre-2026 behaviour, kept for bisection; `mindgrab` shells out to the external `brainchop-mindgrab` and fails if it is absent) `--branch-correction <0|1>` (default 1; ROMEO's global 2π correction at both unwrap stages plus MCPC-3D-S branch selection) `--echo-offset <0|1>` (default 1; the intra-frame per-echo 2π offset) `--border-regularization <0|1>` (default 1; `0` gives the plain global truncation) `--temporal-correction <0|1>` `--phase-offset <mcpc|none>` `--noise-frames <N>`, `-f` (drop N trailing frames) `--weights <romeo|romeo2|romeo3|romeo4|romeo6>` (default `romeo4`; governs both the MCPC-3D-S and the multi-echo unwrap) `--mask <file>` (use this mask verbatim for both unwrapping stages, instead of the default tiered mask; **mutually exclusive with `--mask-mode`**) `--save-intermediates` `--n-cpus <N>`, `-n` `--gz <0|1>`. See also `niimath --medic --help`
   - each of `--mask-mode robustmask`, `--branch-correction 0`, `--echo-offset 0` and `--border-regularization 0` restores the pipeline as it stood immediately before that stage was added, which is how a change in output is attributed to a stage. There is no single combination that returns to the original, because the temporal grouping's switch from magnitude to unwrapped phase has no flag — it corrects what the statistic measures rather than exposing a choice. Note also that the grouping presupposes the global branch has been pinned, so `--branch-correction 0` leaves every frame alone in its group and the temporal correction inactive; niimath says so when that happens
   - a supplied `--mask` is binary, so it becomes the **core tier** with no border ring and `--border-regularization` has nothing to do. `--mask` counts a voxel as inside the brain when its value is **`>= 1`**, not merely non-zero — that is the measured convention of the reference. A fractional probability map is therefore not a mask: threshold it first (`niimath p.nii -thr 0.5 -bin mask.nii`). NaN is treated as outside, and a mask with no voxel `>= 1` is an error rather than an empty result
   - the whole run is held in RAM by design (a 4D `.nii.gz` cannot be seeked, so nothing is gained by streaming it): the work arrays are `nx·ny·nz × frames × (2·echoes + 3) × 4` bytes and that budget is printed at startup; the peak adds one echo pair of input (they are loaded and released one echo at a time) and one echo pair in transit — measured peak for 170 frames × 2 echoes at 76×76×46 is 1.54 GB single-threaded and 2.04 GB at 8 threads writing uncompressed (the reference tool needs 3.50 GB and 4.21 GB for the same run). The default mask and border stages add allocations the printed budget does not cover — a per-frame magnitude mask, a per-thread morphology workspace, and the border filter's own buffers — which is most of the gap between the banner and the 8-thread peak. It grows linearly with frames × echoes, so estimate natively and use `-unwarp` where memory is tight (a WebAssembly build has a 4 GiB ceiling)
   - an optional stdlib-only BIDS wrapper (`medic.py`, in the `medic_bench` repository) discovers multi-echo runs, reads the parameters from their JSON sidecars, and drives `--medic` and `-unwarp` for you
   - a clean-room emulation developed from the published method, black-box measurement of the reference tool, and two non-source documents supplied by the method's author (a theory note on phase-offset ambiguities, and correspondence); no reference source was read at any point. every convention it implements is recorded in the `medic_bench` repository, which also holds the benchmarks and the patent analysis. **No equivalence with the reference tool is claimed.** `-unwarp` does reproduce it closely — fed the reference's own displacement map it matches the reference's corrected images to nrmse 3.5e-5 — and `--medic` now agrees closely over the brain **interior** while differing at the border. Read the two separately: restricted to the core tier the native field maps agree at r 0.995 / p99 3.5 Hz and the displacement maps to 0.02 mm median / 0.3 mm p99 excluding folds, and the corrected images correlate 0.99991 end to end; including the three-voxel border ring the global field p99 is ≈ 60 Hz, because that is where the reference's own field is broadband. Quoting the global figure alone misrepresents both. The reference's brain-mask construction **is** now emulated (valid-tier Dice 0.9998) and so is its border-aware low-rank filtering; its iteration-limited field inversion is the one part deliberately not reproduced
   - please cite Van et al. 2026, *Imaging Neuroscience* 4, [doi:10.1162/IMAG.a.1262](https://doi.org/10.1162/IMAG.a.1262), and the ROMEO reference above for the unwrapping

 - `-skullstrip [-faithful]` : **AFNI-style surface skull stripping** — no template, no mask, no network. Expands a surface outward from inside the head until it wraps the brain (the method of AFNI `3dSkullStrip -no_use_edge`), then applies the resulting mask. Scalar 3D input only, computed in float32. An ordinary chain operation, so further operations may follow. **OFF by default, and 64-bit native only:** build with `SKULLSTRIP=1 make` or `cmake -DENABLE_SKULLSTRIP=ON`. The wasm, tiny and nano targets reject the request explicitly rather than silently dropping it — WebAssembly support is a stated goal of the design but is **not implemented**, so today this is a desktop-only command
   - output convention: in-mask voxels keep their **original intensities** and out-of-mask voxels become the image **minimum**. This is a deliberate divergence from AFNI's default, which rescales intensities
   - **the name is re-used.** It previously aliased the template+mask removal now spelled `-deface`. The new `-skullstrip` needs **no template and no mask** (`-faithful` is its only option) — see the note under `deface` above before migrating an old command line
   - `-skullstrip -faithful` runs the reference deformation kernel instead of the optimised default. Same algorithm and the same quality — the two differ only in how the ray walk is compiled and in whether the surface's node loop runs on several threads — but `-faithful` reproduces the pre-optimisation release **exactly**, which makes it the thing to diff against when checking for a regression. The default is about 1.8x quicker and is what you want otherwise. Across nine varied images on an M4 Pro: 17.4 s with `-faithful`, 9.6 s without. `-faithful` must come before the output filename
   - because the surface converges on a discrete test (an integer count of troubled nodes), a difference far below single-precision rounding can add a whole extra pass. Output is therefore **not** byte-stable between compilers, build systems, or these two kernels; compare masks with Dice, not `cmp`. Within one binary it is exactly reproducible, including across `-p` thread counts — the default kernel spreads the surface's node loop over threads and still gives bit-identical output at any team size, while `-faithful` runs that loop on one thread and so ignores `-p` for the part that dominates its time
   - this is automated research segmentation, **not a diagnostic guarantee**. Inspect the output
 - `--compare <ref>`       : report if images are identical, terminates without saving new image
 - `--bitmap -a name.png`  : mimic fsl slicer (see [niimath-bitmap](https://github.com/rordenlab/niimath-bitmap))
 - `filename.nii`          : mimic fslhd (can also export to a txt file: 'niimath T1.nii 2> T1.txt') report header and terminate without saving new image

## Identical Versus Equivalent Results

This project is designed to provide equivalent results to fslmaths. In most cases, the results are identical, virtually all others are equivalent. The results are not always identical as computations are conducted using floating point representations, where the precise order of instructions can generate small rounding differences. As [Kernighan and Plauger]( https://www.amazon.com/Elements-Programming-Style-Brian-Kernighan/dp/0070341990) note `Floating point numbers are like piles of sand; every time you move one you lose a little sand and pick up a little dirt.` Raw brain imaging data is typically stored as 16-bit integers (and the signal-to-noise is typically a fraction of this dynamic range), whereas niimath uses single (32-bit) or double (64-bit) floating point representations. Therefore, while niimath may generate results that are not identical, the results are intended to be always comparable. For further information on floating point accuracy, suggested readings include [here](https://introcs.cs.princeton.edu/java/91float/) and [here](http://www.freshsources.com/page1/page7/files/Sand-1.pdf).

This project includes the `--compare` argument that allows you to directly the results of niimath and fslmath. A validation repository is also available, which runs hundreds of commands to detect the quality of the output. The validation repository includes two scripts. The `batch.sh` script tests functions that generate identical results. The `close.sh` script conducts tests on functions that provide equivalent but not identical results. For example, for tensor decomposition the vector [1 0 0] is the functionally identical to [-1 0 0] as for fiber tracking the fiber direction ignores vector polarity. When a difference is detected by the `--compare` function, a report is generated allowing the user to determine the equivalence of solutions:

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

Some operations do generate known meaningfully different results. These are listed below, with the rationale for the discrepancy provided:

1. The command "fslmaths inputimg -add 0 outputimg -odt input" can convert a uint8 image float output despite explicit request to retain input type. This occurs if the input image header has a non-unitary scale slope or non-zero intercept. In contrast, niimath retains both the datatype and the intensity scaling parameters.
2. Different versions of fslmaths perform differently for the pass through "fslmaths in out" which is useful for copying files. Old versions will losslessly save in the input datatype, while fslmaths 6.0 converts the data to float. niimath retains the datatype.
3. The fslmaths function `-fillh26` will sometimes fill unconnected regions. An example has been provided to the FSL team. niimath provides the correct solution.
4. The fslmaths `-dilD` function does not do what it claims. It introduces a blurring effect that reduces edge artifacts that plague iterative morphology operations. Unfortunately, this effect is conducted in a consistent order that introduces a spatial shift in signal. In contrast, niimath does the dilation as described. Note there are [better solutions](https://github.com/neurolabusc/niiSmooth) for these functions. The niimath '-edt' operation can also be used for dilation.
6. The fslmaths `-roc` function works differently than described in the help. It appears to ignore voxels near the edge of an image and generates "given object has non-finite elements" if any dimension is less than 12 voxels. When provided with an external noise file, it generates additional columns in the output file that are not described. It does not seem to precisely detect the desired `AROC-thresh`, but samples at different stepped intervals. niimath attempts to emulate the stepped intervals for reporting, but determines the precise cutoff.
7. Be aware that fslmaths help suggests `If you apply a Binary operation (one that takes the current image and a new image together), when one is 3D and the other is 4D, the 3D image is cloned temporally to match the temporal dimensions of the 4D image.` This is not the case for -thr or -uthr: if the second item is 4D, only the first volume is used and the output remains 3D. Particularly odd is uthr: `fslmaths 3D -uthr 4D out` will fill input volume 3D with zeros, regardless of mask values.
8. Perhaps understandably, `fslmaths in1 -rem 0 out` will throw an exception. However, `fslmaths in1 -rem in2 out` will throw an exception if any voxel in the image `in2` is zero. While this seems understandable, niimath provides a description for this error.
9. The fslmaths function `-rem` returns the **integer** modulus remainder. This replicates the C `%` operator. This may be unexpected, e.g. in Python `2.7 % 2` is 0.7, as is Matlab's `mod(2.7, 2)`, as is standard C `fmod`. niimath clones the fslmaths  behavior, but also includes a new function `-mod` to return the modulus fractional remainder.
10. Be aware that fslmaths takes account of whether the image has a negative determinant or not (flipping the first dimension). However, fslstats does not do this, so fslstats coordinates are often misleading. For example, consider an image in RAS orientation, where the command `fslstats tfRAS -x` will give coordinates that are incompatible with fslmath's `tfceS` function. niimath attempts to emulate the behavior of fslmaths for the relevant functions (-index -roi, -tfceS).
11. Neither `-subsamp2` nor `-subsamp2offc` handle anti-aliasing. Be aware that `-subsamp2offc` can exhibit odd edge effects. The problem is simple to describe, for slices in the middle of a volume, and output slice is weighted 50% with the center slice, and 25% for the slice below and the slice above. This makes sense. However, bottom slices (as well as first rows, first columns, last rows, last columns, last slices) the filter weights 75% on the central slice and just 25% on the slice above it. Signal from this 2nd slice is heavily diluted. A better mixture would be 66% edge slice and 33% 2nd slice. This latter solution is used by niimath.
12. fslmaths 6.0.0..6.0.3 were unable to process files where the string ".nii" appears in a folder name. For example, consider the folder "test.niim", the command `fslmaths ~/test.niim/RAS -add 0 tst` will [generate an exception](https://github.com/FCP-INDI/C-PAC/issues/976). niimath will recognize that this is a folder name and not a file extension and work correctly. niimath helped detect this anomaly and it is an example of how a clone can help provide feedback to the developers of the original project.
13. The fslmaths function [`-ztop`](https://github.com/rordenlab/niimath/issues/8) fails to clamp extreme values.

Finally, it is possible that there are some edge cases where niimath fails to replicate fslmath. This is new software, and many of the operations applied by fslmaths are undocumented. If users detect any problems, they are encouraged to generate a Github issue to report the error.

## Superior Performance

Here are some examples of speed up factors you can expect. The sample T1-weighted and resting state data use the [HCP 3T Imaging Protocol](http://protocols.humanconnectome.org/HCP/3T/imaging-protocols.html) sequences. The tests were run on a laptop with a four core (8 thread, 28w) MacOS laptop:

| Command : Seconds (GZ)                                 |  Serial (GZ)  | Parallel (GZ) |
|--------------------------------------------------------|--------------:|--------------:|
| fslmaths rest -s 2.548 out : 270 (424)                 | 5.0x (2.9x)   | 8.6x (6.3x)   |
| fslmaths t1 -kernel boxv 7 -dilM out : 216 (228)       | 245x (41x)    | 225x (72x)    |
| fslmaths rest -Tmean -mul -1 -add rest out : 101 (328) | 2.5x (2.5x)   | 2.8x (4.5x)   |
|  niimath rest -demean out (same output as above)       | 3.5x (3.0x)   | 4.6x (6.2x)   |
| fslmaths rest -bptf 77 8.68 out : 998 (1155)           | 2.0x (2.0x)   | 6.8x (6.7x)   |

Here are the same testson a desktop computer with twelve cores (24 threads, Ryzen 3900X):

| Command : Seconds (GZ)                                 |  Serial (GZ)  | Parallel (GZ) |
|--------------------------------------------------------|--------------:|--------------:|
| fslmaths rest -s 2.548 out : 123 (229)                 | 4.2x (2.4x)   | 9.9x (12.1x)  |
| fslmaths t1 -kernel boxv 7 -dilM out : 156 (159)       | 371x (37x)    | 371x (248x)   |
| fslmaths rest -Tmean -mul -1 -add rest out : 32 (186)  | 1.7x (2.5x)   | 1.8x (7.6x)   |
|  niimath rest -demean out (same output as above)       | 2.6x (2.6x)   | 3.0x (10.8x)  |
| fslmaths rest -bptf 77 8.68 out : 887 (1019)           | 2.6x (2.5x)   | 23x (23.0x)   |

Gaussian smoothing (`-s`/`-dog`/`unsharp`) uses a contiguous vectorizable kernel in every build (native, WASM, and the shared registration pyramid). Neighborhood mean, minimum, maximum, and erosion filters keep their local gathers but evaluate adjacent interior outputs in SIMD lanes, avoiding the non-finite propagation errors of separable running-sum and deque filters.

## Converting voxelwise images to a triangulated mesh

niimath can convert NIfTI images to meshes, suitable for viewing in Surfice, blender, SUMA, FreeSurfer and other tools. The features are based on [nii2mesh](https://github.com/neurolabusc/nii2mesh) and the features are almost identical. However, the order of arguments is different to match the expectations of fslmaths/niimath. So the call `nii2mesh -r 1 bet.nii.gz r100.ply` becomes `niimath bet.nii.gz -mesh -r 1 r100.ply`. The benefit of niimath is that you can apply voxel-based operations before you create your mesh. This allows you to apply morphological operations (`-close`, `-ero`, `-dilM`). As an example, to apply a 4mm Gaussian smooth before creating a mesh, you could run `./niimath mni152.nii.gz -s 4 -mesh -i 122 -l 0 -b 1 b1.ply`. As described on the [nii2mesh](https://github.com/neurolabusc/nii2mesh) page, you can create independent meshes for each area in an atlas using the command:

```
niimath D99_atlas_v2.0_right.nii.gz -mesh -p 0 -s 10 -a D99_v2.0_labels_semicolon.txt ./gii/D99s10roi.gii
```
Both programs allow you to explicitly set the isolevel using the `-i` value, so `-i 128` we render a surface for voxels brighter than 128. One minor difference between the programs is that niimath allows you also request `dark`, `medium` and `bright` using the `-i d`, `-i m` and `-i b` commands respectively. These use Otsu's method, and typically identify pleasing values. Also, if the user does not specify an isolevel be aware that nii2mesh chooses the middle brightness (the midpoint between the darkest and brightest value) while niimath uses the medium Otsu threshold. The latter is more robust to outliers. Here are examples illustrating this usage:

```
niimath bet.nii.gz -mesh -i 128 Isolevel128.gii
niimath bet.nii.gz -mesh -i d darkIsolevel.gii
niimath bet.nii.gz -mesh -i m medIsolevel.gii
niimath bet.nii.gz -mesh -i b brightIsolevel.gii
```

Mesh quality. Marching cubes uses the Lewiner tables by default, which resolve the ambiguous cube configurations against the trilinear interpolant (`-o 1` selects the classic tables). Simplification never breaks the topology of the surface (a link condition rejects any collapse that would create a non-manifold edge). `-q` sets the quality level: `-q 2` (the default) adds a self-intersection guard to both the `-s` smoothing and the simplification and a lossless finishing pass; `-q 1` drops the guards and the lossless finish and runs in about 40% of the time (2.3 s against 5.7 s for a smoothed, simplified MNI surface); `-q 0` is fastest and writes uncompressed mz3. `-r 1` leaves the mesh unsimplified at every quality level. An existing mesh can be processed the same way: `niimath in.mz3 -s 10 -r 0.5 out.mz3` accepts `-r`, `-s`, `-q` and `-v`. With `-v 1` every stage prints a `mesh check` line: components, Euler characteristic and genus, boundary edges and holes, non-manifold edges and vertices, and self-intersecting triangles. A second, half-edge simplifier (`-n 1`) exists as a reference: it stops within one face of the requested count and refuses non-manifold input, at the same geometric fidelity and about 30% more time. It is not compiled by default (`Q2=1 make` or `-DENABLE_QUADRIC2=ON`).

## Creating bitmaps

You can use the `--bitmap` option to visualize the results of any operations. This option has arguments inspired by fsl's slicer, but introduces new features. The [niimath-bitmap](https://github.com/rordenlab/niimath-bitmap) provides examples and documentation.

## WebAssembly

niimath can also be compiled to WebAssembly (Wasm) allowing it to be inserted into web pages and Node.js projects. Here is a [live demo](https://niivue.github.io/niivue-niimath/) with links to source code and instructions.

## License

<!-- codespell-ignore-line --> niimath is licensed under the 2-Clause BSD License. Except where noted, the code was written by Chris Rorden in 2020-2022. The code in `tensor.c` was written by Daniel Glen (2004) from the US National Institutes of Health and is not copyrighted (though it is included here with the permission of the author). The FSL team graciously allowed the text strings (help, warning and error messages) to be copied verbatim. Taylor Hanayik from the FSL group provided pseudo-code for some functions where there is little available documentation. The PolygoniseCube function comes from Cory Bloyd's public domain [Marching Cubes example](http://paulbourke.net/geometry/polygonise/) program described here. The bwlabel.cpp file was written by Jesper Andersson, who has explicitly allowed this to be shared using the BSD 2-Clause license. The [high performance](https://github.com/gaspardpetit/base64) base64.cpp was written by Jouni Malinen and is distributed under the BSD license. The mesh simplification was written by [Sven Forstmann](https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification) and distributed under the MIT license. It was ported from C++ to C by Chris Rorden.  The [radixsort.c](https://github.com/bitshifter/radixsort) was written by Cameron Hart (2014) using the zlib license.

The `-romeo` phase-unwrapping command (`src/romeo.c`) is a C port of [ROMEO.jl](https://github.com/korbinian90/ROMEO.jl) and the [MriResearchTools.jl](https://github.com/korbinian90/MriResearchTools.jl) helpers its command-line app uses, by Korbinian Eckstein, Barbara Dymerska and Simon Robinson, together with the 2π range reduction from the Julia standard library. All are distributed under the MIT license; the upstream copyright and permission notices are preserved verbatim in `src/romeo.LICENSE`. Unlike the GPL module below, `-romeo` is compiled in **by default**, so a standard niimath binary contains this MIT-licensed component — MIT is compatible with the 2-Clause BSD License, so the binary as a whole remains BSD-2-Clause. `ROMEO=0 make` (or `-DENABLE_ROMEO=OFF`) omits it.

The `--medic` and `-unwarp` commands (`src/medic.c`) are original BSD-2-Clause code by the niimath authors: a clean-room emulation of the MEDIC method published by Van et al. (*Imaging Neuroscience* 4, 2026, [doi:10.1162/IMAG.a.1262](https://doi.org/10.1162/IMAG.a.1262)), developed from the paper and from black-box measurement of the reference tool's public executables. No reference implementation, test, build product or debug symbol was read, and no code from it is included; the measurements that fix each convention are recorded in the `medic_bench` repository. Phase unwrapping is performed by the MIT-licensed `-romeo` port described above, so `--medic` requires it (`ROMEO=0` implies `MEDIC=0`); `MEDIC=0 make` or `-DENABLE_MEDIC=OFF` omits MEDIC alone.

The `-moco` and `-stc` commands (`src/moco.c`, `src/stc.c`) are original BSD-2-Clause code by the niimath authors. Both emulate a published AFNI method whose reference implementation is copyrighted by the Medical College of Wisconsin — `3dvolreg`, `mri_3dalign`, `thd_rot3d` and `thd_shear3d` for `-moco`; `3dTshift` and its FFT for `-stc`. Those files were **GPL-2** when this code was written; on 12 May 2026 MCW relicensed its 1994-2000 AFNI code to **CC BY 4.0**, which removes the copyleft bar but introduces attribution and change-notice duties. Those sources were **not** read, translated or paraphrased; they served only as black-box oracles. The clean-room specification is the published method (Cox & Jesmanowicz 1999 for `-moco`; AFNI's published `3dTshift -help` and `-verbose` output for `-stc`) together with measured inputs and outputs, recorded in the `moco_bench` repository (`test/moco_reference_manifest.md` and `test/stc_reference_manifest.md`). The FFT in `stc.c` is original niimath code — a batched Stockham autosort kernel; no FFT implementation was read or adapted. Neither command carries any attribution obligation as a result. A binary containing these commands remains BSD-2-Clause.


The optional `-skullstrip` command (`src/skullstrip.c`) adapts public-domain AFNI code by Robert W. Cox and colleagues (NIMH): spatial normalisation from `thd_brainormalize.c` and `thd_automask.c`, and the surface deformation and touchup stages from `SUMA_BrainWrap.c`, which carries no copyright notice and so falls under AFNI's US-Government-work clause — a US Government work is not copyrightable (17 U.S.C. §105). AFNI states this affirmatively rather than by implication: its `LICENSE.txt` declares the tree a "United States Government Work" apart from a listed set of exceptions and states that "contributions without explicit licensing will be assumed to be entered into the public domain", and its `README.copyright` dates the rule to work after 15 Jan 2001 (the adapted files' first commits are 2001-2004, by NIH authors). AFNI's `SUMA_3dedge3`, which wraps Malandain's GPL-3.0 `Extract_Gradient_Maxima_3D`, is deliberately out of scope, which is why niimath implements only the `-no_use_edge` behaviour. **One item was resolved rather than argued.** The nearest-neighbour index conversion was originally written after reading `THD_3dmm_to_3dind_warn` in AFNI's `thd_coords.c`, which carries an MCW copyright header and predates the public-domain cutoff (GPL-2 at the time; CC BY 4.0 since 12 May 2026). It has been **replaced by a clean-room reimplementation** (`ss_world_to_index`): an implementer who had read neither AFNI nor niimath reproduced a 19,139-row table of measured input/output behaviour, working only from that table, and the resulting masks are byte-identical. The protocol, the table, the sufficiency check and the implementer's derivation account are in the `skullstrip_bench` repository's `clean_room/` directory. The relicense has not changed this: a clean-room result carries no attribution duty, where adapting the CC BY original would. The surface primitives (icosphere, adjacency, intersection testing, rasterisation) are original BSD-2-Clause code. A binary containing this command remains BSD-2-Clause. `-skullstrip` is **off by default** (`SKULLSTRIP=1 make`, `cmake -DENABLE_SKULLSTRIP=ON`) and is absent from every released binary.

The optional `-spm_coreg` and `-spm_deface` commands are the project's entire copyleft payload, carried in the [niimath_gpl](https://github.com/rordenlab/niimath_gpl) submodule at `src/GPL` and enabled only when built with `make GPL=1` / `cmake -DENABLE_GPL=ON` (`-DHAVE_GPL`). They link SPM's `spm_coreg` module, which is **GPL-2 or later**, so distribute a `GPL=1` binary under those terms. The default build contains none of this code and remains BSD-2-Clause. The version string reported by `niimath` ends in ` GPL` or ` BSD` to indicate which applies, and `release_smoke.py` asserts both directions — a BSD binary that still carried `-spm_coreg` would fail the release gate.

## Links

  - [imbibe](https://github.com/jonclayden/imbibe) is a R wrapper for niimath, allowing the performance of tuned code with the convenience of a scripting language.
  - [3dcalc](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html) is AFNI's tool for image arithmetic.
  - [c3d](https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md) provides mathematical functions and format conversion for medical images.
  - [fslmaths](https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro3/index.html) is the inspiration for niimath.

## Citation

  - Rorden C, Webster M, Drake C,  Jenkinson M, Clayden JD, Li N, Hanayik T ([2024](https://apertureneuro.org/article/94384-niimath-and-fslmaths-replication-as-a-method-to-enhance-popular-neuroimaging-tools)) niimath and fslmaths: replication as a method to enhance popular neuroimaging tools. Aperture Neuro.4. doi:10.52294/001c.94384

# niimath

[![Build status](https://ci.appveyor.com/api/projects/status/7o0xp2fgbhadkgn1?svg=true)](https://ci.appveyor.com/project/neurolabusc/niimath)

## About

It is said that `imitation is the sincerest form of flattery`. This project emulates the popular [fslmaths](https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro3/index.html) tool. fslmaths is a `general image calculator` and is not only one of the foundational tools for FSL's brain imaging pipelines (such as [FEAT](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT)), but has also been widely adopted by many tools. This popularity suggests that it fulfills an important niche. While scientists are often encouraged to discover novel solutions, it sometimes seems that replication is undervalued. Here are some specific reasons for creating this tool:

1. While fslmaths is provided for without charge, it is not [open source](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence). This limits its inclusion in other projects, in particular for commercial exploitation.
2. Using an open source license allows niimath to build with open source libraries that the FSL team can not use. Specifically, the CloudFlare zlib provides dramatically faster performance than the public domain library used by fslmaths. n.b. Subsequently, we helped update [CloudFlare zlib](https://github.com/cloudflare/zlib/pull/19) that allows recent FSL releases to use this library,  improving the speed for all FSL tools.
3. Minimal dependencies allow easy distribution, compilation and development. For example, it can be compiled for MacOS, Linux and Windows (fsl can not target Windows).
4. Designed from ground up to optionally use parallel processing (OpenMP and CloudFlare-enhanced [pigz](https://github.com/madler/pigz)).
5. Most programs are developed organically, with new features added as need arises. Cloning an existing tool provides a full specification, which can lead to optimization. niimath uses explicit single and double precision pipelines that allow the compiler to better use advanced instructions (every x86_64 CPU provides SSE, but high level code has trouble optimizing these routines). The result is that modern compilers are able to create operations that are limited by memory bandwidth, obviating the need for [hand tuning](https://github.com/neurolabusc/simd) the code.
6. Developing a robust regression testing dataset has allowed us to discover a few edge cases where fslmaths provides anomalous or unexpected answers (see below). Therefore, this can benefit the popular tool that is being cloned.
7. While the code is completely reverse engineered, the FSL team has been gracious to allow us to copy their error messages and help information. This allows true plug in compatibility. They have also provided pseudo code for poorly documented routines. This will allow the community to better understand the actual algorithms.
8. This project provides an open-source foundation to introduce new features that fill gaps with the current FSL tools (e.g. unsharp, sobel, resize functions). For future releases, Bob Cox has graciously provided permission to use code from [AFNI's](https://afni.nimh.nih.gov) 3dTshift and 3dBandpass tools that provide performance unavailable within [FSL](https://neurostars.org/t/bandpass-filtering-different-outputs-from-fsl-and-nipype-custom-function/824). Including them in this project ensures they work in a familiar manner to other FSL tools (and leverage the same environment variables).

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

The easiest way to build niimath on a Unix computer is to use cmake:

```
git clone https://github.com/rordenlab/niimath.git
cd niimath; mkdir build; cd build; cmake ..
make
```

If you want to enable OpenMP support on macOS, you have to install `libomp` first using `brew install libomp`, and then use `cmake -DOPENMP_XCODE=ON ..` to configure the project in the above commands.

Likewise, if you are compiling on Windows using cmake:

```
git clone https://github.com/rordenlab/niimath.git
cd niimath & mkdir build & cd build & cmake ..
cmake --build .
```
Alternatively, you can compile the software by running the terminal command `make` from the project's `src` folder if you are running Linux (or execute `windows.bat` if you are running Windows):

```
git clone https://github.com/rordenlab/niimath.git
cd niimath/src
make
```

You can also compile this project to Web Assembly so it can be embedded in a web page, as shown in the [live demo](https://niivue.github.io/niivue-niimath/).

```
git clone https://github.com/rordenlab/niimath.git
cd niimath/src
make wasm
```

Advanced users using the `Makefile` may want to run `CF=1 OMP=1 make -j` to make a version that uses OpenMP (parallel processing) and the CloudFlare accelerated compression library. You may need to edit the `Makefile` for your compiler name. On MacOS, the default C compiler is Clang, which has [poor OpenMP](https://github.com/neurolabusc/simd) support. Therefore, MacOS users may want to install the gcc compiler (for example, `brew install gcc@9`).

For Windows, using the cmake method described above is highly recommended. However, you can also compile the project directly from the command line (here without the `-DHAVE_ZLIB` directive, so gz files will not be supported) :

```
cl /Feniimath niimath.c core.c tensor.c bwlabel.c bw.c core32.c core64.c fdr.c meshify.c MarchingCubes.c quadric.c base64.c radixsort.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib -DNII2MESH
```

Simply running `make` in the `src` folder should compile niimath on Linux. This should work regardless of if you use the Clang/LLVM or gcc compiler. However, the resulting executable will only work with specific versions of Linux. If you want to make a universal Linux release you can use [holy-build-box](https://github.com/FooBarWidget/holy-build-box). Be aware that this uses an old version of the gcc compiler (4.8.5), so the resulting performance may not be optimized for your system.

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

niimath provides the same commands as [fslmaths](https://mandymejia.com/fsl-maths-commands/), so you can use it just as you would fslmaths. If you are brave, you can even rename it fslmaths and use it as a drop in replacement. You can also modify your environment variables to unleash advanced features:

 - Just like fslmaths, it uses your [`FSLOUTPUTTYPE` Environment Variable ](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslEnvironmentVariables) to determine output file format. Unix users can specify `export NIFTI_GZ` or `export NIFTI` from the command line or profile to select between compressed (smaller) or uncompressed (faster) results. Windows users can use `set` instead of `export`.
 - To turn on parallel processing and threading, you can either set the environment variable `export AFNI_COMPRESSOR=PIGZ`. If the environment variable `AFNI_COMPRESSOR` does not exist, or is set to any value other than `PIGZ` you will get single threaded compresson.

niimath has a few features not provided by fslmaths:

 - `bandpass <hp> <lp> <tr>`: Butterworth filter, highpass and lowpass in Hz,TR in seconds (zero-phase 2*2nd order filtfilt)
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
 - `resize <X> <Y> <Z> <m>` : grow (>1) or shrink (<1) image. Method <m> (0=nearest,1=linear,2=spline,3=Lanczos,4=Mitchell)\n");
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
 - `--compare <ref>`       : report if images are identical, terminates without saving new image\n");
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

## WebAssembly

niimath can also be compiled to WebAssembly (Wasm) allowing it to be inserted into web pages and Node.js projects. Here is a [live demo](https://niivue.github.io/niivue-niimath/) with links to source code and instructions.

## License

<!-- codespell-ignore-line --> niimath is licensed under the 2-Clause BSD License. Except where noted, the code was written by Chris Rorden in 2020-2022. The code in `tensor.c` was written by Daniel Glen (2004) from the US National Institutes of Health and is not copyrighted (though it is included here with the permission of the author). The FSL team graciously allowed the text strings (help, warning and error messages) to be copied verbatim. The Butterworth Filter Coefficients in `bw.c` are from [Exstrom Labs](http://www.exstrom.com/journal/sigproc/) and the authors provided permission for it to be included in this project under the [LGPL](https://www.gnu.org/licenses/lgpl-3.0.en.html), the file provides additional details. Taylor Hanayik from the FSL group provided pseudo-code for some functions where there is little available documentation. The PolygoniseCube function comes from Cory Bloyd's public domain [Marching Cubes example](http://paulbourke.net/geometry/polygonise/) program described here. The bwlabel.cpp file was written by Jesper Andersson, who has explicitly allowed this to be shared using the BSD 2-Clause license. The [high performance](https://github.com/gaspardpetit/base64) base64.cpp was written by Jouni Malinen and is distributed under the BSD license. The mesh simplification was written by [Sven Forstmann](https://github.com/sp4cerat/Fast-Quadric-Mesh-Simplification) and distributed under the MIT license. It was ported from C++ to C by Chris Rorden.  The [radixsort.c](https://github.com/bitshifter/radixsort) was written by Cameron Hart (2014) using the zlib license.

## Links

  - [imbibe](https://github.com/jonclayden/imbibe) is a R wrapper for niimath, allowing the performance of tuned code with the convenience of a scripting language.
  - [3dcalc](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html) is AFNI's tool for image arithmetic.
  - [c3d](https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md) provides mathematical functions and format conversion for medical images.
  - [fslmaths](https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro3/index.html) is the inspiration for niimath.

## Citation

  - Rorden C, Webster M, Drake C,  Jenkinson M, Clayden JD, Li N, Hanayik T ([2024](https://apertureneuro.org/article/94384-niimath-and-fslmaths-replication-as-a-method-to-enhance-popular-neuroimaging-tools)) niimath and fslmaths: replication as a method to enhance popular neuroimaging tools. Aperture Neuro.4. doi:10.52294/001c.94384

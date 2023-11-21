### Introduction

niimath can be compiled to [WebAssembly (WASM)](https://webassembly.org) so that it can be inserted into JavaScript projects. For example, it is embedded into [NiiVue](https://github.com/niivue/niivue) web pages. It can also be run from the command line using [Node.js](https://nodejs.org/en/). Node.js provides a simple way to test the WASM. However, be aware that the WASM code does not have the full functionality of the natively compiled C code.

### Compiling a minimal Node.js test

The niimath.js project provides a similar command line operation to fslmaths. It uses [NIFTI-Reader-JS](https://github.com/rii-mango/NIFTI-Reader-JS) to load NIfTI images from disk.

 1. Compile C to WASM: `emcc -O2 -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4GB -s WASM=1 -DUSING_WASM -I. core32.c nifti2_wasm.c core.c walloc.c -o funcx.js`
 3. Run the JavaScript example: `node test.js`

**The development branch requires additional files:**
**emcc -O2 -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4GB -s WASM=1 -DUSING_WASM -DNII2MESH -I. core32.c nifti2_wasm.c core.c meshify.c MarchingCubes.c  walloc.c base64.c bwlabel.c quadric.c radixsort.c -o funcx.js**

This project does not load or save images from disk. It simply generates fake images and processes them. It reports the time for each operation:

```
dims: 128x128x128x1 nbyper: 4
input  0 , 1 , 2 , 3 , 0 , 1 ...
   142ms -> ' -dog 2 3 ' -> output  0 , 1 , 0 , 0 , 1 , 1 ...

dims: 128x128x128x1 nbyper: 4
input  0 , 1 , 2 , 3 , 0 , 1 ...
   174ms -> ' -dehaze 3 -dog 2 3 ' -> output  0 , 1 , 0 , 0 , 1 , 1 ...
...
```

### Compiling a fslmaths/niimath Node.js clone

The niimath.js project provides a similar command line operation to fslmaths. It uses [NIFTI-Reader-JS](https://github.com/rii-mango/NIFTI-Reader-JS) to load NIfTI images from disk.

 1. Install dependencies: `npm install nifti-reader-js`
 2. Compile C to WASM: `emcc -O2 -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4GB -s TOTAL_MEMORY=268435456 -s WASM=1 -DUSING_WASM -I. core32.c nifti2_wasm.c core.c walloc.c -o funcx.js`
 3. Run the JavaScript command using operations similar to fslmaths/niimath: `node niimath.js T1.nii -sqr sT1.nii`


### Limitations

The table below shows the supported functions. In general, functions that change the shape of the image (e.g. cropping or shrinking) are not supported and listed as `NO` in the JavaScript (`JS`) column. 

This WASM unit only works on a single input image, so the call `node niimath.js in.nii -add 2 out.nii` will work, while the operation `node niimath.js in.nii -add in2.nii out.nii` will fail.

The WASM code will never provide any test output. The function returns zero on success and some other value on failure. In contrast, the native niimath executable will often generate error messages helping the user determine why the process failed.

WASM memory management is limited. JavaScript protects memory, improving security at the expense of performance. In contrast, C allows users direct access to memory. C code compiled to WASM must bridge this gap. This project uses [walloc](https://github.com/wingo/walloc) to allocate memory. JavaScript does not appreciate memory reallocation, which re-references pointers. However, many niimath commands must temporarily allocate memory for processing. The consequence is that one must compile niimath to WASM with a large TOTAL_MEMORY value to avoid memory reallocation. In other words, you need to pre-allocate the worst case scenario memory usage, rather than growing and shrinking memory as usage demands. This does mean that this code will require a considerable amount of memory, even if you are only viewing small images.

| Function                   | JS  |                                                                   |
|----------------------------|---- |--------------------------------------------------------------------|
| -bandpass [hp] [lp] [tr]  | NO  |  Butterworth filter, highpass and lowpass in Hz,TR in seconds |
| -bptfm [hp] [lp]          | YES |  Same as bptf but does not remove mean (emulates fslmaths < 5.0.7) |
| -bwlabel [conn]           | NO  |  Connected component labelling for non-zero voxels (conn sets neighbors:  6, 18, 26)  |
| -c2h                      | YES |  reverse h2c transform |
| -ceil                     | YES |  round voxels upwards to the nearest integer |
| -crop [tmin] [tsize]      | NO  |  remove volumes, starts with 0 not 1! Inputting -1 for a size will set it to the full range |
| -dehaze [mode]            | YES |  set dark voxels to zero (mode 1..5; higher yields more surviving voxels) |
| -detrend                  | YES |  remove linear trend (and mean) from input |
| -demean                   | YES |  remove average signal across volumes (requires 4D input) |
| -edt                      | YES |  estimate Euler Distance Transform (distance field). Assumes isotropic input |
| -floor                    | YES |  round voxels downwards to the nearest integer |
| -mod                      | YES |  modulus fractional remainder - same as '-rem' but includes fractions |
| -otsu [mode]              | YES |  binarize image using Otsu's method (mode 1..5; higher yields more bright voxels)) |
| -power [exponent]         | YES |  raise the current image by following exponent |
| -h2c                      | YES |  convert CT scans from 'Hounsfield' to 'Cormack' units to emphasize soft tissue contrast |
| -p [threads]              | YES |  set maximum number of parallel threads. DISABLED :  recompile for OpenMP support |
| -resize [X] [Y] [Z] [m]   | YES |  grow (>1) or shrink (<1) image. Method [0=nearest,1=linear,2=spline,3=Lanczos,4=Mitchell] |
| -round                    | YES |  round voxels to the nearest integer |
| -sobel                    | YES |  fast edge detection |
| -sobel_binary             | YES |  sobel creating binary edge |
| -tensor_2lower            | NO  |  convert FSL style upper triangle image to NIfTI standard lower triangle order |
| -tensor_2upper            | NO  |  convert NIfTI standard lower triangle image to FSL style upper triangle order |
| -tensor_decomp_lower      | NO  |  as tensor_decomp except input stores lower diagonal (AFNI, ANTS, Camino convention) |
| -trunc                    | YES |  truncates the decimal value from floating point value and returns integer value |
| -unsharp  [sigma] [scl]   | YES |  edge enhancing unsharp mask (sigma in mm, not voxels; 1.0 is typical for amount (scl)) |
| -dog [sPos] [sNeg]        | YES |  difference of gaussian with zero-crossing edges (positive and negative sigma mm) |
| -dogr [sPos] [sNeg]       | YES |  as dog, without zero-crossing (raw rather than binarized data) |
| -dogx [sPos] [sNeg]       | YES |  as dog, zero-crossing for 2D sagittal slices |
| -dogy [sPos] [sNeg]       | YES |  as dog, zero-crossing for 2D coronal slices |
| -dogz [sPos] [sNeg]       | YES |  as dog, zero-crossing for 2D axial slices |
| --compare [ref]           | YES |  report if images are identical, terminates without saving new image |
| filename.nii              | YES |  mimic fslhd (can also export to a txt file : 'niimath T1.nii 2] T1.txt') report header and terminate without saving new image |
| -add    | YES |  add following input to current image |
| -sub    | YES |  subtract following input from current image |
| -mul    | YES |  multiply current image by following input |
| -div    | YES |  divide current image by following input |
| -rem    | YES |  modulus remainder - divide current image by following input and take remainder |
| -mas    | YES |  use (following image >0) to mask current image |
| -thr    | YES |  use following number to threshold current image (zero anything below the number) |
| -thrp   | YES |  use following percentage (0-100) of ROBUST RANGE to threshold current image (zero anything below the number) |
| -thrP   | YES |  use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold below |
| -uthr   | YES |  use following number to upper-threshold current image (zero anything above the number) |
| -uthrp  | YES |  use following percentage (0-100) of ROBUST RANGE to upper-threshold current image (zero anything above the number) |
| -uthrP  | YES |  use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold above |
| -clamp  | YES |  use following percentage (0-100) of ROBUST RANGE to threshold current image (anything below set to this threshold) |
| -uclamp | YES |  use following percentage (0-100) of ROBUST RANGE to threshold current image (anything above set to this threshold) |
| -max    | YES |  take maximum of following input and current image |
| -min    | YES |  take minimum of following input and current image |
| -seed   | YES |  seed random number generator with following number |
| -restart  | NO  |  replace the current image with input for future processing operations |
| -save  | NO  |  save the current working image to the input filename |
| -inm [mean]  | YES |   (-i i ip.c) intensity normalisation (per 3D volume mean) |
| -ing [mean]  | YES |   (-I i ip.c) intensity normalisation, global 4D mean) |
| -s [sigma]  | YES |  create a gauss kernel of sigma mm and perform mean filtering |
| -exp    | YES |  exponential |
| -log    | YES |  natural logarithm |
| -sin    | YES |  sine function |
| -cos    | YES |  cosine function |
| -tan    | YES |  tangent function |
| -asin   | YES |  arc sine function |
| -acos   | YES |  arc cosine function |
| -atan   | YES |  arc tangent function |
| -sqr    | YES |  square |
| -sqrt   | YES |  square root |
| -recip  | YES |  reciprocal (1/current image) |
| -abs    | YES |  absolute value |
| -bin    | YES |  use (current image >0) to binarise |
| -binv   | YES |  binarise and invert (binarisation and logical inversion) |
| -fillh  | YES |  fill holes in a binary mask (holes are internal - i.e. do not touch the edge of the FOV) |
| -fillh26  | YES |  fill holes using 26 connectivity |
| -index  | YES |  replace each nonzero voxel with a unique (subject to wrapping) index number |
| -grid [value] [spacing]  | YES |  add a 3D grid of intensity [value] with grid spacing [spacing] |
| -edge   | YES |  edge strength |
| -tfce [H] [E] [connectivity] | NO  |  enhance with TFCE, e.g. -tfce 2 0.5 6 (maybe change 6 to 26 for skeletons) |
| -tfceS [H] [E] [connectivity] [X] [Y] [Z] [tfce_thresh] | NO  |  show support area for voxel (X,Y,Z) |
| -nan    | YES |  replace NaNs (improper numbers) with 0 |
| -nanm   | YES |  make NaN (improper number) mask with 1 for NaN voxels, 0 otherwise |
| -rand   | YES |  add uniform noise (range 0 | YES | 1) |
| -randn  | YES |  add Gaussian noise (mean=0 sigma=1) |
| -range  | YES |  set the output calmin/max to full data range |
| -tensor_decomp  | NO  |  convert a 4D (6-timepoint )tensor image into L1,2,3,FA,MD,MO,V1,2,3 (remaining image in pipeline is FA) |
| -kernel 3D  | YES |  3x3x3 box centered on target voxel (set as default kernel) |
| -kernel 2D  | YES |  3x3x1 box centered on target voxel |
| -kernel box    [size]      | YES |  all voxels in a cube of width [size] mm centered on target voxel |
| -kernel boxv   [size]      | YES |  all voxels in a cube of width [size] voxels centered on target voxel, CAUTION | YES |  size should be an odd number |
| -kernel boxv3  [X] [Y] [Z] | YES |  all voxels in a cuboid of dimensions X x Y x Z centered on target voxel, CAUTION | YES |  size should be an odd number |
| -kernel gauss  [sigma]     | YES |  gaussian kernel (sigma in mm, not voxels) |
| -kernel sphere [size]      | YES |  all voxels in a sphere of radius [size] mm centered on target voxel |
| -kernel file   [filename]  | NO  |  use external file as kernel |
| -dilM     | YES |  Mean Dilation of non-zero voxels |
| -dilD     | YES |  Maximum Dilation of non-zero voxels (emulating output of fslmaths 6.0.1, max not modal) |
| -dilF     | YES |  Maximum filtering of all voxels |
| -dilall   | YES |  Apply -dilM repeatedly until the entire FOV is covered |
| -ero      | YES |  Erode by zeroing non-zero voxels when zero voxels found in kernel |
| -eroF     | YES |  Minimum filtering of all voxels |
| -fmedian  | YES |  Median Filtering  |
| -fmean    | YES |  Mean filtering, kernel weighted (conventionally used with gauss kernel) |
| -fmeanu   | YES |  Mean filtering, kernel weighted, un-normalized (gives edge effects) |
| -s [sigma]  | YES |  create a gauss kernel of sigma mm and perform mean filtering |
| -subsamp2   | NO  |  downsamples image by a factor of 2 (keeping new voxels centered on old) |
| -subsamp2offc   | NO  |  downsamples image by a factor of 2 (non-centered) |
| -Tmean    | NO  |  mean across time |
| -Tstd     | NO  |  standard deviation across time |
| -Tmax     | NO  |  max across time |
| -Tmaxn    | NO  |  time index of max across time |
| -Tmin     | NO  |  min across time |
| -Tmedian  | NO  |  median across time |
| -Tperc [percentage]  | NO  |  nth percentile (0-100) of FULL RANGE across time |
| -Tar1     | NO  |  temporal AR(1) coefficient (use -odt float and probably demean first) |
| -pval     | YES |  Nonparametric uncorrected P-value, assuming timepoints are the permutations; first timepoint is actual (unpermuted) stats image |
| -pval0    | YES |  Same as -pval, but treat zeros as missing data |
| -cpval    | YES |  Same as -pval, but gives FWE corrected P-values |
| -ztop     | YES |  Convert Z-stat to (uncorrected) P |
| -ptoz     | YES |  Convert (uncorrected) P to Z |
| -rank     | NO  |  Convert data to ranks (over T dim) |
| -ranknorm | NO  |  Transform to Normal dist via ranks |
| -roi [xmin] [xsize] [ymin] [ysize] [zmin] [zsize] [tmin] [tsize]  | NO  |  zero outside roi (using voxel coordinates). Inputting -1 for a size will set it to the full image extent for that dimension |
| -bptf  [hp_sigma] [lp_sigma]  | NO  |  (-t in ip.c) Bandpass temporal filtering; nonlinear highpass and Gaussian linear lowpass (with sigmas in volumes, not seconds); set either sigma[0 to skip that filter |
| -roc [AROC-thresh] [outfile] [4Dnoiseonly] [truth]  | NO  |  take (normally binary) truth and test c |

### Alternative

One can use emscriptens built in file system and memory support if you compile as follows:

```
emcc -s USE_ZLIB=1 -O1 -DHAVE_ZLIB -DFSLSTYLE -DREJECT_COMPLEX -DNII2MESH -std=gnu99 -o niimathwasm.js niimath.c MarchingCubes.c meshify.c quadric.c base64.c radixsort.c fdr.c bwlabel.c bw.c core.c tensor.c core32.c core64.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib -lm -lz
```
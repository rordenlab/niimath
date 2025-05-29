/* ----------------------------------------------------------------------
 * A clone of fslmaths by Chris Rorden (2020)
 *
 * compile example (consider -pedantic or -Wall):
 *
 * gcc -O3 -lm -lz -o niimath niimath.c  niftilib/nifti1_io.c znzlib/znzlib.c -I./niftilib
 *
 * OpenMP (parallel threading)
 *  gcc-9  -fopenmp -lm  -I./darwin ./darwin/libz.a -DHAVE_ZLIB -o niimath niimath.c  niftilib/nifti1_io.c znzlib/znzlib.c -I./niftilib -I./znzlib
 *
 *----------------------------------------------------------------------
 */

#ifdef _WIN32
	#include <fcntl.h>
	#include <io.h>
#endif
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <nifti2_io.h>
#include "core32.h" //all 32-bit functions
#ifdef HAVE_64BITS
	#include "core64.h" //all 64-bit functions
#endif
#ifdef NII2MESH
	#include <stdbool.h>
	#include "meshtypes.h"
	#include "meshify.h"
	#include "quadric.h"
	#ifdef _MSC_VER
		#include <io.h>
		#define access _access
		#define F_OK    00       /* Test for existence.  */
	#else
		#include <unistd.h>
	#endif
	#include <stdio.h>
#endif

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#if defined(_OPENMP)
   #define kOMPsuf " OpenMP"
#else
   #define kOMPsuf ""
#endif
#if defined(__ICC) || defined(__INTEL_COMPILER)
	#define kCCsuf  " IntelCC" STR(__INTEL_COMPILER)
#elif defined(_MSC_VER)
	#define kCCsuf  " MSC" STR(_MSC_VER)
#elif defined(__clang__)
	#define kCCsuf  " Clang" STR(__clang_major__) "." STR(__clang_minor__) "." STR(__clang_patchlevel__)
#elif defined(__GNUC__) || defined(__GNUG__)
    #define kCCsuf  " GCC" STR(__GNUC__) "." STR(__GNUC_MINOR__) "." STR(__GNUC_PATCHLEVEL__)
#else
	#define kCCsuf " CompilerNA" //unknown compiler!
#endif
#if defined(__APPLE__)
	#define kOS "MacOS"
#elif (defined(__linux) || defined(__linux__))
	#define kOS "Linux"
#else
	#define kOS "Windows"
#endif

#define kMTHdate "v1.0.20250529"
#define kMTHvers kMTHdate kOMPsuf kCCsuf

#ifdef NII2MESH

bool isMz3(const char *fnm) {
	char basenm[768], ext[768] = "";
	strcpy(basenm, fnm);
	strip_ext(basenm); // ~/file.nii -> ~/file
	if (strlen(fnm) > strlen(basenm))
		strcpy(ext, fnm + strlen(basenm));
	return strstr(ext, ".mz3");
}

void read_mz3(const char* filename, vec3d **verts, vec3i **tris, int* nvert, int* ntri) {
	#ifdef _MSC_VER
		#pragma pack(2)
			struct mz3hdr {
				uint16_t SIGNATURE, ATTR;
				uint32_t NFACE, NVERT, NSKIP;
			};
		#pragma pack()
	#else
		struct __attribute__((__packed__)) mz3hdr {
			uint16_t SIGNATURE, ATTR;
			uint32_t NFACE, NVERT, NSKIP;
		};
	#endif
	if( access( filename, F_OK ) != 0 ) {
		printf("Unable to find a mz3 named %s\n", filename);
		exit(EXIT_FAILURE);
	}
	FILE *fp = fopen(filename,"rb");
	struct mz3hdr h;
	size_t bytes_read = fread(&h, sizeof(struct mz3hdr), 1, fp);
	if (bytes_read <= 0) {
		printf("Unable to read %s\n", filename);
		exit(EXIT_FAILURE);
	}
	uint16_t sig = 23117;
	#ifdef HAVE_ZLIB
	if (sig != h.SIGNATURE) {
		fclose(fp);
		gzFile fgz = gzopen(filename, "r");
		if (! fgz) {
			printf("gzopen error %s\n", filename);
			exit(EXIT_FAILURE);
		}
		int bytes_read = gzread(fgz, &h, sizeof(struct mz3hdr));
		if (sig != h.SIGNATURE) {
			gzclose(fgz);
			printf("Unable to read compressed mz3 %s\n", filename);
			exit(EXIT_FAILURE);
		}
		*ntri = (int) h.NFACE;
		*nvert = (int) h.NVERT;
		uint32_t skip = h.NSKIP;
		if (skip > 0) {
			#define GZSKIP_BUFSIZE 4096
				unsigned char buf[GZSKIP_BUFSIZE];
				uint32_t remaining = skip;
				while (remaining > 0) {
						size_t n = (remaining > GZSKIP_BUFSIZE) ? GZSKIP_BUFSIZE : remaining;
						int br = gzread(fgz, buf, (unsigned int)n);
						if (br <= 0) {
								printf("Unable to skip %u bytes in compressed file %s\n", skip, filename);
								exit(EXIT_FAILURE);
						}
						remaining -= br;
				}
		}
		uint32_t tribytes = h.NFACE * sizeof(vec3i);
		*tris = (vec3i *) malloc(tribytes);
		void * imgRaw = (void *) *tris;
		bytes_read = gzread(fgz, imgRaw, tribytes);
		if (bytes_read <= 0) {
			printf("Unable to read compressed triangles %s\n", filename);
			exit(EXIT_FAILURE);
		}
		uint32_t vertbytes32 = h.NVERT * 3 * sizeof(float);
		float * verts32 = (float *) malloc(vertbytes32);
		bytes_read = gzread(fgz, verts32, vertbytes32);
		if (bytes_read <= 0) {
			printf("Unable to read vertices %s\n", filename);
			exit(EXIT_FAILURE);
		}
		*verts = (vec3d *) malloc(h.NVERT * sizeof(vec3d));
		vec3d *vs = *verts;
		int j = 0;
		for (int i = 0; i < h.NVERT; i++) {
			vec3d v;
			v.x = verts32[j++];
			v.y = verts32[j++];
			v.z = verts32[j++];
			vs[i] = v;
		}
		free(verts32);
		gzclose(fgz);
		return;
	}
	#endif
	if (sig != h.SIGNATURE) {
		fclose(fp);
		printf("Unable to read mz3 (unable to read gz compressed) %s\n", filename);
		exit(EXIT_FAILURE);
	}
	*ntri = (int) h.NFACE;
	*nvert = (int) h.NVERT;
	uint32_t skip = h.NSKIP;
	fseek(fp, (int)(sizeof(struct mz3hdr) + skip), SEEK_SET);
	uint32_t tribytes = h.NFACE * sizeof(vec3i);
	*tris = (vec3i *) malloc(tribytes);
	void * imgRaw = (void *) *tris;
	bytes_read = fread(imgRaw, tribytes, 1, fp);
	if (bytes_read <= 0) {
		printf("Unable to read triangles %s\n", filename);
		exit(EXIT_FAILURE);
	}
	uint32_t vertbytes32 = h.NVERT * 3 * sizeof(float);
	float * verts32 = (float *) malloc(vertbytes32);
	bytes_read = fread(verts32, vertbytes32, 1, fp);
	if (bytes_read <= 0) {
		printf("Unable to read vertices %s\n", filename);
		exit(EXIT_FAILURE);
	}
	*verts = (vec3d *) malloc(h.NVERT * sizeof(vec3d));
	vec3d *vs = *verts;
	int j = 0;
	for (int i = 0; i < h.NVERT; i++) {
		vec3d v;
		v.x = verts32[j++];
		v.y = verts32[j++];
		v.z = verts32[j++];
		vs[i] = v;
	}
	free(verts32);
	fclose(fp);
}

int simplify_mz3(const char * innm, const char * outnm, float reduceFraction, bool verbose, int quality) {
	vec3d *pts = NULL;
	vec3i *tris = NULL;
	int ntri, npt;
	read_mz3(innm, &pts, &tris, &npt, &ntri);
	double aggressiveness = 7.0; //7 = default for Simplify.h
	if (quality == 0) //fast
		aggressiveness = 8.0;
	if (quality == 2) //best
		aggressiveness = 5.0;
	int startVert = npt;
	int startTri = ntri;
	int target_count = round((float)ntri * reduceFraction);
	double startTime = clockMsec();
	quadric_simplify_mesh(&pts, &tris, &npt, &ntri, target_count, aggressiveness, verbose, (quality > 1));
	if (verbose)
		printf("simplify vertices %d->%d triangles %d->%d (r = %g): %ld ms\n", startVert, npt, startTri, ntri, (float)ntri / (float) startTri, timediff(startTime, clockMsec()));
	save_mesh(outnm, tris, pts, ntri, npt, (quality > 0));
	free(tris);
	free(pts);
	return EXIT_SUCCESS;
}

int mainMz3(int argc,char **argv) {
	int quality = 1;
	float reduceFraction = 0.25;
	bool verbose = true;
	if (argc > 3) {
		for (int i=2;i<(argc-1);i++) {
			if (strcmp(argv[i],"-q") == 0)
				quality = atoi(argv[i+1]);
			if (strcmp(argv[i],"-r") == 0)
				reduceFraction = atof(argv[i+1]);
			if (strcmp(argv[i],"-v") == 0)
				verbose = atoi(argv[i+1]);
		}
	}
	if ((reduceFraction <= 0.0) || (reduceFraction >= 1.0)) {
		fprintf(stderr,"Mesh reduction factor should be > 0 and < 1.\n");
		return(EXIT_FAILURE);
	}
	return simplify_mz3(argv[1], argv[argc-1], reduceFraction, verbose, quality);
}
#endif // NII2MESH

#ifndef __EMSCRIPTEN__
int show_help( void ) {
	printf("Chris Rorden's niimath version %s (%llu-bit %s)\n",kMTHvers, (unsigned long long) sizeof(size_t)*8, kOS);
    //printf("Chris Rorden's niimath version %s (%llu-bit %s)\n", kMATHvers, (unsigned long long) sizeof(size_t)*8, kOS);
    printf("    Math for NIfTI images inspired by fslmaths without encumbrance problems\n\n");
	printf("\n");
	printf("Usage: niimath [-dt <datatype>] <first_input> [operations and inputs] <output> [-odt <datatype>]\n");
	printf("\n");
	printf("Datatype information:\n");
	printf(" -dt sets the datatype used internally for calculations (default float for all except double images)\n");
	printf(" -odt sets the output datatype ( default is float )\n");
	printf(" Possible datatypes are: char short int float double input input_force\n");
	printf(" ""input"" will set the datatype to that of the original image\n");
	printf("\n");
	printf("New operations: (not in fslmaths)\n");
#ifdef HAVE_BUTTERWORTH
	printf(" -bandpass <hp> <lp> <tr> : Butterworth filter, highpass and lowpass in Hz,TR in seconds (zero-phase 2*2nd order filtfilt)\n");
#endif
	printf(" -bptfm <hp> <lp>         : Same as bptf but does not remove mean (emulates fslmaths < 5.0.7)\n");
#ifdef NII2MESH
	printf(" -bwlabel <conn>          : Connected component labelling for non-zero voxels (conn sets neighbors: 6, 18, 26) \n");
#endif
	printf(" -c2h                     : reverse h2c transform\n");
	printf(" -ceil                    : round voxels upwards to the nearest integer\n");
#ifdef HAVE_CONFORM
	printf(" -ras                     : reorder and flip dimensions to RAS orientation\n");
	printf(" -conform                 : reslice to 1mm size in coronal slice direction with 256^3 voxels\n");
	printf(" -comply <nx> <ny> <nz> <dx> <dy> <dz> <f_high> <isLinear> : conform to axial slice with dx*dy*dzmm size and dx*dy*dz voxels. f_high bright clamping (0.98 for top 2%%). Linear (1) or nearest-neighbor (0)\n");
	printf(" -reslice <target>        : reslice to match image 'target' using linear interpolation\n");
	printf(" -reslice_nn <target>     : reslice to match image 'target' using nearest neighbor interpolation\n");
	printf(" -reslice_mask <mask>     : reslice mask to current image using nearest neighbor; set voxels ≤ 0 in mask to minimum intensity\n");

#endif
	printf(" -close <thr> <dx1> <dx2> : morphological close that binarizes with `thr`, dilates with `dx1` and erodes with `dx2` (fills bubbles with `thr`)\n");
	printf(" -crop <tmin> <tsize>     : remove volumes, starts with 0 not 1! Inputting -1 for a size will set it to the full range\n");
	printf(" -dehaze <mode>           : set dark voxels to zero (mode 1..5; higher yields more surviving voxels)\n");
	printf(" -detrend                 : remove linear trend (and mean) from input\n");
	printf(" -demean                  : remove average signal across volumes (requires 4D input)\n");
	printf(" -dilate <thr> <dx>        : morphological bilate binarizes with `thr`, grows up to distance `dx`\n");
	printf(" -edt                     : estimate Euler Distance Transform (distance field). Assumes isotropic input\n");
	printf(" -erode <thr> <dx>        : morphological erode binarizes with `thr`, shrinks within distance `dx`\n");
	printf(" -floor                   : round voxels downwards to the nearest integer\n");
	printf(" -gz <mode>               : NIfTI gzip mode (0=uncompressed, 1=compressed, else FSL environment; default -1)\n");
	printf(" -h2c                     : convert CT scans from 'Hounsfield' to 'Cormack' units to emphasize soft tissue contrast\n");
#ifdef NII2MESH
	printf("\n");
	printf(" The mesh option has multiple sub-options:\n");
	printf(" -mesh                    : meshify requires 'd'ark, 'm'edium, 'b'right or numeric isosurface ('niimath bet -mesh -i d mesh.gii')\n");
	// We should indent the next few help lines to indicate that they are sub-options of -mesh
	// DO NOT USE TABS!!! Use spaces to indent
	printf("    -i <isovalue>            : 'd'ark, 'm'edium, 'b'right or numeric isosurface\n");
	printf("    -a <atlasFile>           : roi based atlas to mesh\n");
	printf("    -b <fillBubbles>         : fill bubbles\n");
	printf("    -l <onlyLargest>         : only largest\n");
	printf("    -o <originalMC>          : original marching cubes\n");
	printf("    -q <quality>             : quality\n");
	printf("    -s <postSmooth>          : post smooth\n");
	printf("    -r <reduceFraction>      : reduce fraction\n");
	printf("    -v <verbose>             : verbose\n");
	// add -hollow <threshold> <wallThickness> if mesh support enabled
	printf(" -hollow <threshold> <thickness> : hollow out a mesh\n");
#endif
	printf(" -mod                     : modulus fractional remainder - same as '-rem' but includes fractions\n");
	printf(" -otsu <mode>             : binarize image using Otsu's method (mode 1..5; higher yields more bright voxels)\n");
	printf(" -power <exponent>        : raise the current image by following exponent\n");
	printf(" -qform <code>            : set qform_code\n");
	printf(" -sform <code>            : set sform_code\n");
	#if defined(_OPENMP)
	printf(" -p <threads>             : set maximum number of parallel threads (to turn on by default 'export AFNI_COMPRESSOR=PIGZ')\n");
	#else
	printf(" -p <threads>             : set maximum number of parallel threads. DISABLED: recompile for OpenMP support\n");
	#endif
	printf(" -resize <X> <Y> <Z> <m>  : grow (>1) or shrink (<1) image. Method <m> (0=nearest,1=linear,2=spline,3=Lanczos,4=Mitchell)\n");
	printf(" -round                   : round voxels to the nearest integer\n");
	printf(" -sedt                     : estimate signed Euler Distance Transform (distance field). Assumes isotropic input\n");
	printf(" -sobel                   : fast edge detection\n");
	printf(" -sobel_binary            : sobel creating binary edge\n");
	printf(" -tensor_2lower           : convert FSL style upper triangle image to NIfTI standard lower triangle order\n");
	printf(" -tensor_2upper           : convert NIfTI standard lower triangle image to FSL style upper triangle order\n");
	printf(" -tensor_decomp_lower     : as tensor_decomp except input stores lower diagonal (AFNI, ANTS, Camino convention)\n");
	printf(" -trunc                   : truncates the decimal value from floating point value and returns integer value\n");
	printf(" -unsharp  <sigma> <scl>  : edge enhancing unsharp mask (sigma in mm, not voxels [1 is typical]; scl is amount [0.5 medium, 1.0 heavy])\n");
	printf(" -dog <sPos> <sNeg>       : difference of gaussian with zero-crossing edges (positive and negative sigma mm)\n");
	printf(" -dogr <sPos> <sNeg>      : as dog, without zero-crossing (raw rather than binarized data)\n");
	printf(" -dogx <sPos> <sNeg>      : as dog, zero-crossing for 2D sagittal slices\n");
	printf(" -dogy <sPos> <sNeg>      : as dog, zero-crossing for 2D coronal slices\n");
	printf(" -dogz <sPos> <sNeg>      : as dog, zero-crossing for 2D axial slices\n");
	printf(" --compare <ref>          : report if images are identical, terminates without saving new image\n");
	printf(" --compare <theshr> <ref> : report if images are identical, terminates without saving, exits success if difference less than thresh\n");
	printf(" filename.nii             : mimic fslhd (can also export to a txt file: 'niimath T1.nii 2> T1.txt') report header and terminate without saving new image\n");
	printf("\n");
	printf("Binary operations:\n");
	printf("  (some inputs can be either an image or a number)\n");
	printf(" -add <input>             : add following input to current image\n");
	printf(" -sub <input>             : subtract following input from current image\n");
	printf(" -mul <input>             : multiply current image by following input\n");
	printf(" -div <input>             : divide current image by following input\n");
	printf(" -rem <number>            : modulus remainder - divide current image by following input and take remainder\n");
	printf(" -mas <file>              : use (following image>0) to mask current image\n");
	printf(" -thr <number>            : use following number to threshold current image (zero anything below the number)\n");
	printf(" -thrp <input>>           : use following percentage (0-100) of ROBUST RANGE to threshold current image (zero anything below the number)\n");
	printf(" -thrP <input>            : use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold below\n");
	printf(" -uthr <number>           : use following number to upper-threshold current image (zero anything above the number)\n");
	printf(" -uthrp <input>           : use following percentage (0-100) of ROBUST RANGE to upper-threshold current image (zero anything above the number)\n");
	printf(" -uthrP <input>           : use following percentage (0-100) of ROBUST RANGE of non-zero voxels and threshold above\n");
	printf(" -clamp <input>           : use following percentage (0-100) of ROBUST RANGE to threshold current image (anything below set to this threshold)\n");
	printf(" -uclamp <input>          : use following percentage (0-100) of ROBUST RANGE to threshold current image (anything above set to this threshold)\n");
	printf(" -max <input>             : take maximum of following input and current image\n");
	printf(" -min <input>             : take minimum of following input and current image\n");
	printf(" -seed <number>           : seed random number generator with following number\n");
	printf(" -restart <file>          : replace the current image with input for future processing operations\n");
	printf(" -save : save the current working image to the input filename\n");
	printf(" -inm <mean>              :  (-i i ip.c) intensity normalisation (per 3D volume mean)\n");
	printf(" -ing <mean>              :  (-I i ip.c) intensity normalisation, global 4D mean)\n");
	printf(" -s <sigma>               : create a gauss kernel of sigma mm and perform mean filtering\n");
	printf("\n");
	printf("Basic unary operations:\n");
	printf(" -exp   : exponential\n");
	printf(" -log   : natural logarithm\n");
	printf(" -sin   : sine function\n");
	printf(" -cos   : cosine function\n");
	printf(" -tan   : tangent function\n");
	printf(" -asin  : arc sine function\n");
	printf(" -acos  : arc cosine function\n");
	printf(" -atan  : arc tangent function\n");
	printf(" -sqr   : square\n");
	printf(" -sqrt  : square root\n");
	printf(" -recip : reciprocal (1/current image)\n");
	printf(" -abs   : absolute value\n");
	printf(" -bin   : use (current image>0) to binarise\n");
	printf(" -binv  : binarise and invert (binarisation and logical inversion)\n");
	printf(" -fillh : fill holes in a binary mask (holes are internal - i.e. do not touch the edge of the FOV)\n");
	printf(" -fillh26 : fill holes using 26 connectivity\n");
	printf(" -index : replace each nonzero voxel with a unique (subject to wrapping) index number\n");
	printf(" -grid <value> <spacing> : add a 3D grid of intensity <value> with grid spacing <spacing>\n");
	printf(" -edge  : edge strength\n");
	printf(" -tfce <H> <E> <connectivity>: enhance with TFCE, e.g. -tfce 2 0.5 6 (maybe change 6 to 26 for skeletons)\n");
	printf(" -tfceS <H> <E> <connectivity> <X> <Y> <Z> <tfce_thresh>: show support area for voxel (X,Y,Z)\n");
	printf(" -nan   : replace NaNs (improper numbers) with 0\n");
	printf(" -nanm  : make NaN (improper number) mask with 1 for NaN voxels, 0 otherwise\n");
	printf(" -rand  : add uniform noise (range 0:1)\n");
	printf(" -randn : add Gaussian noise (mean=0 sigma=1)\n");
	printf(" -range : set the output calmin/max to full data range\n");
	printf("\n");
	printf("Matrix operations:\n");
	printf(" -tensor_decomp : convert a 4D (6-timepoint )tensor image into L1,2,3,FA,MD,MO,V1,2,3 (remaining image in pipeline is FA)\n");
	printf("\n");
	printf("Kernel operations (set BEFORE filtering operation if desired):\n");
	printf(" -kernel 3D : 3x3x3 box centered on target voxel (set as default kernel)\n");
	printf(" -kernel 2D : 3x3x1 box centered on target voxel\n");
	printf(" -kernel box    <size>     : all voxels in a cube of width <size> mm centered on target voxel\n");
	printf(" -kernel boxv   <size>     : all voxels in a cube of width <size> voxels centered on target voxel, CAUTION: size should be an odd number\n");
	printf(" -kernel boxv3  <X> <Y> <Z>: all voxels in a cuboid of dimensions X x Y x Z centered on target voxel, CAUTION: size should be an odd number\n");
	printf(" -kernel gauss  <sigma>    : gaussian kernel (sigma in mm, not voxels)\n");
	printf(" -kernel sphere <size>     : all voxels in a sphere of radius <size> mm centered on target voxel\n");
	printf(" -kernel file   <filename> : use external file as kernel\n");
	printf("\n");
	printf("Spatial Filtering operations: N.B. all options apart from -s use the default kernel or that _previously_ specified by -kernel\n");
	printf(" -dilM    : Mean Dilation of non-zero voxels\n");
	printf(" -dilD    : Maximum Dilation of non-zero voxels (emulating output of fslmaths 6.0.1, max not modal)\n");
	printf(" -dilF    : Maximum filtering of all voxels\n");
	printf(" -dilall  : Apply -dilM repeatedly until the entire FOV is covered\n");
	printf(" -ero     : Erode by zeroing non-zero voxels when zero voxels found in kernel\n");
	printf(" -eroF    : Minimum filtering of all voxels\n");
	printf(" -fmedian : Median Filtering \n");
	printf(" -fmean   : Mean filtering, kernel weighted (conventionally used with gauss kernel)\n");
	printf(" -fmeanu  : Mean filtering, kernel weighted, un-normalized (gives edge effects)\n");
	printf(" -s <sigma> : create a gauss kernel of sigma mm and perform mean filtering\n");
	printf(" -subsamp2  : downsamples image by a factor of 2 (keeping new voxels centered on old)\n");
	printf(" -subsamp2offc  : downsamples image by a factor of 2 (non-centered)\n");
	printf("\n");
	printf("Dimensionality reduction operations:\n");
	printf("  (the ""T"" can be replaced by X, Y or Z to collapse across a different dimension)\n");
	printf(" -Tmean   : mean across time\n");
	printf(" -Tstd    : standard deviation across time\n");
	printf(" -Tmax    : max across time\n");
	printf(" -Tmaxn   : time index of max across time\n");
	printf(" -Tmin    : min across time\n");
	printf(" -Tmedian : median across time\n");
	printf(" -Tperc <percentage> : nth percentile (0-100) of FULL RANGE across time\n");
	printf(" -Tar1    : temporal AR(1) coefficient (use -odt float and probably demean first)\n");
	printf("\n");
	printf("Basic statistical operations:\n");
	printf(" -pval    : Nonparametric uncorrected P-value, assuming timepoints are the permutations; first timepoint is actual (unpermuted) stats image\n");
	printf(" -pval0   : Same as -pval, but treat zeros as missing data\n");
	printf(" -cpval   : Same as -pval, but gives FWE corrected P-values\n");
	printf(" -ztop    : Convert Z-stat to (uncorrected) P\n");
	printf(" -ptoz    : Convert (uncorrected) P to Z\n");
	printf(" -ztopc    : Convert Z-stat to (uncorrected but clamped) P\n");
	printf(" -ptozc    : Convert (uncorrected but clamped) P to Z\n");
	printf(" -rank    : Convert data to ranks (over T dim)\n");
	printf(" -ranknorm: Transform to Normal dist via ranks\n");
	printf("\n");
	printf("Multi-argument operations:\n");
	printf(" -roi <xmin> <xsize> <ymin> <ysize> <zmin> <zsize> <tmin> <tsize> : zero outside roi (using voxel coordinates). Inputting -1 for a size will set it to the full image extent for that dimension\n");
	printf(" -bptf  <hp_sigma> <lp_sigma> : (-t in ip.c) Bandpass temporal filtering; nonlinear highpass and Gaussian linear lowpass (with sigmas in volumes, not seconds); set either sigma<0 to skip that filter\n");
	printf(" -roc <AROC-thresh> <outfile> [4Dnoiseonly] <truth> : take (normally binary) truth and test current image in ROC analysis against truth. <AROC-thresh> is usually 0.05 and is limit of Area-under-ROC measure FP axis. <outfile> is a text file of the ROC curve (triplets of values: FP TP threshold). If the truth image contains negative voxels these get excluded from all calculations. If <AROC-thresh> is positive then the [4Dnoiseonly] option needs to be set, and the FP rate is determined from this noise-only data, and is set to be the fraction of timepoints where any FP (anywhere) is seen, as found in the noise-only 4d-dataset. This is then controlling the FWE rate. If <AROC-thresh> is negative the FP rate is calculated from the zero-value parts of the <truth> image, this time averaging voxelwise FP rate over all timepoints. In both cases the TP rate is the average fraction of truth=positive voxels correctly found\n");
	printf("\n");
	printf("Combining 4D and 3D images:\n");
	printf(" If you apply a Binary operation (one that takes the current image and a new image together), when one is 3D and the other is 4D,\n");
	printf(" the 3D image is cloned temporally to match the temporal dimensions of the 4D image\n");
	printf("\n");
	printf("e.g. niimath inputVolume -add inputVolume2 output_volume\n");
	printf("     niimath inputVolume -add 2.5 output_volume\n");
	printf("     niimath inputVolume -add 2.5 -mul inputVolume2 output_volume\n");
	printf("\n");
	printf("     niimath 4D_inputVolume -Tmean -mul -1 -add 4D_inputVolume demeaned_4D_inputVolume\n");
    return 0;
}
#endif //ifndef __EMSCRIPTEN__

int main(int argc, char * argv[]) {
	//fslmaths in.nii out.nii changes datatype to float, here we retain (similar to earlier versions of fslmaths)
	//fslmsths in.nii -rem 10 out.nii uses integer modulus not fmod
	//fslmaths robust range not fully described, this emulation is close
	//fslmaths ing/inm are listed as "unary" but should be listed as binary
	//"niimath in.nii" for fslhd style output
	#ifdef _WIN32
		_setmode(_fileno(stdin), _O_BINARY);
	#endif
	if( argc == 2 ) { //special case "niimath img.nii" reports header, a bit like fslhd
		nifti_image * nim = nifti_image_read(argv[1], 0); // issue 5: read header, but not data
		if( nim ) {
			nifti_image_infodump(nim);
			nifti_image_free( nim );
			return 0; //minimal command has input and output: "niimath  in.nii  out.nii"
		}
	}
	#ifndef __EMSCRIPTEN__
	if( argc < 3 ) return show_help(); //minimal command has input and output: "niimath  in.nii  out.nii"
	#endif 
#ifdef NII2MESH
	if (isMz3(argv[1])) {
		return mainMz3(argc, argv);
	}
#endif

	int dtCalc = DT_FLOAT32; //data type for calculation
	int ac = 1;
	if( ! strcmp(argv[ac], "-dt") ) {
		if (! strcmp(argv[ac+1], "double") ) {
			dtCalc = DT_FLOAT64;
			//fprintf(stderr,"'-dt' error: Double calculations not yet supported\n");
			//return 1;
		} else if (strcmp(argv[ac+1], "float") ) {
			fprintf(stderr,"'-dt' error: only 'float' or 'double' calculations supported\n");
			return 1;
		}
		ac += 2;
		if( argc < (ac+2) ) return 1; //insufficient arguments remain
	}
	if (dtCalc == DT_FLOAT32)
		return(main32(argc, argv));
	else
	#ifdef HAVE_64BITS
		return(main64(argc, argv));
	#else
	{
		fprintf(stderr,"'-dt' error: only not compiled for 'double' calculations\n");
		return(EXIT_FAILURE);
	}
	#endif

} //main()

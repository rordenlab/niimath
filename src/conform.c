#include "core.h"
#include "print.h"
#include <limits.h>
#include <math.h>
#include <nifti2_io.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef EMSCRIPTEN
	#define _mm_malloc(size, alignment) malloc(size)
	#define _mm_free(ptr) free(ptr)
#else
	#if defined(_MSC_VER)
		// MSVC (x86 or ARM64)
		#include <malloc.h>  // provides _aligned_malloc / _aligned_free
		#define _mm_malloc(size, alignment) _aligned_malloc(size, alignment)
		#define _mm_free(ptr) _aligned_free(ptr)
	#elif defined(__x86_64__) || defined(__SSE__)
		// GCC/Clang on x86
		#include <immintrin.h>
	#elif defined(__aarch64__) || defined(__arm__) || defined(__ARM_NEON)
		// GCC/Clang on ARM: use custom or fallback
		#include "arm_malloc.h"
	#else
		// fallback
		#define _mm_malloc(size, alignment) malloc(size)
		#define _mm_free(ptr) free(ptr)
	#endif
#endif

// conform.py functions follow
// Python->C port of
//  https://github.com/Deep-MI/FastSurfer/blob/dev/FastSurferCNN/data_loader/conform.py
//  Copyright 2019, AI in Medical Imaging, German Center for Neurodegenerative Diseases (DZNE), Bonn
//   http://www.apache.org/licenses/LICENSE-2.0
//  n.b. did not look at FreeSurfer C code due to licensing restrictions

void scalecrop(float *img, size_t voxnum, float dst_min, float dst_max, float src_min, float scale) {
	// Crop the intensity ranges to specific min and max values.
	for (size_t i = 0; i < voxnum; i++) {
		float val = img[i];
		val = dst_min + scale * (val - src_min);
		val = fmax(val, dst_min);
		val = fmin(val, dst_max);
		img[i] = val;
	}
}

void getscale(float *img, size_t voxnum, float dst_min, float dst_max, float f_low, float f_high, float *src_min, float *scale) {
	*scale = 1.0; // default scale
	*src_min = img[0];
	float src_max = img[0];
	for (size_t i = 0; i < voxnum; i++) {
		src_max = fmax(img[i], src_max);
		*src_min = fmin(img[i], *src_min);
	}
	// count num_nonzero_voxels (nz)
	size_t nz = 0;
	for (size_t i = 0; i < voxnum; i++) {
		if (fabs(img[i]) >= 1e-15) {
			nz++;
		}
	}
	if (*src_min < 0.0) {
		printfx("conform: Input image has value(s) below 0.0 !\n");
	}
	printfx("conform input:    min: %f  max: %f\n", *src_min, src_max);
	if ((f_low <= 0.0) && (f_high >= 1.0)) {
		return;
	}
	// Compute histogram
	size_t histosize = 1000;
	float bin_size = (src_max - *src_min) / histosize;
	int *hist = (int *)calloc(histosize, sizeof(int));
	if (!hist)
		return;
	for (size_t i = 0; i < voxnum; i++) {
		float val = img[i];
		int bin = (int)((val - *src_min) / bin_size);
		bin = fmin(bin, histosize - 1);
		hist[bin]++;
	}
	// Compute cumulative sum
	int *cs = (int *)calloc(histosize, sizeof(int));
	if (!cs) {
		free(hist);
		return;
	}
	cs[0] = hist[0];
	for (size_t i = 1; i < histosize; i++) {
		cs[i] = cs[i - 1] + hist[i];
	}
	// Get lower limit
	size_t nth = (size_t)(f_low * voxnum);
	size_t idx = 0;
	while (idx < histosize) {
		if (cs[idx] >= nth) {
			break;
		}
		idx++;
	}
	*src_min = idx * bin_size + *src_min;
	// Get upper limit
	nth = voxnum - (size_t)((1.0 - f_high) * nz);
	idx = 0;
	while (idx < histosize - 1) {
		if (cs[idx + 1] >= nth) {
			break;
		}
		idx++;
	}
	src_max = idx * bin_size + *src_min;
	// Scale
	if (*src_min != src_max) {
		*scale = (dst_max - dst_min) / (src_max - *src_min);
	} else {
		*scale = 1;
	}
	printfx("Rescale:  min: %f  max: %f  scale: %f\n", *src_min, src_max, *scale);
	free(cs);
	free(hist);
}

void voxelIntensityScale(nifti_image *nim, float f_high) {
	size_t nvox = nim->nx * nim->ny * nim->nz;
	float *img = (float *)nim->data;
	float dst_min = 0;
	float dst_max = 255;
	float f_low = 0.0;
	// n.b. fastsurfer conform.py uses f_high = 0.98
	//   mri_convert command line output reports fhi=0.999
	// float f_high = 0.98; //fsl style 2% robust range
	// float f_high = 0.999; //fastsurfer/mri_convert compatibility
	float src_min;
	float scale;
	getscale(img, nvox, dst_min, dst_max, f_low, f_high, &src_min, &scale);
	scalecrop(img, nvox, dst_min, dst_max, src_min, scale);
	nim->scl_slope = 1.0;
	nim->scl_inter = 0.0;
}

vec4 scaleVec4(vec4 a, float b) {
	// https://glmatrix.net/docs/vec4.js.html#line239
	vec4 out;
	out.v[0] = a.v[0] * b;
	out.v[1] = a.v[1] * b;
	out.v[2] = a.v[2] * b;
	out.v[3] = a.v[3] * b;
	return out;
}

vec4 subtractVec4(vec4 a, vec4 b) {
	// https://glmatrix.net/docs/vec4.js.html#line114
	vec4 out;
	out.v[0] = a.v[0] - b.v[0];
	out.v[1] = a.v[1] - b.v[1];
	out.v[2] = a.v[2] - b.v[2];
	out.v[3] = a.v[3] - b.v[3];
	return out;
}

mat44 scaleMat44(const mat44 a, const vec4 v) {
	// https://glmatrix.net/docs/mat4.js.html#line624
	mat44 out;
	float x = v.v[0];
	float y = v.v[1];
	float z = v.v[2];
	out.m[0][0] = a.m[0][0] * x;
	out.m[0][1] = a.m[0][1] * x;
	out.m[0][2] = a.m[0][2] * x;
	out.m[0][3] = a.m[0][3] * x;
	out.m[1][0] = a.m[1][0] * y;
	out.m[1][1] = a.m[1][1] * y;
	out.m[1][2] = a.m[1][2] * y;
	out.m[1][3] = a.m[1][3] * y;
	out.m[2][0] = a.m[2][0] * z;
	out.m[2][1] = a.m[2][1] * z;
	out.m[2][2] = a.m[2][2] * z;
	out.m[2][3] = a.m[2][3] * z;
	out.m[3][0] = a.m[3][0];
	out.m[3][1] = a.m[3][1];
	out.m[3][2] = a.m[3][2];
	out.m[3][3] = a.m[3][3];
	return out;
}

/*void printMat44(mat44 m) {
	printfx("m=[%g %g %g %g; %g %g %g %g; %g %g %g %g; %g %g %g %g]\n",
		m.m[0][0], m.m[0][1],m.m[0][2], m.m[0][3],
		m.m[1][0], m.m[1][1],m.m[1][2], m.m[1][3],
		m.m[2][0], m.m[2][1],m.m[2][2], m.m[2][3],
		m.m[3][0], m.m[3][1],m.m[3][2], m.m[3][3]
		);
}
void printVec4(vec4 v) {
	printfx("v= [%g, %g, %g, %g]\n", v.v[0], v.v[1], v.v[2], v.v[3]);
}*/

void conformVox2Vox(const int *inDims, mat44 *in_affine, const int outDims[3], const float outMM[3], int toRAS, mat44 *out_affine, mat44 *vox2vox, mat44 *inv_vox2vox) {
	mat44 affine;
	memcpy(&affine, in_affine, sizeof(mat44));
	vec4 half = setVec4(inDims[0] / 2.0f, inDims[1] / 2.0f, inDims[2] / 2.0f);
	vec4 Pxyz_c = nifti_vect44mat44_mul(half, affine);
	vec4 delta = setVec4(outMM[0], outMM[1], outMM[2]);
	mat44 Mdc;
	if (toRAS) {
		LOAD_MAT44(Mdc, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0);
	} else {
		LOAD_MAT44(Mdc, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0);
	}
	mat44 MdcD = scaleMat44(Mdc, delta);
	vec4 vol_center = setVec4(outDims[0], outDims[1], outDims[2]);
	vol_center = nifti_vect44mat44_mul(vol_center, MdcD);
	vol_center = scaleVec4(vol_center, 0.5);
	vec4 translate = subtractVec4(Pxyz_c, vol_center);
	*out_affine = MdcD;
	out_affine->m[0][3] = translate.v[0];
	out_affine->m[1][3] = translate.v[1];
	out_affine->m[2][3] = translate.v[2];
	mat44 inv_out_affine = nifti_mat44_inverse(*out_affine);
	*vox2vox = nifti_mat44_mul(inv_out_affine, affine);
	*inv_vox2vox = nifti_mat44_inverse(*vox2vox);
}

mat44 f642f32mat44(const nifti_dmat44 *dmat) {
	mat44 result;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.m[i][j] = (float)dmat->m[i][j]; // Convert double to float
		}
	}
	return result;
}

void resliceVox2Vox(mat44 *in_affine, mat44 *out_affine, mat44 *vox2vox, mat44 *inv_vox2vox) {
	// Compute vox2vox = inv(out_affine) * in_affine
	mat44 inv_out_affine = nifti_mat44_inverse(*out_affine);
	*vox2vox = nifti_mat44_mul(inv_out_affine, *in_affine);
	*inv_vox2vox = nifti_mat44_inverse(*vox2vox);
}

int doReslice(nifti_image *nim, const float outPixDims[3], const int outDims[3], mat44 out_affine, mat44 inv_vox2vox, int isLinear) {
	if (nim->datatype != DT_FLOAT32) {
		printfx("conform only supports FLOAT32 images\n");
		return EXIT_FAILURE;
	}
	int nvoxIn = nim->nx * nim->ny * MAX(nim->nz, 1);
	// set output header
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			nim->sto_xyz.m[i][j] = out_affine.m[i][j];
			nim->qto_xyz.m[i][j] = out_affine.m[i][j];
		}
	}
	// set output header quaternion
	nifti_dmat44_to_quatern(nim->sto_xyz,
							&nim->quatern_b, &nim->quatern_c, &nim->quatern_d,
							&nim->qoffset_x, &nim->qoffset_y, &nim->qoffset_z,
							&nim->dx, &nim->dy, &nim->dz, &nim->qfac);
	nim->dx = outPixDims[0];
	nim->dy = outPixDims[1];
	nim->dz = outPixDims[2];
	nim->scl_slope = 1.0;
	nim->scl_inter = 0.0;
	nim->cal_min = 0.0;
	nim->cal_max = 0.0;
	// reslice data
	float *in_img = (float *)_mm_malloc(nvoxIn * sizeof(float), 64); // alloc for each volume to allow openmp
	float *raw_img = (float *)nim->data;
	memcpy(in_img, raw_img, nvoxIn * sizeof(float));
	int dimX = nim->nx;
	int dimY = nim->ny;
	int dimZ = nim->nz;
	int dimXY = dimX * dimY;
	// set output
	nim->nx = outDims[0];
	nim->ny = outDims[1];
	nim->nz = outDims[2];
	int nvoxOut = outDims[0] * outDims[1] * outDims[2];
	nim->nvox = nvoxOut;
	free(nim->data);												   // Free the memory allocated with calloc
	float *out_img = (float *)_mm_malloc(nvoxOut * sizeof(float), 64); // output image
	memset(out_img, 0, nim->nvox * sizeof(float));					   // zero array
	nim->data = (void *)out_img;
	if (nim->data == NULL) {
		printfx("conform failed to allocate memory\n");
		return EXIT_FAILURE;
	}
	// find minimimum of input image, used to fill output voxels outside bounding box
	float mn = in_img[0];
	for (int i = 0; i < nvoxIn; i++) {
		mn = MIN(in_img[i], mn);
	}
	for (int i = 0; i < nvoxOut; i++) {
		out_img[i] = mn;
	}
	int outXY = outDims[0] * outDims[1];
	if (isLinear) {
		for (int z = 0; z < outDims[2]; z++) {
			for (int y = 0; y < outDims[1]; y++) {
				float ixYZ = y * inv_vox2vox.m[0][1] + z * inv_vox2vox.m[0][2] + inv_vox2vox.m[0][3];
				float iyYZ = y * inv_vox2vox.m[1][1] + z * inv_vox2vox.m[1][2] + inv_vox2vox.m[1][3];
				float izYZ = y * inv_vox2vox.m[2][1] + z * inv_vox2vox.m[2][2] + inv_vox2vox.m[2][3];
				for (int x = 0; x < outDims[0]; x++) {
					float ix = x * inv_vox2vox.m[0][0] + ixYZ;
					float iy = x * inv_vox2vox.m[1][0] + iyYZ;
					float iz = x * inv_vox2vox.m[2][0] + izYZ;
					int fx = (int)floorf(ix);
					int fy = (int)floorf(iy);
					int fz = (int)floorf(iz);
					int cx = fx + 1;
					int cy = fy + 1;
					int cz = fz + 1;
					int outIdx = x + y * outDims[0] + z * outXY;
					if (fx < 0 || fy < 0 || fz < 0) {
						// out_img[outIdx] = mn;
						continue;
					}
					if (cx >= dimX || cy >= dimY || cz >= dimZ) {
						// out_img[outIdx] = mn;
						continue;
					}
					float dx = ix - fx;
					float dy = iy - fy;
					float dz = iz - fz;
					float wx0 = 1.0f - dx, wx1 = dx;
					float wy0 = 1.0f - dy, wy1 = dy;
					float wz0 = 1.0f - dz, wz1 = dz;
					int i000 = fx + fy * dimX + fz * dimXY;
					float v = 0.0f;
					v += in_img[i000] * wx0 * wy0 * wz0;
					v += in_img[i000 + dimXY] * wx0 * wy0 * wz1;
					v += in_img[i000 + dimX] * wx0 * wy1 * wz0;
					v += in_img[i000 + dimX + dimXY] * wx0 * wy1 * wz1;
					v += in_img[i000 + 1] * wx1 * wy0 * wz0;
					v += in_img[i000 + 1 + dimXY] * wx1 * wy0 * wz1;
					v += in_img[i000 + 1 + dimX] * wx1 * wy1 * wz0;
					v += in_img[i000 + 1 + dimX + dimXY] * wx1 * wy1 * wz1;
					out_img[outIdx] = v;
				}
			}
		}
	} else {
		for (int z = 0; z < outDims[2]; z++) {
			for (int y = 0; y < outDims[1]; y++) {
				float ixYZ = y * inv_vox2vox.m[0][1] + z * inv_vox2vox.m[0][2] + inv_vox2vox.m[0][3];
				float iyYZ = y * inv_vox2vox.m[1][1] + z * inv_vox2vox.m[1][2] + inv_vox2vox.m[1][3];
				float izYZ = y * inv_vox2vox.m[2][1] + z * inv_vox2vox.m[2][2] + inv_vox2vox.m[2][3];
				for (int x = 0; x < outDims[0]; x++) {
					int outIdx = x + y * outDims[0] + z * outXY;
					int ix = roundf(x * inv_vox2vox.m[0][0] + ixYZ);
					int iy = roundf(x * inv_vox2vox.m[1][0] + iyYZ);
					int iz = roundf(x * inv_vox2vox.m[2][0] + izYZ);
					if (ix < 0 || iy < 0 || iz < 0 || ix >= dimX || iy >= dimY || iz >= dimZ) {
						// out_img[outIdx] = mn;
						continue;
					}
					out_img[outIdx] = in_img[ix + iy * dimX + iz * dimXY];
				}
			}
		}
	} // nearest neighbor
	_mm_free(in_img);
	return EXIT_SUCCESS;
}

int reslice(nifti_image *nim, nifti_image *nim2, int isLinear) {
	int nvoxIn = nim->nx * nim->ny * MAX(nim->nz, 1);
	int nVolIn = nim->nvox / nvoxIn;
	int nVoxOut = nim2->nx * nim2->ny * MAX(nim2->nz, 1);
	int nVolOut = nim2->nvox / nVoxOut;
	if ((nVolIn != 1) || (nVolOut != 1)) {
		printfx("conform failed: only for 3D data not 4D time series.\n");
		return EXIT_FAILURE;
	}
	mat44 in_affine = f642f32mat44(&nim->sto_xyz);
	const int outDims[3] = {(int)nim2->nx, (int)nim2->ny, (int)nim2->nz};
	const float outPixDims[3] = {nim2->dx, nim2->dy, nim2->dz};
	mat44 out_affine = f642f32mat44(&nim2->sto_xyz);
	mat44 vox2vox, inv_vox2vox;
	resliceVox2Vox(&in_affine, &out_affine, &vox2vox, &inv_vox2vox);
	return doReslice(nim, outPixDims, outDims, out_affine, inv_vox2vox, isLinear);
}

int conform_core(nifti_image *nim, const int outDims[3], const float outPixDims[3], float f_high, int isLinear, int isRAS) {
	int nvoxIn = nim->nx * nim->ny * MAX(nim->nz, 1);
	int nVol = nim->nvox / nvoxIn;
	if (nVol != 1) {
		printfx("conform failed: only for 3D data not 4D time series.\n");
		return EXIT_FAILURE;
	}
	// normalize voxel brightness 0..255
	if (f_high > 0.0)
		voxelIntensityScale(nim, f_high);
	// estimate spatial transform
	const int inDims[3] = {(int)nim->nx, (int)nim->ny, (int)nim->nz};
	mat44 in_affine = f642f32mat44(&nim->sto_xyz);
	mat44 out_affine, vox2vox, inv_vox2vox;
	conformVox2Vox(inDims, &in_affine, outDims, outPixDims, isRAS, &out_affine, &vox2vox, &inv_vox2vox);
	return doReslice(nim, outPixDims, outDims, out_affine, inv_vox2vox, isLinear);
}

int conform(nifti_image *nim) {
	const int outDims[3] = {256, 256, 256};
	const float outPixDims[3] = {1.0, 1.0, 1.0};
	// n.b. the freesurfer 0.999 can lead little soft tissue dynamic range for 7T T1w images with arterial flow artifacts
	//  mri_convert command line output reports fhi=0.999
	//  https://github.com/Deep-MI/FastSurfer/blob/4557d3bc4d9d54ed908cd222030bc038efc54c2a/FastSurferCNN/data_loader/conform.py#L285
	// float f_high = 0.98; //fsl style 2% robust range
	// float f_high = 0.999; //fastsurfer/mri_convert compatibility
	const float f_high = 0.98;
	return conform_core(nim, outDims, outPixDims, f_high, 1, 0); // 1,0: isLinear(true), isRAS(false)
}

int comply(nifti_image *nim, const int outDims[3], const float outPixDims[3], float f_high, int isLinear) {
	return conform_core(nim, outDims, outPixDims, f_high, isLinear, 1); // 1: isRAS(true)
}

// Helper function to find the index of the maximum absolute value in a row
int find_max_index(float row[3]) {
	int max_index = 0;
	for (int i = 1; i < 3; i++) {
		if (fabs(row[i]) > fabs(row[max_index])) {
			max_index = i;
		}
	}
	return max_index;
}

// Function to determine the order and flips to convert to RAS
void get_ras_order(mat44 m44, int perms[3]) {
	float *rows[3] = {m44.m[0], m44.m[1], m44.m[2]};
	int used[3] = {0, 0, 0}; // To track used indices
	for (int i = 0; i < 3; i++) {
		int max_index = find_max_index(rows[i]);
		if (used[max_index]) {
			printfx("Error: Duplicate max index detected. Spatial transform matrix is invalid.\n");
			return;
		}
		used[max_index] = 1;
		// Determine the sign based on the value
		perms[i] = (rows[i][max_index] > 0) ? (max_index + 1) : -(max_index + 1);
	}
}

void permute_affine(mat44 inAff, mat44 *outAff, int inDims[3], int perms[3]) {
	// Initialize the output affine matrix to zero
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			outAff->m[i][j] = 0.0;
		}
	}
	// Determine the corner of the input volume that will become the new origin
	int corner[3] = {0, 0, 0};
	;
	for (int i = 0; i < 3; i++) {
		if (perms[i] > 0)
			continue;
		int perm = abs(perms[i]) - 1;
		corner[perm] = inDims[perm] - 1;
	}
	// printfx("Corner %d %d %d\n", corner[0], corner[1], corner[2]);
	vec4 newCorner = setVec4(corner[0], corner[1], corner[2]);
	vec4 newOrigin = nifti_vect44mat44_mul(newCorner, inAff);
	// Fill the rotation/scaling part of the output affine
	for (int i = 0; i < 3; i++) {
		int srcAxis = abs(perms[i]) - 1;
		float sign = (perms[i] > 0) ? 1.0 : -1.0;
		for (int j = 0; j < 3; j++) {
			outAff->m[j][i] = sign * inAff.m[j][srcAxis];
		}
	}
	// Set the translation (new origin)
	for (int i = 0; i < 3; i++) {
		outAff->m[i][3] = newOrigin.v[i];
	}
	// Set the homogeneous row
	outAff->m[3][3] = 1.0;
	// printfx("%g %g %g %g\n", outAff->m[0][0], outAff->m[0][1], outAff->m[0][2], outAff->m[0][3]);
	// printfx("%g %g %g %g\n", outAff->m[1][0], outAff->m[1][1], outAff->m[1][2], outAff->m[1][3]);
	// printfx("%g %g %g %g\n", outAff->m[2][0], outAff->m[2][1], outAff->m[2][2], outAff->m[2][3]);
	// printfx("%g %g %g %g\n", outAff->m[3][0], outAff->m[3][1], outAff->m[3][2], outAff->m[3][3]);
}

int toRAS(nifti_image *nim) {
	mat44 in_affine = f642f32mat44(&nim->sto_xyz);
	int perms[3];
	get_ras_order(in_affine, perms);
	if ((perms[0] == 1) && (perms[1] == 2) && (perms[2] == 3)) {
		// data already in RAS
		return EXIT_SUCCESS;
	}
	//
	int inDims[3] = {(int)nim->nx, (int)nim->ny, (int)nim->nz};
	// printfx("RAS Order: [%d, %d, %d]\n", perms[0], perms[1], perms[2]);
	int outDims[3];
	for (int i = 0; i < 3; i++) {
		outDims[i] = inDims[abs(perms[i]) - 1];
	}
	// Iterate over all voxels in the output array
	// Rather than lots of mults, strides allow adds for rows, columns and slices
	// https://mrtrix.readthedocs.io/en/dev/getting_started/image_data.html#strides
	int offs[3] = {1, inDims[0], inDims[0] * inDims[1]}; // offset between column, row, slice
	int inStrides[3] = {offs[abs(perms[0]) - 1], offs[abs(perms[1]) - 1], offs[abs(perms[2]) - 1]};
	int inStarts[3] = {0, 0, 0}; // offset between row, column, slice
	for (int i = 0; i < 3; i++) {
		if (perms[i] > 0)
			continue;
		int dim = inDims[abs(perms[i]) - 1];
		inStarts[i] = (dim - 1) * inStrides[i];
		inStrides[i] = -inStrides[i];
		// printf("flipping dim[%d]\n", i);
	}
	// printfx("Strides: [%d %d %d]\n", inStrides[0], inStrides[1], inStrides[2]);
	// printfx("Starts: [%d %d %d]\n", inStarts[0], inStarts[1], inStarts[2]);
	int nvox3D = inDims[0] * inDims[1] * inDims[2];
	int nVol = nim->nvox / nvox3D;
	float *in_img = (float *)_mm_malloc(nvox3D * sizeof(float), 64);
	float *ras_img = (float *)nim->data;
	int mx = -1;
	int mn = 1;
	for (int v = 0; v < nVol; v++) {					 // transpose each volume separately
		size_t dstIndex = 0;							 // volume offset
		memcpy(in_img, ras_img, nvox3D * sizeof(float)); // dest <- src, sz
		int zOffset = inStarts[2];
		for (int z = 0; z < outDims[2]; z++) {
			int yOffset = inStarts[1];
			for (int y = 0; y < outDims[1]; y++) {
				int xOffset = inStarts[0];
				for (int x = 0; x < outDims[0]; x++) {
					int srcIndex = xOffset + yOffset + zOffset;
					mx = MAX(srcIndex, mx);
					mn = MIN(srcIndex, mn);
					ras_img[dstIndex] = in_img[srcIndex];
					dstIndex++;
					xOffset += inStrides[0];
				} // for x: column
				yOffset += inStrides[1];
			} // for y: row
			zOffset += inStrides[2];
		} // for z slice
		ras_img += nvox3D;
	} // for v : volume
	_mm_free(in_img);
	if ((mn != 0) || (mx != (nvox3D - 1))) {
		printf("ERROR expected %d..%d not 0..%d\n", mn, mx, nvox3D - 1);
		return EXIT_FAILURE;
	}
	// save header
	// printfx("in dims = [%lld %lld %lld]\n", nim->nx, nim->ny, nim->nz);
	// printfx("in pixdims = [%g %g %g]\n", nim->dx, nim->dy, nim->dz);
	nim->nx = outDims[0];
	nim->ny = outDims[1];
	nim->nz = outDims[2];
	nim->dim[1] = nim->nx;
	nim->dim[2] = nim->ny;
	nim->dim[3] = nim->nz;
	float inPixDims[3] = {(float)nim->dx, (float)nim->dy, (float)nim->dz};
	nim->dx = inPixDims[abs(perms[0]) - 1];
	nim->dy = inPixDims[abs(perms[1]) - 1];
	nim->dz = inPixDims[abs(perms[2]) - 1];
	// printfx("out dims = [%lld %lld %lld]\n", nim->nx, nim->ny, nim->nz);
	// printfx("out pixdims = [%g %g %g]\n", nim->dx, nim->dy, nim->dz);
	mat44 out_affine;
	permute_affine(in_affine, &out_affine, inDims, perms);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			nim->sto_xyz.m[i][j] = out_affine.m[i][j];
			nim->qto_xyz.m[i][j] = out_affine.m[i][j];
		}
	}
	return EXIT_SUCCESS;
}

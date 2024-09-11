#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "print.h"
#include <nifti2_io.h>
#include "core.h"
#ifdef EMSCRIPTEN
	#define _mm_malloc(size, alignment) malloc(size)
	#define _mm_free(ptr) free(ptr)
#else
	#ifdef __x86_64__
		#include <immintrin.h>
	#else
		#include "arm_malloc.h"
	#endif
#endif

// conform.py functions follow
// Python->C port of 
//  https://github.com/Deep-MI/FastSurfer/blob/dev/FastSurferCNN/data_loader/conform.py
//  Copyright 2019, AI in Medical Imaging, German Center for Neurodegenerative Diseases (DZNE), Bonn
//   http://www.apache.org/licenses/LICENSE-2.0
//  n.b. did not look at FreeSurfer C code due to licensing restrictions

void scalecrop(float* img, size_t voxnum, float dst_min, float dst_max, float src_min, float scale) {
// Crop the intensity ranges to specific min and max values.
	for (size_t i = 0; i < voxnum; i++) {
		float val = img[i];
		val = dst_min + scale * (val - src_min);
		val = fmax(val, dst_min);
		val = fmin(val, dst_max);
		img[i] = val;
	}
}

void getscale(float* img, size_t voxnum, float dst_min, float dst_max, float f_low, float f_high, float* src_min, float* scale) {
	*scale = 1.0; //default scale
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
	int* hist = (int*)calloc(histosize, sizeof(int));
	if (!hist) return;
	for (size_t i = 0; i < voxnum; i++) {
		float val = img[i];
		int bin = (int)((val - *src_min) / bin_size);
		bin = fmin(bin, histosize - 1);
		hist[bin]++;
	}
	// Compute cumulative sum
	int* cs = (int*)calloc(histosize, sizeof(int));
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
	//n.b. fastsurfer conform.py uses f_high = 0.98
	//  mri_convert command line output reports fhi=0.999
	//float f_high = 0.98; //fsl style 2% robust range
	//float f_high = 0.999; //fastsurfer/mri_convert compatibility
	float src_min;
	float scale;
	getscale(img, nvox, dst_min, dst_max, f_low, f_high, &src_min, &scale);
	scalecrop(img, nvox, dst_min, dst_max, src_min, scale);
	nim->scl_slope = 1.0;
	nim->scl_inter = 0.0;
}

vec4 scaleVec4(vec4 a, float b){
	//https://glmatrix.net/docs/vec4.js.html#line239
	vec4 out;
	out.v[0] = a.v[0] * b;
	out.v[1] = a.v[1] * b;
	out.v[2] = a.v[2] * b;
	out.v[3] = a.v[3] * b;
	return out;
}

vec4 subtractVec4(vec4 a, vec4 b) {
	//https://glmatrix.net/docs/vec4.js.html#line114
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

void conformVox2Vox(const int* inDims, mat44* in_affine, const int outDims[3], const float outMM[3], int toRAS, mat44* out_affine, mat44* vox2vox, mat44* inv_vox2vox) {
	mat44 affine;
	memcpy(&affine, in_affine, sizeof(mat44));
	vec4 half = setVec4(inDims[0] / 2.0f, inDims[1] / 2.0f, inDims[2] / 2.0f);
	vec4 Pxyz_c = nifti_vect44mat44_mul(half, affine );
	vec4 delta = setVec4(outMM[0], outMM[1], outMM[2]);
	mat44 Mdc;
	if (toRAS) {
		LOAD_MAT44(Mdc, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0);
	} else {
		LOAD_MAT44(Mdc, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0);
	}
	vec4 dims = setVec4(outDims[0], outDims[1], outDims[2]);
	mat44 MdcD = scaleMat44(Mdc, delta);
	vec4 vol_center = setVec4(outDims[0], outDims[1], outDims[2]);
	vol_center = nifti_vect44mat44_mul(vol_center, MdcD );
	vol_center = scaleVec4(vol_center, 0.5);
	vec4 translate = subtractVec4(Pxyz_c, vol_center);
	*out_affine = MdcD;
	out_affine->m[0][3] = translate.v[0];
	out_affine->m[1][3] = translate.v[1];
	out_affine->m[2][3] = translate.v[2];
	mat44 inv_out_affine = nifti_mat44_inverse(*out_affine);
	*vox2vox = nifti_mat44_mul( inv_out_affine, affine);
	*inv_vox2vox = nifti_mat44_inverse(*vox2vox);
}

mat44 f642f32mat44(const nifti_dmat44* dmat) {
	mat44 result;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			result.m[i][j] = (float)dmat->m[i][j];  // Convert double to float
		}
	}
	return result;
}

int conform_core(nifti_image *nim, const int outDims[3], const float outPixDims[3], float f_high, int isLinear, int isRAS) {
	int nvoxIn = nim->nx * nim->ny * MAX(nim->nz, 1);
	int nVol = nim->nvox / nvoxIn;
	if (nVol != 1)
		return 1;
	//normalize voxel brightness 0..255
	if (f_high > 0.0)
		voxelIntensityScale(nim, f_high);
	//estimate spatial transform
	const int inDims[3] = {(int)nim->nx, (int)nim->ny, (int)nim->nz};
	mat44 in_affine = f642f32mat44(&nim->sto_xyz);
	mat44 out_affine, vox2vox, inv_vox2vox;
	conformVox2Vox(inDims, &in_affine, outDims, outPixDims, isRAS, &out_affine, &vox2vox, &inv_vox2vox);
	//set output header
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			nim->sto_xyz.m[i][j] = out_affine.m[i][j];
			nim->qto_xyz.m[i][j] = out_affine.m[i][j];
		}
	}
	//set output header quaternion
	nifti_dmat44_to_quatern(nim->sto_xyz ,
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
	//reslice data
	float *in_img = (float *)_mm_malloc(nvoxIn * sizeof(float), 64); //alloc for each volume to allow openmp
	float *raw_img = (float *)nim->data;
	memcpy(in_img, raw_img, nvoxIn * sizeof(float));
	int dimX = nim->nx;
	int dimY = nim->ny;
	int dimZ = nim->nz;
	int dimXY = dimX * dimY;
	//set output
	nim->nx = outDims[0];
	nim->ny = outDims[1];
	nim->nz = outDims[2];
	int nvoxOut = outDims[0] * outDims[1] * outDims[2];
	nim->nvox = nvoxOut;
	free(nim->data);  // Free the memory allocated with calloc
	float *out_img = (float *)_mm_malloc(nvoxOut * sizeof(float), 64); //output image
	memset(out_img, 0, nim->nvox * sizeof(float)); //zero array
	nim->data = (void *) out_img;
	if (nim->data == NULL) {
		printfx("conform failed to allocate memory\n");
		return -1; // Return an error code if allocation fails
	}
	int i = -1;
	//n.b. fastsurfer conform uses linear interpolation: "order = 1"
	// likewise, mri_convert reports "Reslicing using trilinear interpolation"
	if (isLinear) {
		for (int z = 0; z < outDims[2]; z++) {
			for (int y = 0; y < outDims[1]; y++) {
				// loop hoisting
				int ixYZ = y * inv_vox2vox.m[0][1] + z * inv_vox2vox.m[0][2] + inv_vox2vox.m[0][3];
				int iyYZ = y * inv_vox2vox.m[1][1] + z * inv_vox2vox.m[1][2] + inv_vox2vox.m[1][3];
				int izYZ = y * inv_vox2vox.m[2][1] + z * inv_vox2vox.m[2][2] + inv_vox2vox.m[2][3];
				for (int x = 0; x < outDims[0]; x++) {
					float ix = x * inv_vox2vox.m[0][0] + ixYZ;
					float iy = x * inv_vox2vox.m[1][0] + iyYZ;
					float iz = x * inv_vox2vox.m[2][0] + izYZ;
					int fx = floor(ix);
					int fy = floor(iy);
					int fz = floor(iz);
					i++;
					if (fx < 0 || fy < 0 || fz < 0) {
						continue;
					}
					int cx = fx + 1;
					int cy = fy + 1;
					int cz = fz + 1;
					if (cx >= dimX || cy >= dimY || cz >= dimZ) {
						continue;
					}
					// residual fractions
					float rcx = ix - fx;
					float rcy = iy - fy;
					float rcz = iz - fz;
					float rfx = 1.0 - rcx;
					float rfy = 1.0 - rcy;
					float rfz = 1.0 - rcz;
					//floor voxel index for all 3 dimensions
					int fff = fx + fy * dimX + fz * dimXY;
					float vx = 0.0;
					vx += in_img[fff] * rfx * rfy * rfz;
					vx += in_img[fff + dimXY] * rfx * rfy * rcz;
					vx += in_img[fff + dimX] * rfx * rcy * rfz;
					vx += in_img[fff + dimX + dimXY] * rfx * rcy * rcz;
					vx += in_img[fff + 1] * rcx * rfy * rfz;
					vx += in_img[fff + 1 + dimXY] * rcx * rfy * rcz;
					vx += in_img[fff + 1 + dimX] * rcx * rcy * rfz;
					vx += in_img[fff + 1 + dimX + dimXY] * rcx * rcy * rcz;
					out_img[i] = vx;
				} // x
			} // y
		} // x
	} else { //if linear else nearest neighbor
		for (int z = 0; z < outDims[2]; z++) {
			for (int y = 0; y < outDims[1]; y++) {
					// loop hoisting
					int ixYZ = y * inv_vox2vox.m[0][1] + z * inv_vox2vox.m[0][2] + inv_vox2vox.m[0][3];
					int iyYZ = y * inv_vox2vox.m[1][1] + z * inv_vox2vox.m[1][2] + inv_vox2vox.m[1][3];
					int izYZ = y * inv_vox2vox.m[2][1] + z * inv_vox2vox.m[2][2] + inv_vox2vox.m[2][3];
					for (int x = 0; x < outDims[0]; x++) {
						int ix = round(x * inv_vox2vox.m[0][0] + ixYZ);
						int iy = round(x * inv_vox2vox.m[1][0] + iyYZ);
						int iz = round(x * inv_vox2vox.m[2][0] + izYZ);
						i++;
						if (ix < 0 || iy < 0 || iz < 0) {
							continue;
						}
						if (ix >= dimX || iy >= dimY || iz >= dimZ) {
							continue;
						}
						//out_img[i] = in_img[voxidx(ix, iy, iz)];
						out_img[i] = in_img[ix + iy * dimX + iz * dimXY];
					} // z
			} // y
		} // x
	} //nearest neighbor
	_mm_free(in_img);
	return 0;
}

int conform(nifti_image *nim) {
	const int outDims[3] = {256, 256, 256};
	const float outPixDims[3] = {1.0, 1.0, 1.0};
	//n.b. the freesurfer 0.999 can lead little soft tissue dynamic range for 7T T1w images with arterial flow artifacts
	// mri_convert command line output reports fhi=0.999
	// https://github.com/Deep-MI/FastSurfer/blob/4557d3bc4d9d54ed908cd222030bc038efc54c2a/FastSurferCNN/data_loader/conform.py#L285
	//float f_high = 0.98; //fsl style 2% robust range
	//float f_high = 0.999; //fastsurfer/mri_convert compatibility
	const float f_high = 0.98;
	return conform_core(nim, outDims, outPixDims, f_high, 1, 0); //1,0: isLinear(true), isRAS(false)
}

int comply(nifti_image *nim, const int outDims[3], const float outPixDims[3], float f_high, int isLinear) {
	return conform_core(nim, outDims, outPixDims, f_high, isLinear, 1); //1: isRAS(true)
}
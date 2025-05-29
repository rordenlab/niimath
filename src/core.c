#include <ctype.h>
#include <float.h> //FLT_EPSILON
#include <limits.h>
#include <math.h>
#include <nifti2_io.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "core.h"
#include "print.h"

#ifdef NII2MESH
	#include "meshify.h"
	#include "quadric.h"
#endif

#ifdef EMSCRIPTEN
	#define _mm_malloc(size, alignment) malloc(size)
	#define _mm_free(ptr) free(ptr)
#else
	#ifdef __aarch64__
		#include "arm_malloc.h"
	#else
		#include <immintrin.h>
	#endif
#endif

#if defined(_OPENMP) //compile with 'OMP=1 make -j'
	#include <omp.h>
#endif

#define xmemcpy memcpy

#define _USE_MATH_DEFINES //microsoft compiler

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

int nii_otsu(int* H, int nBin, int mode, int *dark, int *mid, int *bright) {
//H: Histogram H[0..nBin-1] with each bin storing nuumber of pixels of this brightness
//nBin: number of bins in histogram, e.g. 256 for H[0..255]
//mode: segment and levels 1: 3/4, 2: 2/3 3: 1/2, 4: 1/3, 5: 1/4
//dark/bright/high: set to threshold
	if (nBin <= 0 || nBin > 32767) { // or 46340
		printfx("nii_otsu error: nBin (%d) is too large for safe indexing\n", nBin);
		return 0;
	}
	int thresh = 0;
	*dark = 0;
	*mid = 0;
	*bright = 0;
	double Sum = 0.0;
	for (int v = 0; v < nBin; v++)
		Sum = Sum + H[v];
	if (Sum <= 0)
		return 0;
	double *P = (double*) malloc(nBin * nBin * sizeof(double));
	double *S = (double*) malloc(nBin * nBin * sizeof(double));
	P[0] = H[0];
	S[0] = H[0];
	for (int v = 1; v < nBin; v++) {
		double Prob = H[v]/Sum;
		P[v] = P[v-1]+Prob;
		S[v] = S[v-1]+(v+1)*Prob;
	}
	for (int u = 1; u < nBin; u++) {
		for (int v = u; v < nBin; v++) {
			P[(u*nBin)+v] = P[v]-P[u-1];
			S[(u*nBin)+v] = S[v]-S[u-1];
		}
	}
	//result is eq 29 from Liao
	for (int u = 0; u < nBin; u++) {
		for (int v = u; v < nBin; v++) {
			if (P[(u*nBin)+v] != 0) //avoid divide by zero errors...
				P[(u*nBin)+v] = (S[(u*nBin)+v]*S[(u*nBin)+v]) / P[(u*nBin)+v];
		}
	}
	if ((mode == 1) || (mode == 5)) {
		int lo = (int)(0.25*nBin);
		int mi = (int)(0.50*nBin);
		int hi = (int)(0.75*nBin);
		double max = P[lo] + P[((lo+1)*nBin)+mi] + P[((mi+1)*nBin)+hi] + P[((hi+1)*nBin)+(nBin-1)];
		for (int l = 0; l < (nBin-3); l++) {
			for (int m = l + 1; m < (nBin-2); m++) {
				for (int h = m + 1; h < (nBin-1); h++) {
					//double v = P[0][l]+P[l+1][h]+P[h+1][nBin-1];
					double v = P[l] + P[((l+1)*nBin)+m] + P[((m+1)*nBin)+h] + P[((h+1)*nBin)+(nBin-1)];
					if (v > max) {
						lo = l;
						mi = m;
						hi = h;
						max = v;
					} //new max
				}//for h -> hi
			} //for m -> mi
		} //for l -> low
		if (mode == 1)
			thresh = hi;
		else
			thresh = lo;
		*dark = lo;
		*mid = mi;
		*bright = hi;
	} else if ((mode == 2) || (mode == 4)) {
		int lo = (int)(0.33*nBin);
		int hi = (int)(0.67*nBin);
		double max = P[lo] + P[((lo+1)*nBin)+hi] + P[((hi+1)*nBin)+nBin-1];
		for (int l = 0; l < (nBin-2); l++) {
			for (int h = l + 1; h < (nBin-1); h++) {
				double v = P[l]+P[((l+1)*nBin)+h]+P[((h+1)*nBin)+nBin-1];
				if (v > max) {
					lo = l;
					hi = h;
					max = v;
				} //new max
			}//for h -> hi
		} //for l -> low
		if (mode == 1)
			thresh = hi;
		else
			thresh = lo;
		*dark = lo;
		*mid = thresh;
		*bright = hi;
	} else { //two levels:
		thresh = (int)(0.25*nBin); //nBin / 2;
		double max = P[thresh]+P[((thresh+1)*nBin)+nBin-1];
		//exhaustively search
		for (int i = 0; i < (nBin-1); i++) {
			double v = P[i]+P[((i+1)*nBin)+nBin-1];
			if (v > max) {
				thresh = i;
				max = v;
			}//new max
		}
		*dark = thresh;
		*mid = thresh;
		*bright = thresh;
	}
	free(P);
	free(S);
	return thresh;
}

int nifti_save(nifti_image *nim, const char *postfix, gzModes gzMode) {
	char extnii[5] = ".nii"; /* modifiable, for possible uppercase */
	char exthdr[5] = ".hdr";
	char extimg[5] = ".img";
	char extgz[5] = ".gz";
	//e.g if current filename is img.nii and postfix is "_FA", save file "img_FA.nii"
	char *fname_in = nim->fname;
	char *iname_in = nim->iname;
	int nifti_type_in = nim->nifti_type;
	char *hname = (char *)calloc(sizeof(char), strlen(nim->fname) + strlen(postfix) + 8);
	char *iname = (char *)calloc(sizeof(char), strlen(nim->fname) + strlen(postfix) + 8);
	//char * hext = (char *)calloc(sizeof(char),8);
	//char * iext = (char *)calloc(sizeof(char),8);
	const char *ext; //input extension
	ext = nifti_find_file_extension(nim->fname);
	strcpy(hname, nim->fname);
	hname[strlen(hname) - strlen(ext)] = 0;
	strcat(hname, postfix);
	strcpy(iname, hname);
	//default extension: .nii
	//strcpy(hext, extnii);
	//strcpy(iext, extnii);
	//read environment
	//export FSLOUTPUTTYPE=NIFTI
	const char *key = "FSLOUTPUTTYPE";
	char nii2Key[2] = "2";
	char gzKey[3] = "GZ";
	char pairKey[5] = "PAIR";
	char *value;
	value = getenv(key);
	//n* has precedence, resolve conflicts between ->dim[*] and ->n*
	nim->dim[1] = nim->nx; //e.g. crop, subsamp2offc
	nim->dim[2] = nim->ny; //e.g. subsamp2offc
	nim->dim[3] = nim->nz; //e.g. subsamp2offc
	nim->dim[4] = nim->nt; //e.g. 4D -> 3D operations like mean
	nim->dim[5] = nim->nu;
	nim->dim[6] = nim->nv;
	nim->dim[7] = nim->nw;
	//d* has precedence, resolve conflicts between ->pixdim[*] and d*
	nim->pixdim[1] = nim->dx;
	nim->pixdim[2] = nim->dy;
	nim->pixdim[3] = nim->dz;
	nim->pixdim[4] = nim->dt;
	//set dime[0]
	int maxDim = 1;
	for (int i = 2; i < 8; i++)
		if (nim->dim[i] > 1)
			maxDim = i;
	nim->dim[0] = maxDim;
	nim->ndim = maxDim;
	//nim->dim[0] = 3;
	//nim->dim[4] = 1;
	int isGz = 0;
	int isNifti2 = 0;
	if ((value != NULL) && strstr(value, nii2Key))
		isNifti2 = 1; //NIFTI2_GZ, NIFTI2_PAIR_GZ, NIFTI_GZ, NIFTI_PAIR_GZ
#ifdef HAVE_ZLIB // if compression is requested, make sure of suffix
	if ((value == NULL) || strstr(value, gzKey))
		isGz = 1; //NIFTI2_GZ, NIFTI2_PAIR_GZ, NIFTI_GZ, NIFTI_PAIR_GZ
	if (gzMode == GZ_FALSE)
		isGz = 0;
	if (gzMode == GZ_TRUE)
		isGz = 1;
#endif
	if ((value != NULL) && strstr(value, pairKey)) {
		strcat(hname, exthdr);
		strcat(iname, extimg);
		if (isNifti2)
			nim->nifti_type = NIFTI_FTYPE_NIFTI2_2;
		else
			nim->nifti_type = NIFTI_FTYPE_NIFTI1_2;

		if (isGz)
			strcat(iname, extgz);
	} else {
		strcat(hname, extnii);
		strcat(iname, extnii);
		if (isNifti2)
			nim->nifti_type = NIFTI_FTYPE_NIFTI2_1;
		else
			nim->nifti_type = NIFTI_FTYPE_NIFTI1_1;
		if (isGz) {
			strcat(hname, extgz);
			strcat(iname, extgz);
		}
	}
	//append extensions...
	nim->fname = hname;
	nim->iname = iname;
	nifti_image_write(nim);
	free(hname);
	if (nim->iname != NULL)
		free(iname);
	//return to input names
	nim->fname = fname_in;
	nim->iname = iname_in;
	nim->nifti_type = nifti_type_in;
	return 0;
}

mat44 xform(nifti_image *nim) {
	if ((nim->sform_code == NIFTI_XFORM_UNKNOWN) && (nim->qform_code == NIFTI_XFORM_UNKNOWN)) {
		mat44 m; //4x4 matrix includes translations
		LOAD_MAT44(m, nim->dx, 0.0, 0.0, 0.0, 0.0, nim->dy, 0.0, 0.0, 0.0, 0.0, nim->dz, 0.0);
		return m;
	}
	nifti_dmat44 AA = nim->sto_xyz;
	if (nim->sform_code < nim->qform_code) //give precedence to SForm, like SPM but unlike VTK tools like ANTs
		AA = nim->qto_xyz; //note qform more constrained than sform: quaternions can not store shears, matrices can
	mat44 m; //4x4 matrix includes translations
	LOAD_MAT44(m, AA.m[0][0], AA.m[0][1], AA.m[0][2], AA.m[0][3],
		AA.m[1][0], AA.m[1][1], AA.m[1][2], AA.m[1][3],
		AA.m[2][0], AA.m[2][1], AA.m[2][2], AA.m[2][3]);
	return m;
}

int neg_determ(nifti_image *nim) {
	//returns -1 for negative determinant, +1 for positive
	mat44 AA = xform(nim);
	mat33 m;
	LOAD_MAT33(m, AA.m[0][0], AA.m[0][1], AA.m[0][2], AA.m[1][0], AA.m[1][1], AA.m[1][2], AA.m[2][0], AA.m[2][1], AA.m[2][2]);
	//printf("determ = %g\n", nifti_mat33_determ(m));
	if (nifti_mat33_determ(m) < 0)
		return 1;
	return 0;
} //report if negative determinant, e.g. we don't want negative volume, eg. "brain volume of -1400cc"

nifti_image *nifti_image_read2(const char *hname, int read_data) {
	//in fslmaths 6.0.1 the commands are different, the first preserves cal_min, cal_max
	// fslmaths in out
	// fslmaths in -add 0 out -odt input
	nifti_image *nim = nifti_image_read(hname, read_data);
	if (nim == NULL)
		exit(134);
	nim->cal_min = 0.0;
	nim->cal_max = 0.0;
	//nim->descrip = '';
	char blank_string[128];
	memset(&blank_string[0], 0, sizeof(blank_string));
	memcpy(nim->descrip, blank_string, 79);
	nim->descrip[79] = '\0';
	strcat(nim->descrip, "6.0.5"); //target fslmaths version
	memcpy(nim->aux_file, blank_string, 23);
	nim->aux_file[23] = '\0';
	memcpy(nim->intent_name, blank_string, 15);
	nim->intent_name[15] = '\0';
	return nim;
}

vec4 setVec4(float x, float y, float z) {
	vec4 v = {{x, y, z, 1.0}};
	return v;
}

vec4 nifti_vect44mat44_mul(vec4 v, mat44 m) { //multiply vector * 4x4matrix
	vec4 vO;
	for (int i = 0; i < 4; i++) { //multiply Pcrs * m
		vO.v[i] = 0;
		for (int j = 0; j < 4; j++)
			vO.v[i] += m.m[i][j] * v.v[j];
	}
	return vO;
}

float vertexDisplacement(float x, float y, float z, mat44 m, mat44 m2) {
	//distance between position of voxel [x,y,z] in space m versus space m2
	vec4 vx = setVec4(x, y, z);
	vec4 pos = nifti_vect44mat44_mul(vx, m);
	vec4 pos2 = nifti_vect44mat44_mul(vx, m2);
	return sqrt(sqr(pos.v[0] - pos2.v[0]));
}

float max_displacement_mm(nifti_image *nim, nifti_image *nim2) {
	//examines each corner of two NIfTI images and returns the max difference in vertex location
	// used to detect if two volumes are aligned
	mat44 m = xform(nim);	//4x4 matrix includes translations
	mat44 m2 = xform(nim2); //4x4 matrix includes translations
	float mx = vertexDisplacement(0, 0, 0, m, m2);
	mx = MAX(mx, vertexDisplacement(nim->nx - 1, 0, 0, m, m2));
	mx = MAX(mx, vertexDisplacement(nim->nx - 1, nim->ny - 1, 0, m, m2));
	mx = MAX(mx, vertexDisplacement(nim->nx - 1, nim->ny - 1, nim->nz - 1, m, m2));
	mx = MAX(mx, vertexDisplacement(nim->nx - 1, 0, nim->nz - 1, m, m2));
	mx = MAX(mx, vertexDisplacement(0, nim->ny - 1, 0, m, m2));
	mx = MAX(mx, vertexDisplacement(0, nim->ny - 1, nim->nz - 1, m, m2));
	mx = MAX(mx, vertexDisplacement(0, 0, nim->nz - 1, m, m2));
	return mx;
}

in_hdr set_input_hdr(nifti_image *nim) {
	//remember input datatype, slope and intercept in case user saves back to this
	in_hdr ihdr;
	ihdr.datatype = nim->datatype;
	ihdr.scl_slope = nim->scl_slope;
	ihdr.scl_inter = nim->scl_inter;
	return ihdr;
}

static inline int32_t clamp_i32(double x) {
	if (x > INT32_MAX) return INT32_MAX;
	if (x < INT32_MIN) return INT32_MIN;
	return (int32_t)round(x);
}

static inline int16_t clamp_i16(double x) {
	if (x > INT16_MAX) return INT16_MAX;
	if (x < INT16_MIN) return INT16_MIN;
	return (int16_t)round(x);
}

static inline uint16_t clamp_u16(double x) {
	if (x > UINT16_MAX) return UINT16_MAX;
	if (x < 0.0) return 0;
	return (uint16_t)round(x);
}

static inline uint8_t clamp_u8(double x) {
	if (x > UINT8_MAX) return UINT8_MAX;
	if (x < 0.0) return 0;
	return (uint8_t)round(x);
}

int nifti_image_change_datatype(nifti_image *nim, int dt, in_hdr *ihdr) {
	//returns -1 on failure, 0 if okay
	bool isRescale = (nim->scl_slope != 1.0) || (nim->scl_inter != 0);
	if ((nim->datatype == dt) && (!isRescale))
		return 0; //no change!
	if (nim->nvox < 1)
		return -1;
	if (nim->scl_slope == 0.0f)
		nim->scl_slope = 1.0;
	float scl = nim->scl_slope;
	float inter = nim->scl_inter;
	//if (ihdr->datatype == dt) { //saving BACK to original format, e.g. int16 converted to float32 for calculations and saved back to int16
	//	nim->scl_slope = ihdr->scl_slope;
	//	nim->scl_inter = ihdr->scl_inter;
	//} else {
		nim->scl_slope = 1.0f;
		nim->scl_inter = 0.0f;
	//}
	int idt = nim->datatype; //input datatype
	double *f64 = (double *)nim->data;
	float *f32 = (float *)nim->data;
	uint64_t *u64 = (uint64_t *)nim->data;
	int64_t *i64 = (int64_t *)nim->data;
	uint32_t *u32 = (uint32_t *)nim->data;
	int32_t *i32 = (int32_t *)nim->data;
	uint16_t *u16 = (uint16_t *)nim->data;
	int16_t *i16 = (int16_t *)nim->data;
	uint8_t *u8 = (uint8_t *)nim->data;
	int8_t *i8 = (int8_t *)nim->data;
	int ok = -1;
	if (dt == DT_FLOAT64) {
		if (idt == DT_FLOAT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				f64[i] = (f64[i] * scl) + inter;
			return 0;
		}
		nim->datatype = DT_FLOAT64;
		nim->nbyper = 8;
		if (idt == DT_UINT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				f64[i] = (u64[i] * scl) + inter;
			return 0;
		}
		if (idt == DT_INT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				f64[i] = (i64[i] * scl) + inter;
			return 0;
		}
		
		
		//following change nbyper
		void *dat = (void *)calloc(1, nim->nvox * sizeof(double));
		double *o64 = (double *)dat;
		if (idt == DT_FLOAT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o64[i] = (f32[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_UINT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o64[i] = (u32[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_INT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o64[i] = (i32[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_UINT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o64[i] = (u16[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_INT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o64[i] = (i16[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_UINT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o64[i] = (u8[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_INT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o64[i] = (i8[i] * scl) + inter;
			ok = 0;
		}
		if (ok == 0) {
			free(nim->data);
			nim->data = dat;
			return 0;
		}
		free(dat);
	} //if (dt == DT_FLOAT64
	if (dt == DT_FLOAT32) {
		if (idt == DT_FLOAT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				f32[i] = (f32[i] * scl) + inter;
			return 0;
		}
		float *o32 = (float *)nim->data;
		nim->datatype = DT_FLOAT32;
		nim->nbyper = 4;
		if (idt == DT_UINT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (u32[i] * scl) + inter;
			return 0;
		}
		if (idt == DT_INT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (i32[i] * scl) + inter;
			return 0;
		}
		//following change nbyper
		void *dat = (void *)calloc(1, nim->nvox * sizeof(float));
		o32 = (float *)dat;
		if (idt == DT_FLOAT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (f64[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_UINT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (u64[i] * scl) + inter;
			return 0;
		}
		if (idt == DT_INT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (i64[i] * scl) + inter;
			return 0;
		}
		if (idt == DT_UINT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (u16[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_INT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (i16[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_UINT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (u8[i] * scl) + inter;
			ok = 0;
		}
		if (idt == DT_INT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = (i8[i] * scl) + inter;
			ok = 0;
		}
		if (ok == 0) {
			free(nim->data);
			nim->data = dat;
			return 0;
		}
		free(dat);
	} //if (dt == DT_FLOAT32)
	if (dt == DT_INT32) {
		if (idt == DT_INT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				i32[i] = clamp_i32((i32[i] * scl) + inter);
			return 0;
		}
		int32_t *o32 = (int32_t *)nim->data;
		nim->datatype = DT_INT32;
		nim->nbyper = 4;
		if (idt == DT_UINT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = clamp_i32(round((u32[i] * scl) + inter));
			return 0;
		}
		if (idt == DT_FLOAT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = clamp_i32(round((f32[i] * scl) + inter));
			return 0;
		}
		//following change nbyper
		void *dat = (void *)calloc(1, nim->nvox * sizeof(int32_t));
		o32 = (int32_t *)dat;
		if (idt == DT_FLOAT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = clamp_i32(round(f64[i] * scl) + inter);
			ok = 0;
		}
		if (idt == DT_UINT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = clamp_i32(round((u16[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_INT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = clamp_i32(round((i16[i] * scl) + inter));
			free(nim->data);
			nim->data = dat;
			ok = 0;
		}
		if (idt == DT_UINT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = clamp_i32(round((u8[i] * scl) + inter));
			free(nim->data);
			nim->data = dat;
			ok = 0;
		}
		if (idt == DT_INT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o32[i] = clamp_i32(round((i8[i] * scl) + inter));
			free(nim->data);
			nim->data = dat;
			ok = 0;
		}
		if (ok == 0) {
			free(nim->data);
			nim->data = dat;
			return 0;
		}
		free(dat);
	} //if (dt == DT_INT32)
	if (dt == DT_INT16) {
		if (idt == DT_INT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				i16[i] = clamp_i16((i16[i] * scl) + inter);
			return 0;
		}
		int16_t *o16 = (int16_t *)nim->data;
		nim->datatype = DT_INT16;
		nim->nbyper = 2;
		if (idt == DT_UINT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_i16(round((u16[i] * scl) + inter));
			return 0;
		}
		//following change nbyper
		void *dat = (void *)calloc(1, nim->nvox * sizeof(int16_t));
		o16 = (int16_t *)dat;
		if (idt == DT_FLOAT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_i16(round(f64[i] * scl) + inter);
			ok = 0;
		}
		if (idt == DT_UINT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_i16(round((u32[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_FLOAT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_i16(round((f32[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_UINT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_i16(round((u8[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_INT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_i16(round((i8[i] * scl) + inter));
			ok = 0;
		}
		if (ok == 0) {
			free(nim->data);
			nim->data = dat;
			return 0;
		}
		free(dat);
	} //if (dt == DT_INT16)
	if (dt == DT_UINT16) {
		if (idt == DT_UINT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				u16[i] = clamp_u16((u16[i] * scl) + inter);
			return 0;
		}
		uint16_t *o16 = (uint16_t *)nim->data;
		nim->datatype = DT_UINT16;
		nim->nbyper = 2;
		if (idt == DT_INT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_u16(round((i16[i] * scl) + inter));
			return 0;
		}
		//following change nbyper
		void *dat = (void *)calloc(1, nim->nvox * sizeof(int16_t));
		o16 = (uint16_t *)dat;
		if (idt == DT_FLOAT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_u16(round(f64[i] * scl) + inter);
			ok = 0;
		}
		if (idt == DT_UINT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_u16(round((u32[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_FLOAT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_u16(round((f32[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_UINT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_u16(round((u8[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_INT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o16[i] = clamp_u16(round((i8[i] * scl) + inter));
			ok = 0;
		}
		if (ok == 0) {
			free(nim->data);
			nim->data = dat;
			return 0;
		}
		free(dat);
	} //if (dt == DT_UINT16)
	if (dt == DT_UINT8) {
		if (idt == DT_UINT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				u8[i] = clamp_u8(round((u8[i] * scl) + inter));
			return 0;
		}
		uint8_t *o8 = (uint8_t *)nim->data;
		nim->datatype = DT_UINT8;
		nim->nbyper = 1;
		if (idt == DT_INT8) {
			for (size_t i = 0; i < nim->nvox; i++)
				o8[i] = clamp_u8(round((i8[i] * scl) + inter));
			return 0;
		}
		//following change nbyper
		void *dat = (void *)calloc(1, nim->nvox * sizeof(uint8_t));
		o8 = (uint8_t *)dat;
		if (idt == DT_FLOAT64) {
			for (size_t i = 0; i < nim->nvox; i++)
				o8[i] = clamp_u8(round(f64[i] * scl) + inter);
			ok = 0;
		}
		if (idt == DT_UINT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o8[i] = clamp_u8(round((u16[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_INT16) {
			for (size_t i = 0; i < nim->nvox; i++)
				o8[i] = clamp_u8(round((i16[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_UINT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o8[i] = clamp_u8(round((u32[i] * scl) + inter));
			ok = 0;
		}
		if (idt == DT_FLOAT32) {
			for (size_t i = 0; i < nim->nvox; i++)
				o8[i] = clamp_u8(round((f32[i] * scl) + inter));
			ok = 0;
		}
		if (ok == 0) {
			free(nim->data);
			nim->data = dat;
			return 0;
		}
		free(dat);
	} //if (dt == DT_UINT16)
	printfx("nifti_image_change_datatype: Unsupported datatype %d -> %d\n", idt, dt);
	return ok;
} //nifti_image_change_datatype()

int *make_kernel_file(nifti_image *nim, int *nkernel, char *fin) {
	nifti_image *nim2 = nifti_image_read(fin, 1);
	if (!nim2) {
		printfx("make_kernel_file: failed to read NIfTI image '%s'\n", fin);
		return NULL;
	}
	int x = nim2->nx;
	int y = nim2->ny;
	int z = nim2->nz;
	int xlo = (int)(-x / 2);
	int ylo = (int)(-y / 2);
	int zlo = (int)(-z / 2);
	in_hdr ihdr = set_input_hdr(nim2);
	if (nifti_image_change_datatype(nim2, DT_FLOAT32, &ihdr) != 0) {
		nifti_image_free(nim2);
		return NULL;
	}
	int n = 0;
	float *f32 = (float *)nim2->data;
	double sum = 0.0;
	for (int i = 0; i < nim2->nvox; i++) {
		if (f32[i] == 0)
			continue;
		sum += fabs(f32[i]);
		n++;
	}
	if ((sum == 0.0) || (n == 0))
		return NULL;
	*nkernel = n;
	int *kernel = (int *)_mm_malloc((n * 4) * sizeof(int), 64); //4 values: offset, xpos, ypos, weight
	//for evenly weighted voxels:
	//int kernelWeight = (int)((double)INT_MAX/(double)n); //requires <limits.h>
	double kernelWeight = (double)INT_MAX / sum;
	int vx = -1;
	int i = 0;
	for (int zi = zlo; zi < (zlo + z); zi++)
		for (int yi = ylo; yi < (ylo + y); yi++)
			for (int xi = xlo; xi < (xlo + x); xi++) {
				vx++;
				if (f32[vx] == 0)
					continue;
				kernel[i] = xi + (yi * nim->nx) + (zi * nim->nx * nim->ny);
				kernel[i + n] = xi; //left-right wrap detection
				kernel[i + n + n] = yi; //anterior-posterior wrap detection
				kernel[i + n + n + n] = (int)(kernelWeight * f32[vx]); //kernel height (weight)
				i++;
			}
	nifti_image_free(nim2);
	return kernel;
} //make_kernel_file()

int *make_kernel_sphere(nifti_image *nim, int *nkernel, double mm) {
	// sphere of radius <size> mm centered on target voxel
	mm = fabs(mm);
	if (mm == 0.0)
		return NULL;
	int x = (2 * floor(mm / nim->dx)) + 1;
	int y = (2 * floor(mm / nim->dy)) + 1;
	int z = (2 * floor(mm / nim->dz)) + 1;
	int xlo = (int)(-x / 2);
	int ylo = (int)(-y / 2);
	int zlo = (int)(-z / 2);
	//first pass: determine number of surviving voxels (n)
	int n = 0;
	for (int zi = zlo; zi < (zlo + z); zi++)
		for (int yi = ylo; yi < (ylo + y); yi++)
			for (int xi = xlo; xi < (xlo + x); xi++) {
				float dx = (xi * nim->dx);
				float dy = (yi * nim->dy);
				float dz = (zi * nim->dz);
				float dist = sqrt(dx * dx + dy * dy + dz * dz);
				if (dist > mm)
					continue;
				n++;
			}
	*nkernel = n;
	int *kernel = (int *)_mm_malloc((n * 4) * sizeof(int), 64); //4 values: offset, xpos, ypos, weight
	int kernelWeight = (int)((double)INT_MAX / (double)n);		//requires <limits.h>
	//second pass: fill surviving voxels
	int i = 0;
	for (int zi = zlo; zi < (zlo + z); zi++)
		for (int yi = ylo; yi < (ylo + y); yi++)
			for (int xi = xlo; xi < (xlo + x); xi++) {
				float dx = (xi * nim->dx);
				float dy = (yi * nim->dy);
				float dz = (zi * nim->dz);
				float dist = sqrt(dx * dx + dy * dy + dz * dz);
				if (dist > mm)
					continue;
				kernel[i] = xi + (yi * nim->nx) + (zi * nim->nx * nim->ny);
				kernel[i + n] = xi; //left-right wrap detection
				kernel[i + n + n] = yi; //anterior-posterior wrap detection
				kernel[i + n + n + n] = kernelWeight; //kernel height
				i++;
			}
	return kernel;
}

#ifdef NII2MESH
int nii2mesh (float * img, nifti_image * nim, int originalMC, float isolevel, float reduceFraction, int preSmooth, bool onlyLargest, bool fillBubbles, int postSmooth, bool verbose, char * outnm, int quality) {
	vec3d *pts = NULL;
	vec3i *tris = NULL;
	int ntri, npt;
	if (nim->datatype != DT_FLOAT32) {
		printfx("'-dt double' does not support mesh\n" );
		return EXIT_FAILURE;
	}
	short dim[3] = {(short)nim->nx, (short)nim->ny, (short)nim->nz};
	if (meshify(img, dim, originalMC, isolevel, &tris, &pts, &ntri, &npt, preSmooth, onlyLargest, fillBubbles, verbose) != EXIT_SUCCESS)
		return EXIT_FAILURE;
	float srow_x[4] = {(float)nim->sto_xyz.m[0][0], (float)nim->sto_xyz.m[0][1], (float)nim->sto_xyz.m[0][2], (float)nim->sto_xyz.m[0][3]} ;
	float srow_y[4] = {(float)nim->sto_xyz.m[1][0], (float)nim->sto_xyz.m[1][1], (float)nim->sto_xyz.m[1][2], (float)nim->sto_xyz.m[1][3]} ;
	float srow_z[4] = {(float)nim->sto_xyz.m[2][0], (float)nim->sto_xyz.m[2][1], (float)nim->sto_xyz.m[2][2], (float)nim->sto_xyz.m[2][3]} ;
	apply_sform(tris, pts, ntri, npt, srow_x, srow_y, srow_z);
	double startTime = clockMsec();
	if (postSmooth > 0) {
		laplacian_smoothHC(pts, tris, npt, ntri, 0.1, 0.5, postSmooth, true);
		if (verbose)
			printfx("post-smooth: %ld ms\n", timediff(startTime, clockMsec()));
		startTime = clockMsec();
	}
	if ((reduceFraction < 1.0) || (quality > 1)) { //lossless for high quality
		double aggressiveness = 7.0; //7 = default for Simplify.h
		if (quality == 0) //fast
			aggressiveness = 8.0;
		if (quality == 2) //best
			aggressiveness = 5.0;
		int startVert = npt;
		int startTri = ntri;
		int target_count = round((float)ntri * reduceFraction);
		quadric_simplify_mesh(&pts, &tris, &npt, &ntri, target_count, aggressiveness, verbose, (quality > 1));
		if (verbose)
			printfx("simplify vertices %d->%d triangles %d->%d (r = %g): %ld ms\n", startVert, npt, startTri, ntri, (float)ntri / (float) startTri, timediff(startTime, clockMsec()));
		startTime = clockMsec();
	}
	save_mesh(outnm, tris, pts, ntri, npt, (quality > 0));
	if (verbose)
		printfx("save to disk: %ld ms\n", timediff(startTime, clockMsec()));
	free(tris);
	free(pts);
	return EXIT_SUCCESS;
}
#endif

int nifti_mesh(nifti_image * nim, float darkThresh, float midThresh, float brightThresh, float imgMax, int arg, int argc, char *argv[]) {
	#ifdef NII2MESH
	#define mxStr 1024
	float isolevel = midThresh;
	float reduceFraction = 0.25;
	int preSmooth = true;
	bool onlyLargest = true;
	bool fillBubbles = false;
	int postSmooth = 0;
	int originalMC = 0;
	int quality = 1;
	bool verbose = true;
	char atlasFilename[mxStr] = "";
	for (int i=arg;i<argc;i++) {
		if (strcmp(argv[i],"-a") == 0)
			strcpy(atlasFilename, argv[i+1]);
		if (strcmp(argv[i],"-b") == 0)
			fillBubbles = atoi(argv[i+1]);
		if (strcmp(argv[i],"-i") == 0) {
			if (strlen(argv[i+1]) < 1) continue;
			if (toupper(argv[i+1][0]) == 'D')
				isolevel = darkThresh;
			else if (toupper(argv[i+1][0]) == 'M')
				isolevel = midThresh;
			else if (toupper(argv[i+1][0]) == 'B')
				isolevel = brightThresh;
			else
				isolevel = atof(argv[i+1]);
		}
		if (strcmp(argv[i],"-l") == 0)
			onlyLargest = atoi(argv[i+1]);
		if (strcmp(argv[i],"-o") == 0)
			originalMC = atoi(argv[i+1]);
		if (strcmp(argv[i],"-p") == 0)
			preSmooth = atoi(argv[i+1]);
		if (strcmp(argv[i],"-q") == 0)
			quality = atoi(argv[i+1]);
		if (strcmp(argv[i],"-s") == 0)
			postSmooth = atoi(argv[i+1]);
		if (strcmp(argv[i],"-r") == 0)
			reduceFraction = atof(argv[i+1]);
		if (strcmp(argv[i],"-v") == 0)
			verbose = atoi(argv[i+1]);
	}
	if (verbose)
		printfx("bubbles=%d isolevel=%g preSmooth=%d quality=%d smooth=%d reduction=%g verbose=%d\n",
			fillBubbles, isolevel, preSmooth, quality, postSmooth, reduceFraction, verbose);
	if (strlen(atlasFilename) > 0) {
		int nLabel = trunc(imgMax);
		int nvox = (nim->nx * nim->ny * nim->nz);
		float * img = (float *)nim->data;
		onlyLargest = false;
		if (imgMax < 1.0) {
			printfx("intensity range not consistent with an indexed atlas (maximum intensity %g)\n", imgMax);
			exit(EXIT_FAILURE);
		}
		char basenm[mxStr], ext[mxStr] = "";
		#define kLabelStrLen 32
		typedef struct  {
			char str[kLabelStrLen];
		} tstr;
		tstr *atlasLabels = (tstr *) malloc((nLabel+1) * sizeof(tstr));
		//We need to use a struct to support MSVC C90, with gcc and clang:
		//  char atlasLabels[nLabel+1][kLabelStrLen];
		for (int i = 0; i <= nLabel; i++)
			snprintf (atlasLabels[i].str, kLabelStrLen-1, "%d", i);
		if (strcmp("1", atlasFilename) != 0) {
			FILE *fp = fopen(atlasFilename,"rt");
			if (fp == NULL) {
				printfx("Unable to find atlas names '%s'\n", atlasFilename);
			} else {
				char str[mxStr], s[mxStr];
				while(fgets(str, mxStr, fp)) {
					strncpy(s, strtok(str,";"), mxStr);
					int i = atoi(s);
					if ((i < 0) || (i > nLabel)) continue;
					strncpy(s, strtok(NULL,";"), mxStr);
					int len = snprintf (atlasLabels[i].str, kLabelStrLen-1, "%s.k%d", s, i);
					if (len < 0) exit(EXIT_FAILURE);
					//remove illegal characters, e.g. 'PACo/Pir' -> 'PACo-Pir'
					if (len < 1) continue;
					for (int j = 0; j < len; j++)
						if ((atlasLabels[i].str[j] < 1) || (atlasLabels[i].str[j] == ' ') || (atlasLabels[i].str[j] == ',') || (atlasLabels[i].str[j] == '/') || (atlasLabels[i].str[j] == '\\') || (atlasLabels[i].str[j] == '%') || (atlasLabels[i].str[j] == '*') || (atlasLabels[i].str[j] == 9) || (atlasLabels[i].str[j] == 10) || (atlasLabels[i].str[j] == 11) || (atlasLabels[i].str[j] == 13))
							atlasLabels[i].str[j] = '-';
				}
				fclose(fp);
			}
		}
		//next, parse name and extension for output files
		strcpy(basenm, argv[argc-1]);
		strip_ext(basenm); // ~/file.nii -> ~/file
		if (strlen(argv[argc-1]) > strlen(basenm))
			strcpy(ext, argv[argc-1] + strlen(basenm));
		#if defined(_OPENMP) //compile with 'OMP=1 make -j'
			int maxNumThreads = omp_get_max_threads();
			printfx("Using %d threads\n", maxNumThreads);
			omp_set_num_threads(maxNumThreads);
		#endif
		int partial_OK, nOK;
		#pragma omp parallel private(partial_OK) shared(nOK)
		{
			partial_OK = 0;
			nOK = 0;
			#pragma omp for
			for (int i = 1; i <= nLabel; i++) {
				printfx("%d/%d\n", i, nLabel);
				float * imgbinary = (float *) malloc(nvox*sizeof(float));
				int n1 = 0;
				float lo = i - 0.5;
				float hi = i + 0.5;
				for (int j = 0; j < nvox; j++) {
					int n = 0;
					if ((img[j] > lo) && (img[j] < hi))
						n = 1;
					imgbinary[j] = n;
					n1 += n;
				}
				if (n1 == 0) {
					printfx("Skipping %d: no voxels with this intensity\n", i);
					continue;
				}
				char outnm[mxStr];
				if (snprintf(outnm,sizeof(outnm),"%s%s%s", basenm, atlasLabels[i].str, ext) < 0) exit(EXIT_FAILURE);
				int reti = nii2mesh(imgbinary, nim, originalMC, 0.5, reduceFraction, preSmooth, onlyLargest, fillBubbles, postSmooth, verbose, outnm, quality);
				if (reti == EXIT_SUCCESS)
					partial_OK ++;
				free(imgbinary);
			} //for nLabel
			#pragma omp critical
			{
				nOK += partial_OK;
			}
		}
		free(atlasLabels);
		printfx("Converted %d regions of interest\n", nOK);
		if (nOK == 0)
			return EXIT_FAILURE;
		return EXIT_SUCCESS;
	} else {
		float * img = (float *)nim->data;
		return nii2mesh (img, nim, originalMC, isolevel, reduceFraction, preSmooth,onlyLargest, fillBubbles, postSmooth, verbose, argv[argc-1], quality);
	}
	#else
		printfx("Not compiled for meshify.\n");
		return EXIT_FAILURE;
	#endif
}

int *make_kernel(nifti_image *nim, int *nkernel, int x, int y, int z) {
	//returns voxels in kernel
	x = MAX(1, x);
	y = MAX(1, y);
	z = MAX(1, z);
	if (((x % 2) == 0) || ((y % 2) == 0) || ((z % 2) == 0))
		printfx("Off-center kernel due to even dimensions.\n");
	int n = x * y * z;
	*nkernel = n;
	int *kernel = (int *)_mm_malloc((n * 4) * sizeof(int), 64); //4 values: offset, xpos, ypos, weight
	int xlo = (int)(-x / 2);
	int ylo = (int)(-y / 2);
	int zlo = (int)(-z / 2);
	int i = 0;
	int kernelWeight = (int)((double)INT_MAX / (double)n); //requires <limits.h>
	for (int zi = zlo; zi < (zlo + z); zi++)
		for (int yi = ylo; yi < (ylo + y); yi++)
			for (int xi = xlo; xi < (xlo + x); xi++) {
				//printf("%d %d %d\n", xi,yi,zi);
				kernel[i] = xi + (yi * nim->nx) + (zi * nim->nx * nim->ny);
				kernel[i + n] = xi; //left-right wrap detection
				kernel[i + n + n] = yi; //anterior-posterior wrap detection
				kernel[i + n + n + n] = kernelWeight; //kernel height
				i++;
			}
	return kernel;
}

//box filter, aka nearest neighbor
#define box_support (0.5)
static double box_filter(double t) {
	if ((t > -0.5) && (t <= 0.5))
		return (1.0);
	return (0.0);
}

//triangle filter, aka linear
#define triangle_support (1.0)
static double triangle_filter(double t) {
	if (t < 0.0)
		t = -t;
	if (t < 1.0)
		return (1.0 - t);
	return (0.0);
}

#define B_spline_support (2.0)
static double B_spline_filter(double t) {
	double tt;
	if (t < 0)
		t = -t;
	if (t < 1) {
		tt = t * t;
		return ((.5 * tt * t) - tt + (2.0 / 3.0));
	} else if (t < 2) {
		t = 2 - t;
		return ((1.0 / 6.0) * (t * t * t));
	}
	return (0.0);
}

static double sinc(double x) {
	x *= M_PI;
	if (x != 0)
		return (sin(x) / x);
	return (1.0);
}

#define Lanczos3_support (3.0)
static double Lanczos3_filter(double t) {
	if (t < 0)
		t = -t;
	if (t < 3.0)
		return (sinc(t) * sinc(t / 3.0));
	return (0.0);
}

#define Mitchell_support (2.0)
#define B (1.0 / 3.0)
#define C (1.0 / 3.0)
static double Mitchell_filter(double t) {
	double tt;
	tt = t * t;
	if (t < 0)
		t = -t;
	if (t < 1.0) {
		t = (((12.0 - 9.0 * B - 6.0 * C) * (t * tt)) + ((-18.0 + 12.0 * B + 6.0 * C) * tt) + (6.0 - 2 * B));
		return (t / 6.0);
	} else if (t < 2.0) {
		t = (((-1.0 * B - 6.0 * C) * (t * tt)) + ((6.0 * B + 30.0 * C) * tt) + ((-12.0 * B - 48.0 * C) * t) + (8.0 * B + 24 * C));
		return (t / 6.0);
	}
	return (0.0);
}

CLIST *createFilter(int srcXsize, int dstXsize, int filterMethod) {
	//Schumacher's resampler in Graphics Gems 3, for improvements see
	//see https://github.com/richgel999/imageresampler/blob/master/resampler.cpp
	// method 1 is the default: linear
	double (*filterf)(double t) = triangle_filter;
	double fwidth = triangle_support;
	if (filterMethod == 0) {
		filterf = box_filter;
		fwidth = box_support;
	} else if (filterMethod == 2) {
		filterf = B_spline_filter;
		fwidth = B_spline_support;
	} else if (filterMethod == 3) {
		filterf = Lanczos3_filter;
		fwidth = Lanczos3_support;
	} else if (filterMethod == 4) {
		filterf = Mitchell_filter;
		fwidth = Mitchell_support;
	}
	CLIST *contrib = (CLIST *)calloc(dstXsize, sizeof(CLIST));
	double xscale = (double)dstXsize / (double)srcXsize;
	double width, fscale, weight;
	double center, left, right;
	int i, j, k, n;
	if (xscale < 1.0) { //image reduction requires anti-aliasing
		width = fwidth / xscale;
		fscale = 1.0 / xscale;
		for (i = 0; i < dstXsize; ++i) {
			contrib[i].n = 0;
			contrib[i].p = (CONTRIB *)calloc((int)(width * 2 + 1), sizeof(CONTRIB));
			center = (double)i / xscale;
			left = ceil(center - width);
			right = floor(center + width);
			for (j = (int)left; j <= (int)right; ++j) {
				weight = center - (double)j;
				weight = (*filterf)(weight / fscale) / fscale;
				if (j < 0) {
					n = -j;
				} else if (j >= srcXsize) {
					n = (srcXsize - j) + srcXsize - 1;
				} else {
					n = j;
				}
				k = contrib[i].n++;
				contrib[i].p[k].pixel = n;
				contrib[i].p[k].weight = weight;
			}
		}
	} else { //if shrink else zoom
		for (i = 0; i < dstXsize; ++i) {
			contrib[i].n = 0;
			contrib[i].p = (CONTRIB *)calloc((int)(fwidth * 2 + 1),
											 sizeof(CONTRIB));
			center = (double)i / xscale;
			left = ceil(center - fwidth);
			right = floor(center + fwidth);
			for (j = (int)left; j <= (int)right; ++j) {
				weight = center - (double)j;
				weight = (*filterf)(weight);
				if (j < 0)
					n = -j;
				else if (j >= srcXsize)
					n = (srcXsize - j) + srcXsize - 1;
				else
					n = j;
				k = contrib[i].n++;
				contrib[i].p[k].pixel = n;
				contrib[i].p[k].weight = weight;
			}
		}
	} //if shrink else zoom
	return contrib;
} //createFilter()

//https://raw.githubusercontent.com/afni/afni/b6a9f7a21c1f3231ff09efbd861f8975ad48e525/src/mri_stats.c
/*******************************************************************/
/****    Given p, return x such that Q(x)=p, for 0 < p < 1.     ****/
/****    Q(x) = 1-P(x) = reversed cdf of N(0,1) variable.       ****/
/*******************************************************************/
//qg and qginv are from AFNI, early AFNI code is GPL later is public domain
//  these are likely early routines, but generic formula
//  if one wants to remove GPL, this should be investigated more carefully

double qg(double x) { return 0.5 * erfc(x / 1.414213562373095); }

double qginv(double p) {
	double dp, dx, dt, ddq, dq;
	int newt; /* not Gingrich, but Isaac */

	dp = (p <= 0.5) ? (p) : (1.0 - p); /* make between 0 and 0.5 */

	if (dp <= 6.1172e-39) { /* cut off at 13 sigma */
		dx = 13.0;
		return ((p <= 0.5) ? (dx) : (-dx));
	}
	/**  Step 1:  use 26.2.23 from Abramowitz and Stegun **/
	dt = sqrt(-2.0 * log(dp));
	dx = dt - ((.010328 * dt + .802853) * dt + 2.515517) / (((.001308 * dt + .189269) * dt + 1.432788) * dt + 1.);
	/**  Step 2:  do 3 Newton steps to improve this (uses the math library erfc function) **/
	for (newt = 0; newt < 3; newt++) {
		dq = 0.5 * erfc(dx / 1.414213562373095) - dp;
		ddq = exp(-0.5 * dx * dx) / 2.506628274631000;
		dx = dx + dq / ddq;
	}
	if (dx > 13.0)
		dx = 13.0;
	return ((p <= 0.5) ? (dx) : (-dx)); /* return with correct sign */
}

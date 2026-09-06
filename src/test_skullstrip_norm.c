// Milestone 1 harness: run ss_normalize() on an image and write the byte working
// volume as a NIfTI on the fixed grid, so it can be compared against AFNI's
// -write_spatnorm output. Not part of the shipped binary.
//
//   cc -O2 -o test_ss_norm test_skullstrip_norm.c skullstrip.c nifti_io.c -lz -lm
//   ./test_ss_norm in.nii.gz out.nii.gz

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nifti_io.h"
#include "skullstrip.h"

static float *to_float(nifti_image *nim, long long nv) {
	float *f = (float *)malloc(sizeof(float) * (size_t)nv);
	long long i;
	double sl = (nim->scl_slope != 0.0) ? nim->scl_slope : 1.0;
	double in = nim->scl_inter;
	if (!f)
		return NULL;
	for (i = 0; i < nv; i++) {
		double v;
		switch (nim->datatype) {
		case DT_UINT8: v = ((unsigned char *)nim->data)[i]; break;
		case DT_INT8: v = ((signed char *)nim->data)[i]; break;
		case DT_INT16: v = ((short *)nim->data)[i]; break;
		case DT_UINT16: v = ((unsigned short *)nim->data)[i]; break;
		case DT_INT32: v = ((int *)nim->data)[i]; break;
		case DT_UINT32: v = ((unsigned int *)nim->data)[i]; break;
		case DT_FLOAT32: v = ((float *)nim->data)[i]; break;
		case DT_FLOAT64: v = ((double *)nim->data)[i]; break;
		default:
			fprintf(stderr, "unsupported datatype %d\n", nim->datatype);
			free(f);
			return NULL;
		}
		f[i] = (float)(v * sl + in);
	}
	return f;
}

int main(int argc, char **argv) {
	nifti_image *nim;
	ss_norm n;
	long long nv, nvox;
	float *f;

	if (argc < 3) {
		fprintf(stderr, "usage: %s in.nii[.gz] out.nii[.gz]\n", argv[0]);
		return 1;
	}
	nim = nifti_image_read(argv[1], 1);
	if (!nim) {
		fprintf(stderr, "cannot read %s\n", argv[1]);
		return 1;
	}
	nv = (long long)nim->nvox;
	f = to_float(nim, nv);
	if (!f)
		return 1;
	free(nim->data);
	nim->data = f;
	nim->datatype = DT_FLOAT32;
	nim->nbyper = 4;
	nim->scl_slope = 1.0;
	nim->scl_inter = 0.0;

	if (ss_normalize(nim, nim->datatype, &n)) {
		fprintf(stderr, "ss_normalize failed\n");
		return 1;
	}
	printf("CM       = [%g %g %g]\n", n.icm, n.jcm, n.kcm);
	printf("ktop     = %d   kbot = %d\n", n.ktop, n.kbot);
	printf("support  = %d\n", n.support);
	printf("clip99   = %d\n", n.clip99);
	printf("a        = [%g %g %g]\n", n.ai, n.aj, n.ak);
	printf("b        = [%g %g %g]\n", n.bi, n.bj, n.bk);
	printf("flip     = [%d %d %d]\n", n.fi, n.fj, n.fk);

	// Round-trip mode: pull the working volume back onto the ORIGINAL grid and write
	// that instead, so alignment with the input can be checked directly.
	if (argc > 3 && strcmp(argv[3], "--restore") == 0) {
		float *back = (float *)malloc(sizeof(float) * (size_t)nv);
		if (!back)
			return 1;
		if (ss_restore(&n, n.vol, nim, back, 0)) {
			fprintf(stderr, "ss_restore failed\n");
			return 1;
		}
		free(nim->data);
		nim->data = back;
		if (nifti_set_filenames(nim, argv[2], 0, 1))
			return 1;
		nifti_image_write(nim);
		ss_norm_free(&n);
		return 0;
	}

	// Reuse the input container: replace its grid with the fixed working grid.
	nvox = (long long)SS_NX * SS_NY * SS_NZ;
	free(nim->data);
	nim->data = malloc(sizeof(float) * (size_t)nvox);
	if (!nim->data)
		return 1;
	for (long long i = 0; i < nvox; i++)
		((float *)nim->data)[i] = (float)n.vol[i];
	nim->ndim = nim->dim[0] = 3;
	nim->nx = nim->dim[1] = SS_NX;
	nim->ny = nim->dim[2] = SS_NY;
	nim->nz = nim->dim[3] = SS_NZ;
	nim->nt = nim->dim[4] = 1;
	nim->nu = nim->dim[5] = 1;
	nim->nv = nim->dim[6] = 1;
	nim->nw = nim->dim[7] = 1;
	nim->nvox = (size_t)nvox;
	nim->dx = nim->pixdim[1] = SS_DXYZ;
	nim->dy = nim->pixdim[2] = SS_DXYZ;
	nim->dz = nim->pixdim[3] = SS_DXYZ;
	// RAI grid: +i runs R->L and +j runs A->P, so world x and y DECREASE with index.
	nim->sform_code = NIFTI_XFORM_SCANNER_ANAT;
	nim->qform_code = 0;
	memset(&nim->sto_xyz, 0, sizeof(nim->sto_xyz));
	nim->sto_xyz.m[0][0] = -SS_DXYZ; nim->sto_xyz.m[0][3] = -SS_XORG;
	nim->sto_xyz.m[1][1] = -SS_DXYZ; nim->sto_xyz.m[1][3] = -SS_YORG;
	nim->sto_xyz.m[2][2] = SS_DXYZ;  nim->sto_xyz.m[2][3] = SS_ZORG;
	nim->sto_xyz.m[3][3] = 1.0;
	if (nifti_set_filenames(nim, argv[2], 0, 1))
		return 1;
	nifti_image_write(nim);
	ss_norm_free(&n);
	return 0;
}

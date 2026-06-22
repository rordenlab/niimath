// niimath --dtifit : linear diffusion tensor fit, emulating FSL's dtifit.
//
// Diffusion model: log(S_i) = log(S0) - b_i * g_i' D g_i, solved per voxel by
// ordinary least squares over all volumes (FSL dtifit's default). The 6-volume
// tensor is then decomposed into FA/MD/L*/V* by niimath's existing, validated
// EIG_tsfunc (src/tensor.c) -- the same eigen code reached by -tensor_decomp.
//
// The tensor *fitting* math (design matrix, log-linear least squares) follows
// AFNI's 3dDWItoDT.c linear path (R.W. Cox / D.R. Glen, public domain). We port
// only the linear solution; AFNI's nonlinear NEWUOA path is intentionally omitted.
//
// ponytail: linear-only, no nonlinear/WLS/kurtosis. Unsupported FSL flags are
// rejected with a clear message rather than silently ignored.

#ifdef HAVE_DTIFIT

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core.h"   // nifti_save, nifti_image_change_datatype, set_input_hdr, gzModes
#include "tensor.h" // EIG_tsfunc
#include "dtifit.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define XFLIP_AUTO (-1)

// Checked allocation: dtifit handles whole-brain DWI buffers that can be large
// enough for allocation to realistically fail; abort loudly rather than risk a
// null dereference or partial-state corruption.
static void *xalloc(size_t n, int zero) {
	void *p = zero ? calloc(1, n) : malloc(n);
	if (!p && n) {
		printf("dtifit: out of memory (requested %zu bytes)\n", n);
		exit(EXIT_FAILURE);
	}
	return p;
}

// Read whitespace-separated floats from a text file. Returns count, fills *out
// (caller frees). Returns -1 on error.
static int read_floats(const char *fname, float **out) {
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("dtifit: unable to open '%s'\n", fname);
		return -1;
	}
	int cap = 256, n = 0;
	float *v = (float *)xalloc(cap * sizeof(float), 0);
	double d;
	while (fscanf(fp, "%lf", &d) == 1) {
		if (n >= cap) {
			cap *= 2;
			v = (float *)realloc(v, cap * sizeof(float));
			if (!v) { printf("dtifit: out of memory reading '%s'\n", fname); exit(EXIT_FAILURE); }
		}
		v[n++] = (float)d;
	}
	fclose(fp);
	*out = v;
	return n;
}

// Invert n x n matrix a (row-major) in place via Gauss-Jordan on the augmented
// [a | I] with partial pivoting. Returns 0 on success, 1 if singular.
static int mat_inv(double *a, int n) {
	double *aug = (double *)xalloc((size_t)n * 2 * n * sizeof(double), 0);
	int w = 2 * n;
	for (int r = 0; r < n; r++) {
		for (int c = 0; c < n; c++) aug[r * w + c] = a[r * n + c];
		for (int c = 0; c < n; c++) aug[r * w + n + c] = (r == c) ? 1.0 : 0.0;
	}
	for (int col = 0; col < n; col++) {
		int best = col;
		double bestv = fabs(aug[col * w + col]);
		for (int r = col + 1; r < n; r++) {
			double v = fabs(aug[r * w + col]);
			if (v > bestv) { bestv = v; best = r; }
		}
		if (bestv < 1e-12) { free(aug); return 1; }
		if (best != col)
			for (int k = 0; k < w; k++) {
				double t = aug[col * w + k]; aug[col * w + k] = aug[best * w + k]; aug[best * w + k] = t;
			}
		double d = aug[col * w + col];
		for (int k = 0; k < w; k++) aug[col * w + k] /= d;
		for (int r = 0; r < n; r++) {
			if (r == col) continue;
			double f = aug[r * w + col];
			for (int k = 0; k < w; k++) aug[r * w + k] -= f * aug[col * w + k];
		}
	}
	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++) a[r * n + c] = aug[r * w + n + c];
	free(aug);
	return 0;
}

static double sform_det3(nifti_image *nim) {
	double (*m)[4] = (nim->sform_code > 0) ? nim->sto_xyz.m : nim->qto_xyz.m;
	return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
	     - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
	     + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Load an image and convert to float32 in place.
static nifti_image *read_f32(const char *fname) {
	nifti_image *nim = nifti_image_read(fname, 1);
	if (!nim) {
		printf("dtifit: failed to read '%s'\n", fname);
		return NULL;
	}
	in_hdr ihdr = set_input_hdr(nim);
	if (nifti_image_change_datatype(nim, DT_FLOAT32, &ihdr) != 0) {
		printf("dtifit: failed to convert '%s' to float32\n", fname);
		nifti_image_free(nim);
		return NULL;
	}
	return nim;
}

int nii_dtifit(int argc, char *argv[]) {
	const char *fdwi = NULL, *fmask = NULL, *fbvec = NULL, *fbval = NULL, *obase = NULL;
	int xflip = XFLIP_AUTO;
	// argv[1] is "--dtifit"; parse the rest with FSL dtifit letters.
	for (int i = 2; i < argc; i++) {
		const char *a = argv[i];
		if (!strcmp(a, "-k") && i + 1 < argc) fdwi = argv[++i];
		else if (!strcmp(a, "-m") && i + 1 < argc) fmask = argv[++i];
		else if (!strcmp(a, "-r") && i + 1 < argc) fbvec = argv[++i];
		else if (!strcmp(a, "-b") && i + 1 < argc) fbval = argv[++i];
		else if (!strcmp(a, "-o") && i + 1 < argc) obase = argv[++i];
		else if (!strcmp(a, "-xflip") && i + 1 < argc) {
			const char *v = argv[++i];
			if (!strcmp(v, "auto")) xflip = XFLIP_AUTO;
			else if (!strcmp(v, "0")) xflip = 0;
			else if (!strcmp(v, "1")) xflip = 1;
			else { printf("dtifit: -xflip must be 0, 1, or auto (got '%s')\n", v); return EXIT_FAILURE; }
		} else {
			printf("dtifit: unsupported option '%s'.\n", a);
			printf("  niimath --dtifit implements FSL dtifit's default linear tensor fit only.\n");
			printf("  Supported: -k <dwi> -r <bvec> -b <bval> -o <base> [-m <mask>] [-xflip 0|1|auto]\n");
			printf("  Unsupported FSL features (e.g. --wls, --sse, --kurt, --gradnonlin, --save_tensor) are not available.\n");
			return EXIT_FAILURE;
		}
	}
	if (!fdwi || !fbvec || !fbval || !obase) {
		printf("dtifit: missing required arguments.\n");
		printf("  usage: niimath --dtifit -k <dwi> -r <bvec> -b <bval> -o <base> [-m <mask>] [-xflip 0|1|auto]\n");
		return EXIT_FAILURE;
	}

	// --- gradient table ---
	float *bval = NULL, *bvec = NULL;
	int nb = read_floats(fbval, &bval);
	int ng = read_floats(fbvec, &bvec);
	if (nb < 1 || ng < 1) { free(bval); free(bvec); return EXIT_FAILURE; }

	nifti_image *nim = read_f32(fdwi);
	if (!nim) { free(bval); free(bvec); return EXIT_FAILURE; }
	int nx = nim->nx, ny = nim->ny, nz = nim->nz;
	size_t nvox3D = (size_t)nx * ny * nz;
	int nvol = (int)(nim->nvox / nvox3D);
	if (nb != nvol || ng != 3 * nvol) {
		printf("dtifit: gradient count mismatch: %d volumes, bval=%d, bvec=%d (need bvec=3*nvol)\n", nvol, nb, ng);
		nifti_image_free(nim); free(bval); free(bvec); return EXIT_FAILURE;
	}
	// bvec stored as 3 rows: x[0..nvol-1], y[..], z[..]
	float *bx = bvec, *by = bvec + nvol, *bz = bvec + 2 * nvol;

	// --- FSL determinant x-flip ---
	double det = sform_det3(nim);
	int doflip = (xflip == 1) || (xflip == XFLIP_AUTO && det > 0.0);
	if (doflip)
		for (int i = 0; i < nvol; i++) bx[i] = -bx[i];
	int nb0 = 0;
	float maxb = 0;
	for (int i = 0; i < nvol; i++) { if (bval[i] < 50.0f) nb0++; if (bval[i] > maxb) maxb = bval[i]; }
	printf("dtifit: %d volumes, %d b0, maxb=%g, sform det=%g, x-flip=%s\n",
	       nvol, nb0, (double)maxb, det, doflip ? "yes" : "no");

	// --- design matrix A (nvol x 7) and pseudo-inverse Pinv = (A'A)^-1 A' (7 x nvol) ---
	const int P = 7;
	double *A = (double *)xalloc((size_t)nvol * P * sizeof(double), 0);
	for (int i = 0; i < nvol; i++) {
		double b = bval[i], gx = bx[i], gy = by[i], gz = bz[i];
		double *r = A + (size_t)i * P;
		r[0] = 1.0;
		r[1] = -b * gx * gx;
		r[2] = -2.0 * b * gx * gy;
		r[3] = -2.0 * b * gx * gz;
		r[4] = -b * gy * gy;
		r[5] = -2.0 * b * gy * gz;
		r[6] = -b * gz * gz;
	}
	double M[49];
	for (int j = 0; j < P; j++)
		for (int k = 0; k < P; k++) {
			double s = 0;
			for (int i = 0; i < nvol; i++) s += A[(size_t)i * P + j] * A[(size_t)i * P + k];
			M[j * P + k] = s;
		}
	if (mat_inv(M, P) != 0) {
		printf("dtifit: design matrix is singular (check bvec/bval)\n");
		nifti_image_free(nim); free(bval); free(bvec); free(A); return EXIT_FAILURE;
	}
	double *Pinv = (double *)xalloc((size_t)P * nvol * sizeof(double), 0); // Pinv[k*nvol + i]
	for (int k = 0; k < P; k++)
		for (int i = 0; i < nvol; i++) {
			double s = 0;
			for (int j = 0; j < P; j++) s += M[k * P + j] * A[(size_t)i * P + j];
			Pinv[(size_t)k * nvol + i] = s;
		}
	// self-check: Pinv * A ~= I_7
	double maxerr = 0;
	for (int a = 0; a < P; a++)
		for (int c = 0; c < P; c++) {
			double s = 0;
			for (int i = 0; i < nvol; i++) s += Pinv[(size_t)a * nvol + i] * A[(size_t)i * P + c];
			double e = fabs(s - (a == c ? 1.0 : 0.0));
			if (e > maxerr) maxerr = e;
		}
	if (maxerr > 1e-4) {
		printf("dtifit: pseudo-inverse self-check failed (max |Pinv*A - I| = %g)\n", maxerr);
		nifti_image_free(nim); free(bval); free(bvec); free(A); free(Pinv); return EXIT_FAILURE;
	}

	// --- per-voxel fit ---
	float *dwi = (float *)nim->data;
	float *mask = NULL;
	nifti_image *nmask = NULL;
	if (fmask) {
		nmask = read_f32(fmask);
		if (!nmask) { nifti_image_free(nim); free(bval); free(bvec); free(A); free(Pinv); return EXIT_FAILURE; }
		if (nmask->nx != nx || nmask->ny != ny || nmask->nz != nz) {
			printf("dtifit: mask %dx%dx%d does not match DWI %dx%dx%d\n", (int)nmask->nx, (int)nmask->ny, (int)nmask->nz, nx, ny, nz);
			nifti_image_free(nim); nifti_image_free(nmask); free(bval); free(bvec); free(A); free(Pinv); return EXIT_FAILURE;
		}
		mask = (float *)nmask->data;
	}
	float *tensor = (float *)xalloc(nvox3D * 6 * sizeof(float), 1); // Dxx,Dxy,Dxz,Dyy,Dyz,Dzz
	float *S0 = (float *)xalloc(nvox3D * sizeof(float), 1);
	#pragma omp parallel
	{
		double *Y = (double *)xalloc((size_t)nvol * sizeof(double), 0); // per-thread log-signal scratch (any nvol)
		#pragma omp for schedule(dynamic, 4096)
		for (size_t v = 0; v < nvox3D; v++) {
			if (mask && mask[v] <= 0.0f) continue;
			double p[7] = {0};
			int ok = 1;
			for (int i = 0; i < nvol; i++) {
				float s = dwi[v + (size_t)i * nvox3D];
				if (s <= 0.0f) { ok = 0; break; }
				Y[i] = log((double)s);
			}
			if (!ok) continue;
			for (int k = 0; k < 7; k++) {
				double acc = 0;
				const double *pk = Pinv + (size_t)k * nvol;
				for (int i = 0; i < nvol; i++) acc += pk[i] * Y[i];
				p[k] = acc;
			}
			S0[v] = (float)exp(p[0]);
			for (int k = 0; k < 6; k++) tensor[v + (size_t)k * nvox3D] = (float)p[k + 1]; // planar layout
		}
		free(Y);
	}
	// release DWI data, keep nim as geometry template for saving
	free(nim->data); nim->data = NULL; dwi = NULL;
	if (nmask) nifti_image_free(nmask);
	free(A); free(Pinv); free(bval); free(bvec);

	// --- decomposition: 16 output volumes ---
	// 0:L1 1:L2 2:L3  3-5:V1  6-8:V2  9-11:V3  12:FA 13:MD 14:S0 15:MO
	const int NV = 16;
	float *out = (float *)xalloc(nvox3D * NV * sizeof(float), 1);
	const double s6 = 3.0 * sqrt(6.0);
	for (size_t v = 0; v < nvox3D; v++) {
		float t[6];
		for (int k = 0; k < 6; k++) t[k] = tensor[v + (size_t)k * nvox3D];
		if (t[0] == 0 && t[3] == 0 && t[5] == 0) continue; // unfitted/masked
		float val[14];
		EIG_tsfunc(0, t, val, 1 /*upper triangle*/);
		for (int k = 0; k < 14; k++) out[v + (size_t)k * nvox3D] = val[k];
		out[v + (size_t)14 * nvox3D] = S0[v];
		// MO = 3*sqrt(6) * abc / (a^2+b^2+c^2)^1.5, deviatoric eigenvalues a,b,c
		double l1 = val[0], l2 = val[1], l3 = val[2], md = (l1 + l2 + l3) / 3.0;
		double a = l1 - md, bb = l2 - md, c = l3 - md;
		double nrm = a * a + bb * bb + c * c;
		out[v + (size_t)15 * nvox3D] = (nrm > 1e-20) ? (float)(s6 * (a * bb * c) / pow(nrm, 1.5)) : 0.0f;
	}
	free(S0);

	// --- save ---
	char outfull[4096];
	snprintf(outfull, sizeof(outfull), "%s.nii.gz", obase);
	if (nifti_set_filenames(nim, outfull, 0, 1)) {
		printf("dtifit: cannot set output name '%s'\n", outfull);
		free(tensor); free(out); nifti_image_free(nim); return EXIT_FAILURE;
	}
	nim->datatype = DT_FLOAT32;
	nim->nbyper = 4;
	nim->scl_slope = 1.0f;
	nim->scl_inter = 0.0f;
	gzModes gz = GZ_ENVIRONMENT;

	// cmin/cmax set viewer display range (cal_min/cal_max); 0,0 lets viewers auto-scale.
	#define SAVE_SCALAR(buf, suffix, cmin, cmax)           \
		do {                                               \
			nim->ndim = 3; nim->nt = nim->nu = nim->nv = nim->nw = 1; \
			nim->nvox = nvox3D; nim->data = (void *)(buf); \
			nim->cal_min = (cmin); nim->cal_max = (cmax);  \
			nifti_save(nim, suffix, gz);                   \
		} while (0)
	#define SAVE_VEC3(buf, suffix)                         \
		do {                                               \
			nim->ndim = 4; nim->nt = 3; nim->nu = nim->nv = nim->nw = 1; \
			nim->nvox = nvox3D * 3; nim->data = (void *)(buf); \
			nim->cal_min = -1; nim->cal_max = 1;           \
			nifti_save(nim, suffix, gz);                   \
		} while (0)

	// tensor (6 vols, FSL upper-triangle order)
	nim->ndim = 4; nim->nt = 6; nim->nu = nim->nv = nim->nw = 1;
	nim->nvox = nvox3D * 6; nim->data = (void *)tensor;
	nim->cal_min = 0; nim->cal_max = 0;
	nifti_save(nim, "_tensor", gz);

	SAVE_SCALAR(out + 0 * nvox3D, "_L1", 0, 0);
	SAVE_SCALAR(out + 1 * nvox3D, "_L2", 0, 0);
	SAVE_SCALAR(out + 2 * nvox3D, "_L3", 0, 0);
	SAVE_VEC3(out + 3 * nvox3D, "_V1");
	SAVE_VEC3(out + 6 * nvox3D, "_V2");
	SAVE_VEC3(out + 9 * nvox3D, "_V3");
	SAVE_SCALAR(out + 12 * nvox3D, "_FA", 0, 1);
	SAVE_SCALAR(out + 13 * nvox3D, "_MD", 0, 0);
	SAVE_SCALAR(out + 14 * nvox3D, "_S0", 0, 0);
	SAVE_SCALAR(out + 15 * nvox3D, "_MO", -1, 1);

	nim->data = NULL; // buffers freed below, not owned by nim
	nifti_image_free(nim);
	free(tensor);
	free(out);
	printf("dtifit: wrote %s_{FA,MD,L1,L2,L3,V1,V2,V3,S0,MO,tensor}\n", obase);
	return EXIT_SUCCESS;
}

#endif // HAVE_DTIFIT

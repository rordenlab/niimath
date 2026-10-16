// confounds.c - fMRIPrep confounds table (-confounds); see confounds.h.
#define _USE_MATH_DEFINES // microsoft compiler: gates M_PI in <math.h>, must precede it
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "confounds.h"
#include "resample.h"
#include "moco.h"
#include "print.h"
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

#define CC_PERIOD 128.0    // CompCor cosine high-pass cutoff (s)
#define CC_VARIANCE 0.5    // a/tCompCor: components explaining 50% of variance
#define CC_CROWN 24        // crownCompCor: fixed component count
#define TCC_TOP 98.0       // tCompCor: voxels at or above the 98th tSTD percentile
#define ACC_THR 0.99       // aCompCor probability threshold in BOLD space
#define FD_RADIUS 50.0     // Power FD head radius (mm)
#define RMS_RADIUS 80.0    // Jenkinson RMSD radius (mm)
#define SPIKE_FD 0.5
#define SPIKE_DVARS 1.5

// ------------------------------------------------------------------ output table
typedef struct { char name[48]; double *v; } cc_col;
typedef struct { cc_col *c; int n, cap, nt; } cc_table;

static double *col_add(cc_table *T, const char *fmt, ...) {
	if (T->n == T->cap) {
		int cap = T->cap ? 2 * T->cap : 64;
		cc_col *c = (cc_col *)realloc(T->c, (size_t)cap * sizeof(cc_col));
		if (!c) return NULL;
		T->c = c; T->cap = cap;
	}
	cc_col *c = &T->c[T->n];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(c->name, sizeof(c->name), fmt, ap);
	va_end(ap);
	c->v = (double *)calloc((size_t)T->nt, sizeof(double));
	if (!c->v) return NULL;
	T->n++;
	return c->v;
}

/* niworkflows ExpandModel '(dd1(x))^^2': x, x_derivative1, x_derivative1_power2, x_power2 */
static int col_expand(cc_table *T, const char *name, const double *x) {
	int nt = T->nt;
	double *b = col_add(T, "%s", name), *d = col_add(T, "%s_derivative1", name);
	double *d2 = col_add(T, "%s_derivative1_power2", name), *p2 = col_add(T, "%s_power2", name);
	if (!b || !d || !d2 || !p2) return 1;
	for (int t = 0; t < nt; t++) {
		b[t] = x[t];
		d[t] = t ? x[t] - x[t - 1] : NAN;
		d2[t] = d[t] * d[t];
		p2[t] = x[t] * x[t];
	}
	return 0;
}

static int table_write(const cc_table *T, const char *path) {
	FILE *f = fopen(path, "w");
	if (!f) { printfx("-confounds: cannot write '%s'\n", path); return 1; }
	for (int i = 0; i < T->n; i++) fprintf(f, "%s%s", T->c[i].name, i + 1 < T->n ? "\t" : "\n");
	for (int t = 0; t < T->nt; t++)
		for (int i = 0; i < T->n; i++) {
			double v = T->c[i].v[t];
			if (isnan(v)) fputs("n/a", f); else fprintf(f, "%.10g", v);
			fputc(i + 1 < T->n ? '\t' : '\n', f);
		}
	return fclose(f) != 0;
}

// ------------------------------------------------------------------ small numerics
static double sel_kth(double *a, int64_t n, int64_t k) {   // 3-way quickselect (reorders a)
	int64_t lo = 0, hi = n - 1;
	while (lo < hi) {
		double pivot = a[lo + ((hi - lo) >> 1)];
		int64_t lt = lo, gt = hi, i = lo;
		while (i <= gt) {
			if (a[i] < pivot) { double t = a[lt]; a[lt++] = a[i]; a[i++] = t; }
			else if (a[i] > pivot) { double t = a[gt]; a[gt--] = a[i]; a[i] = t; }
			else i++;
		}
		if (k < lt) hi = lt - 1; else if (k > gt) lo = gt + 1; else break;
	}
	return a[k];
}

/* numpy.percentile of n values (scratch clobbered): linear (lower = 0) or 'lower' (lower = 1) */
static double percentile(double *a, int64_t n, double pct, int lower) {
	double rank = pct / 100.0 * (double)(n - 1);
	int64_t k = (int64_t)rank;
	double v = sel_kth(a, n, k);
	if (lower || k + 1 >= n) return v;
	double hi = a[k + 1];
	for (int64_t i = k + 2; i < n; i++) if (a[i] < hi) hi = a[i];
	return v + (rank - (double)k) * (hi - v);
}

/* Cyclic Jacobi eigen-decomposition of symmetric A (n x n, destroyed); eigenvalues descending in w,
   eigenvectors in the columns of V. */
static void jacobi_eig(double *A, int n, double *w, double *V) {
	for (int i = 0; i < n * n; i++) V[i] = 0.0;
	double scale = 0.0;   /* the convergence test is relative: Gram entries grow with the voxel count */
	for (int i = 0; i < n; i++) { V[i * n + i] = 1.0; scale += A[i * n + i] * A[i * n + i]; }
	for (int sweep = 0; sweep < 100; sweep++) {
		double off = 0.0;
		for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) off += A[p * n + q] * A[p * n + q];
		if (off <= 1e-24 * scale) break;
		for (int p = 0; p < n; p++)
			for (int q = p + 1; q < n; q++) {
				double apq = A[p * n + q];
				if (fabs(apq) < 1e-300) continue;
				double theta = (A[q * n + q] - A[p * n + p]) / (2.0 * apq);
				double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
				double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
				for (int k = 0; k < n; k++) {   /* A = J^T A J */
					double akp = A[k * n + p], akq = A[k * n + q];
					A[k * n + p] = c * akp - s * akq;
					A[k * n + q] = s * akp + c * akq;
				}
				for (int k = 0; k < n; k++) {
					double apk = A[p * n + k], aqk = A[q * n + k];
					A[p * n + k] = c * apk - s * aqk;
					A[q * n + k] = s * apk + c * aqk;
				}
				for (int k = 0; k < n; k++) {
					double vkp = V[k * n + p], vkq = V[k * n + q];
					V[k * n + p] = c * vkp - s * vkq;
					V[k * n + q] = s * vkp + c * vkq;
				}
			}
	}
	for (int i = 0; i < n; i++) w[i] = A[i * n + i];
	for (int i = 0; i < n; i++) {   /* selection sort, descending, columns follow */
		int m = i;
		for (int j = i + 1; j < n; j++) if (w[j] > w[m]) m = j;
		if (m == i) continue;
		double tw = w[i]; w[i] = w[m]; w[m] = tw;
		for (int k = 0; k < n; k++) { double tv = V[k * n + i]; V[k * n + i] = V[k * n + m]; V[k * n + m] = tv; }
	}
}

// ------------------------------------------------------------------ CompCor
/* nipype _cosine_drift: DCT-II columns k = 1..order-1 (the constant is handled separately). */
static int cosine_order(int n, double tr) {
	int order = (int)floor(2.0 * n * (1.0 / CC_PERIOD) * tr);
	return order < 1 ? 1 : order;
}

static double cosine_basis(int n, int k, int t) {
	return sqrt(2.0 / n) * cos(M_PI / n * (t + 0.5) * k);
}

/* nipype compute_noise_components (pre_filter cosine): voxels of Y (nt x n3) in mask, after `skip`
   leading volumes; cosine-filter + demean, unit-variance per voxel, SVD (via the n x n Gram
   matrix).  Writes `ncomp` (fixed) or the 50%-variance count of columns named prefix%02d. */
static int compcor(cc_table *T, const float *Y, int64_t n3, int nt, int skip, double tr,
                   const uint8_t *mask, int ncomp, const char *prefix) {
	int n = nt - skip, order = cosine_order(n, tr);
	int64_t nv = 0;
	for (int64_t v = 0; v < n3; v++) nv += mask[v] != 0;
	if (nv == 0 || n < 2) return 0;   /* failure_mode NaN with no columns: nothing to report */
	double *M = (double *)malloc((size_t)nv * n * sizeof(double));
	double *C = (double *)calloc((size_t)n * n, sizeof(double));
	double *w = (double *)malloc((size_t)n * sizeof(double)), *V = (double *)malloc((size_t)n * n * sizeof(double));
	int64_t *idx = (int64_t *)malloc((size_t)nv * sizeof(int64_t));
	if (!M || !C || !w || !V || !idx) { free(M); free(C); free(w); free(V); free(idx); return 1; }
	nv = 0;
	for (int64_t v = 0; v < n3; v++) if (mask[v]) idx[nv++] = v;
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int64_t i = 0; i < nv; i++) {
		double *x = M + (size_t)i * n, mean = 0.0;
		for (int t = 0; t < n; t++) { x[t] = Y[(size_t)(t + skip) * n3 + idx[i]]; mean += x[t]; }
		mean /= n;
		for (int t = 0; t < n; t++) x[t] -= mean;
		for (int k = 1; k < order; k++) {   /* orthonormal columns: least squares = projection */
			double b = 0.0;
			for (int t = 0; t < n; t++) b += cosine_basis(n, k, t) * x[t];
			for (int t = 0; t < n; t++) x[t] -= b * cosine_basis(n, k, t);
		}
		double ss = 0.0;
		for (int t = 0; t < n; t++) ss += x[t] * x[t];
		double sd = sqrt(ss / n);
		if (!(sd > 0.0)) sd = 1.0;
		for (int t = 0; t < n; t++) x[t] /= sd;
	}
#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic)
#endif
	for (int a = 0; a < n; a++)
		for (int b = a; b < n; b++) {
			double s = 0.0;
			for (int64_t i = 0; i < nv; i++) s += M[(size_t)i * n + a] * M[(size_t)i * n + b];
			C[a * n + b] = C[b * n + a] = s;
		}
	jacobi_eig(C, n, w, V);
	if (ncomp <= 0) {
		double tot = 0.0, cum = 0.0;
		for (int i = 0; i < n; i++) tot += w[i] > 0 ? w[i] : 0;
		ncomp = n;
		for (int i = 0; i < n; i++) {
			cum += (w[i] > 0 ? w[i] : 0) / tot;
			if (cum >= CC_VARIANCE) { ncomp = i + 1; break; }
		}
	}
	if (ncomp > n) ncomp = n;
	int rc = 0;
	for (int c = 0; c < ncomp && !rc; c++) {
		double *col = col_add(T, "%s%02d", prefix, c);
		if (!col) { rc = 1; break; }
		for (int t = 0; t < n; t++) col[t + skip] = V[t * n + c];   /* dummies stay 0, as nipype pads */
	}
	free(M); free(C); free(w); free(V); free(idx);
	return rc;
}

/* nipype TCompCor mask: within `brain`, quadratic-Legendre-detrended tSTD at or above its 98th percentile. */
static int tcompcor_mask(const float *Y, int64_t n3, int nt, int skip, const uint8_t *brain, uint8_t *out) {
	int n = nt - skip;
	int64_t nv = 0;
	for (int64_t v = 0; v < n3; v++) nv += brain[v] != 0;
	double *sd = (double *)malloc((size_t)(nv ? nv : 1) * sizeof(double));
	double *sc = (double *)malloc((size_t)(nv ? nv : 1) * sizeof(double));
	int64_t *idx = (int64_t *)malloc((size_t)(nv ? nv : 1) * sizeof(int64_t));
	if (!sd || !sc || !idx || n < 3) { free(sd); free(sc); free(idx); return 1; }
	/* design [1, P1, P2] on linspace(-1, 1, n); normal equations solved once */
	double *B = (double *)malloc((size_t)3 * n * sizeof(double)), G[3][3] = {{0}}, Gi[3][3];
	if (!B) { free(sd); free(sc); free(idx); return 1; }
	for (int t = 0; t < n; t++) {
		double x = n > 1 ? -1.0 + 2.0 * t / (n - 1) : 0.0;
		B[t] = 1.0; B[n + t] = x; B[2 * n + t] = 0.5 * (3.0 * x * x - 1.0);
	}
	for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++)
		for (int t = 0; t < n; t++) G[a][b] += B[a * n + t] * B[b * n + t];
	double det = G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1]) - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0])
	           + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]);
	for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
		int a1 = (b + 1) % 3, a2 = (b + 2) % 3, b1 = (a + 1) % 3, b2 = (a + 2) % 3;
		Gi[a][b] = (G[a1][b1] * G[a2][b2] - G[a1][b2] * G[a2][b1]) / det;
	}
	nv = 0;
	for (int64_t v = 0; v < n3; v++) if (brain[v]) idx[nv++] = v;
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int64_t i = 0; i < nv; i++) {
		double r[3] = {0, 0, 0}, beta[3] = {0, 0, 0}, ss = 0.0, mean = 0.0;
		for (int t = 0; t < n; t++) {
			double y = Y[(size_t)(t + skip) * n3 + idx[i]];
			for (int a = 0; a < 3; a++) r[a] += B[a * n + t] * y;
		}
		for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) beta[a] += Gi[a][b] * r[b];
		for (int t = 0; t < n; t++) {
			double e = Y[(size_t)(t + skip) * n3 + idx[i]] - beta[0] - beta[1] * B[n + t] - beta[2] * B[2 * n + t];
			mean += e; ss += e * e;
		}
		mean /= n;
		sd[i] = sqrt(fmax(ss / n - mean * mean, 0.0));
	}
	memcpy(sc, sd, (size_t)nv * sizeof(double));
	double thr = nv ? percentile(sc, nv, TCC_TOP, 0) : 0.0;
	memset(out, 0, (size_t)n3);
	for (int64_t i = 0; i < nv; i++) out[idx[i]] = sd[i] >= thr;
	free(sd); free(sc); free(idx); free(B);
	return 0;
}

// ------------------------------------------------------------------ DVARS (nipype compute_dvars)
static int dvars(const float *Y, int64_t n3, int nt, const uint8_t *brain, double *nstd, double *stdz) {
	int64_t nv = 0;
	for (int64_t v = 0; v < n3; v++) nv += brain[v] != 0;
	double *all = (double *)malloc((size_t)nv * nt * sizeof(double));
	double *x = (double *)malloc((size_t)nv * nt * sizeof(double));
	uint8_t *keep = (uint8_t *)malloc((size_t)(nv ? nv : 1));
	if (!all || !x || !keep || nv == 0 || nt < 2) { free(all); free(x); free(keep); return 1; }
	int64_t i = 0;
	for (int64_t v = 0; v < n3; v++) {
		if (!brain[v]) continue;
		for (int t = 0; t < nt; t++) x[(size_t)i * nt + t] = Y[(size_t)t * n3 + v];
		i++;
	}
	memcpy(all, x, (size_t)nv * nt * sizeof(double));
	double scale = 1000.0 / percentile(all, nv * nt, 50.0, 0);   /* intensity_normalization=1000 */
	int64_t nkeep = 0;
#ifdef _OPENMP
	#pragma omp parallel for reduction(+:nkeep)
#endif
	for (int64_t j = 0; j < nv; j++) {
		double *y = x + (size_t)j * nt, *s = all + (size_t)j * nt;   /* s: this voxel's scratch row */
		for (int t = 0; t < nt; t++) y[t] *= scale;
		memcpy(s, y, (size_t)nt * sizeof(double));
		double q75 = percentile(s, nt, 75.0, 1), q25 = percentile(s, nt, 25.0, 1);
		double func_sd = (q75 - q25) / 1.349;
		keep[j] = func_sd > 0.0;   /* remove_zerovariance, variance_tol 0 */
		if (!keep[j]) continue;
		double mean = 0.0, r0 = 0.0, r1 = 0.0;
		for (int t = 0; t < nt; t++) mean += y[t];
		mean /= nt;
		for (int t = 0; t < nt; t++) r0 += (y[t] - mean) * (y[t] - mean);
		for (int t = 1; t < nt; t++) r1 += (y[t] - mean) * (y[t - 1] - mean);
		double ar1 = r0 > 0 ? r1 / r0 : 0.0;   /* Yule-Walker order 1 (nitime AR_est_YW) */
		s[0] = sqrt(2.0 * (1.0 - ar1)) * func_sd;
		nkeep++;
	}
	double sdsum = 0.0;   /* summed serially: a + reduction on doubles would depend on the team size */
	for (int64_t j = 0; j < nv; j++) if (keep[j]) sdsum += all[(size_t)j * nt];
	free(all);
	double diff_sd_mean = sdsum / (double)nkeep;
	nstd[0] = stdz[0] = NAN;
	for (int t = 1; t < nt; t++) {
		double ss = 0.0;
		for (int64_t j = 0; j < nv; j++) {
			if (!keep[j]) continue;
			double d = x[(size_t)j * nt + t] - x[(size_t)j * nt + t - 1];
			ss += d * d;
		}
		nstd[t] = sqrt(ss / (double)nkeep);
		stdz[t] = nstd[t] / diff_sd_mean;
	}
	free(x); free(keep);
	return 0;
}

// ------------------------------------------------------------------ motion (fMRIPrep FSLMotionParams)
static void mat_mul4(const double A[4][4], const double B[4][4], double C[4][4]) {
	double T[4][4];
	for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) {
		T[r][c] = 0.0;
		for (int k = 0; k < 4; k++) T[r][c] += A[r][k] * B[k][c];
	}
	memcpy(C, T, sizeof(T));
}

static void mat_inv4(const double A[4][4], double B[4][4]) {
	nifti_dmat44 m, i;
	memcpy(m.m, A, sizeof(m.m));
	i = nifti_dmat44_inverse(m);
	memcpy(B, i.m, sizeof(i.m));
}

/* Per-volume FSL matrices (nitransforms FSLLinearTransformArray.from_ras of the RAS pull). */
static int fsl_matrices(const nifti_image *nim, const double *hmc, int nt, double (*F)[4][4]) {
	nifti_dmat44 S = (nim->sform_code > 0) ? nim->sto_xyz : nim->qto_xyz;
	double A[4][4], Ai[4][4], swp[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}, spc[4][4] = {{0}};
	memcpy(A, S.m, sizeof(A));
	mat_inv4(A, Ai);
	double z[3] = {nim->dx, nim->dy, nim->dz};
	for (int a = 0; a < 3; a++) spc[a][a] = z[a];
	spc[3][3] = 1.0;
	double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
	           + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
	if (det > 0) { swp[0][0] = -1.0; swp[0][3] = (nim->nx - 1) * z[0]; }
	double swpi[4][4], spci[4][4], pre[4][4], post[4][4];
	mat_inv4(swp, swpi);
	mat_inv4(spc, spci);
	mat_mul4(spci, swpi, pre);     /* FSL ref -> voxel */
	mat_mul4(A, pre, pre);         /* -> RAS world */
	mat_mul4(spc, Ai, post);       /* RAS world -> scaled voxel */
	mat_mul4(swp, post, post);     /* -> FSL moving */
	for (int t = 0; t < nt; t++) {
		double M[3][4], P[4][4] = {{0}}, W[4][4];
		if (moco_pull_affine(nim, hmc + 6 * t, M)) return 1;
		for (int r = 0; r < 3; r++) for (int c = 0; c < 4; c++) P[r][c] = M[r][c];
		P[3][3] = 1.0;
		mat_mul4(A, P, W);         /* RAS pull: ref world -> moving world */
		mat_mul4(W, Ai, W);
		mat_mul4(post, W, W);
		mat_mul4(W, pre, W);       /* FSL ref -> FSL moving (pull) */
		mat_inv4(W, F[t]);         /* nitransforms stores the inverse */
	}
	return 0;
}

// ------------------------------------------------------------------ driver
static uint8_t *load_mask_on(const char *fn, const nifti_image *grid, const double *pull, double thr) {
	nifti_image *m = nii_read_float32(fn);
	if (!m) return NULL;
	if (nii_reject_5d(m, "-confounds")) { nifti_image_free(m); return NULL; }
	int64_t n3 = grid->nx * grid->ny * grid->nz;
	float *v = NULL;
	if (pull) v = resample_volume(m, grid, pull);
	else if (nii_same_grid(m, grid)) {
		v = (float *)malloc((size_t)n3 * sizeof(float));
		if (v) memcpy(v, m->data, (size_t)n3 * sizeof(float));
	} else printfx("-confounds: '%s' is not on the BOLD grid\n", fn);
	nifti_image_free(m);
	if (!v) return NULL;
	uint8_t *b = (uint8_t *)malloc((size_t)n3);
	for (int64_t i = 0; b && i < n3; i++) b[i] = v[i] > thr;
	free(v);
	return b;
}

static void mean_signal(const float *Y, int64_t n3, int nt, const uint8_t *m, double *out) {
	for (int t = 0; t < nt; t++) {
		double s = 0.0;
		int64_t k = 0;
		for (int64_t v = 0; v < n3; v++) if (m[v]) { s += Y[(size_t)t * n3 + v]; k++; }
		out[t] = k ? s / k : NAN;
	}
}

/* -tsnr: per brain voxel, mean / residual SD over the non-dummy volumes after regressing out an intercept,
   Legendre trends up to degree 1 + floor(duration / 150 s) (AFNI -polort A) and the 6 motion parameters.
   The map borrows nim's header (restored); the JSON carries the brain median for reports. */
static int write_tsnr(nifti_image *nim, const float *Y, int64_t n3, int nt, int skip, double tr, const double *hmc,
                      const uint8_t *brain, const char *path, const char *json) {
	int n = nt - skip, deg = 1 + (int)floor(n * tr / 150.0), p = 1 + deg + 6;
	if (n <= p + 1) { printfx("-tsnr: %d volumes cannot support %d regressors\n", n, p); return 1; }
	double *X = (double *)calloc((size_t)p * n, sizeof(double)), *G = (double *)calloc((size_t)p * p * 2, sizeof(double));
	float *map = (float *)calloc((size_t)n3, sizeof(float));
	int64_t nv = 0;
	for (int64_t v = 0; v < n3; v++) nv += brain[v] != 0;
	double *vals = (double *)malloc((size_t)(nv ? nv : 1) * sizeof(double));
	int64_t *idx = (int64_t *)malloc((size_t)(nv ? nv : 1) * sizeof(int64_t));
	int rc = 1;
	if (!X || !G || !map || !vals || !idx) goto done;
	for (int t = 0; t < n; t++) {   /* columns: Legendre P0..Pdeg on [-1, 1], then motion */
		double x = -1.0 + 2.0 * t / (n - 1), pm = 1.0, pk = x;
		X[t] = 1.0;
		X[n + t] = x;
		for (int k = 1; k < deg; k++) {
			double pn = ((2 * k + 1) * x * pk - k * pm) / (k + 1);
			pm = pk; pk = pn;
			X[(size_t)(k + 1) * n + t] = pn;
		}
		for (int j = 0; j < 6; j++) X[(size_t)(1 + deg + j) * n + t] = hmc[6 * (t + skip) + j];
	}
	int q = 1;   /* drop constant columns past the intercept (a still head has constant motion); dof uses the rest */
	for (int a = 1; a < p; a++) {
		const double *col = X + (size_t)a * n;
		int t = 1;
		while (t < n && col[t] == col[0]) t++;
		if (t == n) continue;
		if (q != a) memmove(X + (size_t)q * n, col, (size_t)n * sizeof(double));
		q++;
	}
	p = q;
	/* Gauss-Jordan inverse of X'X (p <= ~15), [G | I] -> [I | G^-1] */
	for (int a = 0; a < p; a++) {
		for (int b = 0; b < p; b++)
			for (int t = 0; t < n; t++) G[a * 2 * p + b] += X[(size_t)a * n + t] * X[(size_t)b * n + t];
		G[a * 2 * p + p + a] = 1.0;
	}
	for (int c = 0; c < p; c++) {
		int piv = c;
		for (int r = c + 1; r < p; r++) if (fabs(G[r * 2 * p + c]) > fabs(G[piv * 2 * p + c])) piv = r;
		if (fabs(G[piv * 2 * p + c]) < 1e-12) { printfx("-tsnr: singular design (collinear regressors)\n"); goto done; }
		for (int k = 0; k < 2 * p; k++) { double s = G[c * 2 * p + k]; G[c * 2 * p + k] = G[piv * 2 * p + k]; G[piv * 2 * p + k] = s; }
		double d = G[c * 2 * p + c];
		for (int k = 0; k < 2 * p; k++) G[c * 2 * p + k] /= d;
		for (int r = 0; r < p; r++) {
			if (r == c) continue;
			double f = G[r * 2 * p + c];
			for (int k = 0; k < 2 * p; k++) G[r * 2 * p + k] -= f * G[c * 2 * p + k];
		}
	}
	nv = 0;
	for (int64_t v = 0; v < n3; v++) if (brain[v]) idx[nv++] = v;
	int oom = 0;
#pragma omp parallel reduction(| : oom)
	{
	double *r = (double *)malloc((size_t)p * 2 * sizeof(double));
	if (!r) oom = 1;
#pragma omp for
	for (int64_t i = 0; i < nv; i++) {
		if (!r) continue;
		int64_t v = idx[i];
		double *beta = r + p, mean = 0.0, ss = 0.0;
		memset(r, 0, (size_t)p * 2 * sizeof(double));
		for (int t = 0; t < n; t++) {
			double y = Y[(size_t)(t + skip) * n3 + v];
			mean += y;
			for (int a = 0; a < p; a++) r[a] += X[(size_t)a * n + t] * y;
		}
		for (int a = 0; a < p; a++) for (int b = 0; b < p; b++) beta[a] += G[a * 2 * p + p + b] * r[b];
		for (int t = 0; t < n; t++) {
			double e = Y[(size_t)(t + skip) * n3 + v];
			for (int a = 0; a < p; a++) e -= beta[a] * X[(size_t)a * n + t];
			ss += e * e;
		}
		double sd = sqrt(ss / (n - p));
		map[v] = sd > 0.0 ? (float)(mean / n / sd) : 0.0f;
	}
	free(r);
	}
	if (oom) goto done;
	for (int64_t i = 0; i < nv; i++) vals[i] = map[idx[i]];
	double med = nv ? percentile(vals, nv, 50.0, 0) : 0.0;
	nifti_image saved = *nim;
	nim->data = map; nim->nt = nim->dim[4] = 1; nim->nvox = n3; nim->ndim = nim->dim[0] = 3;
	nim->scl_slope = 1.0f; nim->scl_inter = 0.0f;
	nim->fname = nim->iname = (char *)path;
	rc = nifti_save(nim, "", GZ_ENVIRONMENT);
	*nim = saved;
	if (!rc && json) {
		FILE *f = fopen(json, "w");
		rc = !f || fprintf(f, "{\n  \"tsnr_median_brain\": %.2f,\n  \"volumes\": %d,\n  \"polort\": %d,\n"
		                       "  \"regressors\": \"intercept, Legendre 1..%d, 6 motion parameters\",\n  \"regressors_used\": %d\n}\n",
		                   med, n, deg, deg, p) < 0;
		if (f && fclose(f)) rc = 1;
	}
	if (rc) printfx("-tsnr: failed writing '%s'\n", path);
done:
	free(X); free(G); free(map); free(vals); free(idx);
	return rc;
}

int nii_confounds(nifti_image *nim, const confounds_opts *o, const char *tsv) {
	int nt = (nim->ndim > 3 && nim->nt > 1) ? (int)nim->nt : 1;
	int64_t n3 = nim->nx * nim->ny * nim->nz;
	if (nim->datatype != DT_FLOAT32 || nt < 3 || !o->brain || !o->ref || !(o->tr > 0 && o->tr <= 60) ||
	    o->dummy < 0 || o->dummy >= nt - 2) {
		printfx("-confounds: needs a float32 4D series, -tr (seconds, <= 60), -hmc, -ref, -brain and a valid -dummy\n");
		return 1;
	}
	int rc = 1, skip = o->dummy;
	cc_table T = {NULL, 0, 0, nt};
	uint8_t *brain = NULL, *acc[3] = {NULL, NULL, NULL}, *tcc = NULL, *crown = NULL;
	double (*F)[4][4] = NULL;
	double *tmp = (double *)malloc((size_t)nt * sizeof(double)), *par = NULL;
	float *Y = resample_hmc(nim, o->hmc, o->vsm, o->bwarp);   /* motion (+SDC) corrected series on the boldref grid */
	brain = load_mask_on(o->brain, nim, NULL, 0.0);
	if (!Y || !brain || !tmp) goto done;
	double pull[16];
	if (o->acc[0]) {
		nifti_dmat44 f, fi;
		for (int i = 0; i < 16; i++) f.m[i / 4][i % 4] = o->ftm[i];
		fi = nifti_dmat44_inverse(f);   /* boldref world -> T1w world */
		for (int i = 0; i < 16; i++) pull[i] = fi.m[i / 4][i % 4];
		for (int k = 0; k < 3; k++) {
			if (!(acc[k] = load_mask_on(o->acc[k], nim, pull, ACC_THR))) goto done;
			for (int64_t v = 0; v < n3; v++) acc[k][v] &= brain[v];
		}
	}
	/* --- signals --- */
	static const char *sig_name[4] = {"global_signal", "csf", "white_matter", "csf_wm"};
	const uint8_t *sig_mask[4] = {brain, acc[0], acc[1], acc[2]};
	for (int k = 0; k < 4; k++) {
		if (!sig_mask[k]) continue;
		mean_signal(Y, n3, nt, sig_mask[k], tmp);
		if (k == 3) { double *c = col_add(&T, "%s", sig_name[k]); if (!c) goto done; memcpy(c, tmp, nt * sizeof(double)); }
		else if (col_expand(&T, sig_name[k], tmp)) goto done;
	}
	if (!(tcc = (uint8_t *)malloc((size_t)n3)) || tcompcor_mask(Y, n3, nt, skip, brain, tcc)) goto done;
	{
		double *c = col_add(&T, "tcompcor");
		if (!c) goto done;
		mean_signal(Y, n3, nt, tcc, c);
	}
	/* --- DVARS, FD, RMSD --- */
	double *sdv = col_add(&T, "std_dvars"), *dv = col_add(&T, "dvars");
	if (!sdv || !dv || dvars(Y, n3, nt, brain, dv, sdv)) goto done;
	F = (double (*)[4][4])malloc((size_t)nt * sizeof(*F));
	if (!F || fsl_matrices(nim, o->hmc, nt, F)) goto done;
	{
		/* centre of gravity of the boldref, in scaled voxels (ndi.center_of_mass x zooms) */
		nifti_image *ref = nii_read_float32(o->ref);
		if (!ref) goto done;
		if (nii_reject_5d(ref, "-confounds")) { nifti_image_free(ref); goto done; }
		double cog[3] = {0, 0, 0}, mass = 0.0;
		const float *r = (const float *)ref->data;
		for (int64_t v = 0; v < ref->nx * ref->ny * ref->nz; v++) {
			cog[0] += r[v] * (double)(v % ref->nx);
			cog[1] += r[v] * (double)((v / ref->nx) % ref->ny);
			cog[2] += r[v] * (double)(v / (ref->nx * ref->ny));
			mass += r[v];
		}
		double z[3] = {ref->dx, ref->dy, ref->dz};
		for (int a = 0; a < 3; a++) cog[a] = cog[a] / mass * z[a];
		nifti_image_free(ref);
		double *fd = col_add(&T, "framewise_displacement"), *rmsd = col_add(&T, "rmsd");
		if (!fd || !rmsd) goto done;
		if (!(par = (double *)malloc((size_t)nt * 6 * sizeof(double)))) goto done;
		for (int t = 0; t < nt; t++) {
			double R[3][3];
			for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) R[i][j] = F[t][j][i];   /* mats = F[:3,:3]^T */
			for (int i = 0; i < 3; i++) {
				double rc0 = 0.0;
				for (int j = 0; j < 3; j++) rc0 += R[i][j] * cog[j];
				par[t * 6 + i] = F[t][i][3] - rc0 + cog[i];
			}
			/* scipy Rotation.as_euler('XYZ'): R = Rx(a) Ry(b) Rz(c) */
			double sb = R[0][2] > 1 ? 1 : (R[0][2] < -1 ? -1 : R[0][2]);
			par[t * 6 + 3] = atan2(-R[1][2], R[2][2]);
			par[t * 6 + 4] = asin(sb);
			par[t * 6 + 5] = atan2(-R[0][1], R[0][0]);
		}
		double c[3] = {0.5 * (nim->nx - 1) * nim->dx, 0.5 * (nim->ny - 1) * nim->dy, 0.5 * (nim->nz - 1) * nim->dz};
		fd[0] = rmsd[0] = NAN;
		for (int t = 1; t < nt; t++) {
			double s = 0.0;
			for (int k = 0; k < 6; k++) s += fabs(par[t * 6 + k] - par[(t - 1) * 6 + k]) * (k < 3 ? 1.0 : FD_RADIUS);
			fd[t] = s;
			double Pi[4][4], D[4][4];
			mat_inv4(F[t - 1], Pi);
			mat_mul4(F[t], Pi, D);
			for (int i = 0; i < 4; i++) D[i][i] -= 1.0;
			double tt = 0.0, tr = 0.0;
			for (int i = 0; i < 3; i++) {
				double ti = D[i][3];
				for (int j = 0; j < 3; j++) { ti += D[i][j] * c[j]; tr += D[i][j] * D[i][j]; }
				tt += ti * ti;
			}
			rmsd[t] = sqrt(tt + tr * RMS_RADIUS * RMS_RADIUS / 5.0);
		}
		/* --- CompCor --- */
		if (compcor(&T, Y, n3, nt, skip, o->tr, tcc, 0, "t_comp_cor_")) goto done;
		static const char *acc_prefix[3] = {"c_comp_cor_", "w_comp_cor_", "a_comp_cor_"};
		for (int k = 0; k < 3 && acc[0]; k++)
			if (compcor(&T, Y, n3, nt, skip, o->tr, acc[k], 0, acc_prefix[k])) goto done;
		if (o->t1mask) {
			uint8_t *t1m = load_mask_on(o->t1mask, nim, pull, 0.5);
			crown = (uint8_t *)calloc((size_t)n3, 1);
			if (!t1m || !crown) { free(t1m); goto done; }
			for (int64_t v = 0; v < n3; v++) t1m[v] |= brain[v];   /* union */
			int64_t nx = nim->nx, ny = nim->ny, nz = nim->nz;
			for (int64_t v = 0; v < n3; v++) {
				if (t1m[v]) continue;
				int64_t x = v % nx, y = (v / nx) % ny, z = v / (nx * ny), hit = 0;
				for (int dz = -2; dz <= 2 && !hit; dz++)   /* skimage ball(2) dilation */
					for (int dy = -2; dy <= 2 && !hit; dy++)
						for (int dx = -2; dx <= 2 && !hit; dx++) {
							if (dx * dx + dy * dy + dz * dz > 4) continue;
							int64_t xx = x + dx, yy = y + dy, zz = z + dz;
							if (xx < 0 || yy < 0 || zz < 0 || xx >= nx || yy >= ny || zz >= nz) continue;
							hit = t1m[xx + nx * (yy + ny * zz)];
						}
				if (!hit) continue;
				double m = 0.0, ss = 0.0;   /* clip to voxels whose series has variance */
				for (int t = 0; t < nt; t++) m += Y[(size_t)t * n3 + v];
				m /= nt;
				for (int t = 0; t < nt; t++) ss += (Y[(size_t)t * n3 + v] - m) * (Y[(size_t)t * n3 + v] - m);
				crown[v] = ss > 0.0;
			}
			free(t1m);
			if (compcor(&T, Y, n3, nt, skip, o->tr, crown, CC_CROWN, "edge_comp_")) goto done;
		}
		/* --- cosine basis and non-steady-state indicators --- */
		int n = nt - skip, order = cosine_order(n, o->tr);
		for (int k = 1; k < order; k++) {
			double *c2 = col_add(&T, "cosine%02d", k - 1);
			if (!c2) goto done;
			for (int t = 0; t < n; t++) c2[t + skip] = cosine_basis(n, k, t);
		}
		for (int k = 0; k < skip; k++) {
			double *c2 = col_add(&T, "non_steady_state_outlier%02d", k);
			if (!c2) goto done;
			c2[k] = 1.0;
		}
		/* --- motion parameters, expanded --- */
		static const char *mp_name[6] = {"trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"};
		for (int k = 0; k < 6; k++) {
			for (int t = 0; t < nt; t++) tmp[t] = par[t * 6 + k];
			if (col_expand(&T, mp_name[k], tmp)) goto done;
		}
		/* --- spike regressors: FD > 0.5 mm or std DVARS > 1.5 --- */
		int ns = 0;
		for (int t = 0; t < nt; t++) {
			if (!((fd[t] > SPIKE_FD) || sdv[t] > SPIKE_DVARS)) continue;
			double *c2 = col_add(&T, "motion_outlier%02d", ns++);
			if (!c2) goto done;
			c2[t] = 1.0;   /* fd/sdv stay valid: each column owns its storage */
		}
	}
	rc = table_write(&T, tsv);
	if (!rc && o->tsnr) rc = write_tsnr(nim, Y, n3, nt, skip, o->tr, o->hmc, brain, o->tsnr, o->tsnr_json);
done:
	if (rc) printfx("-confounds: failed\n");
	for (int i = 0; i < T.n; i++) free(T.c[i].v);
	free(T.c); free(Y); free(tmp); free(par); free(F); free(brain); free(tcc); free(crown);
	for (int k = 0; k < 3; k++) free(acc[k]);
	return rc;
}

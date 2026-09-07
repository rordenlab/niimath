// refill.c - REFILL dynamic distortion correction (see refill.h)
//
// Port of the MATLAB reference (simon-mri/REFILL-Dynamic-Distortion-Correction, MIT):
//   rescale.m, refill_calc_readout_gradient.m, refill_calc_fms.m, refill_do_dc.m,
//   aspire_unwarp.m, smoothn.m/dctn.m/idctn.m (Garcia 2010, BSD).
// Every stage is validated against that code's own intermediates (test/refill_reference_manifest.md).
//
// Conventions that are load-bearing and invisible in the output:
//   * Everything works in STORED voxel order.  MATLAB dim 1 = NIfTI i (x, fastest), dim 2 = j (y),
//     dim 3 = k (z).  The readout ramp runs along x; the phase-encoding axis is y.
//   * rescale.m runs in double (load_untouch_nii with a volume index returns double), globally
//     over every volume it is handed: (2pi*(x-min))/range + (-pi).
//   * (uw - ramp)/TE is evaluated in SINGLE (double array minus single array is single in MATLAB);
//     the GE (uw2-uw1)/dTE difference is double.  Both are stored float32.
//   * smoothn is an early-stopped iteration, not a converged solution.  Its result is defined by
//     the initial guess (bwdist nearest-neighbour fill with MATLAB's tie rule, then a DCT low-pass
//     keeping the first ceil(n/10) coefficients per axis), RF = 1.75, TolZ = 1e-3, MaxIter = 100.
//     One extra iteration moves the result by ~1e-3 relative, so the iteration count must match.
//   * bwdist ties: the LOWEST column-major linear index among equidistant features (measured).
//     Three per-line brute-force passes (x, then y, then z), each scanning upward with a strict
//     '<', reproduce it exactly; the O(n^2) per line is ~250M ops per volume, not worth an
//     envelope algorithm that would need its own tie handling.
//   * aspire_unwarp is a FORWARD map: sample k sits at k + vsm(k) and the line is re-read on the
//     integer grid by interp1(...,'linear','extrap') after MATLAB sorts the (jittered) grid.  The
//     1e-5*rand jitter is omitted; ties keep stored order.  Interpolation runs in double.
//   * Non-finite tests use magnitude guards (the tree is -ffast-math -fno-finite-math-only).
// Diagnostics go through printfx (stderr): the output image may be going to stdout.
//
// ponytail: matrix (O(n^2)) DCT per axis, OpenMP over volumes.  An FFT-based DCT is the upgrade
// if the 100-volume smoothn stage ever dominates a pipeline; measured in the manifest.

#include "refill.h"
#include "romeo.h"
#include "print.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define RF_PI 3.14159265358979323846
#define RF_INF_I64 ((int64_t)1 << 60)

static int rf_finite_f(float v) { return v >= -FLT_MAX && v <= FLT_MAX; }

refill_opts refill_opts_default(void) {
	refill_opts o;
	memset(&o, 0, sizeof o);
	o.clip_lo = -600.0; o.clip_hi = 2000.0;   // refill_set_defaults: fmthreshl / fmthreshh (7T)
	o.s = 2.0;                                 // data.smoothns
	o.qthresh = 0.5;                           // data.qthresh_epi
	o.echo1 = 1; o.echo2 = 3;                  // data.echoes_to_use (bipolar FLASH)
	return o;
}

/* ================================ I/O helpers ================================ */

/* Read an auxiliary image on the working image's grid as float32 (scaled), with `want`
   volumes (0 = any).  The wrapper's nii_reject_oversize_aux already refused a huge header; the
   INT_MAX re-check after the load is the fail-closed guarantee (.zst has no bounded header read).
   Dimensions must match; the MATLAB reference works in voxel space and never consults a
   transform, so a differing world frame only warns (as -romeo treats its magnitude). */
static nifti_image *rf_read_aux(const nifti_image *nim, const char *fn, const char *op, const char *what, int64_t want) {
	nifti_image *n = nifti_image_read(fn, 1);
	in_hdr ihdr;
	int64_t n3 = (int64_t)nim->nx * nim->ny * nim->nz;
	if (!n) { printfx("%s: failed to read %s '%s'\n", op, what, fn); return NULL; }
	if (n->nu > 1 || n->nv > 1 || n->nw > 1 || (int64_t)n->nvox > INT_MAX) {
		printfx("%s: %s '%s' must be 3D or 4D with at most INT_MAX voxels\n", op, what, fn);
		nifti_image_free(n); return NULL;
	}
	if (n->nx != nim->nx || n->ny != nim->ny || n->nz != nim->nz || (want && n->nvox != want * n3)) {
		printfx("%s: %s '%s' must be %dx%dx%d with %s volume(s)\n", op, what, fn,
			(int)nim->nx, (int)nim->ny, (int)nim->nz, want ? (want == 1 ? "1" : "the input's") : "any number of");
		nifti_image_free(n); return NULL;
	}
	if (!(max_displacement_mm((nifti_image *)nim, n) <= 0.5f))
		printfx("%s: warning: %s '%s' has a different spatial transform (>0.5mm); voxel grids are used as stored\n", op, what, fn);
	ihdr = set_input_hdr(n);
	if (n->datatype != DT_FLOAT32 || (n->scl_slope != 0.0f && n->scl_slope != 1.0f) || n->scl_inter != 0.0f) {
		if (nifti_image_change_datatype(n, DT_FLOAT32, &ihdr) != 0) {
			printfx("%s: failed to convert %s '%s' to float32\n", op, what, fn);
			nifti_image_free(n); return NULL;
		}
	}
	return n;
}

/* Write a float32 companion on the working image's grid.  `dir`+`name` gives an explicit
   <dir>/<name>.nii[.gz] (the -steps dump); otherwise `postfix` is appended to the output name.
   nifti_save reads nim->data and recomputes dim[] from n*. */
static int rf_write(nifti_image *nim, const char *dir, const char *name, const char *postfix,
	const float *vals, int64_t n3, int nt, gzModes gzMode) {
	void *savedata = nim->data;
	char *savefname = nim->fname, *path = NULL;
	int saved_nt = nim->nt, rc;
	int64_t saved_nvox = nim->nvox;
	if (dir) {
		size_t len = strlen(dir) + strlen(name) + 8;
		path = (char *)malloc(len);
		if (!path) return 1;
		snprintf(path, len, "%s/%s.nii", dir, name);
		nim->fname = path;
		postfix = "";
	}
	nim->data = (void *)vals;
	nim->nt = nt; nim->nvox = n3 * nt;
	rc = nifti_save(nim, postfix, gzMode);
	nim->data = savedata; nim->fname = savefname;
	nim->nt = saved_nt; nim->nvox = saved_nvox;
	nim->dim[4] = saved_nt; nim->ndim = nim->dim[0] = (saved_nt > 1) ? 4 : 3;   /* nifti_save rewrote dim[] */
	free(path);
	return rc;
}

static int rf_step(nifti_image *nim, const refill_opts *o, const char *name, const float *vals, int64_t n3, int nt, gzModes gzMode) {
	if (!o->steps) return 0;
	return rf_write(nim, o->steps, name, NULL, vals, n3, nt, gzMode);
}

/* A byte mask as a float32 companion: `postfix` on the output name and/or a -steps dump `name`. */
static int rf_write_u8(nifti_image *nim, const refill_opts *o, const char *name, const char *postfix,
	const uint8_t *m, int64_t n3, gzModes gzMode) {
	float *f = (float *)malloc((size_t)n3 * sizeof(float));
	int64_t i;
	int rc = 0;
	if (!f) return 1;
	for (i = 0; i < n3; i++) f[i] = (float)m[i];
	if (postfix) rc |= rf_write(nim, NULL, NULL, postfix, f, n3, 1, gzMode);
	if (name) rc |= rf_step(nim, o, name, f, n3, 1, gzMode);
	free(f);
	return rc;
}

/* ================================ rescale.m ================================ */

/* (new_range*(mat-old_min))/old_range + new_min over the WHOLE array, in double. */
static void rf_rescale(float *p, int64_t n) {
	double mn = DBL_MAX, mx = -DBL_MAX;
	int64_t i;
	for (i = 0; i < n; i++) { if (p[i] < mn) mn = p[i]; if (p[i] > mx) mx = p[i]; }
	{
		const double range = mx - mn;
		for (i = 0; i < n; i++) p[i] = (float)((2.0 * RF_PI * ((double)p[i] - mn)) / range + (-RF_PI));
	}
}

/* ============================ readout gradient ============================ */

/* refill_calc_readout_gradient, arc=2: half the phase slope along x between the first EPI volume
   and the readout-reversed REFILL volume.  `rg` = angle(exp(i(p2-p1))) (double, both volumes
   rescaled on their own); the complex sum runs in double (MATLAB's single-precision sum order is
   unreproducible; the printed gradient agrees to 6 decimals).  Returns the gradient in rad/voxel. */
static double rf_gradient(const float *m1, const double *rg, int nx, int ny, int nz) {
	double sr = 0.0, si = 0.0;
	int x, y, z;
	for (z = 0; z < nz; z++) for (y = 0; y < ny; y++) {
		const int64_t row = (int64_t)nx * (y + (int64_t)ny * z);
		for (x = 1; x < nx; x++) {
			double a = rg[row + x] - rg[row + x - 1];
			float m = m1[row + x];
			float re = m * (float)cos(a), im = m * (float)sin(a);   /* single .* double complex: MATLAB narrows the double, multiplies in single */
			if (!(rf_finite_f(re) && rf_finite_f(im))) continue;       /* diffMap(isnan)=0 */
			sr += re; si += im;
		}
	}
	return atan2(si, sr) / 2.0;   /* angle / gradient_rescale_factor (shift = 1) */
}

/* ================================ smoothn ================================ */

typedef struct {
	int n;
	double *c;   /* n*n orthonormal DCT-II matrix, c[k*n+i] = w_k cos(pi*(2i+1)k/(2n)) */
} rf_dct;

static int rf_dct_init(rf_dct *d, int n) {
	int k, i;
	d->n = n;
	d->c = (double *)malloc((size_t)n * n * sizeof(double));
	if (!d->c) return 1;
	for (k = 0; k < n; k++) {
		double w = sqrt((k == 0 ? 1.0 : 2.0) / n);
		for (i = 0; i < n; i++) d->c[(int64_t)k * n + i] = w * cos(RF_PI * (2.0 * i + 1.0) * k / (2.0 * n));
	}
	return 0;
}

/* One axis of the N-D DCT (forward) or its inverse (transpose), in place, via a line buffer. */
static void rf_dct_axis(double *a, const rf_dct *d, int axis, int nx, int ny, int nz, int inverse, double *tmp) {
	const int n = d->n;
	const int64_t step = (axis == 0) ? 1 : (axis == 1) ? nx : (int64_t)nx * ny;
	const int nb = (axis == 0) ? ny : nx;
	const int nc = (axis == 2) ? ny : nz;
	const int64_t bstep = (axis == 0) ? nx : 1;
	const int64_t cstep = (axis == 2) ? nx : (int64_t)nx * ny;
	int b, c, k, i;
	for (c = 0; c < nc; c++) for (b = 0; b < nb; b++) {
		double *line = a + b * bstep + c * cstep;
		for (i = 0; i < n; i++) tmp[i] = line[i * step];
		for (k = 0; k < n; k++) {
			double acc = 0.0;
			if (inverse) { for (i = 0; i < n; i++) acc += d->c[(int64_t)i * n + k] * tmp[i]; }
			else         { for (i = 0; i < n; i++) acc += d->c[(int64_t)k * n + i] * tmp[i]; }
			line[k * step] = acc;
		}
	}
}

/* bwdist-equivalent nearest-feature fill: z[missing] = y[nearest finite], MATLAB tie rule.
   Squared distances stay integer, so equality is exact.  Scratch: d, idx (int64, n3 each) and
   lb (int64, 2*max(ny,nz)) for one line's best distance/index. */
static void rf_nnfill(const double *y, const uint8_t *fin, double *z, int nx, int ny, int nz, int64_t *d, int64_t *idx, int64_t *lb) {
	const int64_t n3 = (int64_t)nx * ny * nz;
	int64_t i;
	int x, yy, zz, p;
	for (i = 0; i < n3; i++) { d[i] = fin[i] ? 0 : RF_INF_I64; idx[i] = i; }
	/* x pass: along each line, nearest finite x, lowest x on ties */
	for (zz = 0; zz < nz; zz++) for (yy = 0; yy < ny; yy++) {
		const int64_t row = (int64_t)nx * (yy + (int64_t)ny * zz);
		for (x = 0; x < nx; x++) {
			int64_t best = RF_INF_I64, bi = row + x;
			for (p = 0; p < nx; p++) {
				if (!fin[row + p]) continue;
				{ int64_t c = (int64_t)(x - p) * (x - p); if (c < best) { best = c; bi = row + p; } }
			}
			d[row + x] = best; idx[row + x] = bi;
		}
	}
	/* y pass: combine lines of the x pass, lowest y on ties (carrying that line's lowest x) */
	for (zz = 0; zz < nz; zz++) for (x = 0; x < nx; x++) {
		const int64_t col = x + (int64_t)nx * ny * zz;
		int64_t *bd = lb, *bidx = lb + ny;
		for (yy = 0; yy < ny; yy++) {
			int64_t best = RF_INF_I64, bi = col + (int64_t)nx * yy;
			for (p = 0; p < ny; p++) {
				int64_t dp = d[col + (int64_t)nx * p];
				if (dp >= RF_INF_I64) continue;
				{ int64_t c = dp + (int64_t)(yy - p) * (yy - p); if (c < best) { best = c; bi = idx[col + (int64_t)nx * p]; } }
			}
			bd[yy] = best; bidx[yy] = bi;
		}
		for (yy = 0; yy < ny; yy++) { d[col + (int64_t)nx * yy] = bd[yy]; idx[col + (int64_t)nx * yy] = bidx[yy]; }
	}
	/* z pass, lowest z on ties */
	for (yy = 0; yy < ny; yy++) for (x = 0; x < nx; x++) {
		const int64_t col = x + (int64_t)nx * yy, zstep = (int64_t)nx * ny;
		int64_t *bd = lb, *bidx = lb + nz;
		for (zz = 0; zz < nz; zz++) {
			int64_t best = RF_INF_I64, bi = col + zstep * zz;
			for (p = 0; p < nz; p++) {
				int64_t dp = d[col + zstep * p];
				if (dp >= RF_INF_I64) continue;
				{ int64_t c = dp + (int64_t)(zz - p) * (zz - p); if (c < best) { best = c; bi = idx[col + zstep * p]; } }
			}
			bd[zz] = best; bidx[zz] = bi;
		}
		for (zz = 0; zz < nz; zz++) { d[col + zstep * zz] = bd[zz]; idx[col + zstep * zz] = bidx[zz]; }
	}
	for (i = 0; i < n3; i++) z[i] = fin[i] ? y[i] : y[idx[i]];
}

typedef struct {
	int nx, ny, nz;
	rf_dct dx, dy, dz;
	double *gamma;   /* n3 */
} rf_smooth_ctx;

static void rf_dctn(const rf_smooth_ctx *c, double *a, int inverse, double *tmp) {
	rf_dct_axis(a, &c->dx, 0, c->nx, c->ny, c->nz, inverse, tmp);
	rf_dct_axis(a, &c->dy, 1, c->nx, c->ny, c->nz, inverse, tmp);
	rf_dct_axis(a, &c->dz, 2, c->nx, c->ny, c->nz, inverse, tmp);
}

static void rf_smooth_free(rf_smooth_ctx *c) { free(c->dx.c); free(c->dy.c); free(c->dz.c); free(c->gamma); memset(c, 0, sizeof *c); }

static int rf_smooth_init(rf_smooth_ctx *c, int nx, int ny, int nz, double s) {
	const int64_t n3 = (int64_t)nx * ny * nz;
	int x, y, z;
	memset(c, 0, sizeof *c);
	c->nx = nx; c->ny = ny; c->nz = nz;
	if (rf_dct_init(&c->dx, nx) || rf_dct_init(&c->dy, ny) || rf_dct_init(&c->dz, nz)) { rf_smooth_free(c); return 1; }
	c->gamma = (double *)malloc((size_t)n3 * sizeof(double));
	if (!c->gamma) { rf_smooth_free(c); return 1; }
	/* Lambda = -2*(d - sum_i cos(pi*(k_i)/n_i)) with d = ndims(y) = 3; Gamma = 1/(1+s*Lambda^2) */
	for (z = 0; z < nz; z++) for (y = 0; y < ny; y++) for (x = 0; x < nx; x++) {
		double lam = cos(RF_PI * x / nx) + cos(RF_PI * y / ny) + cos(RF_PI * z / nz);
		lam = -2.0 * (3.0 - lam);
		c->gamma[x + (int64_t)nx * (y + (int64_t)ny * z)] = 1.0 / (1.0 + s * lam * lam);
	}
	return 0;
}

/* smoothn(y, s) for one volume: y (double, n3) with `fin` marking finite entries; result in z.
   Scratch (caller-owned, per thread): dcty, z0 (double n3), idx, dist (int64 n3), tmp (double max n).
   Returns the iteration count, or -1 when no voxel is valid. */
static int rf_smoothn(const rf_smooth_ctx *c, double *y, const uint8_t *fin, double *z,
	double *dcty, double *z0, int64_t *dist, int64_t *idx, int64_t *lb, double *tmp) {
	const int nx = c->nx, ny = c->ny, nz = c->nz;
	const int64_t n3 = (int64_t)nx * ny * nz;
	int64_t i, nfin = 0;
	int isweighted, nit = 0;
	double tol = 1.0, RF;
	for (i = 0; i < n3; i++) nfin += fin[i] ? 1 : 0;
	/* No valid voxel at all: MATLAB's bwdist fill errors here; filling from the masked-OUT
	   values and smoothing them would publish garbage at exit 0. */
	if (nfin == 0) return -1;
	isweighted = (nfin < n3);
	if (isweighted) {
		/* InitialGuess: nearest-neighbour fill, then keep the first ceil(n/10) DCT coefficients per axis */
		int x, yy, zz;
		const int kx = (nx + 9) / 10, ky = (ny + 9) / 10, kz = (nz + 9) / 10;
		rf_nnfill(y, fin, z, nx, ny, nz, dist, idx, lb);
		rf_dctn(c, z, 0, tmp);
		for (zz = 0; zz < nz; zz++) for (yy = 0; yy < ny; yy++) for (x = 0; x < nx; x++)
			if (x >= kx || yy >= ky || zz >= kz) z[x + (int64_t)nx * (yy + (int64_t)ny * zz)] = 0.0;
		rf_dctn(c, z, 1, tmp);
	} else {
		for (i = 0; i < n3; i++) z[i] = 0.0;
	}
	memcpy(z0, z, (size_t)n3 * sizeof(double));
	for (i = 0; i < n3; i++) if (!fin[i]) y[i] = 0.0;
	RF = 1.0 + 0.75 * isweighted;
	while (tol > 1e-3 && nit < 100) {
		double num = 0.0, den = 0.0;
		nit++;
		/* DCTy = dctn(Wtot.*(y-z)+z), W = fin */
		for (i = 0; i < n3; i++) dcty[i] = fin[i] ? y[i] : z[i];
		rf_dctn(c, dcty, 0, tmp);
		for (i = 0; i < n3; i++) dcty[i] *= c->gamma[i];
		rf_dctn(c, dcty, 1, tmp);
		for (i = 0; i < n3; i++) z[i] = RF * dcty[i] + (1.0 - RF) * z[i];
		/* tol = isweighted*norm(z0(:)-z(:))/norm(z(:)) */
		for (i = 0; i < n3; i++) { double dd = z0[i] - z[i]; num += dd * dd; den += z[i] * z[i]; }
		tol = isweighted * sqrt(num) / sqrt(den);
		memcpy(z0, z, (size_t)n3 * sizeof(double));
	}
	return nit;
}

/* Mask + clip + smoothn over `nt` float32 volumes in place.  mask3d (n3 bytes, nonzero = keep)
   broadcasts over volumes.  Prints the per-volume iteration counts when `steps` dumps are on.
   OpenMP over volumes; every volume is independent so the result is thread-count invariant. */
static int rf_mask_smooth(float *fm, int64_t n3, int nt, const uint8_t *mask3d, int nx, int ny, int nz,
	const refill_opts *o, const char *op) {
	rf_smooth_ctx c;
	int nmax = nx > ny ? nx : ny;
	int t, oom = 0, bad = 0;
	if (nz > nmax) nmax = nz;
	if (o->s != 0.0 && rf_smooth_init(&c, nx, ny, nz, o->s)) { printfx("%s: out of memory\n", op); return 1; }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt < omp_get_max_threads() ? nt : omp_get_max_threads())
#endif
	{
		/* -s 0 needs only `fin`; the six full-volume smoothn buffers (~31 MiB per worker on
		   128x128x40) are allocated only when smoothing runs. */
		const int sm = (o->s != 0.0);
		double *y = sm ? (double *)malloc((size_t)n3 * sizeof(double)) : NULL;
		double *z = sm ? (double *)malloc((size_t)n3 * sizeof(double)) : NULL;
		double *dcty = sm ? (double *)malloc((size_t)n3 * sizeof(double)) : NULL;
		double *z0 = sm ? (double *)malloc((size_t)n3 * sizeof(double)) : NULL;
		int64_t *dist = sm ? (int64_t *)malloc((size_t)n3 * sizeof(int64_t)) : NULL;
		int64_t *idx = sm ? (int64_t *)malloc((size_t)n3 * sizeof(int64_t)) : NULL;
		double *tmp = sm ? (double *)malloc((size_t)nmax * sizeof(double)) : NULL;
		int64_t *lb = sm ? (int64_t *)malloc((size_t)2 * nmax * sizeof(int64_t)) : NULL;
		uint8_t *fin = (uint8_t *)malloc((size_t)n3);
		int myoom = !fin || (sm && !(y && z && dcty && z0 && dist && idx && tmp && lb));
		if (myoom) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
			oom = 1;
		}
#ifdef _OPENMP
#pragma omp barrier
#endif
		if (!oom) {
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
			for (t = 0; t < nt; t++) {
				float *v = fm + (int64_t)t * n3;
				int64_t i;
				int nit = 0;
				for (i = 0; i < n3; i++) {
					const float f = v[i];
					/* one_tp(one_tp_copy<lo | one_tp_copy>hi | qmask==0) = NaN; a non-finite input is missing too */
					fin[i] = (rf_finite_f(f) && !(f < o->clip_lo) && !(f > o->clip_hi) && mask3d[i]) ? 1 : 0;
					if (sm) y[i] = f;
				}
				if (o->s != 0.0) {
					nit = rf_smoothn(&c, y, fin, z, dcty, z0, dist, idx, lb, tmp);
					if (nit < 0) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
						bad = 1;
						continue;
					}
					for (i = 0; i < n3; i++) v[i] = (float)z[i];
				} else {
					for (i = 0; i < n3; i++) v[i] = fin[i] ? v[i] : NAN;
				}
				if (o->steps) {
#ifdef _OPENMP
#pragma omp critical
#endif
					printfx("%s: smoothn volume %d: %d iterations\n", op, t + 1, nit);
				}
			}
		}
		free(y); free(z); free(dcty); free(z0); free(dist); free(idx); free(lb); free(tmp); free(fin);
	}
	if (o->s != 0.0) rf_smooth_free(&c);
	if (oom) { printfx("%s: out of memory\n", op); return 1; }
	if (bad) { printfx("%s: a volume has no valid voxel after masking and clipping; refusing to write\n", op); return 1; }
	return 0;
}

/* ================================ ROMEO ================================ */

/* Unwrap `neco` echoes of `phase` (radians, in place) with magnitude `mag`, the way
   refill_calc_fms calls the ROMEO CLI: template 1, -g, robustmask from the magnitude.
   individual = `-i` (EPI).  qmap/mask are optional n3 outputs. */
static int rf_romeo(const char *op, float *phase, const float *mag, int nx, int ny, int nz, int neco,
	const double *TEs, int individual, float *qmap, uint8_t *mask) {
	romeo_opts o = romeo_opts_default();
	o.template_echo = 1;
	o.correctglobal = 1;
	o.individual = individual;
	o.mask_sel = RM_MASK_ROBUST;
	if (romeo_unwrap_frame(phase, mag, neco, nx, ny, nz, neco, TEs, &o, NULL, mask, qmap)) { printfx("%s: ROMEO unwrapping failed\n", op); return 1; }
	return 0;
}

static int rf_common_checks(const nifti_image *nim, const char *op, int *nx, int *ny, int *nz, int *nt, int64_t *n3) {
	if (nim->nu > 1 || nim->nv > 1 || nim->nw > 1 || (int64_t)nim->nvox > INT_MAX) { printfx("%s: the input must be 3D or 4D with at most INT_MAX voxels\n", op); return 1; }
	*nx = nim->nx; *ny = nim->ny; *nz = nim->nz; *n3 = (int64_t)*nx * *ny * *nz;
	*nt = (int)(nim->nvox / *n3);
	if (*ny < 2 || *nx < 2) { printfx("%s: x and y need at least 2 voxels\n", op); return 1; }
	return 0;
}

/* ================================ -refill-gefm ================================ */

int refill_gefm(nifti_image *nim, const char *phasefile, const char *maskfile,
	double te1_ms, double te2_ms, const refill_opts *o, gzModes gzMode) {
	const char *op = "-refill-gefm";
	nifti_image *ph = NULL, *mk = NULL;
	float *p2 = NULL, *m2 = NULL, *fm = NULL;
	uint8_t *mask = NULL;
	int64_t n3, i;
	int nx, ny, nz, ne, rc = 1;
	if (rf_common_checks(nim, op, &nx, &ny, &nz, &ne, &n3)) return 1;
	if (o->echo1 < 1 || o->echo2 < 1 || o->echo1 > ne || o->echo2 > ne || o->echo1 == o->echo2) {
		printfx("%s: echoes %d,%d are not two distinct volumes of a %d-echo input\n", op, o->echo1, o->echo2, ne); return 1;
	}
	ph = rf_read_aux(nim, phasefile, op, "phase", ne);
	mk = ph ? rf_read_aux(nim, maskfile, op, "mask", 1) : NULL;
	if (!mk) goto done;

	/* rescale over ALL echoes (refill_copy_rescale hands rescale.m the whole 4D array) */
	rf_rescale((float *)ph->data, ph->nvox);
	if (rf_step(nim, o, "ge_p", (const float *)ph->data, n3, ne, gzMode)) goto done;

	/* select the two echoes (ROMEO -e [e1,e2] subsets phase AND magnitude) */
	p2 = (float *)malloc((size_t)2 * n3 * sizeof(float));
	m2 = (float *)malloc((size_t)2 * n3 * sizeof(float));
	mask = (uint8_t *)malloc((size_t)n3);
	if (!p2 || !m2 || !mask) { printfx("%s: out of memory\n", op); goto done; }
	memcpy(p2, (float *)ph->data + (int64_t)(o->echo1 - 1) * n3, (size_t)n3 * sizeof(float));
	memcpy(p2 + n3, (float *)ph->data + (int64_t)(o->echo2 - 1) * n3, (size_t)n3 * sizeof(float));
	memcpy(m2, (float *)nim->data + (int64_t)(o->echo1 - 1) * n3, (size_t)n3 * sizeof(float));
	memcpy(m2 + n3, (float *)nim->data + (int64_t)(o->echo2 - 1) * n3, (size_t)n3 * sizeof(float));
	{
		double TEs[2]; TEs[0] = te1_ms; TEs[1] = te2_ms;
		if (rf_romeo(op, p2, m2, nx, ny, nz, 2, TEs, 0, NULL, NULL)) goto done;
	}
	if (rf_step(nim, o, "ge_pd_uw", p2, n3, 2, gzMode)) goto done;

	/* fm = (uw2 - uw1)/dTE, double (both operands are double in MATLAB), stored float32 */
	fm = (float *)nii_malloc((size_t)n3, sizeof(float));
	{
		const double dte = (te2_ms - te1_ms) * 1e-3;
		for (i = 0; i < n3; i++) fm[i] = (float)(((double)p2[n3 + i] - (double)p2[i]) / dte);
	}
	if (rf_step(nim, o, "ge_fm", fm, n3, 1, gzMode)) goto done;
	for (i = 0; i < n3; i++) mask[i] = (((const float *)mk->data)[i] != 0.0f) ? 1 : 0;
	if (rf_mask_smooth(fm, n3, 1, mask, nx, ny, nz, o, op)) goto done;

	free(nim->data);
	nim->data = fm; fm = NULL;
	nim->nt = 1; nim->dim[4] = 1; nim->ndim = 3; nim->dim[0] = 3; nim->nvox = n3;
	rc = 0;
done:
	free(fm); free(p2); free(m2); free(mask);
	if (ph) nifti_image_free(ph);
	if (mk) nifti_image_free(mk);
	return rc;
}

/* ================================ -refill-epifm ================================ */

int refill_epifm(nifti_image *nim, const char *phasefile, const char *refillfile,
	double te_ms, const refill_opts *o, gzModes gzMode) {
	const char *op = "-refill-epifm";
	nifti_image *ph = NULL, *rp = NULL;
	float *ramp = NULL, *qmap = NULL, *fm = NULL;
	double *rg = NULL, *TEs = NULL;
	uint8_t *mask = NULL, *qmask = NULL;
	int64_t n3, i;
	int nx, ny, nz, nt, t, rc = 1;
	if (rf_common_checks(nim, op, &nx, &ny, &nz, &nt, &n3)) return 1;
	ph = rf_read_aux(nim, phasefile, op, "phase", nt);
	rp = ph ? rf_read_aux(nim, refillfile, op, "readout-reversed phase", 0) : NULL;   /* volume 1 is used */
	if (!rp) goto done;

	/* readout gradient from volume 1 of the EPI phase and the REFILL volume, each rescaled alone */
	ramp = (float *)malloc((size_t)n3 * sizeof(float));
	rg = (double *)malloc((size_t)n3 * sizeof(double));
	if (!ramp || !rg) { printfx("%s: out of memory\n", op); goto done; }
	{
		float *p1 = (float *)malloc((size_t)n3 * sizeof(float));
		float *pr = (float *)rp->data;
		double g;
		int x, y, z;
		if (!p1) { printfx("%s: out of memory\n", op); goto done; }
		memcpy(p1, ph->data, (size_t)n3 * sizeof(float));
		rf_rescale(p1, n3); rf_rescale(pr, n3);
		for (i = 0; i < n3; i++) { double a = (double)pr[i] - (double)p1[i]; rg[i] = atan2(sin(a), cos(a)); }   /* angle(exp(i*(p2-p1))) */
		free(p1);
		g = rf_gradient((const float *)nim->data, rg, nx, ny, nz);
		printfx("%s: there is a linear gradient in the readout direction of %2.6f Hz/voxel - removing\n", op, g / (2.0 * RF_PI));
		/* readoutGradient = rValues * gradient: double rValues, single gradient -> single.
		   REFERENCE SIGN, reproduced by default: the code forms angle(exp(i(p_REFILL - p_1))) and
		   then SUBTRACTS the ramp, whereas the paper's Eq. 7 defines phi_G from theta_1 - theta_REFILL.
		   The two differ by a sign, so the reference doubles the residual readout gradient instead
		   of removing it -- measured on sub1 of the validation set, where the x-slope of
		   (EPI - FLASH field map) is 0.78 rad/s/voxel with the reference sign and 0.009 with the
		   paper's.  -ramp-fix applies the paper's sign; see test/refill_reference_manifest.md. */
		if (o->ramp_fix) g = -g;
		for (z = 0; z < nz; z++) for (y = 0; y < ny; y++) for (x = 0; x < nx; x++)
			ramp[x + (int64_t)nx * (y + (int64_t)ny * z)] = (float)((double)(-nx / 2 + x) * (double)(float)g);
		if (o->steps) {
			float *b = (float *)malloc((size_t)n3 * sizeof(float));
			int brc;
			if (!b) { printfx("%s: out of memory\n", op); goto done; }
			for (i = 0; i < n3; i++) b[i] = (float)rg[i];
			brc = rf_step(nim, o, "basis_for_gradient_calc", b, n3, 1, gzMode) || rf_step(nim, o, "ramp", ramp, n3, 1, gzMode);
			free(b);
			if (brc) goto done;
		}
	}

	/* rescale over ALL volumes, unwrap individually (-t epi -i --template 1 -g -q) */
	rf_rescale((float *)ph->data, ph->nvox);
	if (rf_step(nim, o, "epi_p", (const float *)ph->data, n3, nt, gzMode)) goto done;
	TEs = (double *)malloc((size_t)nt * sizeof(double));
	qmap = (float *)malloc((size_t)n3 * sizeof(float));
	mask = (uint8_t *)malloc((size_t)n3);
	qmask = (uint8_t *)malloc((size_t)n3);
	if (!TEs || !qmap || !mask || !qmask) { printfx("%s: out of memory\n", op); goto done; }
	for (t = 0; t < nt; t++) TEs[t] = 1.0;   /* -t epi: ones(neco) */
	if (rf_romeo(op, (float *)ph->data, (const float *)nim->data, nx, ny, nz, nt, TEs, 1, qmap, mask)) goto done;
	if (rf_step(nim, o, "epi_phase_uw", (const float *)ph->data, n3, nt, gzMode)) goto done;

	/* fm = (p_uw - ramp)/TE: double minus single is SINGLE in MATLAB, so this is float arithmetic */
	fm = (float *)nii_malloc((size_t)n3 * (size_t)nt, sizeof(float));
	{
		const float te = (float)(te_ms * 1e-3);
		const float *uw = (const float *)ph->data;
		for (t = 0; t < nt; t++) for (i = 0; i < n3; i++) fm[(int64_t)t * n3 + i] = uw[(int64_t)t * n3 + i] - ramp[i];
		if (rf_step(nim, o, "epi_phase_uw_rc", fm, n3, nt, gzMode)) goto done;
		for (i = 0; i < (int64_t)nt * n3; i++) fm[i] = fm[i] / te;
	}
	if (rf_step(nim, o, "epi_fm", fm, n3, nt, gzMode)) goto done;
	for (i = 0; i < n3; i++) qmask[i] = (qmap[i] > (float)o->qthresh) ? 1 : 0;
	if (rf_write_u8(nim, o, "qmask", NULL, qmask, n3, gzMode)) goto done;
	if (rf_mask_smooth(fm, n3, nt, qmask, nx, ny, nz, o, op)) goto done;
	/* side outputs last, so a failure above leaves no orphaned companions */
	if (rf_step(nim, o, "epi_quality", qmap, n3, 1, gzMode) || rf_write_u8(nim, o, "epi_mask", "_mask", mask, n3, gzMode) ||
		rf_write(nim, NULL, NULL, "_quality", qmap, n3, 1, gzMode)) {
		printfx("%s: failed to write a side output\n", op); goto done;
	}

	free(nim->data);
	nim->data = fm; fm = NULL;
	rc = 0;
done:
	free(fm); free(ramp); free(rg); free(TEs); free(qmap); free(mask); free(qmask);
	if (ph) nifti_image_free(ph);
	if (rp) nifti_image_free(rp);
	return rc;
}

/* ============================ VSM + aspire_unwarp ============================ */

/* refill_do_dc: `if strcmp(epi_PE_dir,'y') vsm=-vsm`; "y-" keeps the sign.  The wrapper has
   already rejected any other direction. */
static int rf_pe_sign(const char *dir) { return (dir[1] == '-') ? 1 : -1; }

/* vsm = (fm/(2pi))*dwell*N_PE (rad/s -> voxels), non-finite fm -> 0, in float like MATLAB's
   single field map times double scalars. */
static void rf_vsm(const float *fm, float *vsm, int64_t n, double dwell, int npe, int sign) {
	const double k = dwell * npe / (2.0 * RF_PI) * sign;
	int64_t i;
	for (i = 0; i < n; i++) { float f = fm[i]; vsm[i] = rf_finite_f(f) ? (float)(f * k) : 0.0f; }
}

/* aspire_unwarp for one line along y (stride `step`): sample k (1-based k+1) sits at k+1+vsm[k];
   interp1(grid, v, 1:ny, 'linear', 'extrap') after a stable sort of the grid.  `pos`/`val`/`ord`
   are scratch of length ny. */
static void rf_unwarp_line(const float *vsm, const float *v, float *out, int ny, int64_t step,
	double *pos, double *val, int *ord) {
	/* vsm/v/out share stride `step`; out may alias v because every read of v is copied to val first */
	int k, j, i;
	for (k = 0; k < ny; k++) { pos[k] = (double)(k + 1) + (double)vsm[k * step]; ord[k] = k; }
	/* stable insertion sort by position (ny is a few hundred at most; comparator-free for WASM) */
	for (i = 1; i < ny; i++) {
		int o = ord[i]; double p = pos[o];
		for (j = i - 1; j >= 0 && pos[ord[j]] > p; j--) ord[j + 1] = ord[j];
		ord[j + 1] = o;
	}
	for (k = 0; k < ny; k++) val[k] = (double)v[ord[k] * step];
	j = 0;
	for (k = 0; k < ny; k++) {
		const double q = (double)(k + 1);
		double x0, x1, t;
		while (j + 1 < ny - 1 && pos[ord[j + 1]] <= q) j++;   /* segment [j, j+1] with pos[j] <= q, clamped to the ends for extrapolation */
		x0 = pos[ord[j]]; x1 = pos[ord[j + 1]];
		t = (x1 != x0) ? (q - x0) / (x1 - x0) : 0.0;
		out[k * step] = (float)(val[j] + t * (val[j + 1] - val[j]));
	}
}

/* Forward-unwarp `nt` volumes of `img` (in place) with the matching volumes of `vsm`
   (nt_vsm == nt, or 1 to broadcast).  OpenMP over (z,x) lines; each line is independent. */
static int rf_unwarp_volumes(float *img, const float *vsm, int nx, int ny, int nz, int nt, int nt_vsm, const char *op) {
	const int64_t n3 = (int64_t)nx * ny * nz;
	int t, oom = 0;
	for (t = 0; t < nt; t++) {
		float *v = img + (int64_t)t * n3;
		const float *s = vsm + (int64_t)(nt_vsm > 1 ? t : 0) * n3;
		int l;
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
			double *pos = (double *)malloc((size_t)ny * sizeof(double));
			double *val = (double *)malloc((size_t)ny * sizeof(double));
			int *ord = (int *)malloc((size_t)ny * sizeof(int));
			if (!(pos && val && ord)) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
				oom = 1;
			}
#ifdef _OPENMP
#pragma omp barrier
#endif
			if (!oom) {
#ifdef _OPENMP
#pragma omp for
#endif
				for (l = 0; l < nz * nx; l++) {
					const int z = l / nx, x = l % nx;
					const int64_t base = x + (int64_t)nx * ny * z;
					rf_unwarp_line(s + base, v + base, v + base, ny, nx, pos, val, ord);
				}
			}
			free(pos); free(val); free(ord);
		}
		if (oom) { printfx("%s: out of memory\n", op); return 1; }
	}
	return 0;
}

/* Median over the finite values of `n` doubles (destroys the buffer); MATLAB's even-count mean
   of the two middle values.  Comparator-free quickselect. */
static double rf_select(double *a, int64_t n, int64_t k) {
	int64_t lo = 0, hi = n - 1;
	while (lo < hi) {
		double piv = a[(lo + hi) / 2];
		int64_t i = lo, j = hi;
		while (i <= j) {
			while (a[i] < piv) i++;
			while (a[j] > piv) j--;
			if (i <= j) { double t = a[i]; a[i] = a[j]; a[j] = t; i++; j--; }
		}
		if (k <= j) hi = j; else if (k >= i) lo = i; else return a[k];
	}
	return a[k];
}
static double rf_median(double *a, int64_t n) {
	double m = rf_select(a, n, n / 2);
	if (n % 2 == 0) { double m2 = rf_select(a, n, n / 2 - 1); return (m + m2) / 2.0; }
	return m;
}

/* ================================ -refill-centre ================================ */

int refill_centre(nifti_image *nim, const char *maskfile, double dwell, const char *dir,
	const refill_opts *o, gzModes gzMode) {
	const char *op = "-refill-centre";
	nifti_image *mk = NULL;
	float *vsm = NULL, *uw = NULL;
	double *vals = NULL;
	int64_t n3, i, nv = 0;
	int nx, ny, nz, nt, sign, t, rc = 1;
	double med;
	if (rf_common_checks(nim, op, &nx, &ny, &nz, &nt, &n3)) return 1;
	sign = rf_pe_sign(dir);
	mk = rf_read_aux(nim, maskfile, op, "mask", 1);
	if (!mk) goto done;
	vsm = (float *)malloc((size_t)nim->nvox * sizeof(float));
	uw = (float *)malloc((size_t)nim->nvox * sizeof(float));
	if (!vsm || !uw) { printfx("%s: out of memory\n", op); goto done; }
	/* preliminary VSM and unwarped field maps, only to find the in-mask median */
	rf_vsm((const float *)nim->data, vsm, nim->nvox, dwell, ny, sign);
	if (rf_step(nim, o, "vsm_prelim", vsm, n3, nt, gzMode)) goto done;
	for (i = 0; i < (int64_t)nim->nvox; i++) { float f = ((const float *)nim->data)[i]; uw[i] = rf_finite_f(f) ? f : 0.0f; }
	if (rf_unwarp_volumes(uw, vsm, nx, ny, nz, nt, nt, op)) goto done;
	/* REFERENCE QUIRK, reproduced deliberately: refill_do_dc.m does `epi_fm(mask~=1)=NaN` with a
	   3D logical mask on the 4D unwarped field maps.  MATLAB applies a logical index that is
	   smaller than the array to the FIRST numel(mask) linear elements, i.e. to VOLUME 1 ONLY, so
	   `median(vector(epi_fm),'omitnan')` pools the masked first volume with EVERY voxel of volumes
	   2..N (measured: sub1 gives -156.71 rad/s this way and -145.39 with all volumes masked).
	   The intent is clearly an in-mask median over all volumes, but the goal here is equivalence
	   with the reference as it runs; see test/refill_reference_manifest.md. */
	vals = (double *)malloc((size_t)nim->nvox * sizeof(double));
	if (!vals) { printfx("%s: out of memory\n", op); goto done; }
	for (t = 0; t < nt; t++) for (i = 0; i < n3; i++) {
		if (t == 0 && ((const float *)mk->data)[i] != 1.0f) continue;
		{ float f = uw[(int64_t)t * n3 + i]; if (rf_finite_f(f)) vals[nv++] = f; }
	}
	if (nv < 1) { printfx("%s: the mask selects no finite voxel\n", op); goto done; }
	med = rf_median(vals, nv);
	printfx("%s: Median in centre of EPI-FM=%4.2f radss-1; removing ... (median=%.9g)\n", op, med, med);
	/* current_fm = current_fm - median (single minus double -> single); NaN -> 0 */
	for (i = 0; i < (int64_t)nim->nvox; i++) {
		float f = ((float *)nim->data)[i];
		((float *)nim->data)[i] = rf_finite_f(f) ? (f - (float)med) : 0.0f;
	}
	rc = 0;
done:
	free(vsm); free(uw); free(vals);
	if (mk) nifti_image_free(mk);
	return rc;
}

/* ================================ -refill-unwarp ================================ */

int refill_unwarp(nifti_image *nim, const char *fmfile, double dwell, const char *dir,
	const refill_opts *o, gzModes gzMode) {
	const char *op = "-refill-unwarp";
	nifti_image *fm = NULL;
	float *vsm = NULL;
	int64_t n3;
	int nx, ny, nz, nt, ntf, sign, rc = 1;
	if (rf_common_checks(nim, op, &nx, &ny, &nz, &nt, &n3)) return 1;
	sign = rf_pe_sign(dir);
	/* a 4D input needs 1 or nt field-map volumes; a 3D input takes volume 1 of any count
	   (aspire_unwarp loops over the IMAGE's volumes) */
	fm = rf_read_aux(nim, fmfile, op, "field map", 0);
	if (!fm) goto done;
	ntf = (int)(fm->nvox / n3);
	if (nt > 1 && ntf != 1 && ntf != nt) {
		printfx("%s: field map '%s' has %d volumes but the input has %d (need 1, or the same count)\n", op, fmfile, ntf, nt); goto done;
	}
	vsm = (float *)malloc((size_t)fm->nvox * sizeof(float));
	if (!vsm) { printfx("%s: out of memory\n", op); goto done; }
	rf_vsm((const float *)fm->data, vsm, fm->nvox, dwell, ny, sign);
	if (rf_step(nim, o, "vsm", vsm, n3, ntf, gzMode)) goto done;
	/* a 3D working image takes field-map volume 1 (aspire_unwarp loops over the IMAGE's volumes) */
	if (rf_unwarp_volumes((float *)nim->data, vsm, nx, ny, nz, nt, (nt == 1) ? 1 : ntf, op)) goto done;
	rc = 0;
done:
	free(vsm);
	if (fm) nifti_image_free(fm);
	return rc;
}

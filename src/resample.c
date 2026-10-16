// resample.c - one-shot resampling (-resample); see resample.h.
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "resample.h"
#include "moco.h"
#include "print.h"

/* voxel -> world in mm, whatever the header's xyz_units: displacements, matrices and -res are mm */
static nifti_dmat44 best_affine(const nifti_image *n) {
	nifti_dmat44 m = (n->sform_code > 0) ? n->sto_xyz : n->qto_xyz;
	double f = xyz_units_to_mm(n->xyz_units);
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 4; c++) m.m[r][c] *= f;
	return m;
}

/* n * k floats, NULL if the byte count overflows size_t (wasm32) */
static float *alloc_f32(int64_t n, size_t k) {
	size_t b;
	return nii_mul_size((size_t)n, k * sizeof(float), &b) ? NULL : (float *)malloc(b);
}

static void apply44(const nifti_dmat44 *m, const double in[3], double out[3]) {
	for (int r = 0; r < 3; r++)
		out[r] = m->m[r][0] * in[0] + m->m[r][1] * in[1] + m->m[r][2] * in[2] + m->m[r][3];
}

static void apply44d(const double *m, const double in[3], double out[3]) {
	for (int r = 0; r < 3; r++)
		out[r] = m[r * 4] * in[0] + m[r * 4 + 1] * in[1] + m[r * 4 + 2] * in[2] + m[r * 4 + 3];
}

/* Trilinear sample of greedy's 3-component field (component-major 5D NIfTI), 0 outside. */
static void warp_vector(const nifti_image *w, const nifti_dmat44 *ijk, const double ras[3], double v[3]) {
	double p[3];
	apply44(ijk, ras, p);
	int64_t nx = w->nx, ny = w->ny, nz = w->nz, n3 = nx * ny * nz;
	const float *f = (const float *)w->data;
	v[0] = v[1] = v[2] = 0.0;
	if (!(p[0] > -1.0 && p[0] < nx && p[1] > -1.0 && p[1] < ny && p[2] > -1.0 && p[2] < nz)) return;   /* also NaN */
	int x0 = (int)floor(p[0]), y0 = (int)floor(p[1]), z0 = (int)floor(p[2]);
	double fx = p[0] - x0, fy = p[1] - y0, fz = p[2] - z0;
	for (int dz = 0; dz < 2; dz++)
		for (int dy = 0; dy < 2; dy++)
			for (int dx = 0; dx < 2; dx++) {
				int64_t x = x0 + dx, y = y0 + dy, z = z0 + dz;
				if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) continue;
				double wgt = (dx ? fx : 1 - fx) * (dy ? fy : 1 - fy) * (dz ? fz : 1 - fz);
				int64_t i = x + nx * (y + ny * z);
				for (int c = 0; c < 3; c++) v[c] += wgt * f[i + c * n3];
			}
}

/* x (RAS mm) moved by greedy's LPS displacement at x */
static void warp_point(const nifti_image *w, const nifti_dmat44 *ijk, double x[3]) {
	double u[3];
	warp_vector(w, ijk, x, u);
	x[0] -= u[0]; x[1] -= u[1]; x[2] += u[2];
}

/* In-place cubic B-spline prefilter of n samples at stride s (Unser; mirror boundary as scipy). */
static void bspline_line(float *c, int64_t n, int64_t s) {
	if (n < 2) return;
	const double z = sqrt(3.0) - 2.0, lambda = 6.0;
	for (int64_t i = 0; i < n; i++) c[i * s] *= lambda;
	int64_t horizon = (int64_t)ceil(log(1e-7) / log(fabs(z)));
	double sum = c[0];
	if (horizon < n) {   /* causal init: truncated sum of the mirror-extended signal */
		double zk = z;
		for (int64_t k = 1; k < horizon; k++) { sum += zk * c[k * s]; zk *= z; }
	} else {   /* short line: the exact mirror sum, whose tail the truncation would drop (Unser 1999) */
		double zn = z, z2n = pow(z, (double)(n - 1));
		sum += z2n * c[(n - 1) * s];
		z2n *= z2n / z;
		for (int64_t k = 1; k < n - 1; k++) { sum += (zn + z2n) * c[k * s]; zn *= z; z2n /= z; }
		sum /= 1.0 - zn * zn;
	}
	c[0] = (float)sum;
	for (int64_t i = 1; i < n; i++) c[i * s] += (float)(z * c[(i - 1) * s]);
	c[(n - 1) * s] = (float)((z / (z * z - 1.0)) * (c[(n - 1) * s] + z * c[(n - 2) * s]));
	for (int64_t i = n - 2; i >= 0; i--) c[i * s] = (float)(z * (c[(i + 1) * s] - c[i * s]));
}

static void bspline_prefilter(float *v, int64_t nx, int64_t ny, int64_t nz) {
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int64_t yz = 0; yz < ny * nz; yz++) bspline_line(v + yz * nx, nx, 1);
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int64_t xz = 0; xz < nx * nz; xz++) bspline_line(v + (xz / nx) * nx * ny + xz % nx, ny, nx);
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int64_t xy = 0; xy < nx * ny; xy++) bspline_line(v + xy, nz, nx * ny);
}

static inline int64_t mirror(int64_t i, int64_t n) {
	if (n == 1) return 0;
	int64_t p = 2 * (n - 1);
	i %= p;
	if (i < 0) i += p;
	return i < n ? i : p - i;
}

static inline void cubic_w(double t, double w[4]) {
	double t2 = t * t, t3 = t2 * t, u = 1.0 - t;
	w[0] = u * u * u / 6.0;
	w[1] = (3.0 * t3 - 6.0 * t2 + 4.0) / 6.0;
	w[2] = (-3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0) / 6.0;
	w[3] = t3 / 6.0;
}

/* Evaluate prefiltered coefficients c at voxel p; 0 outside [0, n-1] (scipy mode='constant'). */
static float bspline_eval(const float *c, const int64_t n[3], const double p[3]) {
	const double tol = 1e-3;
	for (int a = 0; a < 3; a++)
		if (!(p[a] >= -tol && p[a] <= (double)(n[a] - 1) + tol)) return 0.0f;   /* also NaN */
	int64_t i0[3];
	double w[3][4];
	for (int a = 0; a < 3; a++) {
		double f = floor(p[a]);
		i0[a] = (int64_t)f - 1;
		cubic_w(p[a] - f, w[a]);
	}
	double sum = 0.0;
	for (int k = 0; k < 4; k++) {
		int64_t z = mirror(i0[2] + k, n[2]);
		for (int j = 0; j < 4; j++) {
			int64_t y = mirror(i0[1] + j, n[1]);
			double wzy = w[2][k] * w[1][j];
			const float *row = c + n[0] * (y + n[1] * z);
			for (int i = 0; i < 4; i++) sum += wzy * w[0][i] * row[mirror(i0[0] + i, n[0])];
		}
	}
	return (float)sum;
}

/* Trilinear sample of volume `vol` of a float map at voxel p (edge-clamped). */
static double trilin(const nifti_image *m, int vol, const double p[3]) {
	int64_t n[3] = {m->nx, m->ny, m->nz}, i0[3];
	double f[3];
	for (int a = 0; a < 3; a++) {
		double q = p[a] > 0 ? (p[a] < n[a] - 1 ? p[a] : n[a] - 1) : 0;   /* NaN -> 0 */
		i0[a] = (int64_t)floor(q);
		if (i0[a] >= n[a] - 1) i0[a] = n[a] > 1 ? n[a] - 2 : 0;
		f[a] = n[a] > 1 ? q - i0[a] : 0.0;
	}
	const float *d = (const float *)m->data + (size_t)vol * n[0] * n[1] * n[2];
	double s = 0.0;
	for (int c = 0; c < 8; c++) {
		int64_t x = i0[0] + (c & 1), y = i0[1] + ((c >> 1) & 1), z = i0[2] + (c >> 2);
		if (x >= n[0]) x = n[0] - 1;
		if (y >= n[1]) y = n[1] - 1;
		if (z >= n[2]) z = n[2] - 1;
		double w = ((c & 1) ? f[0] : 1 - f[0]) * (((c >> 1) & 1) ? f[1] : 1 - f[1]) * ((c >> 2) ? f[2] : 1 - f[2]);
		s += w * d[x + n[0] * (y + n[1] * z)];
	}
	return s;
}

/* Sample nt volumes of src (float32) at per-voxel src coordinates `coords` (onv x 3), optionally
   through per-volume index-space pulls hmcM (12 doubles, 3x4 row-major) and a voxel-shift map
   (fMRIPrep resample_vol order: motion first, then the shift sampled at the reference position),
   into dst (nt x onv). */
static int pull_series(const nifti_image *src, int nt, const float *coords, int64_t onv,
                       const double (*hmcM)[12], const resample_vsm *vsm, float *dst) {
	int64_t in[3] = {src->nx, src->ny, src->nz}, inv = in[0] * in[1] * in[2];
	float *coef = (float *)malloc((size_t)inv * sizeof(float));
	if (!coef) return 1;
	int per_frame = vsm && vsm->map->nvox > inv;   /* 4D map: one volume per frame */
	for (int t = 0; t < nt; t++) {
		memcpy(coef, (const float *)src->data + (size_t)t * inv, (size_t)inv * sizeof(float));
		bspline_prefilter(coef, in[0], in[1], in[2]);
		const double *M = hmcM ? hmcM[t] : NULL;
		float *o3 = dst + (size_t)t * onv;
		int vt = per_frame ? t : 0;
#ifdef _OPENMP
		#pragma omp parallel for
#endif
		for (int64_t v = 0; v < onv; v++) {
			double p[3] = {coords[v * 3], coords[v * 3 + 1], coords[v * 3 + 2]}, q[3], jac = 1.0;
			if (M) apply44d(M, p, q); else memcpy(q, p, sizeof(q));   /* 3x4 row-major */
			if (vsm) {
				/* a static map is aligned to the reference (sample at p); a per-frame map (MEDIC) lives in
				   that frame's own native position, i.e. where motion has already pulled the voxel (q) */
				double at[3];
				memcpy(at, per_frame ? q : p, sizeof(at));
				q[vsm->axis] += trilin(vsm->map, vt, at);
				if (vsm->jacobian) {   /* 1 + d(vsm)/d(axis), central difference */
					double a[3], b[3];
					memcpy(a, at, sizeof(a)); memcpy(b, at, sizeof(b));
					a[vsm->axis] += 0.5; b[vsm->axis] -= 0.5;
					jac = 1.0 + trilin(vsm->map, vt, a) - trilin(vsm->map, vt, b);
				}
			}
			o3[v] = (float)(jac * bspline_eval(coef, in, q));
		}
	}
	free(coef);
	return 0;
}

static double (*hmc_pulls(const nifti_image *nim, const double *hmc, int nt))[12] {
	double (*M)[12] = (double (*)[12])malloc((size_t)nt * sizeof(*M));
	for (int t = 0; M && t < nt; t++)
		if (moco_pull_affine(nim, hmc + 6 * t, (double (*)[4])M[t])) { free(M); return NULL; }
	return M;
}

float *resample_hmc(const nifti_image *nim, const double *hmc, const resample_vsm *vsm, const nifti_image *bwarp) {
	int64_t n3 = nim->nx * nim->ny * nim->nz;
	int nt = (nim->ndim > 3 && nim->nt > 1) ? (int)nim->nt : 1;
	float *coords = alloc_f32(n3, 3);
	float *dst = alloc_f32(n3, nt);
	double (*M)[12] = hmc ? hmc_pulls(nim, hmc, nt) : NULL;
	int bad = !coords || !dst || (hmc && !M);
	for (int64_t v = 0; !bad && v < n3; v++) {
		coords[v * 3] = (float)(v % nim->nx);
		coords[v * 3 + 1] = (float)((v / nim->nx) % nim->ny);
		coords[v * 3 + 2] = (float)(v / (nim->nx * nim->ny));
	}
	if (!bad && bwarp) {   /* own grid -> world -> warped world -> own grid */
		nifti_dmat44 A = best_affine(nim), Ai = nifti_dmat44_inverse(A), W = nifti_dmat44_inverse(best_affine(bwarp));
#ifdef _OPENMP
		#pragma omp parallel for
#endif
		for (int64_t v = 0; v < n3; v++) {
			double p[3] = {coords[v * 3], coords[v * 3 + 1], coords[v * 3 + 2]}, x[3];
			apply44(&A, p, x);
			warp_point(bwarp, &W, x);
			apply44(&Ai, x, p);
			for (int a = 0; a < 3; a++) coords[v * 3 + a] = (float)p[a];
		}
	}
	if (!bad) bad = pull_series(nim, nt, coords, n3, (const double (*)[12])M, vsm, dst);
	free(coords); free(M);
	if (bad) { free(dst); return NULL; }
	return dst;
}

float *resample_volume(const nifti_image *src, const nifti_image *grid, const double *pull) {
	nifti_dmat44 S = best_affine(grid), ijk = nifti_dmat44_inverse(best_affine(src));
	int64_t n3 = grid->nx * grid->ny * grid->nz;
	float *coords = alloc_f32(n3, 3);
	float *dst = alloc_f32(n3, 1);
	if (!coords || !dst) { free(coords); free(dst); return NULL; }
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int64_t v = 0; v < n3; v++) {
		double i3[3] = {(double)(v % grid->nx), (double)((v / grid->nx) % grid->ny), (double)(v / (grid->nx * grid->ny))};
		double x[3], y[3];
		apply44(&S, i3, x);
		if (pull) { apply44d(pull, x, y); memcpy(x, y, sizeof(x)); }
		apply44(&ijk, x, y);
		for (int a = 0; a < 3; a++) coords[v * 3 + a] = (float)y[a];
	}
	int bad = pull_series(src, 1, coords, n3, NULL, NULL, dst);
	free(coords);
	if (bad) { free(dst); return NULL; }
	return dst;
}

nifti_image *nii_read_float32(const char *fn) {
	nifti_image *n = nifti_image_read(fn, 1);
	if (!n) { printfx("cannot read '%s'\n", fn); return NULL; }
	in_hdr ihdr = set_input_hdr(n);
	if ((n->datatype != DT_FLOAT32 || (n->scl_slope != 0.0f && n->scl_slope != 1.0f) || n->scl_inter != 0.0f)
	    && nifti_image_change_datatype(n, DT_FLOAT32, &ihdr) != 0) {
		printfx("cannot convert '%s' to float32\n", fn);
		nifti_image_free(n);
		return NULL;
	}
	return n;
}

int nii_resample(nifti_image *nim, const resample_opts *o, const char *out, gzModes gz) {
	if (nim->datatype != DT_FLOAT32) {
		printfx("-resample: requires float32 input\n");
		return 1;
	}
	/* output grid: the grid's geometry, optionally re-voxelled isotropically over the same FOV */
	const nifti_image *g = o->grid;
	nifti_dmat44 S = best_affine(g);
	int64_t on[3] = {g->nx, g->ny, g->nz};
	if (o->res > 0.0) {
		for (int a = 0; a < 3; a++) {
			double scale = o->res / sqrt(S.m[0][a] * S.m[0][a] + S.m[1][a] * S.m[1][a] + S.m[2][a] * S.m[2][a]);
			double na = (double)on[a] / scale;
			if (!(na <= INT32_MAX)) na = INT32_MAX;   /* rejected below; keeps lround defined */
			on[a] = (int64_t)lround(na);
			if (on[a] < 1) on[a] = 1;
			/* keep the outer corner (voxel -0.5) fixed: new centre 0 sits half a new voxel inside */
			double shift = 0.5 * (scale - 1.0);
			for (int r = 0; r < 3; r++) {
				S.m[r][3] += S.m[r][a] * shift;
				S.m[r][a] *= scale;
			}
		}
	}
	if ((double)on[0] * on[1] * on[2] > INT32_MAX) {
		printfx("-resample: output grid %lldx%lldx%lld is too large\n", (long long)on[0], (long long)on[1], (long long)on[2]);
		return 1;
	}
	int64_t onv = on[0] * on[1] * on[2];
	int nt = (nim->ndim > 3 && nim->nt > 1) ? (int)nim->nt : 1;
	nifti_dmat44 bold_ijk = nifti_dmat44_inverse(best_affine(nim));
	nifti_dmat44 warp_ijk;
	if (o->warp) warp_ijk = nifti_dmat44_inverse(best_affine(o->warp));
	nifti_dmat44 bwarp_ijk;
	if (o->bwarp) bwarp_ijk = nifti_dmat44_inverse(best_affine(o->bwarp));
	float *ref = alloc_f32(onv, 3);
	float *dst = alloc_f32(onv, nt);
	double (*hmcM)[12] = o->hmc ? hmc_pulls(nim, o->hmc, nt) : NULL;
	if (!ref || !dst || (o->hmc && !hmcM)) {
		free(ref); free(dst); free(hmcM);
		printfx("-resample: out of memory or bad motion parameters\n");
		return 1;
	}
	/* volume-independent part of the pull: output voxel -> input (boldref) voxel */
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (int64_t v = 0; v < onv; v++) {
		double ijk[3] = {(double)(v % on[0]), (double)((v / on[0]) % on[1]), (double)(v / (on[0] * on[1]))};
		double x[3], y[3];
		apply44(&S, ijk, x);
		if (o->warp) {
			double u[3];
			warp_vector(o->warp, &warp_ijk, x, u);
			double q[3] = {x[0] - u[0], x[1] - u[1], x[2] + u[2]};   /* LPS vector -> RAS */
			apply44d(o->aff, q, x);
		}
		if (o->ftm) { apply44d(o->ftm, x, y); memcpy(x, y, sizeof(x)); }
		if (o->bwarp) warp_point(o->bwarp, &bwarp_ijk, x);
		apply44(&bold_ijk, x, y);
		for (int a = 0; a < 3; a++) ref[v * 3 + a] = (float)y[a];
	}
	int bad = pull_series(nim, nt, ref, onv, (const double (*)[12])hmcM, o->vsm, dst);
	free(ref); free(hmcM);
	if (bad) { free(dst); printfx("-resample: out of memory\n"); return 1; }
	/* write: borrow nim's header (TR, units, intent) with the output geometry, then restore it */
	nifti_image saved = *nim;
	nim->nx = nim->dim[1] = on[0]; nim->ny = nim->dim[2] = on[1]; nim->nz = nim->dim[3] = on[2];
	nim->nt = nim->dim[4] = nt; nim->nvox = onv * nt;
	nim->ndim = nim->dim[0] = nt > 1 ? 4 : 3;
	for (int a = 0; a < 3; a++)
		nim->pixdim[a + 1] = sqrt(S.m[0][a] * S.m[0][a] + S.m[1][a] * S.m[1][a] + S.m[2][a] * S.m[2][a]);
	nim->dx = nim->pixdim[1]; nim->dy = nim->pixdim[2]; nim->dz = nim->pixdim[3];
	nim->xyz_units = NIFTI_UNITS_MM;   /* S is mm; time_units is a separate field and survives */
	nim->sto_xyz = S; nim->sto_ijk = nifti_dmat44_inverse(S);
	nim->sform_code = g->sform_code > 0 ? g->sform_code : (g->qform_code > 0 ? g->qform_code : 1);
	nim->qform_code = nim->sform_code;
	nim->qto_xyz = S; nim->qto_ijk = nim->sto_ijk;
	nifti_dmat44_to_quatern(S, &nim->quatern_b, &nim->quatern_c, &nim->quatern_d,
	                        &nim->qoffset_x, &nim->qoffset_y, &nim->qoffset_z, NULL, NULL, NULL, &nim->qfac);
	nim->scl_slope = 1.0f; nim->scl_inter = 0.0f;
	nim->data = dst;
	nim->fname = nim->iname = (char *)out;
	int rc = nifti_save(nim, "", gz);
	*nim = saved;
	free(dst);
	return rc;
}

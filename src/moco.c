// moco.c - rigid-body motion correction for 4D datasets (-moco)
//
// Clean-room implementation of Cox RW & Jesmanowicz A, "Real-Time 3D Image Registration for
// Functional MRI", Magn Reson Med 42:1014-1018 (1999).  AFNI's GPL-2 3dvolreg/3drotate sources
// were not read; they served only as black-box oracles.  See moco.h and
// the moco_bench repository's test/moco_reference_manifest.md.
//
// Structure, following the paper:
//   * a rigid transform is factored into FOUR 3D shears (Appendix), each of which displaces
//     data along one index axis only, so every interpolation is a constant shift of a
//     contiguous row -- an 8-tap Lagrange kernel evaluated ONCE per row (Eq [1]-[5]);
//   * translations fold into the same four passes (Eq [6], [7]);
//   * the six motion parameters are fit by repeated linearization of
//     E(a) = sum_x w(x) [ J(T[a]x) - I(x) ]^2 , with w a smoothed copy of the base.
//
// Conventions below marked "measured" were established against the pinned oracle; see the
// manifest for the experiment behind each one.

#define _USE_MATH_DEFINES // microsoft compiler: gates M_PI in <math.h>, must precede it
#include <math.h>
#include <stddef.h>   /* ptrdiff_t: do NOT rely on <omp.h> to drag this in */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#ifdef _WIN32
	#include <io.h>
	#include <process.h>
	#define getpid _getpid
#else
	#include <unistd.h>
#endif
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "moco.h"
#include "print.h"
#include "core32.h"   /* nifti_smooth_gauss_f32 */

#define MOCO_ORDER 8              // heptic = 7th-order Lagrange = 8 taps (measured: the default)
#define MOCO_HALF  3              // taps span j-3 .. j+4
#define MOCO_PAD   4              // measured: 3dvolreg zero-pads 4 voxels on each of 6 planes
#define MOCO_MAXITE 23            // documented default
#define MOCO_XTHRESH 0.010        // documented default, voxels
#define MOCO_RTHRESH 0.020        // documented default, degrees
#define MOCO_DELTA 0.700          // documented default, voxels (finite-difference step)
#define MOCO_WSIGMA 3.0           // measured: weight blur sigma = 3 voxels (edt_blur sigx/dx)
#define MOCO_WTRUNC 2.5           // measured: edt_blur "sfac 2.5" truncation
#define MOCO_WCLIP 0.025        // measured: weight below 2.5% of its max is zeroed (wtrim)
#define MOCO_EDGING 0.05          // documented: -edging default is 5% of each brick size
#define MOCO_COST_RISE 0.05       // reject a step only if the weighted cost rises by >5% (divergence)

typedef double mo33[3][3];

// ------------------------------------------------------------------ small matrix helpers
static void m_ident(mo33 a) {
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) a[i][j] = (i == j) ? 1.0 : 0.0;
}

static void m_mul(const mo33 a, const mo33 b, mo33 out) {
	mo33 t;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			double s = 0.0;
			for (int k = 0; k < 3; k++) s += a[i][k] * b[k][j];
			t[i][j] = s;
		}
	memcpy(out, t, sizeof(mo33));
}

static double m_det(const mo33 a) {
	return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
	     - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
	     + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

static int m_inv(const mo33 a, mo33 out) {
	double d = m_det(a);
	if (!(fabs(d) > 1e-30)) return 1;
	double r = 1.0 / d;
	mo33 t;
	t[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * r;
	t[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * r;
	t[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * r;
	t[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * r;
	t[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * r;
	t[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * r;
	t[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * r;
	t[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * r;
	t[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * r;
	memcpy(out, t, sizeof(mo33));
	return 0;
}

static void m_vec(const mo33 a, const double v[3], double out[3]) {
	double t[3];
	for (int i = 0; i < 3; i++) t[i] = a[i][0] * v[0] + a[i][1] * v[1] + a[i][2] * v[2];
	memcpy(out, t, sizeof(t));
}

// ------------------------------------------------------------------ rotation from parameters
// Measured convention: AFNI DICOM axes are +x = Left, +y = Posterior, +z = Superior;
// roll rotates about +z (I-S), pitch about +x (R-L), yaw about +y (A-P), right-hand rule, and
// the three compose as R = Ry(yaw) . Rx(pitch) . Rz(roll).
static void moco_rot(double roll_deg, double pitch_deg, double yaw_deg, mo33 out) {
	const double d2r = M_PI / 180.0;
	double cr = cos(roll_deg * d2r), sr = sin(roll_deg * d2r);
	double cp = cos(pitch_deg * d2r), sp = sin(pitch_deg * d2r);
	double cy = cos(yaw_deg * d2r), sy = sin(yaw_deg * d2r);
	mo33 Rz = {{cr, -sr, 0}, {sr, cr, 0}, {0, 0, 1}};
	mo33 Rx = {{1, 0, 0}, {0, cp, -sp}, {0, sp, cp}};
	mo33 Ry = {{cy, 0, sy}, {0, 1, 0}, {-sy, 0, cy}};
	m_mul(Ry, Rx, out);
	m_mul(out, Rz, out);
}

// ------------------------------------------------------------------ four-shear factorization
// A shear along axis `ax` adds alpha*x_u + beta*x_v to coordinate ax, where (u,v) are the other
// two axes in increasing order.
typedef struct {
	int ax;
	double alpha, beta;
} moco_shear;

static void moco_others(int a, int *u, int *v) {
	*u = (a == 0) ? 1 : 0;
	*v = (a == 2) ? 1 : 2;
}

static void shear_mat(const moco_shear *s, mo33 out) {
	int u, v;
	moco_others(s->ax, &u, &v);
	m_ident(out);
	out[s->ax][u] = s->alpha;
	out[s->ax][v] = s->beta;
}

// Solve a 2x2 system.  Eq [4] and Eq [2] can be singular yet consistent (the paper's 2D-rotation
// class); fall back to the minimum-norm least-squares solution there and let the caller's
// reconstruction check accept or reject it.
static int solve2(double m00, double m01, double m10, double m11, double r0, double r1,
                  double *x0, double *x1) {
	double det = m00 * m11 - m01 * m10;
	if (fabs(det) > 1e-11) {
		*x0 = (r0 * m11 - m01 * r1) / det;
		*x1 = (m00 * r1 - r0 * m10) / det;
		return 0;
	}
	// minimum-norm solution of the consistent rank-deficient system via the pseudo-inverse of
	// a 2x2 built from its non-zero row
	double n0 = m00 * m00 + m01 * m01, n1 = m10 * m10 + m11 * m11;
	if (n0 >= n1) {
		if (!(n0 > 1e-24)) return 1;
		*x0 = m00 * r0 / n0;
		*x1 = m01 * r0 / n0;
	} else {
		if (!(n1 > 1e-24)) return 1;
		*x0 = m10 * r1 / n1;
		*x1 = m11 * r1 / n1;
	}
	// consistency check on the row we discarded
	double c0 = m00 * (*x0) + m01 * (*x1) - r0;
	double c1 = m10 * (*x0) + m11 * (*x1) - r1;
	double sc = 1.0 + fabs(r0) + fabs(r1);
	if (fabs(c0) > 1e-9 * sc || fabs(c1) > 1e-9 * sc) return 1;
	return 0;
}

// Canonical factorization A = S3 S1 S2 S3 (the paper's (i,j,k) = (3,1,2) ordering).
static int factor_canonical(const mo33 A, mo33 S[4]) {
	double a2 = A[0][1], a3 = A[0][2];
	double b2 = A[1][1], b3 = A[1][2];
	double c2 = A[2][1], c3 = A[2][2];
	if (!(fabs(b3) > 1e-11)) return 1;
	double al0, be0;
	if (solve2(a2, b2, a3, b3, c2 - (b2 - 1.0) / b3, c3 - 1.0, &al0, &be0)) return 1;   // Eq [4]
	moco_shear s0 = {2, al0, be0};
	shear_mat(&s0, S[0]);
	mo33 inv0, A1;
	if (m_inv(S[0], inv0)) return 1;
	m_mul(inv0, A, A1);
	double al1, be1;
	if (solve2(A1[1][1], A1[2][1], A1[1][2], A1[2][2], A1[0][1], A1[0][2], &al1, &be1)) return 1; // Eq [2]
	moco_shear s1 = {0, al1, be1};
	shear_mat(&s1, S[1]);
	mo33 inv1, A2;
	if (m_inv(S[1], inv1)) return 1;
	m_mul(inv1, A1, A2);
	moco_shear s2 = {1, A2[1][0] - A2[1][2] * A2[2][0], A2[1][2]};
	shear_mat(&s2, S[2]);
	mo33 inv2, A3;
	if (m_inv(S[2], inv2)) return 1;
	m_mul(inv2, A2, A3);
	moco_shear s3 = {2, A3[2][0], A3[2][1]};
	shear_mat(&s3, S[3]);
	// A correct factorization reconstructs A exactly; make that the acceptance test.
	mo33 rec;
	m_mul(S[0], S[1], rec);
	m_mul(rec, S[2], rec);
	m_mul(rec, S[3], rec);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			if (fabs(rec[i][j] - A[i][j]) > 1e-9) return 1;
	return 0;
}

// Permutation sending real axis i -> canonical 2 (z), j -> 0 (x), k -> 1 (y).
static void perm_mat(int i, int j, int k, mo33 P) {
	memset(P, 0, sizeof(mo33));
	P[2][i] = 1.0;
	P[0][j] = 1.0;
	P[1][k] = 1.0;
}

// Try all six axis orderings, keep the one minimising max(|alpha|,|beta|).  Exactly-planar
// rotations tie: the survivors are sign-mirrors, and AFNI takes the paper's Eq [1] form whose
// leading shear is -tan(phi/2) < 0 (measured).
static int moco_factor(const mo33 A, moco_shear out[4]) {
	int isI = 1;
	for (int i = 0; i < 3 && isI; i++)
		for (int j = 0; j < 3; j++)
			if (fabs(A[i][j] - ((i == j) ? 1.0 : 0.0)) > 1e-14) { isI = 0; break; }
	if (isI) {   // paper: "if A = I, no factorization is needed" -- keep passes for translation
		int ax[4] = {0, 1, 2, 0};
		for (int m = 0; m < 4; m++) { out[m].ax = ax[m]; out[m].alpha = out[m].beta = 0.0; }
		return 0;
	}
	int have = 0;
	double best_dist = 0.0;
	int best_lead = 0;
	for (int i = 0; i < 3; i++) {
		int u, v;
		moco_others(i, &u, &v);
		int jk[2][2] = {{u, v}, {v, u}};
		for (int q = 0; q < 2; q++) {
			int j = jk[q][0], k = jk[q][1];
			mo33 P, Pt, tmp, At, S[4];
			perm_mat(i, j, k, P);
			for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) Pt[r][c] = P[c][r];
			m_mul(P, A, tmp);
			m_mul(tmp, Pt, At);
			if (factor_canonical(At, S)) continue;
			int axes[4] = {i, j, k, i};
			moco_shear cand[4];
			int bad = 0;
			for (int m = 0; m < 4; m++) {
				mo33 back;
				m_mul(Pt, S[m], tmp);
				m_mul(tmp, P, back);
				int uu, vv;
				moco_others(axes[m], &uu, &vv);
				cand[m].ax = axes[m];
				cand[m].alpha = back[axes[m]][uu];
				cand[m].beta = back[axes[m]][vv];
				mo33 chk;
				shear_mat(&cand[m], chk);
				for (int r = 0; r < 3 && !bad; r++)
					for (int c = 0; c < 3; c++)
						if (fabs(chk[r][c] - back[r][c]) > 1e-9) { bad = 1; break; }
			}
			if (bad) continue;
			double dist = 0.0;
			for (int m = 0; m < 4; m++) {
				double a = fabs(cand[m].alpha), b = fabs(cand[m].beta);
				if (a > dist) dist = a;
				if (b > dist) dist = b;
			}
			int lead = ((cand[0].alpha + cand[0].beta) < 0.0) ? 0 : 1;
			if (!have || dist < best_dist - 1e-12 ||
			    (fabs(dist - best_dist) <= 1e-12 && lead < best_lead)) {
				have = 1;
				best_dist = dist;
				best_lead = lead;
				memcpy(out, cand, sizeof(cand));
			}
		}
	}
	return have ? 0 : 1;
}

// ------------------------------------------------------------------ 1D Lagrange row shift
// Weights for MOCO_ORDER nodes 0..ORDER-1 evaluated at position t.
static void lagrange_w(double t, double *w) {
	for (int m = 0; m < MOCO_ORDER; m++) {
		double p = 1.0;
		for (int q = 0; q < MOCO_ORDER; q++)
			if (q != m) p *= (t - (double)q) / ((double)m - (double)q);
		w[m] = p;
	}
}

// out[i] = in[i - shift] along a strided row of length n; zero outside.  `shift` is constant
// along the row, so the stencil offset and the eight weights are computed once (the paper's
// efficiency argument for shear-based rotation).
static void shift_row(const float *src, float *dst, int n, ptrdiff_t stride, double shift) {
	if (shift == 0.0) {
		for (int i = 0; i < n; i++) dst[(ptrdiff_t)i * stride] = src[(ptrdiff_t)i * stride];
		return;
	}
	/* (int)floor(ps) is undefined for NaN or an out-of-range magnitude, and i + lo then
	   overflows signed int.  A divergent fit is the only way to get here, so fail closed. */
	if (!(shift > -(double)(n + MOCO_ORDER) && shift < (double)(n + MOCO_ORDER))) {
		for (int i = 0; i < n; i++) dst[(ptrdiff_t)i * stride] = 0.0f;
		return;
	}
	double ps = -shift;
	double fl = floor(ps);
	int j0 = (int)fl;
	double w[MOCO_ORDER];
	lagrange_w((ps - fl) + (double)MOCO_HALF, w);
	int lo = j0 - MOCO_HALF;
	for (int i = 0; i < n; i++) {
		int base = i + lo;
		double acc = 0.0;
		if (base >= 0 && base + MOCO_ORDER - 1 < n) {
			const float *s = src + (ptrdiff_t)base * stride;
			for (int m = 0; m < MOCO_ORDER; m++) acc += w[m] * (double)s[(ptrdiff_t)m * stride];
		} else {
			for (int m = 0; m < MOCO_ORDER; m++) {
				int q = base + m;
				if (q >= 0 && q < n) acc += w[m] * (double)src[(ptrdiff_t)q * stride];
			}
		}
		dst[(ptrdiff_t)i * stride] = (float)acc;
	}
}

// One shear pass: displace along `ax` by alpha*(u-cu) + beta*(v-cv) + delta.
static void shear_pass(float *vol, const int dim[3], const double ctr[3], const moco_shear *sh,
                      double delta, float *rowbuf) {
	int ax = sh->ax, u, v;
	moco_others(ax, &u, &v);
	ptrdiff_t str[3] = {1, dim[0], (ptrdiff_t)dim[0] * dim[1]};
	int n = dim[ax];
	for (int iv = 0; iv < dim[v]; iv++) {
		for (int iu = 0; iu < dim[u]; iu++) {
			double s = sh->alpha * ((double)iu - ctr[u]) + sh->beta * ((double)iv - ctr[v]) + delta;
			ptrdiff_t off = (ptrdiff_t)iu * str[u] + (ptrdiff_t)iv * str[v];
			float *row = vol + off;
			for (int i = 0; i < n; i++) rowbuf[i] = row[(ptrdiff_t)i * str[ax]];
			// rowbuf is contiguous, so the shift runs with unit stride
			shift_row(rowbuf, rowbuf + n, n, 1, s);
			for (int i = 0; i < n; i++) row[(ptrdiff_t)i * str[ax]] = rowbuf[n + i];
		}
	}
}

// ------------------------------------------------------------------ full rigid warp
typedef struct {
	int dim[3];           // unpadded dims
	int pdim[3];          // padded dims
	double ctr[3];        // grid centre of the PADDED volume, index units
	mo33 idx2mm;         // index -> AFNI DICOM mm (linear part)
	mo33 mm2idx;
	size_t pn;            // padded voxel count
} moco_geom;

// out (unpadded, dim) = in (unpadded, dim) resampled through the rigid transform given by the
// six parameters.  scratch must hold pn floats; rowbuf 2*max(pdim) floats.
// (R, s) act about the grid centre: T(x) = R(x - c) + c + s.  Composing T2 after T1 gives
// R = R2 R1, s = R2 s1 + s2.
static void moco_compose(const mo33 R2, const double s2[3], const mo33 R1, const double s1[3],
                         mo33 Ro, double so[3]) {
	mo33 Rt;
	m_mul(R2, R1, Rt);
	double t[3];
	m_vec(R2, s1, t);
	for (int i = 0; i < 3; i++) so[i] = t[i] + s2[i];
	memcpy(Ro, Rt, sizeof(mo33));
}

static void moco_invert(const mo33 R, const double s[3], mo33 Ro, double so[3]) {
	mo33 Rt;
	for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Rt[i][j] = R[j][i];  // orthogonal
	double t[3];
	m_vec(Rt, s, t);
	for (int i = 0; i < 3; i++) so[i] = -t[i];
	memcpy(Ro, Rt, sizeof(mo33));
}

// Inverse of moco_rot for R = Ry(yaw) Rx(pitch) Rz(roll).
static void moco_unrot(const mo33 R, double *roll, double *pitch, double *yaw) {
	const double r2d = 180.0 / M_PI;
	double sp = -R[1][2];
	if (sp > 1.0) sp = 1.0;
	if (sp < -1.0) sp = -1.0;
	*pitch = asin(sp) * r2d;
	*roll = atan2(R[1][0], R[1][1]) * r2d;
	*yaw = atan2(R[0][2], R[2][2]) * r2d;
}

static int moco_warp_rs(const float *in, float *out, const moco_geom *g, const mo33 R,
                        const double s_mm[3], float *pad, float *rowbuf) {
	mo33 tmp, A;
	m_mul(g->mm2idx, R, tmp);
	m_mul(tmp, g->idx2mm, A);
	double s_idx[3];
	m_vec(g->mm2idx, s_mm, s_idx);
	moco_shear sh[4];
	if (moco_factor(A, sh)) return 1;
	// Translation folding, Eq [7]: t0 + S0 t1 + S0 S1 t2 = s, each t_m along its own axis.
	mo33 M, S0, S1;
	shear_mat(&sh[0], S0);
	shear_mat(&sh[1], S1);
	mo33 S0S1;
	m_mul(S0, S1, S0S1);
	for (int r = 0; r < 3; r++) {
		M[r][0] = (r == sh[0].ax) ? 1.0 : 0.0;
		M[r][1] = S0[r][sh[1].ax];
		M[r][2] = S0S1[r][sh[2].ax];
	}
	mo33 Minv;
	if (m_inv(M, Minv)) return 1;
	double d[3];
	m_vec(Minv, s_idx, d);
	// pad
	memset(pad, 0, g->pn * sizeof(float));
	for (int k = 0; k < g->dim[2]; k++)
		for (int j = 0; j < g->dim[1]; j++) {
			const float *sp = in + ((size_t)k * g->dim[1] + j) * g->dim[0];
			float *dp = pad + (((size_t)(k + MOCO_PAD) * g->pdim[1] + (j + MOCO_PAD)) * g->pdim[0]) + MOCO_PAD;
			memcpy(dp, sp, (size_t)g->dim[0] * sizeof(float));
		}
	// out = in o T^-1 with T = T0 T1 T2 T3, so the passes unwind in reverse order.
	for (int m = 3; m >= 0; m--)
		shear_pass(pad, g->pdim, g->ctr, &sh[m], (m < 3) ? d[m] : 0.0, rowbuf);
	// crop
	for (int k = 0; k < g->dim[2]; k++)
		for (int j = 0; j < g->dim[1]; j++) {
			const float *sp = pad + (((size_t)(k + MOCO_PAD) * g->pdim[1] + (j + MOCO_PAD)) * g->pdim[0]) + MOCO_PAD;
			float *dp = out + ((size_t)k * g->dim[1] + j) * g->dim[0];
			memcpy(dp, sp, (size_t)g->dim[0] * sizeof(float));
		}
	return 0;
}

static int moco_warp(const float *in, float *out, const moco_geom *g, const double par[6],
                     float *pad, float *rowbuf) {
	mo33 R;
	moco_rot(par[0], par[1], par[2], R);
	double s_mm[3] = {par[4], par[5], par[3]};      // (dL, dP, dS) -> DICOM (x, y, z)
	return moco_warp_rs(in, out, g, R, s_mm, pad, rowbuf);
}

// ------------------------------------------------------------------ weight image
// Separable Gaussian, sigma in voxels, truncated at MOCO_WTRUNC sigma (measured: edt_blur
// reports sigx = 3 voxels and sfac 2.5).  Normalised, edges renormalised.
static int moco_weight(const float *base, float *wt, const int dim[3], const double vox[3]) {
	size_t n = (size_t)dim[0] * dim[1] * dim[2];
	for (size_t i = 0; i < n; i++) wt[i] = base[i];
	/* Same kernel the reference reports (edt_blur sigx 3 voxels, sfac 2.5): niimath's exported
	   smoother uses exp(-i^2/2 sigma^2) truncated at ceil(2.5 sigma) with edge renormalisation
	   by the truncated mass, which is the identical formula.  Passing 1 mm voxel sizes makes
	   sigma_mm == sigma_voxels. */
	/* edt_blur reports a single sigma in MILLIMETRES derived from dx (sigx = 3*dx on all three
	   datasets), applied to every axis -- so on anisotropic voxels the per-axis sigma in VOXELS
	   differs from 3. */
	{
		double smm = MOCO_WSIGMA * vox[0];
		if (nifti_smooth_gauss_f32(wt, dim[0], dim[1], dim[2], 1, 1.0f, 1.0f, 1.0f,
		                           (float)(smm / vox[0]), (float)(smm / vox[1]),
		                           (float)(smm / vox[2]), (float)MOCO_WTRUNC) != 0)
			return 1;
	}
	// non-negative weight (the paper requires w >= 0)
	for (size_t i = 0; i < n; i++) if (!(wt[i] > 0.0f)) wt[i] = 0.0f;
	/* Documented: -edging zeroes the default weight in a border around the base.  With no
	   -edging and no AFNI_VOLREG_EDGING the default is 5% of each brick size, rounded to
	   nearest -- the reference reports exactly "Edging: x=4 y=4 z=2" for a 76 x 76 x 45 grid. */
	int eg[3];
	for (int i = 0; i < 3; i++) {
		eg[i] = (int)(MOCO_EDGING * (double)dim[i] + 0.5);
		if (eg[i] * 2 >= dim[i]) eg[i] = (dim[i] - 1) / 2;   /* never zero the whole brick */
		if (eg[i] < 0) eg[i] = 0;
	}
	for (int k = 0; k < dim[2]; k++)
		for (int j = 0; j < dim[1]; j++)
			for (int i = 0; i < dim[0]; i++)
				if (i < eg[0] || i >= dim[0] - eg[0] ||
				    j < eg[1] || j >= dim[1] - eg[1] ||
				    k < eg[2] || k >= dim[2] - eg[2])
					wt[((size_t)k * dim[1] + j) * dim[0] + i] = 0.0f;
	/* MEASURED: the reference zeroes weight below MOCO_WCLIP of the maximum and then fits only
	   the bounding box of what survives ("wtrim" in -verbose).  Recovered by inverting the
	   reported trim boxes: 0.0250-0.0255 on all three validation datasets, with the resulting
	   zero-weight fraction matching (21.9 vs 21.7, 18.9 vs 17.0, 30.7 vs 29.9 per cent).  The box
	   itself is only an optimisation -- voxels outside it already carry zero weight. */
	float wmax = 0.0f;
	for (size_t i = 0; i < n; i++) if (wt[i] > wmax) wmax = wt[i];
	float wcut = (float)(MOCO_WCLIP * (double)wmax);
	for (size_t i = 0; i < n; i++) if (wt[i] <= wcut) wt[i] = 0.0f;
	return 0;
}

// ------------------------------------------------------------------ 6x6 Cholesky
static int chol6(const double M[6][6], const double b[6], double x[6]) {
	double L[6][6];
	memset(L, 0, sizeof(L));
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j <= i; j++) {
			double s = M[i][j];
			for (int k = 0; k < j; k++) s -= L[i][k] * L[j][k];
			if (i == j) {
				if (!(s > 1e-20)) return 1;      // singular / not positive definite
				L[i][i] = sqrt(s);
			} else {
				L[i][j] = s / L[j][j];
			}
		}
	}
	double y[6];
	for (int i = 0; i < 6; i++) {
		double s = b[i];
		for (int k = 0; k < i; k++) s -= L[i][k] * y[k];
		y[i] = s / L[i][i];
	}
	for (int i = 5; i >= 0; i--) {
		double s = y[i];
		for (int k = i + 1; k < 6; k++) s -= L[k][i] * x[k];
		x[i] = s / L[i][i];
	}
	return 0;
}

// ------------------------------------------------------------------ geometry from the header
// niimath stores voxel->world in RAS+; AFNI DICOM is +x Left, +y Posterior, +z Superior, so the
// linear part is negated in x and y.
static int moco_geom_init(const nifti_image *nim, moco_geom *g) {
	nifti_dmat44 m = (nim->sform_code > 0 && !(nim->sform_code < nim->qform_code)) ? nim->sto_xyz
	        : (nim->qform_code > 0 ? nim->qto_xyz : nim->sto_xyz);
	double pix[3] = {fabs(nim->dx), fabs(nim->dy), fabs(nim->dz)};
	for (int i = 0; i < 3; i++) if (!(pix[i] > 0.0)) pix[i] = 1.0;
	memset(g->idx2mm, 0, sizeof(mo33));
	if (nim->sform_code <= 0 && nim->qform_code <= 0) {
		for (int i = 0; i < 3; i++) g->idx2mm[i][i] = pix[i];
	} else {
		/* Measured: the reference ignores obliquity and works in the dataset's storage axes --
		   the same geometry `3dinfo -d3` reports.  Snap each index axis to the DICOM axis it is
		   most aligned with (AFNI DICOM is +x Left, +y Posterior, +z Superior, so the RAS+
		   sform's x and y rows are negated) and keep only the sign and the voxel size.  Using
		   the oblique sform here biases the z-coupled parameters by ~20%. */
		double M[3][3];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				M[i][j] = (i < 2 ? -1.0 : 1.0) * (double)m.m[i][j];
		/* Score all six permutations globally: a greedy column-by-column pick can miss the
		   closest signed permutation when the obliquity is strong. */
		static const int PERM[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
		int bestp = -1;
		double bestscore = -1.0;
		for (int q = 0; q < 6; q++) {
			double sc = 0.0;
			for (int j = 0; j < 3; j++) sc += fabs(M[PERM[q][j]][j]);
			if (sc > bestscore) { bestscore = sc; bestp = q; }
		}
		if (bestp < 0) return 1;
		for (int j = 0; j < 3; j++) {
			int i = PERM[bestp][j];
			g->idx2mm[i][j] = (M[i][j] < 0.0 ? -1.0 : 1.0) * pix[j];
		}
	}
	if (m_inv(g->idx2mm, g->mm2idx)) return 1;
	for (int i = 0; i < 3; i++) {
		g->pdim[i] = g->dim[i] + 2 * MOCO_PAD;
		g->ctr[i] = ((double)g->pdim[i] - 1.0) * 0.5;
	}
	g->pn = (size_t)g->pdim[0] * g->pdim[1] * g->pdim[2];
	return 0;
}

// ------------------------------------------------------------------ external reference image
/* Read `fn` as a float32 image and check it can serve as the registration base for `nim`.
   Returns the image (caller frees with nifti_image_free) or NULL after printing why.
   -moco works entirely inside the input's voxel grid -- the four shear passes shift rows of that
   grid -- so an external reference is usable only when it IS that grid.  Anything else would need
   a resampling step whose interpolation is not part of the measured contract, and would leave the
   .1D parameters referring to a grid the caller never supplied; reject instead of resampling. */
static nifti_image *moco_read_ref(const char *fn, nifti_image *nim) {
	{	/* Header-only preflight: reject a malformed or oversized reference BEFORE its payload is
		   decompressed and allocated. */
		nifti_image *h = nifti_image_read(fn, 0);
		int bad = 0;
		if (!h) {
			printfx("-moco: failed to read the header of reference image '%s'\n", fn);
			return NULL;
		}
		if (h->nvox < 1 || h->nx < 1 || h->ny < 1 || h->nz < 1) {
			printfx("-moco: reference image '%s' has invalid dimensions\n", fn);
			bad = 1;
		} else if (h->nu > 1 || h->nv > 1 || h->nw > 1) {
			printfx("-moco: reference image '%s' has more than 4 dimensions\n", fn);
			bad = 1;
		} else if ((int64_t)h->nvox > INT_MAX) {
			printfx("-moco: reference image '%s' exceeds INT_MAX voxels; -moco is not a huge-image-safe operation\n", fn);
			bad = 1;
		}
		nifti_image_free(h);
		if (bad) return NULL;
	}
	nifti_image *ref = nifti_image_read(fn, 1);
	if (!ref) {
		printfx("-moco: failed to read reference image '%s'\n", fn);
		return NULL;
	}
	/* Re-check after the load as well as before it: the preflight is what avoids decompressing a
	   huge payload, this is the fail-closed guarantee if the file changed between the two reads. */
	if (ref->nvox < 1 || ref->nx < 1 || ref->ny < 1 || ref->nz < 1 ||
	    ref->nu > 1 || ref->nv > 1 || ref->nw > 1 || (int64_t)ref->nvox > INT_MAX) {
		printfx("-moco: reference image '%s' changed on disk or has unusable dimensions\n", fn);
		nifti_image_free(ref);
		return NULL;
	}
	if (ref->nx != nim->nx || ref->ny != nim->ny || ref->nz != nim->nz) {
		printfx("-moco: reference image '%s' is not on the same grid as the input "
		        "(%lldx%lldx%lld vs %lldx%lldx%lld). -moco registers within the input voxel grid; "
		        "reslice the reference onto the input first.\n",
		        fn, (long long)ref->nx, (long long)ref->ny, (long long)ref->nz,
		        (long long)nim->nx, (long long)nim->ny, (long long)nim->nz);
		nifti_image_free(ref);
		return NULL;
	}
	/* Same 0.001 mm corner-displacement gate --qc and --medic use for their own same-grid
	   requirement: it covers rotation, scale, origin and xyz_units in one call. */
	float disp = max_displacement_mm(ref, nim);
	if (!(disp <= 0.001f)) {
		printfx("-moco: reference image '%s' has the same dimensions as the input but a different "
		        "voxel-to-world transform (corners differ by up to %g mm). -moco registers within "
		        "the input voxel grid; reslice the reference onto the input first.\n",
		        fn, (double)disp);
		nifti_image_free(ref);
		return NULL;
	}
	in_hdr ihdr = set_input_hdr(ref);
	/* Convert when the stored type is not float32, but ALSO when it IS float32 and carries a
	   non-trivial scl_slope/scl_inter -- otherwise a scaled float32 reference is used raw, which
	   would bias the intensity scale the fit profiles out. */
	if (ref->datatype != DT_FLOAT32 ||
	    (ref->scl_slope != 0.0f && ref->scl_slope != 1.0f) || ref->scl_inter != 0.0f) {
		if (nifti_image_change_datatype(ref, DT_FLOAT32, &ihdr) != 0) {
			printfx("-moco: failed to convert reference image '%s' to float32\n", fn);
			nifti_image_free(ref);
			return NULL;
		}
	}
	if (!ref->data) {
		printfx("-moco: reference image '%s' has no voxel data\n", fn);
		nifti_image_free(ref);
		return NULL;
	}
	if ((ref->ndim > 3) && (ref->nt > 1))
		printfx("-moco: reference image '%s' is 4D; using its volume 0\n", fn);
	return ref;
}

// ------------------------------------------------------------------ entry point
int nii_moco(nifti_image *nim, const char *par_path, int ref_vol, const char *ref_file) {
	if (!nim || nim->datatype != DT_FLOAT32 || !nim->data) {
		printfx("-moco: internal error (expected float32 image)\n");
		return 1;
	}
	int nt = (nim->ndim > 3) ? nim->nt : 1;
	if (nim->ndim > 4) {
		printfx("-moco requires a scalar 4D image (got %lldD)\n", (long long)nim->ndim);
		return 1;
	}
	if (nim->ndim < 4 || nt < 2) {
		printfx("-moco requires a 4D image with more than one volume (got nt = %lld)\n", (long long)nt);
		return 1;
	}
	if (!ref_file && (ref_vol < 0 || ref_vol >= nt)) {
		printfx("-moco -ref %d is outside the input series (it has %lld volumes, 0..%lld)\n",
		        ref_vol, (long long)nt, (long long)(nt - 1));
		return 1;
	}
	moco_geom g;
	g.dim[0] = nim->nx; g.dim[1] = nim->ny; g.dim[2] = nim->nz;
	if (g.dim[0] < 1 || g.dim[1] < 1 || g.dim[2] < 1) {
		printfx("-moco: degenerate spatial dimensions\n");
		return 1;
	}
	if (moco_geom_init(nim, &g)) {
		printfx("-moco: singular voxel-to-world transform\n");
		return 1;
	}
	size_t nvol = (size_t)g.dim[0] * g.dim[1] * g.dim[2];
	/* wasm32 is a 32-bit size_t and the huge gate only bounds the VOXEL count, so these
	   products must be checked before they reach malloc (AGENTS.md: callers pre-compute
	   overflow-prone products with nii_mul_size). */
	size_t nb_out, nb_vol, nb_deriv, nb_pad;
	if (nii_mul_size(nvol, (size_t)nt * sizeof(float), &nb_out) ||
	    nii_mul_size(nvol, sizeof(float), &nb_vol) ||
	    nii_mul_size(nvol, 6 * sizeof(float), &nb_deriv) ||
	    nii_mul_size(g.pn, sizeof(float), &nb_pad)) {
		printfx("-moco: image too large for this build\n");
		return 1;
	}
	const float *img = (const float *)nim->data;

	/* Every cleanup-owned pointer is declared and NULLed here, before any `goto done`, so the
	   single cleanup block can never free an indeterminate pointer. */
	int rc = 0;
	nifti_image *ref = NULL;
	float *wt = NULL, *out = NULL, *deriv = NULL, *pad = NULL;
	float *rowbuf = NULL, *tmpA = NULL, *tmpB = NULL;
	double *par = NULL;
	size_t nb_row = 0;
	int mx = g.pdim[0] > g.pdim[1] ? g.pdim[0] : g.pdim[1];
	if (g.pdim[2] > mx) mx = g.pdim[2];
	if (nii_mul_size((size_t)mx, 2 * sizeof(float), &nb_row)) {
		printfx("-moco: image too large for this build\n");
		return 1;
	}
	/* The registration base.  `base_idx` is the sub-brick of the INPUT that is the base and is
	   therefore copied through unchanged with an all-zero parameter row; an external reference is
	   not a sub-brick of the input, so it is -1 there and every volume gets registered. */
	int base_idx = ref_file ? -1 : ref_vol;
	const float *base;
	if (ref_file) {
		ref = moco_read_ref(ref_file, nim);
		if (!ref) return 1;
		base = (const float *)ref->data;
	} else {
		base = img + (size_t)base_idx * nvol;
	}
	wt = (float *)malloc(nb_vol);
	out = (float *)malloc(nb_out);          /* every voxel is written below; see the memcpy of the
	                                           base sub-brick and the per-volume writes */
	par = (double *)calloc((size_t)nt * 6, sizeof(double));
	deriv = (float *)malloc(nb_deriv);
	pad = (float *)malloc(nb_pad);
	rowbuf = (float *)malloc(nb_row);
	tmpA = (float *)malloc(nb_vol);
	tmpB = (float *)malloc(nb_vol);
	if (!wt || !out || !par || !deriv || !pad || !rowbuf || !tmpA || !tmpB) {
		printfx("-moco: out of memory\n");
		rc = 1;
		goto done;
	}

	// Derivative images of the base: central differences at MOCO_DELTA voxels.  The rotation
	// step is the angle whose displacement at the largest in-plane radius equals MOCO_DELTA
	// voxels, so all six columns of the normal equations carry comparable magnitude.
	double vmm[3] = {0, 0, 0};
	for (int i = 0; i < 3; i++) {
		double s = 0.0;
		for (int r = 0; r < 3; r++) s += g.idx2mm[r][i] * g.idx2mm[r][i];
		vmm[i] = sqrt(s);
	}
	double dmm = MOCO_DELTA * (vmm[0] + vmm[1] + vmm[2]) / 3.0;
	double ext[3];
	for (int i = 0; i < 3; i++) ext[i] = 0.5 * (double)g.dim[i] * vmm[i];
	double lever[3];                       // radius seen by roll(z), pitch(x), yaw(y)
	lever[0] = sqrt(ext[0] * ext[0] + ext[1] * ext[1]);
	lever[1] = sqrt(ext[1] * ext[1] + ext[2] * ext[2]);
	lever[2] = sqrt(ext[0] * ext[0] + ext[2] * ext[2]);
	double step[6];
	for (int p = 0; p < 3; p++)
		step[p] = (lever[p] > 1e-9) ? (dmm / lever[p]) * (180.0 / M_PI) : 1.0;
	/* MEASURED from -verbose: d/dx, d/dy, d/dz use delta = 0.7 * the voxel size of that DICOM
	   axis (1.96184 = 0.7*2.80263 and 1.96 = 0.7*2.8 on bold1), not 0.7 * the mean. */
	double dvox[3] = {0, 0, 0};
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) dvox[i] += fabs(g.idx2mm[i][j]);
	step[3] = MOCO_DELTA * dvox[2];   /* dS -> DICOM z */
	step[4] = MOCO_DELTA * dvox[0];   /* dL -> DICOM x */
	step[5] = MOCO_DELTA * dvox[1];   /* dP -> DICOM y */

	if (moco_weight(base, wt, g.dim, vmm)) {
		printfx("-moco: out of memory building the registration weight\n");
		rc = 1;
		goto done;
	}
	for (int p = 0; p < 6 && !rc; p++) {
		double pp[6] = {0, 0, 0, 0, 0, 0}, pm[6] = {0, 0, 0, 0, 0, 0};
		pp[p] = step[p];
		pm[p] = -step[p];
		if (moco_warp(base, tmpA, &g, pp, pad, rowbuf) ||
		    moco_warp(base, tmpB, &g, pm, pad, rowbuf)) { rc = 1; break; }
		float *d = deriv + (size_t)p * nvol;
		double inv = 1.0 / (2.0 * step[p]);
		for (size_t i = 0; i < nvol; i++) d[i] = (float)(((double)tmpA[i] - (double)tmpB[i]) * inv);
	}
	if (rc) {
		printfx("-moco: shear factorization failed while building derivatives\n");
		goto done;
	}
	// Normal equations are constant across iterations (derivatives are taken at identity).
	double NE[6][6];
	for (int p = 0; p < 6; p++)
		for (int q = p; q < 6; q++) {
			const float *dp = deriv + (size_t)p * nvol, *dq = deriv + (size_t)q * nvol;
			double s = 0.0;
			for (size_t i = 0; i < nvol; i++) {
				double w = (double)wt[i];
				if (!(w > 0.0)) continue;
				s += w * dp[i] * dq[i];
			}
			NE[p][q] = NE[q][p] = s;
		}
	/* A base drawn from the series is copied through unchanged (and keeps its all-zero parameter
	   row).  With an external reference there is no such sub-brick: the loop below writes every
	   volume, so nothing is copied here. */
	if (base_idx >= 0)
		memcpy(out + (size_t)base_idx * nvol, base, nvol * sizeof(float));

	/* A worker that cannot allocate MUST NOT leave its output volume unwritten: `out` would
	   ship whatever was in it and -1Dfile would report "no motion" for that frame, at exit 0. */
	int worker_oom = 0, worker_fail = 0, nfit_failed = 0, first_failed = -1;
#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic) if (nt > 2)
#endif
	for (int t = 0; t < nt; t++) {
		if (t == base_idx) continue;
		float *lpad = (float *)malloc(nb_pad);
		float *lrow = (float *)malloc(nb_row);
		float *lw = (float *)malloc(nb_vol);
		if (!lpad || !lrow || !lw) {
			free(lpad); free(lrow); free(lw);
#ifdef _OPENMP
			#pragma omp atomic write
#endif
			worker_oom = 1;
			continue;
		}
		const float *mov = img + (size_t)t * nvol;
		/* Repeated linearization about IDENTITY: at every iteration the moving volume is
		   resampled through the CURRENT total estimate and compared with the untouched base, so
		   the derivative basis (taken once, at identity) is the correct basis at the point we
		   linearize about, and the fixed point is the true optimum.  Warping the base forward
		   instead leaves a stale basis and converges ~7% short on every parameter. */
		mo33 Rtot, Rbest;
		double stot[3] = {0, 0, 0}, sbest[3] = {0, 0, 0};
		m_ident(Rtot);
		m_ident(Rbest);
		/* Gauss-Newton with a fixed basis can take a step that makes things worse when the
		   linear model is a poor fit (small or near-degenerate volumes).  Keep the best estimate
		   seen and abandon the iteration the first time the weighted cost rises, rather than
		   letting it run away -- an unguarded loop reached 806 mm of "translation" on a 40 mm
		   deep synthetic volume while still exiting 0. */
		int fitted = 0;            /* did any iteration produce a usable step? */
		double best_cost = -1.0;   /* minimum cost seen, with Rbest/sbest its transform */
		double prev_cost = -1.0;   /* previous iteration's cost, for the rise test */
		for (int it = 0; it < MOCO_MAXITE; it++) {
			mo33 Rinv;
			double sinv[3];
			moco_invert(Rtot, stot, Rinv, sinv);
			if (moco_warp_rs(mov, lw, &g, Rinv, sinv, lpad, lrow)) break;
			/* MEASURED: the reference fits a seventh parameter, an intensity SCALE between the
			   base and the moving volume (the leading number on its -verbose "First fit" lines:
			   0.992952 where the measured brightness ratio is 0.9929, 0.862325 where it is
			   0.8679).  Without it, a run whose volumes are dimmer than the base -- bold3 is 13%
			   dimmer -- fits the brightness difference into the motion parameters.  The scale is
			   linear given the transform, so profile it out in closed form (variable projection)
			   rather than carrying a 7th column: for fixed `a` this IS the exact optimum. */
			double sxy = 0.0, sxx = 0.0;
			for (size_t i = 0; i < nvol; i++) {
				/* Skip zero-weight voxels rather than multiplying by 0.0: 0.0 * NaN is NaN, so a
				   single non-finite voxel anywhere -- including the >=21% that the edging border
				   and MOCO_WCLIP have already excluded -- would poison the whole sum and silently
				   abandon the fit for every volume. */
				double w = (double)wt[i];
				if (!(w > 0.0)) continue;
				double m = (double)lw[i];
				sxy += w * (double)base[i] * m;
				sxx += w * m * m;
			}
			/* No information in the weighted region (an all-zero or dropped frame, or a transform
			   that pushed the data out of the FOV): substituting gfac = 1 makes every iteration
			   identical, so neither the divergence guard nor the convergence test can fire and 23
			   identical steps compose into fabricated motion (measured: 16.6 deg of roll for a
			   volume containing no data).  Also reject an absurd or non-positive scale, which is
			   what a contrast-inverted or dropout frame produces. */
			if (!(sxx > 1e-30 * (1.0 + fabs(sxy)))) break;
			double gfac = sxy / sxx;
			if (!(gfac > 0.125 && gfac < 8.0)) break;
			double cost = 0.0;
			for (size_t i = 0; i < nvol; i++) {
				double w = (double)wt[i];
				if (!(w > 0.0)) continue;
				double r = (double)base[i] - gfac * (double)lw[i];
				cost += w * r * r;
			}
			if (!(cost >= 0.0)) break;                       /* non-finite: keep the last good */
			/* Only a CLEAR rise counts as divergence.  Near convergence the cost routinely
			   ticks up slightly without the step being bad, and reverting there costs real
			   accuracy (measured: a 1e-6 tolerance took max rotation error from 0.0086 to
			   0.0703 deg on the reference run, while MOCO_COST_RISE bounds the pathological
			   case just as well). */
			if (prev_cost >= 0.0 && cost > prev_cost * (1.0 + MOCO_COST_RISE)) {
				/* Diverging: fall back to the MINIMUM-cost transform seen, which is not
				   necessarily the immediately preceding one -- several sub-threshold rises can
				   accumulate before a rejection. */
				memcpy(Rtot, Rbest, sizeof(mo33));
				memcpy(stot, sbest, sizeof(stot));
				break;
			}
			prev_cost = cost;
			if (best_cost < 0.0 || cost < best_cost) {
				best_cost = cost;
				memcpy(Rbest, Rtot, sizeof(mo33));
				memcpy(sbest, stot, sizeof(stot));
			}
			double b[6];
			for (int p = 0; p < 6; p++) {
				const float *dp = deriv + (size_t)p * nvol;
				double sum = 0.0;
				for (size_t i = 0; i < nvol; i++) {
					double w = (double)wt[i];
					if (!(w > 0.0)) continue;
					sum += w * dp[i] * ((double)base[i] - gfac * (double)lw[i]);
				}
				b[p] = -sum;
			}
			double dx[6];
			if (chol6(NE, b, dx)) break;
			int wild = 0;
			for (int p = 0; p < 6; p++)
				if (!(fabs(dx[p]) < 1.0e6)) wild = 1;   /* non-finite or divergent */
			if (wild) break;
			fitted = 1;
			mo33 Rd;
			moco_rot(dx[0], dx[1], dx[2], Rd);
			double sd[3] = {dx[4], dx[5], dx[3]};
			moco_compose(Rd, sd, Rtot, stot, Rtot, stot);
			double dr = fabs(dx[0]) > fabs(dx[1]) ? fabs(dx[0]) : fabs(dx[1]);
			if (fabs(dx[2]) > dr) dr = fabs(dx[2]);
			/* dx[3..5] are (dS, dL, dP) = DICOM (z, x, y); the storage axes may also be
			   permuted, so convert the step to VOXELS through mm2idx rather than dividing by
			   vmm[p-3], which silently pairs dS with the x voxel size. */
			double dstep_mm[3] = {dx[4], dx[5], dx[3]};
			double dstep_vox[3];
			m_vec(g.mm2idx, dstep_mm, dstep_vox);
			double dtv = 0.0;
			for (int p = 0; p < 3; p++) {
				double v = fabs(dstep_vox[p]);
				if (v > dtv) dtv = v;
			}
			if (dtv < MOCO_XTHRESH && dr < MOCO_RTHRESH) break;
		}
		/* MEASURED: do NOT substitute the minimum-cost transform here.  Cost is evaluated at
		   the START of each iteration, so best_cost lags one step behind, and at convergence the
		   final (unscored) step is genuinely better -- selecting the minimum instead moved
		   rotation p95 from 0.0067 to 0.0099 deg and max from 0.0086 to 0.0703 against the
		   oracle.  Rbest/sbest are a DIVERGENCE BACKSTOP, not a minimum-cost selector: they are
		   restored only when the cost rises by more than MOCO_COST_RISE. */
		if (!fitted) {
			/* Both updates are shared state written from every worker.  The counter is a plain
			   increment, but `first_failed` is a read-compare-write that `omp atomic` cannot
			   express, so it needs a critical section -- an unguarded minimum here is a data
			   race that can report the wrong volume number in the diagnostic below. */
#ifdef _OPENMP
			#pragma omp atomic update
#endif
			nfit_failed++;
#ifdef _OPENMP
			#pragma omp critical(moco_first_failed)
#endif
			{
				if (first_failed < 0 || t < first_failed) first_failed = t;
			}
		}
		{
			mo33 Rinv;
			double sinv[3];
			moco_invert(Rtot, stot, Rinv, sinv);
			if (moco_warp_rs(mov, out + (size_t)t * nvol, &g, Rinv, sinv, lpad, lrow) == 0) {
				float lo = mov[0], hi = mov[0];
				for (size_t i = 1; i < nvol; i++) { if (mov[i] < lo) lo = mov[i]; if (mov[i] > hi) hi = mov[i]; }
				float *o = out + (size_t)t * nvol;
				for (size_t i = 0; i < nvol; i++) { if (o[i] < lo) o[i] = lo; if (o[i] > hi) o[i] = hi; }
			} else {
				/* Could not resample this volume: that is an internal failure, not a property
				   of the data, so fail the whole run rather than shipping a plausible image. */
				memcpy(out + (size_t)t * nvol, mov, nvol * sizeof(float));
				m_ident(Rinv);
				sinv[0] = sinv[1] = sinv[2] = 0.0;
#ifdef _OPENMP
				#pragma omp atomic write
#endif
				worker_fail = 1;
			}
			/* Measured (manifest 3.2): the file records the CORRECTION that maps this
			   sub-brick back onto the base, i.e. the INVERSE of the fitted motion.  Extract it
			   from the inverse transform rather than negating the parameters, which would only
			   agree to first order. */
			double roll, pitch, yaw;
			moco_unrot(Rinv, &roll, &pitch, &yaw);
			double *pr = par + (size_t)t * 6;
			pr[0] = roll; pr[1] = pitch; pr[2] = yaw;
			pr[3] = sinv[2]; pr[4] = sinv[0]; pr[5] = sinv[1];   /* dS, dL, dP */
		}
		free(lpad); free(lrow); free(lw);
	}

	if (worker_oom) {
		printfx("-moco: out of memory registering one or more volumes\n");
		rc = 1;
		goto done;
	}
	if (worker_fail) {
		printfx("-moco: could not resample one or more volumes (shear factorization failed)\n");
		rc = 1;
		goto done;
	}
	if (nfit_failed > 0) {
		/* Name a frame: a bare count cannot be told apart from genuinely motionless volumes in
		   the .1D, and a pipeline regressing those columns out would use wrong nuisance terms. */
		printfx("-moco: warning - the fit produced no step for %d of %lld volumes (first: volume "
		        "%d); they are reported as zero motion and passed through uncorrected. A "
		        "non-finite voxel or an empty frame is the usual cause.\n",
		        nfit_failed, (long long)(nt - (base_idx >= 0 ? 1 : 0)), first_failed);
	}
	if (par_path) {
		/* Write through a sibling temporary and rename, so a failed run cannot truncate or
		   delete a parameter file the user already had at this path. */
		size_t plen = strlen(par_path);
		char *tmpname = (char *)malloc(plen + 32);   /* ".mocotmp" + pid digits + NUL */
		if (!tmpname) {
			printfx("-moco: out of memory\n");
			rc = 1;
			goto done;
		}
		/* O_EXCL never follows a pre-existing symlink and never truncates another writer's file.
		   The name is exclusive rather than globally unique, so a crashed run plus PID reuse can
		   leave a stale sibling; try a few suffixes before giving up. */
		int fd = -1;
		for (int attempt = 0; attempt < 8 && fd < 0; attempt++) {
			snprintf(tmpname, plen + 32, "%s.mocotmp%ld_%d", par_path, (long)getpid(), attempt);
			/* O_BINARY (Windows only) keeps the .1D byte-identical across platforms: without it
			   the CRT translates every '\n' this file writes into "\r\n". */
#ifdef _WIN32
			fd = open(tmpname, O_WRONLY | O_CREAT | O_EXCL | O_BINARY, 0600);
#else
			fd = open(tmpname, O_WRONLY | O_CREAT | O_EXCL, 0600);
#endif
		}
		FILE *f = (fd < 0) ? NULL : fdopen(fd, "wb");
		if (!f) {
			if (fd >= 0) { close(fd); remove(tmpname); }
			printfx("-moco: cannot write '%s'\n", par_path);
			free(tmpname);
			rc = 1;
			goto done;
		}
		int bad = 0;
		for (int t = 0; t < nt; t++) {
			const double *p = par + (size_t)t * 6;
			if (fprintf(f, "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
			            p[0], p[1], p[2], p[3], p[4], p[5]) < 0) { bad = 1; break; }
		}
		if (fclose(f) != 0) bad = 1;
#ifdef _WIN32
		/* The Windows CRT's rename() fails when the destination exists, which would break every
		   re-run that overwrites an existing -1Dfile. */
		if (!bad) remove(par_path);
#endif
		if (!bad && rename(tmpname, par_path) != 0) bad = 1;
		if (bad) {
			remove(tmpname);
			printfx("-moco: failed writing '%s'\n", par_path);
			free(tmpname);
			rc = 1;
			goto done;
		}
		free(tmpname);
	}
	free(nim->data);
	nim->data = out;
	out = NULL;

done:
	free(deriv); free(pad); free(rowbuf); free(tmpA); free(tmpB);
	free(wt); free(out); free(par);
	nifti_image_free(ref);
	return rc;
}

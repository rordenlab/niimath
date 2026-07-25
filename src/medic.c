// medic.c - MEDIC (Multi-Echo DIstortion Correction) for niimath
//
// Clean-room emulation of the workflow in Van et al., Imaging Neuroscience 4 (2026),
// doi:10.1162/IMAG.a.1262.  No Warpkit implementation, test, build product or debug symbol was
// read.  Conventions the paper does not fix were measured through the public executables; every
// one of them is recorded, with its experiment, in test/medic_reference_manifest.md.  Section
// numbers in the comments below refer to that manifest.
//
// Pipeline (manifest §3.12):
//
//   per frame : readphase rescale -> MCPC-3D-S phase offset -> multi-echo ROMEO unwrap
//               -> magnitude-weighted regression                     [raw native field, Hz]
//   over time : temporal 2*pi correction -> rank-10 truncation       [_fieldmaps_native]
//               -> scalar fixed-point inversion                      [_fieldmaps]
//               -> * -TRT * pixdim_PE                                [_displacementmaps]
//
// Phase unwrapping is romeo.c's in-memory frame API.  MEDIC owns no unwrapping code.
//
// FP policy: this translation unit follows the project default (-ffast-math).  Only romeo.c needs
// strict FP, because only its 8-bit integer edge weights are reassociation-sensitive; the
// regression, SVD and resampling here are ordinary numerics.  See AGENTS.md.

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _MSC_VER
	#include <unistd.h>   /* getpid() for the temporary output prefix */
#else
	#include <process.h>
	#define getpid _getpid
#endif
#ifdef _OPENMP
	#include <omp.h>
#endif

#include "medic.h"
#include "romeo.h"

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif
#define MD_2PI 6.283185307179586476925286766559

#define MD_ERR(...) do { fprintf(stderr, "** --medic: "); fprintf(stderr, __VA_ARGS__); } while (0)
#define MD_ERRW(...) do { fprintf(stderr, "** -unwarp: "); fprintf(stderr, __VA_ARGS__); } while (0)

#define MD_MAX_ECHO 64
#define MD_RANK_DEFAULT 10
#define MD_INVERT_ITERS 64      /* fixed point; converged well before this (manifest §3.4) */
#define MD_INVERT_TOL 1e-6f     /* Hz; early exit */
#define MD_CORR_THRESH 0.98     /* paper §2.1.3 magnitude-correlation grouping */
/* Guards against undefined float->int conversion.  Any displacement or sample position outside
   these bounds (including NaN and +-Inf, which fail the comparisons) is treated as out of FOV. */
#define MD_DISP_LIMIT 1.0e9
#define MD_POS_LIMIT 1.0e9

/* ============================== small helpers ============================== */

static float md_wrapf(double x) {
	/* principal value in (-pi, pi]; remainder() is exact and avoids the catastrophic
	   cancellation a naive fmod-based wrap shows for large |x|. */
	return (float)remainder(x, MD_2PI);
}

static int md_axis_index(const char *s, int *sign_suffix) {
	/* i/x -> 0, j/y -> 1, k/z -> 2, with the polarity reported through `sign_suffix`.
	   The two callers treat that polarity DIFFERENTLY, and both are measured:
	     --medic  HONOURS it -- the reference returns an identical native field for j and j- but a
	              near-negated displacement map (corr -0.918), so the sign drives the inversion.
	     -unwarp  IGNORES it -- by then the sign already lives in the stored map (manifest 3.5),
	              and applying it again would double-correct. */
	int idx = -1;
	size_t n;
	if (sign_suffix) *sign_suffix = 1;   /* defined on EVERY return, including the error ones */
	if (!s) return -1;
	n = strlen(s);
	if (n < 1 || n > 2) return -1;
	if (n == 2 && s[1] != '-' && s[1] != '+') return -1;
	switch (tolower((unsigned char)s[0])) {
		case 'i': case 'x': idx = 0; break;
		case 'j': case 'y': idx = 1; break;
		case 'k': case 'z': idx = 2; break;
		default: return -1;
	}
	if (sign_suffix) *sign_suffix = (n == 2 && s[1] == '-') ? -1 : 1;
	return idx;
}

/* Voxel->world (RAS) 3x3 of an image.
 *
 * Precedence MUST match core.c's xform(), because md_same_grid() validates images through
 * max_displacement_mm() -- which uses xform() -- while md_offset_per_mm() resamples with the
 * matrix returned here.  A different rule in the two places means the grid check can pass on one
 * matrix while the physics runs on another: an image with sform_code=1, qform_code=2 and
 * disagreeing sform/qform was accepted as "same grid" and then corrected in the WRONG DIRECTION.
 * xform()'s rule is: sform, unless sform_code < qform_code, then qform; pixdim if both unknown. */
static void md_xform3(const nifti_image *nim, double A[3][3]) {
	nifti_dmat44 m;
	int r, c;
	if (nim->sform_code == NIFTI_XFORM_UNKNOWN && nim->qform_code == NIFTI_XFORM_UNKNOWN) {
		memset(A, 0, 9 * sizeof(double));
		A[0][0] = nim->dx != 0.0 ? nim->dx : 1.0;
		A[1][1] = nim->dy != 0.0 ? nim->dy : 1.0;
		A[2][2] = nim->dz != 0.0 ? nim->dz : 1.0;
		return;
	}
	m = (nim->sform_code < nim->qform_code) ? nim->qto_xyz : nim->sto_xyz;
	for (r = 0; r < 3; r++) for (c = 0; c < 3; c++) A[r][c] = m.m[r][c];
}

static int md_inv3(const double A[3][3], double Inv[3][3]) {
	double d = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
		- A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
		+ A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
	double id;
	if (!(d > 1e-12 || d < -1e-12)) return 1;
	id = 1.0 / d;
	Inv[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * id;
	Inv[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * id;
	Inv[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * id;
	Inv[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * id;
	Inv[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * id;
	Inv[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * id;
	Inv[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * id;
	Inv[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * id;
	Inv[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * id;
	return 0;
}

/* Per-millimetre voxel-space offset for a scalar displacement map along voxel axis `m`.
 *
 * MEASURED convention (manifest §3.5), not a guess: the displacement is a PHYSICAL vector along
 * the canonical world axis that voxel axis `m` is most aligned with -- not along the image's own
 * column direction.  With u = A[:,m]/|A[:,m]|, w = argmax|u|, sigma = sign(u[w]) and
 * kappa = (-1,-1,+1) (the RAS->LPS sign of world axis w):
 *
 *     delta_RAS = d * sigma * kappa[w] * e_w
 *     s         = A^-1 @ delta_RAS          and   out(v) = in(v + s)
 *
 * Verified on 30 synthetic grid/letter combinations exactly, and on the real oblique demo data at
 * nrmse 3.5e-5 (versus 4.3e-2 for the image's own column direction).
 *
 * The map is in MILLIMETRES by contract, so A is normalised to millimetres via xyz_units_to_mm()
 * before inversion.  A mm-unit header (every real neuroimaging NIfTI) scales by exactly 1.0 and is
 * bit-identical to the un-normalised form; a metre- or micron-unit header would otherwise be off
 * by 1000x. */
static int md_offset_per_mm(const nifti_image *nim, int m, double s_per_mm[3]) {
	double A[3][3], Inv[3][3], u[3], nrm = 0.0, kappa[3] = { -1.0, -1.0, 1.0 };
	double delta[3] = { 0.0, 0.0, 0.0 };
	double unit = xyz_units_to_mm(nim->xyz_units);
	int r, c, w = 0;
	md_xform3(nim, A);
	if (unit != 1.0) for (r = 0; r < 3; r++) for (c = 0; c < 3; c++) A[r][c] *= unit;
	for (r = 0; r < 3; r++) nrm += A[r][m] * A[r][m];
	nrm = sqrt(nrm);
	if (!(nrm > 1e-12)) return 1;
	for (r = 0; r < 3; r++) u[r] = A[r][m] / nrm;
	for (r = 1; r < 3; r++) if (fabs(u[r]) > fabs(u[w])) w = r;
	delta[w] = (u[w] >= 0.0 ? 1.0 : -1.0) * kappa[w];
	if (md_inv3(A, Inv)) return 1;
	for (r = 0; r < 3; r++) s_per_mm[r] = Inv[r][0] * delta[0] + Inv[r][1] * delta[1] + Inv[r][2] * delta[2];
	return 0;
}

/* ============================== resampling ============================== */

/* Lanczos-windowed sinc, radius 5, UNNORMALIZED -- measured to 1.3e-8 against the reference
   impulse response (manifest §3.6).  The weights sum to 0.998746 at a half-voxel offset; that
   0.13 % dip is part of the convention, so do NOT normalise it. */
#define MD_LANCZOS_R 5
static double md_lanczos(double t) {
	double pt, pr;
	if (t <= -(double)MD_LANCZOS_R || t >= (double)MD_LANCZOS_R) return 0.0;
	if (t > -1e-12 && t < 1e-12) return 1.0;
	pt = M_PI * t;
	pr = pt / (double)MD_LANCZOS_R;
	return (sin(pt) / pt) * (sin(pr) / pr);
}

/* out[v] = in[v + s(v)] with a separable 3D Lanczos-5 kernel and ZERO fill outside the FOV.
   `disp` is the scalar map in header length units; `s_per_mm` converts it to a voxel offset.
   in/out are one 3D volume each and must not alias. */
static void md_pull(const float *in, float *out, int nx, int ny, int nz,
	const float *disp, const double s_per_mm[3]) {
	const int64_t nxy = (int64_t)nx * ny;
	int z;
#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for (z = 0; z < nz; z++) {
		int x, y;
		for (y = 0; y < ny; y++) {
			for (x = 0; x < nx; x++) {
				int64_t o = (int64_t)x + (int64_t)y * nx + (int64_t)z * nxy;
				double d = (double)disp[o];
				double px, py, pz;
				double wx[2 * MD_LANCZOS_R], wy[2 * MD_LANCZOS_R], wz[2 * MD_LANCZOS_R];
				int bx, by, bz;
				double acc = 0.0;
				int tx, ty, tz;
				/* A non-finite map value would make floor() non-finite and the cast to int
				   UNDEFINED. Treat such a voxel as fully out of FOV (the documented fill). */
				if (!(d >= -MD_DISP_LIMIT && d <= MD_DISP_LIMIT)) { out[o] = 0.0f; continue; }
				px = x + d * s_per_mm[0]; py = y + d * s_per_mm[1]; pz = z + d * s_per_mm[2];
				if (!(px >= -MD_POS_LIMIT && px <= MD_POS_LIMIT) ||
					!(py >= -MD_POS_LIMIT && py <= MD_POS_LIMIT) ||
					!(pz >= -MD_POS_LIMIT && pz <= MD_POS_LIMIT)) { out[o] = 0.0f; continue; }
				bx = (int)floor(px); by = (int)floor(py); bz = (int)floor(pz);
				/* Fast path: an exactly-zero displacement is the identity for this kernel
				   (sinc vanishes at every nonzero integer), which keeps unshifted background
				   voxels bit-exact and skips 1000 taps for them. */
				if (d == 0.0f && px == (double)x && py == (double)y && pz == (double)z) {
					out[o] = in[o];
					continue;
				}
				for (tx = 0; tx < 2 * MD_LANCZOS_R; tx++) wx[tx] = md_lanczos(px - (bx + tx - (MD_LANCZOS_R - 1)));
				for (ty = 0; ty < 2 * MD_LANCZOS_R; ty++) wy[ty] = md_lanczos(py - (by + ty - (MD_LANCZOS_R - 1)));
				for (tz = 0; tz < 2 * MD_LANCZOS_R; tz++) wz[tz] = md_lanczos(pz - (bz + tz - (MD_LANCZOS_R - 1)));
				for (tz = 0; tz < 2 * MD_LANCZOS_R; tz++) {
					int iz = bz + tz - (MD_LANCZOS_R - 1);
					double az;
					if (iz < 0 || iz >= nz || wz[tz] == 0.0) continue;
					az = 0.0;
					for (ty = 0; ty < 2 * MD_LANCZOS_R; ty++) {
						int iy = by + ty - (MD_LANCZOS_R - 1);
						double ay = 0.0;
						const float *row;
						if (iy < 0 || iy >= ny || wy[ty] == 0.0) continue;
						row = in + (int64_t)iy * nx + (int64_t)iz * nxy;
						for (tx = 0; tx < 2 * MD_LANCZOS_R; tx++) {
							int ix = bx + tx - (MD_LANCZOS_R - 1);
							if (ix < 0 || ix >= nx) continue;
							ay += wx[tx] * (double)row[ix];
						}
						az += wy[ty] * ay;
					}
					acc += wz[tz] * az;
				}
				out[o] = (float)acc;
			}
		}
	}
}

/* ============================== I/O helpers ============================== */

/* Read a NIfTI as float32.  Returns the image (caller frees with nifti_image_free) with
   ->data already converted, or NULL. */
static nifti_image *md_read_f32(const char *fn, const char *what) {
	nifti_image *n;
	in_hdr ihdr;
	{	/* Header-only preflight: reject an oversized or malformed image BEFORE decompressing
		   and allocating its payload, rather than after. */
		nifti_image *h = nifti_image_read(fn, 0);
		int bad = 0;
		if (!h) { MD_ERR("failed to read the header of %s '%s'\n", what, fn); return NULL; }
		if (h->nvox < 1 || h->nx < 1 || h->ny < 1 || h->nz < 1) {
			MD_ERR("%s '%s' has invalid dimensions\n", what, fn); bad = 1;
		} else if (h->nu > 1 || h->nv > 1 || h->nw > 1) {
			MD_ERR("%s '%s' has more than 4 dimensions (5D input is out of scope)\n", what, fn); bad = 1;
		} else if ((int64_t)h->nvox > INT_MAX) {
			MD_ERR("%s '%s' exceeds INT_MAX voxels; --medic is not a huge-image-safe operation\n", what, fn);
			bad = 1;
		}
		nifti_image_free(h);
		if (bad) return NULL;
	}
	n = nifti_image_read(fn, 1);
	if (!n) { MD_ERR("failed to read %s '%s'\n", what, fn); return NULL; }
	/* Re-check after the load as well as before it.  The preflight is what stops us decompressing
	   a huge payload; this is the fail-closed guarantee, and it costs nothing.  Dropping it left a
	   TOCTOU window if the file changed between the two reads, and departed from the project's
	   "re-check at use" pattern (nii_admit_current_op). */
	if (n->nvox < 1 || n->nx < 1 || n->ny < 1 || n->nz < 1 ||
		n->nu > 1 || n->nv > 1 || n->nw > 1 || (int64_t)n->nvox > INT_MAX) {
		MD_ERR("%s '%s' changed on disk or has unusable dimensions\n", what, fn);
		nifti_image_free(n); return NULL;
	}
	ihdr = set_input_hdr(n);
	/* Convert when the stored type is not float32, but ALSO when it IS float32 and carries a
	   non-trivial scl_slope/scl_inter -- otherwise a scaled float32 image is silently used raw.
	   That bit anyone storing a displacement map as float32 with a slope. */
	if (n->datatype != DT_FLOAT32 ||
		(n->scl_slope != 0.0f && n->scl_slope != 1.0f) || n->scl_inter != 0.0f) {
		if (nifti_image_change_datatype(n, DT_FLOAT32, &ihdr) != 0) {
			MD_ERR("failed to convert %s '%s' to float32\n", what, fn);
			nifti_image_free(n); return NULL;
		}
	}
	return n;
}

/* Do two images share a grid?  Dimensions exactly, plus the project's existing world-transform
   metric -- max_displacement_mm() measures true 3D corner displacement AND normalises each header
   to millimetres via xyz_units, so it covers rotation, scale, origin and units in one call.  This
   is the same 0.001 mm gate --qc uses for its own same-grid requirement. */
static int md_same_grid(nifti_image *a, nifti_image *b) {
	return a->nx == b->nx && a->ny == b->ny && a->nz == b->nz &&
		max_displacement_mm(a, b) <= 0.001f;
}

/* ============================== -unwarp ============================== */

int medic_unwarp(nifti_image *nim, const char *mapfile, const char *axis) {
	nifti_image *map = NULL;
	double s_per_mm[3];
	int m, nx, ny, nz, nt, mt, t, rc = 1;
	int64_t n3;
	float *out = NULL;
	const float *in;

	if (!nim || nim->datatype != DT_FLOAT32) { MD_ERRW("internal error: expected a float32 working image\n"); return 1; }
	if (nim->nu > 1 || nim->nv > 1 || nim->nw > 1) { MD_ERRW("input must be 3D or 4D\n"); return 1; }
	if ((int64_t)nim->nvox > INT_MAX) { MD_ERRW("input exceeds INT_MAX voxels; -unwarp is not a huge-image-safe operation\n"); return 1; }
	m = md_axis_index(axis, NULL);
	if (m < 0) { MD_ERRW("axis must be one of i j k x y z (a trailing '-' is accepted and ignored)\n"); return 1; }

	nx = nim->nx; ny = nim->ny; nz = (nim->nz < 1 ? 1 : nim->nz);
	n3 = (int64_t)nx * ny * nz;
	if (n3 < 1 || (int64_t)nim->nvox % n3 != 0) { MD_ERRW("invalid image geometry\n"); return 1; }
	nt = (int)((int64_t)nim->nvox / n3);

	map = md_read_f32(mapfile, "displacement map");
	if (!map) return 1;
	if (!md_same_grid(nim, map)) {
		MD_ERRW("displacement map '%s' does not share the input's grid (dimensions and world transform must match)\n", mapfile);
		goto done;
	}
	mt = (int)((int64_t)map->nvox / n3);
	if (mt != 1 && mt != nt) {
		MD_ERRW("displacement map has %d frame(s); expected 1 (broadcast) or %d (one per input frame)\n", mt, nt);
		goto done;
	}
	if (md_offset_per_mm(map, m, s_per_mm)) { MD_ERRW("displacement map has a singular world transform\n"); goto done; }

	out = (float *)nii_malloc((size_t)nim->nvox, sizeof(float));
	in = (const float *)nim->data;
	for (t = 0; t < nt; t++)
		md_pull(in + (int64_t)t * n3, out + (int64_t)t * n3, nx, ny, nz,
			((const float *)map->data) + (int64_t)(mt == 1 ? 0 : t) * n3, s_per_mm);
	memcpy(nim->data, out, (size_t)nim->nvox * sizeof(float));
	rc = 0;
done:
	free(out);
	nifti_image_free(map);
	return rc;
}

/* ============================== symmetric eigensolver ============================== */

/* Cyclic Jacobi for a small dense symmetric matrix.  T is the number of frames, so this is at
 * most a few hundred on real runs and the O(T^3) cost is invisible next to the unwrapping.
 *
 * ponytail: a T x T Jacobi rotation sweep replaces the plan's block-Lanczos / subspace-iteration
 * solver.  It is ~40 lines, deterministic, needs no convergence tuning, and returns the FULL
 * spectrum so "rank = min(10, T, positive numerical rank)" is a truncation rather than an
 * iterative target.  Upgrade path: if T ever reaches the thousands, swap in a Lanczos solver --
 * the caller only needs the leading k eigenpairs.
 *
 * a[] is T*T row-major and is destroyed; v[] receives the eigenvectors as columns; w[] the
 * eigenvalues.  Returns 0 on success. */
static void md_jacobi_eigh(double *a, double *v, double *w, int n) {
	int i, j, p, q, sweep;
	for (i = 0; i < n; i++) for (j = 0; j < n; j++) v[(size_t)i * n + j] = (i == j) ? 1.0 : 0.0;
	for (sweep = 0; sweep < 100; sweep++) {
		double off = 0.0;
		for (p = 0; p < n; p++) for (q = p + 1; q < n; q++) off += a[(size_t)p * n + q] * a[(size_t)p * n + q];
		if (off <= 1e-30) break;
		for (p = 0; p < n - 1; p++) for (q = p + 1; q < n; q++) {
			double apq = a[(size_t)p * n + q], app, aqq, theta, t, c, s;
			if (fabs(apq) < 1e-300) continue;
			app = a[(size_t)p * n + p];
			aqq = a[(size_t)q * n + q];
			theta = (aqq - app) / (2.0 * apq);
			t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
			c = 1.0 / sqrt(t * t + 1.0);
			s = t * c;
			for (i = 0; i < n; i++) {
				double aip = a[(size_t)i * n + p], aiq = a[(size_t)i * n + q];
				a[(size_t)i * n + p] = c * aip - s * aiq;
				a[(size_t)i * n + q] = s * aip + c * aiq;
			}
			for (i = 0; i < n; i++) {
				double api = a[(size_t)p * n + i], aqi = a[(size_t)q * n + i];
				a[(size_t)p * n + i] = c * api - s * aqi;
				a[(size_t)q * n + i] = s * api + c * aqi;
			}
			for (i = 0; i < n; i++) {
				double vip = v[(size_t)i * n + p], viq = v[(size_t)i * n + q];
				v[(size_t)i * n + p] = c * vip - s * viq;
				v[(size_t)i * n + q] = s * vip + c * viq;
			}
		}
	}
	for (i = 0; i < n; i++) w[i] = a[(size_t)i * n + i];
	/* sort descending by eigenvalue, permuting the eigenvector columns with them */
	for (i = 0; i < n - 1; i++) {
		int best = i;
		for (j = i + 1; j < n; j++) if (w[j] > w[best]) best = j;
		if (best == i) continue;
		{
			double tmp = w[i]; w[i] = w[best]; w[best] = tmp;
			for (p = 0; p < n; p++) {
				tmp = v[(size_t)p * n + i]; v[(size_t)p * n + i] = v[(size_t)p * n + best]; v[(size_t)p * n + best] = tmp;
			}
		}
	}
}

/* Rank-`rank` truncation of the Nvox x T matrix F, in place, UNCENTERED (manifest §3.8).
 *
 * Uses the T x T Gram matrix G = F^T F, whose eigenvectors are the right singular vectors of F.
 * Projecting each voxel's time course onto the leading k of them is exactly the truncated SVD:
 * F_k = F V_k V_k^T.  Memory is O(T^2), independent of the voxel count. */
static int md_lowrank(float *F, int64_t nvox, int T, int rank) {
	double *G = NULL, *V = NULL, *w = NULL, *P = NULL;
	int i, j, k, r, rc = 1;
	int64_t v;
	if (rank <= 0 || T <= 1 || rank >= T) return 0;   /* nothing to truncate */
	G = (double *)calloc((size_t)T * T, sizeof(double));
	V = (double *)malloc((size_t)T * T * sizeof(double));
	w = (double *)malloc((size_t)T * sizeof(double));
	P = (double *)malloc((size_t)T * T * sizeof(double));
	if (!G || !V || !w || !P) { MD_ERR("out of memory in the low-rank filter\n"); goto done; }

	/* G = F^T F, accumulated in double.  One pass over the series, frame-major access. */
	for (i = 0; i < T; i++) {
		for (j = i; j < T; j++) {
			double s = 0.0;
			const float *a = F + (int64_t)i * nvox, *b = F + (int64_t)j * nvox;
			for (v = 0; v < nvox; v++) s += (double)a[v] * (double)b[v];
			G[(size_t)i * T + j] = G[(size_t)j * T + i] = s;
		}
	}
	md_jacobi_eigh(G, V, w, T);

	/* Numerical rank: drop directions that are pure round-off relative to the leading one.
	   The series is known finite here (md_all_finite() gates the caller), so r == 0 can only mean
	   a genuinely zero/degenerate spectrum -- e.g. an all-zero field, whose rank-k truncation is
	   itself zero.  That is a NO-OP, not an error; only a non-finite series is an error, and that
	   is caught before this function is reached. */
	r = rank < T ? rank : T;
	for (k = 0; k < r; k++) if (!(w[k] > w[0] * 1e-24) || !(w[k] > 0.0)) { r = k; break; }
	if (r <= 0 || r >= T) { rc = 0; goto done; }   /* nothing to project onto, or nothing to drop */

	/* P = V_r V_r^T (T x T projector) */
	for (i = 0; i < T; i++) for (j = 0; j < T; j++) {
		double s = 0.0;
		for (k = 0; k < r; k++) s += V[(size_t)i * T + k] * V[(size_t)j * T + k];
		P[(size_t)i * T + j] = s;
	}
	{
		int64_t chunk;
		const int64_t CH = 1 << 16;
		float *tmp = NULL;
		int oom = 0;
		{	/* Check the scratch SIZE ARITHMETIC up front.  The malloc itself still happens per
			   thread inside the region, so an OOM there can leave some chunks unfiltered before
			   the reduction reports it -- hence the error text below says so.  Fail-loud, not
			   atomic. */
			size_t bytes;
			if (nii_mul_size((size_t)(CH < nvox ? CH : nvox), (size_t)T, &bytes) ||
				nii_mul_size(bytes, sizeof(float), &bytes)) {
				MD_ERR("low-rank scratch size overflows this build's address space\n");
				goto done;
			}
		}
#ifdef _OPENMP
		#pragma omp parallel private(tmp) reduction(|:oom)
#endif
		{
			tmp = (float *)malloc((size_t)(CH < nvox ? CH : nvox) * T * sizeof(float));
			if (!tmp) oom = 1;
#ifdef _OPENMP
			#pragma omp for schedule(static)
#endif
			for (chunk = 0; chunk < nvox; chunk += CH) {
				int64_t n = (nvox - chunk < CH) ? (nvox - chunk) : CH, q;
				int ii, jj;
				if (!tmp) continue;
				for (ii = 0; ii < T; ii++) memcpy(tmp + (int64_t)ii * n, F + (int64_t)ii * nvox + chunk, (size_t)n * sizeof(float));
				for (ii = 0; ii < T; ii++) {
					float *dst = F + (int64_t)ii * nvox + chunk;
					for (q = 0; q < n; q++) dst[q] = 0.0f;
					for (jj = 0; jj < T; jj++) {
						double p = P[(size_t)ii * T + jj];
						const float *src = tmp + (int64_t)jj * n;
						if (p == 0.0) continue;
						for (q = 0; q < n; q++) dst[q] += (float)(p * (double)src[q]);
					}
				}
			}
			free(tmp);
		}
		if (oom) { MD_ERR("out of memory in the low-rank filter (the field series may be partly filtered)\n"); goto done; }
	}
	rc = 0;
done:
	free(G); free(V); free(w); free(P);
	return rc;
}

/* ============================== MEDIC stages ============================== */

/* Magnitude guards, not isfinite(): this TU is -ffast-math. */
static int md_all_finite(const float *v, int64_t n) {
	int64_t i;
	for (i = 0; i < n; i++) if (!(v[i] >= -FLT_MAX && v[i] <= FLT_MAX)) return 0;
	return 1;
}


typedef struct {
	int nx, ny, nz, neco, nframe;
	int64_t n3;
	double TEs[MD_MAX_ECHO];      /* milliseconds */
	double trt;                   /* seconds */
	int pe_axis;                  /* 0/1/2 */
	int pe_sign;                  /* +1 for i/j/k, -1 for i-/j-/k-  (see md_invert) */
	int rank;
	int temporal;
	int mcpc;
	int noiseframes;
	int save_intermediates;
	const char *maskfile;         /* --mask: use this mask verbatim instead of robustmask */
	const char *prefix;
	nifti_image *tmpl;            /* header template (phase echo 1) */
} md_ctx;

/* readphase: rescale observed [min,max] of the whole series onto [-pi,pi] (manifest §3.1).
   Matches ROMEO's readphase, which is what the reference uses. */
static void md_rescale_phase(float *p, int64_t n) {
	int64_t i;
	double mn = 0.0, mx = 0.0, slope, inter, span;
	int seen = 0;
	for (i = 0; i < n; i++) {
		double v = (double)p[i];
		if (!isfinite(v)) continue;
		if (!seen) { mn = mx = v; seen = 1; }
		else { if (v < mn) mn = v; if (v > mx) mx = v; }
	}
	if (!seen) return;
	span = mx - mn;
	if (fabs(span - MD_2PI) <= 0.1) return;   /* already radians */
	if (!(span > 0.0)) return;
	slope = MD_2PI / span;
	inter = -M_PI - mn * slope;
	for (i = 0; i < n; i++) p[i] = (float)((double)p[i] * slope + inter);
}

/* MCPC-3D-S monopolar phase offset (paper §2.1.2 Eq. 4).
 *
 * MEASURED to be exactly (manifest §4, residual 4.6e-5 rad against a 9.6e-5 rad quantum):
 *
 *   hip       = m1*m2 * exp(i*(phi2 - phi1))          (Hermitian inner product of echoes 1,2)
 *   d_uw      = ROMEO_unwrap(angle(hip), mag=|hip|)
 *   offset    = wrap( phi1 - TE1/(TE2-TE1) * d_uw )
 *
 * with NO spatial smoothing -- the reference's offset matches the unsmoothed expression exactly.
 * Bipolar acquisition and multi-channel combination are out of scope (plan §11).
 *
 * `phase` (n3*neco) is corrected IN PLACE; `offset_out` may be NULL. */
static int md_mcpc3ds(const md_ctx *c, float *phase, const float *mag, const romeo_opts *ro,
	const uint8_t *mask, float *offset_out) {
	const int64_t n3 = c->n3;
	float *hipp = NULL, *hipm = NULL;
	double dTE = c->TEs[1] - c->TEs[0], k;
	int64_t i;
	int e, rc = 1;
	if (c->neco < 2) return 0;
	if (!(fabs(dTE) > 1e-12)) { MD_ERR("echo times 1 and 2 are identical; MCPC-3D-S needs a nonzero echo spacing\n"); return 1; }
	k = c->TEs[0] / dTE;
	hipp = (float *)malloc((size_t)n3 * sizeof(float));
	hipm = (float *)malloc((size_t)n3 * sizeof(float));
	if (!hipp || !hipm) { MD_ERR("out of memory in MCPC-3D-S\n"); goto done; }
	for (i = 0; i < n3; i++) {
		hipp[i] = md_wrapf((double)phase[n3 + i] - (double)phase[i]);
		hipm[i] = (float)((double)mag[i] * (double)mag[n3 + i]);
	}
	{	/* Single-echo spatial unwrap of the phase difference.
		 *
		 * The weight preset comes from the caller (`--weights`, default romeo4) and is applied
		 * here as well as to the multi-echo unwrap -- the reference uses one preset for both, so
		 * overriding it here would make --weights a half-measure.
		 *
		 * Both choices here are MEASURED, not defaults: with the SHARED mask and romeo4 weights
		 * the resulting offset reproduces the reference's own phase_offset EXACTLY (frac exact
		 * 1.0000, p95 0.0000 rad).  romeo3 / a per-HIP robustmask / nomask all leave 11-18 % of
		 * in-mask voxels on a different 2*pi branch.  See manifest section 4. */
		double te1 = fabs(dTE);
		romeo_opts o = *ro;
		o.nTE = 1; o.TEs[0] = te1; o.te_epi = 0; o.template_echo = 1;
		o.individual = 0; o.correctglobal = 0;
		if (romeo_unwrap_frame(hipp, hipm, 1, c->nx, c->ny, c->nz, 1, &te1, &o, mask, NULL)) {
			MD_ERR("ROMEO failed while unwrapping the MCPC-3D-S phase difference\n");
			goto done;
		}
	}
	for (i = 0; i < n3; i++) {
		float off = md_wrapf((double)phase[i] - k * (double)hipp[i]);
		if (offset_out) offset_out[i] = off;
		for (e = 0; e < c->neco; e++)
			phase[(int64_t)e * n3 + i] = md_wrapf((double)phase[(int64_t)e * n3 + i] - (double)off);
	}
	rc = 0;
done:
	free(hipp); free(hipm);
	return rc;
}

/* Magnitude-weighted regression through the origin (paper §3.4; measured exact, manifest §3.2).
   omega = sum(m^2 t phi) / sum(m^2 t^2); field_Hz = omega / 2pi.  t in SECONDS. */
static void md_regress(const md_ctx *c, const float *phase, const float *mag, float *field) {
	const int64_t n3 = c->n3;
	int64_t i;
	int e;
	for (i = 0; i < n3; i++) {
		double num = 0.0, den = 0.0;
		for (e = 0; e < c->neco; e++) {
			double t = c->TEs[e] * 1e-3;
			double m = (double)mag[(int64_t)e * n3 + i];
			double m2 = m * m;
			num += m2 * t * (double)phase[(int64_t)e * n3 + i];
			den += m2 * t * t;
		}
		field[i] = (den > 0.0) ? (float)(num / (den * MD_2PI)) : 0.0f;
	}
}

/* Temporal 2*pi consistency correction (paper §2.1.3, Eqs. 5-6).
 *
 * Frames are grouped by magnitude correlation >= MD_CORR_THRESH on the FIRST echo; each frame's
 * first-echo unwrapped phase is moved to the 2*pi branch nearest its group mean, and every later
 * echo to the branch predicted by ALL previously corrected echoes.
 *
 * That prediction is the paper's Eq. 6: a through-origin fit over the echoes already corrected,
 *
 *     phi_n_predicted = t_n * sum_{i<n}(phi_i * t_i) / sum_{i<n}(t_i^2)
 *
 * which reduces to phi_1 * t_n/t_1 for n = 2 -- so two-echo data cannot distinguish it from the
 * naive "scale echo 1" form, and three-or-more-echo data can.
 *
 * The group mean is taken from an IMMUTABLE SNAPSHOT: correcting in place while reading would make
 * later frames see already-corrected earlier ones and earlier frames not, i.e. a result that
 * depends on traversal order.
 *
 * The reference's observed behaviour -- a spatially uniform whole-2*pi shift per frame per echo --
 * is reproduced by this per-voxel rounding whenever the discrepancy is itself uniform, which is
 * the case that matters (manifest §3.9).  The grouping threshold is the paper's; it is the one
 * parameter here that the black box could not be made to reveal.
 *
 * `uw` is neco * nframe volumes, echo-major within frame.  `mag1` is the first echo's magnitude
 * series (n3 * nframe). */
static int md_temporal(const md_ctx *c, float *uw, const float *mag1) {
	const int64_t n3 = c->n3;
	const int T = c->nframe;
	double *mu = NULL, *sd = NULL, *corr = NULL;
	float *acc = NULL, *snap = NULL;
	int t, u, e, rc = 1;
	int64_t i;
	size_t bytes;
	if (T < 2) return 0;
	if (nii_mul_size((size_t)T, (size_t)T, &bytes) ||
		nii_mul_size(bytes, sizeof(double), &bytes) ||
		nii_mul_size((size_t)n3, (size_t)T, &bytes) ||
		nii_mul_size(bytes, sizeof(float), &bytes)) {
		MD_ERR("frame count %d is too large for the temporal correction on this build\n", T);
		return 1;
	}
	mu = (double *)malloc((size_t)T * sizeof(double));
	sd = (double *)malloc((size_t)T * sizeof(double));
	corr = (double *)malloc((size_t)T * T * sizeof(double));
	acc = (float *)malloc((size_t)n3 * sizeof(float));
	/* snapshot of every frame's FIRST-echo unwrapped phase, so group means are order-independent */
	snap = (float *)malloc((size_t)n3 * T * sizeof(float));
	if (!mu || !sd || !corr || !acc || !snap) { MD_ERR("out of memory in the temporal correction\n"); goto done; }
	for (t = 0; t < T; t++)
		memcpy(snap + (int64_t)t * n3, uw + ((int64_t)t * c->neco) * n3, (size_t)n3 * sizeof(float));

	for (t = 0; t < T; t++) {
		const float *m = mag1 + (int64_t)t * n3;
		double s = 0.0, s2 = 0.0;
		for (i = 0; i < n3; i++) s += (double)m[i];
		mu[t] = s / (double)n3;
		for (i = 0; i < n3; i++) { double d = (double)m[i] - mu[t]; s2 += d * d; }
		sd[t] = sqrt(s2);
	}
#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic)
#endif
	for (t = 0; t < T; t++) {
		int v;
		for (v = 0; v < T; v++) {
			const float *a = mag1 + (int64_t)t * n3, *b = mag1 + (int64_t)v * n3;
			double s = 0.0;
			int64_t q;
			if (v < t) { corr[(size_t)t * T + v] = 0.0; continue; }   /* filled by symmetry below */
			for (q = 0; q < n3; q++) s += ((double)a[q] - mu[t]) * ((double)b[q] - mu[v]);
			corr[(size_t)t * T + v] = (sd[t] > 0.0 && sd[v] > 0.0) ? s / (sd[t] * sd[v]) : 0.0;
		}
	}
	for (t = 0; t < T; t++) for (u = 0; u < t; u++) corr[(size_t)t * T + u] = corr[(size_t)u * T + t];

	for (t = 0; t < T; t++) {
		int ng = 0;
		for (u = 0; u < T; u++) if (corr[(size_t)t * T + u] >= MD_CORR_THRESH) ng++;
		if (ng < 2) continue;   /* a frame alone in its group has no reference to move toward */
		for (i = 0; i < n3; i++) acc[i] = 0.0f;
		for (u = 0; u < T; u++) {
			const float *p;
			if (corr[(size_t)t * T + u] < MD_CORR_THRESH) continue;
			p = snap + (int64_t)u * n3;   /* snapshot, not the live (partly corrected) series */
			for (i = 0; i < n3; i++) acc[i] += p[i];
		}
		{
			float *p1 = uw + ((int64_t)t * c->neco) * n3;
			for (i = 0; i < n3; i++) {
				double ref = (double)acc[i] / (double)ng;
				double n = nearbyint((ref - (double)p1[i]) / MD_2PI);
				p1[i] = (float)((double)p1[i] + MD_2PI * n);
			}
			for (e = 1; e < c->neco; e++) {
				/* Eq. 6: through-origin fit over the echoes already corrected (0..e-1).
				   denom is voxel-independent, so hoist it out of the voxel loop. */
				float *pe = uw + ((int64_t)t * c->neco + e) * n3;
				double tn = c->TEs[e], denom = 0.0;
				int k;
				for (k = 0; k < e; k++) denom += c->TEs[k] * c->TEs[k];
				if (!(denom > 0.0)) continue;
				for (i = 0; i < n3; i++) {
					double num = 0.0, pred, n;
					for (k = 0; k < e; k++)
						num += (double)uw[((int64_t)t * c->neco + k) * n3 + i] * c->TEs[k];
					pred = tn * num / denom;
					n = nearbyint((pred - (double)pe[i]) / MD_2PI);
					pe[i] = (float)((double)pe[i] + MD_2PI * n);
				}
			}
		}
	}
	rc = 0;
done:
	free(mu); free(sd); free(corr); free(acc); free(snap);
	return rc;
}

/* Scalar displacement inversion along the PE VOXEL axis (manifest §3.4):
 *
 *     f_undistorted(y) = f_native( y + s * f_undistorted(y) * TRT )     [voxels]
 *
 * solved by direct iteration from zero with LINEAR interpolation and edge clamping -- measured to
 * be the reference's sampling (linear p95 0.041 Hz; cubic 2.05; nearest 4.07).  Note the
 * composition runs along the voxel axis, unlike the RESAMPLING in §3.5, which runs along the
 * canonical world axis; both were measured separately.
 *
 * `s` is the phase-encoding POLARITY (+1 for `j`, -1 for `j-`).  It is load-bearing: the reference
 * returns an IDENTICAL native field for j and j-, but an inverted field that differs (corr 0.918)
 * and a displacement map that is very nearly negated (corr -0.918, median ratio -0.973).  Getting
 * it wrong on a `j-` acquisition doubles the distortion instead of correcting it.  Verified for
 * both polarities against the reference at displacement p95 0.045 mm (j) and 0.024 mm (j-). */
static void md_invert(const md_ctx *c, const float *fn, float *fu) {
	const int nx = c->nx, ny = c->ny, nz = c->nz;
	const int m = c->pe_axis;
	const int64_t n3 = c->n3;
	const int64_t stride = (m == 0) ? 1 : ((m == 1) ? nx : (int64_t)nx * ny);
	const int len = (m == 0) ? nx : ((m == 1) ? ny : nz);
	int it;
	int64_t i;
	for (i = 0; i < n3; i++) fu[i] = 0.0f;
	for (it = 0; it < MD_INVERT_ITERS; it++) {
		double worst = 0.0;
		int z;
#ifdef _OPENMP
		#pragma omp parallel for schedule(static) reduction(max:worst)
#endif
		for (z = 0; z < nz; z++) {
			int x, y;
			for (y = 0; y < ny; y++) for (x = 0; x < nx; x++) {
				int64_t o = (int64_t)x + (int64_t)y * nx + (int64_t)z * nx * ny;
				int idx = (m == 0) ? x : ((m == 1) ? y : z);
				int64_t base = o - (int64_t)idx * stride;
				double cur = (double)fu[o];
				double pos;
				int lo;
				double frac, a, b, nv;
				if (!(cur >= -MD_DISP_LIMIT && cur <= MD_DISP_LIMIT)) cur = 0.0; /* NaN/Inf -> no shift */
				pos = (double)idx + (double)c->pe_sign * cur * c->trt;
				if (!(pos >= -MD_POS_LIMIT && pos <= MD_POS_LIMIT)) pos = (double)idx;
				if (len < 2) { fu[o] = fn[o]; continue; }   /* single slice along PE: nothing to interpolate */
				if (pos < 0.0) pos = 0.0;
				if (pos > (double)(len - 1)) pos = (double)(len - 1);
				lo = (int)pos;
				if (lo > len - 2) lo = (len > 1) ? len - 2 : 0;
				frac = pos - (double)lo;
				a = (double)fn[base + (int64_t)lo * stride];
				b = (double)fn[base + (int64_t)((lo + 1 < len) ? lo + 1 : lo) * stride];
				nv = a * (1.0 - frac) + b * frac;
				if (fabs(nv - (double)fu[o]) > worst) worst = fabs(nv - (double)fu[o]);
				fu[o] = (float)nv;
			}
		}
		if (worst < MD_INVERT_TOL) break;
	}
}

/* ============================== output ============================== */

#define MD_PATH_MAX 2048

static int md_write(const md_ctx *c, const char *suffix, const float *vol, int nframe, gzModes gz) {
	nifti_image *n = c->tmpl;
	void *savedata = n->data;
	int saved_nt = n->nt, saved_ndim = n->ndim, saved_dt = n->datatype, saved_nbyper = n->nbyper;
	int64_t saved_nvox = n->nvox;
	float saved_slope = n->scl_slope, saved_inter = n->scl_inter;
	char *saved_fname = n->fname, *saved_iname = n->iname;
	/* n->data is NULL here by design: the template's payload is freed after repacking and only
	   its header is retained.  We swap in our own buffer for the write and restore NULL after. */
	char *fname = NULL;
	int rc;
	size_t nb = (size_t)c->n3 * (size_t)nframe;
	float *buf;
	/* nifti_save() derives the output name by stripping the extension from nim->fname and
	   appending the postfix, so point fname at "<prefix>.nii" for the duration of the write. */
	fname = (char *)malloc(strlen(c->prefix) + 8);
	if (!fname) return 1;
	snprintf(fname, strlen(c->prefix) + 8, "%s.nii", c->prefix);
	buf = (float *)nii_malloc(nb, sizeof(float));
	memcpy(buf, vol, nb * sizeof(float));
	n->fname = fname;
	n->iname = fname;
	n->data = buf;
	n->nt = nframe; n->dim[4] = nframe;
	n->ndim = (nframe > 1) ? 4 : 3; n->dim[0] = n->ndim;
	n->nvox = (int64_t)c->n3 * nframe;
	n->datatype = DT_FLOAT32; n->nbyper = 4;
	/* float32 output, unscaled: the reference stores uint16 + scl (manifest §2) and we
	   deliberately do not, because every gate is stated in physical units. */
	n->scl_slope = 1.0f; n->scl_inter = 0.0f;
	rc = nifti_save(n, suffix, gz) ? 1 : 0;
	n->fname = saved_fname; n->iname = saved_iname;
	free(fname);
	n->data = savedata;
	n->nt = saved_nt; n->dim[4] = saved_nt;
	n->ndim = saved_ndim; n->dim[0] = saved_ndim;
	n->nvox = saved_nvox;
	n->datatype = saved_dt; n->nbyper = saved_nbyper;
	n->scl_slope = saved_slope; n->scl_inter = saved_inter;
	free(buf);
	return rc;
}

static const char *const MD_EXT[3] = { ".nii.gz", ".nii.zst", ".nii" };

/* Write one output under the temporary prefix and report the path actually produced.
 *
 * The extension nifti_save() picks depends on gzMode and FSLOUTPUTTYPE, so it has to be
 * discovered by probing -- which means any STALE candidate must be removed first.  A leftover
 * temporary from a crashed run that happened to reuse this PID would otherwise be found by the
 * probe and renamed into place AS THE RESULT, orphaning the file we just wrote.  PID reuse is
 * routine in containers and HPC schedulers. */
static int md_write_temp(md_ctx *c, const char *tmppfx, const char *suffix,
	const float *buf, int T, gzModes gz, char *out, size_t outsz) {
	int i;
	for (i = 0; i < 3; i++) {
		snprintf(out, outsz, "%s%s%s", tmppfx, suffix, MD_EXT[i]);
		remove(out);
	}
	if (md_write(c, suffix, buf, T, gz)) return 1;
	for (i = 0; i < 3; i++) {
		FILE *f;
		snprintf(out, outsz, "%s%s%s", tmppfx, suffix, MD_EXT[i]);
		f = fopen(out, "rb");
		if (f) { fclose(f); return 0; }
	}
	out[0] = '\0';
	return 1;
}



/* ============================== --medic ============================== */

static void md_usage(void) {
	printf("Usage: niimath --medic --magnitude <e1> [<e2> ...] --phase <e1> [<e2> ...] \\\n");
	printf("                --te-ms <t1,t2,...> --total-readout-time <sec> \\\n");
	printf("                --phase-encoding-direction <i|j|k|i-|j-|k-> --out-prefix <path> [options]\n\n");
	printf("The phase-encoding POLARITY is significant: 'j' and 'j-' give opposite displacement maps\n");
	printf("(the native field map is the same). Take it from the BIDS PhaseEncodingDirection.\n\n");
	printf("Multi-Echo DIstortion Correction: estimates a B0 field map per frame from multi-echo\n");
	printf("phase and converts it to an EPI displacement map.\n\n");
	printf("Options:\n");
	printf("  --rank <N>              low-rank truncation of the field-map series (default %d; 0 disables)\n", MD_RANK_DEFAULT);
	printf("  --temporal-correction <0|1>  temporal 2*pi consistency correction (default 1)\n");
	printf("  --phase-offset <mcpc|none>   MCPC-3D-S phase-offset correction (default mcpc)\n");
	printf("  --noise-frames <N>, -f  drop N trailing frames from the outputs (default 0)\n");
	printf("  --n-cpus <N>, -n        OpenMP threads\n");
	printf("  --gz <0|1>              output compression (default: the FSLOUTPUTTYPE environment)\n");
	printf("  --weights <sel>         ROMEO weight preset: romeo|romeo2|romeo3|romeo4|romeo6 (default romeo4)\n");
	printf("  --mask <file>           use this mask verbatim for both unwrapping stages\n");
	printf("                          (default: ROMEO robustmask of the first echo's magnitude)\n");
	printf("  --save-intermediates    also write per-echo unwrapped phase, the masks, and (when\n");
	printf("                          MCPC-3D-S runs) the estimated phase offset\n\n");
	printf("Outputs: <prefix>_fieldmaps_native (Hz, distorted grid), <prefix>_fieldmaps (Hz,\n");
	printf("undistorted grid), <prefix>_displacementmaps (mm, pull map), all float32.\n\n");
	printf("Emulates the MEDIC workflow of Van et al., Imaging Neuroscience 4 (2026),\n");
	printf("doi:10.1162/IMAG.a.1262. Phase unwrapping is the MIT ROMEO port (Dymerska et al. 2020,\n");
	printf("doi:10.1002/mrm.28563). Clean-room: developed from the paper and black-box measurement,\n");
	printf("see test/medic_reference_manifest.md.\n");
}

/* Strict integer parse: the whole token must be consumed and fit.  atoi() silently accepts
   "8abc" as 8 and "abc" as 0, which turns a typo into a wrong run rather than an error. */
static int md_parse_int(const char *s, long *out) {
	char *end;
	long v;
	if (!s || !*s) return 1;
	errno = 0;
	v = strtol(s, &end, 10);
	if (end == s || *end != '\0' || errno == ERANGE) return 1;
	*out = v;
	return 0;
}

/* Strict double parse, same contract. */
static int md_parse_one_double(const char *s, double *out) {
	char *end;
	double v;
	if (!s || !*s) return 1;
	errno = 0;
	v = strtod(s, &end);
	if (end == s || *end != '\0') return 1;
	*out = v;
	return 0;
}

static int md_parse_doubles(const char *s, double *out, int maxn) {
	int n = 0;
	const char *p = s;
	while (*p) {
		char *end;
		double v;
		while (*p == ',' || *p == ' ' || *p == '[' || *p == ']') p++;
		if (!*p) break;
		v = strtod(p, &end);
		if (end == p) return -1;
		if (n >= maxn) return -2;   /* too many: report, never silently truncate */
		out[n++] = v;
		p = end;
	}
	return n;
}

int nii_medic(int argc, char *argv[]) {
	md_ctx c;
	const char *magf[MD_MAX_ECHO], *phaf[MD_MAX_ECHO];
	int nmag = 0, npha = 0, e, t, ac, rc = EXIT_FAILURE;
	int nTE = 0, have_trt = 0, have_pe = 0;
	nifti_image *ph[MD_MAX_ECHO], *mg[MD_MAX_ECHO];
	float *phase = NULL, *mag = NULL, *fields = NULL, *fu = NULL, *disp = NULL;
	int *frc = NULL;
	romeo_opts ro = romeo_opts_default();
	gzModes gz = GZ_ENVIRONMENT;
	/* MEASURED default (manifest section 4): the reference unwraps with romeo4 weights at BOTH
	   stages.  Against its own intermediates, romeo4 puts 99.76 % of in-mask voxels on the same
	   2*pi branch versus 93.01 % for ROMEO's own romeo3 default, 98.17 % for romeo6 and 87.91 %
	   for romeo2.  Overridable with --weights. */
	ro.weights_sel = RM_W_ROMEO4;
	int64_t n3;
	int Tin = 0, T = 0;

	memset(&c, 0, sizeof c);
	memset(ph, 0, sizeof ph);
	memset(mg, 0, sizeof mg);
	c.rank = MD_RANK_DEFAULT;
	c.temporal = 1;
	c.mcpc = 1;
	c.pe_axis = -1;
	c.pe_sign = 1;

	for (ac = 2; ac < argc; ac++) {
		const char *a = argv[ac];
		if (!strcmp(a, "-h") || !strcmp(a, "--help")) { md_usage(); return EXIT_SUCCESS; }
		else if (!strcmp(a, "--magnitude")) {
			while (ac + 1 < argc && argv[ac + 1][0] != '-') {
				if (nmag >= MD_MAX_ECHO) { MD_ERR("too many magnitude echoes (max %d)\n", MD_MAX_ECHO); goto done; }
				magf[nmag++] = argv[++ac];
			}
		} else if (!strcmp(a, "--phase")) {
			while (ac + 1 < argc && argv[ac + 1][0] != '-') {
				if (npha >= MD_MAX_ECHO) { MD_ERR("too many phase echoes (max %d)\n", MD_MAX_ECHO); goto done; }
				phaf[npha++] = argv[++ac];
			}
		} else if (!strcmp(a, "--te-ms") && ac + 1 < argc) {
			nTE = md_parse_doubles(argv[++ac], c.TEs, MD_MAX_ECHO);
			if (nTE == -2) { MD_ERR("--te-ms lists more than %d echo times\n", MD_MAX_ECHO); goto done; }
			if (nTE < 1) { MD_ERR("could not parse --te-ms '%s'\n", argv[ac]); goto done; }
		} else if (!strcmp(a, "--total-readout-time") && ac + 1 < argc) {
			if (md_parse_one_double(argv[++ac], &c.trt)) { MD_ERR("--total-readout-time '%s' is not a number\n", argv[ac]); goto done; }
			have_trt = 1;
		} else if (!strcmp(a, "--phase-encoding-direction") && ac + 1 < argc) {
			c.pe_axis = md_axis_index(argv[++ac], &c.pe_sign); have_pe = 1;
			if (c.pe_axis < 0) { MD_ERR("--phase-encoding-direction must be one of i j k x y z, optionally with a trailing '-' (the polarity is used: j and j- give opposite displacement maps)\n"); goto done; }
		} else if (!strcmp(a, "--out-prefix") && ac + 1 < argc) {
			c.prefix = argv[++ac];
		} else if (!strcmp(a, "--rank") && ac + 1 < argc) {
			long v;
			if (md_parse_int(argv[++ac], &v) || v < 0 || v > 100000) { MD_ERR("--rank '%s' must be a non-negative integer\n", argv[ac]); goto done; }
			c.rank = (int)v;
		} else if (!strcmp(a, "--temporal-correction") && ac + 1 < argc) {
			long v;
			if (md_parse_int(argv[++ac], &v) || (v != 0 && v != 1)) { MD_ERR("--temporal-correction must be 0 or 1\n"); goto done; }
			c.temporal = (int)v;
		} else if (!strcmp(a, "--phase-offset") && ac + 1 < argc) {
			const char *v = argv[++ac];
			if (!strcmp(v, "mcpc")) c.mcpc = 1;
			else if (!strcmp(v, "none")) c.mcpc = 0;
			else { MD_ERR("--phase-offset must be 'mcpc' or 'none'\n"); goto done; }
		} else if ((!strcmp(a, "--noise-frames") || !strcmp(a, "-f")) && ac + 1 < argc) {
			long v;
			if (md_parse_int(argv[++ac], &v) || v < 0 || v > INT_MAX) { MD_ERR("--noise-frames '%s' must be a non-negative integer\n", argv[ac]); goto done; }
			c.noiseframes = (int)v;
		} else if ((!strcmp(a, "--n-cpus") || !strcmp(a, "-n")) && ac + 1 < argc) {
			long lv;
			int nt;
			if (md_parse_int(argv[++ac], &lv) || lv < 1 || lv > 4096) { MD_ERR("--n-cpus '%s' must be a positive integer\n", argv[ac]); goto done; }
			nt = (int)lv;
#ifdef _OPENMP
			if (nt > 0) omp_set_num_threads(nt);
#else
			(void)nt;
#endif
		} else if (!strcmp(a, "--weights") && ac + 1 < argc) {
			const char *v = argv[++ac];
			if (!strcmp(v, "romeo")) ro.weights_sel = RM_W_ROMEO;
			else if (!strcmp(v, "romeo2")) ro.weights_sel = RM_W_ROMEO2;
			else if (!strcmp(v, "romeo3")) ro.weights_sel = RM_W_ROMEO3;
			else if (!strcmp(v, "romeo4")) ro.weights_sel = RM_W_ROMEO4;
			else if (!strcmp(v, "romeo6")) ro.weights_sel = RM_W_ROMEO6;
			else { MD_ERR("--weights must be one of romeo romeo2 romeo3 romeo4 romeo6\n"); goto done; }
		} else if (!strcmp(a, "--mask") && ac + 1 < argc) {
			c.maskfile = argv[++ac];
		} else if (!strcmp(a, "--save-intermediates")) {
			c.save_intermediates = 1;
		} else if (!strcmp(a, "--gz") && ac + 1 < argc) {
			long v;
			if (md_parse_int(argv[++ac], &v) || (v != 0 && v != 1)) { MD_ERR("--gz must be 0 or 1\n"); goto done; }
			gz = v ? GZ_TRUE : GZ_FALSE;
		} else {
			MD_ERR("unrecognized option '%s' (try --medic --help)\n", a);
			goto done;
		}
	}

	if (!nmag || !npha || !nTE || !have_trt || !have_pe || !c.prefix) {
		MD_ERR("missing a required option\n\n");
		md_usage();
		goto done;
	}
	if (nmag != npha) { MD_ERR("%d magnitude file(s) but %d phase file(s); one of each per echo is required\n", nmag, npha); goto done; }
	if (npha < 2) { MD_ERR("at least two echoes are required (got %d)\n", npha); goto done; }
	if (nTE != npha) { MD_ERR("%d echo time(s) given for %d echo(es)\n", nTE, npha); goto done; }
	for (e = 0; e < nTE; e++) {
		if (!(c.TEs[e] > 0.0) || !isfinite(c.TEs[e])) { MD_ERR("echo times must be finite and positive (--te-ms %g)\n", c.TEs[e]); goto done; }
		if (e && !(c.TEs[e] > c.TEs[e - 1])) { MD_ERR("echo times must be strictly increasing\n"); goto done; }
	}
	if (!(c.trt > 0.0) || !isfinite(c.trt)) { MD_ERR("--total-readout-time must be finite and positive\n"); goto done; }
	if (c.noiseframes < 0) { MD_ERR("--noise-frames must be >= 0\n"); goto done; }
	c.neco = npha;

	for (e = 0; e < c.neco; e++) {
		if (!strcmp(phaf[e], "-") || !strcmp(magf[e], "-")) {
			MD_ERR("stdin is not supported: --medic reads several synchronized inputs\n"); goto done;
		}
		ph[e] = md_read_f32(phaf[e], "phase");
		if (!ph[e]) goto done;
		mg[e] = md_read_f32(magf[e], "magnitude");
		if (!mg[e]) goto done;
		if (!md_same_grid(ph[0], ph[e]) || !md_same_grid(ph[0], mg[e])) {
			MD_ERR("echo %d does not share echo 1's grid (dimensions and world transform must match across all echoes and parts)\n", e + 1);
			goto done;
		}
		if (ph[e]->nvox != ph[0]->nvox || mg[e]->nvox != ph[0]->nvox) {
			MD_ERR("echo %d has a different frame count than echo 1\n", e + 1);
			goto done;
		}
	}
	c.tmpl = ph[0];
	c.nx = ph[0]->nx; c.ny = ph[0]->ny; c.nz = (ph[0]->nz < 1 ? 1 : ph[0]->nz);
	c.n3 = n3 = (int64_t)c.nx * c.ny * c.nz;
	if (n3 < 1 || (int64_t)ph[0]->nvox % n3 != 0) { MD_ERR("invalid image geometry\n"); goto done; }
	Tin = (int)((int64_t)ph[0]->nvox / n3);
	T = Tin - c.noiseframes;
	if (T < 1) { MD_ERR("--noise-frames %d leaves no frames (input has %d)\n", c.noiseframes, Tin); goto done; }
	c.nframe = T;

	{	/* Working set, all resident (plan §5.2 as scoped: in-RAM, documented budget).
		   phase (unwrapped in place) + mag + fields + fu + disp
		   = n3 * T * (2*neco + 3) * 4 bytes.  See medic_plan.md: streaming is a decided
		   non-goal because a 4D .nii.gz cannot be seeked anyway. */
		double gb = (double)n3 * T * (2.0 * c.neco + 3.0) * 4.0 / 1073741824.0;
		fprintf(stderr, "--medic: %dx%dx%d, %d echo(es), %d frame(s); working set ~%.2f GiB\n",
			c.nx, c.ny, c.nz, c.neco, T, gb);
	}

	{	/* Checked: n3 * neco * T * 4 can wrap size_t on a 32-bit / FORCE_INT32_MAX build, where
		   medic.c IS compiled.  Chain EVERY factor through the checked multiply -- an earlier
		   version pre-computed `n3 * neco` unchecked and only then handed it over, so the first
		   product could already have wrapped. */
		size_t big, small;
		if (nii_mul_size((size_t)n3, (size_t)c.neco, &big) ||
			nii_mul_size(big, (size_t)T, &big) ||
			nii_mul_size(big, sizeof(float), &big) ||
			nii_mul_size((size_t)n3, (size_t)T, &small) ||
			nii_mul_size(small, sizeof(float), &small)) {
			MD_ERR("%d frame(s) x %d echo(es) exceeds this build's address space\n", T, c.neco);
			goto done;
		}
	}
	phase = (float *)malloc((size_t)n3 * c.neco * T * sizeof(float));
	mag = (float *)malloc((size_t)n3 * c.neco * T * sizeof(float));
	fields = (float *)malloc((size_t)n3 * T * sizeof(float));
	fu = (float *)malloc((size_t)n3 * T * sizeof(float));
	disp = (float *)malloc((size_t)n3 * T * sizeof(float));
	if (!phase || !mag || !fields || !fu || !disp) { MD_ERR("out of memory allocating the working set\n"); goto done; }

	/* Repack to frame-major, echo-minor and rescale each echo's phase series as readphase does.
	 *
	 * Each echo's inputs are released AS SOON AS they have been repacked, rather than all of them
	 * afterwards.  Holding all 2*neco input payloads alive across the whole repack put the true
	 * peak at roughly (4*neco + 3) series while the banner reported (2*neco + 3); freeing per
	 * echo keeps the overshoot to one echo pair.  Echo 0's phase image is kept: it is the header
	 * template (c.tmpl) for every output. */
	for (e = 0; e < c.neco; e++) {
		float *pe = (float *)ph[e]->data;
		md_rescale_phase(pe, (int64_t)n3 * Tin);
		for (t = 0; t < T; t++) {
			memcpy(phase + ((int64_t)t * c.neco + e) * n3, pe + (int64_t)t * n3, (size_t)n3 * sizeof(float));
			memcpy(mag + ((int64_t)t * c.neco + e) * n3, ((float *)mg[e]->data) + (int64_t)t * n3, (size_t)n3 * sizeof(float));
		}
		nifti_image_free(mg[e]); mg[e] = NULL;
		if (e > 0) { nifti_image_free(ph[e]); ph[e] = NULL; }
		else { free(ph[0]->data); ph[0]->data = NULL; }   /* keep the header, drop the payload */
	}

	/* ---- per-frame: MCPC-3D-S -> ROMEO -> weighted regression ------------------------------- */
	{
		int failed = 0;
		/* Only when MCPC RUNS: with --phase-offset none nothing fills this, and writing it would
		   emit uninitialised heap as if it were an image. */
		float *offs = (c.save_intermediates && c.mcpc) ? (float *)calloc((size_t)n3 * T, sizeof(float)) : NULL;
		if (c.save_intermediates && c.mcpc && !offs) { MD_ERR("out of memory for the phase-offset intermediate\n"); goto done; }
		uint8_t *masks = (uint8_t *)malloc((size_t)n3 * T);
		if (!masks) { free(offs); MD_ERR("out of memory allocating the per-frame masks\n"); goto done; }
		/* ONE mask per frame, shared by the MCPC-3D-S phase-difference unwrap and the multi-echo
		   unwrap, as the reference does (manifest section 4).  Either the user's --mask, used
		   verbatim, or ROMEO's robustmask of that frame's first-echo magnitude. */
		frc = (int *)calloc((size_t)T, sizeof(int));
		if (!frc) { free(offs); free(masks); MD_ERR("out of memory\n"); goto done; }
		if (c.maskfile) {
			nifti_image *mk = md_read_f32(c.maskfile, "mask");
			int64_t q;
			if (!mk) { free(offs); free(masks); goto done; }
			if (!md_same_grid(ph[0], mk) || (int64_t)mk->nvox != n3) {
				MD_ERR("--mask must be a single 3D volume on the input grid\n");
				nifti_image_free(mk); free(offs); free(masks); goto done;
			}
			for (t = 0; t < T; t++)
				for (q = 0; q < n3; q++)
					masks[(int64_t)t * n3 + q] = (((const float *)mk->data)[q] != 0.0f) ? 1 : 0;
			nifti_image_free(mk);
		} else {
#ifdef _OPENMP
			#pragma omp parallel for schedule(dynamic)
#endif
			for (t = 0; t < T; t++)
				frc[t] = romeo_robustmask(mag + (int64_t)t * c.neco * n3, c.nx, c.ny, c.nz,
						masks + (int64_t)t * n3) ? 1 : 0;
			for (t = 0; t < T; t++) if (frc[t]) { failed = 1; break; }   /* leaves t = FIRST failure */
			if (failed) {
				MD_ERR("robustmask failed for frame %d\n", t);
				free(offs); free(masks); free(frc); frc = NULL; goto done;
			}
		}
#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
#endif
		for (t = 0; t < T; t++) {
			/* Unwrapped IN PLACE: the wrapped phase is dead once this frame is unwrapped, so a
			   separate `uw` series would just be a third full copy of the run. */
			float *p = phase + (int64_t)t * c.neco * n3;
			const float *m = mag + (int64_t)t * c.neco * n3;
			const uint8_t *mk = masks + (int64_t)t * n3;
			/* Each frame records its own status: a shared `failed` flag written from several
			   threads is a data race, and reading it to skip work makes the result
			   thread-count-dependent. */
			if (c.mcpc && md_mcpc3ds(&c, p, m, &ro, mk, offs ? offs + (int64_t)t * n3 : NULL)) { frc[t] = 1; continue; }
			if (romeo_unwrap_frame(p, m, c.neco, c.nx, c.ny, c.nz, c.neco, c.TEs, &ro, mk, NULL)) frc[t] = 1;
		}
		for (t = 0; t < T; t++) if (frc[t]) { failed = 1; break; }
		if (failed) {
			MD_ERR("phase unwrapping failed for frame %d\n", t);
			free(offs); free(masks); free(frc); frc = NULL; goto done;
		}
		/* Restrict the unwrapped phase to the mask BEFORE anything reads or writes it.
		 *
		 * MEASURED (manifest 3.7): the reference's per-echo unwrapped phase is nonzero exactly on
		 * mask >= 1.  ROMEO's region growing constrains which voxels it VISITS but leaves the rest
		 * holding their wrapped values, so without this the excluded background carries arbitrary
		 * phase into the regression, the temporal grouping and the SVD basis.  It also has to
		 * happen before --save-intermediates writes the unwrapped phase, or the saved diagnostic
		 * contradicts the pipeline it is supposed to document. */
		{
			int64_t q;
			for (t = 0; t < T; t++) {
				const uint8_t *mk = masks + (int64_t)t * n3;
				for (e = 0; e < c.neco; e++) {
					float *p = phase + ((int64_t)t * c.neco + e) * n3;
					for (q = 0; q < n3; q++) if (!mk[q]) p[q] = 0.0f;
				}
			}
		}
		if (c.save_intermediates) {
			float *tmp = (float *)malloc((size_t)n3 * T * sizeof(float));
			int wrc = 0;
			int64_t q;
			if (!tmp) { free(offs); free(masks); free(frc); frc = NULL; MD_ERR("out of memory writing intermediates\n"); goto done; }
			for (e = 0; e < c.neco; e++) {
				char sfx[64];
				for (t = 0; t < T; t++)
					memcpy(tmp + (int64_t)t * n3, phase + ((int64_t)t * c.neco + e) * n3, (size_t)n3 * sizeof(float));
				snprintf(sfx, sizeof sfx, "_unwrapped_echo-%d", e + 1);
				wrc |= md_write(&c, sfx, tmp, T, gz);
			}
			for (q = 0; q < (int64_t)n3 * T; q++) tmp[q] = (float)masks[q];
			wrc |= md_write(&c, "_masks", tmp, T, gz);
			if (offs) wrc |= md_write(&c, "_phase_offset", offs, T, gz);
			free(tmp);
			/* --save-intermediates is an explicit request; a failure to honour it is an error. */
			if (wrc) { free(offs); free(masks); free(frc); frc = NULL; MD_ERR("failed to write an intermediate\n"); goto done; }
		}
		free(offs); free(masks); free(frc); frc = NULL;
	}

	/* ---- temporal 2*pi correction ------------------------------------------------------------ */
	if (c.temporal) {
		float *mag1 = (float *)malloc((size_t)n3 * T * sizeof(float));
		int trc;
		if (!mag1) { MD_ERR("out of memory in the temporal correction\n"); goto done; }
		for (t = 0; t < T; t++) memcpy(mag1 + (int64_t)t * n3, mag + (int64_t)t * c.neco * n3, (size_t)n3 * sizeof(float));
		trc = md_temporal(&c, phase, mag1);
		free(mag1);
		if (trc) goto done;
	}

	/* ---- regression -> raw native field ------------------------------------------------------ */
#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for (t = 0; t < T; t++) {
		md_regress(&c, phase + (int64_t)t * c.neco * n3, mag + (int64_t)t * c.neco * n3, fields + (int64_t)t * n3);
	}

	/* ---- non-finite gate ---------------------------------------------------------------------
	   Checked HERE, on the field series itself, so it is independent of whether the low-rank
	   filter runs at all.  Gating it inside md_lowrank() made the check frame-count dependent:
	   with the default rank 10 and T <= 10 that function returns early, so a NaN sailed straight
	   through to the outputs.  A non-finite field map is never a usable result. */
	if (!md_all_finite(fields, (int64_t)n3 * T)) {
		MD_ERR("the estimated field maps contain non-finite values (check the input phase and "
			"magnitude for NaN/Inf)\n");
		goto done;
	}

	/* ---- rank-10 truncation ------------------------------------------------------------------ */
	if (c.rank > 0 && T > 1 && md_lowrank(fields, n3, T, c.rank)) goto done;

	/* ---- inversion and displacement ---------------------------------------------------------- */
	{
		double vox;
		{	/* Phase-encoding voxel size in MILLIMETRES (pixdim is in the header's xyz_units). */
			double A[3][3], s = 0.0, unit = xyz_units_to_mm(c.tmpl->xyz_units);
			int r;
			md_xform3(c.tmpl, A);
			for (r = 0; r < 3; r++) s += A[r][c.pe_axis] * A[r][c.pe_axis];
			vox = c.tmpl->pixdim[c.pe_axis + 1];
			if (!(vox > 0.0)) vox = sqrt(s);
			vox *= unit;
			if (!(vox > 0.0)) { MD_ERR("phase-encoding voxel size is zero\n"); goto done; }
		}
#ifdef _OPENMP
		#pragma omp parallel for schedule(static)
#endif
		for (t = 0; t < T; t++) {
			int64_t q;
			md_invert(&c, fields + (int64_t)t * n3, fu + (int64_t)t * n3);
			for (q = 0; q < n3; q++)
				disp[(int64_t)t * n3 + q] =
					(float)(-(double)c.pe_sign * (double)fu[(int64_t)t * n3 + q] * c.trt * vox);
		}
	}

	/* Fail-atomic via SIBLING TEMPORARIES, never by deleting final names.
	 *
	 * An earlier revision wrote straight to the final paths and, on failure, removed all three
	 * output names in every extension -- which destroyed a PREVIOUS run's results that this
	 * invocation had never touched, and could not restore a final file nifti_save had already
	 * truncated.  Now every output is written under a temporary prefix; only once all three have
	 * been written successfully are they renamed into place.  A failure leaves the previous run
	 * untouched and removes only files this invocation created. */
	{
		static const char *const outs[3] = { "_fieldmaps_native", "_fieldmaps", "_displacementmaps" };
		const float *bufs[3];
		char tmppfx[MD_PATH_MAX];
		char made[3][MD_PATH_MAX];
		const char *saved_prefix = c.prefix;
		int k, wrc = 0, nmade = 0;
		bufs[0] = fields; bufs[1] = fu; bufs[2] = disp;
		if ((int)strlen(saved_prefix) + 32 >= MD_PATH_MAX) {
			MD_ERR("--out-prefix is too long\n");
			goto done;
		}
		snprintf(tmppfx, sizeof tmppfx, "%s.medictmp%ld", saved_prefix, (long)getpid());
		c.prefix = tmppfx;
		for (k = 0; k < 3 && !wrc; k++) {
			wrc = md_write_temp(&c, tmppfx, outs[k], bufs[k], T, gz, made[nmade], MD_PATH_MAX);
			if (!wrc) nmade++;
		}
		c.prefix = saved_prefix;
		if (wrc) {
			/* Nothing has been renamed yet, so no final output exists or has been disturbed. */
			int q;
			MD_ERR("failed to write %s%s; previous outputs left untouched\n",
				saved_prefix, outs[k > 0 ? k - 1 : 0]);
			for (q = 0; q < nmade; q++) remove(made[q]);
			goto done;
		}
		/* All three temporaries exist; publish them.
		 *
		 * Per-file rename() is atomic, but three of them are not atomic AS A SET: if the second
		 * fails, the first is already published and cannot be put back.  That window is tiny (same
		 * directory, same filesystem, no I/O) but it is real, so on failure we report EXACTLY
		 * which outputs are in place rather than claiming a clean rollback we cannot perform. */
		for (k = 0; k < 3; k++) {
			char final[MD_PATH_MAX];
			/* The extension sits immediately after "<tmppfx><suffix>" -- computed by offset, not
			   by searching for ".nii", which finds the FIRST occurrence and mangles a legitimate
			   --out-prefix such as "out.nii". */
			const char *ext = made[k] + strlen(tmppfx) + strlen(outs[k]);
			snprintf(final, sizeof final, "%s%s%s", saved_prefix, outs[k], ext);
			if (rename(made[k], final) != 0) {
				int q;
				MD_ERR("failed to move %s into place\n", final);
				if (k > 0) {
					MD_ERR("  %d of 3 outputs were already published and CANNOT be rolled back:\n", k);
					for (q = 0; q < k; q++) MD_ERR("    %s%s%s\n", saved_prefix, outs[q],
						made[q] + strlen(tmppfx) + strlen(outs[q]));
					MD_ERR("  the output set is INCOMPLETE; delete it before reusing this prefix\n");
				}
				for (q = k; q < 3; q++) remove(made[q]);
				goto done;
			}
		}
	}
	rc = EXIT_SUCCESS;
done:
	free(phase); free(mag); free(fields); free(fu); free(disp); free(frc);
	for (e = 0; e < MD_MAX_ECHO; e++) { if (ph[e]) nifti_image_free(ph[e]); if (mg[e]) nifti_image_free(mg[e]); }
	return rc;
}

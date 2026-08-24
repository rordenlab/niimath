// medic.c - MEDIC (Multi-Echo DIstortion Correction) for niimath
//
// Clean-room emulation of the workflow in Van et al., Imaging Neuroscience 4 (2026),
// doi:10.1162/IMAG.a.1262.  No Warpkit implementation, test, build product or debug symbol was
// read.  Conventions the paper does not fix were measured through the public executables; every
// one of them is recorded, with its experiment, in test/medic_reference_manifest.md of the medic_bench repository.  Section
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
#include <sys/stat.h>
#ifndef _MSC_VER
	#include <unistd.h>   /* getpid() for the temporary output prefix */
	#include <sys/wait.h>
#else
	#include <process.h>
	#define getpid _getpid
	#ifndef S_ISREG
		#define S_ISREG(m) (((m) & _S_IFMT) == _S_IFREG)
	#endif
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
/* Convergence tolerance for the inversion fixed point, in Hz.
 *
 * NOT 1e-6: the iterate is float32 and field values run to ~200 Hz, where one ULP is ~1.5e-5 Hz,
 * so a 1e-6 threshold can never be met and the "did not converge" diagnostic fired on every frame
 * of every run -- a warning that is always on is worse than none.  1e-3 Hz is comfortably above
 * float32 resolution and corresponds to ~6e-8 mm of displacement, i.e. six orders of magnitude
 * below the 0.05 mm gate. */
#define MD_INVERT_TOL 1e-3f
#define MD_CORR_THRESH 0.98     /* paper §2.1.3 magnitude-correlation grouping */
/* Ceiling for the DENSE low-rank solver: it forms a T x T Gram matrix and runs cyclic Jacobi on
   it, so cost is O(T^3) and memory O(T^2).  Generous for real fMRI (the target workload is ~600
   frames); it exists so a pathological frame count fails with an explanation rather than wrapping
   an allocation on a 32-bit build. */
#define MD_MAX_FRAMES_DENSE 8192
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
	const float *disp, const double s_per_mm[3], int allow_omp) {
	const int64_t nxy = (int64_t)nx * ny;
	int z;
#ifndef _OPENMP
	(void)allow_omp;   /* the caller's level choice is meaningless without OpenMP */
#endif
#ifdef _OPENMP
	#pragma omp parallel for schedule(static) if (allow_omp)
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

/* Header-only read, for validating geometry before any payload is loaded or any work array is
   allocated.  Keeps the true allocation peak at (2*echoes + 3) work series plus ONE echo pair,
   instead of holding all 2*echoes payloads alive while the work set is allocated. */
static nifti_image *md_read_hdr(const char *fn, const char *what) {
	nifti_image *h = nifti_image_read(fn, 0);
	if (!h) { MD_ERR("failed to read the header of %s '%s'\n", what, fn); return NULL; }
	if (h->nvox < 1 || h->nx < 1 || h->ny < 1 || h->nz < 1) {
		MD_ERR("%s '%s' has invalid dimensions\n", what, fn); nifti_image_free(h); return NULL;
	}
	if (h->nu > 1 || h->nv > 1 || h->nw > 1) {
		MD_ERR("%s '%s' has more than 4 dimensions (5D input is out of scope)\n", what, fn);
		nifti_image_free(h); return NULL;
	}
	if ((int64_t)h->nvox > INT_MAX) {
		MD_ERR("%s '%s' exceeds INT_MAX voxels; --medic is not a huge-image-safe operation\n", what, fn);
		nifti_image_free(h); return NULL;
	}
	return h;
}

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
	int frame_parallel = 0;
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
#ifdef _OPENMP
	frame_parallel = (nt >= omp_get_max_threads());
	#pragma omp parallel for schedule(static) if (frame_parallel)
#endif
	for (t = 0; t < nt; t++)
		md_pull(in + (int64_t)t * n3, out + (int64_t)t * n3, nx, ny, nz,
			((const float *)map->data) + (int64_t)(mt == 1 ? 0 : t) * n3, s_per_mm,
			!frame_parallel);
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
 * EXPERIMENTAL, and deliberately labelled so: rank 10 is the paper's figure and a synthetic probe
 * confirmed the reference truncates at 10, but the reference's real 170-frame output retains a
 * broadband residual past component 10 (manifest 5.2) whose origin is unresolved and has NOT been
 * guessed at here.  This is the least reference-faithful stage; --rank 0 skips it entirely.
 *
 * Uses the T x T Gram matrix G = F^T F, whose eigenvectors are the right singular vectors of F.
 * Projecting each voxel's time course onto the leading k of them is exactly the truncated SVD:
 * F_k = F V_k V_k^T.  Memory is O(T^2), independent of the voxel count. */
static int md_lowrank(float *F, int64_t nvox, int T, int rank) {
	double *G = NULL, *V = NULL, *w = NULL, *P = NULL;
	int i, j, k, r, rc = 1;
	int64_t v;
	if (rank <= 0 || T <= 1 || rank >= T) return 0;   /* nothing to truncate */
	{	/* T*T is formed HERE, and neither the working-set check (which bounds n3*T) nor the
		   temporal check (which does not run under --temporal-correction 0) covers it.  On a
		   32-bit / FORCE_INT32_MAX build a large T wraps the element count to zero and the Gram
		   writes then run off the end. */
		size_t bytes;
		if (T > MD_MAX_FRAMES_DENSE) {
			MD_ERR("the low-rank filter needs a dense %dx%d matrix; %d frames exceeds this "
				"solver's %d-frame ceiling (use --rank 0 to skip it)\n",
				T, T, T, MD_MAX_FRAMES_DENSE);
			return 1;
		}
		if (nii_mul_size((size_t)T, (size_t)T, &bytes) || nii_mul_size(bytes, sizeof(double), &bytes)) {
			MD_ERR("the %dx%d low-rank matrices exceed this build's address space\n", T, T);
			return 1;
		}
	}
	G = (double *)calloc((size_t)T * T, sizeof(double));
	V = (double *)malloc((size_t)T * T * sizeof(double));
	w = (double *)malloc((size_t)T * sizeof(double));
	P = (double *)malloc((size_t)T * T * sizeof(double));
	if (!G || !V || !w || !P) { MD_ERR("out of memory in the low-rank filter\n"); goto done; }

	/* G = F^T F, accumulated in double. */
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


/* ================================ the tiered brain mask ======================================
 *
 * Gap 2 of medic_bench/MEDIC_GAPS.md.  The default mask is a union of two independent estimates
 * -- one from magnitude, one from phase coherence -- cleaned up morphologically and encoded in
 * three tiers: 2 = core, 1 = border ring, 0 = outside.  Everything that asks "is this voxel
 * valid" uses `> 0`; only the border filter distinguishes 1 from 2.
 *
 * The two components are combined precisely because they fail differently: an Otsu magnitude
 * threshold is occasionally too aggressive, while the quality map is permissive but noisy.
 *
 * The author's own note postdates the gap analysis and softens it -- any reasonable brain mask
 * plus a bit extra outside the brain would do, and `--mask-mode mindgrab` is the intended
 * replacement.  This is implemented anyway, and for a diagnostic reason rather than a fidelity
 * one: with matching masks, every remaining field-map difference is attributable to a stage
 * DOWNSTREAM of the mask, which is what makes the later milestones' numbers interpretable.
 *
 * Note `romeo_robustmask` remains the right default for the `-romeo` chain op; this replaces it
 * only inside --medic, and `--mask-mode robustmask` keeps the old behaviour for bisection. */

#define MD_MASK_BINS 256

/* --mask-mode.  TIERED is the default; ROBUSTMASK is the pre-Gap-2 behaviour, kept for
   bisection; MINDGRAB shells out to the external brainchop-mindgrab executable. */
#define MD_MASK_TIERED 0
#define MD_MASK_ROBUST 1
#define MD_MASK_MINDGRAB 2

/* skimage.filters.threshold_otsu at its default 256 bins: histogram the finite values over
   [min, max], maximise the inter-class variance, and return the BIN CENTRE.  Callers threshold
   with a strict `>`, as skimage's own examples do.  Returns NaN if there is nothing to
   threshold. */
static double md_otsu(const float *v, int64_t n) {
	double lo = 0.0, hi = 0.0, w = 0.0, cum = 0.0, cumw = 0.0, tot = 0.0, totw = 0.0;
	double best = -1.0, thr = (double)NAN;
	int64_t h[MD_MASK_BINS], i;
	int b, seen = 0;
	for (b = 0; b < MD_MASK_BINS; b++) h[b] = 0;
	for (i = 0; i < n; i++) {
		double x = (double)v[i];
		if (!isfinite(x)) continue;
		if (!seen) { lo = hi = x; seen = 1; }
		else { if (x < lo) lo = x; if (x > hi) hi = x; }
	}
	if (!seen || !(hi > lo)) return thr;
	w = (hi - lo) / MD_MASK_BINS;
	for (i = 0; i < n; i++) {
		double x = (double)v[i];
		int k;
		if (!isfinite(x)) continue;
		k = (int)((x - lo) / w);
		if (k < 0) k = 0;
		if (k >= MD_MASK_BINS) k = MD_MASK_BINS - 1;
		h[k]++;
	}
	for (b = 0; b < MD_MASK_BINS; b++) {
		double c = lo + (b + 0.5) * w;
		tot += (double)h[b];
		totw += (double)h[b] * c;
	}
	for (b = 0; b < MD_MASK_BINS - 1; b++) {
		double c = lo + (b + 0.5) * w, w1, w2, m1, m2, var;
		cum += (double)h[b];
		cumw += (double)h[b] * c;
		w1 = cum; w2 = tot - cum;
		if (!(w1 > 0.0) || !(w2 > 0.0)) continue;
		m1 = cumw / w1;
		m2 = (totw - cumw) / w2;
		var = w1 * w2 * (m1 - m2) * (m1 - m2);
		if (var > best) { best = var; thr = c; }
	}
	return thr;
}

/* Neighbour offsets for a 6-, 18- or 26-connected structuring element (scipy's
   generate_binary_structure(3, 1|2|3)).  Returns the count. */
static int md_nbr_offsets(int conn, int d[26][3]) {
	int rank = (conn == 6) ? 1 : (conn == 18) ? 2 : 3;
	int n = 0, dx, dy, dz;
	for (dz = -1; dz <= 1; dz++) for (dy = -1; dy <= 1; dy++) for (dx = -1; dx <= 1; dx++) {
		if (!dx && !dy && !dz) continue;
		if (abs(dx) + abs(dy) + abs(dz) > rank) continue;
		d[n][0] = dx; d[n][1] = dy; d[n][2] = dz; n++;
	}
	return n;
}

/* Erosion and dilation use an 18-CONNECTED element -- faces and edges, not corners.
 *
 * Gap 2 does not say which structuring element, and 6-connected (scipy's default) is the obvious
 * guess.  It is wrong: MEASURED against warpkit's own recovered mask on echo2 frame 0, a search
 * over {6, 18, 26} for erosion and dilation, {6, 26} for the connected components and {6, 26} for
 * the hole fill puts 18 far ahead -- valid-tier Dice 0.9997 at 18, against 0.9304 at 6 (13% too
 * small) and 0.9699 at 26 (4% too large).  The component and hole-fill connectivities barely
 * matter by comparison (0.9997 vs 0.9995).  See medic_bench's tools/m2_mask_search.py. */
#define MD_MORPH_CONN 18
#define MD_CC_CONN 26
#define MD_FILL_CONN 26

/* Per-frame scratch for the mask build.  Allocated once per worker rather than per call: the
   mask is built inside the frame loop, and a malloc/free pair per morphological pass would
   dominate it. */
typedef struct { int32_t *lab; int32_t *stk; uint8_t *a, *b; float *q; } md_maskwork;

static void md_maskwork_free(md_maskwork *w) {
	if (!w) return;
	free(w->lab); free(w->stk); free(w->a); free(w->b); free(w->q);
	memset(w, 0, sizeof *w);
}

static int md_maskwork_alloc(md_maskwork *w, int64_t n3) {
	memset(w, 0, sizeof *w);
	/* The flood-fill stack holds voxel indices as int32 -- one entry per voxel is already 4 MB
	   per worker on a typical frame, and --medic is not a huge-image op (it is deliberately
	   absent from kHugeSafeOps), so the narrower index is the right trade. */
	if (n3 > 2147483647LL) return 1;
	w->lab = (int32_t *)malloc((size_t)n3 * sizeof(int32_t));
	w->stk = (int32_t *)malloc((size_t)n3 * sizeof(int32_t));
	w->a = (uint8_t *)malloc((size_t)n3);
	w->b = (uint8_t *)malloc((size_t)n3);
	w->q = (float *)malloc((size_t)n3 * sizeof(float));
	if (!w->lab || !w->stk || !w->a || !w->b || !w->q) { md_maskwork_free(w); return 1; }
	return 0;
}

/* Flood-fill from every voxel of `seed` value `want` in `m`, labelling reachable voxels.
   Iterative with an explicit stack: n3 can be 2e5 and a recursive fill would blow the stack. */
static void md_flood(const uint8_t *m, int nx, int ny, int nz, int conn, uint8_t want,
	int32_t *lab, int32_t *stk, int64_t start, int32_t tag, int64_t *count) {
	int d[26][3], nd = md_nbr_offsets(conn, d), k;
	int64_t sp = 0, n = 0;
	lab[start] = tag; stk[sp++] = (int32_t)start; n = 1;
	while (sp > 0) {
		int64_t v = stk[--sp];
		int x = (int)(v % nx), y = (int)((v / nx) % ny), z = (int)(v / ((int64_t)nx * ny));
		for (k = 0; k < nd; k++) {
			int xx = x + d[k][0], yy = y + d[k][1], zz = z + d[k][2];
			int64_t u;
			if (xx < 0 || yy < 0 || zz < 0 || xx >= nx || yy >= ny || zz >= nz) continue;
			u = xx + (int64_t)nx * (yy + (int64_t)ny * zz);
			if (lab[u] || m[u] != want) continue;
			lab[u] = tag; stk[sp++] = (int32_t)u; n++;
		}
	}
	if (count) *count = n;
}

/* Keep only the largest connected component of the foreground.  A no-op on an empty mask. */
static void md_largest_cc(uint8_t *m, int nx, int ny, int nz, int conn, md_maskwork *w) {
	int64_t n3 = (int64_t)nx * ny * nz, i, best = 0;
	int32_t tag = 0, bestt = 0;
	memset(w->lab, 0, (size_t)n3 * sizeof(int32_t));
	for (i = 0; i < n3; i++) {
		int64_t c = 0;
		if (!m[i] || w->lab[i]) continue;
		md_flood(m, nx, ny, nz, conn, 1, w->lab, w->stk, i, ++tag, &c);
		if (c > best) { best = c; bestt = tag; }
	}
	if (!bestt) return;
	for (i = 0; i < n3; i++) if (m[i] && w->lab[i] != bestt) m[i] = 0;
}

/* Fill holes: background not reachable from the volume border becomes foreground.  `conn` is
   the connectivity of the BACKGROUND, which is what a structuring element means here. */
static void md_fill_holes(uint8_t *m, int nx, int ny, int nz, int conn, md_maskwork *w) {
	int64_t n3 = (int64_t)nx * ny * nz, i;
	int x, y, z;
	memset(w->lab, 0, (size_t)n3 * sizeof(int32_t));
	for (z = 0; z < nz; z++) for (y = 0; y < ny; y++) for (x = 0; x < nx; x++) {
		int64_t v;
		if (x && y && z && x < nx - 1 && y < ny - 1 && z < nz - 1) continue;   /* border only */
		v = x + (int64_t)nx * (y + (int64_t)ny * z);
		if (m[v] || w->lab[v]) continue;
		md_flood(m, nx, ny, nz, conn, 0, w->lab, w->stk, v, 1, NULL);
	}
	for (i = 0; i < n3; i++) if (!m[i] && !w->lab[i]) m[i] = 1;
}

/* Binary erosion or dilation, `iters` passes with a 6-connected structuring element.
   EROSION TREATS THE VOLUME BORDER AS INSIDE -- without that it eats inward from every FOV face
   and the mask loses a shell of brain wherever the head touches the edge of the field of view. */
static void md_morph(uint8_t *m, int nx, int ny, int nz, int iters, int dilate, md_maskwork *w) {
	int64_t n3 = (int64_t)nx * ny * nz, v;
	int d[26][3], nd = md_nbr_offsets(MD_MORPH_CONN, d), it, k, x, y, z;
	for (it = 0; it < iters; it++) {
		memcpy(w->b, m, (size_t)n3);
		for (z = 0; z < nz; z++) for (y = 0; y < ny; y++) for (x = 0; x < nx; x++) {
			int hit = 0;
			v = x + (int64_t)nx * (y + (int64_t)ny * z);
			if (w->b[v] == (dilate ? 1 : 0)) continue;      /* already at the extreme */
			for (k = 0; k < nd && !hit; k++) {
				int xx = x + d[k][0], yy = y + d[k][1], zz = z + d[k][2];
				int64_t u;
				/* Outside the volume: foreground for erosion (so a face of the FOV cannot
				   erode the mask) and background for dilation.  Neither sets `hit`, which is
				   why the two cases collapse to one `continue`. */
				if (xx < 0 || yy < 0 || zz < 0 || xx >= nx || yy >= ny || zz >= nz) continue;
				u = xx + (int64_t)nx * (yy + (int64_t)ny * zz);
				if (dilate ? w->b[u] : !w->b[u]) hit = 1;
			}
			m[v] = (uint8_t)(dilate ? (hit ? 1 : 0) : (hit ? 0 : 1));
		}
	}
}

/* Gap 2 step 2: the magnitude brain mask of ONE echo, WITHOUT any caller-specific extra erosion.
   The extra erosions the other three call sites need (branch scoring 2, temporal grouping 1,
   echo offset 0) are applied by the caller to a copy of this, never folded in here -- eroding
   before the required dilation is a different mask, and keeping the two apart is what stops a
   future caller getting that order wrong. */
static int md_mag_brainmask(const float *mag, int nx, int ny, int nz, md_maskwork *w, uint8_t *out) {
	int64_t n3 = (int64_t)nx * ny * nz, i, nz_count = 0;
	double thr = md_otsu(mag, n3);
	if (!isfinite(thr)) { memset(out, 1, (size_t)n3); return 0; }
	for (i = 0; i < n3; i++) out[i] = ((double)mag[i] > thr) ? 1 : 0;
	md_fill_holes(out, nx, ny, nz, MD_FILL_CONN, w);
	md_morph(out, nx, ny, nz, 2, 0, w);
	md_largest_cc(out, nx, ny, nz, MD_CC_CONN, w);
	md_morph(out, nx, ny, nz, 2, 1, w);
	for (i = 0; i < n3; i++) if (out[i]) nz_count++;
	/* An empty magnitude component is a failed threshold, not a statement that the frame is
	   empty; propagating it would silently produce an all-zero field map at exit 0. */
	if (!nz_count) memset(out, 1, (size_t)n3);
	return 0;
}

/* Gap 2 steps 1, 3 and 4: the tiered mask.  `out` gets 2 (core), 1 (border ring) or 0. */
/* Gap 2's three DERIVED masks -- branch scoring (2 extra erosions), temporal grouping (1) and the
   intra-frame echo offset (0) -- are all the magnitude brain mask with a different final erosion
   count.  Only the first is needed so far; `cons` receives it, or NULL if the caller does not want
   it.  The base result and the extra erosion stay separate on purpose: eroding before the
   dilation md_mag_brainmask ends with is a DIFFERENT mask, and keeping them apart is what stops a
   future caller getting that order wrong. */
#define MD_CONS_ERODE 2

static int md_tiered_mask(const float *phase, const float *mag, int neco, int nx, int ny, int nz,
	const double *TEs, const romeo_opts *ro, md_maskwork *w, uint8_t *out, uint8_t *cons,
	float *qout) {
	int64_t n3 = (int64_t)nx * ny * nz, i, ncore = 0;
	double thr;
	/* 1. voxel quality, from phase alone (all-ones magnitude weights). */
	if (romeo_voxelquality(phase, neco, nx, ny, nz, TEs, ro, w->q)) return 1;
	if (qout) memcpy(qout, w->q, (size_t)n3 * sizeof(float));
	thr = md_otsu(w->q, n3);
	if (isfinite(thr)) {
		for (i = 0; i < n3; i++) w->a[i] = ((double)w->q[i] > thr) ? 1 : 0;
		md_fill_holes(w->a, nx, ny, nz, MD_FILL_CONN, w);
		md_largest_cc(w->a, nx, ny, nz, MD_CC_CONN, w);
	} else memset(w->a, 0, (size_t)n3);
	/* 2. magnitude, shortest echo. */
	if (md_mag_brainmask(mag, nx, ny, nz, w, out)) return 1;
	if (cons) {
		memcpy(cons, out, (size_t)n3);
		md_morph(cons, nx, ny, nz, MD_CONS_ERODE, 0, w);
	}
	/* 3. core = dilate(largest_cc(erode(largest_cc(magnitude | quality), 2)), 2). */
	for (i = 0; i < n3; i++) out[i] = (uint8_t)((out[i] || w->a[i]) ? 1 : 0);
	md_largest_cc(out, nx, ny, nz, MD_CC_CONN, w);
	md_morph(out, nx, ny, nz, 2, 0, w);
	md_largest_cc(out, nx, ny, nz, MD_CC_CONN, w);
	md_morph(out, nx, ny, nz, 2, 1, w);
	for (i = 0; i < n3; i++) if (out[i]) ncore++;
	if (!ncore) return 1;
	/* 4. tiers: dilate the core by 3 and add.  2 = core, 1 = ring, 0 = outside. */
	memcpy(w->a, out, (size_t)n3);
	md_morph(w->a, nx, ny, nz, 3, 1, w);
	for (i = 0; i < n3; i++) out[i] = (uint8_t)(out[i] + w->a[i]);
	return 0;
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
	int mask_mode;                /* MD_MASK_* */
	int branch;                   /* --branch-correction: global 2*pi handling (Gap 3) */
	const char *maskfile;         /* --mask: use this mask verbatim, overriding --mask-mode */
	const char *prefix;
	nifti_image *tmpl;            /* header template (phase echo 1) */
} md_ctx;

/* readphase: rescale the stored phase onto [-pi,pi].
 *
 * The range comes from FRAME 0 ONLY, per echo, combined across echoes by the MODE -- not from
 * the whole series, and there is no already-in-radians short-circuit.  All three were MEASURED
 * against wk-medic --debug, which prints the range it settled on; the probe is
 * tools/m1_rescale_probe.py in medic_bench and its record is manifest section 3.1.
 *
 *   frame 0 only   : a frame-0 range of [-100,100] under later frames spanning [-200,200] is
 *                    reported as 100, and the mirror case (wide frame 0) as 200.
 *   mode, ties low : three echoes at [-100,100],[-100,100],[-50,50] report 100; a two-echo tie
 *                    at [-100,100],[-200,200] reports min -200 and max +100, i.e. each end
 *                    independently takes the SMALLEST of the tied values.
 *   non-finite     : an echo whose frame 0 contains a NaN contributes NaN, which loses every
 *                    tie -- equivalently, non-finite entries are ignored.
 *   no shortcut    : a phantom whose stored phase spans exactly 2.0 rad comes back with the
 *                    same field as its exactly-[-pi,pi] twin (ratio 1.00005; a short-circuit
 *                    would have given 0.3183), so the range is always applied.
 *
 * Why it matters: per-echo whole-series extrema diverge from this when the stored range varies
 * by frame (one outlying frame widens the range and compresses every other) or by echo (the
 * mode is robust to one odd echo; a per-echo min/max is not). */

/* Mode of `n` doubles: most frequent value, ties broken by the smallest, non-finite ignored.
 * Returns NaN when nothing is finite.  n is the echo count, so O(n^2) is free. */
static double md_mode(const double *v, int n) {
	double best = (double)NAN;
	int i, j, bestc = 0;
	for (i = 0; i < n; i++) {
		int c = 0;
		if (!isfinite(v[i])) continue;
		for (j = 0; j < n; j++) if (v[j] == v[i]) c++;
		if (c > bestc || (c == bestc && v[i] < best)) { bestc = c; best = v[i]; }
	}
	return best;
}

static void md_frame0_range(const float *p, int64_t n3, double *mn, double *mx) {
	int64_t i;
	double lo = (double)p[0], hi = lo;
	for (i = 1; i < n3; i++) {
		double v = (double)p[i];
		if (isnan(v)) { lo = hi = (double)NAN; break; }   /* NaN propagates within an echo */
		if (v < lo) lo = v;
		if (v > hi) hi = v;
	}
	*mn = lo; *mx = hi;
}

static void md_apply_phase_scale(float *p, int64_t n, double mn, double mx) {
	double slope = MD_2PI / (mx - mn), inter = -M_PI - mn * slope;
	int64_t i;
	for (i = 0; i < n; i++) p[i] = (float)((double)p[i] * slope + inter);
}

/* ===================== Gap 3(b): global 2*pi branch selection ===============================
 *
 * MCPC-3D-S recovers the coil phase offset by extrapolating the unwrapped echo-1-minus-echo-0
 * phase difference back to t = 0.  The unwrapper only restores turns BETWEEN neighbouring voxels,
 * so the absolute turn count of that difference is undetermined and every integer N gives a
 * self-consistent (offset, field) pair -- a ladder of candidates one wrap of field apart.  ROMEO's
 * global correction picks the rung whose field is smallest, which is right whenever the bulk field
 * sits well inside the half-wrap 1/(2*dTE) and is a coin toss when it does not.  Because the
 * choice is remade on every frame from that frame's own noisy median, a subject whose field sits
 * near the boundary gets ISOLATED FRAMES displaced by exactly one wrap: a time series that is
 * stable to a fraction of a hertz across most frames and a full wrap out on a scattered subset.
 *
 * No statistic computed from the phase can rank the candidates against each other -- for evenly
 * spaced echoes they all predict the same recorded phase.  What CAN be detected is a candidate
 * that is not a rung at all: the offset comes from one unwrapping and the field from a second, so
 * when the first tips, the pipeline pairs an offset from one rung with a field from another, and
 * the line fitted through the offset-removed echoes no longer passes through the origin.  Its
 * INTERCEPT is that signature, and it is a global constant, identical at every echo.
 *
 * So the procedure is a filter followed by a rule, and the order is load-bearing: the intercept
 * DISCARDS candidates the data contradicts, and only among the survivors does the smallest-field
 * prior choose.  The prior can therefore never override the evidence.
 *
 * Specification: MEDIC_GAPS.md Gap 3, and warpkit's notes/phase-offset-ambiguities.pdf section 8
 * (an approved theory reference; no warpkit source was read).  The reference logs its own
 * per-frame decision under --debug -- "branch selection: n=+0 (intercepts={...})" -- which is what
 * this was validated against.
 *
 * Two details are easy to get wrong and both are stated in the note.  The intercept reduction is a
 * PLAIN median, because it estimates one global constant; the tie-break reduction is a median
 * WEIGHTED BY THE SECOND ECHO'S SQUARED MAGNITUDE, because it summarises a wide right-skewed
 * distribution whose tails come from air-tissue interfaces, and the field is read from a
 * difference whose noise is dominated by the weaker echo.  A mean, an M-estimator, an argmax, or
 * m0^2 weighting all flip between consecutive frames. */

/* Median of `n` floats, destructive.  Quickselect, no comparator -- medic.c is in the wasm build
   and emscripten's qsort dispatches its comparator through call_indirect (AGENTS.md). */
static float md_select_kth(float *v, int64_t n, int64_t k) {
	int64_t lo = 0, hi = n - 1;
	while (lo < hi) {
		float p = v[(lo + hi) / 2], t;
		int64_t i = lo, j = hi;
		while (i <= j) {
			while (v[i] < p) i++;
			while (v[j] > p) j--;
			if (i <= j) { t = v[i]; v[i] = v[j]; v[j] = t; i++; j--; }
		}
		if (k <= j) hi = j;
		else if (k >= i) lo = i;
		else return v[k];
	}
	return v[lo];
}

static double md_median(float *v, int64_t n) {
	if (n < 1) return (double)NAN;
	if (n & 1) return (double)md_select_kth(v, n, n / 2);
	{	/* Even n: the average of the two central order statistics.  The second selection runs
		   over the already-partitioned lower half, so it is cheap. */
		double a = (double)md_select_kth(v, n, n / 2);
		double b = (double)md_select_kth(v, n / 2, n / 2 - 1);
		return 0.5 * (a + b);
	}
}

/* Weighted median: the value at which the cumulative weight, over voxels sorted by value, first
   reaches half the total.  Heapsort on (value, weight) pairs -- again comparator-free.  Only
   reached when more than one candidate survives the filter, which is the uncommon case. */
static void md_wsort(float *v, float *w, int64_t n) {
	int64_t start, end;
	for (start = n / 2 - 1; start >= 0; start--) {
		int64_t root = start;
		for (;;) {
			int64_t child = 2 * root + 1, sw;
			if (child >= n) break;
			if (child + 1 < n && v[child] < v[child + 1]) child++;
			if (!(v[root] < v[child])) break;
			{ float t = v[root]; v[root] = v[child]; v[child] = t;
			  t = w[root]; w[root] = w[child]; w[child] = t; }
			sw = child; root = sw;
		}
	}
	for (end = n - 1; end > 0; end--) {
		int64_t root = 0;
		{ float t = v[0]; v[0] = v[end]; v[end] = t; t = w[0]; w[0] = w[end]; w[end] = t; }
		for (;;) {
			int64_t child = 2 * root + 1;
			if (child >= end) break;
			if (child + 1 < end && v[child] < v[child + 1]) child++;
			if (!(v[root] < v[child])) break;
			{ float t = v[root]; v[root] = v[child]; v[child] = t;
			  t = w[root]; w[root] = w[child]; w[child] = t; }
			root = child;
		}
	}
}

static double md_wmedian(float *v, float *w, int64_t n) {
	double tot = 0.0, cum = 0.0;
	int64_t i;
	if (n < 1) return (double)NAN;
	md_wsort(v, w, n);
	for (i = 0; i < n; i++) tot += (double)w[i];
	if (!(tot > 0.0)) return (double)v[n / 2];
	for (i = 0; i < n; i++) { cum += (double)w[i]; if (cum >= 0.5 * tot) return (double)v[i]; }
	return (double)v[n - 1];
}

/* Choose the branch N in {-1, 0, +1} to add to the unwrapped phase difference `duw`.
 *
 * `omega` is the generous unwrapping mask (binary), `omegac` the conservative interior mask used
 * for both reductions.  Returns 0 on success with *nsel set; nonzero on allocation failure.
 * A degenerate or unmeasurable case leaves *nsel = 0, which changes nothing -- a failed fit is
 * not evidence for any candidate. */
static int md_branch_select(const md_ctx *c, const float *phase, const float *mag,
	const float *duw, const romeo_opts *ro, const uint8_t *omega, const uint8_t *omegac,
	int *nsel, double *scores) {
	const int64_t n3 = c->n3;
	const double dTE = c->TEs[1] - c->TEs[0];           /* ms */
	const double t0 = c->TEs[0] * 1e-3, t1 = c->TEs[1] * 1e-3, dt = t1 - t0;   /* s */
	const double k = c->TEs[0] / dTE;
	double delta, sigma, cN[3], best = 0.0;
	float *psi = NULL, *buf = NULL, *fld = NULL, *wgt = NULL;
	int64_t i, nc = 0, j;
	int n, rc = 1, cons[3], ncons = 0;
	romeo_opts o = *ro;
	double TEs2[2];

	*nsel = 0;
	if (scores) scores[0] = scores[1] = scores[2] = (double)NAN;
	if (c->neco < 2 || !(fabs(dTE) > 1e-12)) return 0;
	/* delta = |W(2*pi*k)| is the analytic separation between neighbouring candidates' offsets.
	   When TE0/dTE is an integer the offset does not move at all, the candidates are
	   indistinguishable, and no choice can change the output.  Tested against a tolerance, not
	   exact zero: W(2*pi*k) evaluates to a rounding residue near 1e-16 rather than to zero. */
	delta = fabs((double)md_wrapf(MD_2PI * k));
	if (!(delta > 1e-9)) return 0;

	for (i = 0, nc = 0; i < n3; i++) if (omegac[i]) nc++;
	if (nc < 1) return 0;                                /* no interior to reduce over */

	psi = (float *)malloc((size_t)2 * n3 * sizeof(float));
	buf = (float *)malloc((size_t)nc * sizeof(float));
	fld = (float *)malloc((size_t)3 * nc * sizeof(float));
	wgt = (float *)malloc((size_t)nc * sizeof(float));
	if (!psi || !buf || !fld || !wgt) goto done;

	TEs2[0] = c->TEs[0]; TEs2[1] = c->TEs[1];
	o.nTE = 2; o.TEs[0] = TEs2[0]; o.TEs[1] = TEs2[1];
	o.te_epi = 0; o.template_echo = 1; o.individual = 0;
	o.correctglobal = 1;   /* the second imposition of the global rule; see the note, section 5 */

	for (n = -1; n <= 1; n++) {
		double sl;
		for (i = 0; i < n3; i++) {
			float off = md_wrapf((double)phase[i] - k * ((double)duw[i] + MD_2PI * n));
			psi[i] = md_wrapf((double)phase[i] - (double)off);
			psi[n3 + i] = md_wrapf((double)phase[n3 + i] - (double)off);
		}
		if (romeo_unwrap_frame(psi, mag, 2, c->nx, c->ny, c->nz, 2, TEs2, &o, omega, NULL))
			goto done;
		for (i = 0, j = 0; i < n3; i++) {
			if (!omegac[i]) continue;
			sl = ((double)psi[n3 + i] - (double)psi[i]) / dt;      /* rad/s */
			buf[j] = (float)((double)psi[i] - sl * t0);            /* intercept at t = 0 */
			fld[(n + 1) * nc + j] = (float)(sl / MD_2PI);          /* Hz */
			if (n == -1) wgt[j] = (float)((double)mag[n3 + i] * (double)mag[n3 + i]);
			j++;
		}
		cN[n + 1] = fabs(md_median(buf, nc));
		if (scores) scores[n + 1] = cN[n + 1];
	}

	/* sigma = min(max_N c_N, delta).  Taking the SMALLER of the observed and analytic scales makes
	   the consistency test harder to pass, which biases the procedure toward leaving N = 0 -- the
	   safe direction, since a failed fit is not evidence for any candidate. */
	sigma = cN[0];
	if (cN[1] > sigma) sigma = cN[1];
	if (cN[2] > sigma) sigma = cN[2];
	if (delta < sigma) sigma = delta;
	if (!(sigma > 0.0) || !(sigma <= DBL_MAX)) { rc = 0; goto done; }

	for (n = -1; n <= 1; n++) if (cN[n + 1] < 0.5 * sigma) cons[ncons++] = n;
	if (ncons == 1) *nsel = cons[0];
	else if (ncons > 1) {
		/* A tie: the survivors are indistinguishable to the data, so the prior decides -- the
		   smallest |field|, by the m1^2-weighted median over the interior mask. */
		int b;
		for (b = 0; b < ncons; b++) {
			double f;
			for (j = 0; j < nc; j++) { buf[j] = fld[(cons[b] + 1) * nc + j]; }
			/* md_wmedian sorts BOTH arrays, so the weights must be a scratch copy: reusing the
			   shared `wgt` would leave it permuted for the next survivor. */
			f = fabs(md_wmedian(buf, wgt, nc));
			if (b == 0 || f < best) { best = f; *nsel = cons[b]; }
			for (j = 0; j < nc; j++) wgt[j] = 0.0f;   /* rebuilt below */
			for (i = 0, j = 0; i < n3; i++) if (omegac[i]) { wgt[j] = (float)((double)mag[n3 + i] * (double)mag[n3 + i]); j++; }
		}
	}
	rc = 0;
done:
	free(psi); free(buf); free(fld); free(wgt);
	return rc;
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
	const uint8_t *mask, const uint8_t *omegac, double *blog, float *offset_out) {
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
		o.individual = 0;
		/* Gap 3(a): ROMEO's own global correction, ON at BOTH unwrap calls -- this one and the
		   multi-echo unwrap.  It removes a whole-volume 2*pi offset that is otherwise left
		   wherever the region growing happens to seed, which makes the field's absolute level
		   depend on the mask.  --branch-correction 0 restores the previous behaviour. */
		o.correctglobal = c->branch ? 1 : 0;
		if (romeo_unwrap_frame(hipp, hipm, 1, c->nx, c->ny, c->nz, 1, &te1, &o, mask, NULL)) {
			MD_ERR("ROMEO failed while unwrapping the MCPC-3D-S phase difference\n");
			goto done;
		}
	}
	if (c->branch && omegac) {
		/* Gap 3(b): choose the 2*pi branch of the unwrapped difference before extrapolating it
		   back to t = 0.  Applied here, to `hipp`, so the offset below is computed from the
		   corrected difference exactly as it would have been from an uncorrected one. */
		int nsel = 0;
		if (md_branch_select(c, phase, mag, hipp, ro, mask, omegac, &nsel, blog ? blog + 1 : NULL)) {
			MD_ERR("out of memory during 2*pi branch selection\n");
			goto done;
		}
		if (blog) blog[0] = (double)nsel;
		if (nsel) for (i = 0; i < n3; i++) hipp[i] = (float)((double)hipp[i] + MD_2PI * nsel);
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
 * `uw` and `mag` are neco * nframe volumes, echo-major within frame. */
static int md_temporal(const md_ctx *c, float *uw, const float *mag, const uint8_t *masks) {
	const int64_t n3 = c->n3;
	const int T = c->nframe;
	double *mu = NULL, *sd = NULL, *corr = NULL;
	float *acc = NULL, *snap = NULL, *allacc = NULL;
	int32_t *cnt = NULL, *allcnt = NULL;   /* per-voxel count of frames valid at that voxel */
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
	cnt = (int32_t *)malloc((size_t)n3 * sizeof(int32_t));
	/* snapshot of every frame's FIRST-echo unwrapped phase, so group means are order-independent */
	snap = (float *)malloc((size_t)n3 * T * sizeof(float));
	if (!mu || !sd || !corr || !acc || !snap || !cnt) { MD_ERR("out of memory in the temporal correction\n"); goto done; }
	for (t = 0; t < T; t++)
		memcpy(snap + (int64_t)t * n3, uw + ((int64_t)t * c->neco) * n3, (size_t)n3 * sizeof(float));

	for (t = 0; t < T; t++) {
		const float *m = mag + (int64_t)t * c->neco * n3;
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
			const float *a = mag + (int64_t)t * c->neco * n3;
			const float *b = mag + (int64_t)v * c->neco * n3;
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
		/* Fast path: when EVERY frame is in this frame's group -- the common case, because a
		   quiescent run has all magnitudes correlating well above 0.98 -- the accumulation is
		   identical for every t, so compute it once and reuse it.  EXACT, not an approximation:
		   the same values summed in the same order, hoisted out of the t loop.  Takes the inner
		   work from O(T * group * n3) (~7.7e9 adds on the 170-frame demo) to O(T * n3). */
		if (ng == T) {
			if (!allacc) {
				allacc = (float *)malloc((size_t)n3 * sizeof(float));
				allcnt = (int32_t *)malloc((size_t)n3 * sizeof(int32_t));
				if (!allacc || !allcnt) { MD_ERR("out of memory in the temporal correction\n"); goto done; }
				for (i = 0; i < n3; i++) { allacc[i] = 0.0f; allcnt[i] = 0; }
				for (u = 0; u < T; u++) {
					const float *pu = snap + (int64_t)u * n3;
					const uint8_t *mu2 = masks ? masks + (int64_t)u * n3 : NULL;
					/* `!= 0`, NOT the mask value -- the mask carries TIERS (see the sibling
					   accumulation below). */
					if (mu2) for (i = 0; i < n3; i++) { allacc[i] += pu[i]; allcnt[i] += (mu2[i] != 0); }
					else     for (i = 0; i < n3; i++) { allacc[i] += pu[i]; allcnt[i]++; }
				}
			}
			memcpy(acc, allacc, (size_t)n3 * sizeof(float));
			memcpy(cnt, allcnt, (size_t)n3 * sizeof(int32_t));
		} else {
		/* Accumulate the group mean PER VOXEL over the frames that are valid AT THAT VOXEL.
		 *
		 * Masks are per frame and generally differ between frames, so a fixed group size would
		 * average structural zeros from frames where the voxel is outside the mask into the mean
		 * of frames where it is inside -- biasing the reference and, worse, letting an EXCLUDED
		 * voxel be pushed off zero by a 2*pi correction, silently undoing the mask gating. */
		for (i = 0; i < n3; i++) { acc[i] = 0.0f; cnt[i] = 0; }
		for (u = 0; u < T; u++) {
			const float *p;
			const uint8_t *mu_ = masks ? masks + (int64_t)u * n3 : NULL;
			if (corr[(size_t)t * T + u] < MD_CORR_THRESH) continue;
			p = snap + (int64_t)u * n3;   /* snapshot, not the live (partly corrected) series */
			/* Branchless: the phase is ALREADY zero outside the mask (gated before this
			   function runs), so the sum needs no test -- only the per-voxel valid count does.
			   This inner loop runs T * group_size * n3 times (~7.7e9 on the 170-frame demo), so a
			   per-voxel branch here is worth removing. */
			/* `!= 0`, NOT the mask value: the mask carries TIERS (2 core, 1 border ring), so
			   a bare `cnt[i] += mu_[i]` would count every core voxel twice and halve the group
			   mean. */
			if (mu_) for (i = 0; i < n3; i++) { acc[i] += p[i]; cnt[i] += (mu_[i] != 0); }
			else     for (i = 0; i < n3; i++) { acc[i] += p[i]; cnt[i]++; }
		}
		}
		{
			float *p1 = uw + ((int64_t)t * c->neco) * n3;
			const uint8_t *mt = masks ? masks + (int64_t)t * n3 : NULL;
			for (i = 0; i < n3; i++) {
				double ref, n;
				if (mt && !mt[i]) continue;          /* excluded here: leave the gated zero alone */
				if (cnt[i] < 2) continue;            /* no other valid frame to move toward */
				ref = (double)acc[i] / (double)cnt[i];
				n = nearbyint((ref - (double)p1[i]) / MD_2PI);
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
					if (mt && !mt[i]) continue;      /* keep excluded voxels at their gated zero */
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
	free(mu); free(sd); free(corr); free(acc); free(snap); free(cnt); free(allacc); free(allcnt);
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
static int md_invert(const md_ctx *c, const float *fn, float *fu, int64_t *nfold, int64_t *nunconv,
	int allow_omp) {
	const int nx = c->nx, ny = c->ny, nz = c->nz;
	const int m = c->pe_axis;
	const int64_t n3 = c->n3;
	const int64_t stride = (m == 0) ? 1 : ((m == 1) ? nx : (int64_t)nx * ny);
	const int len = (m == 0) ? nx : ((m == 1) ? ny : nz);
	int it, converged = 0;
	int64_t i, folds = 0, slow = 0;
#ifndef _OPENMP
	(void)allow_omp;
#endif
	for (i = 0; i < n3; i++) fu[i] = 0.0f;
	for (it = 0; it < MD_INVERT_ITERS; it++) {
		double worst = 0.0;
		int z;
		/* On the LAST iteration, count the voxels still moving: a global max is dominated by a
		   handful of oscillating folded voxels, which made every frame look unconverged. */
		const int last = (it == MD_INVERT_ITERS - 1);
		slow = 0;
#ifdef _OPENMP
		#pragma omp parallel for schedule(static) reduction(max:worst) reduction(+:slow) if (allow_omp)
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
				double frac, a, b, nv, delta;
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
				delta = fabs(nv - (double)fu[o]);
				if (delta > worst) worst = delta;
				if (last && delta > (double)MD_INVERT_TOL) slow++;
				fu[o] = (float)nv;
			}
		}
		if (worst < MD_INVERT_TOL) { converged = 1; break; }
	}
	/* Folding detector: where pe_sign*d(field*TRT)/d(PE index) <= -1 the forward map is not
	   monotone, so the inverse is genuinely multi-valued and the fixed point picks one branch.
	   That is where the reference and this implementation disagree most (manifest 3.4 records
	   p99 1.4 mm, max 9.0 mm there), so it is worth reporting rather than hiding. */
	{
		const int64_t stride2 = stride;
		int z;
		for (z = 0; z < nz; z++) {
			int x, y;
			for (y = 0; y < ny; y++) for (x = 0; x < nx; x++) {
				int64_t o = (int64_t)x + (int64_t)y * nx + (int64_t)z * nx * ny;
				int idx = (m == 0) ? x : ((m == 1) ? y : z);
				double dd;
				if (idx + 1 >= len) continue;
				dd = (double)c->pe_sign *
					((double)fu[o + stride2] - (double)fu[o]) * c->trt;
				if (dd <= -1.0) folds++;
			}
		}
	}
	if (nfold) *nfold = folds;
	if (nunconv) *nunconv = converged ? 0 : slow;
	return converged ? 0 : 1;
}

/* ============================== output ============================== */

#define MD_PATH_MAX 2048

static void md_restore_backups(char bak[][MD_PATH_MAX], char final[][MD_PATH_MAX], int n) {
	int i;
	for (i = 0; i < n; i++) {
		if (!bak[i][0]) continue;
		if (rename(bak[i], final[i]) != 0)
			MD_ERR("rollback could not restore %s; its backup remains at %s\n", final[i], bak[i]);
	}
}

static int md_write(const md_ctx *c, const char *suffix, const float *vol, int nframe, gzModes gz) {
	nifti_image *n = c->tmpl;
	void *savedata = n->data;
	int saved_nt = n->nt, saved_ndim = n->ndim, saved_dt = n->datatype, saved_nbyper = n->nbyper;
	int64_t saved_nvox = n->nvox;
	float saved_slope = n->scl_slope, saved_inter = n->scl_inter;
	char *saved_fname = n->fname, *saved_iname = n->iname;
	/* n->data is NULL here by design: the template's payload is freed after repacking and only
	   its header is retained.  We lend the output buffer for the write and restore NULL after. */
	char *fname = NULL;
	int rc;
	/* nifti_save() derives the output name by stripping the extension from nim->fname and
	   appending the postfix, so point fname at "<prefix>.nii" for the duration of the write. */
	fname = (char *)malloc(strlen(c->prefix) + 8);
	if (!fname) return 1;
	snprintf(fname, strlen(c->prefix) + 8, "%s.nii", c->prefix);
	n->fname = fname;
	n->iname = fname;
	/* The writer consumes data synchronously and does not take ownership.  Borrow the resident
	   series directly: copying each 4D output added pure memory bandwidth and an allocation. */
	n->data = (void *)vol;
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
	return rc;
}

/* Run an external executable with an argv vector.  execvp/_spawnvp rather than system(), so the
   output prefix -- which is user text -- never reaches a shell.  Returns 0 only on exit 0. */
static int md_spawn(const char *const argv[]) {
#ifdef _MSC_VER
	intptr_t rc = _spawnvp(_P_WAIT, argv[0], (const char *const *)argv);
	return (rc == 0) ? 0 : 1;
#else
	pid_t pid = fork();
	int st = 0;
	if (pid < 0) return 1;
	if (pid == 0) {
		execvp(argv[0], (char *const *)argv);
		_exit(127);
	}
	if (waitpid(pid, &st, 0) < 0) return 1;
	return (WIFEXITED(st) && WEXITSTATUS(st) == 0) ? 0 : 1;
#endif
}

/* --mask-mode mindgrab: build the tiers from the external brainchop-mindgrab executable.
 *
 * `brainchop-mindgrab in.nii --mask m.nii --border MM` already exposes exactly what this stage
 * needs -- a brain mask, and the same mask grown by a stated distance -- so two invocations
 * differenced give the core and ring tiers directly, with no Otsu, no morphology chain, no
 * connected-component pass and no voxel-quality union.  `--border MM` IS the author's "a bit
 * extra outside the brain".
 *
 * ONE mask for the whole run, from the temporal mean of the shortest echo's magnitude, not one
 * per frame: mindgrab is a neural network and 2*T invocations would dominate the run by two
 * orders of magnitude.  That matches --mask, which is also a single 3D volume.
 *
 * This is a NEW EXTERNAL DEPENDENCY and --medic has none, so it is never a silent fallback: if
 * the executable is missing the run fails with an explanation, and it is not the default. */
static int md_mindgrab_mask(md_ctx *c, const float *mag, int64_t n3, int T, uint8_t *masks) {
	char pre[MD_PATH_MAX], fmag[MD_PATH_MAX], fcore[MD_PATH_MAX], fring[MD_PATH_MAX], border[32];
	float *mean = NULL;
	nifti_image *mc = NULL, *mr = NULL;
	md_ctx tc = *c;
	double mm;
	int64_t i;
	int t, rc = 1;

	/* Ring width: Gap 2 grows the core by 3 voxels with a 6-connected element, so the closest
	   millimetre equivalent is three of the SMALLEST voxel dimension. */
	mm = c->tmpl->dx;
	if (c->tmpl->dy > 0.0 && c->tmpl->dy < mm) mm = c->tmpl->dy;
	if (c->tmpl->dz > 0.0 && c->tmpl->dz < mm) mm = c->tmpl->dz;
	if (!(mm > 0.0)) mm = 1.0;
	snprintf(border, sizeof border, "%g", 3.0 * mm);

	snprintf(pre, sizeof pre, "%s_mindgrab%d", c->prefix, (int)getpid());
	snprintf(fmag, sizeof fmag, "%s_mag.nii", pre);
	snprintf(fcore, sizeof fcore, "%s_core.nii", pre);
	snprintf(fring, sizeof fring, "%s_ring.nii", pre);

	mean = (float *)calloc((size_t)n3, sizeof(float));
	if (!mean) { MD_ERR("out of memory for the mindgrab input\n"); return 1; }
	for (t = 0; t < T; t++)
		for (i = 0; i < n3; i++) mean[i] += mag[(int64_t)t * c->neco * n3 + i];
	for (i = 0; i < n3; i++) mean[i] /= (float)T;

	tc.prefix = pre;
	if (md_write(&tc, "_mag", mean, 1, GZ_FALSE)) {
		MD_ERR("could not write the temporary magnitude '%s' for mindgrab\n", fmag);
		goto done;
	}
	{	/* --border 0 gives the core tier */
		const char *av[8] = { "brainchop-mindgrab", fmag, "--mask", fcore, "--border", "0", NULL, NULL };
		if (md_spawn(av)) {
			MD_ERR("brainchop-mindgrab failed or is not on PATH.  --mask-mode mindgrab needs it "
				"(see https://github.com/neuroneural/brainchop-mindgrab); use --mask-mode tiered, "
				"or run mindgrab yourself and pass the result with --mask.\n");
			goto done;
		}
	}
	{
		const char *av[8] = { "brainchop-mindgrab", fmag, "--mask", fring, "--border", border, NULL, NULL };
		if (md_spawn(av)) { MD_ERR("brainchop-mindgrab failed while growing the border ring\n"); goto done; }
	}
	mc = md_read_f32(fcore, "mindgrab core mask");
	if (!mc) goto done;
	mr = md_read_f32(fring, "mindgrab border mask");
	if (!mr) goto done;
	if (!md_same_grid(c->tmpl, mc) || !md_same_grid(c->tmpl, mr) ||
		(int64_t)mc->nvox != n3 || (int64_t)mr->nvox != n3) {
		MD_ERR("mindgrab returned a mask off the input grid\n");
		goto done;
	}
	for (i = 0; i < n3; i++) {
		const float *vc = (const float *)mc->data, *vr = (const float *)mr->data;
		masks[i] = (uint8_t)(vc[i] >= 0.5f ? 2 : (vr[i] >= 0.5f ? 1 : 0));
	}
	{
		int64_t ncore = 0;
		for (i = 0; i < n3; i++) if (masks[i] == 2) ncore++;
		if (!ncore) { MD_ERR("mindgrab returned an empty brain mask\n"); goto done; }
	}
	for (t = 1; t < T; t++) memcpy(masks + (int64_t)t * n3, masks, (size_t)n3);
	rc = 0;
done:
	nifti_image_free(mc); nifti_image_free(mr);
	free(mean);
	remove(fmag); remove(fcore); remove(fring);
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
	if (md_write(c, suffix, buf, T, gz)) goto fail;
	for (i = 0; i < 3; i++) {
		FILE *f;
		snprintf(out, outsz, "%s%s%s", tmppfx, suffix, MD_EXT[i]);
		f = fopen(out, "rb");
		if (f) { fclose(f); return 0; }
	}
	/* Wrote something we cannot name -- e.g. a .hdr/.img pair.  Fall through and clean up. */
fail:
	/* A short write, a compressor-close error or a disk-full partial leaves a file the caller
	   cannot know about, because it only tracks paths we successfully returned.  Remove ours. */
	for (i = 0; i < 3; i++) {
		snprintf(out, outsz, "%s%s%s", tmppfx, suffix, MD_EXT[i]);
		remove(out);
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
	printf("                          EXPERIMENTAL: rank 10 is what the paper specifies, but the\n");
	printf("                          reference's own output retains a broadband residual past\n");
	printf("                          component 10 whose origin is unresolved, so this stage is the\n");
	printf("                          least reference-faithful part of the pipeline (--rank 0 skips it)\n");
	printf("  --temporal-correction <0|1>  temporal 2*pi consistency correction (default 1)\n");
	printf("  --phase-offset <mcpc|none>   MCPC-3D-S phase-offset correction (default mcpc)\n");
	printf("  --noise-frames <N>, -f  drop N trailing frames from the outputs (default 0)\n");
	printf("  --n-cpus <N>, -n        OpenMP threads\n");
	printf("  --gz <0|1>              output compression (default: the FSLOUTPUTTYPE environment)\n");
	printf("  --weights <sel>         ROMEO weight preset: romeo|romeo2|romeo3|romeo4|romeo6 (default romeo4)\n");
	printf("  --mask-mode <sel>       brain mask: tiered|robustmask|mindgrab (default tiered)\n");
	printf("                          tiered     Otsu magnitude mask unioned with ROMEO's phase\n");
	printf("                                     voxel quality, encoded as 2 core / 1 border ring\n");
	printf("                                     / 0 outside; every stage tests > 0\n");
	printf("                          robustmask ROMEO's robustmask of the first echo's magnitude,\n");
	printf("                                     the pre-tier behaviour, kept for bisection\n");
	printf("                          mindgrab   two calls to the external brainchop-mindgrab, one\n");
	printf("                                     per tier; needs it on PATH and fails if absent\n");
	printf("  --branch-correction <0|1>  global 2*pi handling (default 1): ROMEO's own global\n");
	printf("                          correction at both unwrap stages\n");
	printf("  --mask <file>           use this mask verbatim for both unwrapping stages,\n");
	printf("                          overriding --mask-mode; a supplied mask is the core tier,\n");
	printf("                          so it has no border ring\n");
	printf("  --save-intermediates    also write per-echo unwrapped phase, the masks, and (when\n");
	printf("                          MCPC-3D-S runs) the estimated phase offset\n\n");
	printf("Outputs: <prefix>_fieldmaps_native (Hz, distorted grid), <prefix>_fieldmaps (Hz,\n");
	printf("undistorted grid), <prefix>_displacementmaps (mm, pull map), all float32.\n\n");
	printf("Emulates the MEDIC workflow of Van et al., Imaging Neuroscience 4 (2026),\n");
	printf("doi:10.1162/IMAG.a.1262. Phase unwrapping is the MIT ROMEO port (Dymerska et al. 2020,\n");
	printf("doi:10.1002/mrm.28563). Clean-room: developed from the paper and black-box measurement,\n");
	printf("see test/medic_reference_manifest.md in the medic_bench repository.\n");
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
	double phmin[MD_MAX_ECHO], phmax[MD_MAX_ECHO];   /* frame-0 extrema, per echo */
	int *frc = NULL;
	uint8_t *maskbuf = NULL;   /* per-frame masks, retained through the temporal correction */
	romeo_opts ro = romeo_opts_default();
	gzModes gz = GZ_ENVIRONMENT;
	/* MEASURED default (manifest section 4): the reference unwraps with romeo4 weights at BOTH
	   stages.  Against its own intermediates, romeo4 puts 99.76 % of in-mask voxels on the same
	   2*pi branch versus 93.01 % for ROMEO's own romeo3 default, 98.17 % for romeo6 and 87.91 %
	   for romeo2.  Overridable with --weights. */
	ro.weights_sel = RM_W_ROMEO4;
	int64_t n3;
	int Tin = 0, T = 0;
#ifdef _OPENMP
	int frame_parallel = 0;
#endif

#ifdef _OPENMP
	/* Belt and braces against an OMP_NESTED=true / OMP_MAX_ACTIVE_LEVELS>1 environment: every
	   parallel region below is designed as the single active level. */
	omp_set_max_active_levels(1);
#endif
	memset(&c, 0, sizeof c);
	memset(ph, 0, sizeof ph);
	memset(mg, 0, sizeof mg);
	c.rank = MD_RANK_DEFAULT;
	c.temporal = 1;
	c.branch = 1;
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
		} else if (!strcmp(a, "--branch-correction") && ac + 1 < argc) {
			long v;
			if (md_parse_int(argv[++ac], &v) || (v != 0 && v != 1)) { MD_ERR("--branch-correction must be 0 or 1\n"); goto done; }
			c.branch = (int)v;
		} else if (!strcmp(a, "--mask-mode") && ac + 1 < argc) {
			const char *v = argv[++ac];
			if (!strcmp(v, "tiered")) c.mask_mode = MD_MASK_TIERED;
			else if (!strcmp(v, "robustmask")) c.mask_mode = MD_MASK_ROBUST;
			else if (!strcmp(v, "mindgrab")) c.mask_mode = MD_MASK_MINDGRAB;
			else { MD_ERR("--mask-mode must be tiered, robustmask or mindgrab\n"); goto done; }
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
	{	/* One-file NIfTI is the stated scope.  A PAIR setting makes nifti_save write .hdr/.img,
		   which the output transaction cannot discover or publish -- it used to leave orphaned
		   temporaries and exit 1 AFTER doing all the work.  Refuse before computing anything. */
		const char *fot = getenv("FSLOUTPUTTYPE");
		if (fot && strstr(fot, "PAIR")) {
			MD_ERR("FSLOUTPUTTYPE=%s selects a .hdr/.img pair; --medic writes one-file NIfTI only "
				"(use NIFTI, NIFTI_GZ or NIFTI_ZST)\n", fot);
			goto done;
		}
	}
	if (c.noiseframes < 0) { MD_ERR("--noise-frames must be >= 0\n"); goto done; }
	c.neco = npha;

	for (e = 0; e < c.neco; e++) {
		if (!strcmp(phaf[e], "-") || !strcmp(magf[e], "-")) {
			MD_ERR("stdin is not supported: --medic reads several synchronized inputs\n"); goto done;
		}
		ph[e] = md_read_hdr(phaf[e], "phase");
		if (!ph[e]) goto done;
		mg[e] = md_read_hdr(magf[e], "magnitude");
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
	/* Gap 3(a): ROMEO's global correction is ON for the multi-echo unwrap too, not just for the
	   MCPC-3D-S phase difference.  romeo_opts_default() leaves it off, which is right for the
	   `-romeo` chain op and wrong here. */
	ro.correctglobal = c.branch ? 1 : 0;
#ifdef _OPENMP
	/* Parallelise the per-frame loops whenever there is more than one frame.
	 *
	 * NOT `T >= omp_get_max_threads()`: the inner regions these loops would be yielding to are
	 * too small to compensate.  romeo.c has exactly two OpenMP regions -- rm_calculateweights and
	 * rm_compute_b0 (unused here) -- so romeo_robustmask has NO inner parallelism at all, and the
	 * unwrap keeps only the weight kernel while MCPC and the region growing stay serial.  Gating
	 * on the thread count therefore made every run with fewer frames than cores fall off a cliff:
	 * measured on a 7-frame 64x64x40 case, 0.05 s at 7 threads versus 0.15 s at 8 -- three times
	 * slower for asking for one more thread, and worse on a many-core node, where most task-fMRI
	 * runs have fewer frames than cores.  With `T > 1` the outer region is active whenever it can
	 * do anything, and omp_set_max_active_levels(1) keeps the inner regions from nesting under
	 * it; at T == 1 the outer region is inactive, so the inner ones still get the team. */
	frame_parallel = (T > 1);
#endif

	{	/* Working set, all resident (plan §5.2 as scoped: in-RAM, documented budget).
		   phase (unwrapped in place) + mag + fields + fu + disp
		   = n3 * T * (2*neco + 3) * 4 bytes.  Streaming is a decided
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

	/* Load, rescale and repack ONE ECHO PAIR AT A TIME.
	 *
	 * Geometry was validated from headers alone above, so no payload is resident when the work
	 * arrays are allocated.  An earlier version loaded all 2*neco payloads during validation and
	 * only freed them during the repack -- which is AFTER the work allocation, so the true peak
	 * was ~(4*neco + 3) series while the banner reported (2*neco + 3).  Now the overshoot is one
	 * echo pair regardless of echo count.  Echo 0's phase image is retained, payload dropped, as
	 * the header template (c.tmpl) for every output. */
	for (e = 0; e < c.neco; e++) {
		nifti_image *pi, *mi;
		float *pe;
		if (e > 0) { nifti_image_free(ph[e]); ph[e] = NULL; }   /* header stub no longer needed */
		nifti_image_free(mg[e]); mg[e] = NULL;
		pi = md_read_f32(phaf[e], "phase");
		if (!pi) goto done;
		mi = md_read_f32(magf[e], "magnitude");
		if (!mi) { nifti_image_free(pi); goto done; }
		pe = (float *)pi->data;
		md_frame0_range(pe, n3, &phmin[e], &phmax[e]);   /* rescale happens after the loop */
		for (t = 0; t < T; t++) {
			memcpy(phase + ((int64_t)t * c.neco + e) * n3, pe + (int64_t)t * n3, (size_t)n3 * sizeof(float));
			memcpy(mag + ((int64_t)t * c.neco + e) * n3, ((float *)mi->data) + (int64_t)t * n3, (size_t)n3 * sizeof(float));
		}
		nifti_image_free(mi);
		if (e == 0) { free(pi->data); pi->data = NULL; nifti_image_free(ph[0]); ph[0] = pi; c.tmpl = ph[0]; }
		else nifti_image_free(pi);
	}

	/* One range for the whole run: the mode across echoes of frame 0's extrema.  Applied here
	   rather than per echo above, because the mode cannot be taken until every echo has been
	   seen -- see md_mode / md_frame0_range. */
	{
		double mn = md_mode(phmin, c.neco), mx = md_mode(phmax, c.neco);
		if (!(mx > mn)) {
			MD_ERR("frame 0 of the phase data has no usable range (mode across echoes gave "
				"[%g, %g]); a constant or all-NaN first frame cannot be rescaled\n", mn, mx);
			goto done;
		}
		md_apply_phase_scale(phase, (int64_t)n3 * c.neco * T, mn, mx);
		/* Provenance, and the only direct read-out of the rule: the reference prints the same
		   thing.  Stderr, never stdout -- the output image may be going to stdout. */
		fprintf(stderr, "--medic: phase rescale range [%g, %g] (frame 0, mode across echoes)\n", mn, mx);
	}

	/* ---- per-frame: MCPC-3D-S -> ROMEO -> weighted regression ------------------------------- */
	{
		/* EVERY pointer this block owns is declared and NULLed here, above the first `goto`.
		   Declaring one below a goto that targets a label which frees it is legal C but leaves
		   it INDETERMINATE, and clang then deletes the guards -- see the nifti_bptf incident in
		   AGENTS.md.  Nothing below may add a declaration lower down. */
		int failed = 0;
		float *offs = NULL;      /* phase offset, only when MCPC actually runs */
		float *qual = NULL;      /* ROMEO voxel quality, one of the two mask ingredients */
		uint8_t *masks = NULL;   /* per-frame tiers: 2 core, 1 border ring, 0 outside */
		uint8_t *cons = NULL;    /* per-frame conservative interior mask, for branch scoring */
		double *blog = NULL;     /* per-frame branch decision and its three intercepts */
		int want_qual = c.save_intermediates && !c.maskfile && c.mask_mode == MD_MASK_TIERED;
		if (c.save_intermediates && c.mcpc) {
			/* Only when MCPC RUNS: with --phase-offset none nothing fills this, and writing it
			   would emit uninitialised heap as if it were an image. */
			offs = (float *)calloc((size_t)n3 * T, sizeof(float));
			if (!offs) { MD_ERR("out of memory for the phase-offset intermediate\n"); goto done; }
		}
		if (want_qual) {
			/* The voxel-quality map is one of the two mask ingredients, so it belongs with the
			   mask in the intermediates: it is the only way to see WHICH ingredient moved a
			   boundary. */
			qual = (float *)calloc((size_t)n3 * T, sizeof(float));
			if (!qual) { free(offs); MD_ERR("out of memory for the voxel-quality intermediate\n"); goto done; }
		}
		masks = (uint8_t *)malloc((size_t)n3 * T);
		if (!masks) { free(offs); free(qual); MD_ERR("out of memory allocating the per-frame masks\n"); goto done; }
		if (c.branch && c.mcpc && c.neco > 1) {
			cons = (uint8_t *)malloc((size_t)n3 * T);
			if (!cons) { free(offs); free(qual); free(masks); MD_ERR("out of memory allocating the branch-scoring masks\n"); goto done; }
			if (c.save_intermediates) {
				/* Diagnostic only, and COLLECTED rather than printed from inside the frame loop:
				   printing there would interleave across threads and reorder run to run.  It is
				   read nowhere -- a flag must never gate work that reaches a voxel. */
				blog = (double *)calloc((size_t)4 * T, sizeof(double));
				if (!blog) { free(offs); free(qual); free(masks); free(cons); MD_ERR("out of memory\n"); goto done; }
			}
		}
		/* ONE mask per frame, shared by the MCPC-3D-S phase-difference unwrap and the multi-echo
		   unwrap, as the reference does (manifest section 4).
		 *
		 * The values are TIERS -- 2 core, 1 border ring, 0 outside (Gap 2).  Everything that asks
		 * "is this voxel valid" tests `> 0`; only the border filter distinguishes 1 from 2.  A
		 * user's --mask is binary and maps to the CORE tier, since a supplied mask has no border
		 * ring, so the border filter has nothing to do on that path.
		 *
		 * ROMEO must be handed a BINARY mask: rm_build_ctx forms `magnitude * (float)mask[i]`, so
		 * a tier of 2 would double the magnitude weights over the core. */
		frc = (int *)calloc((size_t)T, sizeof(int));
		if (!frc) { free(offs); free(qual); free(masks); free(cons); free(blog); MD_ERR("out of memory\n"); goto done; }
		if (c.maskfile) {
			nifti_image *mk = md_read_f32(c.maskfile, "mask");
			int64_t q;
			if (!mk) { free(offs); free(qual); free(masks); free(cons); free(blog); goto done; }
			if (!md_same_grid(ph[0], mk) || (int64_t)mk->nvox != n3) {
				MD_ERR("--mask must be a single 3D volume on the input grid\n");
				nifti_image_free(mk); free(offs); free(qual); free(masks); free(cons); free(blog); goto done;
			}
			/* MEASURED contract (manifest 3.7): in-mask is `>= 1`, not merely nonzero.  A
			   fractional probability map is therefore NOT a mask -- threshold it first
			   (`niimath p.nii -thr 0.5 -bin m.nii`).  NaN fails the comparison and is excluded,
			   which is the safe direction. */
			for (t = 0; t < T; t++)
				for (q = 0; q < n3; q++)
					masks[(int64_t)t * n3 + q] = (((const float *)mk->data)[q] >= 1.0f) ? 2 : 0;
			nifti_image_free(mk);
			{	/* An all-fractional or empty mask would silently zero every voxel downstream. */
				int64_t nz = 0;
				for (q = 0; q < n3; q++) if (masks[q]) nz++;
				if (nz == 0) {
					MD_ERR("--mask '%s' has no voxel >= 1 (MEDIC's in-mask test is `>= 1`; "
						"threshold a probability map first, e.g. niimath m.nii -thr 0.5 -bin m_bin.nii)\n",
						c.maskfile);
					free(offs); free(qual); free(masks); free(cons); free(blog); free(frc); frc = NULL; goto done;
				}
			}
		} else if (c.mask_mode == MD_MASK_ROBUST) {
#ifdef _OPENMP
			#pragma omp parallel for schedule(dynamic) if (frame_parallel)
#endif
			for (t = 0; t < T; t++)
				frc[t] = romeo_robustmask(mag + (int64_t)t * c.neco * n3, c.nx, c.ny, c.nz,
						masks + (int64_t)t * n3) ? 1 : 0;
			for (t = 0; t < T; t++) if (frc[t]) { failed = 1; break; }   /* leaves t = FIRST failure */
			if (failed) {
				MD_ERR("robustmask failed for frame %d\n", t);
				free(offs); free(qual); free(masks); free(cons); free(blog); free(frc); frc = NULL; goto done;
			}
			/* Binary robustmask, promoted to the core tier: --mask-mode robustmask exists to
			   bisect against the pre-Gap-2 behaviour, and it has no border ring either. */
			{
				int64_t q;
				for (q = 0; q < (int64_t)n3 * T; q++) masks[q] = (uint8_t)(masks[q] ? 2 : 0);
			}
		} else if (c.mask_mode == MD_MASK_MINDGRAB) {
			if (md_mindgrab_mask(&c, mag, n3, T, masks)) {
				free(offs); free(qual); free(masks); free(cons); free(blog); free(frc); frc = NULL; goto done;
			}
		} else {
#ifdef _OPENMP
			#pragma omp parallel
#endif
			{
				md_maskwork w;
				int tt;
				int oom = md_maskwork_alloc(&w, n3);
#ifdef _OPENMP
				#pragma omp for schedule(dynamic)
#endif
				for (tt = 0; tt < T; tt++) {
					if (oom) { frc[tt] = 1; continue; }
					frc[tt] = md_tiered_mask(phase + (int64_t)tt * c.neco * n3,
							mag + (int64_t)tt * c.neco * n3, c.neco,
							c.nx, c.ny, c.nz, c.TEs, &ro, &w,
							masks + (int64_t)tt * n3,
							cons ? cons + (int64_t)tt * n3 : NULL,
							qual ? qual + (int64_t)tt * n3 : NULL) ? 1 : 0;
				}
				md_maskwork_free(&w);
			}
			for (t = 0; t < T; t++) if (frc[t]) { failed = 1; break; }
			if (failed) {
				MD_ERR("the tiered brain mask failed for frame %d (out of memory, or the frame "
					"has no usable brain component)\n", t);
				free(offs); free(qual); free(masks); free(cons); free(blog); free(frc); frc = NULL; goto done;
			}
		}
		if (cons && (c.maskfile || c.mask_mode != MD_MASK_TIERED)) {
			/* The other mask modes build no magnitude brain mask of their own, so the conservative
			   scoring mask is the supplied one eroded instead.  The specification asks only that
			   it sit well inside the unwrapping mask and exclude edge and low-signal voxels --
			   nothing depends on how either is built. */
#ifdef _OPENMP
			#pragma omp parallel
#endif
			{
				md_maskwork w;
				int tt;
				int oom = md_maskwork_alloc(&w, n3);
#ifdef _OPENMP
				#pragma omp for schedule(dynamic)
#endif
				for (tt = 0; tt < T; tt++) {
					int64_t q;
					if (oom) { frc[tt] = 1; continue; }
					for (q = 0; q < n3; q++)
						cons[(int64_t)tt * n3 + q] = (uint8_t)(masks[(int64_t)tt * n3 + q] ? 1 : 0);
					md_morph(cons + (int64_t)tt * n3, c.nx, c.ny, c.nz, MD_CONS_ERODE, 0, &w);
				}
				md_maskwork_free(&w);
			}
			for (t = 0; t < T; t++) if (frc[t]) { failed = 1; break; }
			if (failed) {
				MD_ERR("out of memory building the branch-scoring mask for frame %d\n", t);
				free(offs); free(qual); free(masks); free(cons); free(blog); free(frc); frc = NULL; goto done;
			}
		}
#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic) if (frame_parallel)
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
			/* ROMEO gets a BINARY mask: rm_build_ctx forms `magnitude * (float)mask[i]`, so a
			   core tier of 2 would silently double the magnitude weights inside the brain. */
			uint8_t *mkb = (uint8_t *)malloc((size_t)n3);
			int64_t q;
			if (!mkb) { frc[t] = 1; continue; }
			for (q = 0; q < n3; q++) mkb[q] = (uint8_t)(mk[q] ? 1 : 0);
			if (c.mcpc && md_mcpc3ds(&c, p, m, &ro, mkb, cons ? cons + (int64_t)t * n3 : NULL,
					blog ? blog + (size_t)4 * t : NULL,
					offs ? offs + (int64_t)t * n3 : NULL)) { free(mkb); frc[t] = 1; continue; }
			if (romeo_unwrap_frame(p, m, c.neco, c.nx, c.ny, c.nz, c.neco, c.TEs, &ro, mkb, NULL)) frc[t] = 1;   /* ro.correctglobal set from --branch-correction above */
			free(mkb);
		}
		if (blog) {
			for (t = 0; t < T; t++)
				fprintf(stderr, "--medic: frame %d branch selection: n=%+d (intercepts "
					"{-1: %.4e, 0: %.4e, +1: %.4e})\n", t, (int)blog[(size_t)4 * t],
					blog[(size_t)4 * t + 1], blog[(size_t)4 * t + 2], blog[(size_t)4 * t + 3]);
		}
		for (t = 0; t < T; t++) if (frc[t]) { failed = 1; break; }
		if (failed) {
			MD_ERR("phase unwrapping failed for frame %d\n", t);
			free(offs); free(qual); free(masks); free(cons); free(blog); free(frc); frc = NULL; goto done;
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
			if (!tmp) { free(offs); free(qual); free(masks); free(cons); free(blog); free(frc); frc = NULL; MD_ERR("out of memory writing intermediates\n"); goto done; }
			for (e = 0; e < c.neco; e++) {
				char sfx[64];
				for (t = 0; t < T; t++)
					memcpy(tmp + (int64_t)t * n3, phase + ((int64_t)t * c.neco + e) * n3, (size_t)n3 * sizeof(float));
				snprintf(sfx, sizeof sfx, "_unwrapped_echo-%d", e + 1);
				wrc |= md_write(&c, sfx, tmp, T, gz);
			}
			for (q = 0; q < (int64_t)n3 * T; q++) tmp[q] = (float)masks[q];
			wrc |= md_write(&c, "_masks", tmp, T, gz);
			if (qual) wrc |= md_write(&c, "_voxelquality", qual, T, gz);
			if (offs) wrc |= md_write(&c, "_phase_offset", offs, T, gz);
			free(tmp);
			/* --save-intermediates is an explicit request; a failure to honour it is an error. */
			if (wrc) { free(offs); free(qual); free(masks); free(cons); free(blog); free(frc); frc = NULL; MD_ERR("failed to write an intermediate\n"); goto done; }
		}
		free(offs); free(qual); free(cons); free(blog); free(frc); frc = NULL;
		maskbuf = masks;   /* retained: md_temporal must know which samples are valid */
	}

	/* ---- temporal 2*pi correction ------------------------------------------------------------ */
	if (c.temporal) {
		if (md_temporal(&c, phase, mag, maskbuf)) goto done;
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
		int64_t inv_folds = 0, inv_unconv = 0;
#ifdef _OPENMP
		const int outer_par = frame_parallel;
#else
		const int outer_par = 0;
#endif
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
		/* ONE active level of parallelism, chosen by frame count.
		 *
		 * md_invert() is itself parallel over slices, so a parallel frame loop around it nested
		 * two regions: with nesting disabled (the default) the inner team collapsed to one thread;
		 * with OMP_NESTED=true it oversubscribed.  For T > 1, parallelise frames and run
		 * md_invert serially: each frame also has useful work without an inner OpenMP region.
		 * For T = 1, leave the outer loop serial and parallelise md_invert over slices.  The
		 * accumulators use reductions either way, so diagnostics are exact and thread-count
		 * independent. */
#ifdef _OPENMP
		#pragma omp parallel for schedule(static) reduction(+:inv_folds) reduction(+:inv_unconv) if (outer_par)
#endif
		for (t = 0; t < T; t++) {
			int64_t q, nf = 0, nu = 0;
			md_invert(&c, fields + (int64_t)t * n3, fu + (int64_t)t * n3, &nf, &nu, !outer_par);
			inv_folds += nf; inv_unconv += nu;
			for (q = 0; q < n3; q++)
				disp[(int64_t)t * n3 + q] =
					(float)(-(double)c.pe_sign * (double)fu[(int64_t)t * n3 + q] * c.trt * vox);
		}
		/* Report, do not hide: an unconverged frame or a folded region is exactly where this
		   implementation and the reference disagree most. */
		if (inv_unconv || inv_folds) {
			double tot = (double)n3 * T;
			fprintf(stderr, "--medic: displacement inversion: %lld voxel(s) (%.3f%%) still moving "
				"by >%g Hz after %d iterations; %lld folded adjacent pair(s) (%.3f%% of voxels) "
				"mark columns where the forward map is not monotone, so the inverse is "
				"multi-valued and the branch chosen is arbitrary\n",
				(long long)inv_unconv, 100.0 * (double)inv_unconv / tot, (double)MD_INVERT_TOL,
				MD_INVERT_ITERS, (long long)inv_folds, 100.0 * (double)inv_folds / tot);
		}
	}

	/* Publish the three outputs as a RECOVERABLE transaction.
	 *
	 * History, because two previous attempts were wrong: writing straight to the final paths and
	 * deleting them on failure destroyed a PREVIOUS run's results; sibling temporaries plus a
	 * plain rename loop fixed that but still published a MIXED set if the second rename failed
	 * (new file 1, old files 2 and 3), which is worse than either a clean old set or a clean new
	 * one because it looks usable.
	 *
	 * So: write all three temporaries; move any existing finals aside to `.medicbak`; rename the
	 * temporaries in; on ANY failure put the backups back and remove our temporaries.  Individual
	 * renames are atomic and same-directory, so the restore path is itself renames.  This is
	 * recoverable rather than atomic -- a crash between two renames still leaves a mixed set --
	 * and it is deliberately NOT labelled atomic. */
	{
		static const char *const outs[3] = { "_fieldmaps_native", "_fieldmaps", "_displacementmaps" };
		const float *bufs[3];
		char tmppfx[MD_PATH_MAX], made[3][MD_PATH_MAX];
		char final[3][MD_PATH_MAX], bak[3][MD_PATH_MAX];
		int had_final[3] = { 0, 0, 0 };
		const char *saved_prefix = c.prefix;
		int k, q, wrc = 0, nmade = 0, npub = 0;
		bufs[0] = fields; bufs[1] = fu; bufs[2] = disp;
		/* Reserve room for the longest suffix we ever append: "<out>.medictmp<pid>.nii.gz". */
		if ((int)strlen(saved_prefix) + 64 >= MD_PATH_MAX) { MD_ERR("--out-prefix is too long\n"); goto done; }
		snprintf(tmppfx, sizeof tmppfx, "%s.medictmp%ld", saved_prefix, (long)getpid());
		c.prefix = tmppfx;
		for (k = 0; k < 3 && !wrc; k++) {
			wrc = md_write_temp(&c, tmppfx, outs[k], bufs[k], T, gz, made[nmade], MD_PATH_MAX);
			if (!wrc) nmade++;
		}
		c.prefix = saved_prefix;
		if (wrc) {
			MD_ERR("failed to write %s%s; existing outputs left untouched\n",
				saved_prefix, outs[k > 0 ? k - 1 : 0]);
			for (q = 0; q < nmade; q++) remove(made[q]);
			goto done;
		}
		for (k = 0; k < 3; k++) {
			/* The extension sits immediately after "<tmppfx><suffix>" -- by offset, never by
			   searching for ".nii", which finds the FIRST occurrence and mangles a legitimate
			   --out-prefix such as "out.nii". */
			const char *ext = made[k] + strlen(tmppfx) + strlen(outs[k]);
			snprintf(final[k], MD_PATH_MAX, "%s%s%s", saved_prefix, outs[k], ext);
			snprintf(bak[k], MD_PATH_MAX, "%s.medicbak%ld", final[k], (long)getpid());
		}
		/* PREFLIGHT the destinations before touching anything.  A path occupied by a directory
		   (or anything that is not a regular file) must be refused, not renamed out of the way --
		   moving a user's directory aside would be a surprising side effect, and a non-empty one
		   cannot then be cleaned up. */
		for (k = 0; k < 3; k++) {
			struct stat st;
			if (stat(final[k], &st) == 0) {
				had_final[k] = 1;
				if (S_ISREG(st.st_mode)) continue;
				MD_ERR("%s exists and is not a regular file; refusing to replace it\n", final[k]);
				for (q = 0; q < 3; q++) remove(made[q]);
				goto done;
			}
			if (errno != ENOENT) {
				MD_ERR("cannot inspect existing output %s; refusing to replace it\n", final[k]);
				for (q = 0; q < 3; q++) remove(made[q]);
				goto done;
			}
		}
		for (k = 0; k < 3; k++) {   /* move existing finals aside (absent is fine) */
			/* Clear any stale backup FIRST, even when there is no final to move aside: a
			   leftover <final>.medicbak<pid> from a crashed run that reused this PID would
			   otherwise sit on disk forever, since nothing later renames or removes it. */
			remove(bak[k]);
			if (!had_final[k]) { bak[k][0] = '\0'; continue; }
			if (rename(final[k], bak[k]) != 0) {
				MD_ERR("failed to preserve existing output %s; aborting publication\n", final[k]);
				md_restore_backups(bak, final, k);
				for (q = 0; q < 3; q++) remove(made[q]);
				goto done;
			}
		}
		for (k = 0; k < 3; k++) {
			remove(final[k]);   /* MSVC rename() will not replace an existing destination */
			if (rename(made[k], final[k]) != 0) {
				MD_ERR("failed to publish %s; rolling back\n", final[k]);
				for (q = 0; q < npub; q++) remove(final[q]);              /* undo our publishes */
				md_restore_backups(bak, final, 3);
				for (q = k; q < 3; q++) remove(made[q]);
				goto done;
			}
			npub++;
		}
		for (k = 0; k < 3; k++) if (bak[k][0]) remove(bak[k]);
	}
	rc = EXIT_SUCCESS;
done:
	free(phase); free(mag); free(fields); free(fu); free(disp); free(frc); free(maskbuf);
	for (e = 0; e < MD_MAX_ECHO; e++) { if (ph[e]) nifti_image_free(ph[e]); if (mg[e]) nifti_image_free(mg[e]); }
	return rc;
}

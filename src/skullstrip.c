// AFNI-style surface skull stripping (-skullstrip). See skullstrip.h for the licence
// boundary, which is the single authoritative statement: normalisation, intensity prep,
// deformation and touchup are ALL adapted from public-domain AFNI (thd_brainormalize.c,
// thd_automask.c, SUMA_BrainWrap.c -- the last carries no copyright notice and is therefore
// a non-copyrightable US Government work). The ONE carve-out is SUMA_3dedge3, which wraps
// Malandain's GPL-3.0 code: do not read or adapt it, and treat anything needing the edge
// volume (-use_edge) as out of scope.

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327   /* MSVC: not in math.h without _USE_MATH_DEFINES */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>
#include "print.h"   // printfx -> stderr; diagnostics must NEVER touch stdout,
                     // which niimath uses for the output image ("-" / -gz 0 -).
#include "skullstrip.h"

#define SS_IJK(i, j, k) ((i) + (j) * nx + (k) * nxy)

// Developer hooks are COMPILE-GATED. In a normal build SS_DEV_HOOKS is undefined and these
// are not read at all, so no environment variable can change a shipped result -- validating
// the syntax was not enough, because a hidden input that alters the segmentation is the
// problem, not a malformed one. Build the benchmark harness with -DSS_DEV_HOOKS=1 to restore
// SS_VAR, SS_NODE_DBG and SS_STAGE_MAX; the manifest's experiments need them and nothing else
// does. (SKULLSTRIP_VERBOSE stays live in every build: it only prints.)
#ifdef SS_DEV_HOOKS
// Never bare atoi(): `atoi("abc")` is 0, so a typo in SS_STAGE_MAX silently truncated the
// stage machine to one pass and degraded the mask at exit 0 with nothing in the output to
// show for it. Garbage or out-of-range keeps the default and says so on stderr.
static int ss_env_int(const char *name, int dflt, int lo, int hi) {
	const char *e = getenv(name);
	char *end = NULL;
	long v;
	if (!e || !*e)
		return dflt;
	errno = 0;
	v = strtol(e, &end, 10);
	if (errno || !end || *end || v < lo || v > hi) {
		fprintf(stderr, "skullstrip: ignoring %s=\"%s\" (want an integer in [%d,%d])\n",
				name, e, lo, hi);
		return dflt;
	}
	return (int)v;
}
#else
#define ss_env_int(name, dflt, lo, hi) (dflt)
#endif

static int ss_verbose(void) {
	static int v = -1;
	if (v < 0) {
		const char *e = getenv("SKULLSTRIP_VERBOSE");
		v = (e && *e && *e != '0') ? 1 : 0;
	}
	return v;
}

// Diagnostics. NOTE the trap this macro sets: its whole body is guarded on ss_verbose(),
// so ANY function call written as an argument is not evaluated in an ordinary run. Only pass
// pure value-computing helpers (ss_mask_count, ss_mean_seg_len); never a call that does work.
// Two pipeline stages were once written as SSV arguments and silently did not run unless
// SKULLSTRIP_VERBOSE was set -- a diagnostic that changed the segmentation.
#define SSV(...)                       \
	do {                               \
		if (ss_verbose())              \
			fprintf(stderr, __VA_ARGS__); \
	} while (0)

static short ss_shortize(float x) {
	// AFNI SHORTIZE: round-to-nearest with clamping into short range.
	int v = (int)((x < 0.0f) ? (x - 0.5f) : (x + 0.5f));
	if (v > 32767)
		v = 32767;
	if (v < -32768)
		v = -32768;
	return (short)v;
}

// ---------------------------------------------------------------------------
// Cliplevels. Adapted from partial_cliplevel() / get_octant_clips() /
// pointclip() in thd_brainormalize.c.
// ---------------------------------------------------------------------------

typedef struct {
	float c000, c100, c010, c110, c001, c101, c011, c111;
	float x0, x1, dxi, y0, y1, dyi, z0, z1, dzi;
	float cmin, cmax;
} ss_clipvec;

#define SS_NHIST 32767

// oom is set to 1 if the histogram allocation fails. It CANNOT be signalled through the return
// value: 1.0f is also AFNI's legitimate result for a region with <=999 positive voxels, so an
// OOM was indistinguishable from a valid answer and ss_octant_clips' c000 < 0 error flag never
// fired. Measured by fault injection: 18 of the faults inside this function exited 0 with a
// mask up to Dice 0.88 from baseline.
static float ss_partial_cliplevel(const short *sar, int nx, int ny, int nz, float mfrac, int *oom,
		int ibot, int itop, int jbot, int jtop, int kbot, int ktop) {
	int nxy = nx * ny;
	int ii, jj, kk, qq, ncut, nold, npos, nhalf, ib, val;
	double dsum;
	int *hist;

	if (mfrac <= 0.0f || mfrac >= 0.99f)
		mfrac = 0.50f;
	if (ibot < 0)
		ibot = 0;
	if (jbot < 0)
		jbot = 0;
	if (kbot < 0)
		kbot = 0;
	if (itop >= nx)
		itop = nx - 1;
	if (jtop >= ny)
		jtop = ny - 1;
	if (ktop >= nz)
		ktop = nz - 1;
	if (itop < ibot || jtop < jbot || ktop < kbot)
		return 1.0f;

	hist = (int *)calloc(SS_NHIST + 1, sizeof(int));
	if (!hist) {
		if (oom) *oom = 1;
		return 1.0f;
	}

	dsum = 0.0;
	npos = 0;
	for (kk = kbot; kk <= ktop; kk++)
		for (jj = jbot; jj <= jtop; jj++)
			for (ii = ibot; ii <= itop; ii++) {
				val = sar[SS_IJK(ii, jj, kk)];
				if (val > 0 && val <= SS_NHIST) {
					hist[val]++;
					dsum += (double)val * (double)val;
					npos++;
				}
			}

	if (npos <= 999) {
		free(hist);
		return 1.0f;
	}

	// Start the cut so it includes the upper 65% of positive voxels, but never
	// below half the RMS.
	qq = (int)(0.65 * npos);
	ib = (int)rint(0.5 * sqrt(dsum / npos));
	for (kk = 0, ii = SS_NHIST - 1; ii >= ib && kk < qq; ii--)
		kk += hist[ii];

	ncut = ii;
	qq = 0;
	do {
		for (npos = 0, ii = ncut; ii < SS_NHIST; ii++)
			npos += hist[ii];
		nhalf = npos / 2;
		for (kk = 0, ii = ncut; ii < SS_NHIST && kk < nhalf; ii++)
			kk += hist[ii];
		nold = ncut;
		ncut = (int)(mfrac * ii);
		qq++;
	} while (qq < 20 && ncut != nold);

	free(hist);
	return (float)ncut;
}

static ss_clipvec ss_octant_clips(const short *sar, int nx, int ny, int nz, float mfrac) {
	ss_clipvec cv;
	int it = nx - 1, jt = ny - 1, kt = nz - 1;
	int ii, jj, kk, ic, jc, kc;
	long long ijk;
	double xcm = 0, ycm = 0, zcm = 0, sum = 0;
	float val;
	int oom = 0;

	memset(&cv, 0, sizeof(cv));
	cv.c000 = -1.0f; // error flag

	for (ijk = 0, kk = 0; kk < nz; kk++)
		for (jj = 0; jj < ny; jj++)
			for (ii = 0; ii < nx; ii++, ijk++) {
				val = (float)sar[ijk];
				if (val <= 0.0f)
					continue;
				sum += val;
				xcm += val * ii;
				ycm += val * jj;
				zcm += val * kk;
			}
	if (sum == 0.0)
		return cv;
	ic = (int)rint(xcm / sum);
	jc = (int)rint(ycm / sum);
	kc = (int)rint(zcm / sum);

	val = 0.5f * ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, 0, it, 0, jt, 0, kt);

	cv.c000 = ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, 0, ic + 2, 0, jc + 2, 0, kc + 2);
	cv.c100 = ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, ic - 2, it, 0, jc + 2, 0, kc + 2);
	cv.c010 = ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, 0, ic + 2, jc - 2, jt, 0, kc + 2);
	cv.c110 = ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, ic - 2, it, jc - 2, jt, 0, kc + 2);
	cv.c001 = ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, 0, ic + 2, 0, jc + 2, kc - 2, kt);
	cv.c101 = ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, ic - 2, it, 0, jc + 2, kc - 2, kt);
	cv.c011 = ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, 0, ic + 2, jc - 2, jt, kc - 2, kt);
	cv.c111 = ss_partial_cliplevel(sar, nx, ny, nz, mfrac, &oom, ic - 2, it, jc - 2, jt, kc - 2, kt);

	if (oom) {
		cv.c000 = -1.0f;   // the caller's existing error flag
		return cv;
	}

	// Don't let any octant's clip level get too small.
	if (cv.c000 < val) cv.c000 = val;
	if (cv.c100 < val) cv.c100 = val;
	if (cv.c010 < val) cv.c010 = val;
	if (cv.c110 < val) cv.c110 = val;
	if (cv.c001 < val) cv.c001 = val;
	if (cv.c101 < val) cv.c101 = val;
	if (cv.c011 < val) cv.c011 = val;
	if (cv.c111 < val) cv.c111 = val;

	{
		float mn = cv.c000, mx = cv.c000;
		float a[8] = {cv.c000, cv.c100, cv.c010, cv.c110, cv.c001, cv.c101, cv.c011, cv.c111};
		for (ii = 1; ii < 8; ii++) {
			if (a[ii] < mn) mn = a[ii];
			if (a[ii] > mx) mx = a[ii];
		}
		cv.cmin = mn;
		cv.cmax = mx;
	}

	cv.x0 = 0.5f * ic; cv.x1 = 0.5f * (ic + it);
	cv.y0 = 0.5f * jc; cv.y1 = 0.5f * (jc + jt);
	cv.z0 = 0.5f * kc; cv.z1 = 0.5f * (kc + kt);
	cv.dxi = (cv.x1 > cv.x0) ? 1.0f / (cv.x1 - cv.x0) : 0.0f;
	cv.dyi = (cv.y1 > cv.y0) ? 1.0f / (cv.y1 - cv.y0) : 0.0f;
	cv.dzi = (cv.z1 > cv.z0) ? 1.0f / (cv.z1 - cv.z0) : 0.0f;
	return cv;
}

// Cliplevel at a point, trilinearly interpolated between octant centres.
static float ss_pointclip(int ii, int jj, int kk, const ss_clipvec *cv) {
	float x1, y1, z1, x0, y0, z0;
	x1 = (ii - cv->x0) * cv->dxi; if (x1 < 0.0f) x1 = 0.0f; else if (x1 > 1.0f) x1 = 1.0f;
	y1 = (jj - cv->y0) * cv->dyi; if (y1 < 0.0f) y1 = 0.0f; else if (y1 > 1.0f) y1 = 1.0f;
	z1 = (kk - cv->z0) * cv->dzi; if (z1 < 0.0f) z1 = 0.0f; else if (z1 > 1.0f) z1 = 1.0f;
	x0 = 1.0f - x1; y0 = 1.0f - y1; z0 = 1.0f - z1;
	return cv->c000 * x0 * y0 * z0 + cv->c100 * x1 * y0 * z0 +
	       cv->c010 * x0 * y1 * z0 + cv->c110 * x1 * y1 * z0 +
	       cv->c001 * x0 * y0 * z1 + cv->c101 * x1 * y0 * z1 +
	       cv->c011 * x0 * y1 * z1 + cv->c111 * x1 * y1 * z1;
}

// ---------------------------------------------------------------------------
// Mask morphology. Adapted from thd_automask.c (public domain).
// ---------------------------------------------------------------------------

static long long ss_mask_count(long long n, const unsigned char *m) {
	long long i, c = 0;
	for (i = 0; i < n; i++)
		if (m[i])
			c++;
	return c;
}

// Keep only clusters of at least csize voxels (6-connected). AFNI clustedit3D.
// Returns 0 on success, -1 on allocation failure. It MUST report: it zeroes clusters as it
// walks, so an OOM midway leaves a partially destroyed mask that still looks plausible.
static int ss_clustedit3D(int nx, int ny, int nz, unsigned char *mmm, long long csize) {
	long long nxy = (long long)nx * ny, nxyz = nxy * nz, ijk, ijk_last = 0;
	int *stack = NULL;
	unsigned char *keep;
	long long nstack_cap = 0;

	if (nx < 1 || ny < 1 || nz < 1 || !mmm || csize < 2)
		return 0;
	keep = (unsigned char *)calloc((size_t)nxyz, 1);
	if (!keep)
		return -1;
	nstack_cap = 4096;
	stack = (int *)malloc(sizeof(int) * (size_t)nstack_cap);
	if (!stack) {
		free(keep);
		return -1;
	}

	while (1) {
		long long nnow;
		for (ijk = ijk_last; ijk < nxyz; ijk++)
			if (mmm[ijk])
				break;
		if (ijk == nxyz)
			break;
		ijk_last = ijk + 1;

		mmm[ijk] = 0;
		nnow = 1;
		stack[0] = (int)ijk;
		for (long long icl = 0; icl < nnow; icl++) {
			long long p = stack[icl];
			int ii = (int)(p % nx), jj = (int)((p % nxy) / nx), kk = (int)(p / nxy);
			int d;
			int nb[6], nnb = 0;
			if (ii > 0) nb[nnb++] = (int)(p - 1);
			if (ii < nx - 1) nb[nnb++] = (int)(p + 1);
			if (jj > 0) nb[nnb++] = (int)(p - nx);
			if (jj < ny - 1) nb[nnb++] = (int)(p + nx);
			if (kk > 0) nb[nnb++] = (int)(p - nxy);
			if (kk < nz - 1) nb[nnb++] = (int)(p + nxy);
			for (d = 0; d < nnb; d++) {
				if (!mmm[nb[d]])
					continue;
				if (nnow == nstack_cap) {
					long long ncap = nstack_cap + 4096 + nstack_cap / 8;
					int *t = (int *)realloc(stack, sizeof(int) * (size_t)ncap);
					if (!t) {
						free(stack);
						free(keep);
						return -1;
					}
					stack = t;
					nstack_cap = ncap;
				}
				mmm[nb[d]] = 0;
				stack[nnow++] = nb[d];
			}
		}
		if (nnow >= csize)
			for (long long q = 0; q < nnow; q++)
				keep[stack[q]] = 1;
	}
	memcpy(mmm, keep, (size_t)nxyz);
	free(stack);
	free(keep);
	return 0;
}

// Keep only the single largest 6-connected cluster. AFNI THD_mask_clust.
// Returns 0 on success, -1 on allocation failure. Same reason as ss_clustedit3D: measured,
// failing its 16 KB stack alloc turned "keep the largest cluster" into a no-op and published
// a mask 46% too large (Dice 0.813) at exit 0.
static int ss_mask_clust(int nx, int ny, int nz, unsigned char *mmm) {
	long long nxy = (long long)nx * ny, nxyz = nxy * nz, ijk, ijk_last = 0;
	long long nbest = 0, nstack_cap = 4096;
	int *stack = (int *)malloc(sizeof(int) * 4096);
	int *best = NULL;

	if (!mmm || !stack) {
		free(stack);
		return -1;
	}
	while (1) {
		long long nnow;
		for (ijk = ijk_last; ijk < nxyz; ijk++)
			if (mmm[ijk])
				break;
		if (ijk == nxyz)
			break;
		ijk_last = ijk + 1;
		mmm[ijk] = 0;
		nnow = 1;
		stack[0] = (int)ijk;
		for (long long icl = 0; icl < nnow; icl++) {
			long long p = stack[icl];
			int ii = (int)(p % nx), jj = (int)((p % nxy) / nx), kk = (int)(p / nxy);
			int nb[6], nnb = 0, d;
			if (ii > 0) nb[nnb++] = (int)(p - 1);
			if (ii < nx - 1) nb[nnb++] = (int)(p + 1);
			if (jj > 0) nb[nnb++] = (int)(p - nx);
			if (jj < ny - 1) nb[nnb++] = (int)(p + nx);
			if (kk > 0) nb[nnb++] = (int)(p - nxy);
			if (kk < nz - 1) nb[nnb++] = (int)(p + nxy);
			for (d = 0; d < nnb; d++) {
				if (!mmm[nb[d]])
					continue;
				if (nnow == nstack_cap) {
					long long ncap = nstack_cap + 4096 + nstack_cap / 8;
					int *t = (int *)realloc(stack, sizeof(int) * (size_t)ncap);
					if (!t) {
						free(stack);
						free(best);
						return -1;
					}
					stack = t;
					nstack_cap = ncap;
				}
				mmm[nb[d]] = 0;
				stack[nnow++] = nb[d];
			}
		}
		if (nnow > nbest) {
			int *t = (int *)realloc(best, sizeof(int) * (size_t)nnow);
			if (!t) {
				free(stack);
				free(best);
				return -1;
			}
			best = t;
			memcpy(best, stack, sizeof(int) * (size_t)nnow);
			nbest = nnow;
		}
	}
	for (long long q = 0; q < nbest; q++)
		mmm[best[q]] = 1;
	free(stack);
	free(best);
	return 0;
}

// Fill a zero voxel that is bracketed by mask within nside along ANY axis.
// AFNI THD_mask_fillin_once, taking nside as a parameter (2 for the normalisation mask, 10 for -fill_hole).
int ss_mask_fillin_once(int nx, int ny, int nz, unsigned char *mmm, int nside) {
	long long nxy = (long long)nx * ny, nxyz = nxy * nz;
	int nsx = (nx - 1) / 2, nsy = (ny - 1) / 2, nsz = (nz - 1) / 2;
	unsigned char *nnn;
	int ii, jj, kk, s;

	if (!mmm || nside <= 0)
		return 0;
	if (nsx > nside) nsx = nside;
	if (nsy > nside) nsy = nside;
	if (nsz > nside) nsz = nside;
	if (nsx == 0 && nsy == 0 && nsz == 0)
		return 0;
	nnn = (unsigned char *)calloc((size_t)nxyz, 1);
	if (!nnn)
		return -1;

	for (kk = nsz; kk < nz - nsz; kk++)
		for (jj = nsy; jj < ny - nsy; jj++)
			for (ii = nsx; ii < nx - nsx; ii++) {
				long long iv = ii + (long long)jj * nx + (long long)kk * nxy;
				int plus, minus;
				if (mmm[iv])
					continue;
				plus = minus = 0;
				for (s = 1; s <= nsx; s++) {
					if (mmm[iv + s]) plus = 1;
					if (mmm[iv - s]) minus = 1;
				}
				if (plus && minus) { nnn[iv] = 1; continue; }
				plus = minus = 0;
				for (s = 1; s <= nsy; s++) {
					if (mmm[iv + (long long)s * nx]) plus = 1;
					if (mmm[iv - (long long)s * nx]) minus = 1;
				}
				if (plus && minus) { nnn[iv] = 1; continue; }
				plus = minus = 0;
				for (s = 1; s <= nsz; s++) {
					if (mmm[iv + (long long)s * nxy]) plus = 1;
					if (mmm[iv - (long long)s * nxy]) minus = 1;
				}
				if (plus && minus) { nnn[iv] = 1; continue; }
			}
	for (long long q = 0; q < nxyz; q++)
		if (nnn[q])
			mmm[q] = 1;
	free(nnn);
	return 0;
}

// Neighbour offsets for NN1 (6), NN2 adds 12, NN3 adds 8. Edge voxels clamp their
// out-of-range neighbour to themselves, matching AFNI.
static int ss_count_nbrs(const unsigned char *m, int nx, int ny, int nz, int ii, int jj, int kk, int NN) {
	long long nxy = (long long)nx * ny;
	long long kz = (long long)kk * nxy, km = (kk == 0) ? kz : kz - nxy, kp = (kk == nz - 1) ? kz : kz + nxy;
	long long jy = (long long)jj * nx, jm = (jj == 0) ? jy : jy - nx, jp = (jj == ny - 1) ? jy : jy + nx;
	int im = (ii == 0) ? 0 : ii - 1, ip = (ii == nx - 1) ? ii : ii + 1;
	int num;
	num = m[ii + jy + km] + m[im + jy + kz] + m[ii + jm + kz] + m[ii + jp + kz] +
	      m[ip + jy + kz] + m[ii + jy + kp];
	if (NN >= 2)
		num += m[im + jy + km] + m[ii + jm + km] + m[ii + jp + km] + m[ip + jy + km] +
		       m[im + jm + kz] + m[im + jp + kz] + m[ip + jm + kz] + m[ip + jp + kz] +
		       m[im + jy + kp] + m[ii + jm + kp] + m[ii + jp + kp] + m[ip + jy + kp];
	if (NN == 3)
		num += m[im + jm + km] + m[im + jp + km] + m[ip + jm + km] + m[ip + jp + km] +
		       m[im + jm + kp] + m[im + jp + kp] + m[ip + jm + kp] + m[ip + jp + kp];
	return num;
}

// AFNI THD_mask_dilate: a zero voxel with >= ndil nonzero neighbours joins the mask.
static int ss_mask_dilate(int nx, int ny, int nz, unsigned char *mmm, int ndil, int NN) {
	long long nxy = (long long)nx * ny, nxyz = nxy * nz;
	unsigned char *nnn = (unsigned char *)calloc((size_t)nxyz, 1);
	int ii, jj, kk;
	if (!nnn)
		return -1;
	if (ndil < 1)
		ndil = 1;
	for (kk = 0; kk < nz; kk++)
		for (jj = 0; jj < ny; jj++)
			for (ii = 0; ii < nx; ii++) {
				long long p = ii + (long long)jj * nx + (long long)kk * nxy;
				if (mmm[p])
					continue;
				if (ss_count_nbrs(mmm, nx, ny, nz, ii, jj, kk, NN) >= ndil)
					nnn[p] = 1;
			}
	for (long long q = 0; q < nxyz; q++)
		if (nnn[q])
			mmm[q] = 1;
	free(nnn);
	return 0;
}

// AFNI THD_mask_erode: a mask voxel lacking a FULL neighbourhood erodes; volume
// voxels on the array face always erode. redilate re-grows afterwards.
static int ss_mask_erode(int nx, int ny, int nz, unsigned char *mmm, int redilate, int NN) {
	long long nxy = (long long)nx * ny, nxyz = nxy * nz;
	unsigned char *nnn = (unsigned char *)calloc((size_t)nxyz, 1);
	int ii, jj, kk;
	int need = (NN >= 3) ? 26 : (NN == 2 ? 18 : 6);
	if (!nnn)
		return -1;
	for (kk = 0; kk < nz; kk++)
		for (jj = 0; jj < ny; jj++)
			for (ii = 0; ii < nx; ii++) {
				long long p = ii + (long long)jj * nx + (long long)kk * nxy;
				if (!mmm[p])
					continue;
				if (ii == 0 || jj == 0 || kk == 0 || ii == nx - 1 || jj == ny - 1 || kk == nz - 1) {
					nnn[p] = 1;
					continue;
				}
				if (ss_count_nbrs(mmm, nx, ny, nz, ii, jj, kk, NN) < need)
					nnn[p] = 1;
			}
	for (long long q = 0; q < nxyz; q++)
		if (nnn[q])
			mmm[q] = 0;
	free(nnn);
	if (redilate)
		return ss_mask_dilate(nx, ny, nz, mmm, 1, NN);
	return 0;
}

// ---------------------------------------------------------------------------
// Cubic warp. Adapted from mri_warp3D_cubic (public domain).
// ---------------------------------------------------------------------------

#define SS_PM1(x) (-(x) * ((x) - 1) * ((x) - 2))
#define SS_P00(x) (3 * ((x) + 1) * ((x) - 1) * ((x) - 2))
#define SS_PP1(x) (-3 * (x) * ((x) + 1) * ((x) - 2))
#define SS_PP2(x) ((x) * ((x) + 1) * ((x) - 1))
#define SS_PFACTOR 4.62962963e-3f // 1/216

static int ss_clampi(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }
// NaN-safe: `!(v >= lo)` is true for NaN, so a NaN maps to lo rather than through the cast.
static float ss_clampf(float v, float lo, float hi) { return !(v >= lo) ? lo : (v > hi ? hi : v); }
static double ss_clampd(double v, double lo, double hi) { return !(v >= lo) ? lo : (v > hi ? hi : v); }

static void ss_warp3D_cubic(const short *src, int nx, int ny, int nz,
		short *dst, int onx, int ony, int onz,
		float ai, float bi, float aj, float bj, float ak, float bk) {
	long long nxy = (long long)nx * ny;
	int oi, oj, ok;
	short smin = 32767, smax = -32768;

	for (long long q = 0; q < nxy * nz; q++) {
		if (src[q] < smin) smin = src[q];
		if (src[q] > smax) smax = src[q];
	}

	for (ok = 0; ok < onz; ok++) {
		// Clamp in FLOAT before the cast, ONLY to make the conversion defined. An
		// out-of-range float->int conversion is UB, and a legal header with a tiny pixdim
		// reaches it (measured: pixdim 1e-9 gives fz ~ -1e11, flagged by UBSan; arm64
		// saturates, wasm32 traps). The bound is +/-1e9, far outside any real grid, because
		// a TIGHT bound would change `dz = fz - kz` and therefore the interpolation weights
		// -- measured, clamping at +/-4 moved three of the five validation masks.
		// ss_clampi below still does the actual index clamping.
		float fz = ss_clampf(ak * ok + bk, -1.0e9f, 1.0e9f);
		int kz = (int)floorf(fz);
		float dz = fz - kz;
		float wzm = SS_PM1(dz), wz0 = SS_P00(dz), wzp = SS_PP1(dz), wzq = SS_PP2(dz);
		for (oj = 0; oj < ony; oj++) {
			float fy = ss_clampf(aj * oj + bj, -1.0e9f, 1.0e9f);
			int jy = (int)floorf(fy);
			float dy = fy - jy;
			float wym = SS_PM1(dy), wy0 = SS_P00(dy), wyp = SS_PP1(dy), wyq = SS_PP2(dy);
			for (oi = 0; oi < onx; oi++) {
				float fx = ss_clampf(ai * oi + bi, -1.0e9f, 1.0e9f);
				int ix = (int)floorf(fx);
				float dx = fx - ix;
				float wxm = SS_PM1(dx), wx0 = SS_P00(dx), wxp = SS_PP1(dx), wxq = SS_PP2(dx);
				float val = 0.0f;
				int dk, dj, di;
				const float wz[4] = {wzm, wz0, wzp, wzq};
				const float wy[4] = {wym, wy0, wyp, wyq};
				const float wx[4] = {wxm, wx0, wxp, wxq};

				// Outside the source volume reads as 0, matching AFNI's default.
				if (fx < -1.0f || fx > nx || fy < -1.0f || fy > ny || fz < -1.0f || fz > nz) {
					dst[(long long)oi + (long long)oj * onx + (long long)ok * onx * ony] = 0;
					continue;
				}
				for (dk = 0; dk < 4; dk++) {
					int kk = ss_clampi(kz - 1 + dk, 0, nz - 1);
					float az = wz[dk];
					if (az == 0.0f)
						continue;
					for (dj = 0; dj < 4; dj++) {
						int jj = ss_clampi(jy - 1 + dj, 0, ny - 1);
						float ay = wy[dj] * az;
						if (ay == 0.0f)
							continue;
						for (di = 0; di < 4; di++) {
							int iii = ss_clampi(ix - 1 + di, 0, nx - 1);
							val += wx[di] * ay * (float)src[iii + (long long)jj * nx + (long long)kk * nxy];
						}
					}
				}
				val *= SS_PFACTOR;
				// Cubic overshoot is clipped to the input data range, as AFNI does.
				if (val < smin) val = smin;
				if (val > smax) val = smax;
				dst[(long long)oi + (long long)oj * onx + (long long)ok * onx * ony] = ss_shortize(val);
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Orientation: snap the world transform to the nearest storage axes and build the
// RAI permutation. AFNI reads its orientation codes from the dataset header; niimath
// derives them from sform/qform the same way -moco does (nearest-axis snapping).
// ---------------------------------------------------------------------------

static int ss_orient_codes(const nifti_image *nim, int *fi, int *fj, int *fk, double *worst_cos) {
	// Column c of the transform gives the world displacement per step of voxel axis c.
	// Assign each voxel axis to the world axis it moves along most.
	double m[3][3];
	int c, r, used[3] = {0, 0, 0};
	int axis_of[3], sign_of[3];
	const nifti_dmat44 *M;

	if (nim->sform_code > 0)
		M = &nim->sto_xyz;
	else if (nim->qform_code > 0)
		M = &nim->qto_xyz;
	else
		return 1;

	double worst = 1.0;

	for (r = 0; r < 3; r++)
		for (c = 0; c < 3; c++)
			m[r][c] = M->m[r][c];

	// Greedy nearest-axis assignment. It is PROVABLY optimal wherever "nearest axis" means
	// anything: rows of an orthogonal matrix are unit vectors, so two columns cannot both have
	// a component > 1/sqrt(2) in the same row, and therefore whenever every voxel axis is
	// within 45 degrees of a distinct world axis that assignment is unique and greedy is
	// forced onto it. (Checked against an exhaustive 6-permutation optimum over 20M random
	// rotations: every disagreement had an axis more than 45 degrees off, the closest at
	// 45.004.) Beyond that both answers are equally arbitrary -- which is why the caller
	// warns; see the obliquity note there.
	for (c = 0; c < 3; c++) {
		double best = -1.0;
		int br = -1;
		for (r = 0; r < 3; r++) {
			double a = fabs(m[r][c]);
			if (!used[r] && a > best) {
				best = a;
				br = r;
			}
		}
		if (br < 0 || best <= 0.0)
			return 1;
		// Record the worst DIRECTION COSINE, i.e. the column's dominant component divided by
		// the column norm. `best` alone is a raw affine magnitude and scales with the voxel
		// size, so comparing it to cos(30) made a cardinal 0.5 mm image look 60 degrees
		// oblique and let a strongly oblique 4 mm image pass unremarked.
		{
			double cn = sqrt(m[0][c] * m[0][c] + m[1][c] * m[1][c] + m[2][c] * m[2][c]);
			double cosang = (cn > 0.0) ? best / cn : 0.0;
			if (cosang < worst)
				worst = cosang;
		}
		used[br] = 1;
		axis_of[c] = br;                       // 0=x(R+),1=y(A+),2=z(S+) in NIfTI world
		sign_of[c] = (m[br][c] > 0) ? 1 : -1;
	}

	// AFNI's RAI target: +i is R->L, +j is A->P, +k is I->S. NIfTI world +x is toward
	// R and +y toward A, so RAI's i and j run OPPOSITE to NIfTI x and y, while k
	// agrees with z.
	{
		int want_sign[3] = {-1, -1, +1}; // for world axis x,y,z to become RAI i,j,k
		int code[3];
		for (c = 0; c < 3; c++) {
			int wa = axis_of[c];
			int s = sign_of[c] * want_sign[wa];
			code[c] = (wa + 1) * s; // +/-1 = R/L run, +/-2 = A/P run, +/-3 = I/S run
		}
		// code[c] says which RAI output axis voxel axis c feeds, and in which direction.
		// AFNI's flip encoding wants, for each of its ii/jj/kk, the SOURCE axis.
		*fi = *fj = *fk = 0;
		for (c = 0; c < 3; c++) {
			int tgt = abs(code[c]) - 1; // 0=i(RAI x),1=j,2=k
			int v = (c + 1) * ((code[c] > 0) ? 1 : -1);
			if (tgt == 0) *fi = v;
			else if (tgt == 1) *fj = v;
			else *fk = v;
		}
		if (*fi == 0 || *fj == 0 || *fk == 0)
			return 1;
	}
	if (worst_cos) *worst_cos = worst;
	return 0;
}

// Permute/flip the source float volume into an RAI-ordered short volume.
static short *ss_to_rai_short(const nifti_image *nim, int in_datatype, int fi, int fj, int fk,
		int *rnx, int *rny, int *rnz, float *rdx, float *rdy, float *rdz) {
	const float *in = (const float *)nim->data;
	int snx = (int)nim->nx, sny = (int)nim->ny, snz = (int)nim->nz;
	int sdim[3] = {snx, sny, snz};
	// Voxel sizes normalised to MILLIMETRES. Everything downstream -- SS_DXYZ, SS_ZHEIGHT,
	// SS_CM_DEPTH, the head-extent clipping, the centre-of-mass depth -- is in mm, so a
	// legal metre- or micron-encoded header would otherwise be treated as a radically
	// different physical size (a 1 m head, or a 0.001 mm one). This is the same class of bug
	// AGENTS.md records as still open for the registration engine; here it costs one
	// multiply, so there is no reason to carry it. mm headers scale by exactly 1.0 and are
	// bit-identical to before.
	// Deliberately local rather than calling core.c's xyz_units_to_mm(): skullstrip.c is a
	// self-contained TU, linked on its own by the mesh selftest, the CMake selftest and the
	// benchmark drivers, and none of those want core.c's dependency graph for three lines.
	// Keep the two in agreement -- same switch, same mm default for unspecified units.
	const float uscl = (nim->xyz_units == NIFTI_UNITS_METER)  ? 1000.0f
	                 : (nim->xyz_units == NIFTI_UNITS_MICRON) ? 0.001f
	                                                          : 1.0f;
	// Spacing comes from PIXDIM, orientation from the coded sform/qform. That is two header
	// authorities, and it is DELIBERATE: it is what the reference does. AFNI's own NIfTI
	// reader (thd_niftiread.c:495-503) sets the dataset's voxel dimensions from
	// `nim->pixdim[1..3]`, falling back to 1 when non-positive, and applies exactly the
	// metre/micron conversion below -- it never uses the affine column norms for spacing.
	// Deriving spacing from the selected transform instead was tried and REJECTED: on the
	// five validation images pixdim and the column norms agree to 1e-7 relative, so it
	// changed three of five masks purely by amplifying float noise through the stage loop's
	// discrete convergence test, and moved Dice against AFNI down (0.9697 -> 0.9685).
	// A header where the two genuinely disagree is malformed; warn rather than silently
	// pick one, since the whole pipeline is a port and pixdim is what the port must follow.
	float sdel[3] = {(float)fabs(nim->dx) * uscl, (float)fabs(nim->dy) * uscl,
	                 (float)fabs(nim->dz) * uscl};
	{
		const nifti_dmat44 *M = (nim->sform_code > 0) ? &nim->sto_xyz
		                      : (nim->qform_code > 0) ? &nim->qto_xyz : NULL;
		const float pd[3] = {(float)fabs(nim->dx), (float)fabs(nim->dy), (float)fabs(nim->dz)};
		if (M)
			for (int c = 0; c < 3; c++) {
				double cn = sqrt((double)M->m[0][c] * M->m[0][c] +
				                 (double)M->m[1][c] * M->m[1][c] +
				                 (double)M->m[2][c] * M->m[2][c]);
				if (cn > 1e-6 && pd[c] > 1e-6 && fabs(cn - pd[c]) / pd[c] > 0.01)
					printfx("skullstrip: WARNING: pixdim[%d] (%g) disagrees with the coded "
							"transform's column norm (%g) by more than 1%%. Orientation comes "
							"from the transform and spacing from pixdim, matching AFNI, so a "
							"header where these differ will be normalised inconsistently.\n",
							c + 1, (double)pd[c], cn);
			}
	}
	int src_of[3] = {abs(fi) - 1, abs(fj) - 1, abs(fk) - 1};
	int rev[3] = {fi < 0, fj < 0, fk < 0};
	int onx = sdim[src_of[0]], ony = sdim[src_of[1]], onz = sdim[src_of[2]];
	long long n = (long long)snx * sny * snz;
	short *out;
	float maxabs = 0.0f, scale;
	int integral = 1;
	long long q;

	for (q = 0; q < n; q++) {
		float v = in[q];
		if (!isfinite(v))
			continue;
		float a = fabsf(v);
		if (a > maxabs)
			maxabs = a;
		if (integral && v != floorf(v))
			integral = 0;
	}
	if (maxabs <= 0.0f)
		return NULL;

	// AFNI branches on the STORED datatype: an integer dataset is copied verbatim, anything
	// else is rescaled to fill the short range. That choice changes the histogram and hence
	// every clip level downstream, so it has to be right.
	//
	// niimath promotes the image to float32 before any op runs, so the caller passes the
	// PRE-PROMOTION datatype in. Inferring it from the values instead -- which is what this
	// did until an audit called it out three times -- is single-voxel fragile: a stored FLOAT
	// image whose samples all happen to be integral took the integer path, and adding 0.5 to
	// one background voxel flipped the whole volume's normalisation (measured: +6.7% retained
	// volume, Dice 0.9659 on coarsely quantised float data).
	//
	// The value inspection survives ONLY as a fallback for in_datatype == DT_NONE, which the
	// shipped dispatch never passes; it exists so the standalone harnesses, which have no
	// niimath header to consult, still behave sensibly.
	{
		int stored_is_int;
		switch (in_datatype) {
			case DT_UINT8: case DT_INT16: case DT_UINT16:
			case DT_INT32: case DT_UINT32: case DT_INT64: case DT_UINT64:
				stored_is_int = 1;
				break;
			case DT_FLOAT32: case DT_FLOAT64:
				stored_is_int = 0;
				break;
			default:   // DT_NONE / anything unexpected: fall back to inspecting the values
				stored_is_int = (integral && maxabs <= 32767.0f);
				break;
		}
		// A scaled integer dataset is not integral once niimath applies scl_slope, and AFNI
		// sees the scaling too, so the verbatim path additionally requires the values to
		// still fit the short range.
		if (stored_is_int && maxabs <= 32767.0f)
			scale = 1.0f;
		else
			scale = 32767.0f / maxabs;
	}

	out = (short *)malloc(sizeof(short) * (size_t)n);
	if (!out)
		return NULL;

	for (int ok = 0; ok < onz; ok++)
		for (int oj = 0; oj < ony; oj++)
			for (int oi = 0; oi < onx; oi++) {
				int o[3] = {oi, oj, ok};
				int s[3];
				long long si, di;
				for (int a = 0; a < 3; a++) {
					int sa = src_of[a];
					s[sa] = rev[a] ? (sdim[sa] - 1 - o[a]) : o[a];
				}
				si = (long long)s[0] + (long long)s[1] * snx + (long long)s[2] * snx * sny;
				di = (long long)oi + (long long)oj * onx + (long long)ok * onx * ony;
				{
					float v = in[si];
					if (!isfinite(v))
						v = 0.0f;
					out[di] = ss_shortize(v * scale);
				}
			}

	*rnx = onx; *rny = ony; *rnz = onz;
	*rdx = sdel[src_of[0]]; *rdy = sdel[src_of[1]]; *rdz = sdel[src_of[2]];
	if (*rdx == 0.0f) *rdx = 1.0f;
	if (*rdy == 0.0f) *rdy = 1.0f;
	if (*rdz == 0.0f) *rdz = 1.0f;
	return out;
}

// ---------------------------------------------------------------------------
// The pipeline. Adapted from mri_brainormalize().
// ---------------------------------------------------------------------------

void ss_norm_free(ss_norm *n) {
	if (!n)
		return;
	free(n->vol);
	memset(n, 0, sizeof(*n));
}

int ss_restore(const ss_norm *n, const unsigned char *wvol, const nifti_image *nim,
		float *dst, int nearest) {
	int snx, sny, snz;
	int sdim[3], src_of[3], rev[3];
	long long nxy = (long long)SS_NX * SS_NY;

	if (!n || !wvol || !nim || !dst || !n->vol)
		return 1;
	if (n->ai == 0.0f || n->aj == 0.0f || n->ak == 0.0f)
		return 1;
	snx = (int)nim->nx; sny = (int)nim->ny; snz = (int)nim->nz;
	sdim[0] = snx; sdim[1] = sny; sdim[2] = snz;
	src_of[0] = abs(n->fi) - 1; src_of[1] = abs(n->fj) - 1; src_of[2] = abs(n->fk) - 1;
	rev[0] = n->fi < 0; rev[1] = n->fj < 0; rev[2] = n->fk < 0;

	// Walk the SOURCE grid in its own storage order, so the output needs no
	// reordering afterwards and the caller's header stays valid untouched.
	for (int sk = 0; sk < snz; sk++)
		for (int sj = 0; sj < sny; sj++)
			for (int si = 0; si < snx; si++) {
				int s[3] = {si, sj, sk};
				float r[3]; // RAI source index
				float o[3]; // working-grid index
				long long di = (long long)si + (long long)sj * snx + (long long)sk * snx * sny;
				for (int a = 0; a < 3; a++) {
					int sa = src_of[a];
					r[a] = rev[a] ? (float)(sdim[sa] - 1 - s[sa]) : (float)s[sa];
				}
				o[0] = (r[0] - n->bi) / n->ai;
				o[1] = (r[1] - n->bj) / n->aj;
				o[2] = (r[2] - n->bk) / n->ak;

				// Clamp before converting. The range check below is unchanged and still does
				// the real work; this only stops the CONVERSION being undefined, which a
				// legal extreme voxel size can otherwise reach (o[] = (r - b)/a with a
				// tiny). The bound is far outside any real grid so no in-range coordinate
				// is touched; NaN maps to -1e9 and is rejected by the same check.
				if (nearest) {
					int oi = (int)lrintf(ss_clampf(o[0], -1.0e9f, 1.0e9f));
					int oj = (int)lrintf(ss_clampf(o[1], -1.0e9f, 1.0e9f));
					int ok = (int)lrintf(ss_clampf(o[2], -1.0e9f, 1.0e9f));
					if (oi < 0 || oi >= SS_NX || oj < 0 || oj >= SS_NY || ok < 0 || ok >= SS_NZ)
						dst[di] = 0.0f;
					else
						dst[di] = (float)wvol[oi + (long long)oj * SS_NX + (long long)ok * nxy];
				} else {
					int oi = (int)floorf(ss_clampf(o[0], -1.0e9f, 1.0e9f));
					int oj = (int)floorf(ss_clampf(o[1], -1.0e9f, 1.0e9f));
					int ok = (int)floorf(ss_clampf(o[2], -1.0e9f, 1.0e9f));
					float fx = o[0] - oi, fy = o[1] - oj, fz = o[2] - ok, acc = 0.0f;
					// The `!(v >= lo)` clamps above give a defined INDEX for a NaN coordinate,
					// but fx/fy/fz are computed from the original o[] and stay NaN, and the
					// range test below is false for NaN. Reject non-finite explicitly, or a
					// NaN weight propagates into the interpolated value. Written as a
					// positive test so NaN fails it.
					if (!(o[0] >= -0.5f && o[0] <= SS_NX - 0.5f &&
							o[1] >= -0.5f && o[1] <= SS_NY - 0.5f &&
							o[2] >= -0.5f && o[2] <= SS_NZ - 0.5f)) {
						dst[di] = 0.0f;
						continue;
					}
					for (int dk = 0; dk < 2; dk++)
						for (int dj = 0; dj < 2; dj++)
							for (int dd = 0; dd < 2; dd++) {
								int xi = ss_clampi(oi + dd, 0, SS_NX - 1);
								int xj = ss_clampi(oj + dj, 0, SS_NY - 1);
								int xk = ss_clampi(ok + dk, 0, SS_NZ - 1);
								float w = (dd ? fx : 1.0f - fx) * (dj ? fy : 1.0f - fy) *
								          (dk ? fz : 1.0f - fz);
								acc += w * (float)wvol[xi + (long long)xj * SS_NX + (long long)xk * nxy];
							}
					dst[di] = acc;
				}
			}
	return 0;
}

int ss_normalize(const nifti_image *nim, int in_datatype, ss_norm *out) {
	short *sar = NULL, *tar = NULL;
	unsigned char *mask = NULL;
	int nx, ny, nz;
	long long nxy, nxyz;
	float dx, dy, dz;
	int fi, fj, fk;
	int ii, jj, kk, ktop, kbot;
	double icm = 0, jcm = 0, kcm = 0, sum = 0;
	long long cnt;
	int rc = 1;

	if (!nim || !out || !nim->data)
		return 1;
	memset(out, 0, sizeof(*out));
	if (nim->datatype != DT_FLOAT32)
		return 1;
	if (nim->nx < 16 || nim->ny < 16 || nim->nz < 16) {
		printfx("skullstrip: needs at least 16 voxels in every dimension\n");
		return 1;
	}

	{
		// Warn on strong obliquity. The nearest-axis assignment is provably optimal only
		// while every voxel axis is within 45 degrees of a distinct world axis (see
		// ss_orient_codes); past that, the permutation is arbitrary and this whole pipeline
		// -- which ignores obliquity and works on storage axes -- silently returns a
		// plausible but MIS-ORIENTED mask. Measured: a 60-degree oblique sform on T1w runs
		// to completion at rc=0. Warn rather than reject, matching -moco's documented
		// "obliquity ignored" convention, so an oblique-but-usable image still works.
		double worst = 1.0;
		if (ss_orient_codes(nim, &fi, &fj, &fk, &worst)) {
			printfx("skullstrip: cannot determine image orientation (no usable sform/qform)\n");
			return 1;
		}
		if (worst < 0.866) // cos(30 degrees)
			printfx("skullstrip: WARNING: strongly oblique transform (worst axis %.1f degrees "
					"off a cardinal direction). This operation ignores obliquity and works on "
					"storage axes, so the mask may be misoriented. Consider -ras first.\n",
					acos(worst < 1.0 ? worst : 1.0) * 180.0 / 3.14159265358979);
	}
	SSV("skullstrip: RAI flip codes = [%d %d %d]\n", fi, fj, fk);

	sar = ss_to_rai_short(nim, in_datatype, fi, fj, fk, &nx, &ny, &nz, &dx, &dy, &dz);
	if (!sar) {
		printfx("skullstrip: degenerate intensity range or out of memory\n");
		return 1;
	}
	nxy = (long long)nx * ny;
	nxyz = nxy * nz;
	SSV("skullstrip: RAI grid %dx%dx%d  d=[%g %g %g]\n", nx, ny, nz, dx, dy, dz);

	// (b) binary mask from spatially varying cliplevels
	{
		ss_clipvec bvec = ss_octant_clips(sar, nx, ny, nz, 0.40f);
		if (bvec.c000 < 0.0f) {
			printfx("skullstrip: could not determine clip levels\n");
			goto done;
		}
		mask = (unsigned char *)malloc((size_t)nxyz);
		if (!mask)
			goto done;
		for (long long ijk = 0, k2 = 0; k2 < nz; k2++)
			for (int j2 = 0; j2 < ny; j2++)
				for (int i2 = 0; i2 < nx; i2++, ijk++)
					mask[ijk] = (sar[ijk] >= ss_pointclip(i2, j2, (int)k2, &bvec));
		// EVERY one of these can fail on allocation, and a silent skip is not a small
		// error: measured, a failed 16 KB stack alloc inside ss_mask_clust turned "keep the
		// largest cluster" into a no-op and published a mask 46% too large at exit 0. They
		// fail OPEN -- the mask gets bigger, i.e. less is stripped -- so nothing downstream
		// notices. Propagate every one.
		if (ss_clustedit3D(nx, ny, nz, mask, (long long)rint(0.02 * (double)nxyz)))
			goto oom;
		SSV("skullstrip: clip mask %lld voxels\n", ss_mask_count(nxyz, mask));
	}

	if (ss_mask_fillin_once(nx, ny, nz, mask, 2) ||
			ss_mask_dilate(nx, ny, nz, mask, 5, 2) ||
			ss_mask_dilate(nx, ny, nz, mask, 5, 2))
		goto oom;

	cnt = ss_mask_count(nxyz, mask);
	SSV("skullstrip: filled mask %lld voxels\n", cnt);
	if (cnt <= 999) {
		printfx("skullstrip: too little signal to form a head mask\n");
		goto done;
	}

	// (c) superior clip: topmost slice with three slices in a row carrying enough
	{
		long long *zc = (long long *)malloc(sizeof(long long) * nz);
		long long z1 = (long long)(0.010 * nxy), z2 = (long long)(0.015 * nxy),
		          z3 = (long long)(0.020 * nxy);
		if (!zc)
			goto done;
		for (kk = nz - 1; kk >= 0; kk--)
			zc[kk] = ss_mask_count(nxy, mask + (long long)kk * nxy);
		for (kk = nz - 1; kk > 2; kk--)
			if (zc[kk] >= z1 && zc[kk - 1] >= z2 && zc[kk - 2] >= z3)
				break;
		free(zc);
		if (kk <= 2) {
			printfx("skullstrip: could not find the top of the head\n");
			goto done;
		}
		ktop = kk;
		if (ktop < nz - 1)
			memset(mask + (long long)(ktop + 1) * nxy, 0, (size_t)(nxy * (nz - 1 - ktop)));
		// Clamp in double before the cast: a legal but absurd pixdim (measured: 1e-9) makes
		// SS_ZHEIGHT/dz ~ 1e11 and the conversion undefined. The slice index is clamped
		// into range immediately below either way.
		jj = (int)ss_clampf(ktop - SS_ZHEIGHT / dz, -1.0e9f, 1.0e9f);
		kbot = jj;
		if (jj >= 0)
			memset(mask, 0, (size_t)(nxy * (jj + 1)));
		SSV("skullstrip: top clip above slice %d, bot clip below slice %d\n", ktop, jj);
	}

	cnt = ss_mask_count(nxyz, mask);
	if (cnt <= 999) {
		printfx("skullstrip: head mask collapsed after clipping\n");
		goto done;
	}
	out->support = (int)cnt;

	// (d) apply mask. NOTE: this does NOT remove negatives, despite what this comment used to
	// claim. The clip level is only positive in the non-degenerate case, and the mask has since
	// been GROWN by ss_mask_fillin_once + two ss_mask_dilate passes, which pull negative voxels
	// back in. The working volume is clamped once after the warp instead -- see that comment.
	for (long long q = 0; q < nxyz; q++)
		if (!mask[q])
			sar[q] = 0;
	free(mask);
	mask = NULL;

	// Centre of mass of the masked image, in RAI source voxel indices. Measured over
	// the top SS_CM_DEPTH mm only (the CMTOP branch -- see the note in skullstrip.h),
	// NOT the whole volume: neck and shoulders would drag the centre inferior.
	{
		int cm_bot = (int)rint(ss_clampd((double)ktop - SS_CM_DEPTH / dz, -1.0e9, 1.0e9));
		if (cm_bot < 0)
			cm_bot = 0;
		for (int k2 = cm_bot; k2 <= ktop; k2++)
			for (jj = 0; jj < ny; jj++)
				for (ii = 0; ii < nx; ii++) {
					long long ijk = ii + (long long)jj * nx + (long long)k2 * nxy;
					float val = (float)sar[ijk];
					sum += val;
					icm += val * ii;
					jcm += val * jj;
					kcm += val * k2;
				}
		SSV("skullstrip: CM measured over slices %d..%d\n", cm_bot, ktop);
	}
	if (sum == 0.0) {
		printfx("skullstrip: empty image after masking\n");
		goto done;
	}
	out->icm = (float)(icm / sum);
	out->jcm = (float)(jcm / sum);
	out->kcm = (float)(kcm / sum);
	out->ktop = ktop;
	out->kbot = kbot;
	SSV("skullstrip: CM = [%g %g %g] (RAI voxels)\n", out->icm, out->jcm, out->kcm);

	// (e) resample onto the fixed grid with the CM placed at (0,20,0) mm
	out->ai = SS_DXYZ / dx;
	out->bi = out->icm - (SS_XCM - SS_XORG) / dx;
	out->aj = SS_DXYZ / dy;
	out->bj = out->jcm - (SS_YCM - SS_YORG) / dy;
	out->ak = SS_DXYZ / dz;
	out->bk = out->kcm - (SS_ZCM - SS_ZORG) / dz;
	SSV("skullstrip: a=[%g %g %g] b=[%g %g %g]\n", out->ai, out->aj, out->ak,
			out->bi, out->bj, out->bk);

	tar = (short *)malloc(sizeof(short) * SS_NX * SS_NY * SS_NZ);
	if (!tar)
		goto done;
	ss_warp3D_cubic(sar, nx, ny, nz, tar, SS_NX, SS_NY, SS_NZ,
			out->ai, out->bi, out->aj, out->bj, out->ak, out->bk);
	// Clamp negatives ONCE, here, where the working volume is established. The head mask is
	// GROWN by ss_mask_fillin_once + two ss_mask_dilate passes AFTER thresholding, so negative
	// voxels get pulled back in and survive; the comment at the threshold step claiming it
	// "also removes negatives, since the clip level is positive" stops being true once that
	// clip level is 0. Downstream, hist[tar[q]]++ (below) indexes a calloc(32768) with a short,
	// so a single negative is a heap underflow WRITE -- measured at 174,812 such writes on a
	// 48^3 int16 fixture whose background is -1. That path was only protected by a clamp nested
	// inside `if (bvec.c000 > 0.0f)`, which the degenerate second-pass octant cliplevel skips.
	// Provably not a behaviour change: ss_octant_clips ignores non-positive voxels and the
	// uniformize branch already maps negatives to 0; verified byte-identical on a real T1.
	for (long long q = 0; q < (long long)SS_NX * SS_NY * SS_NZ; q++)
		if (tar[q] < 0)
			tar[q] = 0;
	free(sar);
	sar = NULL;

	nx = SS_NX; ny = SS_NY; nz = SS_NZ;
	nxy = (long long)nx * ny;
	nxyz = nxy * nz;

	// rescale to partially uniformize: each voxel divided by its local cliplevel
	{
		ss_clipvec bvec = ss_octant_clips(tar, nx, ny, nz, 0.40f);
		// c000 < 0 is the ALLOCATION-FAILURE flag, not "no uniformization needed". Testing
		// `> 0` lumped the two together and silently skipped this whole stage on OOM, then
		// reported success -- the same sentinel read two different ways by two callers.
		if (bvec.c000 < 0.0f) {
			printfx("skullstrip: could not determine clip levels (second pass)\n");
			goto done;
		}
		if (bvec.c000 > 0.0f) {
			for (long long ijk = 0, k2 = 0; k2 < nz; k2++)
				for (jj = 0; jj < ny; jj++)
					for (ii = 0; ii < nx; ii++, ijk++) {
						float bval = ss_pointclip(ii, jj, (int)k2, &bvec);
						float sv = (bval != 0.0f) ? 1000.0f * tar[ijk] / bval : 0.0f;
						short s = ss_shortize(sv);
						tar[ijk] = (s < 0) ? 0 : s;
					}
		}
	}

	// REMASK: threshold at 40% of the median, erode, keep largest cluster, fill holes
	{
		int *hist = (int *)calloc(32768, sizeof(int));
		int sbot, stop2;
		long long nmask, nhalf;
		if (!hist)
			goto done;
		for (long long q = 0; q < nxyz; q++)
			hist[tar[q]]++;
		for (sbot = 1; sbot < 32768 && hist[sbot] == 0; sbot++)
			;
		for (stop2 = 32767; stop2 > sbot && hist[stop2] == 0; stop2--)
			;
		if (sbot < 32768 && stop2 > sbot) {
			int cbot;
			nmask = 0;
			for (ii = sbot; ii <= stop2; ii++)
				nmask += hist[ii];
			nhalf = nmask / 2;
			nmask = 0;
			for (ii = sbot; ii <= stop2 && nmask < nhalf; ii++)
				nmask += hist[ii];
			cbot = (int)(0.40 * ii);
			mask = (unsigned char *)malloc((size_t)nxyz);
			if (!mask) {
				free(hist);
				goto done;
			}
			for (long long q = 0; q < nxyz; q++)
				mask[q] = (tar[q] > cbot);
			// free(hist) before each exit: `hist` is local to this block and the common
			// cleanup at done: does not know about it. These two gotos were added with the
			// OOM propagation and leaked it (caught in audit).
			if (ss_mask_erode(nx, ny, nz, mask, 1, 2) ||
					ss_mask_clust(nx, ny, nz, mask)) {
				free(hist);
				goto oom;
			}
			for (long long q = 0; q < nxyz; q++)
				mask[q] = !mask[q];
			if (ss_mask_clust(nx, ny, nz, mask)) {
				free(hist);
				goto oom;
			}
			for (long long q = 0; q < nxyz; q++)
				mask[q] = !mask[q];
			for (long long q = 0; q < nxyz; q++)
				if (!mask[q])
					tar[q] = 0;
			free(mask);
			mask = NULL;
			SSV("skullstrip: remask threshold %d\n", cbot);
		}
		free(hist);
	}

	// clip the top 1% -- this ALWAYS runs; AFNI's `else` is #if 0'd out
	{
		int *hist = (int *)calloc(32768, sizeof(int));
		long long tot = 0, acc = 0;
		if (!hist)
			goto done;
		for (long long q = 0; q < nxyz; q++)
			hist[tar[q]]++;
		for (ii = 0; ii < 32767; ii++)
			tot += hist[ii];
		tot = (long long)(0.01 * (double)tot);
		for (ii = 32767; ii > 0 && acc < tot; ii--)
			acc += hist[ii];
		out->clip99 = ii;
		for (long long q = 0; q < nxyz; q++)
			if (tar[q] > ii)
				tar[q] = (short)ii;
		free(hist);
		SSV("skullstrip: 99%% clip at %d\n", out->clip99);
	}

	// convert to the byte working volume AFNI actually deforms against
	out->vol = (unsigned char *)malloc((size_t)nxyz);
	if (!out->vol)
		goto done;
	{
		int mx = 0;
		for (long long q = 0; q < nxyz; q++)
			if (tar[q] > mx)
				mx = tar[q];
		if (mx > 255) {
			float fac = 255.0f / mx;
			SSV("skullstrip: scaling to byte by %g\n", fac);
			for (long long q = 0; q < nxyz; q++)
				out->vol[q] = (unsigned char)(fac * tar[q] + 0.49f);
		} else {
			for (long long q = 0; q < nxyz; q++)
				out->vol[q] = (unsigned char)tar[q];
		}
	}

	out->fi = fi; out->fj = fj; out->fk = fk;
	rc = 0;
	goto done;
oom:
	// rc is still 1 here: it is only cleared on the success path above.
	printfx("skullstrip: out of memory during spatial normalisation\n");
done:
	free(sar);
	free(tar);
	free(mask);
	if (rc)
		ss_norm_free(out);
	return rc;
}

// ===========================================================================
// Surface primitives. CLEAN-ROOM -- ordinary computational geometry, written before the
// provenance reversal and owing nothing to SUMA_BrainWrap.c. Gated on analytic answers:
// Euler characteristic, sphere area/volume convergence, exact cube rasterisation.
//
// ONE deliberate exception lives in this section: ss_afni_self_intersect is a faithful port
// of AFNI's SUMA_isSelfIntersect, quirk included, because it decides the retry. It is
// labelled as such at its definition.
// ===========================================================================

void ss_mesh_free(ss_mesh *m) {
	if (!m)
		return;
	free(m->v); free(m->t); free(m->nbr); free(m->nbr_off); free(m->nrm);
	memset(m, 0, sizeof(*m));
}

// Index of the point t/ld along base edge (u,v). The edge is stored once, in the
// canonical low->high direction, so a face traversing it backwards asks for ld-t.
// This is what makes the two faces sharing an edge agree on vertex identity EXACTLY,
// with no coordinate comparison and therefore no tolerance to tune.
static int ss_edge_vertex(int emap[12][12], int off_edge, int ld, int u, int v, int t) {
	int e = emap[u][v];
	int tt = (u < v) ? t : (ld - t);
	return off_edge + e * (ld - 1) + (tt - 1);
}

// Global vertex index of lattice point (i,j) on base face f, i+j <= ld.
static int ss_lattice_index(int bf[20][3], int emap[12][12], int off_edge, int off_face,
		int per_face, int ld, int f, int i, int j) {
	int a = bf[f][0], b = bf[f][1], c = bf[f][2];
	if (i == 0 && j == 0) return a;
	if (i == ld && j == 0) return b;
	if (i == 0 && j == ld) return c;
	if (j == 0) return ss_edge_vertex(emap, off_edge, ld, a, b, i);
	if (i == 0) return ss_edge_vertex(emap, off_edge, ld, a, c, j);
	if (i + j == ld) return ss_edge_vertex(emap, off_edge, ld, b, c, j);
	{
		int idx = 0;
		for (int jj = 1; jj < j; jj++)
			idx += (ld - 1 - jj);
		idx += (i - 1);
		return off_face + f * per_face + idx;
	}
}

// Insertion sort, NOT qsort: these lists are 5-12 entries and there are nv of them.
// emscripten/musl qsort dispatches its comparator through call_indirect, which the
// repo has measured at roughly 100x native -- see the WASM pitfalls in AGENTS.md.
static void ss_isort(int *a, int n) {
	for (int i = 1; i < n; i++) {
		int key = a[i], j = i - 1;
		while (j >= 0 && a[j] > key) {
			a[j + 1] = a[j];
			j--;
		}
		a[j + 1] = key;
	}
}

// Neighbour lists from the triangle list: symmetric and duplicate-free by
// construction (sort + unique per vertex).
static int ss_build_adjacency(ss_mesh *m) {
	int nv = m->nv, nt = m->nt;
	int *deg = (int *)calloc((size_t)nv, sizeof(int));
	int *off = (int *)malloc(sizeof(int) * (size_t)(nv + 1));
	int *tmp = NULL, *fill = NULL;
	long long total = 0;

	if (!deg || !off) {
		free(deg); free(off);
		return 1;
	}
	for (int t = 0; t < nt; t++)
		for (int e = 0; e < 3; e++)
			deg[m->t[3 * t + e]] += 2;
	off[0] = 0;
	for (int i = 0; i < nv; i++)
		off[i + 1] = off[i] + deg[i];
	total = off[nv];
	tmp = (int *)malloc(sizeof(int) * (size_t)total);
	fill = (int *)calloc((size_t)nv, sizeof(int));
	if (!tmp || !fill) {
		free(deg); free(off); free(tmp); free(fill);
		return 1;
	}
	for (int t = 0; t < nt; t++) {
		int a = m->t[3 * t + 0], b = m->t[3 * t + 1], c = m->t[3 * t + 2];
		int tri[3] = {a, b, c};
		for (int e = 0; e < 3; e++) {
			int v0 = tri[e], v1 = tri[(e + 1) % 3], v2 = tri[(e + 2) % 3];
			tmp[off[v0] + fill[v0]++] = v1;
			tmp[off[v0] + fill[v0]++] = v2;
		}
	}
	// Sort + unique each list, compacting in place.
	m->nbr_off = (int *)malloc(sizeof(int) * (size_t)(nv + 1));
	if (!m->nbr_off) {
		free(deg); free(off); free(tmp); free(fill);
		return 1;
	}
	{
		int w = 0;
		m->nbr_off[0] = 0;
		for (int i = 0; i < nv; i++) {
			int lo = off[i], n = fill[i];
			ss_isort(tmp + lo, n);
			for (int k = 0; k < n; k++) {
				if (k > 0 && tmp[lo + k] == tmp[lo + k - 1])
					continue;
				tmp[w++] = tmp[lo + k];
			}
			m->nbr_off[i + 1] = w;
		}
		m->nbr = (int *)malloc(sizeof(int) * (size_t)(w > 0 ? w : 1));
		if (!m->nbr) {
			free(deg); free(off); free(tmp); free(fill);
			return 1;
		}
		memcpy(m->nbr, tmp, sizeof(int) * (size_t)w);
	}
	free(deg); free(off); free(tmp); free(fill);
	return 0;
}

// The 12 icosahedron corners: cyclic permutations of (0, +/-1, +/-phi).
static const float SS_PHI = 1.6180339887498949f;

// Build the base icosahedron with every face wound counter-clockwise seen from
// outside, so the surface normal points away from the centre.
static void ss_base_icosa(float v[12][3], int f[20][3]) {
	const float p = SS_PHI;
	float src[12][3] = {
		{0, 1, p}, {0, -1, p}, {0, 1, -p}, {0, -1, -p},
		{1, p, 0}, {-1, p, 0}, {1, -p, 0}, {-1, -p, 0},
		{p, 0, 1}, {-p, 0, 1}, {p, 0, -1}, {-p, 0, -1}};
	int faces[20][3] = {
		{0, 1, 8}, {0, 8, 4}, {0, 4, 5}, {0, 5, 9}, {0, 9, 1},
		{1, 9, 7}, {1, 7, 6}, {1, 6, 8}, {8, 6, 10}, {8, 10, 4},
		{4, 10, 2}, {4, 2, 5}, {5, 2, 11}, {5, 11, 9}, {9, 11, 7},
		{3, 7, 11}, {3, 11, 2}, {3, 2, 10}, {3, 10, 6}, {3, 6, 7}};
	memcpy(v, src, sizeof(src));
	memcpy(f, faces, sizeof(faces));
}

int ss_icosphere(int ld, float radius, const float centre[3], ss_mesh *m) {
	float bv[12][3];
	int bf[20][3];
	int nv, nt, ne = 30;
	int edge_a[30], edge_b[30], nedge = 0;
	int emap[12][12];
	int ti = 0;
	float cx = centre ? centre[0] : 0.0f, cy = centre ? centre[1] : 0.0f,
	      cz = centre ? centre[2] : 0.0f;

	// ld is bounded as well as positive: nt = 20*ld*ld overflows int at ld >= 10362, and the
	// consistency check below compares two EQUALLY overflowed expressions, so it passes and
	// the mallocs are then sized from a wrapped count. Not reachable from the CLI (the only
	// caller passes SS_LD_NOEDGE) but this is exported in the header and used by the tests.
	if (!m || ld < 1 || ld > 4096 || radius <= 0.0f)
		return 1;
	memset(m, 0, sizeof(*m));
	ss_base_icosa(bv, bf);

	// Canonical edge table.
	for (int a = 0; a < 12; a++)
		for (int b = 0; b < 12; b++)
			emap[a][b] = -1;
	for (int f = 0; f < 20; f++)
		for (int e = 0; e < 3; e++) {
			int a = bf[f][e], b = bf[f][(e + 1) % 3];
			int lo = a < b ? a : b, hi = a < b ? b : a;
			if (emap[lo][hi] < 0) {
				emap[lo][hi] = emap[hi][lo] = nedge;
				edge_a[nedge] = lo;
				edge_b[nedge] = hi;
				nedge++;
			}
		}
	if (nedge != ne)
		return 1;

	// Exact vertex budget: 12 corners + 30*(ld-1) edge points + 20 face interiors.
	nv = 12 + 30 * (ld - 1) + 20 * (ld - 1) * (ld - 2) / 2;
	nt = 20 * ld * ld;
	if (nv != 10 * ld * ld + 2)
		return 1;

	m->v = (float *)malloc(sizeof(float) * 3 * (size_t)nv);
	m->t = (int *)malloc(sizeof(int) * 3 * (size_t)nt);
	m->nrm = (float *)malloc(sizeof(float) * 3 * (size_t)nv);
	if (!m->v || !m->t || !m->nrm) {
		ss_mesh_free(m);
		return 1;
	}
	m->nv = nv;
	m->nt = nt;

	int off_edge = 12, off_face = 12 + 30 * (ld - 1), per_face = (ld - 1) * (ld - 2) / 2;

	// Global index of lattice point (i,j) on face f. i,j >= 0, i+j <= ld.
	// (0,0)=corner a, (ld,0)=corner b, (0,ld)=corner c.
	#define SS_VIDX(f, i, j) ss_lattice_index(bf, emap, off_edge, off_face, per_face, ld, f, i, j)

	for (int f = 0; f < 20; f++) {
		for (int j = 0; j <= ld; j++)
			for (int i = 0; i + j <= ld; i++) {
				int idx = SS_VIDX(f, i, j);
				float wa = (float)(ld - i - j) / ld, wb = (float)i / ld, wc = (float)j / ld;
				const float *a = bv[bf[f][0]], *b = bv[bf[f][1]], *c = bv[bf[f][2]];
				float x = wa * a[0] + wb * b[0] + wc * c[0];
				float y = wa * a[1] + wb * b[1] + wc * c[1];
				float z = wa * a[2] + wb * b[2] + wc * c[2];
				float len = sqrtf(x * x + y * y + z * z);
				if (len <= 0.0f)
					len = 1.0f;
				// Written repeatedly for shared corner/edge points; the value is
				// identical each time because it comes from the same barycentric
				// weights, so there is no order dependence.
				m->v[3 * idx + 0] = cx + radius * x / len;
				m->v[3 * idx + 1] = cy + radius * y / len;
				m->v[3 * idx + 2] = cz + radius * z / len;
			}
		for (int j = 0; j < ld; j++)
			for (int i = 0; i + j < ld; i++) {
				m->t[3 * ti + 0] = SS_VIDX(f, i, j);
				m->t[3 * ti + 1] = SS_VIDX(f, i + 1, j);
				m->t[3 * ti + 2] = SS_VIDX(f, i, j + 1);
				ti++;
				if (i + j < ld - 1) {
					m->t[3 * ti + 0] = SS_VIDX(f, i + 1, j);
					m->t[3 * ti + 1] = SS_VIDX(f, i + 1, j + 1);
					m->t[3 * ti + 2] = SS_VIDX(f, i, j + 1);
					ti++;
				}
			}
	}
	#undef SS_VIDX
	if (ti != nt) {
		ss_mesh_free(m);
		return 1;
	}

	if (ss_build_adjacency(m)) {
		ss_mesh_free(m);
		return 1;
	}
	ss_mesh_normals(m);
	return 0;
}

void ss_mesh_normals(ss_mesh *m) {
	if (!m || !m->nrm)
		return;
	memset(m->nrm, 0, sizeof(float) * 3 * (size_t)m->nv);
	for (int t = 0; t < m->nt; t++) {
		int a = m->t[3 * t + 0], b = m->t[3 * t + 1], c = m->t[3 * t + 2];
		const float *A = m->v + 3 * a, *B = m->v + 3 * b, *C = m->v + 3 * c;
		float u[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
		float w[3] = {C[0] - A[0], C[1] - A[1], C[2] - A[2]};
		float n[3] = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
		              u[0] * w[1] - u[1] * w[0]};
		// UNIT face normals, then a plain mean over the incident faces -- AFNI's
		// SUMA_SurfNorm normalises each face normal BEFORE accumulating it, so every
		// face counts equally regardless of area. Accumulating raw cross products
		// (the usual area weighting) tilts the vertex normal wherever triangle areas
		// are uneven; measured on T1w node 1000 it moved the normal by 3.3e-4, which
		// is enough to move one nearest-neighbour ray sample across a voxel boundary.
		float d = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (d == 0.0f) {
			n[0] = n[1] = n[2] = 1.0f;   // AFNI's degenerate-triangle fallback
		} else {
			n[0] /= d; n[1] /= d; n[2] /= d;
		}
		for (int e = 0; e < 3; e++) {
			int vi = m->t[3 * t + e];
			m->nrm[3 * vi + 0] += n[0];
			m->nrm[3 * vi + 1] += n[1];
			m->nrm[3 * vi + 2] += n[2];
		}
	}
	for (int i = 0; i < m->nv; i++) {
		float *n = m->nrm + 3 * i;
		// Divide by the incident-face count before normalising, as AFNI's SUMA_SurfNorm
		// does. The count comes from the adjacency degree rather than a separate tally: on
		// a closed triangle mesh a vertex has exactly as many incident faces as neighbours,
		// and this is only ever called on one. That matters because the tally used to be a
		// calloc inside a void function -- when it failed this returned early, leaving
		// m->nrm as uninitialised malloc memory while ss_icosphere still reported success
		// (found in audit by forcing the failure). Scaling before normalising cannot change
		// the direction, but it DOES change the rounding: dropping it moved three of the
		// five validation masks, so it stays.
		int deg = m->nbr_off ? (m->nbr_off[i + 1] - m->nbr_off[i]) : 0;
		float len;
		if (deg > 0) {
			n[0] /= deg; n[1] /= deg; n[2] /= deg;
		}
		len = sqrtf(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		if (len > 0.0f) {
			n[0] /= len; n[1] /= len; n[2] /= len;
		}
	}
}

void ss_mesh_smooth(ss_mesh *m, float lambda, float *scratch) {
	if (!m || !scratch)
		return;
	// Double-buffered: every vertex reads the COMPLETE prior state, so the result
	// does not depend on visit order and stays identical under any parallelisation.
	for (int i = 0; i < m->nv; i++) {
		int lo = m->nbr_off[i], hi = m->nbr_off[i + 1], n = hi - lo;
		float cx = 0, cy = 0, cz = 0;
		if (n <= 0) {
			scratch[3 * i + 0] = m->v[3 * i + 0];
			scratch[3 * i + 1] = m->v[3 * i + 1];
			scratch[3 * i + 2] = m->v[3 * i + 2];
			continue;
		}
		for (int k = lo; k < hi; k++) {
			const float *p = m->v + 3 * m->nbr[k];
			cx += p[0]; cy += p[1]; cz += p[2];
		}
		cx /= n; cy /= n; cz /= n;
		scratch[3 * i + 0] = m->v[3 * i + 0] + lambda * (cx - m->v[3 * i + 0]);
		scratch[3 * i + 1] = m->v[3 * i + 1] + lambda * (cy - m->v[3 * i + 1]);
		scratch[3 * i + 2] = m->v[3 * i + 2] + lambda * (cz - m->v[3 * i + 2]);
	}
	memcpy(m->v, scratch, sizeof(float) * 3 * (size_t)m->nv);
}

// --- triangle/triangle intersection -----------------------------------------

static float ss_dot3(const float a[3], const float b[3]) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static void ss_cross3(const float a[3], const float b[3], float o[3]) {
	o[0] = a[1] * b[2] - a[2] * b[1];
	o[1] = a[2] * b[0] - a[0] * b[2];
	o[2] = a[0] * b[1] - a[1] * b[0];
}

// Moller's interval-overlap test, written from the published method (a standard
// result, not AFNI code). Coplanar pairs are reported as NOT intersecting: for a
// closed surface the meaningful failure is a genuine crossing, and treating exact
// coplanarity as a hit makes adjacent flat regions produce false positives.
static int ss_tri_tri_cross(const float *p0, const float *p1, const float *p2,
		const float *q0, const float *q1, const float *q2) {
	float e1[3], e2[3], n1[3], n2[3];
	float du[3], dv[3];
	float dq[3], dp[3];
	float t1[2], t2[2];
	const float EPS = 1e-9f;

	for (int i = 0; i < 3; i++) { e1[i] = p1[i] - p0[i]; e2[i] = p2[i] - p0[i]; }
	ss_cross3(e1, e2, n1);
	for (int i = 0; i < 3; i++) { dq[i] = q0[i] - p0[i]; }
	du[0] = ss_dot3(n1, dq);
	for (int i = 0; i < 3; i++) { dq[i] = q1[i] - p0[i]; }
	du[1] = ss_dot3(n1, dq);
	for (int i = 0; i < 3; i++) { dq[i] = q2[i] - p0[i]; }
	du[2] = ss_dot3(n1, dq);
	if ((du[0] > EPS && du[1] > EPS && du[2] > EPS) ||
			(du[0] < -EPS && du[1] < -EPS && du[2] < -EPS))
		return 0; // q entirely on one side of p's plane
	if (fabsf(du[0]) <= EPS && fabsf(du[1]) <= EPS && fabsf(du[2]) <= EPS)
		return 0; // coplanar

	for (int i = 0; i < 3; i++) { e1[i] = q1[i] - q0[i]; e2[i] = q2[i] - q0[i]; }
	ss_cross3(e1, e2, n2);
	for (int i = 0; i < 3; i++) { dp[i] = p0[i] - q0[i]; }
	dv[0] = ss_dot3(n2, dp);
	for (int i = 0; i < 3; i++) { dp[i] = p1[i] - q0[i]; }
	dv[1] = ss_dot3(n2, dp);
	for (int i = 0; i < 3; i++) { dp[i] = p2[i] - q0[i]; }
	dv[2] = ss_dot3(n2, dp);
	if ((dv[0] > EPS && dv[1] > EPS && dv[2] > EPS) ||
			(dv[0] < -EPS && dv[1] < -EPS && dv[2] < -EPS))
		return 0;

	// Direction of the line of intersection of the two planes; project onto its
	// dominant axis and compare the two intervals.
	{
		float D[3];
		int ax = 0;
		float mx;
		const float *P[3] = {p0, p1, p2}, *Q[3] = {q0, q1, q2};
		ss_cross3(n1, n2, D);
		mx = fabsf(D[0]); ax = 0;
		if (fabsf(D[1]) > mx) { mx = fabsf(D[1]); ax = 1; }
		if (fabsf(D[2]) > mx) { mx = fabsf(D[2]); ax = 2; }
		if (mx <= 0.0f)
			return 0;
		// interval of triangle p on the line
		for (int side = 0; side < 2; side++) {
			const float **T = side ? Q : P;
			const float *dd = side ? dv : du;
			// For side 0 we need p's distances to q's plane (dv), and vice versa.
			const float *dist = side ? du : dv;
			float pv[3];
			int lo = -1;
			float out[2];
			int no = 0;
			(void)dd;
			for (int i = 0; i < 3; i++)
				pv[i] = T[i][ax];
			for (int i = 0; i < 3; i++) {
				int j = (i + 1) % 3;
				float a = dist[i], b = dist[j];
				if ((a > 0 && b < 0) || (a < 0 && b > 0)) {
					float s = a / (a - b);
					if (no < 2)
						out[no++] = pv[i] + s * (pv[j] - pv[i]);
				} else if (fabsf(a) <= EPS) {
					if (no < 2)
						out[no++] = pv[i];
				}
			}
			(void)lo;
			if (no < 2)
				return 0;
			if (out[0] > out[1]) { float tmp = out[0]; out[0] = out[1]; out[1] = tmp; }
			if (side == 0) { t1[0] = out[0]; t1[1] = out[1]; }
			else { t2[0] = out[0]; t2[1] = out[1]; }
		}
		if (t1[1] < t2[0] - EPS || t2[1] < t1[0] - EPS)
			return 0;
	}
	return 1;
}

// AFNI's SUMA_isSelfIntersect, ported FAITHFULLY -- including the part that makes it
// miss most intersections, because it is the retry DECISION and not a quality metric.
//
// AFNI walks every unique edge, ray-casts it against all triangles it does not belong
// to, and then accepts the hit only if the intersection point p satisfies
//     p[a] > ep1[a] && p[a] < ep2[a]   for a = x, y AND z
// with ep1/ep2 the edge's endpoints in stored (lower-node-index-first) order. That test
// can only ever pass when ep2 exceeds ep1 on all three axes, i.e. for edges pointing
// into the (+,+,+) octant -- roughly one edge in eight. AFNI therefore reports "No
// intersections found" on surfaces that genuinely have a handful of folds, and its
// retry fires far less often than a correct test would.
//
// We keep ss_mesh_self_intersections (a correct, order-independent triangle-pair test)
// for REPORTING, and use this one only to decide whether AFNI would have retried.
// Using the correct test for that decision made T1w take 3 expansion attempts where
// AFNI takes 1, which is a behavioural divergence, not a better answer.
//
// StopAt is 1 in the 3dSkullStrip call, so this stops at the first hit.
static long long ss_afni_self_intersect(const ss_mesh *m) {
	const double EPS = 1e-6;
	float lo[3] = {1e30f, 1e30f, 1e30f}, hi[3] = {-1e30f, -1e30f, -1e30f};
	int res, ncell;
	int *cell_cnt = NULL, *cell_off = NULL, *cell_it = NULL, *items = NULL;
	long long hits = 0;

	if (!m || m->nt <= 0)
		return 0;
	for (int i = 0; i < m->nv; i++)
		for (int a = 0; a < 3; a++) {
			float x = m->v[3 * i + a];
			if (x < lo[a]) lo[a] = x;
			if (x > hi[a]) hi[a] = x;
		}
	for (int a = 0; a < 3; a++)
		if (hi[a] <= lo[a]) hi[a] = lo[a] + 1.0f;
	res = (int)(cbrtf((float)m->nt) + 1.0f);
	if (res < 1) res = 1;
	if (res > 64) res = 64;
	ncell = res * res * res;
	cell_cnt = (int *)calloc((size_t)ncell + 1, sizeof(int));
	cell_off = (int *)malloc(sizeof(int) * ((size_t)ncell + 1));
	cell_it = (int *)calloc((size_t)ncell, sizeof(int));
	if (!cell_cnt || !cell_off || !cell_it) {
		free(cell_cnt); free(cell_off); free(cell_it);
		return -1;
	}
	// Cell range of an axis-aligned box. BOTH bounds are clamped, not just c1: a box whose
	// min maps to exactly `res` (mn == hi[a], reachable when a vertex sits on the bounding-box
	// face) leaves c0 == res, and the trailing `c1 < c0` line then RAISES c1 back to that
	// out-of-range value, after which the cell loop writes past the end of items[]/cell_it[].
	// Verified before the fix: a synthetic 220^3 phantom crashed with SIGBUS inside
	// ss_mesh_self_intersections, reproducibly, reached through skullstrip_run.
	// (No // comments inside the macro body -- they would swallow the line continuations.)
	#define SS_BOXCELL(mn, mx, c0, c1)                                                \
		do {                                                                          \
			for (int a = 0; a < 3; a++) {                                             \
				float f0 = (mn[a] - lo[a]) / (hi[a] - lo[a]) * res;                   \
				float f1 = (mx[a] - lo[a]) / (hi[a] - lo[a]) * res;                   \
				c0[a] = (int)f0; c1[a] = (int)f1;                                     \
				if (c0[a] < 0) c0[a] = 0;                                             \
				if (c0[a] > res - 1) c0[a] = res - 1;                                 \
				if (c1[a] > res - 1) c1[a] = res - 1;                                 \
				if (c1[a] < c0[a]) c1[a] = c0[a];                                     \
			}                                                                         \
		} while (0)

	for (int pass = 0; pass < 2; pass++) {
		for (int t = 0; t < m->nt; t++) {
			float mn[3], mx[3];
			int c0[3], c1[3];
			for (int a = 0; a < 3; a++) {
				mn[a] = 1e30f; mx[a] = -1e30f;
				for (int e = 0; e < 3; e++) {
					float x = m->v[3 * m->t[3 * t + e] + a];
					if (x < mn[a]) mn[a] = x;
					if (x > mx[a]) mx[a] = x;
				}
			}
			SS_BOXCELL(mn, mx, c0, c1);
			for (int z = c0[2]; z <= c1[2]; z++)
				for (int y = c0[1]; y <= c1[1]; y++)
					for (int x = c0[0]; x <= c1[0]; x++) {
						int c = x + y * res + z * res * res;
						if (!pass) cell_cnt[c]++;
						else items[cell_off[c] + cell_it[c]++] = t;
					}
		}
		if (!pass) {
			int acc = 0;
			for (int c = 0; c < ncell; c++) { cell_off[c] = acc; acc += cell_cnt[c]; }
			cell_off[ncell] = acc;
			items = (int *)malloc(sizeof(int) * (size_t)(acc > 0 ? acc : 1));
			if (!items) {
				free(cell_cnt); free(cell_off); free(cell_it);
				return -1;
			}
		}
	}

	for (int i = 0; i < m->nv && !hits; i++)
		for (int q = m->nbr_off[i]; q < m->nbr_off[i + 1] && !hits; q++) {
			int j = m->nbr[q];
			const float *ep1, *ep2;
			float mn[3], mx[3];
			int c0[3], c1[3];
			if (j < i)
				continue;   // each undirected edge once, lower index first (SUMA's order)
			ep1 = m->v + 3 * i;
			ep2 = m->v + 3 * j;
			for (int a = 0; a < 3; a++) {
				mn[a] = ep1[a] < ep2[a] ? ep1[a] : ep2[a];
				mx[a] = ep1[a] > ep2[a] ? ep1[a] : ep2[a];
			}
			SS_BOXCELL(mn, mx, c0, c1);
			for (int z = c0[2]; z <= c1[2] && !hits; z++)
				for (int y = c0[1]; y <= c1[1] && !hits; y++)
					for (int x = c0[0]; x <= c1[0] && !hits; x++) {
						int c = x + y * res + z * res * res;
						for (int r = cell_off[c]; r < cell_off[c + 1]; r++) {
							int t = items[r];
							int a0 = m->t[3 * t], a1 = m->t[3 * t + 1], a2 = m->t[3 * t + 2];
							const float *P0, *P1, *P2;
							double e1v[3], e2v[3], pv[3], qv[3], tv[3], det, invd, u, v, tt;
							double px, py, pz;
							if (a0 == i || a1 == i || a2 == i || a0 == j || a1 == j || a2 == j)
								continue;   // AFNI excludes the edge's own two triangles
							P0 = m->v + 3 * a0; P1 = m->v + 3 * a1; P2 = m->v + 3 * a2;
							for (int a = 0; a < 3; a++) {
								e1v[a] = P1[a] - P0[a];
								e2v[a] = P2[a] - P0[a];
								qv[a] = ep2[a] - ep1[a];
							}
							pv[0] = qv[1] * e2v[2] - qv[2] * e2v[1];
							pv[1] = qv[2] * e2v[0] - qv[0] * e2v[2];
							pv[2] = qv[0] * e2v[1] - qv[1] * e2v[0];
							det = e1v[0] * pv[0] + e1v[1] * pv[1] + e1v[2] * pv[2];
							if (det > -1e-12 && det < 1e-12)
								continue;
							invd = 1.0 / det;
							for (int a = 0; a < 3; a++) tv[a] = ep1[a] - P0[a];
							u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * invd;
							if (u <= EPS || u >= 1.0)
								continue;
							{
								double qq[3];
								qq[0] = tv[1] * e1v[2] - tv[2] * e1v[1];
								qq[1] = tv[2] * e1v[0] - tv[0] * e1v[2];
								qq[2] = tv[0] * e1v[1] - tv[1] * e1v[0];
								v = (qv[0] * qq[0] + qv[1] * qq[1] + qv[2] * qq[2]) * invd;
								if (v <= EPS || u + v >= 1.0)
									continue;
								tt = (e2v[0] * qq[0] + e2v[1] * qq[1] + e2v[2] * qq[2]) * invd;
							}
							(void)tt;
							px = P0[0] + u * e1v[0] + v * e2v[0];
							py = P0[1] + u * e1v[1] + v * e2v[1];
							pz = P0[2] + u * e1v[2] + v * e2v[2];
							// AFNI's box test, verbatim: strictly between ep1 and ep2 in the
							// STORED endpoint order, on all three axes.
							if (px > ep1[0] && px < ep2[0] && py > ep1[1] && py < ep2[1] &&
									pz > ep1[2] && pz < ep2[2]) {
								hits++;
								break;
							}
						}
					}
		}
	#undef SS_BOXCELL
	free(cell_cnt); free(cell_off); free(cell_it); free(items);
	return hits;
}

long long ss_mesh_self_intersections(const ss_mesh *m) {
	// Uniform spatial grid over triangle bounding boxes: each triangle is compared
	// only against triangles sharing a cell, so this is not the O(n^2) all-pairs
	// test the plan warns about.
	float lo[3] = {1e30f, 1e30f, 1e30f}, hi[3] = {-1e30f, -1e30f, -1e30f};
	int res, ncell;
	int *cell_cnt = NULL, *cell_off = NULL, *cell_it = NULL, *items = NULL;
	long long hits = 0;

	if (!m || m->nt <= 0)
		return 0;
	for (int i = 0; i < m->nv; i++)
		for (int a = 0; a < 3; a++) {
			float x = m->v[3 * i + a];
			if (x < lo[a]) lo[a] = x;
			if (x > hi[a]) hi[a] = x;
		}
	for (int a = 0; a < 3; a++)
		if (hi[a] <= lo[a])
			hi[a] = lo[a] + 1.0f;
	res = (int)(cbrtf((float)m->nt) + 1.0f);
	if (res < 1) res = 1;
	if (res > 64) res = 64;
	ncell = res * res * res;

	cell_cnt = (int *)calloc((size_t)ncell + 1, sizeof(int));
	cell_off = (int *)malloc(sizeof(int) * ((size_t)ncell + 1));
	cell_it = (int *)calloc((size_t)ncell, sizeof(int));
	if (!cell_cnt || !cell_off || !cell_it) {
		free(cell_cnt); free(cell_off); free(cell_it);
		return -1;
	}

	// Same box->cell mapping as SS_BOXCELL above, including the c0 upper clamp; see the note
	// there for the out-of-range case it prevents.
	#define SS_CELL_RANGE(t, c0, c1)                                                  \
		do {                                                                          \
			for (int a = 0; a < 3; a++) {                                             \
				float mn = 1e30f, mx = -1e30f;                                        \
				for (int e = 0; e < 3; e++) {                                         \
					float x = m->v[3 * m->t[3 * (t) + e] + a];                        \
					if (x < mn) mn = x;                                               \
					if (x > mx) mx = x;                                               \
				}                                                                     \
				c0[a] = (int)((mn - lo[a]) / (hi[a] - lo[a]) * res);                  \
				c1[a] = (int)((mx - lo[a]) / (hi[a] - lo[a]) * res);                  \
				if (c0[a] < 0) c0[a] = 0;                                             \
				if (c0[a] >= res) c0[a] = res - 1;                                    \
				if (c1[a] >= res) c1[a] = res - 1;                                    \
				if (c1[a] < c0[a]) c1[a] = c0[a];                                     \
			}                                                                         \
		} while (0)

	for (int t = 0; t < m->nt; t++) {
		int c0[3], c1[3];
		SS_CELL_RANGE(t, c0, c1);
		for (int k = c0[2]; k <= c1[2]; k++)
			for (int j = c0[1]; j <= c1[1]; j++)
				for (int i = c0[0]; i <= c1[0]; i++)
					cell_cnt[(k * res + j) * res + i]++;
	}
	cell_off[0] = 0;
	for (int c = 0; c < ncell; c++)
		cell_off[c + 1] = cell_off[c] + cell_cnt[c];
	items = (int *)malloc(sizeof(int) * (size_t)(cell_off[ncell] > 0 ? cell_off[ncell] : 1));
	if (!items) {
		free(cell_cnt); free(cell_off); free(cell_it);
		return -1;
	}
	for (int t = 0; t < m->nt; t++) {
		int c0[3], c1[3];
		SS_CELL_RANGE(t, c0, c1);
		for (int k = c0[2]; k <= c1[2]; k++)
			for (int j = c0[1]; j <= c1[1]; j++)
				for (int i = c0[0]; i <= c1[0]; i++) {
					int c = (k * res + j) * res + i;
					items[cell_off[c] + cell_it[c]++] = t;
				}
	}
	#undef SS_CELL_RANGE

	for (int c = 0; c < ncell; c++) {
		int b = cell_off[c], e = cell_off[c + 1];
		for (int x = b; x < e; x++)
			for (int y = x + 1; y < e; y++) {
				int ta = items[x], tb = items[y];
				int share = 0;
				if (ta >= tb)
					continue;
				for (int u = 0; u < 3 && !share; u++)
					for (int w = 0; w < 3; w++)
						if (m->t[3 * ta + u] == m->t[3 * tb + w]) { share = 1; break; }
				if (share)
					continue; // adjacent triangles are not self-intersections
				if (ss_tri_tri_cross(m->v + 3 * m->t[3 * ta + 0], m->v + 3 * m->t[3 * ta + 1],
						m->v + 3 * m->t[3 * ta + 2], m->v + 3 * m->t[3 * tb + 0],
						m->v + 3 * m->t[3 * tb + 1], m->v + 3 * m->t[3 * tb + 2])) {
					// A pair can meet in several cells; count it exactly once, in the
					// lowest-indexed cell where BOTH are registered -- i.e. the min
					// corner of the INTERSECTION of their two cell ranges. Using the
					// combined bounding box instead is wrong: its min corner may be a
					// cell neither triangle occupies, in which case the pair is never
					// counted anywhere and a genuine fold reads as zero intersections.
					int lo_a[3], hi_a[3], lo_b[3], hi_b[3], o0[3];
					int drop = 0, first;
					for (int a = 0; a < 3; a++) {
						float amn = 1e30f, amx = -1e30f, bmn = 1e30f, bmx = -1e30f;
						for (int q = 0; q < 3; q++) {
							float xa = m->v[3 * m->t[3 * ta + q] + a];
							float xb = m->v[3 * m->t[3 * tb + q] + a];
							if (xa < amn) amn = xa;
							if (xa > amx) amx = xa;
							if (xb < bmn) bmn = xb;
							if (xb > bmx) bmx = xb;
						}
						lo_a[a] = (int)((amn - lo[a]) / (hi[a] - lo[a]) * res);
						hi_a[a] = (int)((amx - lo[a]) / (hi[a] - lo[a]) * res);
						lo_b[a] = (int)((bmn - lo[a]) / (hi[a] - lo[a]) * res);
						hi_b[a] = (int)((bmx - lo[a]) / (hi[a] - lo[a]) * res);
						if (lo_a[a] < 0) lo_a[a] = 0;
						if (lo_b[a] < 0) lo_b[a] = 0;
						if (hi_a[a] >= res) hi_a[a] = res - 1;
						if (hi_b[a] >= res) hi_b[a] = res - 1;
						if (hi_a[a] < lo_a[a]) hi_a[a] = lo_a[a];
						if (hi_b[a] < lo_b[a]) hi_b[a] = lo_b[a];
						o0[a] = lo_a[a] > lo_b[a] ? lo_a[a] : lo_b[a];
						{
							int o1 = hi_a[a] < hi_b[a] ? hi_a[a] : hi_b[a];
							if (o1 < o0[a])
								drop = 1;
						}
					}
					first = drop ? c : (o0[2] * res + o0[1]) * res + o0[0];
					if (first == c)
						hits++;
				}
			}
	}
	 free(cell_cnt); free(cell_off); free(cell_it); free(items);
	return hits;
}

// --- closed-surface rasterisation -------------------------------------------
//
// Parity fill along x. For each voxel row (j,k) we find every triangle whose
// projection onto the (y,z) plane contains the row's (y,z) point, take the x of the
// intersection, sort, and fill between alternate crossings.
//
// The correctness hinge is that a row passing exactly along a shared triangle edge
// must be counted ONCE, not twice and not zero times, or the surface leaks. That is
// handled by a half-open rule: a boundary hit counts only for the triangle where the
// corresponding edge is a "top-left" edge in (y,z). Adjacent triangles traverse a
// shared edge in opposite directions, so exactly one of them owns it. No epsilon and
// no ray-direction dependence.

// Which of the two triangles sharing an edge owns a sample lying exactly on it.
// The orientation of this predicate is not free: it must select the LOWER-coordinate
// side, so that (y,z) inclusion matches the lower-inclusive [entry,exit) fill used
// along x. With it inverted the cube's count is right but the region is shifted one
// voxel, which reads as a rasteriser bug and is really an axis-convention mismatch.
static int ss_topleft(float ey, float ez) {
	return (ez < 0.0f) || (ez == 0.0f && ey > 0.0f);
}


static void ss_sort_floats(float *a, int n) {
	for (int i = 1; i < n; i++) {
		float key = a[i];
		int j = i - 1;
		while (j >= 0 && a[j] > key) { a[j + 1] = a[j]; j--; }
		a[j + 1] = key;
	}
}

int ss_mesh_rasterize(const ss_mesh *m, int nx, int ny, int nz, unsigned char *mask) {
	long long nxy = (long long)nx * ny;
	int *cnt = NULL, *off = NULL, *it = NULL, *items = NULL;
	float *xs = NULL;
	int maxrow = 0;

	if (!m || !mask || nx <= 0 || ny <= 0 || nz <= 0)
		return 1;
	memset(mask, 0, (size_t)nxy * nz);

	// Bucket triangles by the (j,k) rows they can touch, so each row only tests
	// triangles that actually span it.
	cnt = (int *)calloc((size_t)ny * nz + 1, sizeof(int));
	off = (int *)malloc(sizeof(int) * ((size_t)ny * nz + 1));
	it = (int *)calloc((size_t)ny * nz, sizeof(int));
	if (!cnt || !off || !it) {
		free(cnt); free(off); free(it);
		return 1;
	}
	#define SS_ROWS(t, j0, j1, k0, k1)                                         \
		do {                                                                   \
			float ymn = 1e30f, ymx = -1e30f, zmn = 1e30f, zmx = -1e30f;        \
			for (int e = 0; e < 3; e++) {                                      \
				const float *P = m->v + 3 * m->t[3 * (t) + e];                 \
				if (P[1] < ymn) ymn = P[1];                                    \
				if (P[1] > ymx) ymx = P[1];                                    \
				if (P[2] < zmn) zmn = P[2];                                    \
				if (P[2] > zmx) zmx = P[2];                                    \
			}                                                                  \
			j0 = (int)ceilf(ymn); j1 = (int)floorf(ymx);                       \
			k0 = (int)ceilf(zmn); k1 = (int)floorf(zmx);                       \
			if (j0 < 0) j0 = 0;                                                \
			if (k0 < 0) k0 = 0;                                                \
			if (j1 > ny - 1) j1 = ny - 1;                                      \
			if (k1 > nz - 1) k1 = nz - 1;                                      \
		} while (0)

	for (int t = 0; t < m->nt; t++) {
		int j0, j1, k0, k1;
		SS_ROWS(t, j0, j1, k0, k1);
		for (int k = k0; k <= k1; k++)
			for (int j = j0; j <= j1; j++)
				cnt[(long long)k * ny + j]++;
	}
	off[0] = 0;
	for (long long r = 0; r < (long long)ny * nz; r++)
		off[r + 1] = off[r] + cnt[r];
	items = (int *)malloc(sizeof(int) * (size_t)(off[(long long)ny * nz] > 0 ? off[(long long)ny * nz] : 1));
	if (!items) {
		free(cnt); free(off); free(it);
		return 1;
	}
	for (int t = 0; t < m->nt; t++) {
		int j0, j1, k0, k1;
		SS_ROWS(t, j0, j1, k0, k1);
		for (int k = k0; k <= k1; k++)
			for (int j = j0; j <= j1; j++) {
				long long r = (long long)k * ny + j;
				items[off[r] + it[r]++] = t;
			}
	}
	#undef SS_ROWS
	for (long long r = 0; r < (long long)ny * nz; r++)
		if (off[r + 1] - off[r] > maxrow)
			maxrow = off[r + 1] - off[r];
	xs = (float *)malloc(sizeof(float) * (size_t)(maxrow > 0 ? maxrow : 1));
	if (!xs) {
		free(cnt); free(off); free(it); free(items);
		return 1;
	}

	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++) {
			long long r = (long long)k * ny + j;
			int b = off[r], e = off[r + 1], n = 0;
			float py = (float)j, pz = (float)k;
			for (int q = b; q < e; q++) {
				int t = items[q];
				const float *A = m->v + 3 * m->t[3 * t + 0];
				const float *B = m->v + 3 * m->t[3 * t + 1];
				const float *C = m->v + 3 * m->t[3 * t + 2];
				// Edge functions in the (y,z) plane, in DOUBLE. Not gold-plating: in
				// float, a face lying exactly on a plane of voxel centres interpolates
				// its crossing to 20.0000019 instead of 20, and the half-open ceil()
				// then leaks one voxel past the surface. In double the products and
				// their sum are exact for representable coordinates, so the crossing
				// lands exactly on the face and the fill is exact.
				double e0 = ((double)B[1] - A[1]) * ((double)pz - A[2]) - ((double)B[2] - A[2]) * ((double)py - A[1]);
				double e1 = ((double)C[1] - B[1]) * ((double)pz - B[2]) - ((double)C[2] - B[2]) * ((double)py - B[1]);
				double e2 = ((double)A[1] - C[1]) * ((double)pz - C[2]) - ((double)A[2] - C[2]) * ((double)py - C[1]);
				double area = ((double)B[1] - A[1]) * ((double)C[2] - A[2]) - ((double)B[2] - A[2]) * ((double)C[1] - A[1]);
				double wsum, xh;
				int inside;
				if (area == 0.0)
					continue; // degenerate in projection: contributes no crossing
				// A back-facing triangle projects with negative area. Normalising the
				// edge functions by negating them ALSO reverses the effective traversal
				// direction, so the tie-break must see reversed edge vectors too.
				// Missing that makes the two oppositely-wound faces of a closed solid
				// disagree about which silhouette edge they own, and rows whose centres
				// lie exactly on a face fill inconsistently.
				{
					float sgn = (area < 0.0) ? -1.0f : 1.0f;
					if (area < 0.0) { e0 = -e0; e1 = -e1; e2 = -e2; }
					inside = 1;
					if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0)
						inside = 0;
					else {
						// Half-open: a hit exactly ON an edge belongs to the triangle for
						// which that edge is top-left. Its neighbour sees the same edge
						// reversed and therefore declines it.
						if (e0 == 0.0 && !ss_topleft(sgn * (B[1] - A[1]), sgn * (B[2] - A[2]))) inside = 0;
						if (e1 == 0.0 && !ss_topleft(sgn * (C[1] - B[1]), sgn * (C[2] - B[2]))) inside = 0;
						if (e2 == 0.0 && !ss_topleft(sgn * (A[1] - C[1]), sgn * (A[2] - C[2]))) inside = 0;
					}
				}
				if (!inside)
					continue;
				// Normalise by the ACTUAL weight sum rather than |area| so the three
				// weights sum to exactly 1; this is what makes a crossing on a planar
				// face come back exactly on the face.
				wsum = e0 + e1 + e2;
				if (wsum == 0.0)
					continue;
				xh = (e1 * (double)A[0] + e2 * (double)B[0] + e0 * (double)C[0]) / wsum;
				if (n < maxrow)
					xs[n++] = (float)xh;
			}
			if (n < 2)
				continue;
			ss_sort_floats(xs, n);
			for (int q = 0; q + 1 < n; q += 2) {
				// HALF-OPEN in x, to match the half-open edge rule used in (y,z).
				// A voxel is inside iff its centre lies in [entry, exit). Mixing a
				// closed x-fill with a half-open row test makes a surface whose faces
				// land exactly on voxel centres resolve differently per axis, which
				// showed up as a cube rasterising to neither 12^3 nor 13^3.
				int i0 = (int)ceilf(xs[q]), i1 = (int)ceilf(xs[q + 1]) - 1;
				if (i0 < 0) i0 = 0;
				if (i1 > nx - 1) i1 = nx - 1;
				for (int i = i0; i <= i1; i++)
					mask[i + (long long)j * nx + (long long)k * nxy] = 1;
			}
		}

	free(cnt); free(off); free(it); free(items); free(xs);
	return 0;
}

// ===========================================================================
// Deformation core. ADAPTED from AFNI's SUMA_BrainWrap.c (SUMA_StretchToFitLeCerveau,
// SUMA_LoadPrepInVol, SUMA_Find_IminImax), a non-copyrightable US Government work; the
// option defaults come from the public block in SUMA_3dSkullStrip.c. The underlying method
// is Smith 2002, "Fast robust automated brain extraction" (HBM 17:143-155), which AFNI's
// own -help names as the algorithm it modifies:
//   tb = (Imax - t2) * SF + t2   with SF (shrink_fac) default 0.6
//
// SUMA_3dedge3 is the one carve-out and is never called -- see the top of this file.
// ===========================================================================

int ss_intensity_stats(const unsigned char *vol, int nx, int ny, int nz, ss_stats *st) {
	long long nxy = (long long)nx * ny, nxyz = nxy * nz;
	long long hist[256], cum, target;
	int i;
	double sx = 0, sy = 0, sz = 0, sw = 0;

	if (!vol || !st)
		return 1;
	memset(st, 0, sizeof(*st));
	memset(hist, 0, sizeof(hist));
	for (long long q = 0; q < nxyz; q++)
		hist[vol[q]]++;

	// t2 / t98: the 2nd and 98th percentiles. BET takes them over the whole image.
	cum = 0; target = (long long)(0.02 * nxyz);
	for (i = 0; i < 256 && cum <= target; i++)
		cum += hist[i];
	st->t2 = (float)(i > 0 ? i - 1 : 0);
	cum = 0; target = (long long)(0.98 * nxyz);
	for (i = 0; i < 256 && cum <= target; i++)
		cum += hist[i];
	st->t98 = (float)(i > 0 ? i - 1 : 0);
	if (st->t98 <= st->t2) {
		// Degenerate contrast: a uniform or near-uniform volume. Fail clearly rather
		// than expanding a surface into a meaningless answer -- this is the gap the
		// Milestone 1 robustness pass recorded and deferred to here.
		printfx("skullstrip: degenerate contrast (2nd and 98th percentiles both %g)\n",
				(double)st->t2);
		return 1;
	}
	st->t = st->t2 + 0.1f * (st->t98 - st->t2);

	// Centre of gravity over voxels above t, with intensity capped at t98 so a few
	// very bright voxels (fat, vessels) cannot drag the centre.
	for (long long k = 0, q = 0; k < nz; k++)
		for (long long j = 0; j < ny; j++)
			for (long long ii = 0; ii < nx; ii++, q++) {
				float v = (float)vol[q];
				if (v <= st->t)
					continue;
				if (v > st->t98)
					v = st->t98;
				sw += v; sx += v * ii; sy += v * j; sz += v * k;
				st->nabove++;
			}
	if (sw <= 0.0 || st->nabove < 1000) {
		printfx("skullstrip: too few voxels above threshold to locate the head\n");
		return 1;
	}
	st->cog[0] = (float)(sx / sw);
	st->cog[1] = (float)(sy / sw);
	st->cog[2] = (float)(sz / sw);
	// Equivalent-sphere radius of the supra-threshold volume (grid is 1 mm isotropic).
	st->radius = (float)cbrt(3.0 * (double)st->nabove / (4.0 * M_PI));

	// tm: median intensity inside a sphere of that radius about the COG, counting
	// only voxels above t (BET's definition).
	{
		long long h2[256], n2 = 0;
		memset(h2, 0, sizeof(h2));
		int i0 = (int)(st->cog[0] - st->radius), i1 = (int)(st->cog[0] + st->radius);
		int j0 = (int)(st->cog[1] - st->radius), j1 = (int)(st->cog[1] + st->radius);
		int k0 = (int)(st->cog[2] - st->radius), k1 = (int)(st->cog[2] + st->radius);
		if (i0 < 0) i0 = 0; if (j0 < 0) j0 = 0; if (k0 < 0) k0 = 0;
		if (i1 > nx - 1) i1 = nx - 1;
		if (j1 > ny - 1) j1 = ny - 1;
		if (k1 > nz - 1) k1 = nz - 1;
		for (int k = k0; k <= k1; k++)
			for (int j = j0; j <= j1; j++)
				for (int ii = i0; ii <= i1; ii++) {
					float dx = ii - st->cog[0], dy = j - st->cog[1], dz = k - st->cog[2];
					if (dx * dx + dy * dy + dz * dz > st->radius * st->radius)
						continue;
					{
						// EVERY voxel inside the sphere, dark ones included. AFNI takes
						// the plain median of the in-sphere values; BET's published
						// definition restricts to voxels above t, and using that here
						// drags tm up (measured: 106 vs AFNI's 94 on T1w). tm caps Imax,
						// which sets tb = (Imax-t2)*lZt + t2, so a tm that is ~13% high
						// inflates every node's threshold in the expansion direction.
						unsigned char v = vol[ii + (long long)j * nx + (long long)k * nxy];
						h2[v]++;
						n2++;
					}
				}
		if (n2 < 1) {
			printfx("skullstrip: no brain-like signal near the centre of gravity\n");
			return 1;
		}
		cum = 0;
		for (i = 0; i < 256; i++) {
			cum += h2[i];
			if (cum >= n2 / 2)
				break;
		}
		st->tm = (float)i;
	}
	// The deformation's image force is f3 = 2*(Imin - tb)/(Imax - t2), and Imax is capped at
	// tm, so its denominator collapses to EXACTLY tm - t2 for every node. A volume with
	// tm == t2 (measured on a hollow-shell phantom: t2=0, t98=255, tm=0 -- legal scalar 3D
	// float32) therefore makes su3 NaN from iteration 0, every vertex goes NaN, and the
	// coordinate clamp lets NaN through because its comparisons are false. The run then
	// spends a full multi-stage expansion to report "produced an empty mask". Checking t98
	// alone is not enough.
	if (st->tm <= st->t2) {
		printfx("skullstrip: degenerate contrast (median inside the head sphere equals the "
				"2nd percentile, %g)\n", (double)st->t2);
		return 1;
	}
	SSV("skullstrip: t2=%g t98=%g t=%g tm=%g cog=[%g %g %g] r=%g nabove=%lld\n",
			(double)st->t2, (double)st->t98, (double)st->t, (double)st->tm,
			(double)st->cog[0], (double)st->cog[1], (double)st->cog[2],
			(double)st->radius, st->nabove);
	return 0;
}

// Continuous coordinate -> integer voxel index on a regular grid, plus an out-of-range
// report. Sampling through this is NEAREST NEIGHBOUR and out-of-grid coordinates CLAMP --
// they do not read as 0; see the conventions at ss_probe_t, where NN is load-bearing
// convention 1. CLEAN-ROOM: derived from a table of 19,139 input/output observations by an
// implementer who had read neither AFNI nor niimath, working only from that table. The kit,
// the table, the derivation account and the protocol are in skullstrip_bench/clean_room/.
//
// It replaced an earlier version written after reading AFNI's THD_3dmm_to_3dind_warn in
// src/thd_coords.c -- a Medical College of Wisconsin file from 1998, outside AFNI's
// post-2001 public-domain rule. Everything else this file adapts is public-domain AFNI that
// was vetted first; that one was reached indirectly, by following a call out of the vetted
// SUMA_BrainWrap.c. Do not "simplify" this back toward the shape of that function.
//
// (thd_coords.c was GPL-2 when this was clean-roomed; MCW relicensed its 1994-2000 AFNI code
// to CC BY 4.0 on 2026-05-12, so the copyleft bar is gone -- but CC BY's attribution and
// change-notice duties are not, which is why the clean-room result remains preferable.)
//
// Three properties the table forced, and each is load-bearing:
//   - the offset is 0.49, not 0.5, in VOXEL units, so an exact half-voxel coordinate falls
//     to the LOWER index and every boundary shifts by 1% of a voxel;
//   - the arithmetic is SINGLE precision -- in double it disagrees on the anisotropic grids;
//   - the conversion TRUNCATES toward zero rather than flooring. Invisible in the index
//     (negatives clamp to 0 either way) but visible in the out-of-range flag, which is
//     asymmetric: a coordinate a whole voxel below the grid still reports in range.
static int ss_world_to_index(float coord, float origin, float delta, int n, int *out_of_range) {
	float v = (coord - origin) / delta;
	int i = (int)(v + 0.49f);   // truncates toward zero; see the note above

	*out_of_range = 0;
	if (i < 0) {
		i = 0;
		*out_of_range = 1;
	} else if (i > n - 1) {
		i = n - 1;
		*out_of_range = 1;
	}
	return i;
}

// The out-of-grid tail of the ray sampler, in its own noinline function ON PURPOSE. Left
// inline, clang if-converts it: it computes the clamped index AND the flag on every step with
// ~20 extra cset/csel/ccmp, so a path taken at most once per ray is paid for by all 37 samples.
// Out of line there is nothing to if-convert and the hot loop keeps only three unsigned
// compares. Used by the SS_FAST kernel; faithful mode keeps the branch inline, because moving
// it perturbs the surrounding float contraction (see skullstrip_kernel.h).
// MSVC has no __attribute__; every GNU extension in this tree carries a _MSC_VER shim (see the
// #pragma pack pair in niimath.c/meshify.c). This was the only unguarded one in the skullstrip
// files, and CMake happily configures ENABLE_SKULLSTRIP=ON for 64-bit MSVC, so it would have
// configured clean and then hard-failed inside cl.exe -- the exact failure mode the pointer-size
// and flavour guards above it were added to prevent. The skullstrip-msvc job in
// skullstrip-build.yml is the only CI that compiles this file with cl.exe -- and it has never
// been observed green, so treat MSVC as unverified.
#ifdef _MSC_VER
#define SS_NOINLINE __declspec(noinline)
#else
#define SS_NOINLINE __attribute__((noinline))
#endif

static SS_NOINLINE void ss_sample_clamp(int nx, int ny, int nz,
		float x, float y, float z, int *i, int *j, int *k, int *out) {
	int oi, oj, ok;
	*i = ss_world_to_index(x, 0.0f, 1.0f, nx, &oi);
	*j = ss_world_to_index(y, 0.0f, 1.0f, ny, &oj);
	*k = ss_world_to_index(z, 0.0f, 1.0f, nz, &ok);
	*out = (oi || oj || ok) ? 1 : 0;   // one flag for the triple, as the reference reports it
}

// Mean of the UNIQUE edge lengths, AFNI's SUMA_MEAN_SEGMENT_LENGTH. Averaging over
// directed neighbour pairs gives the same value (every edge is counted twice), but
// this states the definition the constant it feeds actually has.
static double ss_mean_seg_len(const ss_mesh *m) {
	double acc = 0.0;
	long long cnt = 0;
	for (int i = 0; i < m->nv; i++)
		for (int q = m->nbr_off[i]; q < m->nbr_off[i + 1]; q++) {
			int j = m->nbr[q];
			double dx = m->v[3 * i] - m->v[3 * j], dy = m->v[3 * i + 1] - m->v[3 * j + 1],
			       dz = m->v[3 * i + 2] - m->v[3 * j + 2];
			acc += sqrt(dx * dx + dy * dy + dz * dz);
			cnt++;
		}
	return cnt ? acc / (double)cnt : 1.0;
}

// Mean distance from centre over all vertices (AFNI's SUMA_SO_RADIUS).
static double ss_mesh_radius(const ss_mesh *m, const float c[3]) {
	double acc = 0.0;
	for (int i = 0; i < m->nv; i++) {
		double dx = m->v[3 * i] - c[0], dy = m->v[3 * i + 1] - c[1],
		       dz = m->v[3 * i + 2] - c[2];
		acc += sqrt(dx * dx + dy * dy + dz * dz);
	}
	return m->nv ? acc / (double)m->nv : 0.0;
}

// Total triangle area, for the Stage-1 "is the surface still growing" test.
static double ss_mesh_area(const ss_mesh *m) {
	double acc = 0.0;
	for (int f = 0; f < m->nt; f++) {
		const float *a = m->v + 3 * m->t[3 * f], *b = m->v + 3 * m->t[3 * f + 1],
		            *c = m->v + 3 * m->t[3 * f + 2];
		double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
		double w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
		double cx = u[1] * w[2] - u[2] * w[1], cy = u[2] * w[0] - u[0] * w[2],
		       cz = u[0] * w[1] - u[1] * w[0];
		acc += 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
	}
	return acc;
}

// AFNI's SUMA_WRAP_BRAIN_SMOOTH_NN: `nsmooth` sequential nearest-neighbour geometric
// smoothing passes -- each pass is the UNWEIGHTED mean of the node and its neighbours,
// i.e. lambda = deg/(deg+1) ~ 0.857, NOT a tunable lambda -- followed by a RADIUS
// RESTORATION that scales every node back out by the mean-radius the smoothing lost.
// Dropping the restoration turns a smoothing pass into a shrink, which is why a
// lambda-0.5 pass every 10 iterations is not a substitute for this.
// scratch must hold 6*nv floats.
static void ss_smooth_nn_wrap(ss_mesh *m, int nsmooth, const float c[3], float *scratch) {
	float *a = scratch, *b = scratch + 3 * (size_t)m->nv;
	float *cur = m->v;
	double rref = ss_mesh_radius(m, c), r, dr;

	for (int pass = 0; pass < nsmooth; pass++) {
		float *dst = (pass & 1) ? a : b;
		for (int i = 0; i < m->nv; i++) {
			float sx = cur[3 * i], sy = cur[3 * i + 1], sz = cur[3 * i + 2];
			int deg = m->nbr_off[i + 1] - m->nbr_off[i];
			for (int q = m->nbr_off[i]; q < m->nbr_off[i + 1]; q++) {
				sx += cur[3 * m->nbr[q]];
				sy += cur[3 * m->nbr[q] + 1];
				sz += cur[3 * m->nbr[q] + 2];
			}
			dst[3 * i] = sx / (deg + 1.0f);
			dst[3 * i + 1] = sy / (deg + 1.0f);
			dst[3 * i + 2] = sz / (deg + 1.0f);
		}
		cur = dst;
	}
	if (cur != m->v)
		memcpy(m->v, cur, sizeof(float) * 3 * (size_t)m->nv);

	r = ss_mesh_radius(m, c);
	dr = (rref > 0.0) ? (rref - r) / rref : 0.0;
	for (int i = 0; i < m->nv; i++) {
		float *p = m->v + 3 * i;
		double ux = p[0] - c[0], uy = p[1] - c[1], uz = p[2] - c[2];
		double un = sqrt(ux * ux + uy * uy + uz * uz);
		if (un <= 0.0)
			continue;
		double dn = dr * un + un;
		p[0] = (float)(c[0] + ux / un * dn);
		p[1] = (float)(c[1] + uy / un * dn);
		p[2] = (float)(c[2] + uz / un * dn);
	}
	ss_mesh_normals(m);
}

#define SS_IS_LOWER_ZONE(p, c) (((p)[2] - (c)[2]) < 10.0f)
#define SS_IS_EYE_ZONE(p, c) ((((p)[1] - (c)[1]) < -10.0f) && (((p)[2] - (c)[2]) < 0.0f))

// Mirror of AFNI's SUMA_Find_IminImax: one walk down the inward ray and one up the
// outward ray, recording extremes, the DISTANCES at which they occur, the three Means,
// and the outward profile. The distances are what the touchup conditions are actually
// built on, so a version returning only min/max cannot express them.
//
// Four conventions here are load-bearing and each was a measured divergence:
//   1. sampling is NEAREST NEIGHBOUR (ss_sample_nn), never interpolated;
//   2. BOTH rays start at istep 0, so the value AT the node participates in Imin/Imax
//      and in Means[1] as well as Means[2];
//   3. the walk stops the first time a sample falls outside the grid, but that clamped
//      sample IS still used (AFNI increments istep only when still inside);
//   4. Means are divided by their counts UNCONDITIONALLY, so an empty count yields NaN
//      rather than 0. That difference is not cosmetic: with 0 the eye rule and touchup
//      cond4 both fire (0 < anything), with NaN both compare false. AFNI produces NaN.
typedef struct {
	float mm[2], mmd[2];        // min,max under the node and their distances
	float over[2], overd[2];    // min,max over the node and their distances
	float means[3];             // value at node, mean under, mean over
	float overshish[64];        // outward profile, 1 mm steps
	int n_over;
} ss_probe_t;

// AFNI's SUMA_MAX_SHISH_JUMP: walk the profile until the running max jump exceeds
// thresh_clip, then STOP. i_diffmax is the index BEFORE the jump -- it finds the FIRST
// significant jump, not the global maximum.
static int ss_max_shish_jump(const float *vec, int n, float thresh_clip) {
	float diffmax = 0.0f;
	int i_diffmax = 0;
	for (int k = 1; k < n; k++) {
		float d = vec[k] - vec[k - 1];
		if (d > diffmax) { diffmax = d; i_diffmax = k - 1; }
		if (diffmax > thresh_clip) break;
	}
	return i_diffmax;
}

// ---------------------------------------------------------------------------------------
// The ray walk and the deformation node loop are compiled TWICE, the same way coreFLT.c is
// compiled for float32 and float64: once faithful, once fast. They are not two algorithms --
// the kernel source is one file and the only differences are three things that cost time
// without changing what the code means.
//
// The reason it is a second COMPILATION rather than a runtime `if` is that the differences
// are all codegen, not arithmetic. This surface is knife-edge: the stage loop's convergence
// test is a discrete threshold on an INTEGER troubled-node count, so a sub-ULP change flips
// whether a pass converges and each extra pass adds N_it/2.5 = 100 more iterations. Measured:
// dropping the memset alone, or outlining the clamp alone, or letting OpenMP outline the node
// body -- each on its own moved the benchmark masks by Dice 0.970-0.996, all at IDENTICAL
// quality against AFNI (mean 0.9407 vs 0.9409 over the nine-image set). Nothing is wrong with
// those results; they are simply a different valid landing point. A runtime `if` inside one
// function would perturb BOTH paths and leave neither reproducing the shipped reference.
//
// So: SS_FAST is the default and is 1.8x quicker (17.4 s -> 9.6 s over the nine-image set); `-skullstrip -faithful` selects the other
// specialisation, which is bit-identical to the pre-optimisation release and is what the
// manifest's parity numbers and any regression diff should be taken against. Keep it that
// way -- if the two ever have to diverge in MEANING rather than in codegen, that is a bug.
// ---------------------------------------------------------------------------------------
#define SS_LINKAGE
#define SS_FN(name) name##_faithful
#include "skullstrip_kernel.h"
#undef SS_FN
#undef SS_LINKAGE

#define SS_FAST
#define SS_FN(name) name##_fast
#include "skullstrip_kernel.h"
#undef SS_FN
#undef SS_FAST

// Single dispatch point. `fast` is threaded down from the CLI, not read from the environment:
// a shipped op must not change its output because of a stray shell variable.
int ss_deform_range(const unsigned char *vol, int nx, int ny, int nz, const ss_stats *st,
		ss_mesh *m, float *ztv, const float *stop, int it0, int nit, int niter_total,
		int nnsmooth, float *maxexp_out, int fast) {
	if (fast)
		return ss_deform_range_fast(vol, nx, ny, nz, st, m, ztv, stop, it0, nit,
				niter_total, nnsmooth, maxexp_out);
	return ss_deform_range_faithful(vol, nx, ny, nz, st, m, ztv, stop, it0, nit,
			niter_total, nnsmooth, maxexp_out);
}

// Ventricle avoidance, AFNI's -avoid_vent, DEFAULT ON. Before a single deformation
// iteration runs, every node that is more than r/3 superior OR more than r/2 posterior
// of the centre is pushed radially outward by 1.1*d1 = 22 mm. It is a one-shot
// pre-stretch of the top and back of the starting sphere, not a force.
//
// This was the crack that resolved the manifest's standing contradiction. AFNI's trace
// reports mean segment length l = 3.242030 at it=0 while its node 1000 sits at exactly
// r/2 from the centre -- impossible for a pristine ld=20 icosphere, whose LONGEST edge
// at that radius is 3.032. Applying this pre-stretch to AFNI's own CreateIcosahedron
// output reproduces l = 3.242047. Since l scales BOTH the curvature radius (r = l^2/2|sn|)
// and the image force (su3 = ExpFrac*l*f3), missing it mis-scaled every iteration.
static void ss_prestretch_vent(ss_mesh *m, const ss_stats *st) {
	const float *c = st->cog;
	const float rr = st->radius;
	int ns = 0;
	for (int i = 0; i < m->nv; i++) {
		float *a = m->v + 3 * i;
		double ux, uy, uz, un;
		if (!((a[2] - c[2] > rr / 3.0f) || (a[1] - c[1] > rr / 2.0f)))
			continue;
		ux = a[0] - c[0]; uy = a[1] - c[1]; uz = a[2] - c[2];
		un = sqrt(ux * ux + uy * uy + uz * uz);
		if (un <= 0.0)
			continue;
		a[0] = (float)(c[0] + ux / un * (un + 1.1 * SS_D1_MM));
		a[1] = (float)(c[1] + uy / un * (un + 1.1 * SS_D1_MM));
		a[2] = (float)(c[2] + uz / un * (un + 1.1 * SS_D1_MM));
		ns++;
	}
	ss_mesh_normals(m);
	SSV("skullstrip: ventricle prestretch moved %d of %d nodes, mean seg %.6f\n",
			ns, m->nv, ss_mean_seg_len(m));
}

// ---------------------------------------------------------------------------
// Touchup. ADAPTED from AFNI's SUMA_Suggest_Touchup and SUMA_Reposition_Touchup: -touchup
// exists "to include areas not covered by surface expansion", and the stage loop iterates
// over "troubled nodes" until their count stops changing. The five conditions and the
// node-freezing branch are AFNI's; note in particular that cond2 here is the LIVE form and
// not the superseded one quoted in AFNI's own comment block above it.
// ---------------------------------------------------------------------------

static int ss_suggest_touchup(const unsigned char *vol, int nx, int ny, int nz,
		const ss_stats *st, ss_mesh *m, float *touch, float *stop);

// AFNI's SUMA_Reposition_Touchup: ask SUMA_Suggest_Touchup where the surface fell
// short, then move each troubled node OUTWARD ALONG ITS NORMAL by min(shift, limtouch).
// It is a SINGLE pass, not an iterated relaxation.
//
// Two details are load-bearing:
//   - a FROZEN node (Stop < 0, set by the five-condition detector) is skipped entirely;
//   - in the LOWER zone the shift is averaged with its neighbours' shifts, and if NO
//     neighbour also wants to move the node does not move at all ("only one node wants
//     to move in this hood"). Up top the raw per-node shift is used, because bumpy sulci
//     there are real structure rather than noise.
// AFNI calls this twice: limtouch 6, then Taubin smoothing, then limtouch 2.
int ss_reposition_touchup(const unsigned char *vol, int nx, int ny, int nz,
		const ss_stats *st, ss_mesh *m, float *stop, float limtouch) {
	float *touch = NULL;
	const float *ctr = st->cog;
	int n_troub;

	if (!vol || !st || !m)
		return -1;
	touch = (float *)malloc(sizeof(float) * (size_t)m->nv);
	if (!touch)
		return -1;
	n_troub = ss_suggest_touchup(vol, nx, ny, nz, st, m, touch, stop);
	if (n_troub <= 0) {
		free(touch);
		return 0;
	}
	for (int i = 0; i < m->nv; i++) {
		float *a = m->v + 3 * i;
		const float *n = m->nrm + 3 * i;
		float shft;
		if (stop && stop[i] < 0.0f)
			continue;
		if (!SS_IS_LOWER_ZONE(a, ctr))
			shft = touch[i];
		else {
			shft = touch[i];
			for (int q = m->nbr_off[i]; q < m->nbr_off[i + 1]; q++)
				shft += touch[m->nbr[q]];
			if (shft == touch[i])
				shft = 0.0f;
			else
				shft /= (float)(m->nbr_off[i + 1] - m->nbr_off[i] + 1);
		}
		if (shft != 0.0f) {
			float mv = shft < limtouch ? shft : limtouch;
			a[0] += mv * n[0];
			a[1] += mv * n[1];
			a[2] += mv * n[2];
		}
	}
	ss_mesh_normals(m);
	free(touch);
	return n_troub;
}

// Taubin lambda/mu smoothing with equal neighbour weights: even passes shrink by
// lambda, odd passes re-inflate by mu, so the surface is smoothed without the volume
// loss a plain Laplacian causes. AFNI's -smooth_final, 20 passes at AFNI's constants.
// scratch must hold 6*nv floats.
void ss_taubin_smooth(ss_mesh *m, int niter, float lambda, float mu, float *scratch) {
	float *a = scratch, *b = scratch + 3 * (size_t)m->nv;
	float *cur = m->v;
	for (int it = 0; it < niter; it++) {
		float w = (it & 1) ? mu : lambda;
		float *dst = (it & 1) ? a : b;
		for (int i = 0; i < m->nv; i++) {
			int deg = m->nbr_off[i + 1] - m->nbr_off[i];
			float dx = 0, dy = 0, dz = 0;
			for (int q = m->nbr_off[i]; q < m->nbr_off[i + 1]; q++) {
				dx += cur[3 * m->nbr[q]] - cur[3 * i];
				dy += cur[3 * m->nbr[q] + 1] - cur[3 * i + 1];
				dz += cur[3 * m->nbr[q] + 2] - cur[3 * i + 2];
			}
			if (deg) { dx /= deg; dy /= deg; dz /= deg; }
			dst[3 * i] = cur[3 * i] + w * dx;
			dst[3 * i + 1] = cur[3 * i + 1] + w * dy;
			dst[3 * i + 2] = cur[3 * i + 2] + w * dz;
		}
		cur = dst;
	}
	if (cur != m->v)
		memcpy(m->v, cur, sizeof(float) * 3 * (size_t)m->nv);
	ss_mesh_normals(m);
}




// Per-node "how far short is this node" in mm, written into touch[]; returns the
// count of troubled nodes. AFNI's SUMA_Suggest_Touchup applies five conditions and
// derives the shift from the distance to the minimum ABOVE the node; this is the
// distance-to-brain form of the same idea, kept deliberately short-reach so a node
// cannot leap a sulcus.
static int ss_suggest_touchup(const unsigned char *vol, int nx, int ny, int nz,
		const ss_stats *st, ss_mesh *m, float *touch, float *stop) {
	const int d1 = (int)SS_D1_MM, d4 = 15;
	const float *ctr = st->cog;   // fixed, as everywhere else -- never the running centroid
	int troub = 0;

	for (int i = 0; i < m->nv; i++) {
		const float *p = m->v + 3 * i, *n = m->nrm + 3 * i;
		ss_probe_t pr;
		float tb;
		int c1, c2, c3, c4;

		touch[i] = 0.0f;
		// Pinned to the FAITHFUL probe in both modes, deliberately -- this is not a leftover
		// from the kernel split, so do not "tidy" it into SS_FN(). ss_suggest_touchup probes
		// each node once per STAGE, not once per iteration (~10 calls per run against ~1250),
		// so the fast probe buys nothing measurable here, and pinning it keeps the touchup
		// decision -- which drives the stage-loop convergence test -- identical in both kernels.
		ss_probe_faithful(vol, nx, ny, nz, st, p, n, d1, d4, &pr);
		tb = (pr.mm[1] - st->t2) * 0.5f + st->t2;

		c1 = (pr.overd[0] < pr.mmd[0]) ||
		     (pr.mmd[0] > 1.0f && pr.mm[0] >= pr.over[0]) ||
		     (pr.mm[0] > tb && pr.over[0] < pr.mm[0]);
		// The LIVE cond2, not the version quoted in AFNI's own comment block above
		// it (that older form, `over[1] > mm[1] && over[1] > 0.9*t98 && ...`, is
		// commented out in the source and is NOT what runs).
		c2 = !(pr.over[1] > 1.2f * pr.mm[1] && pr.overd[1] < pr.overd[0]);
		c3 = !(pr.over[0] > 1.2f * pr.means[0]);
		c4 = !(pr.means[2] > pr.means[1]);

		if (!(pr.overd[0] > 0.0f && c1 && c2 && c3))
			continue;

		if (!c4) {
			// Brighter above than below: this is leaking into skull or fat. AFNI
			// FREEZES such a node rather than relaxing it -- Stop[in] = -1 forces
			// su3 = su4 = 0 for the rest of the run. Omitting this is what made the
			// ztv feedback harmful: without the freeze the loop only ever relaxes.
			if (stop)
				stop[i] = -1.0f;
			continue;
		}

		{
			// Outward gradient structure: Down is where the profile first drops
			// hard, Cross where it climbs back; the gap between them is the
			// thickness of the dark band being crossed.
			int Down = 0, Cross = 0;
			float gradthick = 0.0f;
			for (int k = 1; k < pr.n_over; k++) {
				float g = pr.overshish[k] - pr.overshish[k - 1];
				if ((g < -st->tm / 3.0f ||
						(pr.overshish[k] < st->tm / 2.0f && g < 0.0f)) && !Cross)
					Down = k;
				if (Down && g > st->tm / 3.0f && !Cross) {
					Cross = k;
					gradthick = (float)(Cross - Down);
				}
			}
			if (gradthick > 0.0f) {
				if (SS_IS_LOWER_ZONE(p, ctr)) {
					float cap = pr.overd[0] < 4.0f ? pr.overd[0] : 4.0f;
					touch[i] = (float)Down < cap ? (float)Down : cap;
				} else if (gradthick < 3.0f) {
					touch[i] = (float)Down < pr.overd[0] ? (float)Down : pr.overd[0];
				} else {
					touch[i] = pr.overd[0];
				}
			} else {
				// "No big jumps below the belt"; towards the top, go nuts.
				touch[i] = SS_IS_LOWER_ZONE(p, ctr)
				               ? (pr.overd[0] < 4.0f ? pr.overd[0] : 4.0f)
				               : pr.overd[0];
			}
			// Counted in EVERY branch, exactly as AFNI does -- including when the
			// computed shift is 0. The count drives the stage-loop convergence test,
			// so counting only nonzero shifts converges on a different schedule.
			troub++;
		}
	}
	return troub;
}

int ss_build_surface(const unsigned char *vol, int nx, int ny, int nz,
		const ss_stats *st, int ld, int niter, int max_retry,
		ss_mesh *m, long long *n_inter, int *n_tries, int fast) {
	float *ztv = NULL, *touch = NULL, *stop = NULL;
	// AFNI's stage machine, verbatim in structure (SUMA_StretchToFitLeCerveau's
	// do/while). It is NOT "250 iterations then touch up": stage 1 keeps ADDING
	// N_it/2.5 iterations while the surface is still moving (MaxExp > 0.5 mm) and its
	// area is still changing by more than 5%, and stage 2 adds another N_it/2.5 per
	// touchup round. The T1w run reaches ~1250 iterations, which is also what makes
	// AFNI's own runtime make sense.
	int it0 = 0, nit = niter, stage = 0, stage2type = 1, done = 0, keepgoing = 0;
	int past_troub = 0, n_troub = 0, npass = 0, rc = 1;
	double pastarea = 0.0, curarea, darea = 0.0;
	float maxexp = 0.0f;
	const int stage_cap = ss_env_int("SS_STAGE_MAX", 64, 1, 4096);
	// AFNI's self-intersection retry, and note WHERE it sits: the check happens on the
	// surface straight out of the expansion loop, BEFORE any repositioning. AFNI's own
	// comment says the later touchups "might cause some surface intersection, but their
	// effects should be small", so folds introduced after this point are accepted.
	int nnsmooth = SS_NNSMOOTH;
	int retry = (max_retry > 0) ? max_retry : SS_MAX_INTER_ITER;
	long long xs = 0, afni_xs = 0;

	if (!vol || !st || !m)
		return 1;

retry_surface:
	// goto, not return: on every pass after the first, ztv/touch/stop are already allocated
	// and the frees live below this point. rc is still 1 here.
	if (ss_icosphere(ld, 0.5f * st->radius, st->cog, m))
		goto done;
	ss_prestretch_vent(m, st);
	free(ztv); free(touch); free(stop);
	ztv = (float *)malloc(sizeof(float) * (size_t)m->nv);
	touch = (float *)malloc(sizeof(float) * (size_t)m->nv);
	stop = (float *)calloc((size_t)m->nv, sizeof(float));
	if (!ztv || !touch || !stop) {
		free(ztv); free(touch); free(stop); ss_mesh_free(m);
		return 1;
	}
	for (int i = 0; i < m->nv; i++)
		ztv[i] = SS_SHRINK_FAC;
	it0 = 0; nit = niter; stage2type = 1; done = 0; keepgoing = 0;
	past_troub = 0; n_troub = 0; npass = 0; pastarea = 0.0; darea = 0.0;

	do {
		stage = 0;
		if (ss_deform_range(vol, nx, ny, nz, st, m, ztv, stop, it0, nit, niter,
				nnsmooth, &maxexp, fast))
			goto done;
		++stage;
		if (stage == 1) {
			// Still growing? MaxExp is the largest single-node displacement in the LAST
			// iteration of the pass; the area test then asks whether that motion is
			// actually changing the surface or just jitter.
			if (maxexp > 0.5f) {
				if (pastarea == 0.0) {
					pastarea = ss_mesh_area(m);
					keepgoing = 1;
				} else {
					curarea = ss_mesh_area(m);
					darea = (curarea - pastarea) / pastarea;
					keepgoing = (darea < 0.0 ? -darea : darea) > 0.05;
					pastarea = curarea;
				}
				if (keepgoing) {
					it0 = nit;
					nit = nit + (int)(niter / 2.5);
					done = 0;
					SSV("skullstrip: stage 1: MaxExp %.3f darea %.4f -> iters %d..%d\n",
							(double)maxexp, darea, it0, nit);
				} else
					++stage;
			} else
				++stage;
		}
		if (stage == 2) {
			if (stage2type == 1) {
				n_troub = ss_suggest_touchup(vol, nx, ny, nz, st, m, touch, stop);
				if (!n_troub)
					++stage;
				else {
					// Relax the shrink factor exactly where the surface fell short.
					// The cut is PERMANENT -- a node relaxed in one stage stays relaxed.
					for (int i = 0; i < m->nv; i++) {
						if (touch[i] <= 0.0f)
							continue;
						if (touch[i] < 1.0f) ztv[i] *= 0.8f;
						else if (touch[i] < 2.0f) ztv[i] *= 0.7f;
						else if (touch[i] < 3.0f) ztv[i] *= 0.6f;
						else if (touch[i] < 4.0f) ztv[i] *= 0.5f;
						else ztv[i] *= 0.4f;
					}
				}
			} else {
				n_troub = 0;   // stage-2 type 2 is -push_to_edge, out of scope
				++stage;
			}
			it0 = nit;
			nit = nit + (int)(niter / 2.5);
			done = 0;
			if (!past_troub)
				past_troub = n_troub;
			else {
				double dtroub = (double)(past_troub - n_troub) / (double)past_troub;
				if (dtroub > 0.01)
					done = 0;
				else
					++stage;
				past_troub = n_troub;
			}
			SSV("skullstrip: stage 2 type %d: %d troubled -> iters %d..%d\n",
					stage2type, n_troub, it0, nit);
		}
		if (stage > 2) {
			if (stage2type < 2)
				stage2type = 2;
			else
				done = 1;
		}
		// AFNI's own "funding limit", so a surface that never converges still ends.
		if (nit > 8 * stage2type * niter && !done)
			done = 1;
		if (++npass > stage_cap)
			done = 1;
	} while (!done);

	// A NEGATIVE return is an allocation failure, NOT "no folds". Comparing `> 0` made a
	// single failed 37 KB calloc suppress AFNI's retry entirely and publish a folded surface
	// at exit 0 (measured: 7 self-intersections where the retry would have reached 0).
	// These allocations are tiny; if they fail the run is doomed anyway, so fail hard.
	// ONLY the AFNI-compatible predicate runs here, and only it decides the retry. The
	// exhaustive ss_mesh_self_intersections() is a REPORTING measure: running it on every
	// attempt cost a second spatial grid and a second full traversal per retry, and -- worse
	// -- its allocation failure aborted a production run for a number nobody had asked for.
	// It now runs once, at the end, on the mesh actually returned. That also fixes what it
	// described: it used to be measured BEFORE repositioning and then reported as the
	// returned surface's count, while the contract says repositioning may introduce folds.
	afni_xs = ss_afni_self_intersect(m);
	if (afni_xs < 0) {
		printfx("skullstrip: out of memory checking the surface for self-intersection\n");
		goto done;
	}
	if (afni_xs > 0 && retry > 0) {
		SSV("skullstrip: folds by AFNI's own test: %lld, retrying with %d "
				"smoothing passes\n", afni_xs, nnsmooth + SS_NNSMOOTH_STEP);
		nnsmooth += SS_NNSMOOTH_STEP;
		retry--;
		ss_mesh_free(m);
		goto retry_surface;
	}

	// Post-expansion, AFNI's own order: reposition (limit 6 mm) -> Taubin smoothing
	// -> reposition again (limit 2 mm). The two repositions are not redundant: the
	// first closes the large shortfalls, the smoothing removes the spikes it leaves,
	// and the second closes what the smoothing pulled back in.
	{
		float *scratch = (float *)malloc(sizeof(float) * 6 * (size_t)m->nv);
		const int smooth_end = SS_SMOOTH_END;   // AFNI -smooth_final
		if (!scratch)
			goto done;
		{
			// Call, THEN log. Never inline a call into SSV(): the macro guards its whole
			// body on ss_verbose(), so its arguments are not evaluated in an ordinary run.
			// Both of these were written that way, which silently skipped the entire
			// repositioning stage unless SKULLSTRIP_VERBOSE was set -- a diagnostic
			// variable that changed the segmentation. Caught in external audit.
			int t1, t2;
			t1 = ss_reposition_touchup(vol, nx, ny, nz, st, m, stop, 6.0f);
			if (t1 < 0) { free(scratch); goto done; }
			SSV("skullstrip: reposition 1: %d troubled\n", t1);
			if (smooth_end > 0)
				ss_taubin_smooth(m, smooth_end, 0.6307f, -0.6732f, scratch);
			t2 = ss_reposition_touchup(vol, nx, ny, nz, st, m, stop, 2.0f);
			if (t2 < 0) { free(scratch); goto done; }
			SSV("skullstrip: reposition 2: %d troubled\n", t2);
		}
		free(scratch);
	}
	// Measured on the RETURNED mesh, after repositioning, and only when the caller wants it.
	// A failure here is not fatal: it is a diagnostic, so report it as unknown rather than
	// discarding a good surface.
	if (n_inter) {
		xs = ss_mesh_self_intersections(m);
		*n_inter = xs;   // negative means "could not be measured"
	}
	if (n_tries) *n_tries = ((max_retry > 0) ? max_retry : SS_MAX_INTER_ITER) - retry + 1;
	rc = 0;
done:
	free(ztv);
	free(touch);
	free(stop);
	if (rc)
		ss_mesh_free(m);
	return rc;
}

int skullstrip_run(nifti_image *nim, int in_datatype, int fast) {
	ss_norm n;
	ss_stats st;
	ss_mesh m;
	unsigned char *wm = NULL;
	float *src = NULL, *out = NULL;
	long long nvox, wvox = (long long)SS_NX * SS_NY * SS_NZ;
	long long inmask = 0;
	float imin;
	long long xs = 0;
	int tries = 0, rc = 1;

	memset(&m, 0, sizeof(m));
	memset(&n, 0, sizeof(n));
	if (!nim || !nim->data)
		return 1;
	if (nim->datatype != DT_FLOAT32) {
		printfx("skullstrip: unsupported datatype %d (needs float32)\n", nim->datatype);
		return 1;
	}
	nvox = (long long)nim->nx * nim->ny * nim->nz;
	if (nvox <= 0 || (long long)nim->nvox != nvox) {
		printfx("skullstrip: requires a scalar 3D image (not 4D, not oversized)\n");
		return 1;
	}
	SSV("skullstrip: deformation kernel %s\n", fast ? "fast" : "faithful");
	src = (float *)nim->data;

	if (ss_normalize(nim, in_datatype, &n))
		return 1;
	if (ss_intensity_stats(n.vol, SS_NX, SS_NY, SS_NZ, &st))
		goto done;
	if (ss_build_surface(n.vol, SS_NX, SS_NY, SS_NZ, &st, SS_LD_NOEDGE, SS_NITER, 0,
			&m, ss_verbose() ? &xs : NULL, &tries, fast))
		goto done;
	// ONE full-size float buffer, not two. ss_restore writes the binary mask into `out`,
	// and the loop below then rewrites `out` IN PLACE from mask values to source-or-imin.
	// Fail-atomicity is unchanged -- nothing touches the caller's data until the final
	// memcpy -- but peak memory drops by a whole input volume (100 MB on a 192x256x256).
	wm = (unsigned char *)malloc((size_t)wvox);
	out = (float *)malloc(sizeof(float) * (size_t)nvox);
	if (!wm || !out)
		goto done;
	if (ss_mesh_rasterize(&m, SS_NX, SS_NY, SS_NZ, wm))
		goto done;
	// AFNI's -fill_hole (default 10 when touchup is on) is THD_mask_fillin_once, the
	// SAME directional gap-filler the normalisation already uses -- fill a background
	// voxel that has mask on BOTH sides within nside along any one axis. It is not a
	// connected-component hole fill, and the two disagree at the surface.
	if (ss_mask_fillin_once(SS_NX, SS_NY, SS_NZ, wm, 10))
		goto done;
	if (ss_restore(&n, wm, nim, out, 1)) // nearest neighbour: keeps the mask binary
		goto done;

	// Background fill = the image minimum over finite voxels (measured convention).
	imin = 0.0f;
	{
		int seen = 0;
		for (long long q = 0; q < nvox; q++) {
			float v = src[q];
			if (!isfinite(v))
				continue;
			if (!seen || v < imin) {
				imin = v;
				seen = 1;
			}
		}
		if (!seen)
			imin = 0.0f;
	}
	// In place: `out` holds the restored mask on entry and the published values on exit.
	// Read out[q] BEFORE overwriting it.
	for (long long q = 0; q < nvox; q++) {
		if (out[q] > 0.5f) {
			out[q] = src[q];
			inmask++;
		} else {
			out[q] = imin;
		}
	}
	if (inmask < 1) {
		printfx("skullstrip: produced an empty mask\n");
		goto done;
	}
	SSV("skullstrip: %lld voxels retained, background %g, %lld self-intersections\n",
			inmask, (double)imin, xs);

	// Commit only now: everything above could still have failed.
	memcpy(src, out, sizeof(float) * (size_t)nvox);
	rc = 0;

done:
	ss_mesh_free(&m);
	ss_norm_free(&n);
	free(wm);
	free(out);
	return rc;
}

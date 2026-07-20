/* Bias field correction (intensity uniformization)
   Adapted from AFNI's 3dUnifize by RW Cox (public domain)
   Original: https://github.com/afni/afni/blob/master/src/3dUnifize.c

   Faithful to 3dUnifize (matches its output to float rounding on real T1/T2/FLAIR):
   1. Automask (mri_automask_image): gradual (per-octant, trilinear) clip level at clfrac=0.2,
      largest 6-connected cluster, erode, recluster, fill holes, exterior_clip. In VOXELS.
   2. Optionally downsample by 2x ("duplo down", median-7) when nvox > 1e6, for speed
   3. For each voxel, local WM intensity = mean of the 70th-80th percentile within an 18.3-VOXEL
      sphere (0.5x radius in the duplo-downsampled grid)
   4. Upsample back to original resolution ("duplo up", trilinear)
   5. Scale each voxel by 1000/WMI to uniformize white matter; squash extreme highs with tanh
   6. Optional -GM: global gray-matter rescale (mri_GMunifize)
   7. do_mask (on by default in 3dUnifize): re-automask the unified output and zero non-brain */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "unifize.h"

#define PKVAL 1000.0f   /* target WM peak value */
#define PKMID  666.0f   /* target GM value for the -GM global scaling */
#define WMCUT 1300.0f   /* level for WM squashing */
#define WMSCL  200.0f   /* scale for WM squashing */

#define DEFAULT_RAD   18.3f  /* sphere radius in voxels */
#define DEFAULT_PBOT  70.0f  /* bottom percentile */
#define DEFAULT_PTOP  80.0f  /* top percentile */
#define CLFRAC         0.2f  /* automask clip fraction (3dUnifize sets clipfrac=0.2) */

// Comparator-free selection, replacing qsort()+cmp_float. qsort's comparator is a
// per-comparison indirect call — a severe WASM penalty. uf_select is used for the
// image-sized percentile and for per-voxel trimmed means without fully sorting the
// neighborhood.
static float uf_select(float *a, int n, int k) {
	if (n < 1) return 0.0f;
	if (k < 0) k = 0; else if (k >= n) k = n - 1;
	int lo = 0, hi = n - 1;
	while (lo < hi) {
		float pivot = a[lo + ((hi - lo) >> 1)];
		int lt = lo, gt = hi, i = lo;
		while (i <= gt) {
			if (a[i] < pivot) { float t = a[lt]; a[lt++] = a[i]; a[i++] = t; }
			else if (a[i] > pivot) { float t = a[gt]; a[gt--] = a[i]; a[i] = t; }
			else i++;
		}
		if (k < lt) hi = lt - 1; else if (k > gt) lo = gt + 1; else break;
	}
	return a[k];
}

static float thd_cliplevel(const float *data, int nvox, float mfrac);  /* fwd decl */
static float *cliplevel_gradual(const float *data, int nx, int ny, int nz, float mfrac);

/*--- AFNI THD_mask_clust: keep only the largest 6-connected (face) component. ---*/
static void mask_clust(uint8_t *m, int nx, int ny, int nz) {
	int nxy = nx * ny; int nxyz = nxy * nz;
	int *cur = (int *)malloc((size_t)nxyz * sizeof(int));
	int *best = (int *)malloc((size_t)nxyz * sizeof(int));
	if (!cur || !best) { free(cur); free(best); return; }
	int nbest = 0, last = 0;
	for (;;) {
		int seed; for (seed = last; seed < nxyz && !m[seed]; seed++) ;
		if (seed >= nxyz) break;
		last = seed + 1;
		int ncur = 0; m[seed] = 0; cur[ncur++] = seed;
		for (int c = 0; c < ncur; c++) {
			int p = cur[c], i = p % nx, j = (p / nx) % ny, k = p / nxy;
			if (i > 0      && m[p - 1  ]) { m[p - 1  ] = 0; cur[ncur++] = p - 1;   }
			if (i < nx - 1 && m[p + 1  ]) { m[p + 1  ] = 0; cur[ncur++] = p + 1;   }
			if (j > 0      && m[p - nx ]) { m[p - nx ] = 0; cur[ncur++] = p - nx;  }
			if (j < ny - 1 && m[p + nx ]) { m[p + nx ] = 0; cur[ncur++] = p + nx;  }
			if (k > 0      && m[p - nxy]) { m[p - nxy] = 0; cur[ncur++] = p - nxy; }
			if (k < nz - 1 && m[p + nxy]) { m[p + nxy] = 0; cur[ncur++] = p + nxy; }
		}
		if (ncur > nbest) { int *t = best; best = cur; cur = t; nbest = ncur; }
	}
	for (int c = 0; c < nbest; c++) m[best[c]] = 1;
	free(cur); free(best);
}

/*--- 18-neighbour (NN2: 6 faces + 12 edges) count with edge-clamped indices. ---*/
static int nn18_count(const uint8_t *m, int i, int j, int k, int nx, int ny, int nz, int nxy) {
	int kz = k * nxy, km = (k == 0) ? kz : kz - nxy, kp = (k == nz - 1) ? kz : kz + nxy;
	int jy = j * nx,  jm = (j == 0) ? jy : jy - nx,  jp = (j == ny - 1) ? jy : jy + nx;
	int im = (i == 0) ? 0 : i - 1, ip = (i == nx - 1) ? i : i + 1;
	return m[im + jy + km] + m[i + jm + km] + m[i + jy + km] + m[i + jp + km] + m[ip + jy + km]
	     + m[im + jm + kz] + m[im + jy + kz] + m[im + jp + kz]
	     + m[i + jm + kz]                    + m[i + jp + kz]
	     + m[ip + jm + kz] + m[ip + jy + kz] + m[ip + jp + kz]
	     + m[im + jy + kp] + m[i + jm + kp] + m[i + jy + kp] + m[i + jp + kp] + m[ip + jy + kp];
}

/*--- AFNI THD_mask_erodemany: peel npeel layers (a voxel with <17 of its 18 NN2 neighbours set
      is eroded), then redilate from the innermost layer outward. peelcount=1 snaps thin necks. ---*/
static void mask_erodemany(uint8_t *m, int nx, int ny, int nz, int npeel) {
	int nxy = nx * ny; int nxyz = nxy * nz;
	if (npeel < 1 || nxyz < 27) return;
	const int peelthr = 17;
	uint8_t *nnn = (uint8_t *)calloc(nxyz, 1);
	uint8_t *qqq = (uint8_t *)malloc(nxyz);
	if (!nnn || !qqq) { free(nnn); free(qqq); return; }
	for (int pp = 1; pp <= npeel; pp++) {
		for (int k = 0; k < nz; k++)
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++) {
					int v = i + j * nx + k * nxy;
					if (m[v] && nn18_count(m, i, j, k, nx, ny, nz, nxy) < peelthr) nnn[v] = (uint8_t)pp;
				}
		for (int v = 0; v < nxyz; v++) if (nnn[v]) m[v] = 0;
	}
	for (int pp = npeel; pp >= 1; pp--) {
		memset(qqq, 0, nxyz);
		uint8_t bth = (pp == npeel) ? 0 : 1;
		for (int k = 0; k < nz; k++)
			for (int j = 0; j < ny; j++)
				for (int i = 0; i < nx; i++) {
					int v = i + j * nx + k * nxy;
					if (nnn[v] >= pp && !m[v]) qqq[v] = (uint8_t)nn18_count(m, i, j, k, nx, ny, nz, nxy);
				}
		for (int v = 0; v < nxyz; v++) if (qqq[v] > bth) m[v] = 1;
	}
	free(qqq); free(nnn);
}

/*--- AFNI THD_mask_fillin_once: fill an unset voxel if, on any single axis, there is a set
      voxel within nside on BOTH the + and - side. Not an exterior flood. ---*/
static int mask_fillin_once(uint8_t *m, int nx, int ny, int nz, int nside) {
	int nxy = nx * ny; int nxyz = nxy * nz;
	int nsx = (nx - 1) / 2; if (nsx > nside) nsx = nside;
	int nsy = (ny - 1) / 2; if (nsy > nside) nsy = nside;
	int nsz = (nz - 1) / 2; if (nsz > nside) nsz = nside;
	if (nsx == 0 && nsy == 0 && nsz == 0) return 0;
	uint8_t *nnn = (uint8_t *)calloc(nxyz, 1);
	if (!nnn) return 0;
	int nfill = 0;
	for (int k = nsz; k < nz - nsz; k++)
		for (int j = nsy; j < ny - nsy; j++)
			for (int i = nsx; i < nx - nsx; i++) {
				int v = i + j * nx + k * nxy;
				if (m[v]) continue;
				int done = 0;
				if (nsx > 0) {
					int pl = 0, mi = 0;
					for (int l = 1; l <= nsx; l++) if (m[v + l])    { pl = 1; break; }
					if (pl) for (int l = 1; l <= nsx; l++) if (m[v - l]) { mi = 1; break; }
					if (pl && mi) { nnn[v] = 1; nfill++; done = 1; }
				}
				if (!done && nsy > 0) {
					int pl = 0, mi = 0;
					for (int l = 1; l <= nsy; l++) if (m[v + l * nx]) { pl = 1; break; }
					if (pl) for (int l = 1; l <= nsy; l++) if (m[v - l * nx]) { mi = 1; break; }
					if (pl && mi) { nnn[v] = 1; nfill++; done = 1; }
				}
				if (!done && nsz > 0) {
					int pl = 0, mi = 0;
					for (int l = 1; l <= nsz; l++) if (m[v + l * nxy]) { pl = 1; break; }
					if (pl) for (int l = 1; l <= nsz; l++) if (m[v - l * nxy]) { mi = 1; break; }
					if (pl && mi) { nnn[v] = 1; nfill++; }
				}
			}
	if (nfill > 0) for (int v = 0; v < nxyz; v++) if (nnn[v]) m[v] = 1;
	free(nnn);
	return nfill;
}
static int mask_fillin_completely(uint8_t *m, int nx, int ny, int nz, int nside) {
	int nfill = 0, k; do { k = mask_fillin_once(m, nx, ny, nz, nside); nfill += k; } while (k > 0);
	return nfill;
}

/*--- AFNI THD_mask_clip_neighbors: grow the (inverted, non-brain) mask into any below-clip
      voxel adjacent to it — trims the brain's low-intensity exterior shell (exterior_clip). ---*/
static int mask_clip_neighbors(uint8_t *m, const float *mar, int nx, int ny, int nz, float clip_val) {
	int nxy = nx * ny; float tclip = 9999.9f * clip_val;
	int ntot = 0, nnew;
	do {
		nnew = 0;
		for (int k = 1; k < nz - 1; k++) { int k3 = k * nxy;
			for (int j = 1; j < ny - 1; j++) { int j3 = k3 + j * nx;
				for (int i = 1; i < nx - 1; i++) { int p = i + j3;
					if (m[p] || (mar[p] >= clip_val && mar[p] <= tclip)) continue;
					if (m[p-1] || m[p+1] || m[p-nx] || m[p+nx] || m[p-nxy] || m[p+nxy]) { m[p] = 1; nnew++; }
				}
			}
		}
		ntot += nnew;
	} while (nnew > 0);
	return ntot;
}

/*--- AFNI mri_automask_image (3D): gradual clip -> largest cluster -> erode -> recluster ->
      fill holes -> erode -> recluster -> interior-hole fill (largest cluster of complement). ---*/
static uint8_t *compute_automask(const float *data, int nx, int ny, int nz) {
	int nxy = nx * ny; int nvox = nxy * nz;
	uint8_t *m = (uint8_t *)calloc(nvox, 1);
	if (!m) return NULL;
	float *car = cliplevel_gradual(data, nx, ny, nz, CLFRAC);
	int nmm = 0;
	if (car) {
		for (int i = 0; i < nvox; i++) if (data[i] >= car[i]) { m[i] = 1; nmm++; }
		free(car);
	} else {                                   /* gradual failed: fall back to fixed clip */
		float clip = thd_cliplevel(data, nvox, CLFRAC);
		for (int i = 0; i < nvox; i++) if (data[i] >= clip) { m[i] = 1; nmm++; }
	}
	if (nmm == 0 || nx < 2 || ny < 2 || nz < 2) return m;
	mask_clust(m, nx, ny, nz);
	mask_erodemany(m, nx, ny, nz, 1);
	mask_clust(m, nx, ny, nz);
	int ii = mask_fillin_once(m, nx, ny, nz, 1);
	if (ii > 0) { ii = mask_fillin_once(m, nx, ny, nz, 1);
	              if (ii > 0) mask_fillin_once(m, nx, ny, nz, 1); }
	int big = 1;
	int q = (int)lrintf(0.016f * nx); if (q > big) big = q;
	q = (int)lrintf(0.016f * ny); if (q > big) big = q;
	int qz = (int)lrintf(0.016f * nz); if (qz > big) big = qz;
	if (big > 1 || qz > 0) {
		for (int s = 2; s < big; s++) mask_fillin_once(m, nx, ny, nz, s);
		mask_fillin_completely(m, nx, ny, nz, big);
	}
	mask_erodemany(m, nx, ny, nz, 1);
	mask_clust(m, nx, ny, nz);
	/* interior-hole fill via the complement, with exterior_clip=1 (3dUnifize sets it): the
	   inverted (non-brain) mask is grown into below-clip neighbours to trim the low-intensity
	   exterior shell, then re-inverted; if anything was clipped, a final erode+cluster. */
	float clip_val = thd_cliplevel(data, nvox, CLFRAC);
	for (int i = 0; i < nvox; i++) m[i] = !m[i];
	mask_clust(m, nx, ny, nz);
	int jj = mask_clip_neighbors(m, data, nx, ny, nz, clip_val);
	for (int i = 0; i < nvox; i++) m[i] = !m[i];
	if (jj > 0) { mask_erodemany(m, nx, ny, nz, 1); mask_clust(m, nx, ny, nz); }
	return m;
}

/*--- Fast median of 7 values (sorting network) ---*/
#define SWAP(x, y) do { float t_ = x; x = y; y = t_; } while (0)
#define SORT2(a, b) do { if (a > b) SWAP(a, b); } while (0)

static float median7(float *p) {
	SORT2(p[0], p[1]); SORT2(p[4], p[5]); SORT2(p[1], p[2]);
	SORT2(p[5], p[6]); SORT2(p[0], p[1]); SORT2(p[4], p[5]);
	SORT2(p[0], p[4]); SORT2(p[2], p[6]); SORT2(p[1], p[3]);
	SORT2(p[3], p[5]); SORT2(p[1], p[3]); SORT2(p[2], p[3]);
	SORT2(p[3], p[4]); SORT2(p[2], p[3]); return p[3];
}

/*--- Downsample 3D image by 2x using median-of-7 ---*/
static float *duplo_down(const float *data, int nx, int ny, int nz,
	int *onx, int *ony, int *onz) {
	int nxg = nx / 2; if (nxg < 1) nxg = 1;
	int nyg = ny / 2; if (nyg < 1) nyg = 1;
	int nzg = nz / 2; if (nzg < 1) nzg = 1;
	int nxy = nx * ny;
	int nxyg = nxg * nyg;
	float *out = (float *)malloc((size_t)nxg * nyg * nzg * sizeof(float));
	if (!out) return NULL;
	for (int kk = 0; kk < nzg; kk++) {
		int ku = 2 * kk;
		int km = ku - 1; if (km < 0) km = 0;
		int kp = ku + 1; if (kp >= nz) kp = nz - 1;
		for (int jj = 0; jj < nyg; jj++) {
			int ju = 2 * jj;
			int jm = ju - 1; if (jm < 0) jm = 0;
			int jp = ju + 1; if (jp >= ny) jp = ny - 1;
			for (int ii = 0; ii < nxg; ii++) {
				int iu = 2 * ii;
				int im = iu - 1; if (im < 0) im = 0;
				int ip = iu + 1; if (ip >= nx) ip = nx - 1;
				float par[7];
				par[0] = data[iu + ju * nx + ku * nxy];
				par[1] = data[im + ju * nx + ku * nxy];
				par[2] = data[ip + ju * nx + ku * nxy];
				par[3] = data[iu + jm * nx + ku * nxy];
				par[4] = data[iu + jp * nx + ku * nxy];
				par[5] = data[iu + ju * nx + km * nxy];
				par[6] = data[iu + ju * nx + kp * nxy];
				out[ii + jj * nxg + kk * nxyg] = median7(par);
			}
		}
	}
	*onx = nxg; *ony = nyg; *onz = nzg;
	return out;
}

/*--- Upsample 3D image by 2x using trilinear averaging ---*/
static float *duplo_up(const float *data, int nx, int ny, int nz,
	int onx, int ony, int onz) {
	int nxy = nx * ny;
	int nxyo = onx * ony;
	float *out = (float *)calloc((size_t)onx * ony * onz, sizeof(float));
	if (!out) return NULL;
	for (int kk = 0; kk < onz; kk++) {
		int km = kk / 2, kp = km;
		if (km >= nz) km = kp = nz - 1;
		if (kk % 2) { kp = km + 1; if (kp >= nz) kp = nz - 1; }
		for (int jj = 0; jj < ony; jj++) {
			int jm = jj / 2, jp = jm;
			if (jm >= ny) jm = jp = ny - 1;
			if (jj % 2) { jp = jm + 1; if (jp >= ny) jp = ny - 1; }
			for (int ii = 0; ii < onx; ii++) {
				int im = ii / 2, ip = im;
				if (im >= nx) im = ip = nx - 1;
				if (ii % 2) { ip = im + 1; if (ip >= nx) ip = nx - 1; }
				out[ii + jj * onx + kk * nxyo] = 0.125f * (
					data[im + jm * nx + km * nxy] + data[ip + jm * nx + km * nxy] +
					data[im + jp * nx + km * nxy] + data[ip + jp * nx + km * nxy] +
					data[im + jm * nx + kp * nxy] + data[ip + jm * nx + kp * nxy] +
					data[im + jp * nx + kp * nxy] + data[ip + jp * nx + kp * nxy]);
			}
		}
	}
	return out;
}

/*--- Sphere neighborhood offsets ---*/
typedef struct { int di, dj, dk; } SphOff;

static SphOff *build_sphere(float radius, int *nout) {
	int r = (int)ceilf(radius);
	float r2 = radius * radius;
	int count = 0;
	for (int dk = -r; dk <= r; dk++)
		for (int dj = -r; dj <= r; dj++)
			for (int di = -r; di <= r; di++)
				if ((float)(di * di + dj * dj + dk * dk) <= r2) count++;
	SphOff *off = (SphOff *)malloc(count * sizeof(SphOff));
	if (!off) { *nout = 0; return NULL; }
	int idx = 0;
	for (int dk = -r; dk <= r; dk++)
		for (int dj = -r; dj <= r; dj++)
			for (int di = -r; di <= r; di++)
				if ((float)(di * di + dj * dj + dk * dk) <= r2) {
					off[idx].di = di;
					off[idx].dj = dj;
					off[idx].dk = dk;
					idx++;
				}
	*nout = count;
	return off;
}

/*--- Core: compute local percentile mean within sphere ---*/
static float *local_percmean(const float *data, const uint8_t *mask,
	int nx, int ny, int nz, float sph_radius, float pbot, float ptop) {
	int nxy = nx * ny, nxyz = nxy * nz;
	int nsph;
	SphOff *sph = build_sphere(sph_radius, &nsph);
	if (!sph) return NULL;
	float *out = (float *)calloc(nxyz, sizeof(float));
	if (!out) { free(sph); return NULL; }
	/* Per-voxel independent output -> OpenMP-parallel with a thread-local neighbourhood
	   buffer (AFNI's mri_local_percmean is likewise OMP-ized; this is the slow step). The
	   result is identical to serial. */
	/* `oom` is reduced across the team so a per-thread buffer failure propagates as a real
	   error (return NULL) instead of a silent all-zero result. EVERY thread must still encounter
	   the `omp for` worksharing construct (OpenMP requires uniform encounter, else deadlock at the
	   implicit barrier), so a thread that failed to allocate does not skip the loop — it enters it
	   and its iterations no-op (leaving out[v]=0) while `oom` carries the failure out. */
	int oom = 0;
#ifdef _OPENMP
	#pragma omp parallel reduction(|| : oom)
#endif
	{
		float *nbar = (float *)malloc((size_t)(nsph > 0 ? nsph : 1) * sizeof(float));
		if (!nbar) oom = 1;
#ifdef _OPENMP
		#pragma omp for schedule(dynamic, 4096)
#endif
		for (int v = 0; v < nxyz; v++) {
			if (!nbar) continue;   /* this thread OOM'd; failure carried by `oom` */
			int ii = v % nx, kk = v / nxy, jj = (v - kk * nxy) / nx;
			int ncount = 0;
			for (int s = 0; s < nsph; s++) {
				int ni = ii + sph[s].di, nj = jj + sph[s].dj, nk = kk + sph[s].dk;
				if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
					continue;
				int idx = ni + nj * nx + nk * nxy;
				if (!mask[idx]) continue;
				nbar[ncount++] = data[idx];
			}
			float val = 0.0f;
			if (ncount >= 2) {
				int q1 = (int)(0.01f * pbot * (ncount - 1));
				int q2 = (int)(0.01f * ptop * (ncount - 1));
				if (q2 > ncount - 1) q2 = ncount - 1;
				/* Trimmed mean of the [q1,q2] percentile band. Two O(n) quickselects
				   isolate that band into nbar[q1..q2] (select q2 over all, then q1 over
				   the lower part); we only sum it, so it needn't be sorted (summing in
				   partition rather than sorted order changes the estimate by <1 float ULP). */
				uf_select(nbar, ncount, q2);
				uf_select(nbar, q2 + 1, q1);
				for (int qq = q1; qq <= q2; qq++) val += nbar[qq];
				val /= (q2 - q1 + 1.0f);
			} else if (ncount == 1) {
				val = nbar[0];
			}
			out[v] = val;
		}
		free(nbar);
	}
	if (oom) { free(out); free(sph); return NULL; }
	free(sph);
	return out;
}

/*--- Median matching AFNI qmed_float: for even n, the average of the two middle
      order statistics (uf_select partitions so a[0..mid-1] <= a[mid]). Modifies `a`. ---*/
static float uf_qmed(float *a, int n) {
	if (n <= 0) return 0.0f;
	if (n == 1) return a[0];
	if (n == 2) return 0.5f * (a[0] + a[1]);
	int mid = n / 2;
	float m = uf_select(a, n, mid);
	if (n & 1) return m;
	float lo = a[0];
	for (int i = 1; i < mid; i++) if (a[i] > lo) lo = a[i];
	return 0.5f * (m + lo);
}

/*--- AFNI THD_cliplevel (float branch, thd_cliplevel.c): a histogram/median clip level
      to separate brain from background. `mfrac` scales the median-of-suprathreshold to the
      cut. Ignores negatives. Returns 0 when too few positive voxels. ---*/
static float thd_cliplevel(const float *data, int nvox, float mfrac) {
	if (mfrac <= 0.0f || mfrac >= 0.99f) mfrac = 0.50f;
	const int nhist = 10000;
	float fac = data[0];                       /* mri_max */
	for (int i = 1; i < nvox; i++) if (data[i] > fac) fac = data[i];
	if (fac < 1.0e-30f) return 0.0f;
	double sfac = (double)nhist / fac;
	int *hist = (int *)calloc((size_t)nhist + 1, sizeof(int));
	if (!hist) return 0.0f;
	double dsum = 0.0; int npos = 0;
	for (int i = 0; i < nvox; i++) {
		if (data[i] > 0.0f) {
			int kk = (int)(sfac * data[i] + 0.499);
			if (kk <= nhist) { hist[kk]++; dsum += (double)kk * (double)kk; npos++; }
		}
	}
	if (npos <= 222) { free(hist); return 0.0f; }
	/* start including the upper 65% of positive voxels */
	int qq = (int)(0.65f * npos);
	int ib = (int)rint(0.5 * sqrt(dsum / npos));
	int kk = 0, ii;
	for (ii = nhist - 1; ii >= ib && kk < qq; ii--) kk += hist[ii];
	/* median-adjustment: cut = mfrac * (median of the values above the cut) */
	int ncut = ii, nold;
	qq = 0;
	do {
		npos = 0;
		for (ii = ncut; ii < nhist; ii++) npos += hist[ii];
		int nhalf = npos / 2;
		kk = 0;
		for (ii = ncut; ii < nhist && kk < nhalf; ii++) kk += hist[ii];
		nold = ncut;
		ncut = (int)(mfrac * ii);
		qq++;
	} while (qq < 66 && ncut != nold);
	free(hist);
	double fclip = ncut / sfac;
	if (fclip > 1.0e38) fclip = 1.0e38;
	return (float)fclip;
}

/*--- Clip level of a rectangular sub-box [xa..xb, ya..yb, za..zb] (inclusive). ---*/
static float cliplevel_partial(const float *data, int nx, int nxy,
	int xa, int xb, int ya, int yb, int za, int zb, float mfrac) {
	int n = (xb - xa + 1) * (yb - ya + 1) * (zb - za + 1);
	if (n < 1) return 0.0f;
	float *sub = (float *)malloc((size_t)n * sizeof(float));
	if (!sub) return 0.0f;
	int t = 0;
	for (int k = za; k <= zb; k++)
		for (int j = ya; j <= yb; j++)
			for (int i = xa; i <= xb; i++) sub[t++] = data[i + j * nx + k * nxy];
	float c = thd_cliplevel(sub, n, mfrac);
	free(sub);
	return c;
}

/*--- AFNI THD_cliplevel_gradual: 8 octant clip levels about the intensity center-of-mass,
      trilinearly interpolated to a per-voxel clip image (gradualize=1, the automask default). ---*/
static float *cliplevel_gradual(const float *data, int nx, int ny, int nz, float mfrac) {
	int nxy = nx * ny; int nvox = nxy * nz;
	if (nvox < 1) return NULL;
	double sx = 0, sy = 0, sz = 0, sw = 0;
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++)
			for (int i = 0; i < nx; i++) {
				float w = data[i + j * nx + k * nxy];
				if (w > 0.0f) { sx += (double)w * i; sy += (double)w * j; sz += (double)w * k; sw += w; }
			}
	if (sw <= 0.0) return NULL;
	int it = nx - 1, jt = ny - 1, kt = nz - 1;
	int ic = (int)lrint(sx / sw), jc = (int)lrint(sy / sw), kc = (int)lrint(sz / sw);
	if (ic < 0) ic = 0; else if (ic > it) ic = it;
	if (jc < 0) jc = 0; else if (jc > jt) jc = jt;
	if (kc < 0) kc = 0; else if (kc > kt) kc = kt;
	float val = 0.333f * thd_cliplevel(data, nvox, mfrac);
	int di = (int)lrintf(0.01f * nx); if (di < 1) di = 1;
	int dj = (int)lrintf(0.01f * ny); if (dj < 1) dj = 1;
	int dk = (int)lrintf(0.01f * nz); if (dk < 1) dk = 1;
	int icm = ic - di; if (icm < 0) icm = 0; int icp = ic + di; if (icp > it) icp = it;
	int jcm = jc - dj; if (jcm < 0) jcm = 0; int jcp = jc + dj; if (jcp > jt) jcp = jt;
	int kcm = kc - dk; if (kcm < 0) kcm = 0; int kcp = kc + dk; if (kcp > kt) kcp = kt;
	float c000 = cliplevel_partial(data, nx, nxy, 0, icp, 0, jcp, 0, kcp, mfrac);
	float c100 = cliplevel_partial(data, nx, nxy, icm, it, 0, jcp, 0, kcp, mfrac);
	float c010 = cliplevel_partial(data, nx, nxy, 0, icp, jcm, jt, 0, kcp, mfrac);
	float c110 = cliplevel_partial(data, nx, nxy, icm, it, jcm, jt, 0, kcp, mfrac);
	float c001 = cliplevel_partial(data, nx, nxy, 0, icp, 0, jcp, kcm, kt, mfrac);
	float c101 = cliplevel_partial(data, nx, nxy, icm, it, 0, jcp, kcm, kt, mfrac);
	float c011 = cliplevel_partial(data, nx, nxy, 0, icp, jcm, jt, kcm, kt, mfrac);
	float c111 = cliplevel_partial(data, nx, nxy, icm, it, jcm, jt, kcm, kt, mfrac);
	if (c000 < val) c000 = val; if (c100 < val) c100 = val;
	if (c010 < val) c010 = val; if (c110 < val) c110 = val;
	if (c001 < val) c001 = val; if (c101 < val) c101 = val;
	if (c011 < val) c011 = val; if (c111 < val) c111 = val;
	float x0 = 0.5f * ic, x1 = 0.5f * (ic + it); float dxi = (x1 > x0) ? 1.0f / (x1 - x0) : 0.0f;
	float y0 = 0.5f * jc, y1 = 0.5f * (jc + jt); float dyi = (y1 > y0) ? 1.0f / (y1 - y0) : 0.0f;
	float z0 = 0.5f * kc, z1 = 0.5f * (kc + kt); float dzi = (z1 > z0) ? 1.0f / (z1 - z0) : 0.0f;
	float *car = (float *)malloc((size_t)nvox * sizeof(float));
	if (!car) return NULL;
	for (int k = 0; k < nz; k++) {
		float zt = (k - z0) * dzi; if (zt < 0) zt = 0; else if (zt > 1) zt = 1; float z0w = 1 - zt;
		for (int j = 0; j < ny; j++) {
			float yt = (j - y0) * dyi; if (yt < 0) yt = 0; else if (yt > 1) yt = 1; float y0w = 1 - yt;
			for (int i = 0; i < nx; i++) {
				float xt = (i - x0) * dxi; if (xt < 0) xt = 0; else if (xt > 1) xt = 1; float x0w = 1 - xt;
				car[i + j * nx + k * nxy] =
					c000 * x0w * y0w * z0w + c100 * xt * y0w * z0w +
					c010 * x0w * yt * z0w  + c110 * xt * yt * z0w +
					c001 * x0w * y0w * zt  + c101 * xt * y0w * zt +
					c011 * x0w * yt * zt   + c111 * xt * yt * zt;
			}
		}
	}
	return car;
}

/*--- AFNI mri_GMunifize (3dUnifize.c): after WM is unified to PKVAL, globally rescale so a
      'typical' GM intensity lands at PKMID, sharpening WM-GM contrast. Operates in place on
      the WM-unified image. No-op (leaves `data` unchanged) when too few voxels qualify. ---*/
static void gm_unifize(float *data, int nvox) {
	/* upper cutoff: reflect the median of super-peak (>PKVAL) values below PKVAL */
	int npval = 0;
	for (int i = 0; i < nvox; i++) if (data[i] > PKVAL) npval++;
	if (npval < 111) return;                      /* 1/6 of being beastly bad */
	float *pval = (float *)malloc((size_t)npval * sizeof(float));
	if (!pval) return;
	for (int i = 0, j = 0; i < nvox; i++) if (data[i] > PKVAL) pval[j++] = data[i];
	float pupper = uf_qmed(pval, npval);
	free(pval);
	pupper = PKVAL - 1.987654321f * (pupper - PKVAL);
	/* lower cutoff from the auto-clip level */
	float plower = thd_cliplevel(data, nvox, 0.4321f);
	/* median of the intermediate 'GM' voxels in [plower, pupper] */
	npval = 0;
	for (int i = 0; i < nvox; i++) if (data[i] >= plower && data[i] <= pupper) npval++;
	if (npval < 111) return;
	pval = (float *)malloc((size_t)npval * sizeof(float));
	if (!pval) return;
	for (int i = 0, j = 0; i < nvox; i++) if (data[i] >= plower && data[i] <= pupper) pval[j++] = data[i];
	float pmid = uf_qmed(pval, npval);
	free(pval);
	/* global linear scale: put pmid at PKMID while keeping WM's PKVAL fixed */
	float pfac = (PKVAL - PKMID) / (PKVAL - pmid);
	plower *= 0.333f;
	for (int i = 0; i < nvox; i++) {
		if (data[i] >= plower) {
			data[i] = pfac * (data[i] - PKVAL) + PKVAL;
			if (data[i] < 0.0f) data[i] = 0.0f;
		} else {
			data[i] = 0.0f;
		}
	}
}

/*--- Main unifize function ---*/
int unifize_image(float *data, int nx, int ny, int nz, float dx, float dy, float dz, int do_gm) {
	int nvox = nx * ny * nz;
	if (nvox < 100) return 1;
	(void)dx; (void)dy; (void)dz;
	/* AFNI's -Urad is in VOXELS (a fixed 18.3-voxel ball via MCW_spheremask(1,1,1,vrad)),
	   NOT millimetres — so use it directly, independent of voxel size, to match 3dUnifize. */
	float vrad = DEFAULT_RAD;
	int do_duplo = (nvox > 1000000);
	/* Step 1: Create working copy and automask */
	float *work = (float *)malloc(nvox * sizeof(float));
	if (!work) return 1;
	memcpy(work, data, nvox * sizeof(float));
	uint8_t *mask = compute_automask(work, nx, ny, nz);
	if (!mask) { free(work); return 1; }
	for (int i = 0; i < nvox; i++)
		if (!mask[i]) work[i] = 0.0f;
	free(mask);
	/* Step 2: Compute local WM intensity map */
	float *wmi;
	if (do_duplo) {
		int dnx, dny, dnz;
		float *down = duplo_down(work, nx, ny, nz, &dnx, &dny, &dnz);
		free(work);
		if (!down) return 1;
		int dnvox = dnx * dny * dnz;
		uint8_t *dmask = (uint8_t *)malloc(dnvox);
		if (!dmask) { free(down); return 1; }
		for (int i = 0; i < dnvox; i++) dmask[i] = (down[i] != 0.0f) ? 1 : 0;
		float *dwmi = local_percmean(down, dmask, dnx, dny, dnz,
			0.5f * vrad + 0.001f, DEFAULT_PBOT, DEFAULT_PTOP);
		free(down); free(dmask);
		if (!dwmi) return 1;
		int onx = 2 * dnx + (nx % 2);
		int ony = 2 * dny + (ny % 2);
		int onz = 2 * dnz + (nz % 2);
		wmi = duplo_up(dwmi, dnx, dny, dnz, onx, ony, onz);
		free(dwmi);
		if (!wmi) return 1;
	} else {
		uint8_t *mask2 = (uint8_t *)malloc(nvox);
		if (!mask2) { free(work); return 1; }
		for (int i = 0; i < nvox; i++) mask2[i] = (work[i] != 0.0f) ? 1 : 0;
		wmi = local_percmean(work, mask2, nx, ny, nz,
			vrad, DEFAULT_PBOT, DEFAULT_PTOP);
		free(work); free(mask2);
		if (!wmi) return 1;
	}
	/* Step 3: Scale input by 1000/WMI and squash extremes */
	for (int i = 0; i < nvox; i++) {
		float scale = (wmi[i] <= 0.0f) ? 0.0f : PKVAL / wmi[i];
		data[i] *= scale;
		if (data[i] > WMCUT)
			data[i] = WMCUT + WMSCL * tanhf((data[i] - WMCUT) / WMSCL);
	}
	free(wmi);
	/* Step 4 (optional -GM): global GM scaling on the WM-unified image (AFNI applies
	   mri_GMunifize after the WM unifize, before the final mask). */
	if (do_gm) gm_unifize(data, nvox);
	/* Step 5 (do_mask, on by default in 3dUnifize): automask the FINAL output and zero
	   non-brain — this is what keeps the result brain-only instead of the sphere-dilated
	   whole-head extent. The mask is recomputed on the unified output, not the input. */
	uint8_t *omask = compute_automask(data, nx, ny, nz);
	if (!omask) return 1;   /* fail closed: never emit a partially-unifized (unmasked) image */
	for (int i = 0; i < nvox; i++) if (!omask[i]) data[i] = 0.0f;
	free(omask);
	return 0;
}

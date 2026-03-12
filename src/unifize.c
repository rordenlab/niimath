/* Bias field correction (intensity uniformization)
   Adapted from AFNI's 3dUnifize by RW Cox (public domain)
   Original: https://github.com/afni/afni/blob/master/src/3dUnifize.c

   Algorithm:
   1. Automask the volume (simple clip-level threshold)
   2. Optionally downsample by 2x ("duplo down") for speed
   3. For each voxel, compute local WM intensity as the mean of
      the 70th-80th percentile of values within an 18.3-voxel sphere
   4. Upsample back to original resolution ("duplo up")
   5. Scale each voxel by 1000/WMI to uniformize white matter
   6. Squash extreme high values with tanh */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "unifize.h"

#define PKVAL 1000.0f   /* target WM peak value */
#define WMCUT 1300.0f   /* level for WM squashing */
#define WMSCL  200.0f   /* scale for WM squashing */

#define DEFAULT_RAD   18.3f  /* sphere radius in voxels */
#define DEFAULT_PBOT  70.0f  /* bottom percentile */
#define DEFAULT_PTOP  80.0f  /* top percentile */

static int cmp_float(const void *a, const void *b) {
	float fa = *(const float *)a, fb = *(const float *)b;
	return (fa > fb) - (fa < fb);
}

/*--- Automask: threshold at fraction of robust maximum ---*/
static uint8_t *compute_automask(const float *data, int nvox) {
	uint8_t *mask = (uint8_t *)calloc(nvox, sizeof(uint8_t));
	if (!mask) return NULL;
	int npos = 0;
	for (int i = 0; i < nvox; i++)
		if (data[i] > 0.0f) npos++;
	if (npos < 100) return mask;
	float *vals = (float *)malloc(npos * sizeof(float));
	if (!vals) { free(mask); return NULL; }
	for (int i = 0, j = 0; i < nvox; i++)
		if (data[i] > 0.0f) vals[j++] = data[i];
	qsort(vals, npos, sizeof(float), cmp_float);
	float p98 = vals[(int)(0.98f * (npos - 1))];
	float thresh = 0.10f * p98;
	free(vals);
	for (int i = 0; i < nvox; i++)
		mask[i] = (data[i] > thresh) ? 1 : 0;
	return mask;
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
	float *nbar = (float *)malloc(nsph * sizeof(float));
	float *out = (float *)calloc(nxyz, sizeof(float));
	if (!nbar || !out) {
		free(sph); free(nbar); free(out);
		return NULL;
	}
	for (int kk = 0; kk < nz; kk++) {
		for (int jj = 0; jj < ny; jj++) {
			for (int ii = 0; ii < nx; ii++) {
				int ncount = 0;
				for (int s = 0; s < nsph; s++) {
					int ni = ii + sph[s].di;
					int nj = jj + sph[s].dj;
					int nk = kk + sph[s].dk;
					if (ni < 0 || ni >= nx || nj < 0 || nj >= ny || nk < 0 || nk >= nz)
						continue;
					int idx = ni + nj * nx + nk * nxy;
					if (!mask[idx]) continue;
					nbar[ncount++] = data[idx];
				}
				float val = 0.0f;
				if (ncount >= 2) {
					qsort(nbar, ncount, sizeof(float), cmp_float);
					int q1 = (int)(0.01f * pbot * (ncount - 1));
					int q2 = (int)(0.01f * ptop * (ncount - 1));
					if (q2 > ncount - 1) q2 = ncount - 1;
					for (int qq = q1; qq <= q2; qq++) val += nbar[qq];
					val /= (q2 - q1 + 1.0f);
				} else if (ncount == 1) {
					val = nbar[0];
				}
				out[ii + jj * nx + kk * nxy] = val;
			}
		}
	}
	free(sph);
	free(nbar);
	return out;
}

/*--- Main unifize function ---*/
int unifize_image(float *data, int nx, int ny, int nz, float dx, float dy, float dz) {
	int nvox = nx * ny * nz;
	if (nvox < 100) return 1;
	/* Scale radius from mm to voxels using average voxel size */
	float voxsize = cbrtf(fabsf(dx * dy * dz));
	if (voxsize < 0.01f) voxsize = 1.0f;
	float vrad = DEFAULT_RAD / voxsize;
	if (vrad < 4.14f) vrad = 4.14f;
	int do_duplo = (nvox > 1000000);
	/* Step 1: Create working copy and automask */
	float *work = (float *)malloc(nvox * sizeof(float));
	if (!work) return 1;
	memcpy(work, data, nvox * sizeof(float));
	uint8_t *mask = compute_automask(work, nvox);
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
	return 0;
}

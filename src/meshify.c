#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER

#else
#include <unistd.h>
#endif
#include <time.h>
#ifdef HAVE_ZLIB
#include <zlib.h>
#ifdef HAVE_JSON
#include "cJSON.h"
#endif
#endif
#include "meshify.h"
#ifdef HAVE_FORMATS
#include "base64.h" //required for GIfTI
#endif
#include "bwlabel.h"
#include "meshtypes.h"
#include "radixsort.h"
#ifdef USE_CLASSIC_CUBES
#include "oldcubes.h"
#else
#include "MarchingCubes.h"
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

static double sqr(double x) {
	return x * x;
}

static double dx(vec3d p0, vec3d p1) {
	return sqrt(sqr(p0.x - p1.x) + sqr(p0.y - p1.y) + sqr(p0.z - p1.z));
}

static int unify_vertices(vec3d **inpt, vec3i *tris, int npt, int ntri, bool verbose) {
	double startTime = clockMsec();
	vec3d *pts = *inpt;
	float *dx_in = (float *)malloc(npt * sizeof(float));
	float *dx_out = (float *)malloc(npt * sizeof(float));
	uint32_t *idx_in = (uint32_t *)malloc(npt * sizeof(uint32_t));
	uint32_t *idx_out = (uint32_t *)malloc(npt * sizeof(uint32_t));
	int *old2new = (int *)malloc(npt * sizeof(int));
	if (!dx_in || !dx_out || !idx_in || !idx_out || !old2new) {
		free(dx_in); free(dx_out); free(idx_in); free(idx_out); free(old2new);
		return npt;
	}
	vec3d ref = pts[0];
	for (int i = 0; i < npt; i++) {
		dx_in[i] = dx(ref, pts[i]);
		idx_in[i] = i;
		old2new[i] = -1;
	}
	radix11sort_f32(dx_in, dx_out, idx_in, idx_out, npt);
	free(dx_in);
	free(idx_in);
	const float tol = 0.00001f;
	int nnew = 0;
	for (int i = 0; i < npt; i++) {
		int vi = idx_out[i];
		if (old2new[vi] >= 0)
			continue;
		vec3d pti = pts[vi];
		float di = dx_out[i];
		old2new[vi] = nnew;
		for (int j = i + 1; j < npt; j++) {
			if ((dx_out[j] - di) >= tol)
				break;
			int vj = idx_out[j];
			if (old2new[vj] < 0 && dx(pti, pts[vj]) < tol)
				old2new[vj] = nnew;
		}
		nnew++;
	}
	free(dx_out);
	free(idx_out);
	if (npt == nnew) {
		if (verbose)
			printf("Unify vertices found no shared vertices\n");
		free(old2new);
		return npt;
	}
	for (int i = 0; i < ntri; i++) {
		tris[i].x = old2new[tris[i].x];
		tris[i].y = old2new[tris[i].y];
		tris[i].z = old2new[tris[i].z];
	}
	vec3d *newpts = (vec3d *)calloc(nnew, sizeof(vec3d));
	char *written = (char *)calloc(nnew, sizeof(char));
	if (!newpts || !written) {
		free(newpts); free(written); free(dx_out); free(idx_out); free(old2new);
		return npt;
	}
	for (int i = 0; i < npt; i++) {
		int ni = old2new[i];
		if (!written[ni]) {
			newpts[ni] = pts[i];
			written[ni] = 1;
		}
	}
	free(written);
	free(*inpt);
	*inpt = newpts;
	free(old2new);
	if (verbose)
		printf("vertex welding %d -> %d: %ld ms\n", npt, nnew, timediff(startTime, clockMsec()));
	return nnew;
}


#ifndef FLT_EPSILON
#define FLT_EPSILON 1.19209290e-07F // float
// #define DBL_EPSILON 2.2204460492503131e-16 // double
#endif

static int remove_degenerate_triangles(vec3d *pts, vec3i **intris, int ntri, bool verbose) {
	// reduces the number of triangles, number of vertices unchanged
	double startTime = clockMsec();
	vec3i *tris = *intris;
	int *isdegenerate = (int *)malloc(ntri * sizeof(int));
	if (!isdegenerate) return ntri;
	int ndegenerate = 0;
	for (int i = 0; i < ntri; i++) {
		// sorted lengths a ≥ b ≥ c
		isdegenerate[i] = 0;
		double l = dx(pts[tris[i].x], pts[tris[i].y]);
		double m = dx(pts[tris[i].x], pts[tris[i].z]);
		double n = dx(pts[tris[i].y], pts[tris[i].z]);
		double c = fmin(fmin(l, m), n);
		double a = fmax(fmax(l, m), n);
		double b = l + m + n - a - c;
		if ((c - (a - b)) <= 0.0) {
			isdegenerate[i] = 1;
			ndegenerate++;
			continue;
		}
		if (tris[i].x == tris[i].y || tris[i].x == tris[i].z || tris[i].y == tris[i].z) {
				isdegenerate[i] = 1;
				ndegenerate++;
				continue;
		}
#define REQUIRE_SIGNIFICANT_AREA
#ifdef REQUIRE_SIGNIFICANT_AREA
		// use Heron’s Formula to eliminate triangles of tiny area
		//  see Kahan: Miscalculating Area and Angles of a Needle-like Triangle
		//  https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
		// n.b. scale area improves watertight check
		double area4 = 0.25 * sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)));
		if (area4 < 1e-10 * (a * b)) {
			isdegenerate[i] = 1;
			ndegenerate++;
		}
#endif // REQUIRE_SIGNIFICANT_AREA
	}
	if (ndegenerate == 0) {
		free(isdegenerate);
		return ntri;
	}
	int newtri = ntri - ndegenerate;
	vec3i *oldtris = (vec3i *)malloc(ntri * sizeof(vec3i));
	if (!oldtris) { free(isdegenerate); return ntri; }
	for (int i = 0; i < ntri; i++)
		oldtris[i] = tris[i];
	free(*intris);
	*intris = (vec3i *)malloc(newtri * sizeof(vec3i));
	if (!*intris) { free(oldtris); free(isdegenerate); *intris = oldtris; return ntri; }
	tris = *intris;
	int j = 0;
	for (int i = 0; i < ntri; i++) {
		if (isdegenerate[i])
			continue;
		tris[j] = oldtris[i];
		j++;
	}
	free(oldtris);
	free(isdegenerate);
	if (verbose)
		printf("remove degenerate triangles %d -> %d: %ld ms\n", ntri, newtri, timediff(startTime, clockMsec()));
	return newtri;
}

static int quick_smooth(float *img, int nx, int ny, int nz) {
	if ((nx < 5) || (ny < 5) || (nz < 5))
		return EXIT_FAILURE;
	int nvox = nx * ny * nz;
	float *tmp = (float *)malloc(nvox * sizeof(float));
	if (!tmp) return EXIT_FAILURE;
#define kwid 2
#define k0 0.45
#define k1 0.225
#define k2 0.05
	int nxy = nx * ny;
	int nxy2 = nxy * 2;
	int nx2 = nx * 2;
	// smooth column direction
	memcpy(tmp, img, nvox * sizeof(float)); // dst,src,n
	for (int z = 0; z < nz; z++) {
		for (int y = 0; y < ny; y++) {
			int zy = (y * nx) + (z * nxy);
			for (int x = kwid; x < (nx - kwid); x++) {
				int v = zy + x;
				tmp[v] = (img[v - 2] * k2) + (img[v - 1] * k1) + (img[v] * k0) + (img[v + 1] * k1) + (img[v + 2] * k2);
			}
		}
	}
	// smooth row direction:
	memcpy(img, tmp, nvox * sizeof(float)); // dst,src,n
	for (int z = 0; z < nz; z++) {
		for (int x = 0; x < nx; x++) {
			int xz = x + (z * nxy);
			for (int y = kwid; y < (ny - kwid); y++) {
				int v = xz + (y * nx);
				tmp[v] = (img[v - nx2] * k2) + (img[v - nx] * k1) + (img[v] * k0) + (img[v + nx] * k1) + (img[v + nx2] * k2);
			}
		}
	}
	memcpy(img, tmp, nvox * sizeof(float)); // dst,src,n
	for (int y = 0; y < ny; y++) {
		for (int x = 0; x < nx; x++) {
			int yx = (y * nx) + x;
			for (int z = kwid; z < (nz - kwid); z++) {
				int v = yx + (z * nxy);
				img[v] = (tmp[v - nxy2] * k2) + (tmp[v - nxy] * k1) + (tmp[v] * k0) + (tmp[v + nxy] * k1) + (tmp[v + nxy2] * k2);
			}
		}
	}
	free(tmp);
	return EXIT_SUCCESS;
}

static void dilate(float *img, size_t dim[3], bool is26) {
	int nx = dim[0];
	int ny = dim[1];
	int nz = dim[2];
	int nxy = nx * ny;
	int nvox = nx * ny * nz;
	uint8_t *mask = (uint8_t *)malloc(nvox * sizeof(uint8_t));
	if (!mask) return;
	memset(mask, 0, nvox * sizeof(uint8_t));
	int numk = 6;
	if (is26)
		numk = 26;
	int32_t *k = (int32_t *)malloc(numk * sizeof(int32_t)); // queue with untested seed
	if (!k) { free(mask); return; }
	if (is26) {
		int j = 0;
		for (int z = -1; z <= 1; z++)
			for (int y = -1; y <= 1; y++)
				for (int x = -1; x <= 1; x++) {
					if ((x == 0) && (y == 0) && (z == 0))
						continue;
					k[j] = x + (y * nx) + (z * nx * ny);
					j++;
				} // for x
	} else {			// if 26 neighbors else 6..
		k[0] = nx * ny; // up
		k[1] = -k[0];	// down
		k[2] = nx;		// anterior
		k[3] = -k[2];	// posterior
		k[4] = 1;		// left
		k[5] = -1;
	}
	for (int z = 1; z < (nz - 1); z++) {
		for (int y = 1; y < (ny - 1); y++) {
			size_t iyz = +(z * nxy) + (y * nx);
			for (int x = 1; x < (nx - 1); x++) {
				size_t vx = iyz + x;
				for (int n = 1; n < numk; n++) {
					if (img[vx + k[n]] > 0)
						mask[vx] = 1;
				} // check all neighbors
			} // x
		} // y
	} // z
	for (int v = 1; v < nvox; v++)
		if (mask[v] > 0)
			img[v] = 1;
	free(mask);
	free(k);
}

double clockMsec(void) { // return milliseconds since midnight
#ifdef _MSC_VER
	clock_t t = clock();
	return (double)((double)t) / (CLOCKS_PER_SEC / 1000.0);
#else
#ifdef __MINGW32__ // issue 4
	time_t seconds_since_midnight = time(NULL) % 86400;
	return seconds_since_midnight;
#else
	struct timespec _t;
	clock_gettime(CLOCK_MONOTONIC, &_t);
	return _t.tv_sec * 1000.0 + (_t.tv_nsec / 1.0e6);
#endif
#endif
}

long timediff(double startTimeMsec, double endTimeMsec) {
	return round(endTimeMsec - startTimeMsec);
}

int meshify(float *img, short dim[3], int originalMC, float isolevel, vec3i **t, vec3d **p, int *nt, int *np, bool preSmooth, bool onlyLargest, bool fillBubbles, bool verbose) {
	// img: input volume
	// hdr: nifti header
	// isolevel: air/surface threshold
	// t: triangle indices e.g. [0,1,3] indicates triangle composed of vertices 0,1,3
	// p: 3D points, aka vertices
	// nt: number of triangles, aka faces
	// np: number of points
	int NX = dim[0];
	int NY = dim[1];
	int NZ = dim[2];
	int nvox = NX * NY * NZ;
	// preSmooth: Gaussian blur to soften image
	if (preSmooth) {
		double startTime = clockMsec();
		quick_smooth(img, NX, NY, NZ);
		if (verbose)
			printf("pre-smooth: %ld ms\n", timediff(startTime, clockMsec()));
	}
	// determine image intensity range - ensure isolvel will detect edge
	float mx = img[0];
	float mn = mx;
	for (int i = 0; i < nvox; i++) {
		mx = fmaxf(mx, img[i]);
		mn = fminf(mn, img[i]);
	}
	if (mn == mx) {
		printf("Error: No variability in image intensity.\n");
		return EXIT_FAILURE;
	}
	if ((isolevel <= mn) || (isolevel > mx)) {
		isolevel = 0.5 * (mn + mx);
		printf("Suggested isolevel out of range. Intensity range %g..%g, setting isolevel to %g\n", mn, mx, isolevel);
	}
	if (verbose)
		printf("intensity range %g..%g, isolevel %g\n", mn, mx, isolevel);
	//(optional) fill bubbles and only extract largest contiguous object
	if ((onlyLargest) || (fillBubbles)) {
		double startTime = clockMsec();
		float *mask = (float *)malloc(nvox * sizeof(float));
		if (!mask) return EXIT_FAILURE;
		size_t dim[3] = {(size_t)NX, (size_t)NY, (size_t)NZ};
		memset(mask, 0, nvox * sizeof(float));
		for (int i = 0; i < nvox; i++)
			if (img[i] >= isolevel)
				mask[i] = 1;
		bwlabel(mask, 18, dim, onlyLargest, fillBubbles);
		if (fillBubbles) {
			for (int i = 0; i < nvox; i++)
				if (mask[i] != 0)
					img[i] = fmax(img[i], isolevel);
		}
		if (onlyLargest) {
			dilate(mask, dim, true); // expand by one voxel to preserve subvoxel edges
			for (int i = 0; i < nvox; i++)
				if (mask[i] == 0)
					img[i] = mn;
		}
		free(mask);
		if (verbose)
			printf("voxel clustering (largest cluster, bubbles): %ld ms\n", timediff(startTime, clockMsec()));
	}
	// edge darken
	float edgeMax = 0.75 * (mn + isolevel);
	int vx = 0;
	int lo[3] = {NX, NY, NZ};
	int hi[3] = {0, 0, 0};
	for (int z = 0; z < NZ; z++) // darken edges
		for (int y = 0; y < NY; y++)
			for (int x = 0; x < NX; x++) {
				if (img[vx] >= isolevel) {
					lo[0] = MIN(x, lo[0]);
					lo[1] = MIN(y, lo[1]);
					lo[2] = MIN(z, lo[2]);
					hi[0] = MAX(x, hi[0]);
					hi[1] = MAX(y, hi[1]);
					hi[2] = MAX(z, hi[2]);
				}
				if ((x == 0) || (y == 0) || (z == 0) || (x == (NX - 1)) || (y == (NY - 1)) || (z == (NZ - 1)))
					img[vx] = fminf(edgeMax, img[vx]);
				vx++;
			}
	// printf("Bounding box for bright voxels: %d..%d %d..%d %d..%d\n", lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]);
	for (int i = 0; i < 3; i++) {
		lo[i] = MAX(lo[i] - 1, 0);
		hi[i] = MIN(hi[i] + 2, dim[i]);
	}
	double startTimeMC = clockMsec();
	vec3d *pts = NULL;
	vec3i *tris = NULL;
	int ntri;
	int npt;
	if (marchingCubes(img, dim, lo, hi, originalMC, isolevel, &pts, &tris, &npt, &ntri) != EXIT_SUCCESS)
		return EXIT_FAILURE;
	if (verbose)
		printf("marching cubes (%dx%dx%d): %ld ms\n", NX, NY, NZ, timediff(startTimeMC, clockMsec()));
	npt = unify_vertices(&pts, tris, npt, ntri, verbose);
	if (npt < 3)
		return EXIT_FAILURE;
	ntri = remove_degenerate_triangles(pts, &tris, ntri, verbose);
	*t = tris;
	*p = pts;
	*nt = ntri;
	*np = npt;
	return EXIT_SUCCESS;
}

/* ================================ mesh quality report =======================================
 * Topology from an edge hash (holes = boundary edges, non-manifold = edges with >2 faces, vertex
 * umbrellas that split into more than one fan), connectivity and genus from Euler's formula, and
 * self-intersections by Moller's 1997 triangle-triangle test on a uniform grid.  No dependencies.
 * Counts TRIANGLES that intersect, not pairs, so grid double-visits need no dedupe. */

typedef struct { uint64_t key; int32_t n; int32_t t0, t1; } mc_edge;   /* t0/t1: first two faces */

static uint64_t mc_ekey(int a, int b) {
	return a < b ? ((uint64_t)a << 32) | (uint32_t)b : ((uint64_t)b << 32) | (uint32_t)a;
}

static mc_edge *mc_efind(mc_edge *h, size_t cap, uint64_t key) {   /* open addressing, key 0 = empty */
	size_t i = (size_t)((key * 0x9E3779B97F4A7C15ull) >> 20) & (cap - 1);
	while (h[i].key && h[i].key != key) i = (i + 1) & (cap - 1);
	return h + i;
}

static int mc_find(int *uf, int i) { while (uf[i] != i) { uf[i] = uf[uf[i]]; i = uf[i]; } return i; }
static void mc_union(int *uf, int a, int b) { a = mc_find(uf, a); b = mc_find(uf, b); if (a != b) uf[a] = b; }

/* Moller: do triangles (v0,v1,v2) and (u0,u1,u2) intersect?  Coplanar pairs are ignored --
   ponytail: a coplanar overlap on a marching-cubes surface is a degenerate fold, vanishingly rare,
   and the 2D case is a page of code for it. */
static void mc_sub(const double *a, const double *b, double *o) { o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2]; }
static void mc_cross(const double *a, const double *b, double *o) {
	o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0]; }
static double mc_dot(const double *a, const double *b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

/* interval of a triangle's projection onto the intersection line, given signed distances d */
static int mc_interval(double p0, double p1, double p2, double d0, double d1, double d2, double *lo, double *hi) {
	double a, b;
	if (d0 * d1 > 0.0) { a = p2 + (p0 - p2) * d2 / (d2 - d0); b = p2 + (p1 - p2) * d2 / (d2 - d1); }
	else if (d0 * d2 > 0.0) { a = p1 + (p0 - p1) * d1 / (d1 - d0); b = p1 + (p2 - p1) * d1 / (d1 - d2); }
	else if (d1 * d2 > 0.0 || d0 != 0.0) { a = p0 + (p1 - p0) * d0 / (d0 - d1); b = p0 + (p2 - p0) * d0 / (d0 - d2); }
	else if (d1 != 0.0) { a = p1 + (p0 - p1) * d1 / (d1 - d0); b = p1 + (p2 - p1) * d1 / (d1 - d2); }
	else if (d2 != 0.0) { a = p2 + (p0 - p2) * d2 / (d2 - d0); b = p2 + (p1 - p2) * d2 / (d2 - d1); }
	else return 0;   /* coplanar */
	*lo = MIN(a, b); *hi = MAX(a, b);
	return 1;
}

int mesh_tri_tri(const double *v0, const double *v1, const double *v2,
	const double *u0, const double *u1, const double *u2) {
	double e1[3], e2[3], n1[3], n2[3], d[3], dv[3], du[3], lo1, hi1, lo2, hi2, mx;
	int i;
	mc_sub(v1, v0, e1); mc_sub(v2, v0, e2); mc_cross(e1, e2, n1);
	{ double d1 = -mc_dot(n1, v0);
	  dv[0] = mc_dot(n1, u0) + d1; dv[1] = mc_dot(n1, u1) + d1; dv[2] = mc_dot(n1, u2) + d1; }
	if ((dv[0] > 0 && dv[1] > 0 && dv[2] > 0) || (dv[0] < 0 && dv[1] < 0 && dv[2] < 0)) return 0;
	mc_sub(u1, u0, e1); mc_sub(u2, u0, e2); mc_cross(e1, e2, n2);
	{ double d2 = -mc_dot(n2, u0);
	  du[0] = mc_dot(n2, v0) + d2; du[1] = mc_dot(n2, v1) + d2; du[2] = mc_dot(n2, v2) + d2; }
	if ((du[0] > 0 && du[1] > 0 && du[2] > 0) || (du[0] < 0 && du[1] < 0 && du[2] < 0)) return 0;
	mc_cross(n1, n2, d);
	mx = fabs(d[0]); i = 0;
	if (fabs(d[1]) > mx) { mx = fabs(d[1]); i = 1; }
	if (fabs(d[2]) > mx) i = 2;
	if (!mc_interval(v0[i], v1[i], v2[i], du[0], du[1], du[2], &lo1, &hi1)) return 0;
	if (!mc_interval(u0[i], u1[i], u2[i], dv[0], dv[1], dv[2], &lo2, &hi2)) return 0;
	return hi1 >= lo2 && hi2 >= lo1;
}

/* Self-intersecting triangles: uniform grid over triangle AABBs, Moller tri-tri on cell-mates that
   share no vertex (coplanar overlap is deliberately not counted).  hit[t] = 1 for each crossing
   triangle.  Returns their count, or -1 when the grid did not fit in memory. */
int mesh_self_intersections(vec3i *tris, vec3d *pts, int ntri, int npt, uint8_t *hit) {
	int nself = 0;
	if (npt < 1) return 0;
	memset(hit, 0, (size_t)ntri);
	double lo[3] = { pts[0].x, pts[0].y, pts[0].z }, hi[3] = { lo[0], lo[1], lo[2] }, elen = 0.0, cell;
	for (int k = 0; k < 3; k++) if (!(fabs(lo[k]) <= DBL_MAX)) return -1;
	int g[3], *cnt = NULL, *cell_tri = NULL;
	int64_t ncell;
	for (int i = 1; i < npt; i++) {
		double p[3] = { pts[i].x, pts[i].y, pts[i].z };
		for (int k = 0; k < 3; k++) { if (!(fabs(p[k]) <= DBL_MAX)) return -1; lo[k] = MIN(lo[k], p[k]); hi[k] = MAX(hi[k], p[k]); }   /* NaN/Inf: no grid */
	}
	for (int t = 0; t < ntri; t++) elen += dx(pts[tris[t].x], pts[tris[t].y]);
	cell = ntri ? 2.0 * elen / ntri : 1.0;
	if (!(cell > 0.0)) cell = 1.0;
	for (int k = 0; k < 3; k++) { double c = (hi[k] - lo[k]) / 256 + 1e-9; if (c > cell) cell = c; }   /* <= 257^3 cells */
	for (int k = 0; k < 3; k++) g[k] = (int)((hi[k] - lo[k]) / cell) + 1;
	ncell = (int64_t)g[0] * g[1] * g[2];
	cnt = (int *)calloc((size_t)ncell + 1, sizeof(int));
	if (!cnt) return -1;
	{
		#define MC_CELLS(t, c0, c1) do { \
			int _v[3] = { tris[t].x, tris[t].y, tris[t].z }; \
			double _p[3][3] = { { pts[_v[0]].x, pts[_v[0]].y, pts[_v[0]].z }, { pts[_v[1]].x, pts[_v[1]].y, pts[_v[1]].z }, { pts[_v[2]].x, pts[_v[2]].y, pts[_v[2]].z } }; \
			for (int _k = 0; _k < 3; _k++) { \
				double _a = MIN(_p[0][_k], MIN(_p[1][_k], _p[2][_k])), _b = MAX(_p[0][_k], MAX(_p[1][_k], _p[2][_k])); \
				c0[_k] = MAX(0, MIN(g[_k] - 1, (int)((_a - lo[_k]) / cell))); c1[_k] = MAX(0, MIN(g[_k] - 1, (int)((_b - lo[_k]) / cell))); } } while (0)
		int64_t total = 0;
		for (int t = 0; t < ntri; t++) {
			int c0[3], c1[3];
			MC_CELLS(t, c0, c1);
			total += (int64_t)(c1[0] - c0[0] + 1) * (c1[1] - c0[1] + 1) * (c1[2] - c0[2] + 1);
			for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++)
				cnt[1 + x + g[0] * (y + (int64_t)g[1] * z)]++;
		}
		if (total > INT_MAX) { free(cnt); return -1; }   /* only a hostile mesh: faces spanning the whole grid */
		for (int64_t i = 0; i < ncell; i++) cnt[i + 1] += cnt[i];
		cell_tri = (int *)malloc((size_t)cnt[ncell] * sizeof(int));
		int *fill = (int *)calloc((size_t)ncell, sizeof(int));
		if (!cell_tri || !fill) { free(fill); free(cnt); free(cell_tri); return -1; }
		{
			for (int t = 0; t < ntri; t++) {
				int c0[3], c1[3];
				MC_CELLS(t, c0, c1);
				for (int z = c0[2]; z <= c1[2]; z++) for (int y = c0[1]; y <= c1[1]; y++) for (int x = c0[0]; x <= c1[0]; x++) {
					int64_t c = x + g[0] * (y + (int64_t)g[1] * z);
					cell_tri[cnt[c] + fill[c]++] = t;
				}
			}
			/* -mesh -a runs nii2mesh per label inside its own parallel region; this one then
			   runs serially on the encountering thread (nested parallelism is inactive) */
			#ifdef _OPENMP
			#pragma omp parallel for schedule(dynamic, 4096)
			#endif
			for (int64_t c = 0; c < ncell; c++) {
				for (int i = cnt[c]; i < cnt[c + 1]; i++) for (int j = i + 1; j < cnt[c + 1]; j++) {
					int a = cell_tri[i], b = cell_tri[j];
					int va[3] = { tris[a].x, tris[a].y, tris[a].z }, vb[3] = { tris[b].x, tris[b].y, tris[b].z }, adj = 0;
					for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) if (va[p] == vb[q]) adj = 1;
					if (adj) continue;
					if (mesh_tri_tri(&pts[va[0]].x, &pts[va[1]].x, &pts[va[2]].x, &pts[vb[0]].x, &pts[vb[1]].x, &pts[vb[2]].x)) {
						#ifdef _OPENMP
						#pragma omp atomic write
						#endif
						hit[a] = 1;
						#ifdef _OPENMP
						#pragma omp atomic write
						#endif
						hit[b] = 1;
					}
				}
			}
		}
		free(fill);
	}
	free(cnt); free(cell_tri);
	for (int t = 0; t < ntri; t++) nself += hit[t];
	return nself;
}

int mesh_report(vec3i *tris, vec3d *pts, int ntri, int npt, const char *label) {
	double t0 = clockMsec();
	size_t cap = 1;
	while (cap < (size_t)ntri * 6) cap <<= 1;
	mc_edge *h = (mc_edge *)calloc(cap, sizeof(mc_edge));
	int *uf = (int *)malloc((size_t)npt * sizeof(int));
	int *deg = (int *)calloc((size_t)npt + 1, sizeof(int));
	int *vt = (int *)malloc((size_t)ntri * 3 * sizeof(int));   /* vertex -> incident triangles, CSR */
	uint8_t *hit = (uint8_t *)calloc((size_t)ntri, 1);
	int nedge = 0, nbound = 0, nnonman = 0, nnmv = 0, ncomp = 0, nloops = 0, nself;
	if (!h || !uf || !deg || !vt || !hit) { free(h); free(uf); free(deg); free(vt); free(hit); return 1; }
	for (int i = 0; i < npt; i++) uf[i] = i;
	for (int t = 0; t < ntri; t++) {
		int v[3] = { tris[t].x, tris[t].y, tris[t].z };
		for (int k = 0; k < 3; k++) {
			mc_edge *e = mc_efind(h, cap, mc_ekey(v[k], v[(k + 1) % 3]) + 1);   /* +1: 0 is empty */
			if (!e->key) { e->key = mc_ekey(v[k], v[(k + 1) % 3]) + 1; e->t0 = t; e->t1 = -1; nedge++; }
			else if (e->t1 < 0) e->t1 = t;
			e->n++;
			mc_union(uf, v[k], v[(k + 1) % 3]);
			deg[v[k] + 1]++;
		}
	}
	for (size_t i = 0; i < cap; i++) if (h[i].key) { if (h[i].n == 1) nbound++; else if (h[i].n > 2) nnonman++; }
	for (int i = 0; i < npt; i++) if (deg[i + 1] && mc_find(uf, i) == i) ncomp++;
	/* boundary loops: union-find over boundary edges' endpoints, on a fresh forest.  Reported as
	   holes; two loops touching at one vertex count once (that vertex is also non-manifold). */
	for (int i = 0; i < npt; i++) uf[i] = i;
	for (size_t i = 0; i < cap; i++) if (h[i].key && h[i].n == 1) {
		uint64_t k = h[i].key - 1;
		mc_union(uf, (int)(k >> 32), (int)(uint32_t)k);
	}
	{	uint8_t *onb = (uint8_t *)calloc((size_t)npt, 1);
		if (onb) {
			for (size_t i = 0; i < cap; i++) if (h[i].key && h[i].n == 1) {
				uint64_t k = h[i].key - 1; onb[k >> 32] = 1; onb[(uint32_t)k] = 1; }
			for (int i = 0; i < npt; i++) if (onb[i] && mc_find(uf, i) == i) nloops++;
			free(onb);
		}
	}
	/* non-manifold vertices: the incident triangles must form ONE fan joined by edges through v */
	for (int i = 0; i < npt; i++) deg[i + 1] += deg[i];
	{	int *fill = (int *)calloc((size_t)npt, sizeof(int));
		if (fill) {
			for (int t = 0; t < ntri; t++) {
				int v[3] = { tris[t].x, tris[t].y, tris[t].z };
				for (int k = 0; k < 3; k++) vt[deg[v[k]] + fill[v[k]]++] = t;
			}
			for (int v = 0; v < npt; v++) {
				int n = deg[v + 1] - deg[v], *tv = vt + deg[v], comps = n;
				if (n < 2 || n > npt) continue;   /* uf is the scratch forest; n > npt only on duplicate-face garbage */
				for (int a = 0; a < n; a++) uf[a] = a;   /* uf reused as a tiny local forest */
				for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) {
					int ta = tv[a], tb = tv[b], va[3] = { tris[ta].x, tris[ta].y, tris[ta].z }, vb[3] = { tris[tb].x, tris[tb].y, tris[tb].z }, shared = 0;
					for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) if (va[p] != v && va[p] == vb[q]) shared = 1;
					if (shared && mc_find(uf, a) != mc_find(uf, b)) { uf[mc_find(uf, a)] = mc_find(uf, b); comps--; }
				}
				if (comps > 1) nnmv++;
			}
			free(fill);
		}
	}
	nself = mesh_self_intersections(tris, pts, ntri, npt, hit);
	{	int euler = npt - nedge + ntri;
		int genus2 = 2 * ncomp - nloops - euler;   /* 2G, from V - E + F = 2C - 2G - B */
		char self[16] = "n/a", genus[16] = "n/a";   /* self: the grid did not fit; genus: only meaningful on a manifold */
		if (nself >= 0) snprintf(self, sizeof self, "%d", nself);
		if (!nnonman && !nnmv) snprintf(genus, sizeof genus, "%g", genus2 / 2.0);
		printf("%s: V=%d E=%d F=%d components=%d euler=%d genus=%s boundary_edges=%d holes=%d nonmanifold_edges=%d nonmanifold_vertices=%d self_intersecting=%s (%ld ms)\n",
			label, npt, nedge, ntri, ncomp, euler, genus, nbound, nloops, nnonman, nnmv, self, timediff(t0, clockMsec()));
	}
	free(h); free(uf); free(deg); free(vt); free(hit);
	return 0;
}

static bool littleEndianPlatform() {
	uint32_t value = 1;
	return (*((char *)&value) == 1);
}

static void swap_4bytes(size_t n, void *ar) { // 4 bytes at a time
	size_t ii;
	unsigned char *cp0 = (unsigned char *)ar, *cp1, *cp2;
	unsigned char tval;
	for (ii = 0; ii < n; ii++) {
		cp1 = cp0;
		cp2 = cp0 + 3;
		tval = *cp1;
		*cp1 = *cp2;
		*cp2 = tval;
		cp1++;
		cp2--;
		tval = *cp1;
		*cp1 = *cp2;
		*cp2 = tval;
		cp0 += 4;
	}
	return;
}

#ifdef HAVE_ZLIB
#ifdef HAVE_JSON

enum TZipMethod { zmZlib,
				  zmGzip,
				  zmBase64,
				  zmLzip,
				  zmLzma,
				  zmLz4,
				  zmLz4hc };

static int zmat_run(const size_t inputsize, unsigned char *inputstr, size_t *outputsize, unsigned char **outputbuf, const int zipid, int *ret, const int iscompress) {
	z_stream zs;
	size_t buflen[2] = {0};
	*outputbuf = NULL;
	zs.zalloc = Z_NULL;
	zs.zfree = Z_NULL;
	zs.opaque = Z_NULL;
	if (inputsize == 0)
		return -1;
	if (iscompress) {
		/** perform compression or encoding   */
		if (zipid == zmBase64) {
			/** base64 encoding  */
			*outputbuf = base64_encode((const unsigned char *)inputstr, inputsize, outputsize);
		} else if (zipid == zmZlib) {
			/** zlib (.zip) or gzip (.gz) compression  */
			if (deflateInit(&zs, (iscompress > 0) ? Z_DEFAULT_COMPRESSION : (-iscompress)) != Z_OK)
				return -2;
			buflen[0] = deflateBound(&zs, inputsize);
			*outputbuf = (unsigned char *)malloc(buflen[0]);
			if (!*outputbuf) { deflateEnd(&zs); return -1; }
			zs.avail_in = inputsize;			 /* size of input, string + terminator*/
			zs.next_in = (Bytef *)inputstr;		 /* input char array*/
			zs.avail_out = buflen[0];			 /* size of output*/
			zs.next_out = (Bytef *)(*outputbuf); /*(Bytef *)(); // output char array*/
			*ret = deflate(&zs, Z_FINISH);
			*outputsize = zs.total_out;
			if (*ret != Z_STREAM_END && *ret != Z_OK)
				return -3;
			deflateEnd(&zs);
		} else {
			return -7;
		}
	} else {
		/** perform decompression or decoding */
		if (zipid == zmBase64) {
			/** base64 decoding  */
			*outputbuf = base64_decode((const unsigned char *)inputstr, inputsize, outputsize);
		} else if (zipid == zmZlib) {
			/** zlib (.zip) or gzip (.gz) decompression */
			int count = 1;
			if (zipid == zmZlib)
				if (inflateInit(&zs) != Z_OK)
					return -2;
			buflen[0] = inputsize * 20;
			*outputbuf = (unsigned char *)malloc(buflen[0]);
			if (!*outputbuf) { inflateEnd(&zs); return -1; }
			zs.avail_in = inputsize;			 /* size of input, string + terminator*/
			zs.next_in = inputstr;				 /* input char array*/
			zs.avail_out = buflen[0];			 /* size of output*/
			zs.next_out = (Bytef *)(*outputbuf); /*(Bytef *)(); // output char array*/
			while ((*ret = inflate(&zs, Z_SYNC_FLUSH)) != Z_STREAM_END && count <= 10) {
				unsigned char *tmp = (unsigned char *)realloc(*outputbuf, (buflen[0] << count));
				if (!tmp) { free(*outputbuf); *outputbuf = NULL; inflateEnd(&zs); return -1; }
				*outputbuf = tmp;
				zs.next_out = (Bytef *)(*outputbuf + (buflen[0] << (count - 1)));
				zs.avail_out = (buflen[0] << (count - 1)); /* size of output*/
				count++;
			}
			*outputsize = zs.total_out;

			if (*ret != Z_STREAM_END && *ret != Z_OK)
				return -3;
			inflateEnd(&zs);
		} else {
			return -7;
		}
	}
	return 0;
}

static int save_jmsh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
	FILE *fp;
	cJSON *root = NULL, *hdr = NULL, *node = NULL, *face = NULL;
	char *jsonstr = NULL;
	int dim[2] = {0, 3}, len[2] = {1, 0};
	size_t compressedbytes, totalbytes;
	unsigned char *compressed = NULL, *buf = NULL;
	int ret = 0, status = 0;
	root = cJSON_CreateObject();
	cJSON_AddItemToObject(root, "_DataInfo_", hdr = cJSON_CreateObject());
	cJSON_AddStringToObject(hdr, "JMeshVersion", "0.5");
	cJSON_AddStringToObject(hdr, "Comment", "Created by nii2mesh");
	cJSON_AddItemToObject(root, "MeshVertex3", node = cJSON_CreateObject());
	cJSON_AddStringToObject(node, "_ArrayType_", "double");
	dim[0] = npt;
	cJSON_AddItemToObject(node, "_ArraySize_", cJSON_CreateIntArray(dim, 2));
	cJSON_AddStringToObject(node, "_ArrayZipType_", "zlib");
	len[1] = dim[0] * dim[1];
	cJSON_AddItemToObject(node, "_ArrayZipSize_", cJSON_CreateIntArray(len, 2));
	totalbytes = dim[0] * dim[1] * sizeof(pts[0].x);
	unsigned int *val = (unsigned int *)malloc(totalbytes);
	if (!val) { fclose(fp); return EXIT_FAILURE; }
	memcpy(val, &(tris[0].x), totalbytes);
	for (int i = 0; i < len[1]; i++)
		val[i]++;
	ret = zmat_run(totalbytes, (unsigned char *)val, &compressedbytes, (unsigned char **)&compressed, zmZlib, &status, 1);
	free(val);
	if (!ret) {
		ret = zmat_run(compressedbytes, compressed, &totalbytes, (unsigned char **)&buf, zmBase64, &status, 1);
		cJSON_AddStringToObject(node, "_ArrayZipData_", (char *)buf);
	}
	if (compressed) {
		free(compressed);
		compressed = NULL;
	}
	if (buf) {
		free(buf);
		buf = NULL;
	}
	cJSON_AddItemToObject(root, "MeshTri3", face = cJSON_CreateObject());
	cJSON_AddStringToObject(face, "_ArrayType_", "uint32");
	dim[0] = ntri;
	cJSON_AddItemToObject(face, "_ArraySize_", cJSON_CreateIntArray(dim, 2));
	cJSON_AddStringToObject(face, "_ArrayZipType_", "zlib");
	len[1] = dim[0] * dim[1];
	cJSON_AddItemToObject(face, "_ArrayZipSize_", cJSON_CreateIntArray(len, 2));
	totalbytes = dim[0] * dim[1] * sizeof(tris[0].x);
	ret = zmat_run(totalbytes, (unsigned char *)&(tris[0].x), &compressedbytes, (unsigned char **)&compressed, zmZlib, &status, 1);
	if (!ret) {
		ret = zmat_run(compressedbytes, compressed, &totalbytes, (unsigned char **)&buf, zmBase64, &status, 1);
		cJSON_AddStringToObject(face, "_ArrayZipData_", (char *)buf);
	}
	if (compressed)
		free(compressed);
	if (buf)
		free(buf);
	jsonstr = cJSON_Print(root);
	if (jsonstr == NULL)
		return EXIT_FAILURE;
	fp = fopen(fnm, "wt");
	if (fp == NULL)
		return EXIT_FAILURE;
	fprintf(fp, "%s\n", jsonstr);
	fclose(fp);
	if (jsonstr)
		free(jsonstr);
	if (root)
		cJSON_Delete(root);
	return EXIT_SUCCESS;
}
#endif // HAVE_JSON
#endif // HAVE_ZLIB

typedef struct {
	float x, y, z;
} vec3s; // single precision (float32)

static vec3s vec3d2vec4s(vec3d v) {
	return (vec3s){.x = (float)v.x, .y = (float)v.y, .z = (float)v.z};
} // convert float64 to float32

int save_mz3(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz, float *perVertexScalar) {
	// https://github.com/neurolabusc/surf-ice/tree/master/mz3
#ifdef _MSC_VER
#pragma pack(2)
	struct mz3hdr {
		uint16_t SIGNATURE, ATTR;
		uint32_t NFACE, NVERT, NSKIP;
	};
#pragma pack()
#else
	struct __attribute__((__packed__)) mz3hdr {
		uint16_t SIGNATURE, ATTR;
		uint32_t NFACE, NVERT, NSKIP;
	};
#endif
	bool hasPerVertexScalar = (perVertexScalar != NULL);
	struct mz3hdr h;
	h.SIGNATURE = 0x5A4D;
	h.ATTR = 3; // isFACE +1 isVERT +2
	if (hasPerVertexScalar)
		h.ATTR += 8; // isSCALAR +8
	h.NFACE = ntri;
	h.NVERT = npt;
	h.NSKIP = 0;
	if (!littleEndianPlatform())
		swap_4bytes(3, &h.NFACE);
	FILE *fp;
#ifdef HAVE_ZLIB
	gzFile fgz;
	if (isGz) {
		fgz = gzopen(fnm, "w");
		if (!fgz)
			return EXIT_FAILURE;
		gzwrite(fgz, &h, sizeof(struct mz3hdr));
	} else
#endif
	{
		fp = fopen(fnm, "wb");
		if (fp == NULL)
			return EXIT_FAILURE;
		fwrite(&h, sizeof(struct mz3hdr), 1, fp);
	}
	if (!littleEndianPlatform()) {
		vec3i *trisSwap = (vec3i *)malloc(ntri * sizeof(vec3i));
		if (!trisSwap) goto mz3_fail;
		for (int i = 0; i < ntri; i++)
			trisSwap[i] = tris[i];
		swap_4bytes(ntri * 3, trisSwap);
#ifdef HAVE_ZLIB
		if (isGz)
			gzwrite(fgz, trisSwap, ntri * sizeof(vec3i));
		else
#endif
			fwrite(trisSwap, ntri * sizeof(vec3i), 1, fp);
		free(trisSwap);
	} else {
#ifdef HAVE_ZLIB
		if (isGz)
			gzwrite(fgz, tris, ntri * sizeof(vec3i));
		else
#endif
			fwrite(tris, ntri * sizeof(vec3i), 1, fp);
	}
	vec3s *pts32 = (vec3s *)malloc(npt * sizeof(vec3s));
	if (!pts32) goto mz3_fail;
	for (int i = 0; i < npt; i++) // double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	if (!littleEndianPlatform())
		swap_4bytes(npt * 3, pts32);
#ifdef HAVE_ZLIB
	if (isGz) {
		gzwrite(fgz, pts32, npt * sizeof(vec3s));
		gzclose(fgz);
	} else
#endif
	{
		fwrite(pts32, npt * sizeof(vec3s), 1, fp);
		if (hasPerVertexScalar) {
			fwrite(perVertexScalar, npt * sizeof(float), 1, fp);
		}
		fclose(fp);
	}
	free(pts32);
	return EXIT_SUCCESS;
mz3_fail:
#ifdef HAVE_ZLIB
	if (isGz) gzclose(fgz); else
#endif
	fclose(fp);
	return EXIT_FAILURE;
}

#ifdef HAVE_FORMATS

static int save_freesurfer(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
	// FreeSurfer Triangle Surface Binary Format http://www.grahamwideman.com/gw/brain/fs/surfacefileformats.htm
	uint8_t magic[3] = {0xFF, 0xFF, 0xFE};
	FILE *fp = fopen(fnm, "wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fwrite(magic, 3, 1, fp);
	time_t t = time(NULL);
	char s[128] = "";
	struct tm *tm = localtime(&t);
	strftime(s, sizeof(s), "created by niimath on %c\n\n", tm);
	fwrite(s, strlen(s), 1, fp);
	int32_t VertexCount = npt;
	int32_t FaceCount = ntri;
	if (littleEndianPlatform()) {
		swap_4bytes(1, &VertexCount);
		swap_4bytes(1, &FaceCount);
	}
	fwrite(&VertexCount, sizeof(int32_t), 1, fp);
	fwrite(&FaceCount, sizeof(int32_t), 1, fp);
	vec3s *pts32 = (vec3s *)malloc(npt * sizeof(vec3s));
	if (!pts32) { fclose(fp); return EXIT_FAILURE; }
	for (int i = 0; i < npt; i++) // double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	if (littleEndianPlatform())
		swap_4bytes(npt * 3, pts32);
	fwrite(pts32, npt * sizeof(vec3s), 1, fp);
	free(pts32);
	if (littleEndianPlatform()) {
		vec3i *trisSwap = (vec3i *)malloc(ntri * sizeof(vec3i));
		if (!trisSwap) { fclose(fp); return EXIT_FAILURE; }
		for (int i = 0; i < ntri; i++)
			trisSwap[i] = tris[i];
		swap_4bytes(ntri * 3, trisSwap);
		fwrite(trisSwap, ntri * sizeof(vec3i), 1, fp);
		free(trisSwap);
	} else
		fwrite(tris, ntri * sizeof(vec3i), 1, fp);
	fclose(fp);
	return EXIT_SUCCESS;
}

static int save_json(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
	FILE *fp = fopen(fnm, "w");
	if (fp == NULL)
		return EXIT_FAILURE;
	fprintf(fp, "{\n");
	fprintf(fp, "\t\"_DataInfo_\":{\n\t\t\"JMeshVersion\":\"0.5\",\n\t\t\"Comment\":\"Created by nii2mesh\"\n\t},\n");
	fprintf(fp, "\t\"MeshVertex3\":[\n");
	for (int i = 0; i < npt; i++)
		fprintf(fp, "[%g,\t%g,\t%g],\n", pts[i].x, pts[i].y, pts[i].z);
	fprintf(fp, "\t],\n\t\"MeshTri3\":[\n");
	for (int i = 0; i < ntri; i++)
		fprintf(fp, "[%d,\t%d,\t%d],\n", tris[i].x + 1, tris[i].y + 1, tris[i].z + 1);
	fprintf(fp, "\t]\n}\n");
	fclose(fp);
	return EXIT_SUCCESS;
}

static int save_off(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
	FILE *fp = fopen(fnm, "w");
	if (fp == NULL)
		return EXIT_FAILURE;
	fprintf(fp, "OFF\n%d\t%d\t0\n", npt, ntri);
	for (int i = 0; i < npt; i++)
		fprintf(fp, "%g %g %g\n", pts[i].x, pts[i].y, pts[i].z);
	for (int i = 0; i < ntri; i++)
		fprintf(fp, "%d %d %d\n", tris[i].x + 1, tris[i].y + 1, tris[i].z + 1);
	fclose(fp);
	return EXIT_SUCCESS;
}

static int save_obj(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
	FILE *fp = fopen(fnm, "w");
	if (fp == NULL)
		return EXIT_FAILURE;
	for (int i = 0; i < npt; i++)
		fprintf(fp, "v %g %g %g\n", pts[i].x, pts[i].y, pts[i].z);
	for (int i = 0; i < ntri; i++)
		fprintf(fp, "f %d %d %d\n", tris[i].x + 1, tris[i].y + 1, tris[i].z + 1);
	fclose(fp);
	return EXIT_SUCCESS;
}

static int save_stl(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
// binary STL http://paulbourke.net/dataformats/stl/
// n.b. like other tools, ignores formal restriction that all adjacent facets must share two common vertices.
// n.b. does not write normal
#ifdef _MSC_VER
#pragma pack(2)
	typedef struct {
		vec3s norm, pts[3];
		uint16_t spacer;
	} tfacet;
#pragma pack()
#else
	typedef struct __attribute__((__packed__)) {
		vec3s norm, pts[3];
		uint16_t spacer;
	} tfacet;
#endif
	FILE *fp = fopen(fnm, "wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	uint8_t hdr[80] = {0};
	fwrite(hdr, 80, 1, fp);
	int32_t nf = ntri;
	fwrite(&nf, sizeof(int32_t), 1, fp);
	tfacet *facets = (tfacet *)malloc(ntri * sizeof(tfacet));
	if (!facets) { fclose(fp); return EXIT_FAILURE; }
	vec3s n0 = (vec3s){.x = 0.0, .y = 0.0, .z = 0.0};
	for (int i = 0; i < ntri; i++) { // double->single precision
		facets[i].norm = n0;
		facets[i].pts[0] = vec3d2vec4s(pts[tris[i].x]);
		facets[i].pts[1] = vec3d2vec4s(pts[tris[i].y]);
		facets[i].pts[2] = vec3d2vec4s(pts[tris[i].z]);
		facets[i].spacer = 0;
	}
	fwrite(facets, ntri * sizeof(tfacet), 1, fp);
	free(facets);
	fclose(fp);
	return EXIT_SUCCESS;
}

static int save_ply(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
#ifdef _MSC_VER
#pragma pack(1)
	typedef struct {
		uint8_t n;
		int32_t x, y, z;
	} vec1b3i;
#pragma pack()
#else
	typedef struct __attribute__((__packed__)) {
		uint8_t n;
		int32_t x, y, z;
	} vec1b3i;
#endif
	FILE *fp = fopen(fnm, "wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fputs("ply\n", fp);
	if (littleEndianPlatform())
		fputs("format binary_little_endian 1.0\n", fp);
	else
		fputs("format binary_big_endian 1.0\n", fp);
	fputs("comment niimath\n", fp);
	char vpts[80];
	sprintf(vpts, "element vertex %d\n", npt);
	fwrite(vpts, strlen(vpts), 1, fp);
	fputs("property float x\n", fp);
	fputs("property float y\n", fp);
	fputs("property float z\n", fp);
	char vfc[80];
	sprintf(vfc, "element face %d\n", ntri);
	fwrite(vfc, strlen(vfc), 1, fp);
	fputs("property list uchar int vertex_indices\n", fp);
	fputs("end_header\n", fp);
	vec3s *pts32 = (vec3s *)malloc(npt * sizeof(vec3s));
	if (!pts32) { fclose(fp); return EXIT_FAILURE; }
	for (int i = 0; i < npt; i++) { // double->single precision
		pts32[i].x = pts[i].x;
		pts32[i].y = pts[i].y;
		pts32[i].z = pts[i].z;
	}
	fwrite(pts32, npt * sizeof(vec3s), 1, fp);
	free(pts32);
	vec1b3i *tris4 = (vec1b3i *)malloc(ntri * sizeof(vec1b3i));
	if (!tris4) { fclose(fp); return EXIT_FAILURE; }
	for (int i = 0; i < ntri; i++) { // double->single precision
		tris4[i].n = 3;
		tris4[i].x = tris[i].x;
		tris4[i].y = tris[i].y;
		tris4[i].z = tris[i].z;
	}
	fwrite(tris4, ntri * sizeof(vec1b3i), 1, fp);
	free(tris4);
	fclose(fp);
	return EXIT_SUCCESS;
}

static int save_gii(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz) {
	// https://www.nitrc.org/projects/gifti/
	// https://stackoverflow.com/questions/342409/how-do-i-base64-encode-decode-in-c
	FILE *fp = fopen(fnm, "wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fputs("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n", fp);
	fputs("<!DOCTYPE GIFTI SYSTEM \"http://www.nitrc.org/frs/download.php/115/gifti.dtd\">\n", fp);
	fputs("<GIFTI Version=\"1.0\"  NumberOfDataArrays=\"2\">\n", fp);
	fputs("   <MetaData>\n", fp);
	fputs("      <MD>\n", fp);
	fputs("          <Name><![CDATA[nii2mesh-version]]></Name>\n", fp);
	fputs("          <Value><![CDATA[nii2mesh, 1 Jan 2022]]></Value>\n", fp);
	fputs("      </MD>\n", fp);
	fputs("   </MetaData>\n", fp);
	fputs("   <LabelTable/>\n", fp);
	fputs("   <DataArray  ArrayIndexingOrder=\"RowMajorOrder\"\n", fp);
	fputs("               DataType=\"NIFTI_TYPE_INT32\"\n", fp);
	char fc[80];
	sprintf(fc, "               Dim0=\"%d\"\n", ntri);
	fwrite(fc, strlen(fc), 1, fp);
	fputs("               Dim1=\"3\"\n", fp);
	fputs("               Dimensionality=\"2\"\n", fp);
#ifdef HAVE_ZLIB
	if (isGz)
		fputs("               Encoding=\"GZipBase64Binary\"\n", fp);
	else
#endif
		fputs("               Encoding=\"Base64Binary\"\n", fp);
	if (littleEndianPlatform())
		fputs("               Endian=\"LittleEndian\"\n", fp);
	else
		fputs("               Endian=\"BigEndian\"\n", fp);
	fputs("               ExternalFileName=\"\"\n", fp);
	fputs("               ExternalFileOffset=\"\"\n", fp);
	fputs("               Intent=\"NIFTI_INTENT_TRIANGLE\">\n", fp);
	fputs("      <MetaData>\n", fp);
	fputs("      </MetaData>\n", fp);
	fputs("      <Data>", fp);
	size_t out_len;
	unsigned char *fcs;
#ifdef HAVE_ZLIB
	if (isGz) {
		unsigned long srcLen = ntri * sizeof(vec3i);
		uLongf destLen = compressBound(srcLen);
		unsigned char *ostream = (unsigned char *)malloc(destLen);
		if (!ostream) { fclose(fp); return EXIT_FAILURE; }
		int res = compress(ostream, &destLen, (const unsigned char *)tris, srcLen);
		if (res != Z_OK)
			printf("Compression error\n");
		fcs = base64_encode(ostream, destLen, &out_len);
		free(ostream);
	} else
#endif
		fcs = base64_encode((const unsigned char *)tris, ntri * sizeof(vec3i), &out_len);
	fwrite(fcs, out_len, 1, fp);
	free(fcs);
	fputs("</Data>\n", fp);
	fputs("   </DataArray>\n", fp);
	fputs("   <DataArray  ArrayIndexingOrder=\"RowMajorOrder\"\n", fp);
	fputs("               DataType=\"NIFTI_TYPE_FLOAT32\"\n", fp);
	char vpts[80];
	sprintf(vpts, "               Dim0=\"%d\"\n", npt);
	fwrite(vpts, strlen(vpts), 1, fp);
	fputs("               Dim1=\"3\"\n", fp);
	fputs("               Dimensionality=\"2\"\n", fp);
#ifdef HAVE_ZLIB
	if (isGz)
		fputs("               Encoding=\"GZipBase64Binary\"\n", fp);
	else
#endif
		fputs("               Encoding=\"Base64Binary\"\n", fp);
	if (littleEndianPlatform())
		fputs("               Endian=\"LittleEndian\"\n", fp);
	else
		fputs("               Endian=\"BigEndian\"\n", fp);
	fputs("               ExternalFileName=\"\"\n", fp);
	fputs("               ExternalFileOffset=\"\"\n", fp);
	fputs("               Intent=\"NIFTI_INTENT_POINTSET\">\n", fp);
	fputs("      <MetaData>\n", fp);
	fputs("      </MetaData>\n", fp);
	fputs("      <CoordinateSystemTransformMatrix>\n", fp);
	fputs("         <DataSpace><![CDATA[NIFTI_XFORM_UNKNOWN]]></DataSpace>\n", fp);
	fputs("         <TransformedSpace><![CDATA[NIFTI_XFORM_UNKNOWN]]></TransformedSpace>\n", fp);
	fputs("         <MatrixData>1.000000 0.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 0.000000 1.000000 </MatrixData>\n", fp);
	fputs("      </CoordinateSystemTransformMatrix>\n", fp);
	fputs("      <Data>", fp);
	vec3s *pts32 = (vec3s *)malloc(npt * sizeof(vec3s));
	if (!pts32) { fclose(fp); return EXIT_FAILURE; }
	for (int i = 0; i < npt; i++) // double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	unsigned char *vts;
#ifdef HAVE_ZLIB
	if (isGz) {
		unsigned long srcLen = npt * sizeof(vec3s);
		uLongf destLen = compressBound(srcLen);
		unsigned char *ostream = (unsigned char *)malloc(destLen);
		if (!ostream) { free(pts32); fclose(fp); return EXIT_FAILURE; }
		int res = compress(ostream, &destLen, (const unsigned char *)pts32, srcLen);
		if (res != Z_OK)
			printf("Compression error\n");
		vts = base64_encode(ostream, destLen, &out_len);
		free(ostream);
	} else
#endif
		vts = base64_encode((const unsigned char *)pts32, npt * sizeof(vec3s), &out_len);
	free(pts32);
	fwrite(vts, out_len, 1, fp);
	free(vts);
	fputs("</Data>\n", fp);
	fputs("   </DataArray>\n", fp);
	fputs("</GIFTI>\n", fp);
	fclose(fp);
	return EXIT_SUCCESS;
}

static int save_vtk(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt) {
	typedef struct {
		uint32_t n, x, y, z;
	} vec4i;
	FILE *fp = fopen(fnm, "wb");
	if (fp == NULL)
		return EXIT_FAILURE;
	fputs("# vtk DataFile Version 3.0\n", fp);
	fputs("this file was written using niimath\n", fp);
	fputs("BINARY\n", fp);
	fputs("DATASET POLYDATA\n", fp);
	char vpts[80];
	sprintf(vpts, "POINTS %d float\n", npt);
	fwrite(vpts, strlen(vpts), 1, fp);
	vec3s *pts32 = (vec3s *)malloc(npt * sizeof(vec3s));
	if (!pts32) { fclose(fp); return EXIT_FAILURE; }
	for (int i = 0; i < npt; i++) // double->single precision
		pts32[i] = vec3d2vec4s(pts[i]);
	if (littleEndianPlatform())
		swap_4bytes(3 * npt, pts32);
	fwrite(pts32, npt * sizeof(vec3s), 1, fp);
	free(pts32);
	char vfac[80];
	sprintf(vfac, "POLYGONS %d %d\n", ntri, ntri * 4);
	fwrite(vfac, strlen(vfac), 1, fp);
	vec4i *tris4 = (vec4i *)malloc(ntri * sizeof(vec4i));
	if (!tris4) { fclose(fp); return EXIT_FAILURE; }
	for (int i = 0; i < ntri; i++) { // double->single precision
		tris4[i].n = 3;
		tris4[i].x = tris[i].x;
		tris4[i].y = tris[i].y;
		tris4[i].z = tris[i].z;
	}
	if (littleEndianPlatform())
		swap_4bytes(4 * ntri, tris4);
	fwrite(tris4, ntri * sizeof(vec4i), 1, fp);
	free(tris4);
	fclose(fp);
	return EXIT_SUCCESS;
}

#endif // HAVE_FORMATS

void strip_ext(char *fname) {
	char *end = fname + strlen(fname);
	while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
		--end;
	}
	if ((end > fname && *end == '.') &&
		(*(end - 1) != '\\' && *(end - 1) != '/')) {
		*end = '\0';
	}
}

int save_mesh(const char *fnm, vec3i *tris, vec3d *pts, int ntri, int npt, bool isGz) {
	char basenm[768], ext[768] = "";
	size_t fnmlen = strlen(fnm);
	if (fnmlen >= sizeof(basenm) - 5) {
		fprintf(stderr, "** filename too long (max %d chars): %s\n", (int)(sizeof(basenm) - 6), fnm);
		return EXIT_FAILURE;
	}
	strncpy(basenm, fnm, sizeof(basenm) - 1);
	basenm[sizeof(basenm) - 1] = '\0';
	strip_ext(basenm); // ~/file.nii -> ~/file
	if (fnmlen > strlen(basenm))
		strncpy(ext, fnm + strlen(basenm), sizeof(ext) - 1);
	ext[sizeof(ext) - 1] = '\0';
	if (strstr(ext, ".mz3"))
		return save_mz3(fnm, tris, pts, ntri, npt, isGz, NULL);
#ifdef HAVE_FORMATS
	else if (strstr(ext, ".gii"))
		return save_gii(fnm, tris, pts, ntri, npt, isGz);
	else if ((strstr(ext, ".inflated")) || (strstr(ext, ".pial")))
		return save_freesurfer(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".json"))
		return save_json(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".off"))
		return save_off(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".obj"))
		return save_obj(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".ply"))
		return save_ply(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".stl"))
		return save_stl(fnm, tris, pts, ntri, npt);
	else if (strstr(ext, ".vtk"))
		return save_vtk(fnm, tris, pts, ntri, npt);
#endif // HAVE_FORMATS
#ifdef HAVE_ZLIB
#ifdef HAVE_JSON
	else if (strstr(ext, ".jmsh"))
		return save_jmsh(fnm, tris, pts, ntri, npt);
#endif // HAVE_JSON
#endif // HAVE_ZLIB
	snprintf(basenm, sizeof(basenm), "%s.mz3", fnm);
	return save_mz3(basenm, tris, pts, ntri, npt, isGz, NULL);
}

static double sform(vec3d p, float srow[4]) {
	return (p.x * srow[0]) + (p.y * srow[1]) + (p.z * srow[2]) + srow[3];
}

void apply_sform(vec3i *t, vec3d *p, int nt, int np, float srow_x[4], float srow_y[4], float srow_z[4]) {
	for (int i = 0; i < np; i++) {
		vec3d v = p[i];
		p[i].x = sform(v, srow_x);
		p[i].y = sform(v, srow_y);
		p[i].z = sform(v, srow_z);
	}
	// detect determinant
	vec3d p0;
	p0.x = srow_x[0] + srow_x[1] + srow_x[2];
	p0.y = srow_y[0] + srow_y[1] + srow_y[2];
	p0.z = srow_z[0] + srow_z[1] + srow_z[2];
	float det = p0.x * p0.y * p0.z;
	if (det >= 0.0)
		return; // positive volume
	// negative volume: we need to reverse the triangle winding
	for (int i = 0; i < nt; i++) {
		vec3i f = t[i];
		t[i].x = f.y;
		t[i].y = f.x;
	}
}

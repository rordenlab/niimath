/* niimath -n4: C port of ITK/ANTs N4 bias-field correction.
 *
 *  Derived from the Insight Toolkit (Copyright 1999-2019 Insight Software Consortium,
 *  Copyright 2020-present NumFOCUS) and ANTs (Copyright the ANTs contributors).
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file
 *  except in compliance with the License. You may obtain a copy of the License at
 *      https://www.apache.org/licenses/LICENSE-2.0.txt  (also LICENSE in this directory)
 *  Unless required by applicable law or agreed to in writing, software distributed under the
 *  License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 *  either express or implied. See the License for the specific language governing
 *  permissions and limitations under the License. See NOTICE in this directory.
 */
// B-spline scattered-data fitting, lattice refinement and grid evaluation, ported from ITK's
// BSplineScatteredDataPointSetToImageFilter / BSplineControlPointImageFilter (cubic, open ends).
#include "n4.h"
#include "n4_priv.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define N4_ORDER N4_SPLINE_ORDER

// itk::BSplineKernelFunction<3> (double). The volatile 6 keeps -freciprocal-math from turning the
// division into a multiply by an inexact 1/6.
static double n4_cubic(double u) {
	static volatile double six = 6.0;
	double a = fabs(u);
	if (a < 1.0) {
		double s = a * a;
		return (4.0 - 6.0 * s + 3.0 * s * a) / six;
	}
	if (a < 2.0) {
		double s = a * a;
		return (8.0 - 12.0 * a + 6.0 * s - s * a) / six;
	}
	return 0.0;
}

int n4_bspline_fit(const float *resid, const uint8_t *incl, const float *conf, const int m[3],
                   const double sp[3], const int nc[3], float *phi) {
	const size_t plane = (size_t)nc[0] * nc[1], nlat = plane * nc[2];
	// Every voxel of a shrunk z slice shares one z span, so a slice touches only 4 control planes.
	// Each slice accumulates its own 4-plane partial and the partials are added in slice order:
	// the result does not depend on the thread count (plan_n4.md M7), though it rounds
	// differently from ITK's point order.
	float *omega = (float *)calloc(nlat, sizeof(float)), *delta = (float *)calloc(nlat, sizeof(float));
	float *pom = (float *)calloc(4 * plane * m[2], sizeof(float)), *pde = (float *)calloc(4 * plane * m[2], sizeof(float));
	unsigned *zspan = (unsigned *)calloc(m[2], sizeof(unsigned));
	int bad = !omega || !delta || !pom || !pde || !zspan;
	float r[3], eps[3];
	for (int d = 0; d < 3; d++) { // ITK: RealType float, spacing double
		volatile double rd = (double)(float)(nc[d] - N4_ORDER) / ((double)(float)(m[d] - 1) * sp[d]);
		r[d] = (float)rd;
		eps[d] = (float)(r[d] * sp[d] * (double)1e-3f);
	}
	if (!bad) {
#pragma omp parallel for schedule(dynamic) reduction(| : bad)
		for (int z = 0; z < m[2]; z++) {
			float *om = pom + 4 * plane * z, *de = pde + 4 * plane * z;
			size_t k = (size_t)z * m[0] * m[1];
			for (int y = 0; y < m[1]; y++)
				for (int x = 0; x < m[0]; x++, k++) {
					if (!incl[k]) continue;
					// ITK's point is float(origin + index*spacing) in a float point set; origin taken as 0
					// ponytail: origin omitted (ITK's float point rounding depends on it, ~1e-7 relative)
					const int ix[3] = {x, y, z};
					float p[3];
					unsigned span[3];
					double kern[3][4];
					int outside = 0;
					for (int d = 0; d < 3; d++) {
						float pt = (float)(ix[d] * sp[d]);
						const float nspans = (float)(nc[d] - N4_ORDER);
						p[d] = (float)((double)pt * r[d]);
						if (fabsf(p[d] - nspans) <= eps[d]) p[d] = nspans - eps[d];
						if (p[d] < 0.0f && fabsf(p[d]) <= eps[d]) p[d] = 0.0f;
						if (p[d] < 0.0f || p[d] >= nspans) outside = 1; // ITK throws: outside the domain
						span[d] = outside ? 0 : (unsigned)p[d];
						for (int j = 0; j < 4; j++) {
							float u = (float)((float)(p[d] - (float)span[d] - (float)j) + 0.5 * (float)(N4_ORDER - 1));
							kern[d][j] = n4_cubic(u);
						}
					}
					if (outside) {
						bad = 1;
						continue;
					}
					zspan[z] = span[2];
					// weights of the 4x4x4 neighbourhood, x fastest, as ITK's iterator visits it
					float B[64], w2Sum = 0.0f;
					for (int j = 0, n = 0; j < 4; j++)
						for (int i = 0; i < 4; i++)
							for (int h = 0; h < 4; h++, n++) {
								float b = 1.0f;
								b = (float)(b * kern[0][h]);
								b = (float)(b * kern[1][i]);
								b = (float)(b * kern[2][j]);
								B[n] = b;
								w2Sum += b * b;
							}
					const float wc = conf ? conf[k] : 1.0f;
					for (int j = 0, n = 0; j < 4; j++)
						for (int i = 0; i < 4; i++)
							for (int h = 0; h < 4; h++, n++) {
								size_t l = ((size_t)j * nc[1] + (span[1] + i)) * nc[0] + (span[0] + h);
								const float t = B[n];
								om[l] = om[l] + wc * t * t;
								de[l] = de[l] + resid[k] * n4_fdiv(t * t * t * wc, w2Sum);
							}
				}
		}
	}
	if (!bad)
		for (int z = 0; z < m[2]; z++) { // slice order, so every control point sums its terms in z order
			const float *om = pom + 4 * plane * z, *de = pde + 4 * plane * z;
			for (size_t l = 0; l < 4 * plane && zspan[z] * plane + l < nlat; l++) {
				omega[zspan[z] * plane + l] += om[l];
				delta[zspan[z] * plane + l] += de[l];
			}
		}
	for (size_t l = 0; l < nlat && !bad; l++) {
		float v = 0.0f;
		if (fabsf(omega[l]) > 0.1f * 1.1920929e-7f) { // !itk::Math::AlmostEquals(omega, 0)
			v = n4_fdiv(delta[l], omega[l]);
			if (!isfinite(v)) v = 0.0f;
		}
		phi[l] = v;
	}
	free(omega); free(delta); free(pom); free(pde); free(zspan);
	return bad;
}

// Cubic subdivision in ITK's indexing (lattice origin one spacing before the domain): fine point
// 2p is the midpoint rule (c[p] + c[p+1]) / 2, fine point 2p+1 the vertex rule
// (c[p] + 6 c[p+1] + c[p+2]) / 8; ITK solves these by SVD, the masks are exact.
void n4_bspline_refine(const float *phi, const int nc[3], float *out) {
	static const float M[2][4] = {{0.5f, 0.5f, 0.0f, 0.0f}, {0.125f, 0.75f, 0.125f, 0.0f}};
	int nf[3];
	for (int d = 0; d < 3; d++) nf[d] = 2 * nc[d] - N4_ORDER;
	for (int z = 0; z < nf[2]; z++)
		for (int y = 0; y < nf[1]; y++)
			for (int x = 0; x < nf[0]; x++) {
				float sum = 0.0f;
				for (int k = 0; k < 4; k++) {
					int cz = z / 2 + k;
					if (cz >= nc[2] || M[z & 1][k] == 0.0f) continue;
					for (int j = 0; j < 4; j++) {
						int cy = y / 2 + j;
						if (cy >= nc[1] || M[y & 1][j] == 0.0f) continue;
						for (int i = 0; i < 4; i++) {
							int cx = x / 2 + i;
							if (cx >= nc[0] || M[x & 1][i] == 0.0f) continue;
							sum += M[x & 1][i] * M[y & 1][j] * M[z & 1][k] *
							       phi[((size_t)cz * nc[1] + cy) * nc[0] + cx];
						}
					}
				}
				out[((size_t)z * nf[1] + y) * nf[0] + x] = sum;
			}
}

// Parametric coordinate of grid index i on an axis of m samples over `spans` spans, with ITK's
// end clamps; the base-spline weights of the four control points it touches go in w.
static int n4_basis(int i, int m, int spans, float eps, unsigned *first, float w[4]) {
	float u = n4_fdiv((float)spans * (float)i, (float)(m - 1)); // n4_core guarantees m >= 2
	if (fabsf(u - (float)spans) <= eps) u = (float)spans - eps;
	if (u < 0.0f && fabsf(u) <= eps) u = 0.0f;
	if (u < 0.0f || u >= (float)spans) return 1;
	*first = (unsigned)u;
	for (int j = 0; j < 4; j++) w[j] = (float)n4_cubic((float)(u - (float)(*first + j) + 0.5 * (N4_ORDER - 1)));
	return 0;
}

int n4_bspline_eval(const float *phi, const int nc[3], const int m[3], const double sp[3], float *out) {
	const int nx = nc[0], ny = nc[1];
	// per-axis basis tables: first control point and four weights for every grid index
	unsigned *first[3] = {NULL, NULL, NULL};
	float *w[3] = {NULL, NULL, NULL};
	int rc = 0;
	for (int d = 0; d < 3 && !rc; d++) {
		first[d] = (unsigned *)malloc(m[d] * sizeof(unsigned));
		w[d] = (float *)malloc(4 * (size_t)m[d] * sizeof(float));
		if (!first[d] || !w[d]) { rc = 1; break; }
		float r = (float)((double)(float)(nc[d] - N4_ORDER) / ((double)(float)(m[d] - 1) * sp[d]));
		float eps = (float)(r * sp[d] * (double)1e-3f);
		for (int i = 0; i < m[d] && !rc; i++) rc = n4_basis(i, m[d], nc[d] - N4_ORDER, eps, &first[d][i], &w[d][4 * i]);
	}
	// collapse z for every lattice (x,y), then y for every lattice x, then x (ITK's order); slices
	// are independent, so the result does not depend on the thread count
	if (!rc) {
#pragma omp parallel
		{
			float *c2 = (float *)malloc((size_t)nx * ny * sizeof(float)), *c1 = (float *)malloc(nx * sizeof(float));
			if (!c2 || !c1) {
#pragma omp atomic write
				rc = 1;
			}
#pragma omp for schedule(static)
			for (int z = 0; z < m[2]; z++) {
				if (!c2 || !c1) continue;
				const float *wz = &w[2][4 * z];
				for (int l = 0; l < nx * ny; l++) {
					float v = 0.0f;
					for (int k = 0; k < 4; k++) v += phi[(size_t)(first[2][z] + k) * nx * ny + l] * wz[k];
					c2[l] = v;
				}
				size_t o = (size_t)z * m[1] * m[0];
				for (int y = 0; y < m[1]; y++) {
					const float *wy = &w[1][4 * y];
					for (int l = 0; l < nx; l++) {
						float v = 0.0f;
						for (int j = 0; j < 4; j++) v += c2[(first[1][y] + j) * nx + l] * wy[j];
						c1[l] = v;
					}
					for (int x = 0; x < m[0]; x++, o++) {
						const float *wx = &w[0][4 * x], *c = &c1[first[0][x]];
						out[o] = c[0] * wx[0] + c[1] * wx[1] + c[2] * wx[2] + c[3] * wx[3];
					}
				}
			}
			free(c2); free(c1);
		}
	}
	for (int d = 0; d < 3; d++) { free(first[d]); free(w[d]); }
	return rc;
}

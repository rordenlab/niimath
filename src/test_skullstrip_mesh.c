// Milestone 2 gate: every surface primitive checked against an analytic answer.
// No AFNI dependency -- these are the tests the plan lists as committed with niimath.
//
//   cc -O2 -DHAVE_ZLIB -o t test_skullstrip_mesh.c skullstrip.c nifti_io.c -lz -lm

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "skullstrip.h"

static int fails = 0;
static void ok(int cond, const char *what, const char *detail) {
	printf("  [%s] %s%s%s\n", cond ? "PASS" : "FAIL", what, detail ? " -- " : "", detail ? detail : "");
	if (!cond)
		fails++;
}
static void okf(int cond, const char *what, const char *fmt, double a, double b) {
	char buf[256];
	snprintf(buf, sizeof buf, fmt, a, b);
	ok(cond, what, buf);
}

// ---- topology ---------------------------------------------------------------

static void test_topology(int ld) {
	ss_mesh m;
	char tag[64];
	snprintf(tag, sizeof tag, "ld=%d", ld);
	if (ss_icosphere(ld, 100.0f, NULL, &m)) {
		printf("  [FAIL] icosphere ld=%d construction\n", ld);
		fails++;
		return;
	}
	printf("ld=%d: nv=%d nt=%d\n", ld, m.nv, m.nt);
	okf(m.nv == 10 * ld * ld + 2, "vertex count == 10*ld^2+2", "%.0f vs %.0f",
			(double)m.nv, (double)(10 * ld * ld + 2));
	okf(m.nt == 20 * ld * ld, "triangle count == 20*ld^2", "%.0f vs %.0f",
			(double)m.nt, (double)(20 * ld * ld));

	// Every undirected edge must have exactly two incident faces, and every directed
	// edge exactly one -- that is consistent winding AND closedness in one check.
	{
		long long ne = 0, bad_dir = 0;
		int *ea = malloc(sizeof(int) * 3 * m.nt), *eb = malloc(sizeof(int) * 3 * m.nt);
		for (int t = 0; t < m.nt; t++)
			for (int e = 0; e < 3; e++) {
				ea[3 * t + e] = m.t[3 * t + e];
				eb[3 * t + e] = m.t[3 * t + (e + 1) % 3];
			}
		ne = 3 * m.nt;
		// count directed duplicates via per-vertex adjacency scan
		for (long long i = 0; i < ne; i++) {
			int a = ea[i], b = eb[i], fwd = 0, rev = 0;
			for (long long q = 0; q < ne; q++) {
				if (ea[q] == a && eb[q] == b) fwd++;
				if (ea[q] == b && eb[q] == a) rev++;
			}
			if (fwd != 1 || rev != 1) { bad_dir++; break; }
		}
		ok(bad_dir == 0, "each directed edge appears exactly once (closed + consistently wound)", NULL);
		free(ea); free(eb);
	}

	// Euler characteristic V - E + F == 2. E = 3F/2 for a closed triangulation.
	{
		int E = 3 * m.nt / 2;
		int chi = m.nv - E + m.nt;
		okf(chi == 2, "Euler characteristic == 2", "chi=%.0f (E=%.0f)", (double)chi, (double)E);
	}

	// Adjacency symmetric and duplicate-free.
	{
		int asym = 0, dup = 0;
		for (int i = 0; i < m.nv; i++) {
			for (int k = m.nbr_off[i]; k < m.nbr_off[i + 1]; k++) {
				int j = m.nbr[k], found = 0;
				if (k > m.nbr_off[i] && m.nbr[k] == m.nbr[k - 1]) dup++;
				if (j == i) dup++;
				for (int q = m.nbr_off[j]; q < m.nbr_off[j + 1]; q++)
					if (m.nbr[q] == i) { found = 1; break; }
				if (!found) asym++;
			}
		}
		ok(asym == 0 && dup == 0, "adjacency symmetric and duplicate-free", NULL);
	}

	// Degrees: exactly 12 vertices of degree 5, the rest degree 6.
	{
		int d5 = 0, d6 = 0, other = 0;
		for (int i = 0; i < m.nv; i++) {
			int d = m.nbr_off[i + 1] - m.nbr_off[i];
			if (d == 5) d5++; else if (d == 6) d6++; else other++;
		}
		okf(d5 == 12 && other == 0, "exactly 12 degree-5 vertices, rest degree-6",
				"d5=%.0f other=%.0f", (double)d5, (double)other);
	}

	// Normals point outward: for a sphere centred at the origin, n.v > 0 everywhere.
	{
		double worst = 1e30;
		for (int i = 0; i < m.nv; i++) {
			const float *v = m.v + 3 * i, *n = m.nrm + 3 * i;
			double d = (v[0] * n[0] + v[1] * n[1] + v[2] * n[2]) /
			           sqrt(v[0] * (double)v[0] + v[1] * (double)v[1] + v[2] * (double)v[2]);
			if (d < worst) worst = d;
		}
		okf(worst > 0.99, "vertex normals point outward", "min cos = %.5f (want > %.2f)", worst, 0.99);
	}
	ss_mesh_free(&m);
}

// ---- analytic convergence ---------------------------------------------------

static void mesh_area_volume(const ss_mesh *m, double *area, double *vol) {
	double A = 0, V = 0;
	for (int t = 0; t < m->nt; t++) {
		const float *a = m->v + 3 * m->t[3 * t + 0];
		const float *b = m->v + 3 * m->t[3 * t + 1];
		const float *c = m->v + 3 * m->t[3 * t + 2];
		double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
		double w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
		double n[3] = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]};
		A += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		// signed volume via the divergence theorem
		V += (a[0] * n[0] + a[1] * n[1] + a[2] * n[2]) / 6.0;
	}
	*area = A;
	*vol = V;
}

static void test_convergence(void) {
	const double R = 100.0;
	const double Aexact = 4.0 * M_PI * R * R, Vexact = 4.0 / 3.0 * M_PI * R * R * R;
	double prevA = 1e30, prevV = 1e30;
	int mono = 1;
	printf("sphere convergence (R=100):\n");
	for (int ld = 4; ld <= 32; ld *= 2) {
		ss_mesh m;
		double A, V;
		if (ss_icosphere(ld, (float)R, NULL, &m))
			continue;
		mesh_area_volume(&m, &A, &V);
		printf("   ld=%-3d area err %8.5f%%   volume err %8.5f%%   (V sign %s)\n", ld,
				100 * (A - Aexact) / Aexact, 100 * (V - Vexact) / Vexact, V > 0 ? "+" : "-");
		if (fabs(A - Aexact) / Aexact > prevA || fabs(V - Vexact) / Vexact > prevV)
			mono = 0;
		prevA = fabs(A - Aexact) / Aexact;
		prevV = fabs(V - Vexact) / Vexact;
		if (ld == 32) {
			okf(prevA < 0.002, "area converges to 4*pi*R^2", "err %.5f%% < %.2f%%", 100 * prevA, 0.2);
			okf(prevV < 0.003, "volume converges to 4/3*pi*R^3", "err %.5f%% < %.2f%%", 100 * prevV, 0.3);
			ok(V > 0, "signed volume positive => outward winding", NULL);
		}
		ss_mesh_free(&m);
	}
	ok(mono, "error decreases monotonically with density", NULL);
}

// ---- rasterisation ----------------------------------------------------------

// An axis-aligned cube as 12 triangles, outward wound, spanning [x0,x1]^3.
static int make_cube(ss_mesh *m, float x0, float x1) {
	static const int f[12][3] = {
		{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},
		{2,3,7},{2,7,6},{1,2,6},{1,6,5},{0,4,7},{0,7,3}};
	memset(m, 0, sizeof(*m));
	m->nv = 8; m->nt = 12;
	m->v = malloc(sizeof(float) * 24);
	m->t = malloc(sizeof(int) * 36);
	m->nrm = malloc(sizeof(float) * 24);
	if (!m->v || !m->t || !m->nrm) return 1;
	{
		float c[8][3] = {{x0,x0,x0},{x1,x0,x0},{x1,x1,x0},{x0,x1,x0},
		                 {x0,x0,x1},{x1,x0,x1},{x1,x1,x1},{x0,x1,x1}};
		memcpy(m->v, c, sizeof(c));
	}
	memcpy(m->t, f, sizeof(f));
	m->nbr_off = calloc(9, sizeof(int));
	m->nbr = calloc(1, sizeof(int));
	if (!m->nbr_off || !m->nbr) return 1;
	ss_mesh_normals(m);
	return 0;
}

static void test_raster(void) {
	const int N = 32;
	unsigned char *mask = malloc((size_t)N * N * N);
	printf("rasterisation:\n");

	// Cube from 8.0 to 20.0. Faces land EXACTLY on voxel centres, which is the
	// ambiguous case. The pinned convention is half-open in every axis: a voxel is
	// inside iff its centre lies in [lo, hi). So centres 8..19 are inside => 12^3.
	{
		ss_mesh c;
		long long cnt = 0;
		if (make_cube(&c, 8.0f, 20.0f) || ss_mesh_rasterize(&c, N, N, N, mask)) {
			ok(0, "cube rasterisation ran", NULL);
		} else {
			for (long long q = 0; q < (long long)N * N * N; q++) cnt += mask[q] ? 1 : 0;
			okf(cnt == 12LL * 12 * 12, "axis-aligned cube exact voxel count (half-open)",
					"%.0f vs %.0f", (double)cnt, (double)(12 * 12 * 12));
			// and exactly the right voxels, not merely the right number
			{
				int wrong = 0;
				for (int k = 0; k < N; k++) for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
					int want = (i >= 8 && i < 20 && j >= 8 && j < 20 && k >= 8 && k < 20);
					if ((mask[i + j * N + (long long)k * N * N] != 0) != want) wrong++;
				}
				okf(wrong == 0, "cube occupies exactly the right voxels", "%.0f wrong", (double)wrong, 0.0);
			}
		}
		ss_mesh_free(&c);
	}

	// Translated by a half voxel: no centre lies on a face, so the half-open rule is
	// not exercised and the answer is unambiguous: centres 9..20 => 12^3.
	{
		ss_mesh c;
		long long cnt = 0;
		if (!make_cube(&c, 8.5f, 20.5f) && !ss_mesh_rasterize(&c, N, N, N, mask)) {
			for (long long q = 0; q < (long long)N * N * N; q++) cnt += mask[q] ? 1 : 0;
			okf(cnt == 12LL * 12 * 12, "half-voxel translated cube exact count",
					"%.0f vs %.0f", (double)cnt, (double)(12 * 12 * 12));
		} else ok(0, "translated cube ran", NULL);
		ss_mesh_free(&c);
	}

	// Sphere volume convergence through the rasteriser.
	{
		ss_mesh s;
		float ctr[3] = {16.0f, 16.0f, 16.0f};
		long long cnt = 0;
		double exact = 4.0 / 3.0 * M_PI * 10.0 * 10.0 * 10.0;
		if (!ss_icosphere(24, 10.0f, ctr, &s) && !ss_mesh_rasterize(&s, N, N, N, mask)) {
			for (long long q = 0; q < (long long)N * N * N; q++) cnt += mask[q] ? 1 : 0;
			okf(fabs(cnt - exact) / exact < 0.02, "rasterised sphere volume within 2% of analytic",
					"%.0f vs %.1f", (double)cnt, exact);
		} else ok(0, "sphere raster ran", NULL);
		ss_mesh_free(&s);
	}

	// Reversing triangle enumeration order must not change the mask -- that is the
	// watertightness/half-open-edge property, and it is where a naive rasteriser leaks.
	{
		ss_mesh s, r;
		float ctr[3] = {16.0f, 16.0f, 16.0f};
		unsigned char *m2 = malloc((size_t)N * N * N);
		if (!ss_icosphere(16, 11.3f, ctr, &s)) {
			memset(&r, 0, sizeof r);
			r.nv = s.nv; r.nt = s.nt;
			r.v = malloc(sizeof(float) * 3 * s.nv);
			r.t = malloc(sizeof(int) * 3 * s.nt);
			r.nrm = malloc(sizeof(float) * 3 * s.nv);
			r.nbr_off = calloc((size_t)s.nv + 1, sizeof(int));
			r.nbr = calloc(1, sizeof(int));
			memcpy(r.v, s.v, sizeof(float) * 3 * s.nv);
			memcpy(r.nrm, s.nrm, sizeof(float) * 3 * s.nv);
			for (int t = 0; t < s.nt; t++)
				memcpy(r.t + 3 * (s.nt - 1 - t), s.t + 3 * t, sizeof(int) * 3);
			ss_mesh_rasterize(&s, N, N, N, mask);
			ss_mesh_rasterize(&r, N, N, N, m2);
			ok(memcmp(mask, m2, (size_t)N * N * N) == 0,
					"mask independent of triangle enumeration order", NULL);
			ss_mesh_free(&r);
		}
		ss_mesh_free(&s);
		free(m2);
	}

	// No interior holes in a rasterised closed surface.
	{
		ss_mesh s;
		float ctr[3] = {16.0f, 16.0f, 16.0f};
		int holes = 0;
		if (!ss_icosphere(20, 12.0f, ctr, &s) && !ss_mesh_rasterize(&s, N, N, N, mask)) {
			for (int k = 1; k < N - 1; k++) for (int j = 1; j < N - 1; j++) for (int i = 1; i < N - 1; i++) {
				long long p = i + j * N + (long long)k * N * N;
				if (mask[p]) continue;
				if (mask[p-1] && mask[p+1] && mask[p-N] && mask[p+N] &&
						mask[p-(long long)N*N] && mask[p+(long long)N*N]) holes++;
			}
			okf(holes == 0, "no interior holes after rasterisation", "%.0f holes", (double)holes, 0.0);
		} else ok(0, "hole check ran", NULL);
		ss_mesh_free(&s);
	}
	free(mask);
}

// ---- self-intersection ------------------------------------------------------

static void test_intersect(void) {
	printf("self-intersection:\n");
	{
		ss_mesh s;
		if (!ss_icosphere(16, 50.0f, NULL, &s)) {
			long long h = ss_mesh_self_intersections(&s);
			okf(h == 0, "clean sphere reports zero intersections", "%.0f found", (double)h, 0.0);
			ss_mesh_free(&s);
		}
	}
	{
		// Fold one vertex far through the surface: that must be detected.
		ss_mesh s;
		if (!ss_icosphere(16, 50.0f, NULL, &s)) {
			long long h0 = ss_mesh_self_intersections(&s);
			int vi = s.nv / 2;
			s.v[3 * vi + 0] *= -1.6f; s.v[3 * vi + 1] *= -1.6f; s.v[3 * vi + 2] *= -1.6f;
			long long h = ss_mesh_self_intersections(&s);
			okf(h0 == 0 && h > 0, "folded vertex is detected", "before %.0f, after %.0f",
					(double)h0, (double)h);
			ss_mesh_free(&s);
		}
	}
	{
		// Two triangles sharing an edge, bent: adjacency must not be a false positive.
		ss_mesh m;
		memset(&m, 0, sizeof m);
		m.nv = 4; m.nt = 2;
		m.v = malloc(sizeof(float) * 12);
		m.t = malloc(sizeof(int) * 6);
		m.nrm = malloc(sizeof(float) * 12);
		m.nbr_off = calloc(5, sizeof(int));
		m.nbr = calloc(1, sizeof(int));
		{
			float v[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
			int t[2][3] = {{0,1,2},{0,2,3}};
			memcpy(m.v, v, sizeof v);
			memcpy(m.t, t, sizeof t);
		}
		{
			long long h = ss_mesh_self_intersections(&m);
			okf(h == 0, "edge-sharing triangles are not a false positive", "%.0f found", (double)h, 0.0);
		}
		ss_mesh_free(&m);
	}
}

// ---- smoothing --------------------------------------------------------------

static void test_smooth(void) {
	ss_mesh m;
	printf("smoothing:\n");
	if (ss_icosphere(12, 40.0f, NULL, &m))
		return;
	{
		float *scratch = malloc(sizeof(float) * 3 * m.nv);
		double r0 = 0, r1 = 0;
		for (int i = 0; i < m.nv; i++)
			r0 += sqrt(m.v[3*i]*(double)m.v[3*i] + m.v[3*i+1]*(double)m.v[3*i+1] + m.v[3*i+2]*(double)m.v[3*i+2]);
		r0 /= m.nv;
		for (int it = 0; it < 5; it++)
			ss_mesh_smooth(&m, 0.5f, scratch);
		for (int i = 0; i < m.nv; i++)
			r1 += sqrt(m.v[3*i]*(double)m.v[3*i] + m.v[3*i+1]*(double)m.v[3*i+1] + m.v[3*i+2]*(double)m.v[3*i+2]);
		r1 /= m.nv;
		// A sphere is already the smooth fixed shape: smoothing shrinks it slightly
		// but must not distort it, so the radius spread must stay tiny.
		{
			double var = 0;
			for (int i = 0; i < m.nv; i++) {
				double r = sqrt(m.v[3*i]*(double)m.v[3*i] + m.v[3*i+1]*(double)m.v[3*i+1] + m.v[3*i+2]*(double)m.v[3*i+2]);
				var += (r - r1) * (r - r1);
			}
			var = sqrt(var / m.nv);
			okf(var / r1 < 0.005, "smoothing preserves sphericity", "rel sd %.6f < %.3f", var / r1, 0.005);
			okf(r1 < r0, "smoothing shrinks (expected for neighbour averaging)", "%.3f -> %.3f", r0, r1);
		}
		free(scratch);
	}
	ss_mesh_free(&m);
}

int main(void) {
	printf("== Milestone 2: surface primitives ==\n\n");
	printf("topology:\n");
	test_topology(1);
	test_topology(3);
	test_topology(20); // the -no_use_edge density
	printf("\n");
	test_convergence();
	printf("\n");
	test_raster();
	printf("\n");
	test_intersect();
	printf("\n");
	test_smooth();
	printf("\n%s (%d failure%s)\n", fails ? "FAILED" : "ALL PASSED", fails, fails == 1 ? "" : "s");
	return fails ? 1 : 0;
}

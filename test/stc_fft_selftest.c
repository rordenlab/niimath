// Isolation test for the batched Stockham FFT primitive in src/stc.c.
//
//   cc -O2 -std=gnu99 -Isrc -o /tmp/stcfft test/stc_fft_selftest.c -lm && /tmp/stcfft
//
// Checks the forward transform against a naive long-double DFT and the inverse round trip at
// EVERY length stc_fft_len() can return for nt = 5..4000 (35 distinct lengths, 10 through 4096),
// over random / impulse / constant / Nyquist inputs, with DIFFERENT data in each SIMD lane so
// cross-lane contamination cannot hide.  Also re-checks stc_fft_len() against the boundary cases
// measured from the oracle in moco_bench's test/stc_reference_manifest.md.  Pass an argument for per-length
// output.  Runs in a few seconds.  `make test` builds and runs it via the `stc-test` target --
// it needs its own compile step, which is the only reason it is not inside release_smoke.py.
//
// It reaches the static kernel by #include-ing stc.c, so nii_stc() is never called and the one
// core.c symbol stc.c references is stubbed below rather than linked.
#include <stdio.h>

int nii_mul_size(size_t a, size_t b, size_t *out);
int nii_mul_size(size_t a, size_t b, size_t *out) {
	if (a && b > (size_t)-1 / a) return 1;
	*out = a * b;
	return 0;
}

#include "stc.c"

static int fails = 0;

static void naive_dft(const double *xr, const double *xi, double *yr, double *yi,
                      int64_t n, double sgn) {
	for (int64_t k = 0; k < n; k++) {
		long double sr = 0.0L, si = 0.0L;
		for (int64_t t = 0; t < n; t++) {
			long double th = sgn * 2.0L * 3.14159265358979323846264338327950288L
			                 * (long double)((t * k) % n) / (long double)n;
			long double c = cosl(th), s = sinl(th);
			sr += (long double)xr[t] * c - (long double)xi[t] * s;
			si += (long double)xr[t] * s + (long double)xi[t] * c;
		}
		yr[k] = (double)sr;
		yi[k] = (double)si;
	}
}

// Run one batched transform on `n`-point data, lane l carrying series l.
static void run(int64_t n, double *lane_r[STC_LANES], double *lane_i[STC_LANES],
                double *out_r[STC_LANES], double *out_i[STC_LANES], double sgn) {
	stc_plan p;
	memset(&p, 0, sizeof(p));
	p.n = n;
	p.nrad = stc_factor(n, p.rad, 40);
	if (!p.nrad) { printf("FAIL: cannot factor %lld\n", (long long)n); fails++; return; }
	p.wr = (double *)malloc(n * sizeof(double));
	p.wi = (double *)malloc(n * sizeof(double));
	for (int64_t j = 0; j < n; j++) {
		double th = 2.0 * M_PI * (double)j / (double)n;
		p.wr[j] = cos(th);
		p.wi[j] = sin(th);
	}
	double *buf = (double *)malloc((size_t)4 * n * STC_LANES * sizeof(double));
	double *ar = buf, *ai = buf + n * STC_LANES;
	double *br = buf + 2 * n * STC_LANES, *bi = buf + 3 * n * STC_LANES;
	for (int64_t t = 0; t < n; t++)
		for (int l = 0; l < STC_LANES; l++) {
			ar[t * STC_LANES + l] = lane_r[l][t];
			ai[t * STC_LANES + l] = lane_i[l][t];
		}
	double *fr, *fi;
	stc_fft(ar, ai, br, bi, n, p.rad, p.nrad, p.wr, p.wi, sgn, &fr, &fi);
	for (int64_t t = 0; t < n; t++)
		for (int l = 0; l < STC_LANES; l++) {
			out_r[l][t] = fr[t * STC_LANES + l];
			out_i[l][t] = fi[t * STC_LANES + l];
		}
	free(buf);
	free(p.wr);
	free(p.wi);
}

static unsigned long rs = 12345;
static double rnd(void) {
	rs = rs * 6364136223846793005UL + 1442695040888963407UL;
	return (double)((rs >> 11) & 0x1FFFFF) / 1048576.0 - 1.0;
}

static void check_len(int64_t n, int verbose) {
	double *lr[STC_LANES], *li[STC_LANES], *or_[STC_LANES], *oi[STC_LANES];
	double *nr = (double *)malloc(n * sizeof(double)), *ni = (double *)malloc(n * sizeof(double));
	for (int l = 0; l < STC_LANES; l++) {
		lr[l] = (double *)malloc(n * sizeof(double));
		li[l] = (double *)malloc(n * sizeof(double));
		or_[l] = (double *)malloc(n * sizeof(double));
		oi[l] = (double *)malloc(n * sizeof(double));
	}
	double worst_f = 0.0, worst_rt = 0.0;
	for (int trial = 0; trial < 4; trial++) {
		double scale = 0.0;
		for (int l = 0; l < STC_LANES; l++)
			for (int64_t t = 0; t < n; t++) {
				if (trial == 0) { lr[l][t] = rnd(); li[l][t] = rnd(); }
				else if (trial == 1) { lr[l][t] = (t == (l * 3) % n) ? 1.0 : 0.0; li[l][t] = 0.0; }
				else if (trial == 2) { lr[l][t] = 1.0 + l; li[l][t] = 0.0; }
				else { lr[l][t] = ((t & 1) ? -1.0 : 1.0) * (1 + l); li[l][t] = 0.0; }
				double m = fabs(lr[l][t]) > fabs(li[l][t]) ? fabs(lr[l][t]) : fabs(li[l][t]);
				if (m > scale) scale = m;
			}
		double *save_r[STC_LANES], *save_i[STC_LANES];
		for (int l = 0; l < STC_LANES; l++) {
			save_r[l] = (double *)malloc(n * sizeof(double));
			save_i[l] = (double *)malloc(n * sizeof(double));
			memcpy(save_r[l], lr[l], n * sizeof(double));
			memcpy(save_i[l], li[l], n * sizeof(double));
		}
		run(n, lr, li, or_, oi, -1.0);
		// forward vs naive long-double DFT, per lane (also proves lanes do not cross-talk)
		for (int l = 0; l < STC_LANES; l++) {
			naive_dft(save_r[l], save_i[l], nr, ni, n, -1.0);
			for (int64_t k = 0; k < n; k++) {
				double e = fabs(or_[l][k] - nr[k]) + fabs(oi[l][k] - ni[k]);
				double rel = e / (scale * (double)n);
				if (rel > worst_f) worst_f = rel;
			}
		}
		// inverse round trip
		double *rr[STC_LANES], *ri[STC_LANES];
		for (int l = 0; l < STC_LANES; l++) {
			rr[l] = (double *)malloc(n * sizeof(double));
			ri[l] = (double *)malloc(n * sizeof(double));
		}
		run(n, or_, oi, rr, ri, +1.0);
		for (int l = 0; l < STC_LANES; l++) {
			for (int64_t t = 0; t < n; t++) {
				double e = fabs(rr[l][t] / (double)n - save_r[l][t])
				         + fabs(ri[l][t] / (double)n - save_i[l][t]);
				double rel = e / scale;
				if (rel > worst_rt) worst_rt = rel;
			}
			free(rr[l]); free(ri[l]); free(save_r[l]); free(save_i[l]);
		}
	}
	const double TOL_F = 1e-13, TOL_R = 1e-13;
	int bad = (worst_f > TOL_F) || (worst_rt > TOL_R);
	if (bad || verbose)
		printf("%s n=%-6lld fwd=%.2e roundtrip=%.2e\n", bad ? "FAIL" : "  ok",
		       (long long)n, worst_f, worst_rt);
	if (bad) fails++;
	for (int l = 0; l < STC_LANES; l++) { free(lr[l]); free(li[l]); free(or_[l]); free(oi[l]); }
	free(nr); free(ni);
}

int main(int argc, char **argv) {
	int verbose = (argc > 1);
	// every length stc_fft_len() can return, for nt = 5 .. 4000
	int64_t prev = 0;
	int nlen = 0;
	for (int64_t nt = 5; nt <= 4000; nt++) {
		int64_t n = stc_fft_len(nt);
		if (n == prev) continue;
		prev = n;
		nlen++;
		check_len(n, verbose);
	}
	printf("lengths exercised: %d\n", nlen);
	// the length rule itself, against the measured boundary cases in the manifest
	struct { int64_t nt, n; } known[] = {
		{5, 10}, {6, 10}, {7, 12}, {8, 12}, {9, 16}, {10, 16}, {11, 16}, {12, 16}, {13, 20},
		{16, 20}, {29, 40}, {61, 80}, {64, 80}, {67, 80}, {125, 160}, {131, 160}, {254, 320},
		{476, 480}, {478, 512}, {508, 512}};
	for (size_t i = 0; i < sizeof(known) / sizeof(known[0]); i++) {
		int64_t got = stc_fft_len(known[i].nt);
		if (got != known[i].n) {
			printf("FAIL: stc_fft_len(%lld) = %lld, expected %lld\n",
			       (long long)known[i].nt, (long long)got, (long long)known[i].n);
			fails++;
		}
	}
	printf(fails ? "FAILURES: %d\n" : "all FFT checks passed (%d failures)\n", fails);
	return fails ? 1 : 0;
}

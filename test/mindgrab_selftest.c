// Closed-form self-test for -mindgrab. Run by `make test` when BRAINCHOP=1; it needs neither
// brainchop, tinygrad, AFNI nor any test data. It #includes mindgrab.c to reach the static
// kernels, the same way test/stc_fft_selftest.c reaches stc.c's FFT, and overrides MG_DIM so
// the real convolution runs on a 32^3 grid (2 MB) instead of the production 256^3 (2 GB).
//
// The full-volume comparison against the brainchop-cli oracle is a separate, external step:
// test/mindgrab_parity.c.
#define MG_DIM 32
#include "../src/mindgrab.c"

#include <float.h>

static int failures = 0;

static void check(int cond, const char *what) {
	if (!cond) {
		printf("FAIL: %s\n", what);
		failures++;
	}
}

static void checkf(double got, double want, double tol, const char *what) {
	double d = fabs(got - want);
	if (!(d <= tol)) {
		printf("FAIL: %s (got %.9g, want %.9g, |d| %.3g > %.3g)\n", what, got, want, d, tol);
		failures++;
	}
}

static long ulp_dist(float a, float b) {
	int32_t ai, bi;
	memcpy(&ai, &a, 4);
	memcpy(&bi, &b, 4);
	if (ai < 0) ai = (int32_t)0x80000000 - ai;
	if (bi < 0) bi = (int32_t)0x80000000 - bi;
	long d = (long)ai - (long)bi;
	return d < 0 ? -d : d;
}

// ---------------------------------------------------------------------------- tanh and GELU

static void test_tanh(void) {
	long worst = 0;
	float worst_at = 0.0f;
	for (int i = -240000; i <= 240000; i++) {
		float x = (float)i * 5e-5f; // -12 .. 12
		long d = ulp_dist(mg_tanhf(x), tanhf(x));
		if (d > worst) {
			worst = d;
			worst_at = x;
		}
	}
	printf("  worst tanh ULP %ld at x=%g\n", worst, worst_at);
	check(worst <= 8, "mg_tanhf within 8 float32 ULP of libm over [-12,12]");
	// In the clamped tail the rational lands up to 2 float32 ULP below 1.0f -- libm's own
	// tanhf(9) is 1 ULP below 1.0f, so this is a 1.2e-7 absolute error on a value the GELU
	// multiplies by 0.5x. Inside the budget, but not exact: assert the bound, not equality.
	check(mg_tanhf(20.0f) <= 1.0f && ulp_dist(mg_tanhf(20.0f), 1.0f) <= 2,
		  "mg_tanhf saturates to +1 within 2 ULP");
	check(mg_tanhf(-20.0f) >= -1.0f && ulp_dist(mg_tanhf(-20.0f), -1.0f) <= 2,
		  "mg_tanhf saturates to -1 within 2 ULP");
	check(mg_tanhf(0.0f) == 0.0f, "mg_tanhf(0) == 0");
}

static void test_gelu(void) {
	// Pinned from tinygrad's Tensor.gelu() -- the tanh approximation, NOT exact erf. Exact erf
	// gives -0.045500 / -0.154286 / 0.345731 / 1.954500 and would move boundary voxels.
	static const float xs[5] = {-2.0f, -0.5f, 0.0f, 0.5f, 2.0f};
	static const double want[5] = {-0.04540232, -0.15428598, 0.0, 0.345714, 1.9545976};
	for (int i = 0; i < 5; i++)
		checkf(mindgrab_gelu(xs[i]), want[i], 1e-6, "mindgrab_gelu matches tinygrad");
	// and against the defining formula in double precision
	for (int i = -60; i <= 60; i++) {
		double x = i * 0.1;
		double u = 0.7978845608028654 * (x + 0.044715 * x * x * x);
		double w = 0.5 * x * (1.0 + tanh(u));
		checkf(mindgrab_gelu((float)x), w, 3e-6 + 1e-6 * fabs(w), "gelu closed form");
	}
}

// -------------------------------------------------------------------------------- quantiles

static void test_quantile(void) {
	// numpy's default "linear" rule: h = q*(n-1); floor/frac; lerp with the >=0.5 branch.
	int64_t cum[256];
	// 100 values: 40 zeros, 60 tens -> sorted a[0..39]=0, a[40..99]=10
	for (int i = 0; i < 256; i++)
		cum[i] = (i < 10) ? 40 : 100;
	// q = 0.5 -> h = 49.5, a[49]=10, a[50]=10 -> 10
	checkf(mg_quantile(cum, 100, 0.5), 10.0, 1e-12, "quantile inside a run");
	// q = 0.39 -> h = 38.61, a[38]=0, a[39]=0 -> 0
	checkf(mg_quantile(cum, 100, 0.39), 0.0, 1e-12, "quantile below the step");
	// q = 39.0/99 -> h = 39.0 exactly, t = 0 -> a[39] = 0
	checkf(mg_quantile(cum, 100, 39.0 / 99.0), 0.0, 1e-9, "quantile exactly on the step");
	// q = 39.25/99 -> h = 39.25, t = 0.25 < 0.5 -> a + (b-a)t = 0 + 10*0.25
	checkf(mg_quantile(cum, 100, 39.25 / 99.0), 2.5, 1e-9, "quantile lerp, t < 0.5");
	// q = 39.75/99 -> t = 0.75 >= 0.5 -> b - (b-a)(1-t) = 10 - 10*0.25
	checkf(mg_quantile(cum, 100, 39.75 / 99.0), 7.5, 1e-9, "quantile lerp, t >= 0.5");
}

// ---------------------------------------------------------------------------- weight table

static void test_weights(void) {
	check(MINDGRAB_NPARAM == 146237, "146237 parameters");
	check(MINDGRAB_LOFF(0) == 0, "layer 0 offset");
	check(MINDGRAB_LOFF(1) == 405, "layer 1 offset");
	check(MINDGRAB_CLS_OFF == 405 + 24 * 6075, "classifier offset");
	check(MINDGRAB_BIAS_OFF == MINDGRAB_CLS_OFF + 30, "bias offset");
	double s = 0.0, q = 0.0;
	for (int i = 0; i < MINDGRAB_NPARAM; i++) {
		s += mindgrab_weights[i];
		q += (double)mindgrab_weights[i] * mindgrab_weights[i];
	}
	// Pinned against the SHA-256-verified checkpoint (see src/mindgrab.LICENSE): a truncated or
	// re-exported table changes these, a reordered one changes only q.
	checkf(s, -1322.458802328659, 1e-6, "weight sum");
	checkf(q, 5363.111318030771, 1e-3, "weight sum of squares");
	// the padded 16th output column must be exactly zero, or GroupNorm would normalise garbage
	for (int l = 0; l < MINDGRAB_NHIDDEN; l++) {
		int cin = 0;
		float *w = mg_pack_weights(l, &cin);
		check(w != NULL, "pack weights");
		if (!w) return;
		check(cin == (l == 0 ? 1 : MG_C), "packed cin");
		int zero_pad = 1, same = 1;
		const float *raw = mindgrab_weights + MINDGRAB_LOFF(l);
		for (int t = 0; t < 27; t++)
			for (int ic = 0; ic < cin; ic++) {
				for (int oc = MG_C; oc < MG_CS; oc++)
					if (w[(t * cin + ic) * MG_CS + oc] != 0.0f) zero_pad = 0;
				for (int oc = 0; oc < MG_C; oc++)
					if (w[(t * cin + ic) * MG_CS + oc] != raw[(t * cin + ic) * MG_C + oc]) same = 0;
			}
		check(zero_pad, "packed weights zero-pad the 16th output channel");
		check(same, "packed weights preserve the table values");
		free(w);
	}
}

// ------------------------------------------------------------------------- convolution shape

// An impulse at one voxel must reproduce the kernel, mirrored, at the dilated offsets -- this is
// the check that would catch the reversed spatial axes (the tensor's (D,H,W) are the NIfTI
// (x,y,z)) and any padding/dilation slip at the volume faces.
static void test_conv_impulse(int dil, int cx, int cy, int cz, const char *where) {
	float *src = (float *)calloc(MG_NVOX * MG_CS, sizeof(float));
	float *dst = (float *)malloc(MG_NVOX * MG_CS * sizeof(float));
	float *w = (float *)calloc((size_t)27 * MG_C * MG_CS, sizeof(float));
	if (!src || !dst || !w) {
		printf("FAIL: impulse alloc\n");
		failures++;
		free(src);
		free(dst);
		free(w);
		return;
	}
	const int ic = 3; // arbitrary input channel carrying the impulse
	src[(((size_t)cz * MG_DIM + cy) * MG_DIM + cx) * MG_CS + ic] = 1.0f;
	for (int t = 0; t < 27; t++)
		for (int i = 0; i < MG_C; i++)
			for (int oc = 0; oc < MG_C; oc++)
				w[((size_t)t * MG_C + i) * MG_CS + oc] = (float)(t * 100 + i * 10 + oc + 1);
	mg_conv(src, dst, w, dil);
	size_t wrong = 0;
	for (int z = 0; z < MG_DIM; z++)
		for (int y = 0; y < MG_DIM; y++)
			for (int x = 0; x < MG_DIM; x++) {
				// output voxel v reads src at v + (t-1)*dil, so the impulse lands on v where
				// v = centre - (t-1)*dil
				int tz = 1 - (z - cz) / dil, ty = 1 - (y - cy) / dil, tx = 1 - (x - cx) / dil;
				int hit = ((z - cz) % dil == 0) && ((y - cy) % dil == 0) && ((x - cx) % dil == 0) &&
						  tz >= 0 && tz < 3 && ty >= 0 && ty < 3 && tx >= 0 && tx < 3;
				const float *o = dst + ((((size_t)z * MG_DIM + y) * MG_DIM) + x) * MG_CS;
				for (int oc = 0; oc < MG_CS; oc++) {
					float want = 0.0f;
					if (hit && oc < MG_C)
						want = w[(((size_t)(tz * 9 + ty * 3 + tx)) * MG_C + ic) * MG_CS + oc];
					if (o[oc] != want) wrong++;
				}
			}
	if (wrong) printf("FAIL: conv impulse %s dil=%d: %zu wrong values\n", where, dil, wrong);
	check(wrong == 0, "conv impulse response");
	free(src);
	free(dst);
	free(w);
}

// -------------------------------------------------------------------- GroupNorm, GELU, argmax

static void test_norm_gelu(void) {
	float *buf = (float *)calloc(MG_NVOX * MG_CS, sizeof(float));
	double *partial = (double *)malloc((size_t)MG_DIM * 2 * MG_C * sizeof(double));
	if (!buf || !partial) {
		printf("FAIL: norm alloc\n");
		failures++;
		free(buf);
		free(partial);
		return;
	}
	// channel 0: exactly two values, +-1 about a mean of 5 -> var 1
	// channel 1: constant -> var 0, so every output is gelu(0) == 0
	for (size_t v = 0; v < MG_NVOX; v++) {
		buf[v * MG_CS + 0] = (v & 1) ? 6.0f : 4.0f;
		buf[v * MG_CS + 1] = 7.0f;
	}
	mg_norm_gelu(buf, partial);
	double scale = 1.0 / sqrt(1.0 + MG_EPS);
	checkf(buf[0 * MG_CS + 0], mindgrab_gelu((float)(-1.0 * scale)), 1e-6, "groupnorm even voxel");
	checkf(buf[1 * MG_CS + 0], mindgrab_gelu((float)(1.0 * scale)), 1e-6, "groupnorm odd voxel");
	checkf(buf[0 * MG_CS + 1], 0.0, 0.0, "groupnorm constant channel -> 0");
	checkf(buf[0 * MG_CS + 2], 0.0, 0.0, "groupnorm empty channel stays 0");
	check(buf[0 * MG_CS + MG_C] == 0.0f, "padded channel stays 0");
	free(buf);
	free(partial);
}

static void test_argmax_tie(void) {
	// Ties go to class 0, matching numpy's argmax. Build an activation whose two logits are
	// identical by construction: zero activations make both logits 0 (the reference drops the
	// classifier bias, so there is nothing to break the tie).
	float *buf = (float *)calloc(MG_NVOX * MG_CS, sizeof(float));
	float *mask = (float *)malloc(MG_NVOX * sizeof(float));
	if (!buf || !mask) {
		printf("FAIL: argmax alloc\n");
		failures++;
		free(buf);
		free(mask);
		return;
	}
	mg_classify(buf, mask);
	size_t fg = 0;
	for (size_t v = 0; v < MG_NVOX; v++)
		if (mask[v] != 0.0f) fg++;
	check(fg == 0, "argmax tie resolves to class 0 (background)");
	free(buf);
	free(mask);
}

static void test_largest_component(void) {
	float *mask = (float *)calloc(MG_NVOX, sizeof(float));
	if (!mask) {
		printf("FAIL: component alloc\n");
		failures++;
		return;
	}
	// Two equal two-voxel 26-connected components. The reference keeps the first maximum,
	// which is also bwlabel's scan-order tie rule.
#define TVOX(x, y, z) (((size_t)(z) * MG_DIM + (y)) * MG_DIM + (x))
	mask[TVOX(1, 1, 1)] = 1.0f;
	mask[TVOX(2, 2, 2)] = 1.0f;
	mask[TVOX(20, 20, 20)] = 1.0f;
	mask[TVOX(21, 21, 21)] = 1.0f;
	size_t dim[3] = {MG_DIM, MG_DIM, MG_DIM};
	check(bwlabel(mask, 26, dim, true, false) == 1, "largest component returns one label");
	check(mask[TVOX(1, 1, 1)] == 1.0f && mask[TVOX(2, 2, 2)] == 1.0f,
		  "largest component keeps the first equal-size cluster");
	check(mask[TVOX(20, 20, 20)] == 0.0f && mask[TVOX(21, 21, 21)] == 0.0f,
		  "largest component removes the later equal-size cluster");
#undef TVOX
	free(mask);
}

int main(void) {
	printf("mindgrab selftest (MG_DIM=%d)\n", MG_DIM);
	test_tanh();
	test_gelu();
	test_quantile();
	test_weights();
	for (int d = 1; d <= 16; d *= 2) {
		test_conv_impulse(d, MG_DIM / 2, MG_DIM / 2, MG_DIM / 2, "interior");
		test_conv_impulse(d, 0, 0, 0, "corner");
		test_conv_impulse(d, MG_DIM - 1, MG_DIM - 1, MG_DIM - 1, "far corner");
	}
	test_norm_gelu();
	test_argmax_tie();
	test_largest_component();
	if (failures) {
		printf("mindgrab selftest: %d FAILURE(S)\n", failures);
		return 1;
	}
	printf("mindgrab selftest: all checks passed\n");
	return 0;
}

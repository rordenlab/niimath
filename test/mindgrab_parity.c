// Layer-by-layer parity harness for -mindgrab against dumps from the brainchop-cli
// reference (scripts/dump_mindgrab_oracle.py). Developer tool, not part of
// any build; it #includes mindgrab.c to reach the static kernels, the same way
// test/stc_fft_selftest.c reaches stc.c's FFT.
//
//   cc -O3 -ffast-math -fno-finite-math-only -I../src \
//      -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
//      -o mindgrab_parity mindgrab_parity.c ../src/mindgrab_weights.c ../src/bwlabel.c \
//      -L/opt/homebrew/opt/libomp/lib -lomp -lm
//   ./mindgrab_parity <oracle_dir>
//
// The oracle directory holds conformed.u8, normalized.f32, layerNN.f32 (channel-last,
// stride 15), logits.f32 and argmax.u8.
#include "../src/mindgrab.c"

#include <time.h>

static int failures = 0;

static void *slurp(const char *dir, const char *name, size_t want) {
	char path[1024];
	snprintf(path, sizeof(path), "%s/%s", dir, name);
	FILE *f = fopen(path, "rb");
	if (!f) return NULL;
	void *p = malloc(want);
	if (!p) {
		fclose(f);
		return NULL;
	}
	size_t got = fread(p, 1, want, f);
	fclose(f);
	if (got != want) {
		free(p);
		return NULL;
	}
	return p;
}

static double now(void) {
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}

// Max |a-b| between our padded buffer and an oracle dump of stride 15.
static double cmp_pad(const char *label, const float *ours, const float *ref, int nchan) {
	double mx = 0.0, sum = 0.0;
	size_t at = 0;
	for (size_t v = 0; v < MG_NVOX; v++)
		for (int c = 0; c < nchan; c++) {
			double d = fabs((double)ours[v * MG_CS + c] - (double)ref[v * nchan + c]);
			if (!isfinite(d)) d = INFINITY;
			sum += d * d;
			if (d > mx) {
				mx = d;
				at = v;
			}
		}
	printf("  %-10s max|d| %.6e  rms %.3e  (worst voxel %zu)\n", label, mx,
		   sqrt(sum / (double)(MG_NVOX * nchan)), at);
	return mx;
}

int main(int argc, char **argv) {
	if (argc < 2) {
		fprintf(stderr, "usage: %s <oracle_dir>\n", argv[0]);
		return 2;
	}
	const char *dir = argv[1];
	unsigned char *conformed = (unsigned char *)slurp(dir, "conformed.u8", MG_NVOX);
	if (!conformed) {
		fprintf(stderr, "cannot read %s/conformed.u8\n", dir);
		return 2;
	}
	if (argc > 2 && strcmp(argv[2], "bench") == 0) { // timing only, no oracle I/O
		float *m = (float *)malloc(MG_NVOX * sizeof(float));
		int reps = (argc > 3) ? atoi(argv[3]) : 1;
		for (int r = 0; r < reps; r++) {
			double t = now();
			if (mindgrab_segment(conformed, m)) return 2;
			size_t fg = 0;
			for (size_t i = 0; i < MG_NVOX; i++)
				if (m[i] != 0.0f) fg++;
			printf("run %d: %.2f s  fg %zu\n", r, now() - t, fg);
			fflush(stdout);
		}
		return 0;
	}
	float *a = (float *)malloc(MG_NVOX * MG_CS * sizeof(float));
	float *b = (float *)malloc(MG_NVOX * MG_CS * sizeof(float));
	double *partial = (double *)malloc((size_t)MG_DIM * 2 * MG_C * sizeof(double));
	float *wpack[MINDGRAB_NHIDDEN];
	memset(wpack, 0, sizeof(wpack));
	int wcin;
	int packed = 1;
	for (int i = 0; i < MINDGRAB_NHIDDEN; i++) {
		wpack[i] = mg_pack_weights(i, &wcin);
		if (!wpack[i]) {
			packed = 0;
			break;
		}
	}
	if (!a || !b || !partial || !packed) {
		fprintf(stderr, "out of memory\n");
		return 2;
	}

	double t0 = now();
	mg_qnormalize(conformed, b);
	float *nref = (float *)slurp(dir, "normalized.f32", MG_NVOX * sizeof(float));
	if (!nref) {
		fprintf(stderr, "cannot read %s/normalized.f32\n", dir);
		failures++;
	} else {
		size_t bad = 0;
		double mx = 0;
		for (size_t i = 0; i < MG_NVOX; i++) {
			if (b[i] != nref[i]) bad++;
			double d = fabs((double)b[i] - (double)nref[i]);
			if (!isfinite(d)) d = INFINITY;
			if (d > mx) mx = d;
		}
		printf("normalized: %zu/%zu voxels differ, max|d| %.3e\n", bad, MG_NVOX, mx);
		if (mx > 2e-7) failures++;
		free(nref);
	}

	mg_conv_first(b, a, wpack[0], mg_dilation[0]);
	mg_norm_gelu(a, partial);
	for (int l = 0; l < MINDGRAB_NHIDDEN; l++) {
		if (l > 0) {
			float *src = (l & 1) ? a : b;
			float *dst = (l & 1) ? b : a;
			mg_conv(src, dst, wpack[l], mg_dilation[l]);
			mg_norm_gelu(dst, partial);
		}
		char name[64];
		snprintf(name, sizeof(name), "layer%02d.f32", l);
		float *ref = (float *)slurp(dir, name, MG_NVOX * MG_C * sizeof(float));
		int required = (l == 0 || l == 1 || l == 5 || l == 24);
		if (!ref && required) {
			fprintf(stderr, "cannot read required %s/%s\n", dir, name);
			failures++;
		} else if (ref) {
			printf("layer %02d (dilation %2d)\n", l, mg_dilation[l]);
			if (cmp_pad(name, (l & 1) ? b : a, ref, MG_C) > 1e-4) failures++;
			free(ref);
		}
	}
	const float *last = ((MINDGRAB_NHIDDEN - 1) & 1) ? b : a;

	// logits, before argmax
	float *lref = (float *)slurp(dir, "logits.f32", MG_NVOX * 2 * sizeof(float));
	if (!lref) {
		fprintf(stderr, "cannot read %s/logits.f32\n", dir);
		failures++;
	} else {
		const float *w = mindgrab_weights + MINDGRAB_CLS_OFF;
		double mx = 0.0;
		double minmargin = 1e30;
		for (size_t i = 0; i < MG_NVOX; i++) {
			const float *x = last + i * MG_CS;
			float l0 = 0.0f, l1 = 0.0f; // the reference drops the classifier bias
			for (int c = 0; c < MG_C; c++) {
				l0 += x[c] * w[c * 2];
				l1 += x[c] * w[c * 2 + 1];
			}
			double d0 = fabs((double)l0 - (double)lref[i * 2]);
			double d1 = fabs((double)l1 - (double)lref[i * 2 + 1]);
			if (!isfinite(d0)) d0 = INFINITY;
			if (!isfinite(d1)) d1 = INFINITY;
			if (d0 > mx) mx = d0;
			if (d1 > mx) mx = d1;
			double m = fabs((double)lref[i * 2] - (double)lref[i * 2 + 1]);
			if (m < minmargin) minmargin = m;
		}
		printf("logits: max|d| %.6e   oracle min margin %.6e\n", mx, minmargin);
		if (mx > 1e-4) failures++;
		free(lref);
	}

	float *mask = (float *)malloc(MG_NVOX * sizeof(float));
	if (!mask) {
		fprintf(stderr, "out of memory\n");
		return 2;
	}
	mg_classify(last, mask);
	unsigned char *aref = (unsigned char *)slurp(dir, "argmax.u8", MG_NVOX);
	if (!aref) {
		fprintf(stderr, "cannot read %s/argmax.u8\n", dir);
		failures++;
	} else {
		size_t bad = 0, fg = 0;
		for (size_t i = 0; i < MG_NVOX; i++) {
			if ((mask[i] != 0.0f) != (aref[i] != 0)) bad++;
			if (mask[i] != 0.0f) fg++;
		}
		printf("argmax: %zu/%zu voxels differ (ours fg %zu)\n", bad, MG_NVOX, fg);
		if (bad != 0) failures++;
		free(aref);
	}
	size_t dim[3] = {MG_DIM, MG_DIM, MG_DIM};
	bwlabel(mask, 26, dim, true, false);
	size_t fg = 0;
	for (size_t i = 0; i < MG_NVOX; i++)
		if (mask[i] != 0.0f) fg++;
	printf("largest component: %zu voxels\n", fg);
	printf("elapsed %.2f s\n", now() - t0);

	char out[1024];
	snprintf(out, sizeof(out), "%s/mask_c.u8", dir);
	FILE *f = fopen(out, "wb");
	if (f) {
		unsigned char *m8 = (unsigned char *)malloc(MG_NVOX);
		if (m8) {
			for (size_t i = 0; i < MG_NVOX; i++)
				m8[i] = mask[i] != 0.0f;
			if (fwrite(m8, 1, MG_NVOX, f) != MG_NVOX) failures++;
		} else {
			failures++;
		}
		fclose(f);
		free(m8);
	}
	for (int i = 0; i < MINDGRAB_NHIDDEN; i++) free(wpack[i]);
	free(mask);
	free(partial);
	free(a);
	free(b);
	free(conformed);
	if (failures) fprintf(stderr, "mindgrab parity: %d FAILURE(S)\n", failures);
	return failures ? 1 : 0;
}

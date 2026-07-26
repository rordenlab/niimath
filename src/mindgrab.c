// mindgrab.c - MindGrab skull stripping (-mindgrab)
//
// Clean-room C implementation of the inference brainchop-cli performs with tinygrad:
// quantile normalisation, a 25-block dilated MeshNet (Fedorov et al. 2017), a 2-class
// argmax and the largest 26-connected component. Only the trained weights are upstream
// (MIT, see src/mindgrab.LICENSE and scripts/export_mindgrab.py); no runtime is linked.
//
// The reference is brainchop-cli 0.1.24 with tinygrad; it was used as a black-box oracle,
// layer by layer. Conventions that silently corrupt the result if they are guessed rather
// than measured, each verified against that oracle:
//
//   * SPATIAL AXES ARE REVERSED. brainchop feeds the network volume.transpose((2,1,0)),
//     so the tensor's (D,H,W) are the NIfTI (x,y,z) and a torch kernel tap (kd,kh,kw) is
//     the (dx,dy,dz) offset. scripts/export_mindgrab.py bakes the transpose into the
//     emitted table; this file consumes taps in x-fastest order.
//   * GELU IS THE TANH APPROXIMATION, not exact erf. Exact erf differs in the 4th decimal
//     -- enough to move decision-boundary voxels.
//   * GroupNorm has num_groups == num_channels == 15 and affine=False, i.e. per-channel
//     normalisation over the whole 256^3 volume with eps 1e-5, no scale and no shift.
//     tinygrad's layernorm normalises by mean((x-mean)^2), the BIASED variance.
//   * Hidden convolutions have NO bias, and the reference does not apply the checkpoint's
//     final-classifier bias (details at mg_classify()).
//   * Dilation cycles 16, 8, 4, 2, 1 five times with padding == dilation, so every layer
//     preserves 256^3.
//   * qnormalize's quantiles are numpy's default linear interpolation on the sorted
//     array. The input is uint8, so a 256-bin histogram reproduces them exactly -- but
//     only with numpy's lerp branch (t < 0.5 ? a+(b-a)t : b-(b-a)(1-t)) and float32
//     rounding of the result, which NEP 50 keeps float32 throughout.
//   * argmax ties go to class 0 (numpy returns the FIRST maximum). No tie was observed on
//     real data; the reference's minimum |logit0-logit1| was 3.0e-5.
//
// Accuracy budget: the reference's own METAL and CPU backends differ by up to 3.1e-5 in
// the logits yet agree on every one of the 16.7 M argmax voxels, the smallest margin
// being 3.0e-5. That is the tolerance this port has to live inside, and what
// test/mindgrab_parity.c measures.

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mindgrab.h"
#include "mindgrab_weights.h"
#include "bwlabel.h"
#include "print.h"

// Grid size. test/mindgrab_selftest.c overrides it so the real kernels can be exercised on a
// 32^3 grid (2 MB instead of 2 GB); nothing in a normal build ever defines it.
#ifndef MG_DIM
#define MG_DIM MINDGRAB_DIM
#endif
#define MG_NVOX ((size_t)MG_DIM * MG_DIM * MG_DIM)
// mg_conv_impl's x-segmentation assumes the largest dilation fits twice across the grid: segment
// 0 keeps the +dil tap for every x in [0, dil), which is only in range when dil <= MG_DIM - dil.
// The largest dilation is 16, so MG_DIM must be >= 32 -- production 256 and the selftest's 32
// (exactly tight) both qualify, but the seam above would otherwise let a smaller override read
// and write past the row with segment 1 inverted and hiding it. C99 negative-array-size assert.
typedef char mg_dim_is_large_enough[(MG_DIM >= 32) ? 1 : -1];
#define MG_C MINDGRAB_CHAN
// Channel stride is padded 15 -> 16 so every activation vector is exactly 64 bytes wide and the
// inner accumulator is exactly four NEON/SSE registers. (The buffers come from plain malloc, so
// the alignment is malloc's 16 bytes -- the win is the shape, not the alignment.) Lane 15 stays
// zero: the padded weight column is zero, so the convolution never writes it and GroupNorm maps
// 0 to 0.
#define MG_CS 16
#define MG_EPS 1e-5

// tanh-GELU is evaluated 6.3e9 times per volume, so libm tanhf would dominate the runtime.
// Odd rational x*P(x^2)/Q(x^2), fitted by scripts/fit_mindgrab_tanh.py; <= 5 float32 ULP
// over the whole range, and tanhf itself rounds to +-1 beyond |x| = 9.011.
#define MG_TANH_LIM 9.0f
#define MG_TANH_P1 0.99999994f
#define MG_TANH_P3 0.130791619f
#define MG_TANH_P5 0.00309865153f
#define MG_TANH_P7 1.11017944e-05f
#define MG_TANH_P9 -2.00294448e-08f
#define MG_TANH_P11 5.19263868e-11f
#define MG_TANH_P13 -8.29118133e-14f
#define MG_TANH_Q0 1.0f
#define MG_TANH_Q2 0.464124829f
#define MG_TANH_Q4 0.024473751f
#define MG_TANH_Q6 0.000253859733f

static const int mg_dilation[MINDGRAB_NHIDDEN] = {
	16, 8, 4, 2, 1, 16, 8, 4, 2, 1, 16, 8, 4, 2, 1, 16, 8, 4, 2, 1, 16, 8, 4, 2, 1};

static inline float mg_tanhf(float x) {
	if (x > MG_TANH_LIM) x = MG_TANH_LIM;
	if (x < -MG_TANH_LIM) x = -MG_TANH_LIM;
	float u = x * x;
	float p = MG_TANH_P13;
	p = p * u + MG_TANH_P11;
	p = p * u + MG_TANH_P9;
	p = p * u + MG_TANH_P7;
	p = p * u + MG_TANH_P5;
	p = p * u + MG_TANH_P3;
	p = p * u + MG_TANH_P1;
	float q = MG_TANH_Q6;
	q = q * u + MG_TANH_Q4;
	q = q * u + MG_TANH_Q2;
	q = q * u + MG_TANH_Q0;
	return (x * p) / q;
}

// 0.5x(1 + tanh(sqrt(2/pi)(x + 0.044715 x^3))), grouped as tinygrad groups it. static: the two
// test harnesses #include this .c and so reach it directly; nothing else outside needs it.
static float mindgrab_gelu(float x) {
	float u = 0.797884583f * (x + 0.044715f * (x * x * x));
	return (0.5f * x) * (1.0f + mg_tanhf(u));
}

// ---------------------------------------------------------------- quantile normalisation

// numpy's default "linear" quantile of the conformed uint8 volume, from a 256-bin
// histogram: exact, and O(nvox) instead of a 16.7 M element sort.
static double mg_quantile(const int64_t *cum, size_t n, double q) {
	double vi = q * (double)(n - 1);
	double fi = floor(vi);
	double t = vi - fi;
	size_t i = (size_t)fi;
	int a = 0, b = 0;
	// smallest value whose cumulative count exceeds the rank
	while (a < 255 && cum[a] <= (int64_t)i) a++;
	b = a;
	while (b < 255 && cum[b] <= (int64_t)(i + 1)) b++;
	if (t < 0.5) return (double)a + (double)(b - a) * t;
	return (double)b - (double)(b - a) * (1.0 - t);
}

static void mg_qnormalize(const unsigned char *vol, float *dst) {
	int64_t hist[256], cum[256];
	memset(hist, 0, sizeof(hist));
	for (size_t i = 0; i < MG_NVOX; i++)
		hist[vol[i]]++;
	int64_t run = 0;
	for (int i = 0; i < 256; i++) {
		run += hist[i];
		cum[i] = run;
	}
	float qlo = (float)mg_quantile(cum, MG_NVOX, 0.02);
	float qhi = (float)mg_quantile(cum, MG_NVOX, 0.98);
	float den = (qhi - qlo) + 1e-3f;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int z = 0; z < MG_DIM; z++) {
		size_t lo = (size_t)z * MG_DIM * MG_DIM;
		for (size_t i = lo; i < lo + (size_t)MG_DIM * MG_DIM; i++) {
			float v = ((float)vol[i] - qlo) / den;
			if (!(v > 0.0f)) v = 0.0f;
			if (v > 1.0f) v = 1.0f;
			dst[i] = v;
		}
	}
}

// ---------------------------------------------------------------------------- the network

// Per-layer padded weights: [27 taps][cin][MG_CS], tap = tz*9 + ty*3 + tx (tx fastest),
// out-channel contiguous. Column MG_C..MG_CS-1 is zero padding.
static float *mg_pack_weights(int layer, int *cin_out) {
	int cin = (layer == 0) ? 1 : MG_C;
	size_t n = (size_t)27 * cin * MG_CS;
	float *w = (float *)calloc(n, sizeof(float));
	if (!w) return NULL;
	const float *src = mindgrab_weights + MINDGRAB_LOFF(layer);
	for (int t = 0; t < 27; t++)
		for (int ic = 0; ic < cin; ic++)
			for (int oc = 0; oc < MG_C; oc++)
				w[((size_t)t * cin + ic) * MG_CS + oc] = src[((size_t)t * cin + ic) * MG_C + oc];
	*cin_out = cin;
	return w;
}

// One dilated 3x3x3 convolution, src channel-last with stride `sstride`, dst with MG_CS.
//
// Two things here are load-bearing for speed, both measured on an M4 Pro at dilation 1:
//
//   * `restrict` and a COMPILE-TIME cin. Without them clang cannot prove dst does not
//     alias the weights, keeps nothing in vector registers, and runs at 27 GFLOP/s.
//   * MG_NB output voxels per inner block. With one voxel the loop issues one 128-bit
//     weight load per FMA and is load-port bound; sharing each loaded weight vector across
//     four x-neighbours makes it 1 load : 4 FMA. 27 -> 37 -> 110 GFLOP/s per core for the
//     three steps. MG_NB=8 is SLOWER (32 accumulator registers exceed the ARM64 file and
//     everything spills again) -- do not raise it without re-measuring.
//
// The three x segments exist so the inner loop carries no bounds test: with
// padding == dilation, x in [dil, 255-dil] has all three x-taps in range, and the two
// margins drop one each.
#define MG_NB 4

static inline void mg_conv_impl(const float *restrict src, const int sstride, const int cin,
								float *restrict dst, const float *restrict w, const int dil) {
	const size_t wlayer = (size_t)cin * MG_CS;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int z = 0; z < MG_DIM; z++) {
		for (int y = 0; y < MG_DIM; y++) {
			ptrdiff_t rowoff[9];
			int rowtap[9];
			int nrow = 0;
			for (int tz = 0; tz < 3; tz++) {
				int zz = z + (tz - 1) * dil;
				if (zz < 0 || zz >= MG_DIM) continue;
				for (int ty = 0; ty < 3; ty++) {
					int yy = y + (ty - 1) * dil;
					if (yy < 0 || yy >= MG_DIM) continue;
					rowoff[nrow] = ((ptrdiff_t)zz * MG_DIM + yy) * MG_DIM;
					rowtap[nrow] = (tz * 3 + ty) * 3;
					nrow++;
				}
			}
			size_t dstrow = ((size_t)z * MG_DIM + y) * MG_DIM;
			int xseg[4] = {0, dil, MG_DIM - dil, MG_DIM};
			for (int seg = 0; seg < 3; seg++) {
				int txlo = (seg == 0) ? 1 : 0;
				int txhi = (seg == 2) ? 1 : 2;
				ptrdiff_t tvox[27];
				const float *tw[27];
				int ntap = 0;
				for (int r = 0; r < nrow; r++)
					for (int tx = txlo; tx <= txhi; tx++) {
						tvox[ntap] = rowoff[r] + (tx - 1) * dil;
						tw[ntap] = w + (size_t)(rowtap[r] + tx) * wlayer;
						ntap++;
					}
				int x = xseg[seg];
				const int xend = xseg[seg + 1];
				for (; x + MG_NB <= xend; x += MG_NB) {
					float acc[MG_NB][MG_CS];
					for (int k = 0; k < MG_NB; k++)
						for (int c = 0; c < MG_CS; c++)
							acc[k][c] = 0.0f;
					for (int t = 0; t < ntap; t++) {
						const float *xin = src + (size_t)(tvox[t] + x) * sstride;
						const float *ww = tw[t];
						for (int ic = 0; ic < cin; ic++) {
							const float *wv = ww + (size_t)ic * MG_CS;
							float s[MG_NB];
							for (int k = 0; k < MG_NB; k++)
								s[k] = xin[k * sstride + ic];
							for (int c = 0; c < MG_CS; c++) {
								float wo = wv[c];
								for (int k = 0; k < MG_NB; k++)
									acc[k][c] += s[k] * wo;
							}
						}
					}
					float *o = dst + (dstrow + x) * MG_CS;
					for (int k = 0; k < MG_NB; k++)
						for (int c = 0; c < MG_CS; c++)
							o[k * MG_CS + c] = acc[k][c];
				}
				for (; x < xend; x++) { // margins can be 1 or 2 voxels wide
					float acc[MG_CS];
					for (int c = 0; c < MG_CS; c++)
						acc[c] = 0.0f;
					for (int t = 0; t < ntap; t++) {
						const float *xin = src + (size_t)(tvox[t] + x) * sstride;
						const float *ww = tw[t];
						for (int ic = 0; ic < cin; ic++) {
							float xi = xin[ic];
							for (int oc = 0; oc < MG_CS; oc++)
								acc[oc] += xi * ww[(size_t)ic * MG_CS + oc];
						}
					}
					float *o = dst + (dstrow + x) * MG_CS;
					for (int c = 0; c < MG_CS; c++)
						o[c] = acc[c];
				}
			}
		}
	}
}

// Two specialisations so cin and the source stride are literal constants inside the kernel.
static void mg_conv_first(const float *restrict src, float *restrict dst,
						  const float *restrict w, int dil) {
	mg_conv_impl(src, 1, 1, dst, w, dil);
}

static void mg_conv(const float *restrict src, float *restrict dst,
					const float *restrict w, int dil) {
	mg_conv_impl(src, MG_CS, MG_C, dst, w, dil);
}

// GroupNorm(num_groups == num_channels, affine=False) followed by tanh-GELU, in place.
// The two moments are accumulated per z-slice in double and combined in slice order, so
// the result does not depend on the thread count.
static void mg_norm_gelu(float *buf, double *partial) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int z = 0; z < MG_DIM; z++) {
		double s[MG_C], q[MG_C];
		for (int c = 0; c < MG_C; c++) {
			s[c] = 0.0;
			q[c] = 0.0;
		}
		const float *p = buf + (size_t)z * MG_DIM * MG_DIM * MG_CS;
		for (int i = 0; i < MG_DIM * MG_DIM; i++) {
			for (int c = 0; c < MG_C; c++) {
				double v = p[(size_t)i * MG_CS + c];
				s[c] += v;
				q[c] += v * v;
			}
		}
		for (int c = 0; c < MG_C; c++) {
			partial[(size_t)z * 2 * MG_C + c] = s[c];
			partial[(size_t)z * 2 * MG_C + MG_C + c] = q[c];
		}
	}
	double mean[MG_C], scale[MG_C];
	for (int c = 0; c < MG_C; c++) {
		double s = 0.0, q = 0.0;
		for (int z = 0; z < MG_DIM; z++) {
			s += partial[(size_t)z * 2 * MG_C + c];
			q += partial[(size_t)z * 2 * MG_C + MG_C + c];
		}
		mean[c] = s / (double)MG_NVOX;
		double var = q / (double)MG_NVOX - mean[c] * mean[c];
		if (!(var > 0.0)) var = 0.0;
		scale[c] = 1.0 / sqrt(var + MG_EPS);
	}
	float fmean[MG_C], fscale[MG_C];
	for (int c = 0; c < MG_C; c++) {
		fmean[c] = (float)mean[c];
		fscale[c] = (float)scale[c];
	}
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int z = 0; z < MG_DIM; z++) {
		float *p = buf + (size_t)z * MG_DIM * MG_DIM * MG_CS;
		for (int i = 0; i < MG_DIM * MG_DIM; i++) {
			float *v = p + (size_t)i * MG_CS;
			for (int c = 0; c < MG_C; c++)
				v[c] = mindgrab_gelu((v[c] - fmean[c]) * fscale[c]);
		}
	}
}

// 1x1x1 classifier + argmax. Ties go to class 0, matching numpy's argmax.
//
// THE CLASSIFIER BIAS IS DELIBERATELY NOT APPLIED. model.pth carries model.25.bias, but the
// reference never loads it: brainchop builds its tinygrad Conv2d from model.json, which has
// no "bias" key, so bias defaults to False -- and tiny_meshnet.convert_keys then maps the
// checkpoint onto the model by ZIPPING the two key lists, which silently drops the 27th
// (bias) entry. The shipped browser export (net_mindgrab.safetensors) likewise contains 26
// weight tensors and no bias, so this is the model as published, not an accident of the CLI.
// Applying the bias moves the logits by up to 1.53 and changes thousands of voxels. The two
// values stay in mindgrab_weights[] (MINDGRAB_BIAS_OFF) so the table remains a faithful copy
// of the checkpoint.
static void mg_classify(const float *src, float *mask) {
	const float *w = mindgrab_weights + MINDGRAB_CLS_OFF;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int z = 0; z < MG_DIM; z++) {
		size_t lo = (size_t)z * MG_DIM * MG_DIM;
		for (size_t i = lo; i < lo + (size_t)MG_DIM * MG_DIM; i++) {
			const float *x = src + i * MG_CS;
			float l0 = 0.0f, l1 = 0.0f;
			for (int c = 0; c < MG_C; c++) {
				l0 += x[c] * w[c * 2];
				l1 += x[c] * w[c * 2 + 1];
			}
			mask[i] = (l1 > l0) ? 1.0f : 0.0f;
		}
	}
}

int mindgrab_segment(const unsigned char *conformed, float *mask) {
	if (!conformed || !mask) return 1;
	float *wpack[MINDGRAB_NHIDDEN];
	int wcin; // scratch: mg_pack_weights reports cin, which the caller does not need
	int nw;
	float *a = NULL, *b = NULL;
	double *partial = NULL;
	int rc = 1;
	for (nw = 0; nw < MINDGRAB_NHIDDEN; nw++) {
		wpack[nw] = mg_pack_weights(nw, &wcin);
		if (!wpack[nw]) goto done;
	}
	a = (float *)malloc(MG_NVOX * MG_CS * sizeof(float));
	b = (float *)malloc(MG_NVOX * MG_CS * sizeof(float));
	partial = (double *)malloc((size_t)MG_DIM * 2 * MG_C * sizeof(double));
	if (!a || !b || !partial) {
		printfx("mindgrab: out of memory (needs about %.1f GB)\n",
				2.0 * (double)MG_NVOX * MG_CS * sizeof(float) / 1e9);
		goto done;
	}
	// The normalised input is one channel, so it fits in the head of b. Layer 0 READS b and
	// writes a; b is not overwritten until layer 1, by which time it has been consumed. No
	// third buffer.
	mg_qnormalize(conformed, b);
	mg_conv_first(b, a, wpack[0], mg_dilation[0]);
	mg_norm_gelu(a, partial);
	for (int l = 1; l < MINDGRAB_NHIDDEN; l++) {
		float *src = (l & 1) ? a : b;
		float *dst = (l & 1) ? b : a;
		mg_conv(src, dst, wpack[l], mg_dilation[l]);
		mg_norm_gelu(dst, partial);
	}
	// layer 0 wrote a; layer l wrote (l odd ? b : a), so the last hidden layer left its
	// output in the buffer selected by (MINDGRAB_NHIDDEN - 1)
	mg_classify(((MINDGRAB_NHIDDEN - 1) & 1) ? b : a, mask);
	free(a);
	free(b);
	a = b = NULL;
	// keep only the largest 26-connected cluster, as the reference does with -bwlabel 26
	{
		size_t dim[3] = {MG_DIM, MG_DIM, MG_DIM};
		bwlabel(mask, 26, dim, true, false);
	}
	// An empty mask is not an allocation failure, so it cannot be reported through rc without
	// changing the exit status of a run that "worked" -- but the caller will then blank every
	// voxel to the image minimum and hand back an all-background image. Say so, loudly: silence
	// here reads as success and the output looks like a legitimate strip of a very dark scan.
	{
		size_t fg = 0;
		for (size_t i = 0; i < MG_NVOX; i++)
			if (mask[i] != 0.0f) fg++;
		if (fg == 0)
			printfx("mindgrab: WARNING - the network found no brain; the output will be blank. "
					"Check that the input is a head image with usable contrast.\n");
	}
	rc = 0;
done:
	for (int i = 0; i < nw; i++)
		free(wpack[i]);
	free(a);
	free(b);
	free(partial);
	return rc;
}

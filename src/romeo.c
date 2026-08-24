/*----------------------------------------------------------------------------------------------
 * romeo.c — ROMEO phase unwrapping (`-romeo`) for niimath.
 *
 * PROVENANCE.  Faithful C port of:
 *   ROMEO.jl            v1.4.0, git 60d83fbb69669560d227c252dcd844afbb6648e1  (MIT)
 *                       Korbinian Eckstein, Barbara Dymerska, Simon Robinson
 *                       https://github.com/korbinian90/ROMEO.jl
 *   MriResearchTools.jl v3.5.0                                              (MIT)
 *                       Korbinian Eckstein
 *                       https://github.com/korbinian90/MriResearchTools.jl
 *                       (robustmask, gaussiansmooth3d box filter, sample/quantile,
 *                        readphase rescaling)
 * plus the Float64 2*pi range reduction from the Julia standard library
 *   base/math.jl + base/special/rem_pio2.jl                                 (MIT)
 *   which is itself derived from FDLIBM __ieee754_rem_pio2 (Sun Microsystems, 1993).
 * Upstream copyright and permission notices are preserved verbatim in src/romeo.LICENSE.
 *
 * Algorithm citation: Dymerska, B., Eckstein, K., Bachrata, B., Siow, B., Trattnig, S.,
 * Shmueli, K., Robinson, S.D., 2020. "Phase Unwrapping with a Rapid Opensource Minimum
 * Spanning TreE AlgOrithm (ROMEO)." Magn Reson Med. doi:10.1002/mrm.28563
 *
 * ---------------------------------------------------------------------------------------------
 * WHY THIS FILE IS COMPILED STRICT-FP (-fno-fast-math -ffp-contract=off)
 *
 * ROMEO is a minimum-spanning-tree region grow whose traversal order is decided by 8-bit
 * integer edge weights produced by `rescale(w) = max(round(Int,(1-w)*255), 1)`.  A single
 * reassociated or FMA-contracted expression can move one weight from bin 137 to bin 138,
 * change the spanning tree, and shift an entire connected region by exactly 2*pi.  The rest of
 * niimath is built whole-program -ffast-math; this translation unit must not be.
 *
 * That is not a theoretical worry - it was MEASURED, building this same source three ways and
 * comparing against the pinned Julia oracle. The middle column is the 76x76x46 validation volume
 * (phase0/mag0, -t 16.8); the right column is the FULL parity suite (test/romeo_compare.py (medic_bench repo)
 * --weights-all: 4 real + 11 synthetic cases, 10 weight selections).
 *
 * NOTE: the pass counts below are AS MEASURED WHEN THIS EXPERIMENT WAS RUN, when the suite had
 * 422 checks. It has since grown to 602 (the six B0 weighting modes are now compared byte-exactly).
 * The counts are left at their measured values rather than rewritten, because they record an
 * experiment, not a current claim -- any re-measurement of the FP policy must re-run the suite and
 * restate them:
 *
 *   FP policy for romeo.o              e0 weight bytes differing   FULL parity suite
 *   ---------------------------------  -------------------------   ---------------------------
 *   -fno-fast-math -ffp-contract=off         0 / 797088            422/422 pass   (SHIPPED)
 *   -fno-fast-math -ffp-contract=fast        0 / 797088            9 FAIL: the readphase rescale
 *     (FMA contraction only)                                       on the line_x/plane_xy
 *                                                                  fixtures differs by 1 float
 *                                                                  ULP, and the pre-rescale
 *                                                                  weights drift past the 4 ULP
 *                                                                  limit (me: 5 ULP)
 *   -ffast-math -fno-finite-math-only      360 / 797088            66 voxels off by >=1 FULL
 *     (the repository-wide default)                                2*pi wrap, max|diff| 12.57 rad
 *
 * The -ffast-math failure mode is NOT float32 rounding: reassociation pushes a weight just past
 * 1.0, the `0 <= w <= 1` guard in rescale() then returns bin 0, and the edge DISAPPEARS from the
 * graph (largest observed bin deviation: 252). Multi-echo output stops matching even at
 * --compare 1e-4.
 *
 * FMA contraction is bit-identical on the validation volume but NOT on the full corpus (row 2) -
 * an earlier version of this comment claimed otherwise on the strength of the single-volume
 * measurement alone. It is also unsafe in principle: it would silently break the Dekker/
 * Cody-Waite double-double arithmetic in the 2*pi range reduction below, where Julia fuses ONLY
 * at its explicit muladd sites (mirrored here as explicit fma() calls). And it buys nothing
 * measurable - the unwrap is 0.02 s on the validation volume (0.07 s wall including gzip
 * output). -fno-fast-math does not inhibit SIMD auto-vectorisation of the non-reduction loops,
 * so no vectorisation is given up; only reassociation is. Change ROMEO_STRICT_FP in
 * src/Makefile to revisit - and re-run the FULL suite, not one volume.
 *
 * PLATFORM CAVEAT on the "bit-identical" claim: on Linux/x86 with gcc, the whole-program
 * -ffast-math on the LINK line pulls in crtfastmath.o, which sets MXCSR FTZ/DAZ process-wide -
 * including inside this strict-FP object. ROMEO's operands (weights in [0,1], phase ~[-pi,pi],
 * magnitudes ~1e3) never reach the denormal range, so the exposure is nil, but every measurement
 * above was made on arm64 macOS. Re-run test/romeo_compare.py (medic_bench repo) on a Linux gcc release build
 * before claiming bit-identity there. Do NOT "fix" this by removing -ffast-math from the link.
 *
 * ---------------------------------------------------------------------------------------------
 * NUMERIC TYPE AUDIT.  Julia's promotion rules are NOT uniform across the
 * six weight terms; each row was confirmed with typeof(...) under the pinned environment above.
 * The C column is what this file implements.  Getting one row wrong is the single most likely
 * cause of a near-miss.
 *
 *   Julia expression                                          Julia type   C type here
 *   --------------------------------------------------------  -----------  ---------------------
 *   gamma(x::Float32)             (single-wrap fold)           Float32      float   rm_gamma_f
 *   gamma(x::Float64)                                          Float64      double  rm_gamma_d
 *   phasecoherence  = 1-abs(g/pi)  (Irrational pi -> Float32)  Float32      float
 *   phasegradientcoherence         (TEs are Float64)           Float64      double
 *   phaselinearity(P,i,j,k) normal                             Float32      float  (rm_pl.f64=0)
 *   phaselinearity(P,i,j,k) isnan -> 0.5                       Float64      double (rm_pl.f64=1)
 *   phaselinearity(P,i,j)  interior product of two Float32     Float32      float
 *   phaselinearity(P,i,j)  boundary fallback 0.9               Float64      double
 *   phaselinearity(P,i,j)  interior, one operand widened       Float64      double
 *   magcoherence = (small/big)^2                               Float32      float
 *   magweight, magweight2   (maxmag is Float64)                Float64      double
 *   getweight accumulator   (weight = 1.0)                     Float64      double
 *   0.1 + 0.9x  factors     (Float32 x widened first)          Float64      double
 *   maxmag = quantile(mag[isfinite], 0.95)                     Float64      double
 *   robustmask quantiles q05/q15/q8/q99                        Float64      double
 *   robustmask high_intensity/noise/threshold (mean of F32)    Float32      float
 *   unwrapvoxel(new,old) = new - 2pi*round((new-old)/2pi)      Float64      double -> stored f32
 *   rem2pi(x::Float32, RoundNearest)                           Float32      (float)rem2pi(f64)
 *   box filter running sum (boxfilterline!)                    Float32      float
 *
 * DELIBERATE DEVIATION, one place only: Julia's `mean(::Vector{Float32})` reduces through
 * `sum`, i.e. pairwise blocks of 1024 whose sequential base case is `@simd`, so LLVM reduces it
 * with 8 vector accumulators on the oracle machine.  That order is a property of the oracle's
 * codegen, not a portable contract (an x86 Julia would differ).  robustmask's two means are
 * therefore accumulated here in double and rounded once to float: more accurate than the
 * oracle, and within the plan's <=1e-6 relative scalar tolerance (observed <=1 float ULP,
 * ~9e-8 relative).  The mask stages themselves are still required to match exactly; a threshold
 * this close only reorders voxels that sit inside a 1-ULP band, and none do on the corpus.
 *
 * ROUNDING.  Julia's round(Int,x) is round-half-to-EVEN (0.5->0, 1.5->2, 2.5->2).  C's
 * round()/lround() are half-away-from-zero and MUST NOT be used: nearbyint()/rint() under the
 * default FE_TONEAREST are the correct spelling.  This governs rescale() (the weight bin) and
 * unwrapvoxel() (the number of 2*pi wraps) — both load-bearing.
 *
 * INDEXING.  Julia's checkbounds(Bool, A, i) on a LINEAR index only tests 1 <= i <= length(A),
 * so a neighbour lookup may cross a row or plane and still be "in bounds".  ROMEO relies on
 * that, then zeroes weights[1,end,:,:] / [2,:,end,:] / [3,:,:,end] so the invalid directed
 * edges never enter the queue.  The same linear semantics appear in phaselinearity's h/k and in
 * unwrapedge's `oo`.  The formulas below are ported literally; do NOT "fix" them into Cartesian
 * bounds checks.  Voxel indices are kept 1-BASED inside the algorithm to match the reference.
 *--------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include "romeo.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define RM_NBINS 256
#define RM_PI_F32 3.14159274101257324f    /* Float32(pi)  = 0x1.921fb6p+1 */
#define RM_2PI_F32 6.28318548202514648f   /* Float32(2pi) = 0x1.921fb6p+2 */
#define RM_PI_F64 3.14159265358979323846  /* Float64(pi) */
#define RM_2PI_F64 6.28318530717958647692 /* Float64(2pi) */

/* MSVC has no strtok_r; its strtok_s takes the identical 3 arguments. romeo.c IS compiled on
   Windows (src/CMakeLists.txt adds it to ADDITIONAL_SRCS, and both AppVeyor and
   release-binaries.yml build with MSVC), so this is a hard build break without the alias. */
#if defined(_MSC_VER) && !defined(strtok_r)
#define strtok_r strtok_s
#endif

#define RM_ERR(...) do { fprintf(stderr, "** -romeo: "); fprintf(stderr, __VA_ARGS__); } while (0)

/* ============================================================================================
 * 1. Julia-compatible numeric primitives
 * ==========================================================================================*/

/* Julia's max/min propagate NaN; C's fmax/fmin return the non-NaN operand. */
static float rm_maxf(float a, float b) { if (isnan(a)) return a; if (isnan(b)) return b; return (b > a) ? b : a; }
static double rm_maxd(double a, double b) { if (isnan(a)) return a; if (isnan(b)) return b; return (b > a) ? b : a; }
static double rm_mind(double a, double b) { if (isnan(a)) return a; if (isnan(b)) return b; return (b < a) ? b : a; }

/* gamma: fold by at most ONE wrap.  The Julia comparison is against the exact irrational pi;
   no float or double lies strictly between Float64(pi) and the true pi, so comparing the
   widened operand against RM_PI_F64 reproduces `x > pi` / `x < -pi` exactly. */
static float rm_gamma_f(float x) {
	if ((double)x < -RM_PI_F64) return x + RM_2PI_F32;
	if ((double)x > RM_PI_F64) return x - RM_2PI_F32;
	return x;
}
static double rm_gamma_d(double x) {
	if (x < -RM_PI_F64) return x + RM_2PI_F64;  /* typeof(x)(2pi) == Float64(2pi) */
	if (x > RM_PI_F64) return x - RM_2PI_F64;
	return x;
}

/* -------------------------------------------------------------------------------------------
 * rem2pi(x, RoundNearest): port of Julia base/math.jl + base/special/rem_pio2.jl.
 * Julia uses an INFINITELY precise 2*pi, so libc remainder(x, 6.283185307179586) is NOT
 * equivalent (it differs by ~n*2.4e-16).  Every muladd below is a real fma(); every other
 * operation must not be contracted, which is why this TU is built -ffp-contract=off.
 * -----------------------------------------------------------------------------------------*/

typedef struct { double hi, lo; } rm_dd;

static double rm_add22condh(double xh, double xl, double yh, double yl) {
	double r = xh + yh;
	double s = (fabs(xh) > fabs(yh)) ? ((((xh - r) + yh) + yl) + xl) : ((((yh - r) + xh) + xl) + yl);
	return r + s;
}

static uint32_t rm_highword(double x) {
	uint64_t u;
	memcpy(&u, &x, sizeof u);
	return (uint32_t)(u >> 32);
}

static int rm_cw2c(double x, double fn, int n, rm_dd *y) {
	const double pio2_1 = 1.57079632673412561417e+00;
	const double pio2_1t = 6.07710050650619224932e-11;
	double z = fma(-fn, pio2_1, x);
	double y1 = fma(-fn, pio2_1t, z);
	double y2 = fma(-fn, pio2_1t, (z - y1));
	y->hi = y1; y->lo = y2;
	return n;
}

static int rm_cwext(double x, uint32_t xhp, rm_dd *y) {
	const double pio2_1 = 1.57079632673412561417e+00;
	const double pio2_1t = 6.07710050650619224932e-11;
	const double pio2_2 = 6.07710050630396597660e-11;
	const double pio2_2t = 2.02226624879595063154e-21;
	const double pio2_3 = 2.02226624871116645580e-21;
	const double pio2_3t = 8.47842766036889956997e-32;
	double fn = nearbyint(x * (2.0 / RM_PI_F64));
	double r = fma(-fn, pio2_1, x);
	double w = fn * pio2_1t;
	int32_t j = (int32_t)(xhp >> 20);
	double y1 = r - w;
	uint32_t high = rm_highword(y1);
	int32_t i = j - (int32_t)((high >> 20) & 0x7ff);
	if (i > 16) {
		double t = r;
		w = fn * pio2_2;
		r = t - w;
		w = fma(fn, pio2_2t, -((t - r) - w));
		y1 = r - w;
		high = rm_highword(y1);
		i = j - (int32_t)((high >> 20) & 0x7ff);
		if (i > 49) {
			t = r;
			w = fn * pio2_3;
			r = t - w;
			w = fma(fn, pio2_3t, -((t - r) - w));
			y1 = r - w;
		}
	}
	y->hi = y1;
	y->lo = (r - y1) - w;
	return (int)(int64_t)fn; /* unsafe_trunc(Int, fn) */
}

/* --- 128-bit helpers for Payne-Hanek.  Written portably (no __int128) so every platform,
   including MSVC and wasm, produces the identical reduction. --------------------------------*/
typedef struct { uint64_t hi, lo; } rm_u128;

static rm_u128 rm_u128_add(rm_u128 a, rm_u128 b) {
	rm_u128 r;
	r.lo = a.lo + b.lo;
	r.hi = a.hi + b.hi + (r.lo < a.lo ? 1u : 0u);
	return r;
}
static rm_u128 rm_u128_mul(uint64_t a, uint64_t b) { /* full 128-bit product */
	uint64_t a0 = a & 0xffffffffu, a1 = a >> 32;
	uint64_t b0 = b & 0xffffffffu, b1 = b >> 32;
	uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
	uint64_t mid = (p00 >> 32) + (p01 & 0xffffffffu) + (p10 & 0xffffffffu);
	rm_u128 r;
	r.lo = (p00 & 0xffffffffu) | (mid << 32);
	r.hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
	return r;
}
static rm_u128 rm_u128_shl(rm_u128 a, int s) {
	rm_u128 r;
	if (s == 0) return a;
	if (s >= 128) { r.hi = 0; r.lo = 0; return r; }
	if (s >= 64) { r.hi = a.lo << (s - 64); r.lo = 0; return r; }
	r.hi = (a.hi << s) | (a.lo >> (64 - s));
	r.lo = a.lo << s;
	return r;
}
static rm_u128 rm_u128_shr(rm_u128 a, int s) { /* logical */
	rm_u128 r;
	if (s == 0) return a;
	if (s >= 128) { r.hi = 0; r.lo = 0; return r; }
	if (s >= 64) { r.lo = a.hi >> (s - 64); r.hi = 0; return r; }
	r.lo = (a.lo >> s) | (a.hi << (64 - s));
	r.hi = a.hi >> s;
	return r;
}
static rm_u128 rm_u128_neg(rm_u128 a) {
	rm_u128 r;
	r.lo = ~a.lo + 1u;
	r.hi = ~a.hi + (r.lo == 0 ? 1u : 0u);
	return r;
}
static int rm_u128_is_zero(rm_u128 a) { return a.hi == 0 && a.lo == 0; }
static int rm_u128_top_set_bit(rm_u128 a) { /* 1-based position of the highest set bit */
	int n = 0;
	uint64_t v;
	if (a.hi) { n = 64; v = a.hi; } else { v = a.lo; }
	while (v) { v >>= 1; n++; }
	return n;
}
static rm_u128 rm_u128_sub(rm_u128 a, rm_u128 b) { return rm_u128_add(a, rm_u128_neg(b)); }

/* fromfraction(f::Int128) — `f` is passed as its two's-complement bit pattern. */
static void rm_fromfraction(rm_u128 f, double *z1o, double *z2o) {
	uint64_t s;
	rm_u128 x, m1shift, x2;
	int n1, n2;
	uint64_t m1, d1, m2, d2, bits;
	double z1, z2;
	if (rm_u128_is_zero(f)) { *z1o = 0.0; *z2o = 0.0; return; }
	s = (f.hi >> 63) ? (UINT64_C(1) << 63) : 0;
	x = (f.hi >> 63) ? rm_u128_neg(f) : f;                 /* abs(f) */
	n1 = rm_u128_top_set_bit(x);
	{
		rm_u128 t = (n1 >= 26) ? rm_u128_shr(x, n1 - 26) : rm_u128_shl(x, 26 - n1);
		m1 = t.lo << 27;
	}
	d1 = ((uint64_t)(int64_t)(n1 - 128 + 1021)) << 52;
	bits = s | (d1 + m1);
	memcpy(&z1, &bits, sizeof z1);
	{
		rm_u128 mm; mm.hi = 0; mm.lo = m1;
		m1shift = (n1 >= 53) ? rm_u128_shl(mm, n1 - 53) : rm_u128_shr(mm, 53 - n1);
	}
	x2 = rm_u128_sub(x, m1shift);
	if (rm_u128_is_zero(x2)) { *z1o = z1; *z2o = 0.0; return; }
	n2 = rm_u128_top_set_bit(x2);
	{
		rm_u128 t = (n2 >= 53) ? rm_u128_shr(x2, n2 - 53) : rm_u128_shl(x2, 53 - n2);
		m2 = t.lo;
	}
	d2 = ((uint64_t)(int64_t)(n2 - 128 + 1021)) << 52;
	bits = s | (d2 + m2);
	memcpy(&z2, &bits, sizeof z2);
	*z1o = z1; *z2o = z2;
}

static const uint64_t RM_INV_2PI[19] = {
	UINT64_C(0x28be60db9391054a), UINT64_C(0x7f09d5f47d4d3770), UINT64_C(0x36d8a5664f10e410),
	UINT64_C(0x7f9458eaf7aef158), UINT64_C(0x6dc91b8e909374b8), UINT64_C(0x01924bba82746487),
	UINT64_C(0x3f877ac72c4a69cf), UINT64_C(0xba208d7d4baed121), UINT64_C(0x3a671c09ad17df90),
	UINT64_C(0x4e64758e60d4ce7d), UINT64_C(0x272117e2ef7e4a0e), UINT64_C(0xc7fe25fff7816603),
	UINT64_C(0xfbcbc462d6829b47), UINT64_C(0xdb4d9fb3c9f2c26d), UINT64_C(0xd3d18fd9a797fa8b),
	UINT64_C(0x5d49eeb1faf97c5e), UINT64_C(0xcf41ce7de294a4ba), UINT64_C(0x9afed7ec47e35742),
	UINT64_C(0x1580cc11bf1edaea)
};

static int rm_paynehanek(double x, rm_dd *y) {
	uint64_t u, X;
	int raw_exponent, k, idx, shift, q;
	uint64_t a1, a2, a3;
	rm_u128 w1, w2, w3, w, f;
	double z_hi, z_lo, y_hi, y_lo;
	const double pio2 = 1.5707963267948966;
	const double pio2_hi = 1.5707963407039642;
	const double pio2_lo = -1.3909067614167116e-8;
	memcpy(&u, &x, sizeof u);
	X = (u & UINT64_C(0x000fffffffffffff)) | (UINT64_C(1) << 52);
	raw_exponent = (int)((u >> 52) & 0x7ff);
	k = raw_exponent - 1023 - 52;
	idx = k >> 6;                 /* arithmetic shift, matches Julia's k >> 6 */
	/* Julia writes `k - (idx << 6)`, but idx is NEGATIVE here for every |x| just above the
	   Payne-Hanek threshold (k ~ -32 -> idx == -1), and left-shifting a negative int is UB in C
	   (caught by UBSan). idx * 64 is the same value and well defined; idx stays in [-1, 15]. */
	shift = k - idx * 64;
	if (shift == 0) {
		a1 = (idx + 0 < 0) ? 0 : RM_INV_2PI[idx + 0];
		a2 = RM_INV_2PI[idx + 1];
		a3 = RM_INV_2PI[idx + 2];
	} else {
		a1 = ((idx < 0) ? 0 : (RM_INV_2PI[idx] << shift)) | (RM_INV_2PI[idx + 1] >> (64 - shift));
		a2 = (RM_INV_2PI[idx + 1] << shift) | (RM_INV_2PI[idx + 2] >> (64 - shift));
		a3 = (RM_INV_2PI[idx + 2] << shift) | (RM_INV_2PI[idx + 3] >> (64 - shift));
	}
	w1.hi = X * a1; w1.lo = 0;                 /* UInt128(X*a1) << 64 (overflow -> integer part) */
	w2 = rm_u128_mul(X, a2);
	w3 = rm_u128_shr(rm_u128_mul(X, a3), 64);
	w = rm_u128_add(rm_u128_add(w1, w2), w3);
	if (x < 0.0) w = rm_u128_neg(w);           /* flipsign(w, x) */
	q = (int)((((int64_t)(rm_u128_shr(w, 125).lo)) + 1) >> 1);
	f = rm_u128_shl(w, 2);
	rm_fromfraction(f, &z_hi, &z_lo);
	y_hi = (z_hi + z_lo) * pio2;
	y_lo = (((z_hi * pio2_hi - y_hi) + z_hi * pio2_lo) + z_lo * pio2_hi) + z_lo * pio2_lo;
	y->hi = y_hi; y->lo = y_lo;
	return q;
}

static int rm_rem_pio2_kernel(double x, rm_dd *y) {
	uint32_t xhp = rm_highword(x) & 0x7fffffffu;
	if (xhp <= 0x400f6a7au) {
		if ((xhp & 0xfffffu) == 0x921fbu) return rm_cwext(x, xhp, y);
		if (xhp <= 0x4002d97cu) return (x > 0.0) ? rm_cw2c(x, 1.0, 1, y) : rm_cw2c(x, -1.0, -1, y);
		return (x > 0.0) ? rm_cw2c(x, 2.0, 2, y) : rm_cw2c(x, -2.0, -2, y);
	}
	if (xhp <= 0x401c463bu) {
		if (xhp <= 0x4015fdbcu) {
			if (xhp == 0x4012d97cu) return rm_cwext(x, xhp, y);
			return (x > 0.0) ? rm_cw2c(x, 3.0, 3, y) : rm_cw2c(x, -3.0, -3, y);
		}
		if (xhp == 0x401921fbu) return rm_cwext(x, xhp, y);
		return (x > 0.0) ? rm_cw2c(x, 4.0, 4, y) : rm_cw2c(x, -4.0, -4, y);
	}
	if (xhp < 0x413921fbu) return rm_cwext(x, xhp, y);
	return rm_paynehanek(x, y);
}

static double rm_rem2pi_d(double x) {
	rm_dd y;
	int n;
	const double pi1o2_h = 1.5707963267948966, pi1o2_l = 6.123233995736766e-17;
	const double pi2o2_h = 3.141592653589793, pi2o2_l = 1.2246467991473532e-16;
	if (isnan(x)) return x;
	if (isinf(x)) return NAN;
	/* Julia: `abs(x) < pi` against the exact irrational.  Float64(pi) itself IS less than the
	   true pi, so the C spelling is <=, not <. */
	if (fabs(x) <= RM_PI_F64) return x;
	n = rm_rem_pio2_kernel(x, &y);
	if ((n & 1) == 0) {
		if ((n & 2) == 2) {
			if (y.hi <= 0) return rm_add22condh(y.hi, y.lo, pi2o2_h, pi2o2_l);
			return rm_add22condh(y.hi, y.lo, -pi2o2_h, -pi2o2_l);
		}
		return y.hi + y.lo;
	}
	if ((n & 2) == 2) return rm_add22condh(y.hi, y.lo, -pi1o2_h, -pi1o2_l);
	return rm_add22condh(y.hi, y.lo, pi1o2_h, pi1o2_l);
}

/* Base: rem2pi(x::Float32, r) = Float32(rem2pi(Float64(x), r)) */
static float rm_rem2pi_f(float x) { return (float)rm_rem2pi_d((double)x); }

/* unwrapvoxel(new, old) = new - 2pi*round((new-old)/2pi).  `new` is always Float32 (an element
   of the working array); `old` is Float32 in the spatial path and Float64 in the temporal path
   (refvalue = wrapped[...] .* (TE_i/TE_ref) widens to Float64), so the subtraction happens in a
   different width in each case.  Both variants are needed. */
static float rm_unwrapvoxel_ff(float nw, float od) {
	float d = nw - od;                       /* Float32 subtraction */
	double q = nearbyint((double)d / RM_2PI_F64);
	return (float)((double)nw - RM_2PI_F64 * q);
}
static float rm_unwrapvoxel_fd(float nw, double od) {
	double d = (double)nw - od;              /* Float64 subtraction */
	double q = nearbyint(d / RM_2PI_F64);
	return (float)((double)nw - RM_2PI_F64 * q);
}

/* from: 1 best, 0 worst.  to: 1 best, NBINS-1 worst, 0 = invalid (never queued). */
static uint8_t rm_rescale(double w) {
	if (0.0 <= w && w <= 1.0) {
		double r = nearbyint((1.0 - w) * (double)(RM_NBINS - 1));
		int64_t ri = (int64_t)r;
		return (uint8_t)(ri < 1 ? 1 : ri);
	}
	return 0;
}

/* ============================================================================================
 * 2. Order statistics: comparator-free quickselect (the WASM qsort prohibition applies),
 *    Julia's type-7 quantile, and Base.middle-based median.
 * ==========================================================================================*/

/* Select the k-th smallest (0-based) of v[0..n), permuting v.  Median-of-three pivot with a
   deterministic fallback; no function-pointer comparator. */
static float rm_select_kth_f(float *v, int64_t n, int64_t k) {
	int64_t lo = 0, hi = n - 1;
	while (lo < hi) {
		int64_t i = lo, j = hi, mid = lo + ((hi - lo) >> 1);
		float a = v[lo], b = v[mid], c = v[hi], p;
		p = (a < b) ? ((b < c) ? b : ((a < c) ? c : a)) : ((a < c) ? a : ((b < c) ? c : b));
		while (i <= j) {
			while (v[i] < p) i++;
			while (v[j] > p) j--;
			if (i <= j) { float t = v[i]; v[i] = v[j]; v[j] = t; i++; j--; }
		}
		if (k <= j) hi = j;
		else if (k >= i) lo = i;
		else return v[k];
	}
	return v[lo];
}
static double rm_select_kth_d(double *v, int64_t n, int64_t k) {
	int64_t lo = 0, hi = n - 1;
	while (lo < hi) {
		int64_t i = lo, j = hi, mid = lo + ((hi - lo) >> 1);
		double a = v[lo], b = v[mid], c = v[hi], p;
		p = (a < b) ? ((b < c) ? b : ((a < c) ? c : a)) : ((a < c) ? a : ((b < c) ? c : b));
		while (i <= j) {
			while (v[i] < p) i++;
			while (v[j] > p) j--;
			if (i <= j) { double t = v[i]; v[i] = v[j]; v[j] = t; i++; j--; }
		}
		if (k <= j) hi = j;
		else if (k >= i) lo = i;
		else return v[k];
	}
	return v[lo];
}

/* isapprox(a::Float32, b::Float32) with Julia's defaults: atol=0, rtol=sqrt(eps(Float32)). */
static int rm_isapprox_f(float a, float b) {
	const float rtol = 3.4526698e-4f; /* sqrt(eps(Float32)) */
	float ma = fabsf(a), mb = fabsf(b);
	float m = (mb > ma) ? mb : ma;
	return fabsf(a - b) <= rtol * m;
}

/* Statistics.quantile(v::Vector{Float32}, p) — type 7 (alpha=beta=1), returns Float64.
   Destroys the ordering of `v`.  n must be > 0.
   PINNED-VERSION GOTCHA: the pinned oracle environment (see test/romeo_oracle.jl in the medic_bench repo)
   resolves the REGISTRY package
   Statistics v1.11.1, not the copy bundled with Julia's own stdlib tree.  1.11.1 computes
   `aleph = n*p + m` with a plain multiply-add; the newer stdlib copy uses `fma(n, p, m)`.  The
   two differ in the last bits (for n=265696, p=0.95: 252411.24999999997 vs 252411.25), which
   moves gamma just off 1/4 and changes maxmag in the 11th digit.  Match the pinned package. */
static double rm_quantile7(float *v, int64_t n, double p) {
	double m, aleph, gam;
	int64_t j;
	float a, b;
	m = 1.0 + p * (1.0 - 1.0 - 1.0);            /* alpha + p*(one-alpha-beta) with alpha=beta=1 */
	aleph = (double)n * p + m;
	j = (int64_t)aleph;                          /* trunc toward zero; aleph >= 0 here */
	if (j < 1) j = 1;
	if (j > n - 1) j = n - 1;
	if (j < 1) j = 1;                            /* n == 1 */
	gam = aleph - (double)j;
	if (gam < 0.0) gam = 0.0;
	if (gam > 1.0) gam = 1.0;
	if (n == 1) { a = v[0]; b = v[0]; }
	else {
		a = rm_select_kth_f(v, n, j - 1);
		b = rm_select_kth_f(v, n, j);            /* v is partially ordered; still correct */
	}
	if (isfinite(a) && isfinite(b) && rm_isapprox_f(a, b))
		return (double)a + gam * (double)(b - a);
	return (1.0 - gam) * (double)a + gam * (double)b;
}

/* Statistics.median: sorted middle, Base.middle(x,y) = x/2 + y/2 for even counts. */
static double rm_median_d(double *v, int64_t n) {
	if (n == 1) return v[0];
	if (n % 2 == 1) return rm_select_kth_d(v, n, n / 2);
	{
		double a = rm_select_kth_d(v, n, n / 2 - 1);
		double b = rm_select_kth_d(v, n, n / 2);
		return a / 2.0 + b / 2.0;
	}
}

/* ============================================================================================
 * 3. MriResearchTools.sample / approxextrema
 *
 * sample(I; n=1e5) takes len=ceil(sqrt(min(n,length))) blocks of len contiguous elements whose
 * starts are round.(Int, range(0, length-len; length=len)), then keeps the finite ones (falling
 * back to the whole array when that is empty).  The block starts are computed here with EXACT
 * integer arithmetic and round-half-to-even, which reproduces Julia's TwicePrecision range
 * without inheriting its floating-point spelling.
 * ==========================================================================================*/

static int64_t rm_sample_len(int64_t n) {
	double nn = (1.0e5 < (double)n) ? 1.0e5 : (double)n;
	double s = sqrt(nn);
	int64_t len = (int64_t)ceil(s);
	if (len < 1) len = 1;
	if (len > n) len = n;
	return len;
}

/* round(Int, (i-1)*(stop)/(len-1)) with stop = n-len, computed exactly (half-to-even). */
static int64_t rm_range_round(int64_t i /*0-based*/, int64_t stop, int64_t len) {
	int64_t num, den, q, r;
	if (len <= 1) return 0;
	num = i * stop;
	den = len - 1;
	q = num / den;
	r = num - q * den;
	if (2 * r > den) return q + 1;
	if (2 * r == den) return (q % 2 == 0) ? q : q + 1;
	return q;
}

/* Gather the sample of v[0..n) into out (capacity >= len*len). Returns 0 when every sampled
   value was non-finite; the caller then runs the whole-array fallback into a larger buffer. */
static int64_t rm_sample_f32(const float *v, int64_t n, float *out) {
	int64_t len = rm_sample_len(n), stop = n - len, b, t, m = 0;
	if (stop < 0) stop = 0;
	for (b = 0; b < len; b++) {
		int64_t s = rm_range_round(b, stop, len);
		for (t = 1; t <= len; t++) {
			int64_t idx = s + t;              /* 1-based Julia index */
			if (idx >= 1 && idx <= n) {
				float x = v[idx - 1];
				if (isfinite(x)) out[m++] = x;
			}
		}
	}
	return m;
}

/* sample(), allocating only what is needed: len*len (~100k) normally, the full array only for
   the fallback. *buf is malloc'd and owned by the caller. Returns the count, or -1 on failure. */
static int64_t rm_sample_alloc(const float *v, int64_t n, float **buf) {
	int64_t len = rm_sample_len(n), cap = len * len, m;
	float *s = (float *)malloc((size_t)cap * sizeof(float));
	if (!s) return -1;
	m = rm_sample_f32(v, n, s);
	if (m == 0) {   /* filter(isfinite, I) over the WHOLE array */
		float *g = (float *)realloc(s, (size_t)n * sizeof(float));
		if (!g) { free(s); return -1; }
		s = g;
		for (int64_t i = 0; i < n; i++) if (isfinite(v[i])) s[m++] = v[i];
	}
	*buf = s;
	return m;
}

/* approxextrema(I) = extrema(sample(I)), falling back to extrema(I) when the sample is flat. */
static int rm_approxextrema(const float *v, int64_t n, float *mn, float *mx) {
	int64_t m, i;
	float *s = NULL;
	float lo, hi;
	m = rm_sample_alloc(v, n, &s);
	if (m < 1) { free(s); return 1; }
	lo = s[0]; hi = s[0];
	for (i = 1; i < m; i++) { if (s[i] < lo) lo = s[i]; if (s[i] > hi) hi = s[i]; }
	free(s);
	if (lo == hi) {
		lo = v[0]; hi = v[0];
		for (i = 1; i < n; i++) { if (v[i] < lo) lo = v[i]; if (v[i] > hi) hi = v[i]; }
	}
	*mn = lo; *mx = hi;
	return 0;
}

/* ============================================================================================
 * 4. Weights (ROMEO.jl src/weights.jl)
 * ==========================================================================================*/

typedef struct { double v; int f64; } rm_pl; /* value plus "was widened to Float64" */

typedef struct {
	const float *P;        /* wrapped phase, template echo (1-based access via P[i-1]) */
	const float *P2;       /* phase2 (Float32) or NULL */
	const double *P2d;     /* phase2 (Float64, temporal-uncertain path) or NULL */
	double TE1, TE2;
	const float *M;        /* magnitude .* mask, or NULL */
	double maxmag;
	const uint8_t *mask;   /* NULL = all true */
	int flags[6];
	int nx, ny, nz;
	int64_t n;             /* nx*ny*nz */
} rm_wctx;

/* phaselinearity(P,i,j,k) — 1-based indices. */
static rm_pl rm_pl3(const float *P, int64_t i, int64_t j, int64_t k) {
	rm_pl out;
	float t = (P[i - 1] - 2.0f * P[j - 1]) + P[k - 1];
	float r = rm_rem2pi_f(t);
	float pl = rm_maxf(0.0f, 1.0f - fabsf(r / 2.0f));
	if (isnan(pl)) { out.v = 0.5; out.f64 = 1; return out; }
	out.v = (double)pl; out.f64 = 0;
	return out;
}

/* phaselinearity(P,i,j) — the linear h/k formula is deliberate (see INDEXING note). */
static rm_pl rm_pl2(const float *P, int64_t i, int64_t j, int64_t n) {
	rm_pl out, a, b;
	int64_t neighbor = j - i, h = i - neighbor, k = j + neighbor;
	if (0 < h && k <= n) {
		a = rm_pl3(P, h, i, j);
		b = rm_pl3(P, i, j, k);
		if (!a.f64 && !b.f64) { out.v = (double)((float)a.v * (float)b.v); out.f64 = 0; }
		else { out.v = a.v * b.v; out.f64 = 1; }
		return out;
	}
	out.v = 0.9; out.f64 = 1;
	return out;
}

static double rm_getweight(const rm_wctx *c, int64_t i, int64_t j) {
	double weight = 1.0;
	if (c->flags[0]) {
		float pc = 1.0f - fabsf(rm_gamma_f(c->P[i - 1] - c->P[j - 1]) / RM_PI_F32);
		weight *= (0.1 + 0.9 * (double)pc);
	}
	if (c->flags[1]) {
		double g1 = (double)rm_gamma_f(c->P[i - 1] - c->P[j - 1]);
		double g2;
		if (c->P2d) g2 = rm_gamma_d(c->P2d[i - 1] - c->P2d[j - 1]);
		else g2 = (double)rm_gamma_f(c->P2[i - 1] - c->P2[j - 1]);
		weight *= (0.1 + 0.9 * rm_maxd(0.0, 1.0 - fabs(g1 - g2 * c->TE1 / c->TE2)));
	}
	if (c->flags[2]) {
		rm_pl pl = rm_pl2(c->P, i, j, c->n);
		weight *= (0.1 + 0.9 * pl.v);
	}
	if (c->M) {
		float mi = c->M[i - 1], mj = c->M[j - 1], small, big;
		/* Base.minmax propagates NaN to BOTH endpoints — `minmax(5f0, NaN32) == (NaN32, NaN32)`,
		   verified in the pinned Julia. A plain `mj < mi` swap leaves one endpoint finite, which
		   changes magweight/magweight2 (flags 5/6) and can retain an edge Julia drops. */
		if (isnan(mi) || isnan(mj)) { small = big = mi + mj; }
		else if (mj < mi) { small = mj; big = mi; }
		else { small = mi; big = mj; }
		if (c->flags[3]) {
			float q = small / big;
			weight *= (0.1 + 0.9 * (double)(q * q));
		}
		if (c->flags[4])
			weight *= (0.1 + 0.9 * (0.5 + 0.5 * rm_mind(1.0, (double)small / (0.5 * c->maxmag))));
		if (c->flags[5])
			weight *= (0.1 + 0.9 * (0.5 + 0.5 * rm_mind(1.0, (0.5 * c->maxmag) / (double)big)));
	}
	return weight;
}

enum { RM_WOUT_U8 = 0, RM_WOUT_F32, RM_WOUT_F64 };

/* calculateweights_romeo: weights[3, nx, ny, nz], then zero the three non-existent planes. */
static void rm_calculateweights(const rm_wctx *c, int outmode, void *out) {
	int64_t n = c->n;
	int64_t stride[3];
	int dim;
	stride[0] = 1; stride[1] = c->nx; stride[2] = (int64_t)c->nx * c->ny;
	memset(out, 0, (size_t)3 * (size_t)n * (outmode == RM_WOUT_U8 ? 1 : (outmode == RM_WOUT_F32 ? 4 : 8)));
	for (dim = 0; dim < 3; dim++) {
		int64_t nb = stride[dim];
		int64_t i0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (i0 = 0; i0 < n; i0++) {
			int64_t I = i0 + 1, J = I + nb;
			double w;
			if (c->mask && !c->mask[i0]) continue;
			if (J > n) continue;
			w = rm_getweight(c, I, J);
			if (outmode == RM_WOUT_U8) ((uint8_t *)out)[3 * i0 + dim] = rm_rescale(w);
			else if (outmode == RM_WOUT_F32) ((float *)out)[3 * i0 + dim] = (float)w;
			else ((double *)out)[3 * i0 + dim] = w;
		}
	}
	/* these edges do not exist */
	{
		int64_t x, y, z;
		int64_t nx = c->nx, ny = c->ny, nz = c->nz;
		for (z = 0; z < nz; z++) for (y = 0; y < ny; y++) {
			int64_t i0 = (nx - 1) + nx * (y + ny * z);
			if (outmode == RM_WOUT_U8) ((uint8_t *)out)[3 * i0 + 0] = 0;
			else if (outmode == RM_WOUT_F32) ((float *)out)[3 * i0 + 0] = 0.0f;
			else ((double *)out)[3 * i0 + 0] = 0.0;
		}
		for (z = 0; z < nz; z++) for (x = 0; x < nx; x++) {
			int64_t i0 = x + nx * ((ny - 1) + ny * z);
			if (outmode == RM_WOUT_U8) ((uint8_t *)out)[3 * i0 + 1] = 0;
			else if (outmode == RM_WOUT_F32) ((float *)out)[3 * i0 + 1] = 0.0f;
			else ((double *)out)[3 * i0 + 1] = 0.0;
		}
		for (y = 0; y < ny; y++) for (x = 0; x < nx; x++) {
			int64_t i0 = x + nx * (y + ny * (nz - 1));
			if (outmode == RM_WOUT_U8) ((uint8_t *)out)[3 * i0 + 2] = 0;
			else if (outmode == RM_WOUT_F32) ((float *)out)[3 * i0 + 2] = 0.0f;
			else ((double *)out)[3 * i0 + 2] = 0.0;
		}
	}
}

/* updateflags: no magnitude disables 4..6; no phase2/TEs pair disables 2. */
static void rm_updateflags(int *flags, int have_p2, int have_TEs, int have_mag) {
	if (!have_mag) { flags[3] = 0; flags[4] = 0; flags[5] = 0; }
	if (!have_p2 || !have_TEs) flags[1] = 0;
}

static void rm_flags_from_sel(int sel, int have_mag, int *flags) {
	int i;
	for (i = 0; i < 6; i++) flags[i] = 0;
	switch (sel) {
	case RM_W_ROMEO:  /* resolved by the app: romeo3 with magnitude, romeo4 without */
		if (have_mag) { flags[0] = flags[1] = flags[3] = 1; }
		else { flags[0] = flags[1] = flags[2] = flags[3] = 1; }
		break;
	case RM_W_ROMEO2: flags[0] = flags[3] = 1; break;
	case RM_W_ROMEO3: flags[0] = flags[1] = flags[3] = 1; break;
	case RM_W_ROMEO4: flags[0] = flags[1] = flags[2] = flags[3] = 1; break;
	case RM_W_ROMEO6: for (i = 0; i < 6; i++) flags[i] = 1; break;
	default: break;
	}
}

/* ============================================================================================
 * 5. voxelquality (ROMEO.jl src/voxelquality.jl)
 * ==========================================================================================*/

static int rm_voxelquality(const rm_wctx *c, float *qmap) {
	int64_t n = c->n, i, x, y, z;
	int64_t nx = c->nx, ny = c->ny, nz = c->nz;
	float *w = (float *)malloc((size_t)3 * (size_t)n * sizeof(float));
	if (!w) return 1;
	rm_calculateweights(c, RM_WOUT_F32, w);
	for (i = 0; i < n; i++) qmap[i] = (w[3 * i] + w[3 * i + 1]) + w[3 * i + 2];
	/* qmap[2:end,:,:] .+= weights[1,1:end-1,:,:] (and the two other directions), in this order */
	for (z = 0; z < nz; z++) for (y = 0; y < ny; y++) for (x = 1; x < nx; x++) {
		int64_t i0 = x + nx * (y + ny * z);
		qmap[i0] += w[3 * (i0 - 1) + 0];
	}
	for (z = 0; z < nz; z++) for (y = 1; y < ny; y++) for (x = 0; x < nx; x++) {
		int64_t i0 = x + nx * (y + ny * z);
		qmap[i0] += w[3 * (i0 - nx) + 1];
	}
	for (z = 1; z < nz; z++) for (y = 0; y < ny; y++) for (x = 0; x < nx; x++) {
		int64_t i0 = x + nx * (y + ny * z);
		qmap[i0] += w[3 * (i0 - nx * ny) + 2];
	}
	for (i = 0; i < n; i++) qmap[i] = qmap[i] / 6.0f;
	free(w);
	return 0;
}

/* ============================================================================================
 * 5b. B0 field map (MriResearchTools romeofunctions.jl: calculateB0_unwrapped / get_B0_snr)
 *
 *   B0  = ((1000/2pi) * sum(phase./TEs .* w; dims=4)) ./ sum(w; dims=4),  non-finite -> 0
 *   snr = sum(mag .* w; dims=4) ./ sum(w; dims=4)
 *
 * Note the association: upstream writes `(1000 / 2pi) * sum(...) ./ sum(...)`, and `*` and `./`
 * are the same precedence and left-associative, so the constant multiplies the NUMERATOR before
 * the division — not the quotient afterwards.
 *
 * WIDTHS: the weight is Float64 for every mode EXCEPT `mag` with a real (Float32) magnitude,
 * where `mag .* weight` and both sums stay Float32.  Verified per mode against the pinned Julia.
 * Without a magnitude the app substitutes a Float64 exp(-TE/20) decay, which also makes the
 * `mag` mode Float64 - hence the `w32` flag depends on BOTH the mode and have_mag.
 * ==========================================================================================*/

/* `m32` is 1 when the magnitude is a real Float32 image: upstream's `mag .* mag .* TEs .* TEs`
   rounds the SQUARE in Float32 and only then widens. Rounding it in double drifts the B0 map by
   ~1.5e-5 Hz and the SNR by ~2e-3. The magnitude-free fallback is Float64 throughout. */
static double rm_b0_weight_d(int mode, double m, double te, int m32, double sim) {
	switch (mode) {
	case RM_B0_PHASE_VAR: return (m32 ? (double)((float)m * (float)m) : m * m) * te * te;
	case RM_B0_AVERAGE: return 1.0;
	case RM_B0_TES: return te;
	case RM_B0_MAG: return m;
	case RM_B0_SIMULATED_MAG: return sim * te;
	default: return m * te;   /* RM_B0_PHASE_SNR */
	}
}

/* `mag` may be NULL: the app then uses a voxel-independent exp(-TE/20) T2* decay. */
static void rm_compute_b0(const float *phase, const float *mag, int64_t n3, int neco,
	const double *TEs, int mode, float *b0, float *snr) {
	const int w32 = (mode == RM_B0_MAG) && (mag != NULL);
	double synth[ROMEO_MAX_TE];   /* exp(-TE/20): echo-only, so hoist it out of the voxel loop */
	int64_t i;
	int e;
	for (e = 0; e < neco; e++) synth[e] = exp(-TEs[e] / 20.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(e)
#endif
	for (i = 0; i < n3; i++) {
		double num = 0.0, den = 0.0, snum = 0.0;   /* num is Float64 in every mode */
		float denf = 0.0f, snumf = 0.0f;
		for (e = 0; e < neco; e++) {
			double te = TEs[e];
			double m = mag ? (double)mag[(int64_t)e * n3 + i] : synth[e];
			double p = (double)phase[(int64_t)e * n3 + i];
			if (w32) {
				float w = (float)m;                /* weight == the Float32 magnitude itself */
				denf += w;
				snumf += (float)m * w;
				num += p / te * (double)w;
			} else {
				double w = rm_b0_weight_d(mode, m, te, mag != NULL, synth[e]);
				den += w;
				snum += m * w;
				num += p / te * w;
			}
		}
		{
			double d = w32 ? (double)denf : den;
			double v = (1000.0 / RM_2PI_F64) * num / d;
			b0[i] = (fabs(v) <= DBL_MAX) ? (float)v : 0.0f;   /* B0[.!isfinite.(B0)] .= 0 */
			snr[i] = w32 ? (snumf / denf) : (float)(snum / den);
		}
	}
}

/* ============================================================================================
 * 6. robustmask (MriResearchTools masking.jl + smoothing.jl)
 * ==========================================================================================*/

/* boxfilterline!(line, boxsize, q) — running average with the asymmetric edge normalisation
   `line[i] = lsum/(r+i)` at the leading edge.  The running sum is Float32 and the middle-part
   update is `lsum += (line[i+r] - popfirst!(q))`, i.e. the difference is formed FIRST; that is
   not the same rounding as `(lsum + line[i+r]) - qold`. */
static void rm_boxfilterline(float *line, int64_t len, int boxsize, float *q) {
	int r = boxsize / 2;
	int64_t i;
	int qh = 0, qn = 0; /* ring buffer head + count, capacity boxsize */
	float lsum = 0.0f;
	for (i = 0; i < r; i++) { lsum += line[i]; q[(qh + qn) % boxsize] = line[i]; qn++; }
	for (i = 1; i <= r + 1; i++) {
		float v = line[i + r - 1];
		lsum += v;
		q[(qh + qn) % boxsize] = v; qn++;
		line[i - 1] = lsum / (float)(r + i);
	}
	for (i = r + 2; i <= len - r; i++) {
		float v = line[i + r - 1];
		float old = q[qh]; qh = (qh + 1) % boxsize; qn--;
		lsum += (v - old);
		q[(qh + qn) % boxsize] = v; qn++;
		line[i - 1] = lsum / (float)boxsize;
	}
	for (i = len - r + 1; i <= len; i++) {
		float old = q[qh]; qh = (qh + 1) % boxsize; qn--;
		lsum -= old;
		line[i - 1] = lsum / (float)(r + len - i + 1);
	}
}

/* checkboxsizes!: force odd, and clamp to sz/2 (also forced odd). */
static int rm_checkboxsize(int bs, int64_t sz) {
	if (bs % 2 == 0) bs += 1;
	if ((double)bs > (double)sz / 2.0) {
		int64_t val = sz / 2;
		if (val % 2 == 0) val += 1;
		bs = (int)val;
	}
	return bs;
}

/* gaussiansmooth3d(image; nbox, boxsizes) with mask=nothing, weight=nothing: nbox passes of the
   running-average box filter over dims 1..3.  This is NOT a Gaussian convolution; substituting
   niimath's nifti_smooth_gauss changes the mask. */
static int rm_boxsmooth3d(float *img, int nx, int ny, int nz, int nbox, const int *boxes) {
	int64_t sz[3];
	int bs[3][4];
	int ibox, dim, d, k;
	int maxbs = 1;
	float *line = NULL, *q = NULL;
	sz[0] = nx; sz[1] = ny; sz[2] = nz;
	if (nbox > 4) return 1;
	for (d = 0; d < 3; d++) for (k = 0; k < nbox; k++) {
		bs[d][k] = rm_checkboxsize(boxes[k], sz[d]);
		if (bs[d][k] > maxbs) maxbs = bs[d][k];
	}
	{
		int64_t maxlen = sz[0] > sz[1] ? sz[0] : sz[1];
		if (sz[2] > maxlen) maxlen = sz[2];
		line = (float *)malloc((size_t)maxlen * sizeof(float));
		q = (float *)malloc((size_t)maxbs * sizeof(float));
		if (!line || !q) { free(line); free(q); return 1; }
	}
	for (ibox = 0; ibox < nbox; ibox++) for (dim = 0; dim < 3; dim++) {
		int boxsize = bs[dim][ibox];
		int64_t len = sz[dim], stride, outer0, outer1, a, b, t;
		if (len == 1 || boxsize < 3) continue;
		stride = (dim == 0) ? 1 : ((dim == 1) ? nx : (int64_t)nx * ny);
		outer0 = (dim == 0) ? ny : ((dim == 1) ? nx : nx);
		outer1 = (dim == 0) ? nz : ((dim == 1) ? nz : ny);
		for (b = 0; b < outer1; b++) for (a = 0; a < outer0; a++) {
			int64_t base;
			if (dim == 0) base = nx * (a + ny * b);
			else if (dim == 1) base = a + (int64_t)nx * ny * b;
			else base = a + nx * b;
			for (t = 0; t < len; t++) line[t] = img[base + t * stride];
			rm_boxfilterline(line, len, boxsize, q);
			for (t = 0; t < len; t++) img[base + t * stride] = line[t];
		}
	}
	free(line); free(q);
	return 0;
}

/* fill_holes(mask) = .!imfill(.!mask, (1, length/20)) with 6-connectivity: any connected
   component of the COMPLEMENT whose size is in [1, length/20] becomes part of the mask. */
static int rm_fill_holes(uint8_t *mask, int nx, int ny, int nz) {
	int64_t n = (int64_t)nx * ny * nz, i;
	double maxhole = (double)n / 20.0;
	int32_t *lab = NULL;
	int32_t *stack = NULL, *sizes = NULL;   /* n <= INT_MAX, so 32-bit halves this scratch */
	int32_t nlab = 0, sp;
	if (1.0 > maxhole) {
		RM_ERR("robustmask needs at least 20 voxels (upstream fill_holes passes (1, n/20) to imfill, which rejects n < 20)\n");
		return 1;
	}
	lab = (int32_t *)calloc((size_t)n, sizeof(int32_t));
	stack = (int32_t *)malloc((size_t)n * sizeof(int32_t));
	sizes = (int32_t *)malloc((size_t)(n + 1) * sizeof(int32_t));
	if (!lab || !stack || !sizes) { free(lab); free(stack); free(sizes); return 1; }
	for (i = 0; i < n; i++) {
		if (mask[i] || lab[i]) continue;
		nlab++;
		sizes[nlab] = 0;
		sp = 0; stack[sp++] = (int32_t)i; lab[i] = nlab;
		while (sp > 0) {
			int64_t v = stack[--sp];
			int64_t z = v / ((int64_t)nx * ny), rem = v % ((int64_t)nx * ny);
			int64_t y = rem / nx, x = rem % nx;
			sizes[nlab]++;
			if (x > 0 && !mask[v - 1] && !lab[v - 1]) { lab[v - 1] = nlab; stack[sp++] = (int32_t)(v - 1); }
			if (x < nx - 1 && !mask[v + 1] && !lab[v + 1]) { lab[v + 1] = nlab; stack[sp++] = (int32_t)(v + 1); }
			if (y > 0 && !mask[v - nx] && !lab[v - nx]) { lab[v - nx] = nlab; stack[sp++] = (int32_t)(v - nx); }
			if (y < ny - 1 && !mask[v + nx] && !lab[v + nx]) { lab[v + nx] = nlab; stack[sp++] = (int32_t)(v + nx); }
			if (z > 0 && !mask[v - (int64_t)nx * ny] && !lab[v - (int64_t)nx * ny]) { lab[v - (int64_t)nx * ny] = nlab; stack[sp++] = (int32_t)(v - (int64_t)nx * ny); }
			if (z < nz - 1 && !mask[v + (int64_t)nx * ny] && !lab[v + (int64_t)nx * ny]) { lab[v + (int64_t)nx * ny] = nlab; stack[sp++] = (int32_t)(v + (int64_t)nx * ny); }
		}
	}
	for (i = 0; i < n; i++) {
		if (!mask[i]) {
			int32_t c = sizes[lab[i]];
			if (1 <= c && (double)c <= maxhole) mask[i] = 1; /* hole filled */
		}
	}
	free(lab); free(stack); free(sizes);
	return 0;
}

typedef struct {           /* observable robustmask intermediates, for the parity dump */
	uint8_t *s1, *s2, *s3, *s4;
	float *sm1, *sm2;
	double q05, q15, q8, q99;
	float high_intensity, noise, threshold;
	int noise_stage;
	int64_t sample_len;
} rm_mask_stages;

/* Frees EVERY stage including s4. The two callers that keep the final mask detach it first
   (`mask = stages.s4; stages.s4 = NULL;`), so this stays symmetric and every error path is a
   single call rather than a free-plus-free-s4 pair. */
static void rm_mask_stages_free(rm_mask_stages *s) {
	free(s->s1); free(s->s2); free(s->s3); free(s->sm1); free(s->sm2); free(s->s4);
	s->s1 = s->s2 = s->s3 = s->s4 = NULL; s->sm1 = s->sm2 = NULL;
}

/* robustmask(weight; factor=1, threshold=nothing).  Returns 0 and fills stages->s4 (owned by
   the caller) on success.  `have_thr` selects the -k qualitymask path. */
#define RM_DROP(p) do { if (!keep_stages) { free(p); (p) = NULL; } } while (0)
static int rm_robustmask(const float *weight, int nx, int ny, int nz,
	int have_thr, double thr_in, int keep_stages, rm_mask_stages *st) {
	int64_t n = (int64_t)nx * ny * nz, i, m;
	float *s = NULL;
	float threshold;
	memset(st, 0, sizeof *st);
	if (!have_thr) {
		double q05, q15, q8, q99, acc;
		int64_t cnt;
		float *tmp = NULL;
		m = rm_sample_alloc(weight, n, &s);
		if (m < 1) { free(s); RM_ERR("magnitude has no finite voxels\n"); return 1; }
		tmp = (float *)malloc((size_t)m * sizeof(float));
		if (!tmp) { free(s); return 1; }
		st->sample_len = m;
		memcpy(tmp, s, (size_t)m * sizeof(float)); q05 = rm_quantile7(tmp, m, 0.05);
		memcpy(tmp, s, (size_t)m * sizeof(float)); q15 = rm_quantile7(tmp, m, 0.15);
		memcpy(tmp, s, (size_t)m * sizeof(float)); q8 = rm_quantile7(tmp, m, 0.80);
		memcpy(tmp, s, (size_t)m * sizeof(float)); q99 = rm_quantile7(tmp, m, 0.99);
		free(tmp);
		st->q05 = q05; st->q15 = q15; st->q8 = q8; st->q99 = q99;
		acc = 0.0; cnt = 0;
		for (i = 0; i < m; i++) if (q8 <= (double)s[i] && (double)s[i] <= q99) { acc += (double)s[i]; cnt++; }
		st->high_intensity = cnt ? (float)(acc / (double)cnt) : NAN;
		acc = 0.0; cnt = 0;
		for (i = 0; i < m; i++) if ((double)s[i] <= q15) { acc += (double)s[i]; cnt++; }
		st->noise = cnt ? (float)(acc / (double)cnt) : NAN;
		st->noise_stage = 1;
		if (st->noise > st->high_intensity / 10.0f) {
			acc = 0.0; cnt = 0;
			for (i = 0; i < m; i++) if ((double)s[i] <= q05) { acc += (double)s[i]; cnt++; }
			st->noise = cnt ? (float)(acc / (double)cnt) : NAN;
			st->noise_stage = 2;
			if (st->noise > st->high_intensity / 10.0f) { st->noise = 0.0f; st->noise_stage = 3; }
		}
		threshold = rm_maxf(5.0f * st->noise, st->high_intensity / 5.0f);
		st->threshold = threshold;
		free(s);
	} else {
		threshold = (float)thr_in;
		st->threshold = threshold;
	}

	/* Allocated as each stage is reached, not all six up front: without -romeo-dump the peak is
	   then ~2 live buffers instead of 6 (12 bytes/voxel -> ~5). */
	st->s1 = (uint8_t *)malloc((size_t)n);
	st->sm1 = (float *)malloc((size_t)n * sizeof(float));
	if (!st->s1 || !st->sm1) { rm_mask_stages_free(st); return 1; }
	for (i = 0; i < n; i++) st->s1[i] = (weight[i] > threshold) ? 1 : 0;
	for (i = 0; i < n; i++) st->sm1[i] = (float)st->s1[i];
	{
		int boxes1[1] = { 5 };
		if (rm_boxsmooth3d(st->sm1, nx, ny, nz, 1, boxes1)) { rm_mask_stages_free(st); return 1; }
	}
	st->s2 = (uint8_t *)malloc((size_t)n);
	if (!st->s2) { rm_mask_stages_free(st); return 1; }
	for (i = 0; i < n; i++) st->s2[i] = ((double)st->sm1[i] > 0.4) ? 1 : 0;
	RM_DROP(st->sm1);
	RM_DROP(st->s1);
	st->s3 = (uint8_t *)malloc((size_t)n);
	if (!st->s3) { rm_mask_stages_free(st); return 1; }
	memcpy(st->s3, st->s2, (size_t)n);
	RM_DROP(st->s2);
	if (rm_fill_holes(st->s3, nx, ny, nz)) { rm_mask_stages_free(st); return 1; }
	st->sm2 = (float *)malloc((size_t)n * sizeof(float));
	if (!st->sm2) { rm_mask_stages_free(st); return 1; }
	for (i = 0; i < n; i++) st->sm2[i] = (float)st->s3[i];
	RM_DROP(st->s3);
	{
		int boxes2[2] = { 3, 3 };
		if (rm_boxsmooth3d(st->sm2, nx, ny, nz, 2, boxes2)) { rm_mask_stages_free(st); return 1; }
	}
	st->s4 = (uint8_t *)malloc((size_t)n);
	if (!st->s4) { rm_mask_stages_free(st); return 1; }
	for (i = 0; i < n; i++) st->s4[i] = ((double)st->sm2[i] > 0.6) ? 1 : 0;
	RM_DROP(st->sm2);
	return 0;
}
#undef RM_DROP

/* ============================================================================================
 * 7. Priority queues (ROMEO.jl src/priorityqueue.jl)
 *
 * PQueue is a bucket queue with LIFO semantics inside a bin: enqueue! pushes and dequeue! pops
 * from the END.  q.min is lowered by enqueue! and advanced past empty bins by dequeue!.  A heap
 * or a FIFO would build a different (still valid) spanning tree and therefore assign different
 * 2*pi multiples in ambiguous regions.  Ported literally.
 * ==========================================================================================*/

typedef struct {
	int64_t **bin;
	int64_t *len, *cap;
	int nbins;
	int min;   /* 1-based bin index; nbins+1 == empty */
	int oom;   /* sticky: a failed enqueue would silently drop an edge from the spanning tree,
	              leaving a whole region unwrapped WRONG with a zero exit status. Fail closed. */
} rm_pq;

static int rm_pq_init(rm_pq *q, int nbins) {
	q->nbins = nbins;
	q->min = nbins + 1;
	q->oom = 0;
	q->bin = (int64_t **)calloc((size_t)nbins + 1, sizeof(int64_t *));
	q->len = (int64_t *)calloc((size_t)nbins + 1, sizeof(int64_t));
	q->cap = (int64_t *)calloc((size_t)nbins + 1, sizeof(int64_t));
	if (!q->bin || !q->len || !q->cap) { free(q->bin); free(q->len); free(q->cap); memset(q, 0, sizeof *q); return 1; }
	return 0;
}
static void rm_pq_free(rm_pq *q) {
	int i;
	if (q->bin) for (i = 1; i <= q->nbins; i++) free(q->bin[i]);
	free(q->bin); free(q->len); free(q->cap);
	memset(q, 0, sizeof *q);
}
static int rm_pq_isempty(const rm_pq *q) { return q->min > q->nbins; }
static int rm_pq_enqueue(rm_pq *q, int64_t item, int w) {
	if (w < 1 || w > q->nbins) { q->oom = 1; return 1; }
	if (q->len[w] == q->cap[w]) {
		int64_t nc = q->cap[w] ? q->cap[w] * 2 : 64;
		size_t bytes;
		int64_t *nb;
		/* The queue holds DUPLICATE edge insertions, so its power-of-two growth is not bounded by
		   the top-level n3 guard; on a 32-bit target nc*8 can wrap to 0 and realloc(p,0) may
		   return non-NULL, letting the store below run through a zero-sized allocation. */
		if (nii_mul_size((size_t)nc, sizeof(int64_t), &bytes)) { q->oom = 1; return 1; }
		nb = (int64_t *)realloc(q->bin[w], bytes);
		if (!nb) { q->oom = 1; return 1; }
		q->bin[w] = nb; q->cap[w] = nc;
	}
	q->bin[w][q->len[w]++] = item;
	if (w < q->min) q->min = w;
	return 0;
}
static int64_t rm_pq_dequeue(rm_pq *q) {
	int64_t e = q->bin[q->min][--q->len[q->min]];
	while (q->min <= q->nbins && q->len[q->min] == 0) q->min++;
	return e;
}

/* getseedqueue + findseed!, collapsed to one scan.
 *
 * Upstream builds a 3*NBINS bucket queue over sum(weights; dims=1) -- with every ZERO weight
 * first substituted by 255, so a voxel with non-existent edges sorts as WORST rather than best
 * -- inserts every voxel in ascending linear index, and pops from the END of the lowest
 * non-empty bin. With maxseeds capped at 1 (the only value this build accepts) the queue is
 * dequeued exactly ONCE against an all-zero `visited`, so the result is simply: the smallest
 * substituted weight sum, ties broken toward the HIGHEST linear index. That is one O(n) scan
 * instead of an int64_t[n3] plus bucket metadata (128 MiB at 256^3).
 *
 * If -max-seeds > 1 is ever ported, the bucket queue must come back: later seeds depend on the
 * ordering of the remainder and on which voxels have since been visited. */
static int64_t rm_find_seed(const uint8_t *w, int64_t n) {
	int64_t i, best = 0;
	int bestsum = 3 * 255 + 1;
	for (i = 0; i < n; i++) {
		int a = w[3 * i] ? w[3 * i] : 255;
		int b = w[3 * i + 1] ? w[3 * i + 1] : 255;
		int c = w[3 * i + 2] ? w[3 * i + 2] : 255;
		int sum = a + b + c;
		if (sum <= bestsum) { bestsum = sum; best = i + 1; }   /* <= : highest index wins ties */
	}
	return best;   /* 1-based voxel index; n >= 1 so always found */
}

/* ============================================================================================
 * 8. grow_region_unwrap! (ROMEO.jl src/algorithm.jl + src/seed.jl)
 * ==========================================================================================*/

static int64_t rm_getedgeindex(int64_t leftvoxel, int dim /*1-based*/) { return dim + 3 * (leftvoxel - 1); }
static int rm_getdimfromedge(int64_t edge) { return (int)((edge - 1) % 3) + 1; }
static int64_t rm_getfirstvoxfromedge(int64_t edge) { return (edge - 1) / 3 + 1; }

typedef struct {
	float *wrapped;
	const uint8_t *weights;
	uint8_t *visited;
	int64_t n;
	int64_t stride[3];
	double wrap_addition;
	/* multi-echo seed correction */
	const float *phase2;
	double TE1, TE2;
	int have_p2;
} rm_grow;

/* getnewedge(v, notvisited, stridelist, i) */
static int64_t rm_getnewedge(const rm_grow *g, int64_t v, int i) {
	int iDim = (i + 1) / 2;                    /* div(i+1,2), 1-based dim */
	int64_t nb = g->stride[iDim - 1];
	if (i % 2 == 0) {
		int64_t t = v + nb;
		if (t >= 1 && t <= g->n && g->visited[t - 1] == 0) return rm_getedgeindex(v, iDim);
		return 0;
	}
	{
		int64_t t = v - nb;
		if (t >= 1 && t <= g->n && g->visited[t - 1] == 0) return rm_getedgeindex(t, iDim);
	}
	return 0;
}

static void rm_seedcorrection(rm_grow *g, int64_t vox) {
	if (g->have_p2) {
		double best = INFINITY;
		int offset = 0, off1, off2;
		for (off1 = -2; off1 <= 2; off1++) for (off2 = -1; off2 <= 1; off2++) {
			double diff = fabs(((double)g->wrapped[vox - 1] + RM_2PI_F64 * off1) / g->TE1
				- ((double)g->phase2[vox - 1] + RM_2PI_F64 * off2) / g->TE2);
			diff += (double)(abs(off1) + abs(off2)) / 100.0;
			if (diff < best) { best = diff; offset = off1; }
		}
		g->wrapped[vox - 1] = (float)((double)g->wrapped[vox - 1] + RM_2PI_F64 * offset);
	} else {
		g->wrapped[vox - 1] = rm_rem2pi_f(g->wrapped[vox - 1]);
	}
}

/* unwrapedge!.  SUBTLE, and load-bearing: upstream initialises `d = 0` as an *Int*, and only the
   two threshold branches assign the Float64 `x`/`-x`; the pass-through branch assigns the Float32
   difference `v`.  So `wrapped[oldvox] + d` is a Float32 addition when d is Int 0 or Float32 v,
   and a Float64 addition when d is +/-x.  The two feed DIFFERENT unwrapvoxel methods: the Float32
   one rounds `new - old` to Float32 before dividing by 2*pi, which can pick a different wrap count
   than the Float64 subtraction near an odd multiple of pi.  The value of `old` is identical in
   both; only the width of the subtraction differs.  Do not collapse these two paths. */
static void rm_unwrapedge(rm_grow *g, int64_t oldvox, int64_t newvox) {
	int64_t oo = 2 * oldvox - newvox;
	double x = g->wrap_addition, dd = 0.0;
	float df = 0.0f;
	int d_is_f64 = 0;
	if (oo >= 1 && oo <= g->n && g->visited[oo - 1] != 0) {
		float v = g->wrapped[oldvox - 1] - g->wrapped[oo - 1];
		if ((double)v < -x) { dd = -x; d_is_f64 = 1; }
		else if ((double)v > x) { dd = x; d_is_f64 = 1; }
		else { df = v; d_is_f64 = 0; }
	}
	if (d_is_f64)
		g->wrapped[newvox - 1] = rm_unwrapvoxel_fd(g->wrapped[newvox - 1], (double)g->wrapped[oldvox - 1] + dd);
	else
		g->wrapped[newvox - 1] = rm_unwrapvoxel_ff(g->wrapped[newvox - 1], g->wrapped[oldvox - 1] + df);
}

/* new_seed_thresh = NBINS - div(NBINS - sum(seed_weights)/3, 2).  sum/3 is Float64 and
   div(::Float64, 2) truncates toward zero.  ONE definition: the -romeo-dump manifest reports
   the same value, and two copies of this expression would drift. */
static double rm_seed_thresh(int w1, int w2, int w3) {
	double t = (double)RM_NBINS - (double)(w1 + w2 + w3) / 3.0;
	return (double)RM_NBINS - trunc(t / 2.0);
}

/* Returns the new seed threshold, or 255 when no unvisited voxel remains. */
static double rm_addseed(rm_grow *g, rm_pq *pq, int64_t *seeds, int *nseeds) {
	int64_t seed = rm_find_seed(g->weights, g->n);
	int i;
	if (seed == 0 || g->visited[seed - 1] != 0) return 255.0;
	for (i = 1; i <= 6; i++) {
		int64_t e = rm_getnewedge(g, seed, i);
		if (e != 0 && g->weights[e - 1] > 0) rm_pq_enqueue(pq, e, g->weights[e - 1]);
	}
	rm_seedcorrection(g, seed);
	seeds[*nseeds] = seed;
	(*nseeds)++;
	g->visited[seed - 1] = (uint8_t)(*nseeds);
	return rm_seed_thresh((int)g->weights[rm_getedgeindex(seed, 1) - 1],
		(int)g->weights[rm_getedgeindex(seed, 2) - 1],
		(int)g->weights[rm_getedgeindex(seed, 3) - 1]);
}

/* grow_region_unwrap!.  maxseeds is capped at 255 upstream; only 1 is supported here (the
   experimental multi-seed/region-merging path is not ported). `pq` may already hold seed edges
   (the temporal-uncertain re-entry), in which case no seed is created. */
static int rm_grow_region(rm_grow *g, rm_pq *pq, int seeded_externally, int maxseeds) {
	int64_t seeds[256];
	int nseeds = 0;
	double new_seed_thresh = 256.0;
	int seeded = 0;
	if (rm_pq_isempty(pq)) {
		if (seeded_externally) return 1;
		new_seed_thresh = rm_addseed(g, pq, seeds, &nseeds);
		seeded = 1;
		if (pq->oom) return 1;
	}
	while (!rm_pq_isempty(pq)) {
		int64_t edge, oldvox, newvox, vox, neighbor;
		int dim, i;
		if (seeded && nseeds < maxseeds && (double)pq->min > new_seed_thresh)
			new_seed_thresh = rm_addseed(g, pq, seeds, &nseeds);
		edge = rm_pq_dequeue(pq);
		dim = rm_getdimfromedge(edge);
		vox = rm_getfirstvoxfromedge(edge);
		neighbor = vox + g->stride[dim - 1];
		if (g->visited[neighbor - 1] == 0) { oldvox = vox; newvox = neighbor; }
		else { oldvox = neighbor; newvox = vox; }
		if (g->visited[newvox - 1] == 0) {
			rm_unwrapedge(g, oldvox, newvox);
			g->visited[newvox - 1] = g->visited[oldvox - 1];
			for (i = 1; i <= 6; i++) {
				int64_t e = rm_getnewedge(g, newvox, i);
				if (e != 0 && g->weights[e - 1] > 0) rm_pq_enqueue(pq, e, g->weights[e - 1]);
			}
			if (pq->oom) return 1;
		}
	}
	return pq->oom ? 1 : 0;   /* fail closed on a dropped edge */
}

/* ============================================================================================
 * 9. unwrap! (ROMEO.jl src/unwrapping.jl)
 * ==========================================================================================*/

/* correctglobal: wrapped .-= 2pi*median(round.(filter(isfinite, wrapped[mask]) ./ 2pi))
 *
 * DELIBERATE DIVERGENCE (safer, and documented rather than reproduced): when nothing inside the
 * mask is finite, Julia reaches `median([])` and throws. Here the correction is simply skipped
 * and the phase is returned unchanged. Erroring out on an all-NaN volume adds nothing a caller
 * can act on, and `-g` is a cosmetic global offset, not part of the unwrap. Identical on any
 * input with at least one finite in-mask voxel, i.e. every real image. */
static int rm_correctglobal(float *w, int64_t n, const uint8_t *mask) {
	int64_t i, m = 0;
	double med;
	double *v = (double *)malloc((size_t)n * sizeof(double));
	if (!v) return 1;
	for (i = 0; i < n; i++) {
		if (mask && !mask[i]) continue;
		if (!isfinite(w[i])) continue;
		v[m++] = nearbyint((double)w[i] / RM_2PI_F64);
	}
	if (m == 0) { free(v); return 0; }
	med = rm_median_d(v, m);
	free(v);
	for (i = 0; i < n; i++) w[i] = (float)((double)w[i] - RM_2PI_F64 * med);
	return 0;
}

/* 3D unwrap on an already-computed weight array. */
static int rm_unwrap3d(float *wrapped, const uint8_t *weights, int nx, int ny, int nz,
	const float *phase2, double TE1, double TE2, int have_p2,
	double wrap_addition, int maxseeds, uint8_t *visited_out) {
	int64_t n = (int64_t)nx * ny * nz, i;
	rm_grow g;
	rm_pq pq;
	int rc;
	uint8_t *visited = (uint8_t *)calloc((size_t)n, 1);
	int64_t wsum = 0;
	if (!visited) return 1;
	for (i = 0; i < 3 * n; i++) wsum += weights[i];
	if (wsum == 0) { free(visited); RM_ERR("unwrap weights are all zero\n"); return 1; }
	g.wrapped = wrapped; g.weights = weights; g.visited = visited; g.n = n;
	g.stride[0] = 1; g.stride[1] = nx; g.stride[2] = (int64_t)nx * ny;
	g.wrap_addition = wrap_addition;
	g.phase2 = phase2; g.TE1 = TE1; g.TE2 = TE2; g.have_p2 = have_p2;
	if (rm_pq_init(&pq, RM_NBINS)) { free(visited); return 1; }
	rc = rm_grow_region(&g, &pq, 0, maxseeds);
	rm_pq_free(&pq);
	if (visited_out) memcpy(visited_out, visited, (size_t)n);
	free(visited);
	return rc;
}

/* ============================================================================================
 * 10. NIfTI helpers, options, parser and runner
 * ==========================================================================================*/

static int rm_dump(const char *dir, const char *name, const void *p, size_t elem, int64_t n) {
	char path[2048];
	FILE *f;
	if (!dir) return 0;
	if (snprintf(path, sizeof path, "%s/%s", dir, name) >= (int)sizeof path) {
		/* Truncation would silently alias two dumps onto one filename (c_qmap_1 / c_qmap_2). */
		RM_ERR("-romeo-dump directory path is too long for '%s'\n", name);
		return 1;
	}
	f = fopen(path, "wb");
	if (!f) { RM_ERR("cannot write dump '%s'\n", path); return 1; }
	if (n > 0 && fwrite(p, elem, (size_t)n, f) != (size_t)n) { fclose(f); RM_ERR("short write to '%s'\n", path); return 1; }
	fclose(f);
	return 0;
}

/* ------------------------------------------------------------------------------------------
 * Primitive self-dump (test hook).  The 2*pi range reduction, gamma, rescale and unwrapvoxel are
 * exercised on a FIXED input table that is duplicated verbatim in test/romeo_oracle.jl (medic_bench repo), so the
 * two sides can be byte-compared.  This is the ONLY coverage of the Payne-Hanek branch: real
 * phase data never reaches |x| >= 2^20*pi/2, but -no-phase-rescale lets a caller feed already
 * unwrapped phase, and phaselinearity would then take arbitrarily large arguments.
 * KEEP THE TABLES IN SYNC WITH test/romeo_oracle.jl (medic_bench repo).
 * ----------------------------------------------------------------------------------------*/

static const double RM_PRIM_D[] = {
	0.0, 1.0, -1.0, 3.141592653589793, -3.141592653589793,
	1.5707963267948966, -1.5707963267948966, 4.71238898038469, 6.283185307179586, -6.283185307179586,
	2.5, -2.5, 3.5, -3.5, 5.0, -5.0, 7.0, -7.0,
	100.0, -100.0, 1000.0, 100000.0, 1000000.0, 1600000.0,
	/* >= 2^20*pi/2 == the Payne-Hanek branch */
	1650000.0, -1650000.0, 1.0e7, -1.0e7, 1.0e10, -1.0e10, 1.0e15, 1.0e20, 1.0e30, 1.0e100,
	0.5, -0.5, 1.5, 2.5000000000000004,
	/* Subnormal inputs. NOTE, correcting an earlier over-claim in this file: these are NOT an
	   FTZ/DAZ probe. Every one takes the |x| <= pi pass-through, and ROMEO's real operands
	   (weights in [0,1], phase ~[-pi,pi], magnitudes ~1e3) never reach the denormal range at
	   all -- which is precisely why a gcc -ffast-math link's process-wide MXCSR FTZ/DAZ has nil
	   exposure here. They pin the accepted input DOMAIN, not the FP mode. */
	5e-324, -5e-324, 1e-310, 2.2250738585072011e-308
};
static const float RM_PRIM_F[] = {
	0.0f, 1.0f, -1.0f, 3.1415925f, 3.1415927f, -3.1415925f, -3.1415927f,
	3.2f, -3.2f, 6.2831855f, -6.2831855f, 9.42477f, 12.566371f,
	1.0e-30f, -1.0e-30f, 100.0f, 1.0e5f, 1.0e7f, 1.0e10f, -1.0e10f, 1.0e20f, 4.0f, -4.0f,
	1.4012984643e-45f, -1.4012984643e-45f, 1.1754942107e-38f   /* subnormal / smallest normal */
};
static const double RM_PRIM_W[] = { /* rescale() inputs, incl. exact bin boundaries */
	0.0, 1.0, 0.5, 0.5 / 255.0, 1.5 / 255.0, 2.5 / 255.0,
	1.0 - 0.5 / 255.0, 1.0 - 1.5 / 255.0, 1.0000000000000002, -1.0e-17,
	0.9980392156862745, 0.99607843137254903, 0.25, 0.75
};
static const float RM_PRIM_UV[][2] = { /* (new, old) pairs for unwrapvoxel */
	{ 1.0f, 1.0f }, { 3.0f, -3.0f }, { 0.0f, 3.1415927f }, { 0.0f, -3.1415927f },
	{ 1.0f, 7.2831855f }, { -1.0f, 100.0f }, { 0.0f, 1.0e10f }, { 3.1415927f, -3.1415927f },
	{ 2.0f, 2.0f + 3.1415927f }, { -2.0f, -2.0f - 3.1415927f }
};

static int rm_dump_primitives(const char *dir) {
	size_t i;
	int drc = 0;
	size_t nd = sizeof RM_PRIM_D / sizeof RM_PRIM_D[0];
	size_t nf = sizeof RM_PRIM_F / sizeof RM_PRIM_F[0];
	size_t nw = sizeof RM_PRIM_W / sizeof RM_PRIM_W[0];
	size_t nu = sizeof RM_PRIM_UV / sizeof RM_PRIM_UV[0];
	double *od = (double *)malloc(nd * sizeof(double));
	float *of = (float *)malloc(nf * sizeof(float) * 2);
	uint8_t *ow = (uint8_t *)malloc(nw);
	float *ou = (float *)malloc(nu * sizeof(float) * 2);
	if (!od || !of || !ow || !ou) { free(od); free(of); free(ow); free(ou); return 1; }
	for (i = 0; i < nd; i++) od[i] = rm_rem2pi_d(RM_PRIM_D[i]);
	for (i = 0; i < nf; i++) { of[i] = rm_rem2pi_f(RM_PRIM_F[i]); of[nf + i] = rm_gamma_f(RM_PRIM_F[i]); }
	for (i = 0; i < nw; i++) ow[i] = rm_rescale(RM_PRIM_W[i]);
	for (i = 0; i < nu; i++) {
		ou[i] = rm_unwrapvoxel_ff(RM_PRIM_UV[i][0], RM_PRIM_UV[i][1]);
		ou[nu + i] = rm_unwrapvoxel_fd(RM_PRIM_UV[i][0], (double)RM_PRIM_UV[i][1]);
	}
	drc |= rm_dump(dir, "c_prim_rem2pi64.f64", od, sizeof(double), (int64_t)nd);
	drc |= rm_dump(dir, "c_prim_rem2pi32_gamma.f32", of, sizeof(float), (int64_t)(2 * nf));
	drc |= rm_dump(dir, "c_prim_rescale.u8", ow, 1, (int64_t)nw);
	drc |= rm_dump(dir, "c_prim_unwrapvoxel.f32", ou, sizeof(float), (int64_t)(2 * nu));
	free(od); free(of); free(ow); free(ou);
	return drc;
}

romeo_opts romeo_opts_default(void) {
	romeo_opts o;
	memset(&o, 0, sizeof o);
	o.weights_sel = RM_W_ROMEO;
	o.mask_sel = RM_MASK_ROBUST;
	o.template_echo = 1;
	o.maxseeds = 1;
	o.qmask_thresh = 0.1;
	o.b0_name = "B0";
	o.b0_weighting = RM_B0_PHASE_SNR;
	return o;
}

/* strtod that must consume the whole token */
static int rm_parse_double(const char *s, double *out) {
	char *end = NULL;
	double v;
	if (!s || !*s) return 1;
	v = strtod(s, &end);
	if (!end || *end != '\0') return 1;
	*out = v;
	return 0;
}
static int rm_parse_int(const char *s, long *out) {
	char *end = NULL;
	long v;
	if (!s || !*s) return 1;
	v = strtol(s, &end, 10);
	if (!end || *end != '\0') return 1;
	*out = v;
	return 0;
}

/* -t: "16.8", "16.8,38.56", "[16.8,38.56]", "[16.8 38.56]", "epi", "epi 5.3".  No eval(). */
static int rm_parse_TEs(const char *s, romeo_opts *o) {
	char buf[4096];
	size_t i, j = 0, len;
	const char *tok;
	char *save = NULL;
	if (!s) return 1;
	len = strlen(s);
	if (len + 1 > sizeof buf) return 1;
	for (i = 0; i < len; i++) {
		char c = s[i];
		if (c == '[' || c == ']' || c == '(' || c == ')') continue;
		buf[j++] = (c == ',' || c == ';') ? ' ' : c;
	}
	buf[j] = '\0';
	o->nTE = 0;
	for (tok = strtok_r(buf, " \t", &save); tok; tok = strtok_r(NULL, " \t", &save)) {
		double v;
		if (o->nTE >= ROMEO_MAX_TE) return 1;
		if (rm_parse_double(tok, &v)) return 1;
		o->TEs[o->nTE++] = v;
	}
	return o->nTE > 0 ? 0 : 1;
}

/* -w: romeo | romeo2 | romeo3 | romeo4 | romeo6 | <flag bits, e.g. 1010> */
static int rm_parse_weights(const char *s, romeo_opts *o) {
	size_t i, len;
	if (!s) return 1;
	if (!strcmp(s, "romeo")) { o->weights_sel = RM_W_ROMEO; return 0; }
	if (!strcmp(s, "romeo2")) { o->weights_sel = RM_W_ROMEO2; return 0; }
	if (!strcmp(s, "romeo3")) { o->weights_sel = RM_W_ROMEO3; return 0; }
	if (!strcmp(s, "romeo4")) { o->weights_sel = RM_W_ROMEO4; return 0; }
	if (!strcmp(s, "romeo6")) { o->weights_sel = RM_W_ROMEO6; return 0; }
	if (!strcmp(s, "bestpath")) {
		RM_ERR("-w bestpath (Abdul-Rahman weights) is not implemented in this build\n");
		return 1;
	}
	len = strlen(s);
	if (len < 1 || len > 6) {
		RM_ERR("unknown -w '%s' (romeo|romeo2|romeo3|romeo4|romeo6|<up to 6 bits, e.g. 1010>)\n", s);
		return 1;
	}
	for (i = 0; i < len; i++) if (s[i] != '0' && s[i] != '1') {
		RM_ERR("unknown -w '%s' (romeo|romeo2|romeo3|romeo4|romeo6|<up to 6 bits, e.g. 1010>)\n", s);
		return 1;
	}
	o->weights_sel = RM_W_FLAGS;
	for (i = 0; i < 6; i++) o->flags[i] = 0;
	for (i = 0; i < len; i++) o->flags[i] = (s[i] == '1');
	return 0;
}

int romeo_parse_subopts(int *pac, int argc, char *argv[], romeo_opts *o, const char *cmd) {
	int ac = *pac;
	(void)cmd;  /* reserved: disambiguates the RM_ERR prefix if a second caller is ever added */
	while (ac < argc) {
		const char *a = argv[ac];
		if (!strcmp(a, "-t")) {
			if (ac + 1 >= argc) { RM_ERR("-t requires echo time(s)\n"); return 1; }
			if (!strcmp(argv[ac + 1], "epi")) {
				o->te_epi = 1;
				o->nTE = 1; o->TEs[0] = 1.0;
				ac += 2;
				if (ac < argc) {
					double v;
					if (!rm_parse_double(argv[ac], &v)) { o->TEs[0] = v; ac++; }
				}
				continue;
			}
			if (rm_parse_TEs(argv[ac + 1], o)) { RM_ERR("cannot parse -t '%s'\n", argv[ac + 1]); return 1; }
			ac += 2;
			continue;
		}
		if (!strcmp(a, "-k")) {
			const char *spec;
			if (ac + 1 >= argc) { RM_ERR("-k requires nomask|robustmask|qualitymask [thr]|<file>\n"); return 1; }
			spec = argv[ac + 1];
			ac += 2;
			if (!strcmp(spec, "nomask")) o->mask_sel = RM_MASK_NONE;
			else if (!strcmp(spec, "robustmask")) o->mask_sel = RM_MASK_ROBUST;
			else if (!strcmp(spec, "qualitymask")) {
				o->mask_sel = RM_MASK_QUALITY;
				if (ac < argc) {
					double v;
					if (!rm_parse_double(argv[ac], &v)) { o->qmask_thresh = v; ac++; }
				}
			} else {
				/* set_mask! treats any token that names an existing file as a mask and errors
				   otherwise, with a hint when the token looks like a bare qualitymask threshold. */
				FILE *probe = fopen(spec, "rb");
				if (!probe) {
					double v;
					if (!rm_parse_double(spec, &v))
						RM_ERR("masking option '%s' is undefined (Maybe '-k qualitymask %s' was meant?)\n", spec, spec);
					else
						RM_ERR("masking option '%s' is undefined (nomask | robustmask | qualitymask [thr] | <mask file>)\n", spec);
					return 1;
				}
				fclose(probe);
				o->mask_sel = RM_MASK_FILE;
				o->mask_file = spec;
			}
			continue;
		}
		if (!strcmp(a, "-w")) {
			if (ac + 1 >= argc) { RM_ERR("-w requires a weight specification\n"); return 1; }
			if (rm_parse_weights(argv[ac + 1], o)) return 1;
			ac += 2;
			continue;
		}
		if (!strcmp(a, "-template")) {
			long v;
			/* Range-check BEFORE narrowing to int: a value above INT_MAX wraps NEGATIVE, and a
			   negative echo index reaches `phase + (template-1)*n3` as a wild pointer (SIGSEGV),
			   while e.g. 4294967297 truncates to 1 and silently unwraps the WRONG echo. */
			if (ac + 1 >= argc || rm_parse_int(argv[ac + 1], &v) || v < 1 || v > INT_MAX) {
				RM_ERR("-template requires a positive integer echo index (1..%d)\n", INT_MAX);
				return 1;
			}
			o->template_echo = (int)v;
			ac += 2;
			continue;
		}
		if (!strcmp(a, "-max-seeds")) {
			long v;
			if (ac + 1 >= argc || rm_parse_int(argv[ac + 1], &v) || v < 1) { RM_ERR("-max-seeds requires a positive integer\n"); return 1; }
			if (v != 1) { RM_ERR("-max-seeds > 1 (with -merge-regions/-correct-regions) is not implemented in this build\n"); return 1; }
			o->maxseeds = (int)v;
			ac += 2;
			continue;
		}
		if (!strcmp(a, "-wrap-addition")) {
			double v;
			if (ac + 1 >= argc || rm_parse_double(argv[ac + 1], &v)) { RM_ERR("-wrap-addition requires a number\n"); return 1; }
			if (v != 0.0) { RM_ERR("-wrap-addition != 0 is not implemented in this build\n"); return 1; }
			o->wrap_addition = v;
			ac += 2;
			continue;
		}
		if (!strcmp(a, "-temporal-uncertain-unwrapping")) {
			double v;
			ac++;
			o->temporal_uncertain = 0.5;
			if (ac < argc && !rm_parse_double(argv[ac], &v)) { o->temporal_uncertain = v; ac++; }
			continue;
		}
		if (!strcmp(a, "-g")) { o->correctglobal = 1; ac++; continue; }
		if (!strcmp(a, "-i")) { o->individual = 1; ac++; continue; }
		if (!strcmp(a, "-v")) { o->verbose = 1; ac++; continue; }
		if (!strcmp(a, "-q")) { o->write_quality = 1; ac++; continue; }
		if (!strcmp(a, "-Q")) { o->write_quality_all = 1; ac++; continue; }
		if (!strcmp(a, "-B")) {
			o->compute_b0 = 1;
			ac++;
			/* nargs='?' upstream: take the next token as the output stem only when it is not
			   another option. The main parser already removed the trailing output filename. */
			if (ac < argc && argv[ac][0] != '-') {
				/* The stem becomes a nifti_save POSTFIX (<out>_<stem>, <out>_<stem>_snr), so it
				   must not truncate, must not carry a path or extension, and must not collide
				   with a side output this command already writes. */
				static const char *const reserved[] = { "mask", "quality", "quality_1", "quality_2",
					"quality_3", "quality_4", "quality_5", "quality_6" };
				const char *nm = argv[ac];
				size_t k, len = strlen(nm);
				if (len < 1 || len > 24) {
					RM_ERR("-B name '%s' must be 1-24 characters\n", nm);
					return 1;
				}
				for (k = 0; k < len; k++) {
					if (!((nm[k] >= 'A' && nm[k] <= 'Z') || (nm[k] >= 'a' && nm[k] <= 'z') ||
						  (nm[k] >= '0' && nm[k] <= '9') || nm[k] == '_' || nm[k] == '-')) {
						RM_ERR("-B name '%s' may only contain letters, digits, '_' and '-' (it is a filename POSTFIX on the output, not a path)\n", nm);
						return 1;
					}
				}
				for (k = 0; k < sizeof reserved / sizeof reserved[0]; k++)
					if (!strcmp(nm, reserved[k])) {
						RM_ERR("-B name '%s' would overwrite the <out>_%s side output; choose another\n", nm, nm);
						return 1;
					}
				o->b0_name = nm;
				ac++;
			}
			continue;
		}
		if (!strcmp(a, "-B0-phase-weighting")) {
			const char *m;
			if (ac + 1 >= argc) { RM_ERR("-B0-phase-weighting requires a mode\n"); return 1; }
			m = argv[ac + 1];
			if (!strcmp(m, "phase_snr")) o->b0_weighting = RM_B0_PHASE_SNR;
			else if (!strcmp(m, "phase_var")) o->b0_weighting = RM_B0_PHASE_VAR;
			else if (!strcmp(m, "average")) o->b0_weighting = RM_B0_AVERAGE;
			else if (!strcmp(m, "TEs")) o->b0_weighting = RM_B0_TES;
			else if (!strcmp(m, "mag")) o->b0_weighting = RM_B0_MAG;
			else if (!strcmp(m, "simulated_mag")) o->b0_weighting = RM_B0_SIMULATED_MAG;
			else {
				RM_ERR("the phase weighting option '%s' is not defined (phase_snr|phase_var|average|TEs|mag|simulated_mag)\n", m);
				return 1;
			}
			ac += 2;
			continue;
		}
		if (!strcmp(a, "-no-mask-out")) { o->no_mask_out = 1; ac++; continue; }
		if (!strcmp(a, "-no-phase-rescale") || !strcmp(a, "-no-rescale")) { o->no_phase_rescale = 1; ac++; continue; }
		if (!strcmp(a, "-romeo-dump")) {
			if (ac + 1 >= argc) { RM_ERR("-romeo-dump requires a directory\n"); return 1; }
			o->dumpdir = argv[ac + 1];
			ac += 2;
			continue;
		}
		/* Options that exist upstream but are deliberately not ported yet: reject explicitly
		   rather than letting them fall through to niimath as an unknown operation. */
		if (!strcmp(a, "-u") || !strcmp(a, "-e") || !strcmp(a, "-threshold") ||
			!strcmp(a, "-merge-regions") || !strcmp(a, "-correct-regions") ||
			!strcmp(a, "-fix-ge-phase")) {
			RM_ERR("'%s' is a ROMEO option that is not implemented in this build\n", a);
			return 1;
		}
		break; /* first unrecognized token: back off so niimath sees the next chain operation */
	}
	*pac = ac;
	return 0;
}

/* ---- auxiliary image loading -------------------------------------------------------------- */

/* Read a NIfTI into float32.
     RM_RD_SCALED   apply slope/intercept in FLOAT arithmetic, matching NIfTI.jl's getindex
     RM_RD_RAW      the STORED values (Julia's `.raw`, used by readphase's second branch)
     RM_RD_NONZERO  1.0 where the STORED value is nonzero, else 0.0 -- the mask rule
                    (`niread(f).raw .!= 0`) evaluated in the ORIGINAL width. Narrowing a float64
                    mask to float32 first would turn e.g. 1e-300 into a false zero. */
enum { RM_RD_SCALED = 0, RM_RD_RAW, RM_RD_NONZERO };
static int rm_read_f32(const char *fn, int raw, float **out, int *nx, int *ny, int *nz, int *nvol) {
	nifti_image *n = nifti_image_read(fn, 1);
	int64_t i, nv;
	float *d = NULL;
	float scl, inter;
	if (!n) { RM_ERR("failed to read '%s'\n", fn); return 1; }
	if (n->nx < 1 || n->ny < 1 || n->nz < 1 || n->nvox < 1) { nifti_image_free(n); RM_ERR("'%s' has invalid dimensions\n", fn); return 1; }
	if (n->nu > 1 || n->nv > 1 || n->nw > 1) {   /* 5D multi-channel is out of scope, not "volume 1" */
		nifti_image_free(n);
		RM_ERR("'%s' has more than 4 dimensions (5D multi-channel input is not supported)\n", fn);
		return 1;
	}
	nv = (int64_t)n->nvox;
	if (nv > INT_MAX) { nifti_image_free(n); RM_ERR("'%s' exceeds INT_MAX voxels\n", fn); return 1; }
	{	/* checked multiply: on a 32-bit/wasm target nv*4 can wrap size_t and under-allocate
		   while the int64 conversion loop below still writes all nv elements */
		size_t bytes;
		if (nii_mul_size((size_t)nv, sizeof(float), &bytes)) {
			nifti_image_free(n);
			RM_ERR("'%s' is too large for this build's address space\n", fn);
			return 1;
		}
		d = (float *)malloc(bytes);
	}
	if (!d) { nifti_image_free(n); return 1; }
	scl = (n->scl_slope == 0.0f) ? 1.0f : n->scl_slope;
	inter = n->scl_inter;
	if (raw != RM_RD_SCALED) { scl = 1.0f; inter = 0.0f; }
#define RM_CVT(T) do { const T *p = (const T *)n->data; \
	if (raw == RM_RD_NONZERO) { for (i = 0; i < nv; i++) d[i] = (p[i] != 0) ? 1.0f : 0.0f; } \
	else { for (i = 0; i < nv; i++) d[i] = (float)p[i] * scl + inter; } } while (0)
	switch (n->datatype) {
	case DT_UINT8: RM_CVT(uint8_t); break;
	case DT_INT8: RM_CVT(int8_t); break;
	case DT_INT16: RM_CVT(int16_t); break;
	case DT_UINT16: RM_CVT(uint16_t); break;
	case DT_INT32: RM_CVT(int32_t); break;
	case DT_UINT32: RM_CVT(uint32_t); break;
	case DT_INT64: RM_CVT(int64_t); break;
	case DT_UINT64: RM_CVT(uint64_t); break;
	case DT_FLOAT32: RM_CVT(float); break;
	case DT_FLOAT64: RM_CVT(double); break;
	default: {
		int dt = n->datatype;   /* read BEFORE the free: nifti_image_free(n) invalidates n */
		free(d);
		nifti_image_free(n);
		RM_ERR("'%s' has an unsupported datatype (%d)\n", fn, dt);
		return 1;
	}
	}
#undef RM_CVT
	*out = d;
	*nx = n->nx; *ny = n->ny; *nz = (n->nz < 1 ? 1 : n->nz);
	*nvol = (int)(nv / ((int64_t)(*nx) * (*ny) * (*nz)));
	if (*nvol < 1) *nvol = 1;
	nifti_image_free(n);
	return 0;
}

/* ---- side outputs -------------------------------------------------------------------------- */

static int rm_save_side(nifti_image *nim, const char *postfix, const float *vals, int64_t n3, gzModes gzMode) {
	/* Save a 3D float32 companion on the working image's grid without disturbing nim->data. */
	void *savedata = nim->data;
	int saved_nt = nim->nt, saved_ndim = nim->ndim;
	int64_t saved_nvox = nim->nvox;
	int rc;
	float *buf = (float *)nii_malloc((size_t)n3, sizeof(float));
	memcpy(buf, vals, (size_t)n3 * sizeof(float));
	nim->data = buf;
	nim->nt = 1; nim->dim[4] = 1;
	nim->ndim = 3; nim->dim[0] = 3;
	nim->nvox = n3;
	rc = nifti_save(nim, postfix, gzMode);
	free(buf);
	nim->data = savedata;
	nim->nt = saved_nt; nim->dim[4] = saved_nt;
	nim->ndim = saved_ndim; nim->dim[0] = saved_ndim;
	nim->nvox = saved_nvox;
	return rc;
}

/* Build a weight context, reducing a 4D input to its template echo exactly as ROMEO's 4D
   calculateweights overload does (phase2 = echo p2ref, TEs = [TE_template, TE_p2ref], magnitude
   = volume `template`).  `magmasked` (n3 floats, caller-owned) receives mag .* mask when a mask
   is present, mirroring parsekwargs.  maxmag is the Float64 0.95 quantile of the finite
   (masked) magnitude.  Returns 0 on success. */
static int rm_build_ctx(rm_wctx *c, const float *phase, const float *mag, int magvol,
	const uint8_t *mask, float *magmasked, const double *TEs, int neco, int template_echo,
	int p2ref, int nx, int ny, int nz, const int *flags) {
	int64_t n3 = (int64_t)nx * ny * nz, i;
	int mte = (magvol > 1) ? template_echo : 1;   /* `size(args[:mag],4) > 1` guard upstream */
	memset(c, 0, sizeof *c);
	c->nx = nx; c->ny = ny; c->nz = nz; c->n = n3;
	c->P = phase + (int64_t)(template_echo - 1) * n3;
	c->P2 = (neco > 1) ? (phase + (int64_t)(p2ref - 1) * n3) : NULL;
	c->TE1 = (neco > 1) ? TEs[template_echo - 1] : 1.0;
	c->TE2 = (neco > 1) ? TEs[p2ref - 1] : 1.0;
	c->mask = mask;
	if (mag) {
		const float *magt = mag + (int64_t)(mte - 1) * n3;
		if (mask && magmasked) {
			for (i = 0; i < n3; i++) magmasked[i] = magt[i] * (float)mask[i];
			c->M = magmasked;
		} else c->M = magt;
		{
			float *tmp = (float *)malloc((size_t)n3 * sizeof(float));
			int64_t m = 0;
			if (!tmp) return 1;
			for (i = 0; i < n3; i++) if (isfinite(c->M[i])) tmp[m++] = c->M[i];
			c->maxmag = (m > 0) ? rm_quantile7(tmp, m, 0.95) : 0.0;
			free(tmp);
		}
	}
	memcpy(c->flags, flags, 6 * sizeof(int));
	rm_updateflags(c->flags, c->P2 != NULL, 1, c->M != NULL);
	return 0;
}

/* ---- shared unwrapping core ------------------------------------------------------------------
 *
 * Everything between "phase is in radians in memory" and "the phase is unwrapped": mask
 * selection, weight calculation and the three unwrapping dispatches.  Both entry points run this
 * ONE body, so the file-oriented `-romeo` and the in-memory frame API used by `--medic` cannot
 * drift.  `romeo_run` additionally owns loading, readphase rescaling, the parity dump and the side
 * outputs; `romeo_unwrap_frame` owns none of those.
 *
 * The struct carries the outputs the caller still needs afterwards (mask, weights, flags) and, for
 * `romeo_run` only, the parity-dump plumbing.  On failure rm_core_run returns non-zero having left
 * the struct in a state rm_core_free can release. */
typedef struct {
	/* caller-supplied */
	float *phase;         /* in/out, n3 * neco, radians */
	const uint8_t *mask_in; /* optional: use this mask verbatim instead of selecting one */
	const float *mag;     /* borrowed, n3 * magvol, or NULL */
	int magvol;
	int nx, ny, nz, neco;
	const double *TEs;
	const romeo_opts *o;
	/* derived / owned */
	int64_t n3;
	uint8_t *mask;
	uint8_t *weights;
	uint8_t *visited;     /* dump only */
	float *magmasked;
	rm_mask_stages stages;
	int flags[6];
	int template_echo, p2ref;
	double TE1, TE2;
	int have_mag;
	/* romeo_run-only parity dump */
	const char *dump;
	FILE *manifest;
	int drc;
} rm_core;

static void rm_core_free(rm_core *c) {
	rm_mask_stages_free(&c->stages);
	free(c->mask); free(c->weights); free(c->magmasked); free(c->visited);
	c->mask = NULL; c->weights = NULL; c->magmasked = NULL; c->visited = NULL;
}

static int rm_core_run(rm_core *c) {
	const romeo_opts *o = c->o;
	const int nx = c->nx, ny = c->ny, nz = c->nz, neco = c->neco;
	const int64_t n3 = c->n3;
	const double *TEs = c->TEs;
	const char *dump = c->dump;
	float *phase = c->phase;
	int64_t i;

	c->template_echo = o->template_echo;
	c->p2ref = 0;
	c->TE1 = c->TE2 = 1.0;

	/* ---- weight selection (needed before the mask: -k qualitymask uses it) ------------------ */
	rm_flags_from_sel(o->weights_sel, c->have_mag, c->flags);
	if (o->weights_sel == RM_W_FLAGS) for (i = 0; i < 6; i++) c->flags[i] = o->flags[i];
	if (neco > 1) {
		c->p2ref = (c->template_echo == 1) ? 2 : c->template_echo - 1;
		c->TE1 = TEs[c->template_echo - 1];
		c->TE2 = TEs[c->p2ref - 1];
	}
	if (c->have_mag) {
		c->magmasked = (float *)malloc((size_t)n3 * sizeof(float));
		if (!c->magmasked) return 1;
	}

	/* ---- mask ------------------------------------------------------------------------------ */
	if (c->mask_in) {
		/* Caller supplied the mask (used by --medic, which shares ONE mask across the MCPC-3D-S
		   phase-difference unwrap and the multi-echo unwrap, as the reference does). */
		c->mask = (uint8_t *)malloc((size_t)n3);
		if (!c->mask) return 1;
		memcpy(c->mask, c->mask_in, (size_t)n3);
	} else if (o->mask_sel == RM_MASK_ROBUST && !c->have_mag) {
		/* load_data_and_resolve_args!: robustmask without a magnitude degrades to nomask */
		fprintf(stderr, " + -romeo: robustmask was chosen but no magnitude is available. No mask is used!\n");
	} else if (o->mask_sel == RM_MASK_ROBUST) {
		int te = c->template_echo < c->magvol ? c->template_echo : c->magvol;
		if (rm_robustmask(c->mag + (int64_t)(te - 1) * n3, nx, ny, nz, 0, 0.0, dump != NULL, &c->stages)) return 1;
		c->mask = c->stages.s4; c->stages.s4 = NULL;
	} else if (o->mask_sel == RM_MASK_QUALITY) {
		/* set_mask!: qmap = voxelquality(phase; get_keyargs(...)) — computed on the still-WRAPPED
		   phase and WITHOUT a mask (data["mask"] does not exist yet), then robustmask(qmap; threshold).
		   voxelquality's own 4D overload defaults p2ref to 2 regardless of `template`. */
		rm_wctx qc;
		float *qmap = (float *)malloc((size_t)n3 * sizeof(float));
		if (!qmap) return 1;
		if (rm_build_ctx(&qc, phase, c->have_mag ? c->mag : NULL, c->magvol, NULL, c->magmasked, TEs, neco,
				c->template_echo, 2, nx, ny, nz, c->flags)) { free(qmap); return 1; }
		if (rm_voxelquality(&qc, qmap)) { free(qmap); return 1; }
		if (dump) c->drc |= rm_dump(dump, "c_qmap_wrapped.f32", qmap, sizeof(float), n3);
		if (rm_robustmask(qmap, nx, ny, nz, 1, o->qmask_thresh, dump != NULL, &c->stages)) { free(qmap); return 1; }
		free(qmap);
		c->mask = c->stages.s4; c->stages.s4 = NULL;
	} else if (o->mask_sel == RM_MASK_FILE) {
		float *mv = NULL;
		int mnx, mny, mnz, mnv;
		if (rm_read_f32(o->mask_file, RM_RD_NONZERO, &mv, &mnx, &mny, &mnz, &mnv)) return 1;
		if (mnx != nx || mny != ny || mnz != nz || mnv != 1) {
			free(mv); RM_ERR("mask dimensions do not match the phase\n"); return 1;
		}
		c->mask = (uint8_t *)malloc((size_t)n3);
		if (!c->mask) { free(mv); return 1; }
		for (i = 0; i < n3; i++) c->mask[i] = (mv[i] != 0.0f) ? 1 : 0; /* RM_RD_NONZERO already applied the raw test */
		free(mv);
	}

	/* ---- weights ---------------------------------------------------------------------------- */
	{
		rm_wctx wc;
		if (rm_build_ctx(&wc, phase, c->have_mag ? c->mag : NULL, c->magvol, c->mask, c->magmasked, TEs, neco,
				c->template_echo, neco > 1 ? c->p2ref : 1, nx, ny, nz, c->flags)) return 1;
		c->weights = (uint8_t *)malloc((size_t)3 * (size_t)n3);
		if (!c->weights) return 1;
		rm_calculateweights(&wc, RM_WOUT_U8, c->weights);

		if (dump) {
			char path[2048];
			snprintf(path, sizeof path, "%s/c_manifest.txt", dump);
			c->manifest = fopen(path, "w");
			c->drc |= rm_dump_primitives(dump);
			if (rm_dump(dump, "c_weights.u8", c->weights, 1, 3 * n3)) return 1;
			{
				double *wd = (double *)malloc((size_t)3 * (size_t)n3 * sizeof(double));
				if (wd) {
					rm_calculateweights(&wc, RM_WOUT_F64, wd);
					c->drc |= rm_dump(dump, "c_weights_prerescale.f64", wd, sizeof(double), 3 * n3);
					free(wd);
				}
			}
			if (c->manifest) {
				/* Seed scalars: the plan's M5 gate names them explicitly, so make them
				   directly comparable with the oracle manifest rather than implied by the
				   (bit-exact) weights they are derived from. */
				{
					int64_t sd = rm_find_seed(c->weights, n3);
					fprintf(c->manifest, "seed_index %lld\n", (long long)sd);
					if (sd > 0) {
						int w1 = c->weights[rm_getedgeindex(sd, 1) - 1];
						int w2 = c->weights[rm_getedgeindex(sd, 2) - 1];
						int w3 = c->weights[rm_getedgeindex(sd, 3) - 1];
						fprintf(c->manifest, "seed_w1 %d\nseed_w2 %d\nseed_w3 %d\n", w1, w2, w3);
						fprintf(c->manifest, "new_seed_thresh %.17g\n", rm_seed_thresh(w1, w2, w3));
					}
				}
				fprintf(c->manifest, "flags_active %d%d%d%d%d%d\n",
					wc.flags[0], wc.flags[1], wc.flags[2], wc.flags[3], wc.flags[4], wc.flags[5]);
				if (c->have_mag) fprintf(c->manifest, "maxmag %.17g\n", wc.maxmag);
				if (c->stages.s1) {
					fprintf(c->manifest, "rm_sample_len %lld\n", (long long)c->stages.sample_len);
					fprintf(c->manifest, "rm_q05 %.17g\nrm_q15 %.17g\nrm_q8 %.17g\nrm_q99 %.17g\n",
						c->stages.q05, c->stages.q15, c->stages.q8, c->stages.q99);
					fprintf(c->manifest, "rm_high_intensity %.9g\nrm_noise %.9g\nrm_noise_stage %d\nrm_threshold %.9g\n",
						(double)c->stages.high_intensity, (double)c->stages.noise, c->stages.noise_stage, (double)c->stages.threshold);
				}
			}
			if (c->stages.s1) {
				c->drc |= rm_dump(dump, "c_mask_s1_thresh.u8", c->stages.s1, 1, n3);
				c->drc |= rm_dump(dump, "c_mask_sm1.f32", c->stages.sm1, sizeof(float), n3);
				c->drc |= rm_dump(dump, "c_mask_s2_smooth1.u8", c->stages.s2, 1, n3);
				c->drc |= rm_dump(dump, "c_mask_s3_fill.u8", c->stages.s3, 1, n3);
				c->drc |= rm_dump(dump, "c_mask_sm2.f32", c->stages.sm2, sizeof(float), n3);
			}
			if (c->mask) c->drc |= rm_dump(dump, "c_mask_s4_final.u8", c->mask, 1, n3);
			if (c->drc) { RM_ERR("one or more -romeo-dump writes failed\n"); return 1; }
		}
	}

	if (o->verbose)
		fprintf(stderr, " + -romeo: weights %d%d%d%d%d%d, mask=%s, template echo %d\n",
			c->flags[0], c->flags[1], c->flags[2], c->flags[3], c->flags[4], c->flags[5],
			c->mask ? (o->mask_sel == RM_MASK_FILE ? "file" : (o->mask_sel == RM_MASK_QUALITY ? "qualitymask" : "robustmask")) : "none",
			c->template_echo);

	/* ---- unwrap ----------------------------------------------------------------------------- */
	if (dump) {   /* the outer copy feeds the parity dump only; rm_unwrap3d owns its working set */
		c->visited = (uint8_t *)calloc((size_t)n3, 1);
		if (!c->visited) return 1;
	}
	if (neco == 1) {
		if (rm_unwrap3d(phase, c->weights, nx, ny, nz, NULL, c->TE1, c->TE2, 0,
				o->wrap_addition, o->maxseeds, c->visited)) return 1;
		if (o->correctglobal && rm_correctglobal(phase, n3, c->mask)) return 1;
	} else if (o->individual) {
		/* unwrap_individual!: each echo is unwrapped spatially with its own weights, using the
		   PREVIOUS echo (echo 2 for echo 1) as phase2.  Echoes are processed in ascending order,
		   so when echo i>1 is unwrapped its reference echo i-1 is ALREADY unwrapped — that is
		   upstream behaviour (Threads.@threads with a shared Dict; the oracle pins 1 thread). */
		int ie;
		for (ie = 1; ie <= neco; ie++) {
			int e2 = (ie == 1) ? 2 : ie - 1;
			rm_wctx ic;
			float *p2copy = (float *)malloc((size_t)n3 * sizeof(float));
			uint8_t *w2 = NULL;
			if (!p2copy) return 1;
			memcpy(p2copy, phase + (int64_t)(e2 - 1) * n3, (size_t)n3 * sizeof(float));
			if (rm_build_ctx(&ic, phase, c->have_mag ? c->mag : NULL, c->magvol, c->mask, c->magmasked, TEs, neco,
					ie, e2, nx, ny, nz, c->flags)) { free(p2copy); return 1; }
			w2 = (uint8_t *)malloc((size_t)3 * (size_t)n3);
			if (!w2) { free(p2copy); return 1; }
			rm_calculateweights(&ic, RM_WOUT_U8, w2);
			if (rm_unwrap3d(phase + (int64_t)(ie - 1) * n3, w2, nx, ny, nz, p2copy,
					ic.TE1, ic.TE2, 1, o->wrap_addition, o->maxseeds, c->visited)) {
				free(p2copy); free(w2); return 1;
			}
			if (o->correctglobal && rm_correctglobal(phase + (int64_t)(ie - 1) * n3, n3, c->mask)) {
				free(p2copy); free(w2); return 1;
			}
			free(p2copy); free(w2);
		}
		if (o->correctglobal) {
			/* correct_multi_echo_wraps!
			 *
			 * DELIBERATE DIVERGENCE (safer): upstream filters the reference and current echoes
			 * INDEPENDENTLY before subtracting them, so a NaN present in only one echo either
			 * throws on a length mismatch or silently pairs mismatched voxels. Here a voxel
			 * contributes only when BOTH echoes are finite at that voxel, which is what the
			 * expression means. Identical whenever the two echoes share a finite mask, i.e.
			 * every real image. */
			int ie2;
			double *v = (double *)malloc((size_t)n3 * sizeof(double));
			if (!v) return 1;
			for (ie2 = 2; ie2 <= neco; ie2++) {
				int iref = ie2 - 1;
				double fac = TEs[ie2 - 1] / TEs[iref - 1], nwraps;
				const float *pr = phase + (int64_t)(iref - 1) * n3;
				float *pe = phase + (int64_t)(ie2 - 1) * n3;
				int64_t m = 0;
				for (i = 0; i < n3; i++) {
					if (c->mask && !c->mask[i]) continue;
					if (!isfinite(pr[i]) || !isfinite(pe[i])) continue;
					v[m++] = nearbyint(((double)pr[i] * fac - (double)pe[i]) / RM_2PI_F64);
				}
				if (m == 0) continue;
				nwraps = rm_median_d(v, m);
				for (i = 0; i < n3; i++) pe[i] = (float)((double)pe[i] + RM_2PI_F64 * nwraps);
			}
			free(v);
		}
	} else {
		/* 4D: spatially unwrap the template echo, then propagate temporally. */
		float *tpl = phase + (int64_t)(c->template_echo - 1) * n3;
		const float *p2 = phase + (int64_t)(c->p2ref - 1) * n3;
		float *p2copy = (float *)malloc((size_t)n3 * sizeof(float));
		int order_i;
		if (!p2copy) return 1;
		memcpy(p2copy, p2, (size_t)n3 * sizeof(float)); /* args[:phase2] is a COPY taken up front */
		if (rm_unwrap3d(tpl, c->weights, nx, ny, nz, p2copy, c->TE1, c->TE2, 1,
				o->wrap_addition, o->maxseeds, c->visited)) { free(p2copy); return 1; }
		free(p2copy);
		if (o->correctglobal && rm_correctglobal(tpl, n3, c->mask)) return 1;
		for (order_i = 0; order_i < neco - 1; order_i++) {
			/* iteration order: (template-1):-1:1, then (template+1):neco */
			int ieco = (order_i < c->template_echo - 1) ? (c->template_echo - 1 - order_i) : (order_i + 2);
			int iref = (ieco < c->template_echo) ? ieco + 1 : ieco - 1;
			double fac = TEs[ieco - 1] / TEs[iref - 1];
			float *w = phase + (int64_t)(ieco - 1) * n3;
			const float *r = phase + (int64_t)(iref - 1) * n3;
			double *refvalue = (double *)malloc((size_t)n3 * sizeof(double));
			if (!refvalue) return 1;
			for (i = 0; i < n3; i++) refvalue[i] = (double)r[i] * fac;
			for (i = 0; i < n3; i++) w[i] = rm_unwrapvoxel_fd(w[i], refvalue[i]);
			if (o->temporal_uncertain > 0.0) {
				/* temporal_uncertain_unwrapping!: spatially re-unwrap low-quality voxels */
				rm_wctx qc;
				float *qual = (float *)malloc((size_t)n3 * sizeof(float));
				float *halfw = (float *)malloc((size_t)n3 * sizeof(float));
				double *halfr = (double *)malloc((size_t)n3 * sizeof(double));
				uint8_t *vis = (uint8_t *)malloc((size_t)n3);
				int any = 0, all = 1;
				if (!qual || !halfw || !halfr || !vis) { free(qual); free(halfw); free(halfr); free(vis); free(refvalue); return 1; }
				for (i = 0; i < n3; i++) { halfw[i] = w[i] / 2.0f; halfr[i] = refvalue[i] / 2.0; }
				memset(&qc, 0, sizeof qc);
				qc.P = halfw; qc.P2d = halfr; qc.TE1 = 1.0; qc.TE2 = 1.0;
				qc.nx = nx; qc.ny = ny; qc.nz = nz; qc.n = n3;
				qc.flags[0] = 1; qc.flags[1] = 1; qc.flags[2] = 1; /* :romeo, no mag -> 4..6 off */
				if (rm_voxelquality(&qc, qual)) { free(qual); free(halfw); free(halfr); free(vis); free(refvalue); return 1; }
				for (i = 0; i < n3; i++) vis[i] = ((double)qual[i] > o->temporal_uncertain) ? 1 : 0;
				for (i = 0; i < n3; i++) {
					int inmask;
					if (c->mask) inmask = c->mask[i] != 0;
					else {
						int64_t s = (int64_t)c->weights[3 * i] + c->weights[3 * i + 1] + c->weights[3 * i + 2];
						inmask = (s < 100);
					}
					if (!inmask) vis[i] = 1;
					if (vis[i]) any = 1; else all = 0;
				}
				if (any && !all) {
					rm_grow g;
					rm_pq pq;
					int64_t stride[3];
					int dim;
					stride[0] = 1; stride[1] = nx; stride[2] = (int64_t)nx * ny;
					g.wrapped = w; g.weights = c->weights; g.visited = vis; g.n = n3;
					g.stride[0] = stride[0]; g.stride[1] = stride[1]; g.stride[2] = stride[2];
					g.wrap_addition = o->wrap_addition;
					g.phase2 = NULL; g.TE1 = c->TE1; g.TE2 = c->TE2; g.have_p2 = 0;
					if (rm_pq_init(&pq, RM_NBINS)) { free(qual); free(halfw); free(halfr); free(vis); free(refvalue); return 1; }
					for (dim = 1; dim <= 3; dim++) {
						int64_t I;
						for (I = 1; I <= n3; I++) {
							int64_t J = I + stride[dim - 1];
							if (J > n3) continue;
							if ((int)vis[I - 1] + (int)vis[J - 1] == 1) {
								int64_t ed = rm_getedgeindex(I, dim);
								if (c->weights[ed - 1] != 0) rm_pq_enqueue(&pq, ed, c->weights[ed - 1]);
							}
						}
					}
					if (rm_grow_region(&g, &pq, 1, o->maxseeds)) {
						rm_pq_free(&pq);
						free(qual); free(halfw); free(halfr); free(vis); free(refvalue);
						RM_ERR("out of memory during temporal-uncertain re-unwrapping\n");
						return 1;
					}
					rm_pq_free(&pq);
				}
				free(qual); free(halfw); free(halfr); free(vis);
			}
			free(refvalue);
		}
	}
	return 0;
}

/* ---- in-memory frame API (used by --medic) ---------------------------------------------------
 *
 * No file I/O, no side outputs, no nifti_image.  `phase` is caller-owned, echo-major
 * (n3 floats per echo) and ALREADY in radians -- readphase rescaling and any phase-offset
 * correction are the caller's business.  It is unwrapped IN PLACE.  `mag` is caller-owned and may
 * be NULL.  `mask_out`, when non-NULL, receives the n3-byte mask (all zero if the options select
 * no mask).  Options that only make sense for the CLI (dump, side outputs, rescaling) are ignored.
 */
int romeo_unwrap_frame(float *phase, const float *mag, int magvol,
	int nx, int ny, int nz, int neco, const double *TEs,
	const romeo_opts *o, const uint8_t *mask_in, uint8_t *mask_out) {
	rm_core c;
	int rc;
	if (!phase || !o || nx < 1 || ny < 1 || nz < 1 || neco < 1 || !TEs) return 1;
	memset(&c, 0, sizeof c);
	c.phase = phase;
	c.mag = mag; c.magvol = mag ? magvol : 0;
	c.nx = nx; c.ny = ny; c.nz = nz; c.neco = neco;
	c.n3 = (int64_t)nx * ny * nz;
	c.TEs = TEs; c.o = o;
	c.mask_in = mask_in;
	c.have_mag = (mag != NULL);
	if (o->template_echo < 1 || o->template_echo > neco) return 1;
	if (c.have_mag && magvol < neco) return 1;
	rc = rm_core_run(&c);
	if (mask_out) {
		if (!rc && c.mask) memcpy(mask_out, c.mask, (size_t)c.n3);
		else memset(mask_out, 0, (size_t)c.n3);
	}
	rm_core_free(&c);
	return rc;
}

/* Compute ROMEO's robustmask from a magnitude volume, for callers that need the mask on its own
   (--medic shares ONE mask across its two unwrapping calls).  `mask` is caller-owned, nx*ny*nz
   bytes.  Returns 0 on success. */
int romeo_robustmask(const float *mag, int nx, int ny, int nz, uint8_t *mask) {
	rm_mask_stages st;
	int64_t n3 = (int64_t)nx * ny * nz;
	memset(&st, 0, sizeof st);
	if (!mag || !mask || nx < 1 || ny < 1 || nz < 1) return 1;
	if (rm_robustmask(mag, nx, ny, nz, 0, 0.0, 0, &st)) { rm_mask_stages_free(&st); return 1; }
	if (!st.s4) { rm_mask_stages_free(&st); return 1; }
	memcpy(mask, st.s4, (size_t)n3);
	rm_mask_stages_free(&st);
	return 0;
}

/* Compute ROMEO's voxel-quality map from WRAPPED phase alone -- all-ones magnitude weights, so
   the result reflects phase coherence only.  This is `voxelquality(phase; TEs, ...)` with no
   `mag` keyword, i.e. the same call ROMEO makes internally for -romeo's quality mask except that
   the magnitude is deliberately not supplied.  --medic's tiered mask unions it with an Otsu
   magnitude mask precisely because the two fail differently: the magnitude mask is occasionally
   too aggressive, the quality mask permissive but noisy.

   `phase` is neco wrapped volumes of nx*ny*nz; `qmap` is caller-owned, nx*ny*nz floats.
   Returns 0 on success. */
int romeo_voxelquality(const float *phase, int neco, int nx, int ny, int nz,
	const double *TEs, const romeo_opts *o, float *qmap) {
	rm_wctx qc;
	int flags[6];
	if (!phase || !qmap || !TEs || !o || neco < 1 || nx < 1 || ny < 1 || nz < 1) return 1;
	/* have_mag = 0: romeo/romeo4 resolve to the magnitude-free flag set, which is what
	   "quality from phase alone" means. */
	rm_flags_from_sel(o->weights_sel, 0, flags);
	if (o->weights_sel == RM_W_FLAGS) memcpy(flags, o->flags, sizeof flags);
	/* p2ref = 2 regardless of the template echo: voxelquality's own 4D overload does that. */
	if (rm_build_ctx(&qc, phase, NULL, 0, NULL, NULL, TEs, neco,
			o->template_echo > 0 ? o->template_echo : 1, 2, nx, ny, nz, flags)) return 1;
	return rm_voxelquality(&qc, qmap);
}

/* ---- the runner ---------------------------------------------------------------------------- */

int romeo_run(nifti_image *nim, const char *magfile, const char *phasefile,
	const in_hdr *ihdr, int is_first_op, const romeo_opts *o, gzModes gzMode) {
	int nx, ny, nz, neco, ret = 1;
	int64_t n3, i;
	float *phase = NULL;     /* owned working copy, [n3 * neco] */
	float *mag = NULL;       /* owned, [n3 * magvol] or NULL */
	int magvol = 0;
	uint8_t *mask = NULL;    /* borrowed from core (core owns it) */
	double *TEs = NULL;
	int have_mag = 0;
	int template_echo;
	const char *dump = o->dumpdir;
	rm_core core;
	int drc = 0;   /* accumulated -romeo-dump status: a PARTIAL parity dump must not read as complete */

	memset(&core, 0, sizeof core);

	if (nim->datatype != DT_FLOAT32) { RM_ERR("internal error: expected float32 working image\n"); return 1; }
	if (nim->nu > 1 || nim->nv > 1 || nim->nw > 1) { RM_ERR("input must be 3D or 4D (echoes on dim 4)\n"); return 1; }
	nx = nim->nx; ny = nim->ny; nz = (nim->nz < 1 ? 1 : nim->nz);
	n3 = (int64_t)nx * ny * nz;
	if (n3 < 1 || (int64_t)nim->nvox % n3 != 0) { RM_ERR("invalid image geometry\n"); return 1; }
	neco = (int)((int64_t)nim->nvox / n3);
	if (neco < 1) { RM_ERR("invalid image geometry\n"); return 1; }
	if ((int64_t)nim->nvox > INT_MAX) { RM_ERR("images with more than INT_MAX voxels are not supported\n"); return 1; }
	{	/* The weight arrays are 3 per voxel and up to 8 bytes wide. On a 32-bit/wasm target
		   (FORCE_INT32_MAX, where romeo.c IS linked) those products could wrap size_t and
		   under-allocate, so use the project's checked multiply rather than trusting them. */
		size_t chk;
		if (nii_mul_size((size_t)n3, 3 * sizeof(double), &chk) ||
			nii_mul_size((size_t)nim->nvox, sizeof(float), &chk)) {
			RM_ERR("image is too large for this build's address space\n");
			return 1;
		}
	}

	if (o->compute_b0 && o->nTE == 0) {
		RM_ERR("echo times are required for B0 calculation (-B needs -t)\n");
		return 1;
	}
	template_echo = o->template_echo;
	if (template_echo < 1 || template_echo > neco) {   /* lower bound too: never index behind the buffer */
		RM_ERR("-template %d is out of range (the image has %d echo(es))\n", template_echo, neco);
		return 1;
	}

	/* ---- phase: readphase rescale branch ---------------------------------------------- */
	phase = (float *)malloc((size_t)nim->nvox * sizeof(float));
	if (!phase) goto done;
	memcpy(phase, nim->data, (size_t)nim->nvox * sizeof(float));
	if (!o->no_phase_rescale) {
		float mn, mx;
		if (rm_approxextrema(phase, (int64_t)nim->nvox, &mn, &mx)) { RM_ERR("phase has no finite voxels\n"); goto done; }
		if (!(fabs((double)(mx - mn) - RM_2PI_F64) <= 0.1)) {
			float *rawv = NULL;
			int rnx, rny, rnz, rnv;
			float rmn, rmx, slope, inter;
			if (!is_first_op) {
				RM_ERR("phase rescaling requires -romeo to be the first operation (its input range is %g, not 2*pi); use -no-phase-rescale to keep the current values\n",
					(double)(mx - mn));
				goto done;
			}
			if (!phasefile) {
				RM_ERR("phase rescaling needs the unscaled stored values, which are unavailable for stdin input; use -no-phase-rescale\n");
				goto done;
			}
			if (rm_read_f32(phasefile, RM_RD_RAW, &rawv, &rnx, &rny, &rnz, &rnv)) goto done;
			if ((int64_t)rnx * rny * rnz * rnv != (int64_t)nim->nvox) {
				free(rawv); RM_ERR("phase file no longer matches the working image\n"); goto done;
			}
			if (rm_approxextrema(rawv, (int64_t)nim->nvox, &rmn, &rmx)) { free(rawv); RM_ERR("phase has no finite voxels\n"); goto done; }
			if (fabs((double)(rmx - rmn) - RM_2PI_F64) <= 0.1) { slope = 1.0f; inter = 0.0f; }
			else {
				slope = (float)(RM_2PI_F64 / (double)(rmx - rmn));
				inter = (float)(-RM_PI_F64 - (double)(rmn * slope));
			}
			for (i = 0; i < (int64_t)nim->nvox; i++) phase[i] = rawv[i] * slope + inter;
			free(rawv);
		}
	}
	if (dump) { if (rm_dump(dump, "c_phase_rescaled.f32", phase, sizeof(float), (int64_t)nim->nvox)) goto done; }
	if (o->verbose) fprintf(stderr, " + -romeo: %dx%dx%d, %d echo(es); phase loaded\n", nx, ny, nz, neco);

	/* ---- echo times ---------------------------------------------------------------------- */
	TEs = (double *)malloc((size_t)(neco > 1 ? neco : 1) * sizeof(double));
	if (!TEs) goto done;
	if (o->nTE == 0) {
		if (neco > 1) { RM_ERR("multi-echo data requires echo times: -t '[te1,te2,...]'\n"); goto done; }
		TEs[0] = 1.0;
	} else if (o->te_epi) {
		for (i = 0; i < neco; i++) TEs[i] = o->TEs[0];
	} else if (o->nTE == neco) {
		for (i = 0; i < neco; i++) TEs[i] = o->TEs[i];
	} else if (neco == 1 && o->nTE == 1) {
		TEs[0] = o->TEs[0];
	} else if (neco == 1 && o->nTE > 1) {
		RM_ERR("%d echo times given for a single-volume image; supply one (-t %g) or a 4D phase\n",
			o->nTE, o->TEs[0]);
		goto done;
	} else {
		RM_ERR("%d echo time(s) given for %d echo(es) in the data\n", o->nTE, neco);
		goto done;
	}

	/* ---- magnitude ------------------------------------------------------------------------ */
	if (magfile && strcmp(magfile, "none") != 0) {
		int mnx, mny, mnz;
		if (rm_read_f32(magfile, RM_RD_SCALED, &mag, &mnx, &mny, &mnz, &magvol)) goto done;
		if (mnx != nx || mny != ny || mnz != nz) {
			RM_ERR("magnitude dimensions (%dx%dx%d) do not match the phase (%dx%dx%d)\n", mnx, mny, mnz, nx, ny, nz);
			goto done;
		}
		if (magvol < neco) {
			RM_ERR("magnitude has %d volume(s) but %d echo(es) are unwrapped\n", magvol, neco);
			goto done;
		}
		have_mag = 1;
	}

	/* ---- mask, weights and unwrapping: the shared core -------------------------------------- */
	core.phase = phase;
	core.mag = have_mag ? mag : NULL;
	core.magvol = magvol;
	core.nx = nx; core.ny = ny; core.nz = nz; core.neco = neco;
	core.n3 = n3;
	core.TEs = TEs; core.o = o;
	core.have_mag = have_mag;
	core.dump = dump;
	if (rm_core_run(&core)) { drc |= core.drc; goto done; }
	drc |= core.drc;
	mask = core.mask;

	/* ---- side outputs ------------------------------------------------------------------------ */
	{
		int want_side = (mask && !o->no_mask_out) || o->compute_b0 || o->write_quality ||
			o->write_quality_all || (dump != NULL);
		float *tmp = want_side ? (float *)malloc((size_t)n3 * sizeof(float)) : NULL;
		int save_rc = 0;
		if (want_side && !tmp) goto done;
		if (mask && !o->no_mask_out) {
			for (i = 0; i < n3; i++) tmp[i] = (float)mask[i];
			save_rc |= rm_save_side(nim, "_mask", tmp, n3, gzMode);
		}
		if (o->compute_b0) {
			float *snr = (float *)malloc((size_t)n3 * sizeof(float));
			char pf[64];
			if (!snr) { free(tmp); goto done; }
			/* ROMEO's own multi-echo -B silently enables MCPC-3D-S monopolar phase-offset
			   correction, which is out of scope here (plan section 7). Say so rather than imply
			   full CLI equivalence. */
			if (neco > 1)
				fprintf(stderr, " + -romeo -B: B0 computed WITHOUT MCPC-3D-S phase-offset correction, which ROMEO's own multi-echo -B applies; the maps are comparable to `romeo --compute-B0 --phase-offset-correction off`\n");
			if (!have_mag && neco > 1)
				fprintf(stderr, " + -romeo -B: B0 frequency estimation without magnitude might result in poor handling of noise in later echoes!\n");
			rm_compute_b0(phase, have_mag ? mag : NULL, n3, neco, TEs, o->b0_weighting, tmp, snr);
			snprintf(pf, sizeof pf, "_%s", o->b0_name);
			save_rc |= rm_save_side(nim, pf, tmp, n3, gzMode);
			snprintf(pf, sizeof pf, "_%s_snr", o->b0_name);
			save_rc |= rm_save_side(nim, pf, snr, n3, gzMode);
			if (dump) {
				drc |= rm_dump(dump, "c_b0.f32", tmp, sizeof(float), n3);
				drc |= rm_dump(dump, "c_b0_snr.f32", snr, sizeof(float), n3);
			}
			free(snr);
		}
		if (o->write_quality || o->write_quality_all || dump) {
			rm_wctx c;
			int qi;
			/* write_qualitymap runs AFTER unwrapping, on data["phase"], and voxelquality's own 4D
			   overload defaults p2ref to 2 regardless of `template` (it does NOT use the
			   template-1 rule that unwrap! applies). */
			if (rm_build_ctx(&c, phase, have_mag ? mag : NULL, magvol, mask, core.magmasked, TEs, neco,
					template_echo, 2, nx, ny, nz, core.flags)) { free(tmp); goto done; }
			if (o->write_quality || dump) {
				if (rm_voxelquality(&c, tmp)) { free(tmp); goto done; }
				if (dump) drc |= rm_dump(dump, "c_qmap.f32", tmp, sizeof(float), n3);
				if (o->write_quality) save_rc |= rm_save_side(nim, "_quality", tmp, n3, gzMode);
			}
			if (o->write_quality_all || dump) {
				for (qi = 0; qi < 6; qi++) {
					rm_wctx c1 = c;
					int64_t x, y, z;
					int allone = 1;
					char pf[32], nm[64];
					memset(c1.flags, 0, sizeof c1.flags);
					c1.flags[qi] = 1;
					rm_updateflags(c1.flags, c1.P2 != NULL, 1, c1.M != NULL);
					if (rm_voxelquality(&c1, tmp)) { free(tmp); goto done; }
					if (dump) { snprintf(nm, sizeof nm, "c_qmap_%d.f32", qi + 1); drc |= rm_dump(dump, nm, tmp, sizeof(float), n3); }
					if (!o->write_quality_all) continue;
					/* `all(qm[1:end-1,1:end-1,1:end-1] .== 1.0)` — an empty range makes all() true,
					   so a singleton dimension skips the map, exactly as upstream. */
					for (z = 0; z + 1 < nz && allone; z++) for (y = 0; y + 1 < ny && allone; y++) for (x = 0; x + 1 < nx && allone; x++)
						if (tmp[x + nx * (y + ny * z)] != 1.0f) allone = 0;
					if (allone) continue;
					snprintf(pf, sizeof pf, "_quality_%d", qi + 1);
					save_rc |= rm_save_side(nim, pf, tmp, n3, gzMode);
				}
			}
		}
		free(tmp);
		if (save_rc) { RM_ERR("failed to write a side output\n"); goto done; }
	}

	if (dump) {
		drc |= rm_dump(dump, "c_visited.u8", core.visited, 1, n3);
		drc |= rm_dump(dump, "c_unwrapped.f32", phase, sizeof(float), (int64_t)nim->nvox);
		if (drc) { RM_ERR("one or more -romeo-dump writes failed\n"); goto done; }
	}

	if (o->verbose) fprintf(stderr, " + -romeo: unwrapping finished\n");
	memcpy(nim->data, phase, (size_t)nim->nvox * sizeof(float));
	ret = 0;

done:
	if (core.manifest) fclose(core.manifest);
	rm_core_free(&core);   /* owns mask, weights, magmasked, visited and the mask stages */
	free(phase); free(mag); free(TEs);
	(void)ihdr;
	return ret;
}

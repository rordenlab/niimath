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
 * That is not a theoretical worry — it was MEASURED on the 76x76x46 validation volume
 * (phase0/mag0, -t 16.8), building this same source three ways and comparing against the pinned
 * Julia oracle:
 *
 *   FP policy for romeo.o                          weight bytes differing   unwrapped
 *   ---------------------------------------------  ----------------------   -------------------
 *   -fno-fast-math -ffp-contract=off  (SHIPPED)          0 / 797088         bit-identical
 *   -fno-fast-math -ffp-contract=fast (FMA only)         0 / 797088         bit-identical
 *   -ffast-math -fno-finite-math-only (repo default)   360 / 797088         66 voxels off by
 *                                                                           >=1 full 2*pi wrap,
 *                                                                           max|diff| 12.57 rad
 *
 * The failure mode is NOT float32 rounding: reassociation pushes a weight just past 1.0, the
 * `0 <= w <= 1` guard in rescale() then returns bin 0, and the edge DISAPPEARS from the graph
 * (largest observed bin deviation: 252).  Multi-echo output stops matching even at --compare
 * 1e-4.  Contraction alone measured clean, but it is left off because (a) it buys nothing —
 * median runtime is 0.07 s either way, the op is I/O-bound — and (b) it would silently break the
 * Dekker/Cody-Waite double-double arithmetic in the 2*pi range reduction below, where Julia
 * fuses ONLY at its explicit muladd sites (mirrored here as explicit fma() calls).  -fno-fast-math
 * does not inhibit SIMD auto-vectorisation of the non-reduction loops, so no vectorisation is
 * given up; only reassociation is.  Change ROMEO_STRICT_FP in src/Makefile to revisit.
 *
 * ---------------------------------------------------------------------------------------------
 * NUMERIC TYPE AUDIT (romeo_plan.md §3.1).  Julia's promotion rules are NOT uniform across the
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
#include <float.h>
#include "romeo.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define RM_NBINS 256
#define RM_PI_F32 3.14159274101257324f    /* Float32(pi)  = 0x1.921fb6p+1 */
#define RM_2PI_F32 6.28318548202514648f   /* Float32(2pi) = 0x1.921fb6p+2 */
#define RM_PI_F64 3.14159265358979323846  /* Float64(pi) */
#define RM_2PI_F64 6.28318530717958647692 /* Float64(2pi) */

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
	shift = k - (idx << 6);
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
   PINNED-VERSION GOTCHA: the environment in romeo_plan.md §2.1 resolves the REGISTRY package
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

/* Gather the sample of v[0..n) into out (caller supplies capacity >= len*len). */
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
	if (m == 0) {
		for (b = 0; b < n; b++) if (isfinite(v[b])) out[m++] = v[b];
	}
	return m;
}
static int64_t rm_sample_capacity(int64_t n) { int64_t len = rm_sample_len(n); return len * len > n ? len * len : n; }

/* approxextrema(I) = extrema(sample(I)), falling back to extrema(I) when the sample is flat. */
static int rm_approxextrema(const float *v, int64_t n, float *mn, float *mx) {
	int64_t cap = rm_sample_capacity(n), m, i;
	float *s = (float *)malloc((size_t)cap * sizeof(float));
	float lo, hi;
	if (!s) return 1;
	m = rm_sample_f32(v, n, s);
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
		if (mj < mi) { small = mj; big = mi; } else { small = mi; big = mj; } /* Base.minmax */
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
	int64_t *stack = NULL, *sizes = NULL;
	int64_t nlab = 0, sp;
	if (1.0 > maxhole) {
		RM_ERR("robustmask needs at least 20 voxels (upstream fill_holes passes (1, n/20) to imfill, which rejects n < 20)\n");
		return 1;
	}
	lab = (int32_t *)calloc((size_t)n, sizeof(int32_t));
	stack = (int64_t *)malloc((size_t)n * sizeof(int64_t));
	sizes = (int64_t *)malloc((size_t)(n + 1) * sizeof(int64_t));
	if (!lab || !stack || !sizes) { free(lab); free(stack); free(sizes); return 1; }
	for (i = 0; i < n; i++) {
		if (mask[i] || lab[i]) continue;
		nlab++;
		sizes[nlab] = 0;
		sp = 0; stack[sp++] = i; lab[i] = (int32_t)nlab;
		while (sp > 0) {
			int64_t v = stack[--sp];
			int64_t z = v / ((int64_t)nx * ny), rem = v % ((int64_t)nx * ny);
			int64_t y = rem / nx, x = rem % nx;
			sizes[nlab]++;
			if (x > 0 && !mask[v - 1] && !lab[v - 1]) { lab[v - 1] = (int32_t)nlab; stack[sp++] = v - 1; }
			if (x < nx - 1 && !mask[v + 1] && !lab[v + 1]) { lab[v + 1] = (int32_t)nlab; stack[sp++] = v + 1; }
			if (y > 0 && !mask[v - nx] && !lab[v - nx]) { lab[v - nx] = (int32_t)nlab; stack[sp++] = v - nx; }
			if (y < ny - 1 && !mask[v + nx] && !lab[v + nx]) { lab[v + nx] = (int32_t)nlab; stack[sp++] = v + nx; }
			if (z > 0 && !mask[v - (int64_t)nx * ny] && !lab[v - (int64_t)nx * ny]) { lab[v - (int64_t)nx * ny] = (int32_t)nlab; stack[sp++] = v - (int64_t)nx * ny; }
			if (z < nz - 1 && !mask[v + (int64_t)nx * ny] && !lab[v + (int64_t)nx * ny]) { lab[v + (int64_t)nx * ny] = (int32_t)nlab; stack[sp++] = v + (int64_t)nx * ny; }
		}
	}
	for (i = 0; i < n; i++) {
		if (!mask[i]) {
			int64_t c = sizes[lab[i]];
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

static void rm_mask_stages_free(rm_mask_stages *s) {
	free(s->s1); free(s->s2); free(s->s3); free(s->sm1); free(s->sm2);
	s->s1 = s->s2 = s->s3 = NULL; s->sm1 = s->sm2 = NULL;
	/* s4 is handed to the caller */
}

/* robustmask(weight; factor=1, threshold=nothing).  Returns 0 and fills stages->s4 (owned by
   the caller) on success.  `have_thr` selects the -k qualitymask path. */
static int rm_robustmask(const float *weight, int nx, int ny, int nz,
	int have_thr, double thr_in, rm_mask_stages *st) {
	int64_t n = (int64_t)nx * ny * nz, i, m;
	float *s = NULL;
	float threshold;
	memset(st, 0, sizeof *st);
	if (!have_thr) {
		int64_t cap = rm_sample_capacity(n);
		double q05, q15, q8, q99, acc;
		int64_t cnt;
		float *tmp = NULL;
		s = (float *)malloc((size_t)cap * sizeof(float));
		tmp = (float *)malloc((size_t)cap * sizeof(float));
		if (!s || !tmp) { free(s); free(tmp); return 1; }
		m = rm_sample_f32(weight, n, s);
		if (m < 1) { free(s); free(tmp); RM_ERR("magnitude has no finite voxels\n"); return 1; }
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

	st->s1 = (uint8_t *)malloc((size_t)n);
	st->sm1 = (float *)malloc((size_t)n * sizeof(float));
	st->s2 = (uint8_t *)malloc((size_t)n);
	st->s3 = (uint8_t *)malloc((size_t)n);
	st->sm2 = (float *)malloc((size_t)n * sizeof(float));
	st->s4 = (uint8_t *)malloc((size_t)n);
	if (!st->s1 || !st->sm1 || !st->s2 || !st->s3 || !st->sm2 || !st->s4) {
		rm_mask_stages_free(st); free(st->s4); st->s4 = NULL; return 1;
	}
	for (i = 0; i < n; i++) st->s1[i] = (weight[i] > threshold) ? 1 : 0;
	for (i = 0; i < n; i++) st->sm1[i] = (float)st->s1[i];
	{
		int boxes1[1] = { 5 };
		if (rm_boxsmooth3d(st->sm1, nx, ny, nz, 1, boxes1)) { rm_mask_stages_free(st); free(st->s4); st->s4 = NULL; return 1; }
	}
	for (i = 0; i < n; i++) st->s2[i] = ((double)st->sm1[i] > 0.4) ? 1 : 0;
	memcpy(st->s3, st->s2, (size_t)n);
	if (rm_fill_holes(st->s3, nx, ny, nz)) { rm_mask_stages_free(st); free(st->s4); st->s4 = NULL; return 1; }
	for (i = 0; i < n; i++) st->sm2[i] = (float)st->s3[i];
	{
		int boxes2[2] = { 3, 3 };
		if (rm_boxsmooth3d(st->sm2, nx, ny, nz, 2, boxes2)) { rm_mask_stages_free(st); free(st->s4); st->s4 = NULL; return 1; }
	}
	for (i = 0; i < n; i++) st->s4[i] = ((double)st->sm2[i] > 0.6) ? 1 : 0;
	return 0;
}

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
} rm_pq;

static int rm_pq_init(rm_pq *q, int nbins) {
	q->nbins = nbins;
	q->min = nbins + 1;
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
	if (w < 1 || w > q->nbins) return 1;
	if (q->len[w] == q->cap[w]) {
		int64_t nc = q->cap[w] ? q->cap[w] * 2 : 64;
		int64_t *nb = (int64_t *)realloc(q->bin[w], (size_t)nc * sizeof(int64_t));
		if (!nb) return 1;
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

/* The seed queue is built once from sum(weights; dims=1) in ascending linear index and only
   ever dequeued, so a counting sort reproduces the bucket layout exactly: within a bin the
   entries are ascending and pop-from-the-end yields the HIGHEST linear index first. */
typedef struct {
	int64_t *items;   /* concatenated bins */
	int64_t *start;   /* bin offsets, 1..nbins+1 */
	int64_t *len;
	int nbins;
	int min;
} rm_seedq;

static void rm_seedq_free(rm_seedq *s) { free(s->items); free(s->start); free(s->len); memset(s, 0, sizeof *s); }

static int rm_seedq_build(rm_seedq *s, const uint8_t *w, int64_t n) {
	int64_t i;
	int b;
	int64_t *fill = NULL;
	s->nbins = 3 * RM_NBINS;
	s->items = (int64_t *)malloc((size_t)n * sizeof(int64_t));
	s->start = (int64_t *)calloc((size_t)s->nbins + 2, sizeof(int64_t));
	s->len = (int64_t *)calloc((size_t)s->nbins + 2, sizeof(int64_t));
	fill = (int64_t *)calloc((size_t)s->nbins + 2, sizeof(int64_t));
	if (!s->items || !s->start || !s->len || !fill) { free(fill); rm_seedq_free(s); return 1; }
	for (i = 0; i < n; i++) {
		int a = w[3 * i] ? w[3 * i] : 255;
		int b1 = w[3 * i + 1] ? w[3 * i + 1] : 255;
		int c = w[3 * i + 2] ? w[3 * i + 2] : 255;
		s->len[a + b1 + c]++;
	}
	s->start[1] = 0;
	for (b = 1; b <= s->nbins; b++) s->start[b + 1] = s->start[b] + s->len[b];
	for (i = 0; i < n; i++) {
		int a = w[3 * i] ? w[3 * i] : 255;
		int b1 = w[3 * i + 1] ? w[3 * i + 1] : 255;
		int c = w[3 * i + 2] ? w[3 * i + 2] : 255;
		int bb = a + b1 + c;
		s->items[s->start[bb] + fill[bb]++] = i + 1;   /* 1-based voxel index, ascending */
	}
	free(fill);
	s->min = s->nbins + 1;
	for (b = 1; b <= s->nbins; b++) if (s->len[b]) { s->min = b; break; }
	return 0;
}

static int64_t rm_seedq_dequeue(rm_seedq *s) {
	int64_t e = s->items[s->start[s->min] + (--s->len[s->min])];
	while (s->min <= s->nbins && s->len[s->min] == 0) s->min++;
	return e;
}
static int rm_seedq_isempty(const rm_seedq *s) { return s->min > s->nbins; }

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

/* Returns the new seed threshold, or 255 when no unvisited voxel remains. */
static double rm_addseed(rm_grow *g, rm_seedq *sq, rm_pq *pq, int64_t *seeds, int *nseeds) {
	int64_t seed = 0;
	int i;
	while (!rm_seedq_isempty(sq)) {
		int64_t ind = rm_seedq_dequeue(sq);
		if (g->visited[ind - 1] == 0) { seed = ind; break; }
	}
	if (seed == 0) return 255.0;
	for (i = 1; i <= 6; i++) {
		int64_t e = rm_getnewedge(g, seed, i);
		if (e != 0 && g->weights[e - 1] > 0) rm_pq_enqueue(pq, e, g->weights[e - 1]);
	}
	rm_seedcorrection(g, seed);
	seeds[*nseeds] = seed;
	(*nseeds)++;
	g->visited[seed - 1] = (uint8_t)(*nseeds);
	{
		int sum = (int)g->weights[rm_getedgeindex(seed, 1) - 1]
			+ (int)g->weights[rm_getedgeindex(seed, 2) - 1]
			+ (int)g->weights[rm_getedgeindex(seed, 3) - 1];
		/* NBINS - div(NBINS - sum/3, 2): sum/3 is Float64, div(::Float64,2) truncates toward 0 */
		double t = (double)RM_NBINS - (double)sum / 3.0;
		return (double)RM_NBINS - trunc(t / 2.0);
	}
}

/* grow_region_unwrap!.  maxseeds is capped at 255 upstream; only 1 is supported here (the
   experimental multi-seed/region-merging path is not ported). `pq` may already hold seed edges
   (the temporal-uncertain re-entry), in which case no seed is created. */
static int rm_grow_region(rm_grow *g, rm_pq *pq, rm_seedq *sq, int maxseeds, uint8_t *out_visited) {
	int64_t seeds[256];
	int nseeds = 0;
	double new_seed_thresh = 256.0;
	int seeded = 0;
	if (rm_pq_isempty(pq)) {
		if (!sq) return 1;
		new_seed_thresh = rm_addseed(g, sq, pq, seeds, &nseeds);
		seeded = 1;
	}
	while (!rm_pq_isempty(pq)) {
		int64_t edge, oldvox, newvox, vox, neighbor;
		int dim, i;
		if (seeded && nseeds < maxseeds && (double)pq->min > new_seed_thresh)
			new_seed_thresh = rm_addseed(g, sq, pq, seeds, &nseeds);
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
		}
	}
	if (out_visited && out_visited != g->visited) memcpy(out_visited, g->visited, (size_t)g->n);
	return 0;
}

/* ============================================================================================
 * 9. unwrap! (ROMEO.jl src/unwrapping.jl)
 * ==========================================================================================*/

/* correctglobal: wrapped .-= 2pi*median(round.(filter(isfinite, wrapped[mask]) ./ 2pi)) */
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
	rm_seedq sq;
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
	if (rm_seedq_build(&sq, weights, n)) { rm_pq_free(&pq); free(visited); return 1; }
	rc = rm_grow_region(&g, &pq, &sq, maxseeds, NULL);
	rm_pq_free(&pq);
	rm_seedq_free(&sq);
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
	snprintf(path, sizeof path, "%s/%s", dir, name);
	f = fopen(path, "wb");
	if (!f) { RM_ERR("cannot write dump '%s'\n", path); return 1; }
	if (n > 0 && fwrite(p, elem, (size_t)n, f) != (size_t)n) { fclose(f); RM_ERR("short write to '%s'\n", path); return 1; }
	fclose(f);
	return 0;
}

romeo_opts romeo_opts_default(void) {
	romeo_opts o;
	memset(&o, 0, sizeof o);
	o.weights_sel = RM_W_ROMEO;
	o.mask_sel = RM_MASK_ROBUST;
	o.template_echo = 1;
	o.maxseeds = 1;
	o.qmask_thresh = 0.1;
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
	(void)cmd;
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
					if (!rm_parse_double(argv[ac], &v)) { o->qmask_thresh = v; o->qmask_thresh_set = 1; ac++; }
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
			if (ac + 1 >= argc || rm_parse_int(argv[ac + 1], &v) || v < 1) { RM_ERR("-template requires a positive integer\n"); return 1; }
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
			!strcmp(a, "-fix-ge-phase") || !strcmp(a, "-B") || !strcmp(a, "-B0-phase-weighting")) {
			RM_ERR("'%s' is a ROMEO option that is not implemented in this build\n", a);
			return 1;
		}
		break; /* first unrecognized token: back off so niimath sees the next chain operation */
	}
	*pac = ac;
	return 0;
}

/* ---- auxiliary image loading -------------------------------------------------------------- */

/* Read a NIfTI into float32.  raw!=0 returns the STORED values (Julia's `.raw`, used for masks
   and for readphase's second branch); otherwise slope/intercept are applied in FLOAT arithmetic,
   matching NIfTI.jl's getindex (raw*scl_slope + scl_inter, all Float32). */
static int rm_read_f32(const char *fn, int raw, float **out, int *nx, int *ny, int *nz, int *nvol) {
	nifti_image *n = nifti_image_read(fn, 1);
	int64_t i, nv;
	float *d = NULL;
	float scl, inter;
	if (!n) { RM_ERR("failed to read '%s'\n", fn); return 1; }
	if (n->nx < 1 || n->ny < 1 || n->nz < 1 || n->nvox < 1) { nifti_image_free(n); RM_ERR("'%s' has invalid dimensions\n", fn); return 1; }
	nv = (int64_t)n->nvox;
	if (nv > INT_MAX) { nifti_image_free(n); RM_ERR("'%s' exceeds INT_MAX voxels\n", fn); return 1; }
	d = (float *)malloc((size_t)nv * sizeof(float));
	if (!d) { nifti_image_free(n); return 1; }
	scl = (n->scl_slope == 0.0f) ? 1.0f : n->scl_slope;
	inter = n->scl_inter;
	if (raw) { scl = 1.0f; inter = 0.0f; }
#define RM_CVT(T) do { const T *p = (const T *)n->data; for (i = 0; i < nv; i++) d[i] = (float)p[i] * scl + inter; } while (0)
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
	default:
		free(d); nifti_image_free(n);
		RM_ERR("'%s' has an unsupported datatype (%d)\n", fn, n->datatype);
		return 1;
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

/* ---- the runner ---------------------------------------------------------------------------- */

int romeo_run(nifti_image *nim, const char *magfile, const char *phasefile,
	const in_hdr *ihdr, int is_first_op, const romeo_opts *o, gzModes gzMode) {
	int nx, ny, nz, neco, ret = 1;
	int64_t n3, i;
	float *phase = NULL;     /* owned working copy, [n3 * neco] */
	float *mag = NULL;       /* owned, [n3 * magvol] or NULL */
	int magvol = 0;
	uint8_t *mask = NULL;    /* owned or NULL */
	uint8_t *weights = NULL;
	float *magmasked = NULL;
	uint8_t *visited = NULL;
	rm_mask_stages stages;
	int flags[6];
	int template_echo, p2ref = 0;
	double TE1 = 1.0, TE2 = 1.0;
	double *TEs = NULL;
	int have_mag = 0;
	const char *dump = o->dumpdir;
	FILE *manifest = NULL;

	memset(&stages, 0, sizeof stages);

	if (nim->datatype != DT_FLOAT32) { RM_ERR("internal error: expected float32 working image\n"); return 1; }
	if (nim->nu > 1 || nim->nv > 1 || nim->nw > 1) { RM_ERR("input must be 3D or 4D (echoes on dim 4)\n"); return 1; }
	nx = nim->nx; ny = nim->ny; nz = (nim->nz < 1 ? 1 : nim->nz);
	n3 = (int64_t)nx * ny * nz;
	if (n3 < 1 || (int64_t)nim->nvox % n3 != 0) { RM_ERR("invalid image geometry\n"); return 1; }
	neco = (int)((int64_t)nim->nvox / n3);
	if (neco < 1) { RM_ERR("invalid image geometry\n"); return 1; }
	if ((int64_t)nim->nvox > INT_MAX) { RM_ERR("images with more than INT_MAX voxels are not supported\n"); return 1; }

	template_echo = o->template_echo;
	if (template_echo > neco) { RM_ERR("-template %d exceeds the %d echo(es) present\n", template_echo, neco); return 1; }

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
			if (rm_read_f32(phasefile, 1, &rawv, &rnx, &rny, &rnz, &rnv)) goto done;
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
	} else if (neco == 1 && o->nTE >= 1) {
		TEs[0] = o->TEs[0];
	} else {
		RM_ERR("%d echo time(s) given for %d echo(es) in the data\n", o->nTE, neco);
		goto done;
	}

	/* ---- magnitude ------------------------------------------------------------------------ */
	if (magfile && strcmp(magfile, "none") != 0) {
		int mnx, mny, mnz;
		if (rm_read_f32(magfile, 0, &mag, &mnx, &mny, &mnz, &magvol)) goto done;
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

	/* ---- weight selection (needed before the mask: -k qualitymask uses it) ------------------ */
	rm_flags_from_sel(o->weights_sel, have_mag, flags);
	if (o->weights_sel == RM_W_FLAGS) for (i = 0; i < 6; i++) flags[i] = o->flags[i];
	if (neco > 1) {
		p2ref = (template_echo == 1) ? 2 : template_echo - 1;
		TE1 = TEs[template_echo - 1];
		TE2 = TEs[p2ref - 1];
	}
	if (have_mag) {
		magmasked = (float *)malloc((size_t)n3 * sizeof(float));
		if (!magmasked) goto done;
	}

	/* ---- mask ------------------------------------------------------------------------------ */
	if (o->mask_sel == RM_MASK_ROBUST && !have_mag) {
		/* load_data_and_resolve_args!: robustmask without a magnitude degrades to nomask */
		fprintf(stderr, " + -romeo: robustmask was chosen but no magnitude is available. No mask is used!\n");
	} else if (o->mask_sel == RM_MASK_ROBUST) {
		int te = template_echo < magvol ? template_echo : magvol;
		if (rm_robustmask(mag + (int64_t)(te - 1) * n3, nx, ny, nz, 0, 0.0, &stages)) goto done;
		mask = stages.s4; stages.s4 = NULL;
	} else if (o->mask_sel == RM_MASK_QUALITY) {
		/* set_mask!: qmap = voxelquality(phase; get_keyargs(...)) — computed on the still-WRAPPED
		   phase and WITHOUT a mask (data["mask"] does not exist yet), then robustmask(qmap; threshold).
		   voxelquality's own 4D overload defaults p2ref to 2 regardless of `template`. */
		rm_wctx qc;
		float *qmap = (float *)malloc((size_t)n3 * sizeof(float));
		if (!qmap) goto done;
		if (rm_build_ctx(&qc, phase, have_mag ? mag : NULL, magvol, NULL, magmasked, TEs, neco,
				template_echo, 2, nx, ny, nz, flags)) { free(qmap); goto done; }
		if (rm_voxelquality(&qc, qmap)) { free(qmap); goto done; }
		if (dump) rm_dump(dump, "c_qmap_wrapped.f32", qmap, sizeof(float), n3);
		if (rm_robustmask(qmap, nx, ny, nz, 1, o->qmask_thresh, &stages)) { free(qmap); goto done; }
		free(qmap);
		mask = stages.s4; stages.s4 = NULL;
	} else if (o->mask_sel == RM_MASK_FILE) {
		float *mv = NULL;
		int mnx, mny, mnz, mnv;
		if (rm_read_f32(o->mask_file, 1, &mv, &mnx, &mny, &mnz, &mnv)) goto done;
		if (mnx != nx || mny != ny || mnz != nz || mnv != 1) {
			free(mv); RM_ERR("mask dimensions do not match the phase\n"); goto done;
		}
		mask = (uint8_t *)malloc((size_t)n3);
		if (!mask) { free(mv); goto done; }
		for (i = 0; i < n3; i++) mask[i] = (mv[i] != 0.0f) ? 1 : 0; /* raw stored values, per niread(...).raw .!= 0 */
		free(mv);
	}

	/* ---- weights ---------------------------------------------------------------------------- */
	{
		rm_wctx c;
		if (rm_build_ctx(&c, phase, have_mag ? mag : NULL, magvol, mask, magmasked, TEs, neco,
				template_echo, neco > 1 ? p2ref : 1, nx, ny, nz, flags)) goto done;
		weights = (uint8_t *)malloc((size_t)3 * (size_t)n3);
		if (!weights) goto done;
		rm_calculateweights(&c, RM_WOUT_U8, weights);

		if (dump) {
			char path[2048];
			snprintf(path, sizeof path, "%s/c_manifest.txt", dump);
			manifest = fopen(path, "w");
			if (rm_dump(dump, "c_weights.u8", weights, 1, 3 * n3)) goto done;
			{
				double *wd = (double *)malloc((size_t)3 * (size_t)n3 * sizeof(double));
				if (wd) {
					rm_calculateweights(&c, RM_WOUT_F64, wd);
					rm_dump(dump, "c_weights_prerescale.f64", wd, sizeof(double), 3 * n3);
					free(wd);
				}
			}
			if (manifest) {
				fprintf(manifest, "flags_active %d%d%d%d%d%d\n",
					c.flags[0], c.flags[1], c.flags[2], c.flags[3], c.flags[4], c.flags[5]);
				if (have_mag) fprintf(manifest, "maxmag %.17g\n", c.maxmag);
				if (stages.s1) {
					fprintf(manifest, "rm_sample_len %lld\n", (long long)stages.sample_len);
					fprintf(manifest, "rm_q05 %.17g\nrm_q15 %.17g\nrm_q8 %.17g\nrm_q99 %.17g\n",
						stages.q05, stages.q15, stages.q8, stages.q99);
					fprintf(manifest, "rm_high_intensity %.9g\nrm_noise %.9g\nrm_noise_stage %d\nrm_threshold %.9g\n",
						(double)stages.high_intensity, (double)stages.noise, stages.noise_stage, (double)stages.threshold);
				}
			}
			if (stages.s1) {
				rm_dump(dump, "c_mask_s1_thresh.u8", stages.s1, 1, n3);
				rm_dump(dump, "c_mask_sm1.f32", stages.sm1, sizeof(float), n3);
				rm_dump(dump, "c_mask_s2_smooth1.u8", stages.s2, 1, n3);
				rm_dump(dump, "c_mask_s3_fill.u8", stages.s3, 1, n3);
				rm_dump(dump, "c_mask_sm2.f32", stages.sm2, sizeof(float), n3);
			}
			if (mask) rm_dump(dump, "c_mask_s4_final.u8", mask, 1, n3);
		}
	}

	/* ---- unwrap ----------------------------------------------------------------------------- */
	visited = (uint8_t *)calloc((size_t)n3, 1);
	if (!visited) goto done;
	if (neco == 1) {
		if (rm_unwrap3d(phase, weights, nx, ny, nz, NULL, TE1, TE2, 0,
				o->wrap_addition, o->maxseeds, visited)) goto done;
		if (o->correctglobal && rm_correctglobal(phase, n3, mask)) goto done;
	} else if (o->individual) {
		/* unwrap_individual!: each echo is unwrapped spatially with its own weights, using the
		   PREVIOUS echo (echo 2 for echo 1) as phase2.  Echoes are processed in ascending order,
		   so when echo i>1 is unwrapped its reference echo i-1 is ALREADY unwrapped — that is
		   upstream behaviour (Threads.@threads with a shared Dict; the oracle pins 1 thread). */
		int ie;
		for (ie = 1; ie <= neco; ie++) {
			int e2 = (ie == 1) ? 2 : ie - 1;
			rm_wctx c;
			float *p2copy = (float *)malloc((size_t)n3 * sizeof(float));
			uint8_t *w2 = NULL;
			if (!p2copy) goto done;
			memcpy(p2copy, phase + (int64_t)(e2 - 1) * n3, (size_t)n3 * sizeof(float));
			if (rm_build_ctx(&c, phase, have_mag ? mag : NULL, magvol, mask, magmasked, TEs, neco,
					ie, e2, nx, ny, nz, flags)) { free(p2copy); goto done; }
			w2 = (uint8_t *)malloc((size_t)3 * (size_t)n3);
			if (!w2) { free(p2copy); goto done; }
			rm_calculateweights(&c, RM_WOUT_U8, w2);
			if (rm_unwrap3d(phase + (int64_t)(ie - 1) * n3, w2, nx, ny, nz, p2copy,
					c.TE1, c.TE2, 1, o->wrap_addition, o->maxseeds, visited)) {
				free(p2copy); free(w2); goto done;
			}
			if (o->correctglobal && rm_correctglobal(phase + (int64_t)(ie - 1) * n3, n3, mask)) {
				free(p2copy); free(w2); goto done;
			}
			free(p2copy); free(w2);
		}
		if (o->correctglobal) {
			/* correct_multi_echo_wraps! */
			int ie2;
			double *v = (double *)malloc((size_t)n3 * sizeof(double));
			if (!v) goto done;
			for (ie2 = 2; ie2 <= neco; ie2++) {
				int iref = ie2 - 1;
				double fac = TEs[ie2 - 1] / TEs[iref - 1], nwraps;
				const float *pr = phase + (int64_t)(iref - 1) * n3;
				float *pe = phase + (int64_t)(ie2 - 1) * n3;
				int64_t m = 0;
				for (i = 0; i < n3; i++) {
					if (mask && !mask[i]) continue;
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
		float *tpl = phase + (int64_t)(template_echo - 1) * n3;
		const float *p2 = phase + (int64_t)(p2ref - 1) * n3;
		float *p2copy = (float *)malloc((size_t)n3 * sizeof(float));
		int order_i;
		if (!p2copy) goto done;
		memcpy(p2copy, p2, (size_t)n3 * sizeof(float)); /* args[:phase2] is a COPY taken up front */
		if (rm_unwrap3d(tpl, weights, nx, ny, nz, p2copy, TE1, TE2, 1,
				o->wrap_addition, o->maxseeds, visited)) { free(p2copy); goto done; }
		free(p2copy);
		if (o->correctglobal && rm_correctglobal(tpl, n3, mask)) goto done;
		for (order_i = 0; order_i < neco - 1; order_i++) {
			/* iteration order: (template-1):-1:1, then (template+1):neco */
			int ieco = (order_i < template_echo - 1) ? (template_echo - 1 - order_i) : (order_i + 2);
			int iref = (ieco < template_echo) ? ieco + 1 : ieco - 1;
			double fac = TEs[ieco - 1] / TEs[iref - 1];
			float *w = phase + (int64_t)(ieco - 1) * n3;
			const float *r = phase + (int64_t)(iref - 1) * n3;
			double *refvalue = (double *)malloc((size_t)n3 * sizeof(double));
			if (!refvalue) goto done;
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
				if (!qual || !halfw || !halfr || !vis) { free(qual); free(halfw); free(halfr); free(vis); free(refvalue); goto done; }
				for (i = 0; i < n3; i++) { halfw[i] = w[i] / 2.0f; halfr[i] = refvalue[i] / 2.0; }
				memset(&qc, 0, sizeof qc);
				qc.P = halfw; qc.P2d = halfr; qc.TE1 = 1.0; qc.TE2 = 1.0;
				qc.nx = nx; qc.ny = ny; qc.nz = nz; qc.n = n3;
				qc.flags[0] = 1; qc.flags[1] = 1; qc.flags[2] = 1; /* :romeo, no mag -> 4..6 off */
				if (rm_voxelquality(&qc, qual)) { free(qual); free(halfw); free(halfr); free(vis); free(refvalue); goto done; }
				for (i = 0; i < n3; i++) vis[i] = ((double)qual[i] > o->temporal_uncertain) ? 1 : 0;
				for (i = 0; i < n3; i++) {
					int inmask;
					if (mask) inmask = mask[i] != 0;
					else {
						int64_t s = (int64_t)weights[3 * i] + weights[3 * i + 1] + weights[3 * i + 2];
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
					g.wrapped = w; g.weights = weights; g.visited = vis; g.n = n3;
					g.stride[0] = stride[0]; g.stride[1] = stride[1]; g.stride[2] = stride[2];
					g.wrap_addition = o->wrap_addition;
					g.phase2 = NULL; g.TE1 = TE1; g.TE2 = TE2; g.have_p2 = 0;
					if (rm_pq_init(&pq, RM_NBINS)) { free(qual); free(halfw); free(halfr); free(vis); free(refvalue); goto done; }
					for (dim = 1; dim <= 3; dim++) {
						int64_t I;
						for (I = 1; I <= n3; I++) {
							int64_t J = I + stride[dim - 1];
							if (J > n3) continue;
							if ((int)vis[I - 1] + (int)vis[J - 1] == 1) {
								int64_t ed = rm_getedgeindex(I, dim);
								if (weights[ed - 1] != 0) rm_pq_enqueue(&pq, ed, weights[ed - 1]);
							}
						}
					}
					rm_grow_region(&g, &pq, NULL, o->maxseeds, NULL);
					rm_pq_free(&pq);
				}
				free(qual); free(halfw); free(halfr); free(vis);
			}
			free(refvalue);
		}
	}

	/* ---- side outputs ------------------------------------------------------------------------ */
	{
		float *tmp = (float *)malloc((size_t)n3 * sizeof(float));
		int save_rc = 0;
		if (!tmp) goto done;
		if (mask && !o->no_mask_out) {
			for (i = 0; i < n3; i++) tmp[i] = (float)mask[i];
			save_rc |= rm_save_side(nim, "_mask", tmp, n3, gzMode);
		}
		if (o->write_quality || o->write_quality_all || dump) {
			rm_wctx c;
			int qi;
			/* write_qualitymap runs AFTER unwrapping, on data["phase"], and voxelquality's own 4D
			   overload defaults p2ref to 2 regardless of `template` (it does NOT use the
			   template-1 rule that unwrap! applies). */
			if (rm_build_ctx(&c, phase, have_mag ? mag : NULL, magvol, mask, magmasked, TEs, neco,
					template_echo, 2, nx, ny, nz, flags)) { free(tmp); goto done; }
			if (o->write_quality || dump) {
				if (rm_voxelquality(&c, tmp)) { free(tmp); goto done; }
				if (dump) rm_dump(dump, "c_qmap.f32", tmp, sizeof(float), n3);
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
					if (dump) { snprintf(nm, sizeof nm, "c_qmap_%d.f32", qi + 1); rm_dump(dump, nm, tmp, sizeof(float), n3); }
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
		rm_dump(dump, "c_visited.u8", visited, 1, n3);
		rm_dump(dump, "c_unwrapped.f32", phase, sizeof(float), (int64_t)nim->nvox);
	}

	memcpy(nim->data, phase, (size_t)nim->nvox * sizeof(float));
	ret = 0;

done:
	if (manifest) fclose(manifest);
	rm_mask_stages_free(&stages);
	free(stages.s4);
	free(phase); free(mag); free(mask); free(weights); free(magmasked); free(visited); free(TEs);
	(void)ihdr;
	return ret;
}

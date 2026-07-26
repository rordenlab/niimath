// stc.c - slice-time correction for 4D datasets (-stc)
//
// Clean-room implementation of AFNI 3dTshift's default Fourier method.  The published help
// supplies the method name and the "detrend -> interpolate -> retrend" outline; AFNI's GPL-2
// 3dTshift.c, its shifting engine and its FFT were NOT read, translated or paraphrased and
// served only as a black-box oracle.  See stc.h and the moco_bench repository's test/stc_reference_manifest.md, which
// records the experiment behind every convention marked "measured" below.
//
// The measured algorithm, per slice z and per voxel time series x[0..nt-1]:
//   s = (t[z] - tzero) / TR                                       (fractional sample shift)
//   |s| < STC_SKIP  ->  the slice is copied verbatim              (measured)
//   otherwise:
//     1. remove the least-squares line a + b*i fitted over i = 0..nt-1;
//     2. zero-pad the residual to N = stc_fft_len(nt);
//     3. multiply DFT bin k by exp(-2*pi*i*k'*s/N), k' folded to [-N/2, N/2), with the
//        Nyquist bin scaled by the REAL factor cos(pi*s);
//     4. inverse transform, keep the first nt real samples;
//     5. clip to the residual's own [min, max]                    (measured)
//     6. add the line back;
//     7. clip to the ORIGINAL series' [min, max]                  (measured)
// Both clips are required and neither subsumes the other; see the manifest.
//
// Performance shape.  The transform is written batched, because the workload is millions of
// independent transforms of identical length rather than one big transform: the lane index is
// the FASTEST axis (SoA), so plain scalar C vectorises.  On top of that, each lane carries TWO
// real series packed as (real, imaginary) of one complex transform -- legitimate here only
// because the frequency-domain filter is exactly Hermitian, so the inverse transform's real and
// imaginary parts are the two shifted series.  Together one batched transform corrects
// 2 * STC_LANES voxels.

#define _USE_MATH_DEFINES // microsoft compiler
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "stc.h"
#include "print.h"

// Measured: AFNI leaves a slice untouched when the fractional shift is below 0.001.
#define STC_SKIP 0.001
// Measured: AFNI's default Fourier path rejects series shorter than 5 samples.
#define STC_MIN_NT 5
// Voxel pairs transformed together.  Lane index is the fastest array axis.
#define STC_LANES 4
#define STC_PAIRS (2 * STC_LANES)
// Bound on the transform length, so the twiddle table and scratch stay sane on wasm32.
#define STC_MAX_N 4194304

// ------------------------------------------------------------------ transform length
// Measured: N is the smallest EVEN integer of the form 2^p * 3^q * 5^r with q,r in {0,1}
// satisfying N >= nt + 4.  The "even" qualifier is only visible at nt = 9..11, where AFNI
// returns 16 rather than 15; the factor 15 itself is real (nt = 476 -> 480 = 15 * 32).
static int64_t stc_fft_len(int64_t nt) {
	const int64_t want = nt + 4;
	int64_t best = 0;
	static const int64_t base[4] = {1, 3, 5, 15};
	for (int b = 0; b < 4; b++) {
		int64_t n = base[b];
		/* Test the candidate BEFORE doubling: checking `n > STC_MAX_N` afterwards would let a
		   candidate sitting exactly on the cap double past it and be returned. */
		while (n < want || (n & 1)) {
			if (n > STC_MAX_N / 2) { n = 0; break; }
			n *= 2;
		}
		if (n > STC_MAX_N) n = 0;
		if (n && (best == 0 || n < best)) best = n;
	}
	return best;
}

// Factor into the radices the butterflies below implement.  Returns the count, or 0 if the
// length is not of the supported form.
static int stc_factor(int64_t n, int *rad, int maxrad) {
	int k = 0;
	static const int r5[3] = {5, 3, 2};
	for (int i = 0; i < 3; i++) {
		while (n % r5[i] == 0) {
			if (k >= maxrad) return 0;
			rad[k++] = r5[i];
			n /= r5[i];
		}
	}
	return (n == 1) ? k : 0;
}

// ------------------------------------------------------------------ batched Stockham FFT
// PROVENANCE: this kernel is original niimath code.  No FFT implementation was read or adapted
// for it -- not kissfft, muFFT or pffft (evaluated as dependencies and rejected: none can
// express AFNI's length set with the batching this workload wants), and certainly not AFNI's own
// GPL-2 csfft.  What IS borrowed is the ALGORITHM, which is classical published mathematics:
// Stockham's autosort FFT (Stockham 1966; the mixed-radix index map is textbook, e.g. Van Loan,
// "Computational Frameworks for the FFT", SIAM 1992), plus the standard radix-3 and radix-5
// butterflies, which are just the roots-of-unity identities written out.  The index map below
// was re-derived from the Cooley-Tukey identity and validated against a naive long-double DFT at
// every supported length before this code shipped.
//
// Stockham autosort: no bit/digit-reversal pass, natural order in and out, at the cost of a
// ping-pong buffer.  Derivation of the index map (L = radices already applied, M = N / L):
//   A_{L'}[j*M' + m + M*q] = sum_a  A_L[(j + L'*a)*M + m] * w_{M'}^{a*m} * w_r^{a*q}
// with L' = L/r and M' = M*r, starting from L = N, M = 1 (the input in natural order).
//
// Every element is a vector of STC_LANES doubles; `sgn` is -1 for the forward transform and +1
// for the inverse (which is left UNNORMALISED -- the caller divides by N).
static void stc_fft(double *ar, double *ai, double *br, double *bi, int64_t n,
                    const int *rad, int nrad, const double *wr, const double *wi, double sgn,
                    double **outr, double **outi) {
	const double S3 = 0.86602540378443864676;   // sin(2*pi/3)
	const double C51 = 0.30901699437494742410;  // cos(2*pi/5)
	const double S51 = 0.95105651629515357212;  // sin(2*pi/5)
	const double C52 = -0.80901699437494742410; // cos(4*pi/5)
	const double S52 = 0.58778525229247312917;  // sin(4*pi/5)
	int64_t L = n, M = 1;
	for (int st = 0; st < nrad; st++) {
		const int r = rad[st];
		const int64_t Lp = L / r, Mp = M * r, tw = n / Mp;
		for (int64_t jp = 0; jp < Lp; jp++) {
			for (int64_t m = 0; m < M; m++) {
				double ur[5][STC_LANES], ui[5][STC_LANES];
				for (int a = 0; a < r; a++) {
					const double *sr = ar + ((jp + Lp * a) * M + m) * STC_LANES;
					const double *si = ai + ((jp + Lp * a) * M + m) * STC_LANES;
					if (a == 0 || m == 0) {
						for (int l = 0; l < STC_LANES; l++) { ur[a][l] = sr[l]; ui[a][l] = si[l]; }
					} else {
						// twiddle exp(sgn * 2*pi*i * a*m / Mp); the table holds
						// (cos, sin)(2*pi*j/n) so the sign is one multiply on the imaginary part
						const int64_t idx = (a * m * tw) % n;
						const double cw = wr[idx], sw = sgn * wi[idx];
						for (int l = 0; l < STC_LANES; l++) {
							ur[a][l] = sr[l] * cw - si[l] * sw;
							ui[a][l] = sr[l] * sw + si[l] * cw;
						}
					}
				}
				double *d0 = br + (jp * Mp + m) * STC_LANES;
				double *e0 = bi + (jp * Mp + m) * STC_LANES;
				const int64_t stride = M * STC_LANES;
				if (r == 2) {
					for (int l = 0; l < STC_LANES; l++) {
						d0[l] = ur[0][l] + ur[1][l];
						e0[l] = ui[0][l] + ui[1][l];
						d0[stride + l] = ur[0][l] - ur[1][l];
						e0[stride + l] = ui[0][l] - ui[1][l];
					}
				} else if (r == 3) {
					const double sv = sgn * S3;
					for (int l = 0; l < STC_LANES; l++) {
						const double pr = ur[1][l] + ur[2][l], pi = ui[1][l] + ui[2][l];
						const double mr = ur[1][l] - ur[2][l], mi = ui[1][l] - ui[2][l];
						const double cr = ur[0][l] - 0.5 * pr, ci = ui[0][l] - 0.5 * pi;
						d0[l] = ur[0][l] + pr;
						e0[l] = ui[0][l] + pi;
						d0[stride + l] = cr - sv * mi;
						e0[stride + l] = ci + sv * mr;
						d0[2 * stride + l] = cr + sv * mi;
						e0[2 * stride + l] = ci - sv * mr;
					}
				} else { // r == 5
					for (int l = 0; l < STC_LANES; l++) {
						const double t1r = ur[1][l] + ur[4][l], t1i = ui[1][l] + ui[4][l];
						const double t2r = ur[2][l] + ur[3][l], t2i = ui[2][l] + ui[3][l];
						const double t3r = ur[1][l] - ur[4][l], t3i = ui[1][l] - ui[4][l];
						const double t4r = ur[2][l] - ur[3][l], t4i = ui[2][l] - ui[3][l];
						const double a1r = ur[0][l] + C51 * t1r + C52 * t2r;
						const double a1i = ui[0][l] + C51 * t1i + C52 * t2i;
						const double b1r = sgn * (S51 * t3r + S52 * t4r);
						const double b1i = sgn * (S51 * t3i + S52 * t4i);
						const double a2r = ur[0][l] + C52 * t1r + C51 * t2r;
						const double a2i = ui[0][l] + C52 * t1i + C51 * t2i;
						const double b2r = sgn * (S52 * t3r - S51 * t4r);
						const double b2i = sgn * (S52 * t3i - S51 * t4i);
						d0[l] = ur[0][l] + t1r + t2r;
						e0[l] = ui[0][l] + t1i + t2i;
						d0[stride + l] = a1r - b1i;
						e0[stride + l] = a1i + b1r;
						d0[2 * stride + l] = a2r - b2i;
						e0[2 * stride + l] = a2i + b2r;
						d0[3 * stride + l] = a2r + b2i;
						e0[3 * stride + l] = a2i - b2r;
						d0[4 * stride + l] = a1r + b1i;
						e0[4 * stride + l] = a1i - b1r;
					}
				}
			}
		}
		double *t;
		t = ar; ar = br; br = t;
		t = ai; ai = bi; bi = t;
		L = Lp; M = Mp;
	}
	*outr = ar;
	*outi = ai;
}

// ------------------------------------------------------------------ per-image state
typedef struct {
	int64_t n;            // transform length
	int rad[40];
	int nrad;
	double *wr, *wi;      // (cos, sin)(2*pi*j/n), j = 0..n-1
	double *phr, *phi;    // Hermitian shift filter for the current slice
} stc_plan;

// Fill the frequency-domain shift filter for fractional shift `s`.  Built symmetrically so the
// filter is EXACTLY Hermitian: that is what makes the two-real-series-per-complex-transform
// packing legal, and it is also the measured Nyquist convention (a real cos(pi*s) factor).
static void stc_filter(stc_plan *p, double s) {
	const int64_t n = p->n, h = n / 2;
	p->phr[0] = 1.0;
	p->phi[0] = 0.0;
	for (int64_t k = 1; k < h; k++) {
		const double th = 2.0 * M_PI * (double)k * s / (double)n;
		const double c = cos(th), sn = sin(th);
		p->phr[k] = c;
		p->phi[k] = -sn;
		p->phr[n - k] = c;
		p->phi[n - k] = sn;
	}
	p->phr[h] = cos(M_PI * s);
	p->phi[h] = 0.0;
}

// ------------------------------------------------------------------ one batch of STC_PAIRS series
// `src`/`dst` address the first sample of the batch's first voxel; consecutive voxels are one
// element apart and consecutive time points are `tstride` elements apart.
static void stc_batch(const float *src, float *dst, int64_t nvox_in_batch, int64_t nt,
                      int64_t tstride, const stc_plan *p, double *buf) {
	const int64_t n = p->n;
	double *ar = buf, *ai = buf + n * STC_LANES;
	double *br = buf + 2 * n * STC_LANES, *bi = buf + 3 * n * STC_LANES;
	double mean[STC_PAIRS], slope[STC_PAIRS], lo[STC_PAIRS], hi[STC_PAIRS];
	double xmin[STC_PAIRS], xmax[STC_PAIRS];
	int bad[STC_PAIRS];
	const double half = 0.5 * (double)(nt - 1);
	const double sdd = (double)nt * ((double)nt * (double)nt - 1.0) / 12.0; // sum (i - half)^2

	/* Clears all four buffers, though strictly only the transform INPUT (ar/ai) needs it: br/bi
	   are pure scratch that the first Stockham stage overwrites completely.  MEASURED and left
	   alone deliberately -- narrowing this to the pad region [nt, n) of ar/ai only moved the
	   254-volume serial kernel from 812 ms to 826 ms, i.e. nothing, because these buffers are
	   L1-resident.  Do not re-optimise it: the narrow form buys no time and makes correctness
	   depend on a non-local property of the FFT's first stage. */
	memset(ar, 0, (size_t)4 * n * STC_LANES * sizeof(double));
	for (int v = 0; v < STC_PAIRS; v++) {
		mean[v] = slope[v] = 0.0;
		lo[v] = hi[v] = xmin[v] = xmax[v] = 0.0;
		bad[v] = (v >= nvox_in_batch);
	}
	// Pass 1: mean, slope, original range; flag any series that is not entirely finite.
	for (int v = 0; v < nvox_in_batch; v++) {
		const float *q = src + v;
		double sum = 0.0, sxy = 0.0, mn = (double)q[0], mx = (double)q[0];
		int ok = 1;
		for (int64_t i = 0; i < nt; i++) {
			const double xv = (double)q[i * tstride];
			if (!(xv >= -DBL_MAX && xv <= DBL_MAX)) { ok = 0; break; }
			sum += xv;
			sxy += ((double)i - half) * xv;
			if (xv < mn) mn = xv;
			if (xv > mx) mx = xv;
		}
		if (!ok) { bad[v] = 1; continue; }
		mean[v] = sum / (double)nt;
		slope[v] = sxy / sdd;
		xmin[v] = mn;
		xmax[v] = mx;
	}
	// Pass 2: detrend into the transform buffer, series 2l real / 2l+1 imaginary of lane l.
	for (int v = 0; v < STC_PAIRS; v++) {
		const int lane = v >> 1;
		double *dstv = (v & 1) ? ai : ar;
		/* A bad lane slot (past the end of the batch, or a non-finite series) still shares a
		   transform with up to seven good ones, so it must stay zero -- which the memset above
		   already guarantees. */
		if (bad[v]) continue;
		const float *q = src + v;
		/* Seed from the first residual, not from 0.  A zero-mean series must straddle zero, so
		   the two agree on real data -- but a constant series detrends to a row of like-signed
		   round-off, and seeding at 0 would widen its clip band by that round-off instead of
		   reproducing the reference exactly. */
		double mn = 0.0, mx = 0.0;
		for (int64_t i = 0; i < nt; i++) {
			const double d = (double)q[i * tstride] - (mean[v] + slope[v] * ((double)i - half));
			dstv[i * STC_LANES + lane] = d;
			if (i == 0) { mn = mx = d; }
			else if (d < mn) mn = d;
			else if (d > mx) mx = d;
		}
		lo[v] = mn;
		hi[v] = mx;
	}
	double *fr, *fi;
	stc_fft(ar, ai, br, bi, n, p->rad, p->nrad, p->wr, p->wi, -1.0, &fr, &fi);
	for (int64_t k = 0; k < n; k++) {
		const double c = p->phr[k], sn = p->phi[k];
		double *xr = fr + k * STC_LANES, *xi = fi + k * STC_LANES;
		for (int l = 0; l < STC_LANES; l++) {
			const double re = xr[l], im = xi[l];
			xr[l] = re * c - im * sn;
			xi[l] = re * sn + im * c;
		}
	}
	// Reuse whichever pair of buffers the forward pass did not end on as the inverse scratch.
	double *sr = (fr == ar) ? br : ar, *si = (fi == ai) ? bi : ai;
	double *gr, *gi;
	stc_fft(fr, fi, sr, si, n, p->rad, p->nrad, p->wr, p->wi, +1.0, &gr, &gi);
	const double inv = 1.0 / (double)n;
	for (int v = 0; v < nvox_in_batch; v++) {
		float *o = dst + v;
		if (bad[v]) {
			// A Fourier shift is global: with a non-finite sample present no output sample of
			// this series is defined.  niimath policy (the manifest); AFNI cannot be an oracle
			// here because its reader zeroes non-finite values before 3dTshift sees them.
			// This applies only to slices that are actually transformed -- a slice below the
			// skip threshold is memcpy'd by the caller and keeps its non-finite samples as they
			// were, which is deliberate: the verbatim copy is the measured behaviour.
			for (int64_t i = 0; i < nt; i++) o[i * tstride] = (float)NAN;
			continue;
		}
		const int lane = v >> 1;
		const double *g = ((v & 1) ? gi : gr) + lane;
		for (int64_t i = 0; i < nt; i++) {
			double y = g[i * STC_LANES] * inv;
			if (y < lo[v]) y = lo[v];
			else if (y > hi[v]) y = hi[v];
			y += mean[v] + slope[v] * ((double)i - half);
			if (y < xmin[v]) y = xmin[v];
			else if (y > xmax[v]) y = xmax[v];
			o[i * tstride] = (float)y;
		}
	}
}

// ------------------------------------------------------------------ header time units
// Returns the seconds-per-header-time-unit scale, or 0 for a missing/unsupported unit.
static double stc_time_scale(int units) {
	switch (units) {
		case NIFTI_UNITS_SEC: return 1.0;
		case NIFTI_UNITS_MSEC: return 1.0e-3;
		case NIFTI_UNITS_USEC: return 1.0e-6;
		default: return 0.0;
	}
}

// ------------------------------------------------------------------ entry point
int nii_stc(nifti_image *nim, const double *times, int ntimes, int have_tzero, double tzero) {
	if (!nim || nim->datatype != DT_FLOAT32 || !nim->data || !times) {
		printfx("-stc: internal error (expected float32 image)\n");
		return 1;
	}
	if (nim->ndim > 4) {
		printfx("-stc requires a scalar 4D image (got %lldD)\n", (long long)nim->ndim);
		return 1;
	}
	const int64_t nt = (nim->ndim > 3) ? nim->nt : 1;
	if (nim->ndim < 4 || nt < STC_MIN_NT) {
		printfx("-stc requires a 4D image with at least %d volumes (got nt = %lld)\n",
		        STC_MIN_NT, (long long)nt);
		return 1;
	}
	const int64_t nx = nim->nx, ny = nim->ny, nz = nim->nz;
	if (nx < 1 || ny < 1 || nz < 1) {
		printfx("-stc: degenerate spatial dimensions\n");
		return 1;
	}
	if (ntimes != (int)nz) {
		printfx("-stc: --slicetiming has %d value%s but the image has %lld slices\n",
		        ntimes, ntimes == 1 ? "" : "s", (long long)nz);
		return 1;
	}
	const double tscale = stc_time_scale(nim->time_units);
	if (tscale == 0.0) {
		printfx("-stc: the header has no usable temporal unit (xyzt_units time field = %d); "
		        "-stc will not assume seconds\n", nim->time_units);
		return 1;
	}
	const double tr = nim->dt * tscale;
	if (!(tr > 0.0 && tr <= DBL_MAX)) {
		printfx("-stc: the header TR is not a finite positive time (pixdim[4] = %g)\n", nim->dt);
		return 1;
	}
	double tmin = times[0], tmax = times[0], tsum = 0.0;
	for (int z = 0; z < ntimes; z++) {
		const double v = times[z];
		if (!(v >= -DBL_MAX && v <= DBL_MAX)) {
			printfx("-stc: slice time %d is not finite\n", z);
			return 1;
		}
		if (v < 0.0 || v > tr) {
			printfx("-stc: slice time %d is %g s, outside [0, TR] with TR = %g s\n", z, v, tr);
			return 1;
		}
		if (v < tmin) tmin = v;
		if (v > tmax) tmax = v;
		tsum += v;
	}
	if (!have_tzero) tzero = tsum / (double)ntimes;
	if (!(tzero >= -DBL_MAX && tzero <= DBL_MAX)) {
		printfx("-stc: -tzero is not finite\n");
		return 1;
	}
	if (tzero < tmin || tzero > tmax) {
		printfx("-stc: -tzero %g s is outside the slice-time range [%g, %g] s\n", tzero, tmin, tmax);
		return 1;
	}

	const int64_t n = stc_fft_len(nt);
	stc_plan plan;
	memset(&plan, 0, sizeof(plan));
	plan.n = n;
	plan.nrad = n ? stc_factor(n, plan.rad, (int)(sizeof(plan.rad) / sizeof(plan.rad[0]))) : 0;
	if (!n || !plan.nrad) {
		printfx("-stc: no supported transform length for nt = %lld\n", (long long)nt);
		return 1;
	}

	const int64_t nxy = nx * ny;
	/* Every product goes through the checked multiply, including the intermediates: writing
	   nii_mul_size(a * b, c * d, ...) would leave a*b and c*d themselves unguarded while reading
	   as if the whole expression were covered. */
	size_t nb_out, nb_tw, nb_buf, nvox_all, nslice;
	if (nii_mul_size((size_t)nxy, (size_t)nz, &nslice) ||
	    nii_mul_size(nslice, (size_t)nt, &nvox_all) ||
	    nii_mul_size(nvox_all, sizeof(float), &nb_out) ||
	    nii_mul_size((size_t)n, sizeof(double), &nb_tw) ||
	    nii_mul_size((size_t)n, (size_t)(4 * STC_LANES) * sizeof(double), &nb_buf)) {
		printfx("-stc: image too large for this build\n");
		return 1;
	}
	int nthread = 1;
#ifdef _OPENMP
	nthread = omp_get_max_threads();
	if (nthread < 1) nthread = 1;
#endif
	/* Never reserve more scratch slices than there is work. A thin slice (small nxy) yields a
	   handful of batches, and reserving one 4*n*STC_LANES buffer per available core would then be
	   almost entirely untouched: at n = 4194304 that is half a gigabyte per core for, possibly,
	   a single batch. The num_threads clause below pins the team to this same bound, which is
	   what makes indexing the scratch by omp_get_thread_num() safe. */
	const int64_t maxbatch = (nxy + STC_PAIRS - 1) / STC_PAIRS;
	if ((int64_t)nthread > maxbatch) nthread = (int)maxbatch;
	if (nthread < 1) nthread = 1;
	size_t nb_scratch;
	if (nii_mul_size(nb_buf, (size_t)nthread, &nb_scratch)) {
		printfx("-stc: image too large for this build\n");
		return 1;
	}

	int rc = 0;
	float *out = (float *)malloc(nb_out);
	double *scratch = (double *)malloc(nb_scratch);
	plan.wr = (double *)malloc(nb_tw);
	plan.wi = (double *)malloc(nb_tw);
	plan.phr = (double *)malloc(nb_tw);
	plan.phi = (double *)malloc(nb_tw);
	if (!out || !scratch || !plan.wr || !plan.wi || !plan.phr || !plan.phi) {
		printfx("-stc: out of memory\n");
		rc = 1;
		goto done;
	}
	for (int64_t j = 0; j < n; j++) {
		const double th = 2.0 * M_PI * (double)j / (double)n;
		plan.wr[j] = cos(th);
		plan.wi[j] = sin(th);
	}

	{
		const float *img = (const float *)nim->data;
		const int64_t nvol = nxy * nz;
		for (int64_t z = 0; z < nz; z++) {
			const double s = (times[z] - tzero) / tr;
			const int64_t off = z * nxy;
			if (fabs(s) < STC_SKIP) {
				// Measured: AFNI copies such a slice verbatim.  Copy per volume so the
				// output is bit-identical to the input for this slice.
				for (int64_t t = 0; t < nt; t++)
					memcpy(out + t * nvol + off, img + t * nvol + off,
					       (size_t)nxy * sizeof(float));
				continue;
			}
			stc_filter(&plan, s);
			const int64_t nbatch = (nxy + STC_PAIRS - 1) / STC_PAIRS;
#ifdef _OPENMP
			/* num_threads caps the team at exactly the number of scratch slices allocated above,
			   so omp_get_thread_num() is always a valid index. A team may be SMALLER than this
			   (the runtime is allowed to shrink it) but never larger. */
			#pragma omp parallel for schedule(static) num_threads(nthread)
#endif
			for (int64_t b = 0; b < nbatch; b++) {
				int tid = 0;
#ifdef _OPENMP
				tid = omp_get_thread_num();
#endif
				const int64_t v0 = b * STC_PAIRS;
				int64_t cnt = nxy - v0;
				if (cnt > STC_PAIRS) cnt = STC_PAIRS;
				stc_batch(img + off + v0, out + off + v0, cnt, nt, nvol, &plan,
				          scratch + (size_t)tid * (size_t)n * 4 * STC_LANES);
			}
		}
		free(nim->data);
		nim->data = out;
		out = NULL;
		// Measured: the oracle moves the time origin to the common time point.  toffset lives
		// in the header's own time units, so convert back out of seconds.
		nim->toffset = tzero / tscale;
	}

done:
	free(out);
	free(scratch);
	free(plan.wr);
	free(plan.wi);
	free(plan.phr);
	free(plan.phi);
	return rc;
}

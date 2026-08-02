// fmap.c - B0 fieldmap EPI distortion correction (-fugue)
//
// Clean-room; see fmap.h for the provenance statement and the licence boundary, and the
// fmap_bench repository's test/fmap_reference_manifest.md (section M1) for the experiment behind
// every value marked "measured" below.  FSL's sources were never read; its executables were used
// only as a black-box oracle.
//
// The measured apply stage, in full:
//
//   s(v)   = fmap(v) / (2*pi) * dwell * N_axis                     [voxels]
//   out(v) = linear_interp_1d(in, v + sigma * s(v) * e_axis)       sigma = +1, or -1 for "a-"
//
// Four things here are easy to get wrong and each was established by a separate experiment:
//
//   1. The interpolation is 1D LINEAR along the unwarp axis, not cubic, spline or Lanczos.  An
//      impulse displaced by half a voxel comes back as exactly two taps of 0.5; by a quarter, as
//      0.75/0.25.  No wider kernel can do that.  (This is why -fugue does not reuse -unwarp's
//      md_pull, which is a 3D Lanczos-5 pull -- see the note at the bottom of this file.)
//   2. The shift is sampled at the OUTPUT voxel.  This is a pure pull with no inversion, so mass
//      is NOT conserved under compression: an impulse through a ramp field comes back summing to
//      893.33, not 1000, and the reference agrees with the pull model on that deficit to six
//      significant figures.  Anything that conserves mass here is a different algorithm.
//   3. N_axis is the image dimension ALONG THE UNWARP AXIS, not a fixed axis.
//   4. A zero fieldmap does not mean zero shift; it means NO DATA.  The shift field is
//      extrapolated over unsupported voxels before use.  Omitting this is not an edge-case
//      refinement -- it changes a third of the voxels in a real image and drops agreement with
//      the reference from r = 0.99967 to r = 0.9940.

#define _USE_MATH_DEFINES // microsoft compiler
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif
/* No <omp.h>: this file uses only #pragma omp, which needs no header. */

#include "fmap.h"
#include "print.h"
#ifdef HAVE_ROMEO
#include "romeo.h" // in-memory phase unwrapping for -fmapprep
#endif

#define FM_ERR(...) do { printfx("-fugue: " __VA_ARGS__); } while (0)
#define FP_ERR(...) do { printfx("-fmapprep: " __VA_ARGS__); } while (0)

// A voxel whose |fieldmap| falls below this carries no data.  Measured bound on the reference's
// own threshold: (3.278e-07, 0.6237] rad/s -- six orders of magnitude of daylight between the
// largest value it discards and the smallest it keeps.  At 1e-6 rad/s the implied shift is
// 6.5e-9 voxels, so nothing physical is discarded by sitting in the middle of that gap.
#define FMAP_EPS 1e-6

/* ============================== argument parsing ============================== */

// Accepts the reference's x/y/z and the BIDS spellings i/j/k, with an optional trailing '-'.
// Returns 0 on success and writes the axis index and the sign.
static int fm_parse_dir(const char *dir, int *axis, double *sigma) {
	size_t len;
	char a;
	if (!dir) return 1;
	len = strlen(dir);
	if (len < 1 || len > 2) return 1;
	if (len == 2) {
		if (dir[1] != '-' && dir[1] != '+') return 1;
		*sigma = (dir[1] == '-') ? -1.0 : 1.0;
	} else
		*sigma = 1.0;
	a = dir[0];
	if (a == 'x' || a == 'i') *axis = 0;
	else if (a == 'y' || a == 'j') *axis = 1;
	else if (a == 'z' || a == 'k') *axis = 2;
	else return 1;
	return 0;
}

/* ============================== fieldmap input ============================== */

// Read a NIfTI as float32.  Header-only preflight first, so a malformed or oversized image is
// rejected before its payload is decompressed and allocated.
/* `op` and `what` are threaded through because this reader is shared: -fugue passes
   ("-fugue", "fieldmap") and -fmapprep passes ("-fmapprep", "magnitude").  Routing its
   messages through FM_ERR instead told a -fmapprep user that "-fugue" could not read their
   "fieldmap", naming an operation absent from their command line and the wrong file role. */
static nifti_image *fm_read_f32(const char *fn, const char *op, const char *what) {
	nifti_image *n;
	in_hdr ihdr;
	{
		nifti_image *h = nifti_image_read(fn, 0);
		int bad = 0;
		if (!h) { printfx("%s: failed to read the header of %s '%s'\n", op, what, fn); return NULL; }
		if (h->nvox < 1 || h->nx < 1 || h->ny < 1 || h->nz < 1) {
			printfx("%s: %s '%s' has invalid dimensions\n", op, what, fn); bad = 1;
		} else if (h->nt > 1 || h->nu > 1 || h->nv > 1 || h->nw > 1) {
			printfx("%s: %s '%s' must be a single 3D volume\n", op, what, fn); bad = 1;
		} else if ((int64_t)h->nvox > INT_MAX) {
			printfx("%s: %s '%s' exceeds INT_MAX voxels; %s is not a huge-image-safe operation\n",
				op, what, fn, op);
			bad = 1;
		}
		nifti_image_free(h);
		if (bad) return NULL;
	}
	n = nifti_image_read(fn, 1);
	if (!n) { printfx("%s: failed to read %s '%s'\n", op, what, fn); return NULL; }
	/* Re-check after the load as well as before it: the preflight is what stops us decompressing
	   a huge payload, this is the fail-closed guarantee, and it costs nothing. */
	if (n->nvox < 1 || n->nx < 1 || n->ny < 1 || n->nz < 1 ||
		n->nt > 1 || n->nu > 1 || n->nv > 1 || n->nw > 1 || (int64_t)n->nvox > INT_MAX) {
		printfx("%s: %s '%s' changed on disk or has unusable dimensions\n", op, what, fn);
		nifti_image_free(n); return NULL;
	}
	ihdr = set_input_hdr(n);
	/* Convert when the stored type is not float32, but ALSO when it IS float32 and carries a
	   non-trivial scl_slope/scl_inter -- otherwise a scaled float32 fieldmap is used raw, which
	   silently rescales every shift in the image. */
	if (n->datatype != DT_FLOAT32 ||
		(n->scl_slope != 0.0f && n->scl_slope != 1.0f) || n->scl_inter != 0.0f) {
		if (nifti_image_change_datatype(n, DT_FLOAT32, &ihdr) != 0) {
			printfx("%s: failed to convert %s '%s' to float32\n", op, what, fn);
			nifti_image_free(n); return NULL;
		}
	}
	return n;
}

/* ============================== shift field ============================== */

// Fill one line of the shift field, in place, over the voxels the fieldmap does not support.
//
// Measured: beyond the ends of the supported run the reference REPLICATES the nearest supported
// value, and a line with no supported voxel at all stays at exactly 0 -- which is what shows the
// fill to be per-line and one-dimensional rather than a 3D diffusion.
//
// Interior gaps are a DELIBERATE DIVERGENCE.  The reference produces a value strictly between a
// nearest-valid fill and linear interpolation, and no single rule tested reproduces both probe
// geometries (a full-slab gap needs ~1 diffusion sweep, a single-column gap ~3-4).  Rather than
// approximate an unidentified scheme, this interpolates linearly: the smoothest interpolant
// consistent with the supported data along the distortion axis, defensible on its own terms.
static void fm_fill_line(float *s, const unsigned char *ok, int64_t n, int64_t step) {
	int64_t i, k, prev = -1;
	for (i = 0; i < n; i++) {
		if (!ok[i]) continue;
		if (prev < 0) {
			for (k = 0; k < i; k++) s[k * step] = s[i * step];      // leading run: replicate
		} else if (i > prev + 1) {
			const double a = s[prev * step], b = s[i * step];       // interior gap: linear
			const double d = (double)(i - prev);
			for (k = prev + 1; k < i; k++)
				s[k * step] = (float)(a + (b - a) * ((double)(k - prev) / d));
		}
		prev = i;
	}
	if (prev < 0) {
		for (k = 0; k < n; k++) s[k * step] = 0.0f;                 // no data anywhere on this line
	} else {
		for (k = prev + 1; k < n; k++) s[k * step] = s[prev * step]; // trailing run: replicate
	}
}

/* ============================== the operation ============================== */

int fmap_unwarp(nifti_image *nim, const char *fmapfile, double dwell, const char *unwarpdir) {
	nifti_image *fm = NULL;
	float *shift = NULL;
	const float *fdat;
	double sigma = 1.0, scale;
	int axis = 1, rc = 1;
	int64_t nx, ny, nz, n3, nt, i;
	int64_t dim[3], str[3], step, nline, nax, bstep, cstep, nb, nc;
	int oom = 0;

	if (!nim || nim->datatype != DT_FLOAT32) {
		FM_ERR("internal error: expected a float32 working image\n"); return 1;
	}
	if (nim->nu > 1 || nim->nv > 1 || nim->nw > 1) { FM_ERR("input must be 3D or 4D\n"); return 1; }
	if ((int64_t)nim->nvox > INT_MAX) {
		FM_ERR("input exceeds INT_MAX voxels; -fugue is not a huge-image-safe operation\n"); return 1;
	}
	if (fm_parse_dir(unwarpdir, &axis, &sigma)) {
		FM_ERR("unwarpdir must be one of x y z (or i j k), with an optional trailing '-'; got '%s'\n",
			unwarpdir ? unwarpdir : "(null)");
		return 1;
	}
	/* Guard the DOUBLE before it reaches any conversion or division.  A NaN dwell would otherwise
	   produce a NaN shift field and a silently all-NaN output. */
	if (!(dwell > 0.0) || !(dwell <= DBL_MAX)) {
		FM_ERR("dwell must be a positive, finite echo spacing in seconds; got %g\n", dwell); return 1;
	}

	nx = nim->nx; ny = nim->ny; nz = (nim->nz < 1 ? 1 : nim->nz);
	n3 = nx * ny * nz;
	if (n3 < 1 || (int64_t)nim->nvox % n3 != 0) { FM_ERR("invalid image geometry\n"); return 1; }
	nt = (int64_t)nim->nvox / n3;

	fm = fm_read_f32(fmapfile, "-fugue", "fieldmap");
	if (!fm) return 1;
	if (fm->nx != nx || fm->ny != ny || fm->nz != nz || max_displacement_mm(nim, fm) > 0.001f) {
		FM_ERR("fieldmap '%s' does not share the input's grid (dimensions and world transform must match)\n",
			fmapfile);
		goto done;
	}
	fdat = (const float *)fm->data;
	for (i = 0; i < n3; i++) {
		/* Magnitude guard rather than isfinite(): the whole program is built -ffast-math, under
		   which isfinite() is not reliable.  This is the project's standing idiom. */
		if (!(fdat[i] >= -FLT_MAX && fdat[i] <= FLT_MAX)) {
			FM_ERR("fieldmap '%s' contains a non-finite value; refusing to write a corrupted image\n",
				fmapfile);
			goto done;
		}
	}

	dim[0] = nx; dim[1] = ny; dim[2] = nz;
	str[0] = 1; str[1] = nx; str[2] = nx * ny;
	step = str[axis];
	nax = dim[axis];
	/* The two axes that are not the unwarp axis, in increasing order, enumerate the lines. */
	{
		int b = (axis == 0) ? 1 : 0;
		int c = (axis == 2) ? 1 : 2;
		bstep = str[b]; cstep = str[c]; nb = dim[b]; nc = dim[c];
	}
	nline = nb * nc;

	shift = (float *)malloc((size_t)n3 * sizeof(float));
	if (!shift) { FM_ERR("out of memory allocating the shift field\n"); goto done; }

	/* Voxels, not millimetres, and N is the dimension along the unwarp axis (measured). */
	scale = dwell * (double)nax / (2.0 * M_PI);
	for (i = 0; i < n3; i++) shift[i] = (float)((double)fdat[i] * scale);

#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for (int64_t L = 0; L < nline; L++) {
		unsigned char ok[512];
		unsigned char *heap = NULL;
		unsigned char *sup = ok;
		const int64_t base = (L % nb) * bstep + (L / nb) * cstep;
		int64_t k;
		if (nax > (int64_t)sizeof(ok)) {
			heap = (unsigned char *)malloc((size_t)nax);
			if (!heap) {
#ifdef _OPENMP
				#pragma omp atomic write
#endif
				oom = 1;
				continue;
			}
			sup = heap;
		}
		for (k = 0; k < nax; k++)
			sup[k] = (unsigned char)(fabs((double)fdat[base + k * step]) > FMAP_EPS);
		fm_fill_line(shift + base, sup, nax, step);
		free(heap);
	}
	if (oom) { FM_ERR("out of memory building the shift field\n"); goto done; }

	/* Resample IN PLACE, one line at a time.  Each line is gathered into a thread-local buffer,
	   rewritten, and stored back, so the whole operation needs no second copy of the image -- only
	   the shift volume above.  Every line is independent and nothing is reduced across lines, so
	   the result is byte-identical whatever the thread count. */
#ifdef _OPENMP
	#pragma omp parallel
#endif
	{
		double stackbuf[512];
		double *heap = NULL;
		double *buf = stackbuf;
		if (nax > (int64_t)(sizeof(stackbuf) / sizeof(stackbuf[0]))) {
			heap = (double *)malloc((size_t)nax * sizeof(double));
			if (!heap) {
#ifdef _OPENMP
				#pragma omp atomic write
#endif
				oom = 1;
			} else
				buf = heap;
		}
		/* EVERY thread must reach the same decision about the worksharing loop below.  Branching
		   on a THREAD-LOCAL allocation result would let one thread skip an `omp for` that the
		   others enter -- the team then encounters different worksharing regions, which is
		   undefined behaviour and in practice hangs at the implicit barrier.  The barrier
		   publishes `oom` (it implies a flush) so the test that follows is team-wide. */
#ifdef _OPENMP
		#pragma omp barrier
#endif
		if (!oom) {
			float *img = (float *)nim->data;
			int64_t tt;
			for (tt = 0; tt < nt; tt++) {
				float *vol = img + tt * n3;
#ifdef _OPENMP
				#pragma omp for schedule(static)
#endif
				for (int64_t L = 0; L < nline; L++) {
					const int64_t base = (L % nb) * bstep + (L / nb) * cstep;
					const float *sl = shift + base;
					int64_t k;
					for (k = 0; k < nax; k++) buf[k] = (double)vol[base + k * step];
					for (k = 0; k < nax; k++) {
						const double src = (double)k + sigma * (double)sl[k * step];
						double fl, w, v;
						int64_t j;
						/* Bound the coordinate BEFORE the cast.  A shift far outside the volume
						   would otherwise make (int64_t)floor(src) undefined; the guard also
						   supplies the measured out-of-FOV behaviour, which is a 0 fill. */
						if (!(src > -1.0 && src < (double)nax)) { vol[base + k * step] = 0.0f; continue; }
						fl = floor(src);
						j = (int64_t)fl;
						w = src - fl;
						v = 0.0;
						if (j >= 0) v += buf[j] * (1.0 - w);
						if (j + 1 < nax) v += buf[j + 1] * w;
						vol[base + k * step] = (float)v;
					}
				}
			}
		}
		free(heap);
	}
	if (oom) {
		/* The image is now partially rewritten, so it cannot be handed back.  The op-loop caller
		   frees `nim` without saving when a non-zero status is returned, so reporting the failure
		   is sufficient; there is no snapshot to restore. */
		FM_ERR("out of memory resampling; the working image is no longer valid\n");
		goto done;
	}
	rc = 0;
done:
	free(shift);
	nifti_image_free(fm);
	return rc;
}

/* ============================== -fmapprep ============================== */

#ifdef HAVE_ROMEO

// Comparator-free selection of the k-th smallest, in place.  Deliberately not qsort: emscripten's
// qsort dispatches the comparator through call_indirect once per comparison, which is ~100x slower
// than native and is a standing rule in this tree for anything touching a per-voxel array.
static float fm_select_kth(float *a, int64_t n, int64_t k) {
	int64_t lo = 0, hi = n - 1;
	while (lo < hi) {
		/* Median-of-three pivot: keeps the already-sorted and reverse-sorted cases -- both of
		   which a masked, spatially coherent fieldmap can approximate -- off the O(n^2) path. */
		int64_t mid = lo + (hi - lo) / 2, i = lo, j = hi;
		float p;
		if (a[mid] < a[lo]) { float t = a[mid]; a[mid] = a[lo]; a[lo] = t; }
		if (a[hi] < a[lo])  { float t = a[hi];  a[hi]  = a[lo]; a[lo] = t; }
		if (a[hi] < a[mid]) { float t = a[hi];  a[hi]  = a[mid]; a[mid] = t; }
		p = a[mid];
		while (i <= j) {
			while (a[i] < p) i++;
			while (a[j] > p) j--;
			if (i <= j) { float t = a[i]; a[i] = a[j]; a[j] = t; i++; j--; }
		}
		if (k <= j) hi = j;
		else if (k >= i) lo = i;
		else return a[k];
	}
	return a[lo];
}

// Maximum branch-correction passes.  Measured on the benchmark: 4 passes reach a fixed point
// (30, 6, 5, 3 voxels moved, then 0).  The cap only bounds a pathological input; convergence is
// the normal exit.
#define FMAP_DEBRANCH_PASSES 10

// Move voxels that ROMEO left on the wrong 2*pi branch onto the branch their neighbours agree on.
//
// This is a DELIBERATE DIVERGENCE from the reference, which applies no post-processing at all
// (measured -- a single-voxel spike passes through fsl_prepare_fieldmap intact).  It exists
// because niimath unwraps with ROMEO rather than PRELUDE, and on this data ROMEO leaves ~30
// isolated voxels a full 2*pi out: measured by their distance from their own 6-neighbour median,
// ROMEO produced 30 such outliers where PRELUDE produced 2.  Each one becomes a ~17-voxel shift
// error in the corrected EPI, so they are visible artefacts, not rounding.
//
// It cannot corrupt a genuine field value: the correction is always an INTEGER multiple of
// 2*pi/deltaTE, i.e. a projection onto the set of unwrappings consistent with the same wrapped
// phase.  A voxel only moves when it sits more than half a wrap from the median of its in-mask
// 6-neighbours, and a true field gradient that steep is beyond what a phase DIFFERENCE can
// represent anyway -- it would already be aliased.
//
// Jacobi, not Gauss-Seidel: every pass reads a snapshot, so the result does not depend on
// traversal order and is byte-identical across thread counts.
static int fm_debranch(float *f, const unsigned char *mask, int64_t nx, int64_t ny, int64_t nz,
	double wrap) {
	float *snap = NULL;
	int64_t n3 = nx * ny * nz, x, y, z, pass;
	if (!(wrap > 0.0) || !(wrap <= DBL_MAX)) return 0;   /* nothing sane to project onto */
	snap = (float *)malloc((size_t)n3 * sizeof(float));
	if (!snap) return 1;
	for (pass = 0; pass < FMAP_DEBRANCH_PASSES; pass++) {
		int64_t moved = 0;
		memcpy(snap, f, (size_t)n3 * sizeof(float));
		for (z = 0; z < nz; z++) for (y = 0; y < ny; y++) for (x = 0; x < nx; x++) {
			const int64_t o = x + y * nx + z * nx * ny;
			double v[6], nm, k;
			int n = 0, a, b;
			if (!mask[o]) continue;
			if (x > 0      && mask[o - 1])          v[n++] = snap[o - 1];
			if (x < nx - 1 && mask[o + 1])          v[n++] = snap[o + 1];
			if (y > 0      && mask[o - nx])         v[n++] = snap[o - nx];
			if (y < ny - 1 && mask[o + nx])         v[n++] = snap[o + nx];
			if (z > 0      && mask[o - nx * ny])    v[n++] = snap[o - nx * ny];
			if (z < nz - 1 && mask[o + nx * ny])    v[n++] = snap[o + nx * ny];
			if (n < 3) continue;   /* too few neighbours to out-vote the voxel itself */
			/* Insertion sort of at most 6 values: small fixed N, so no comparator indirection
			   (emscripten's qsort would dispatch through call_indirect per comparison). */
			for (a = 1; a < n; a++) {
				double t = v[a];
				for (b = a - 1; b >= 0 && v[b] > t; b--) v[b + 1] = v[b];
				v[b + 1] = t;
			}
			nm = v[n / 2];   /* upper central value, matching the demedian convention */
			k = floor(((double)snap[o] - nm) / wrap + 0.5);
			if (k != 0.0) { f[o] = (float)((double)snap[o] - k * wrap); moved++; }
		}
		if (moved == 0) break;
	}
	free(snap);
	return 0;
}

int fmap_prepare(nifti_image *nim, const char *magfile, double delta_te_ms, int debranch) {
	nifti_image *mg = NULL;
	unsigned char *mask = NULL;
	float *sel = NULL;
	float *ph;
	const float *mag;
	double lo = 0.0, hi = 0.0, scale, offset, inv_te, median;
	int64_t nx, ny, nz, n3, i, nmask = 0;
	int rc = 1;
	romeo_opts o;

	if (!nim || nim->datatype != DT_FLOAT32) {
		printfx("-fmapprep: internal error: expected a float32 working image\n"); return 1;
	}
	if (nim->nt > 1 || nim->nu > 1 || nim->nv > 1 || nim->nw > 1) {
		printfx("-fmapprep: the phase-difference input must be a single 3D volume\n"); return 1;
	}
	if ((int64_t)nim->nvox > INT_MAX) {
		printfx("-fmapprep: input exceeds INT_MAX voxels; -fmapprep is not a huge-image-safe operation\n"); return 1;
	}
	if (!(delta_te_ms > 0.0) || !(delta_te_ms <= DBL_MAX)) {
		printfx("-fmapprep: deltaTE must be a positive, finite echo time difference in milliseconds; got %g\n",
			delta_te_ms);
		return 1;
	}
	nx = nim->nx; ny = nim->ny; nz = (nim->nz < 1 ? 1 : nim->nz);
	n3 = nx * ny * nz;
	if (n3 < 1 || (int64_t)nim->nvox != n3) { printfx("-fmapprep: invalid image geometry\n"); return 1; }

	mg = fm_read_f32(magfile, "-fmapprep", "magnitude");
	if (!mg) return 1;
	if (mg->nx != nx || mg->ny != ny || mg->nz != nz || max_displacement_mm(nim, mg) > 0.001f) {
		printfx("-fmapprep: magnitude '%s' does not share the phase image's grid (dimensions and world transform must match)\n",
			magfile);
		goto done;
	}
	mag = (const float *)mg->data;
	ph = (float *)nim->data;

	mask = (unsigned char *)malloc((size_t)n3);
	if (!mask) { printfx("-fmapprep: out of memory allocating the mask\n"); goto done; }
	for (i = 0; i < n3; i++) {
		/* Measured: the mask is the magnitude's nonzero support, used verbatim.  Magnitude-guard
		   rather than isfinite(): the program is built -ffast-math. */
		mask[i] = (unsigned char)(mag[i] >= -FLT_MAX && mag[i] <= FLT_MAX && mag[i] != 0.0f);
		nmask += mask[i];
	}
	if (nmask < 1) {
		printfx("-fmapprep: magnitude '%s' has no nonzero voxels, so there is no brain to build a fieldmap over\n",
			magfile);
		goto done;
	}

	/* Rescale the stored phase so its observed range spans exactly one 2*pi period, which is what
	   a wrapped phase difference is.  The reference instead demands its input pre-scaled onto a
	   fixed 0..4096 and multiplies by 2*pi/4096; doing it from the observed range is equivalent on
	   such an image and additionally accepts the other encodings a converter may emit (0..4095,
	   -4096..4094, and radians already).  The absolute offset is irrelevant -- the demedian at the
	   end removes any constant -- so only the SPAN has to be right. */
	for (i = 0; i < n3; i++) {
		if (!(ph[i] >= -FLT_MAX && ph[i] <= FLT_MAX)) {
			printfx("-fmapprep: the phase image contains a non-finite value\n");
			goto done;
		}
		if (i == 0 || ph[i] < lo) lo = ph[i];
		if (i == 0 || ph[i] > hi) hi = ph[i];
	}
	if (!(hi > lo)) {
		printfx("-fmapprep: the phase image is constant (%g everywhere); it carries no field information\n", lo);
		goto done;
	}
	scale = (2.0 * M_PI) / (hi - lo);
	offset = -M_PI - lo * scale;
	for (i = 0; i < n3; i++) ph[i] = (float)((double)ph[i] * scale + offset);

	o = romeo_opts_default();
	o.nTE = 1;
	o.TEs[0] = delta_te_ms;
	o.te_epi = 0;
	o.template_echo = 1;
	o.individual = 0;
	o.correctglobal = 0;
	o.no_phase_rescale = 1;   /* already in radians, above */
	/* mask_in is supplied, so ROMEO's own robustmask is bypassed entirely and the mask stays
	   byte-identical to the caller's brain extraction -- the point of taking a magnitude at all. */
	if (romeo_unwrap_frame(ph, mag, 1, (int)nx, (int)ny, (int)nz, 1, o.TEs, &o, mask, NULL)) {
		printfx("-fmapprep: ROMEO failed to unwrap the phase difference\n");
		goto done;
	}

	inv_te = 1000.0 / delta_te_ms;   /* deltaTE arrives in milliseconds; the field is rad/s */
	for (i = 0; i < n3; i++) ph[i] = (float)((double)ph[i] * inv_te);

	/* Applied HERE, in fmap.c, strictly AFTER romeo_unwrap_frame() has returned -- romeo.c is not
	   touched and neither --medic nor -romeo can see this.  `-no-debranch` turns it off so the raw
	   ROMEO field can be obtained for parity tracing against a reference MEDIC implementation. */
	if (debranch && fm_debranch(ph, mask, nx, ny, nz, 2.0 * M_PI * inv_te)) {
		printfx("-fmapprep: out of memory correcting 2*pi branch outliers\n");
		goto done;
	}

	/* Measured: subtract the median over the mask, taking the UPPER of the two central values for
	   an even population.  Averaging the two central values is wrong, and not subtly so -- on a
	   field whose mask splits into two equal populations it errs by half the field's range. */
	sel = (float *)malloc((size_t)nmask * sizeof(float));
	if (!sel) { printfx("-fmapprep: out of memory computing the median\n"); goto done; }
	{
		int64_t j = 0;
		for (i = 0; i < n3; i++) if (mask[i]) sel[j++] = ph[i];
		median = (double)fm_select_kth(sel, nmask, nmask / 2);
	}
	for (i = 0; i < n3; i++) ph[i] = mask[i] ? (float)((double)ph[i] - median) : 0.0f;

	rc = 0;
done:
	free(sel);
	free(mask);
	nifti_image_free(mg);
	return rc;
}

#endif // HAVE_ROMEO

/* Why this does not reuse medic_unwarp's md_pull, despite the plan's preference for sharing:
   md_pull is a 3D Lanczos-5 pull driven by a displacement map in MILLIMETRES with an arbitrary
   world-space direction.  -fugue is a 1D LINEAR pull driven by a shift in VOXELS along a single
   storage axis, with a support-aware extrapolation step md_pull has no notion of.  The two share
   a loop shape and nothing else; unifying them would mean threading a mode flag that changes the
   kernel, its dimensionality and its units, which is two functions wearing one name.  Recorded as
   a deviation in the manifest rather than done quietly. */

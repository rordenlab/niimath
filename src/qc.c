// niimath --qc : MRIQC-style anatomical quality metrics from a T1 image + segmentation.
//
// Given a T1-weighted intensity image and a matching integer tissue segmentation
// with the label values that denote CSF and white matter (or, with --pve, soft
// partial-volume maps; below), this computes a subset of MRIQC's anatomical Image
// Quality Metrics (IQMs); --air adds the background (air) metrics. The AGENTS.md
// qc.c entry records the computable-vs-blocked categorisation and the gotchas.
//
// Metric formulas follow MRIQC (mriqc/qc/anatomical.py, mriqc/interfaces/
// anatomical.py; Esteban et al., PLOS ONE 2017). This hard-segmentation variant
// deliberately uses unrounded intensities and NumPy-linear quantiles; MRIQC's
// soft-PVM summary path rounds intensities and uses weighted quantiles.
//   CJV     = (mad_WM + mad_GM) / |median_WM - median_GM|            (median, MAD)
//   CNR     = |median_WM - median_GM| / sqrt(s_bg^2 + s_WM^2 + s_GM^2)  (median, stdv)
//   SNR_t   = median_t / (stdv_t * sqrt(n/(n-1)))                    (median, stdv)
//   WM2MAX  = median_WM / P99.95(image)
//   EFC     = entropy-focus criterion over non-zero voxels
// We lack the air term s_bg, so CNR is emitted as `cnr_noair` (s_bg = 0): valid
// as a relative contrast measure but NOT comparable to MRIQC's normative values.
//
// MRIQC computes per-tissue stats as partial-volume-WEIGHTED statistics over soft
// pvms; with a hard label map the weights are 0/1, i.e. ordinary statistics over
// the masked voxels (the binary special case). --erode 1 erodes tissue masks one
// voxel before the intensity statistics; off by default, as on ABIDE it moved WM/GM
// stdv and MAD away from MRIQC (0.71x -> 1.0x without it). ICV fractions and absolute
// volumes always use the FULL (un-eroded) tissue extent.
//
// --pve <csf> <gm> <wm> takes partial-volume fraction maps instead and weights each
// voxel by its fraction, as MRIQC's summary_stats does: mean, stdv and the
// quantiles are weighted, n = sum of weights, and MAD and kurtosis use the voxels
// above half the tissue's peak fraction. The weighted quantile keeps the hard path's
// NumPy-linear definition (element k sits at rank sum(w before k), the target at
// p*(W-1)), so 0/1 fractions reproduce --seg --erode 0 exactly. No erosion: the
// fractions already model the boundary.
//
// Intensities follow MRIQC unless --raw: IntensityClip (10th..99.9th percentiles of the
// median-filtered positive voxels, shifted and scaled to 0..255, rounded), N4 weighted by the
// brain (HAVE_N4; ANTs defaults), then Harmonize (x 1000 / median of the 18-connected-eroded
// WM: label, or fraction >= 0.9). Raw voxels <= 0 are MRIQC's rotation mask: excluded from
// EFC and from FBER's background.
//
// --mask <brain> (e.g. brainchop-mindgrab, a SynthStrip-like brain+CSF mask) counts every mask
// voxel that is not GM or WM as CSF (label 0 -> CSF; with --pve the CSF fraction takes the
// remainder 1 - csf - gm - wm). MRIQC segments inside SynthStrip's mask, so its CSF includes
// the extra-cerebral CSF a brain parcellation leaves unlabeled (icvs_csf ~0.27 vs ~0.02).

#ifdef HAVE_QC

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core.h" // set_input_hdr, nifti_image_change_datatype
#ifdef HAVE_N4
#include "n4/n4.h" // N4 bias correction in the MRIQC intensity pipeline
#endif
#if defined(HAVE_ALLINEATE) && defined(HAVE_CONFORM)
#define QC_HAVE_AIR
#include "core32.h" // nii_qc_air_masks_f32
#endif

// statsmodels.robust.scale.mad default Gaussian-consistency constant: the
// normalised MAD is median(|x-center|) / MAD_C (== 1.4826 * raw MAD).
#define MAD_C 0.6744897501960817
// Minimum voxels a tissue must retain (after any erosion fallback) to trust its
// statistics; below this the tissue's metrics are reported as nan.
#define QC_MIN_VOX 100
// WM2MAX percentile (MRIQC uses the 99.95th percentile of the whole image).
#define QC_WM2MAX_PCTL 0.9995

typedef struct {
	double mean, stdv, median, mad, p05, p95, kurt;
	double n;    // voxels (or summed fractions) used for the statistics
	double nraw; // full (un-eroded) tissue voxel count or summed fractions, for ICV/volume
	int ok;    // 1 if n >= QC_MIN_VOX and stats are finite
} TissueStats;

static int qc_finite(double v) {
	return v >= -DBL_MAX && v <= DBL_MAX;
}

// ---- comparator-free quickselect (WASM-safe; project rule: no qsort+comparator
// on large/per-voxel arrays). Hoare partition with the middle element as pivot;
// deterministic, in-place, reorders `a`. ----
static float qc_select_kth(float *a, int n, int k) {
	int lo = 0, hi = n - 1;
	while (lo < hi) {
		float pivot = a[lo + (hi - lo) / 2];
		int i = lo, j = hi;
		while (i <= j) {
			while (a[i] < pivot) i++;
			while (a[j] > pivot) j--;
			if (i <= j) {
				float t = a[i]; a[i] = a[j]; a[j] = t;
				i++; j--;
			}
		}
		if (k <= j) hi = j;
		else if (k >= i) lo = i;
		else break;
	}
	return a[k];
}

// numpy-default ('linear') percentile of a[0..n-1], destroying order. p in [0,1].
static double qc_pctl_inplace(float *a, int n, double p) {
	if (n <= 0) return NAN;
	if (n == 1) return a[0];
	double r = p * (n - 1);
	int lo = (int)floor(r);
	if (lo < 0) lo = 0;
	if (lo > n - 1) lo = n - 1;
	double frac = r - lo;
	float vlo = qc_select_kth(a, n, lo);
	if (frac <= 0.0 || lo >= n - 1) return vlo;
	// After selecting index lo, a[lo+1..n-1] are all >= vlo; the next order
	// statistic is their minimum.
	float vhi = a[lo + 1];
	for (int i = lo + 2; i < n; i++)
		if (a[i] < vhi) vhi = a[i];
	return (double)vlo + frac * ((double)vhi - (double)vlo);
}

// Compute the per-tissue estimators selected by qc_plan.md. `vals` holds the
// (post-erosion) intensities; `scratch` is reusable and at least `n` floats.
static void qc_stats(const float *vals, int n, double nraw, float *scratch, TissueStats *st) {
	memset(st, 0, sizeof(*st));
	st->n = n;
	st->nraw = nraw;
	st->median = st->mad = st->mean = st->stdv = st->p05 = st->p95 = st->kurt = NAN;
	if (n < QC_MIN_VOX) { st->ok = 0; return; }
	// mean / population stdv (matches statsmodels DescrStatsW.std, ddof=0)
	double sum = 0.0;
	for (int i = 0; i < n; i++) sum += vals[i];
	double mean = sum / n;
	double m2 = 0.0, m4 = 0.0;
	for (int i = 0; i < n; i++) {
		double d = vals[i] - mean;
		double d2 = d * d;
		m2 += d2;
		m4 += d2 * d2;
	}
	m2 /= n; m4 /= n;
	st->mean = mean;
	st->stdv = sqrt(m2);
	// scipy.stats.kurtosis default: Fisher (excess), bias=True
	st->kurt = (m2 > 0.0) ? (m4 / (m2 * m2) - 3.0) : NAN;
	// Selection mutates order but preserves the multiset, so one copy serves all
	// three quantiles.
	memcpy(scratch, vals, (size_t)n * sizeof(float));
	st->median = qc_pctl_inplace(scratch, n, 0.50);
	st->p05 = qc_pctl_inplace(scratch, n, 0.05);
	st->p95 = qc_pctl_inplace(scratch, n, 0.95);
	// normalised MAD about the median
	for (int i = 0; i < n; i++) scratch[i] = (float)fabs(vals[i] - st->median);
	double rawmad = qc_pctl_inplace(scratch, n, 0.50);
	st->mad = rawmad / MAD_C;
	st->ok = qc_finite(st->median) && qc_finite(st->stdv) && st->stdv >= 0.0;
}

// IEEE bits of f, remapped so unsigned order is float order (memcpy: no aliasing UB).
static uint32_t qc_sort_key(float f) {
	uint32_t u;
	memcpy(&u, &f, sizeof u);
	return u ^ ((u >> 31) ? 0xFFFFFFFFu : 0x80000000u);
}

// Sort v[] ascending carrying w[] along: LSD radix on qc_sort_key (comparator-free, as
// the quickselect above). tv/tw are n-float scratch.
static void qc_sort_pairs(float *v, float *w, int n, float *tv, float *tw) {
	for (int shift = 0; shift < 32; shift += 8) {
		int count[257] = {0};
		for (int i = 0; i < n; i++) count[((qc_sort_key(v[i]) >> shift) & 0xFF) + 1]++;
		for (int b = 0; b < 256; b++) count[b + 1] += count[b];
		for (int i = 0; i < n; i++) {
			int o = count[(qc_sort_key(v[i]) >> shift) & 0xFF]++;
			tv[o] = v[i];
			tw[o] = w[i];
		}
		memcpy(v, tv, (size_t)n * sizeof(float));
		memcpy(w, tw, (size_t)n * sizeof(float));
	}
}

// Weighted NumPy-linear quantile of sorted v (weights w, total W); see the header.
static double qc_wpctl(const float *v, const float *w, int n, double W, double p) {
	double t = p * (W - 1.0), r = 0.0;
	if (t <= 0.0) return v[0];
	for (int k = 0; k < n - 1; k++) {
		double r1 = r + w[k];
		if (t < r1) return (double)v[k] + (t - r) / w[k] * ((double)v[k + 1] - (double)v[k]);
		r = r1;
	}
	return v[n - 1];
}

// Fraction-weighted tissue statistics over the n voxels with w > 0 (vals/w in voxel
// order, both reordered on return). Moments are summed in voxel order, as qc_stats
// does, so the binary case matches it bit for bit.
static void qc_wstats(float *vals, float *w, int n, float *tv, float *tw, TissueStats *st) {
	memset(st, 0, sizeof(*st));
	st->median = st->mad = st->mean = st->stdv = st->p05 = st->p95 = st->kurt = NAN;
	double W = 0.0, sum = 0.0, wmax = 0.0;
	for (int i = 0; i < n; i++) {
		W += w[i];
		sum += (double)w[i] * vals[i];
		if (w[i] > wmax) wmax = w[i];
	}
	st->n = st->nraw = W;
	if (W < QC_MIN_VOX) { st->ok = 0; return; }
	double mean = sum / W, m2 = 0.0;
	for (int i = 0; i < n; i++) {
		double d = vals[i] - mean;
		m2 += w[i] * d * d;
	}
	st->mean = mean;
	st->stdv = sqrt(m2 / W);
	// kurtosis and MAD: unweighted, over the core voxels (w > half the peak fraction)
	int nc = 0;
	double csum = 0.0;
	for (int i = 0; i < n; i++)
		if (w[i] > 0.5 * wmax) { tv[nc++] = vals[i]; csum += vals[i]; }
	double cmean = csum / nc, c2 = 0.0, c4 = 0.0;
	for (int i = 0; i < nc; i++) {
		double d = tv[i] - cmean, d2 = d * d;
		c2 += d2;
		c4 += d2 * d2;
	}
	c2 /= nc; c4 /= nc;
	st->kurt = (c2 > 0.0) ? (c4 / (c2 * c2) - 3.0) : NAN;
	qc_sort_pairs(vals, w, n, tw + n, tw); // tw holds 2n floats: temp weights, then temp values
	st->median = qc_wpctl(vals, w, n, W, 0.50);
	st->p05 = qc_wpctl(vals, w, n, W, 0.05);
	st->p95 = qc_wpctl(vals, w, n, W, 0.95);
	for (int i = 0; i < nc; i++) tv[i] = (float)fabs(tv[i] - st->median);
	st->mad = qc_pctl_inplace(tv, nc, 0.50) / MAD_C;
	st->ok = qc_finite(st->median) && qc_finite(st->stdv);
}

// Load an image and convert in place while respecting scl_slope/scl_inter.
static nifti_image *qc_read_as(const char *fname, int datatype, const char *type_name) {
	nifti_image *nim = nifti_image_read(fname, 1);
	if (!nim) { printf("qc: failed to read '%s'\n", fname); return NULL; }
	in_hdr ihdr = set_input_hdr(nim);
	if (nifti_image_change_datatype(nim, datatype, &ihdr) != 0) {
		printf("qc: failed to convert '%s' to %s\n", fname, type_name);
		nifti_image_free(nim);
		return NULL;
	}
	return nim;
}

// Parse a comma-separated list of integer labels into set[] (returns count, -1 on error).
static int qc_parse_labels(const char *s, int *set, int maxn) {
	int n = 0;
	const char *p = s;
	while (1) {
		while (*p == ' ' || *p == '\t') p++;
		if (!*p) break;
		char *end;
		errno = 0;
		long v = strtol(p, &end, 10);
		if (end == p || errno == ERANGE || v < INT_MIN || v > INT_MAX) {
			printf("qc: invalid integer in label list '%s'\n", s);
			return -1;
		}
		if (v == 0) {
			printf("qc: label 0 is reserved for background and cannot be CSF or WM\n");
			return -1;
		}
		int duplicate = 0;
		for (int i = 0; i < n; i++) if (set[i] == (int)v) duplicate = 1;
		if (!duplicate) {
			if (n >= maxn) { printf("qc: too many labels in '%s'\n", s); return -1; }
			set[n++] = (int)v;
		}
		p = end;
		while (*p == ' ' || *p == '\t') p++;
		if (!*p) break;
		if (*p != ',') { printf("qc: invalid label list '%s'\n", s); return -1; }
		p++;
		while (*p == ' ' || *p == '\t') p++;
		if (!*p) { printf("qc: trailing comma in label list '%s'\n", s); return -1; }
	}
	if (n == 0) { printf("qc: empty label list\n"); return -1; }
	for (int i = 1; i < n; i++) {
		int v = set[i], j = i - 1;
		while (j >= 0 && set[j] > v) { set[j + 1] = set[j]; j--; }
		set[j + 1] = v;
	}
	return n;
}

static int qc_in_set(int v, const int *set, int n) {
	int lo = 0, hi = n - 1;
	while (lo <= hi) {
		int mid = lo + (hi - lo) / 2;
		if (set[mid] == v) return 1;
		if (set[mid] < v) lo = mid + 1;
		else hi = mid - 1;
	}
	return 0;
}

// One-voxel 6-connected erosion (erode=1) or dilation (erode=0) of one class in the
// compact class map, into a binary mask. Border voxels erode away.
static void qc_morph6(const uint8_t *classes, uint8_t tissue, uint8_t *out,
                      int nx, int ny, int nz, int erode) {
	size_t nxy = (size_t)nx * ny;
	for (int z = 0; z < nz; z++)
		for (int y = 0; y < ny; y++)
			for (int x = 0; x < nx; x++) {
				size_t idx = (size_t)z * nxy + (size_t)y * nx + x;
				int self = classes[idx] == tissue;
				if (erode && !self) { out[idx] = 0; continue; }
				int xm = x > 0 && classes[idx - 1] == tissue, xp = x < nx - 1 && classes[idx + 1] == tissue;
				int ym = y > 0 && classes[idx - nx] == tissue, yp = y < ny - 1 && classes[idx + nx] == tissue;
				int zm = z > 0 && classes[idx - nxy] == tissue, zp = z < nz - 1 && classes[idx + nxy] == tissue;
				out[idx] = erode ? (xm && xp && ym && yp && zm && zp) : (self || xm || xp || ym || yp || zm || zp);
			}
}

// Gather one class (or the binary eroded mask) into caller-owned storage.
static int qc_gather(const float *t1, const uint8_t *map, uint8_t selected,
                     size_t nvox, float *out) {
	int n = 0;
	for (size_t i = 0; i < nvox; i++)
		if (map[i] == selected) out[n++] = t1[i];
	return n;
}

// One-voxel erosion with the 18-connected element (scipy generate_binary_structure(3, 2)),
// border voxels eroding away as in scipy's binary_erosion.
static void qc_erode18(const uint8_t *m, uint8_t v, uint8_t *out, int nx, int ny, int nz) {
	size_t nxy = (size_t)nx * ny;
	memset(out, 0, (size_t)nx * ny * nz);
	for (int z = 1; z < nz - 1; z++)
		for (int y = 1; y < ny - 1; y++)
			for (int x = 1; x < nx - 1; x++) {
				size_t i = (size_t)z * nxy + (size_t)y * nx + x;
				int keep = m[i] == v;
				for (int dz = -1; keep && dz <= 1; dz++)
					for (int dy = -1; keep && dy <= 1; dy++)
						for (int dx = -1; keep && dx <= 1; dx++)
							if (abs(dx) + abs(dy) + abs(dz) < 3 && m[i + dz * (long)nxy + dy * (long)nx + dx] != v) keep = 0;
				out[i] = keep;
			}
}

// niworkflows IntensityClip(p_min=10, p_max=99.9) as MRIQC applies it before N4.
// Thresholds come from a radius-3 ball median filter of the positive voxels.
// ponytail: filter sampled on every 2nd voxel per axis (percentiles only); full grid if it ever matters.
static int qc_clip(float *t1, int nx, int ny, int nz, float *scratch, double *lo, double *hi) {
	size_t nxy = (size_t)nx * ny, nvox = nxy * nz;
	float ball[123];
	int n = 0;
	for (int z = 0; z < nz; z += 2)
		for (int y = 0; y < ny; y += 2)
			for (int x = 0; x < nx; x += 2) {
				int k = 0;
				for (int dz = -3; dz <= 3; dz++)
					for (int dy = -3; dy <= 3; dy++)
						for (int dx = -3; dx <= 3; dx++) {
							if (dx * dx + dy * dy + dz * dz > 9) continue;
							// scipy's default 'reflect' border mode
							int xx = x + dx, yy = y + dy, zz = z + dz;
							xx = xx < 0 ? -xx - 1 : xx >= nx ? 2 * nx - xx - 1 : xx;
							yy = yy < 0 ? -yy - 1 : yy >= ny ? 2 * ny - yy - 1 : yy;
							zz = zz < 0 ? -zz - 1 : zz >= nz ? 2 * nz - zz - 1 : zz;
							ball[k++] = t1[(size_t)zz * nxy + (size_t)yy * nx + xx];
						}
				float med = qc_select_kth(ball, k, k / 2);
				if (med > 0.0f) scratch[n++] = med;
			}
	if (n < 2) return 1;
	*hi = qc_pctl_inplace(scratch, n, 0.999);
	*lo = qc_pctl_inplace(scratch, n, 0.10);
	if (!(*hi > *lo)) return 1;
	double s = 255.0 / (*hi - *lo);
	for (size_t i = 0; i < nvox; i++) {
		double v = t1[i] < *lo ? *lo : t1[i] > *hi ? *hi : t1[i];
		t1[i] = (float)round((v - *lo) * s);
	}
	return 0;
}

// MRIQC Harmonize: scale so the eroded WM median is 1000 (un-eroded WM below 30 voxels).
static double qc_harmonize(float *t1, const uint8_t *map, uint8_t v, int nx, int ny, int nz, float *scratch) {
	size_t nvox = (size_t)nx * ny * nz;
	uint8_t *ero = (uint8_t *)malloc(nvox);
	if (!ero) return NAN;
	qc_erode18(map, v, ero, nx, ny, nz);
	int n = qc_gather(t1, ero, 1, nvox, scratch);
	free(ero);
	if (n < 30) n = qc_gather(t1, map, v, nvox, scratch);
	double med = qc_pctl_inplace(scratch, n, 0.5);
	if (!(med > 0.0)) return NAN;
	for (size_t i = 0; i < nvox; i++) t1[i] = (float)(t1[i] * (1000.0 / med));
	return med;
}

typedef struct {
	double lo, hi, wm_median; // clip range; WM median before harmonizing
	int n4;                   // N4 ran
} QcPrep;

// MRIQC's intensity pipeline: clip, N4 (when built; weighted by the brain, map != 0, correcting the
// whole image as MRIQC's weighted N4 does), harmonize on the voxels where map == wm.
static int qc_prep(float *t1, const uint8_t *map, uint8_t wm, int nx, int ny, int nz, float dx, float dy,
                   float dz, float *scratch, QcPrep *prep) {
	if (qc_clip(t1, nx, ny, nz, scratch, &prep->lo, &prep->hi)) {
		printf("qc: too few positive voxels to clip intensities (use --raw)\n");
		return 1;
	}
	prep->n4 = 0;
#ifdef HAVE_N4
	N4Params n4p;
	n4_default_params(&n4p); // ANTs defaults, shrink 4 (as MRIQC runs N4)
	if (nx < 2 * n4p.shrink || ny < 2 * n4p.shrink || nz < 2 * n4p.shrink)
		printf("qc: image too small for N4 at shrink %d; skipping bias correction\n", n4p.shrink);
	else {
		size_t nvox = (size_t)nx * ny * nz;
		for (size_t i = 0; i < nvox; i++) scratch[i] = map[i] ? 1.0f : 0.0f; // N4 confidence
		if (n4_correct(t1, NULL, scratch, nx, ny, nz, dx, dy, dz, &n4p, NULL)) {
			printf("qc: N4 bias correction failed (use --raw)\n");
			return 1;
		}
		prep->n4 = 1;
	}
#else
	(void)dx; (void)dy; (void)dz;
#endif
	prep->wm_median = qc_harmonize(t1, map, wm, nx, ny, nz, scratch);
	if (!qc_finite(prep->wm_median)) {
		printf("qc: no white matter to harmonize on (no --wm label voxels, or no PVE WM fraction >= 0.9); use --raw\n");
		return 1;
	}
	return 0;
}

// Fill one tissue, falling back to its raw class if erosion leaves too few voxels.
static void qc_tissue(const char *name, const float *t1, const uint8_t *classes,
                      uint8_t tissue, double nraw, uint8_t *eroded, size_t nvox,
                      int nx, int ny, int nz, int do_erode, float *vals,
                      float *scratch, TissueStats *st) {
	const uint8_t *use = classes;
	uint8_t selected = tissue;
	int n = 0;
	if (do_erode) {
		qc_morph6(classes, tissue, eroded, nx, ny, nz, 1);
		n = qc_gather(t1, eroded, 1, nvox, vals);
		if (n >= QC_MIN_VOX) { use = eroded; selected = 1; }
		else printf("qc: %s has %d voxels after erosion (< %d); using un-eroded mask for its stats\n",
		            name, n, QC_MIN_VOX);
	}
	if (use == classes) n = qc_gather(t1, use, selected, nvox, vals);
	qc_stats(vals, n, nraw, scratch, st);
	if (!st->ok)
		printf("qc: %s has %d usable voxels (< %d); its metrics reported as nan\n", name, n, QC_MIN_VOX);
}

// Gather one fraction map's w > 0 voxels, then its weighted statistics.
static void qc_pve_tissue(const char *name, const float *t1, const float *frac, size_t nvox,
                          float *vals, float *w, float *tv, float *tw, TissueStats *st) {
	int n = 0;
	for (size_t i = 0; i < nvox; i++)
		if (frac[i] > 0.0f) { vals[n] = t1[i]; w[n++] = frac[i]; }
	qc_wstats(vals, w, n, tv, tw, st);
	if (!st->ok)
		printf("qc: %s fractions sum to %.1f (< %d); its metrics reported as nan\n", name, st->n, QC_MIN_VOX);
}

// ---- air ("hat") metrics: MRIQC's ArtifactMask, on the RAS-canonical image ----
// The hat is the head-free air above MRIQC's landmarks, glabella [0,90,-14] and inion
// [0,-120,-14] (template mm), mapped to voxels and applied as MRIQC's two axis-aligned slabs:
// slices below the inion's, and the anterior half below the glabella's, count as head. (A
// tilted plane through both landmarks seemed more exact but took 9 % more air, and brighter.)
// Artifacts are hat voxels brighter than QC_AIR_ZSCORE MADs of the hat's positive voxels,
// outside the 10 % shell nearest the head, after a 6-connected opening; the background
// statistics are taken over air = hat - artifacts.
// The Otsu head mask is dilated QC_HEAD_DILATE voxels to match MRIQC's gradient-based one
// (sigma 4 smoothing puts a halo round the scalp): 0.86 -> 0.90 Dice, FBER 1.4x -> ~1.03x.
// ponytail: dilation in voxels; MRIQC's sigma scales with zooms, so mm-scale it if anisotropic data matters.
#define QC_HEAD_DILATE 4
// 1 / sqrt(2 / (4 - pi)): Dietrich's Rayleigh correction (mriqc DIETRICH_FACTOR).
#define QC_DIETRICH 0.6551364
#define QC_AIR_ZSCORE 10.0

typedef struct {
	int present;
	double qi_1, fber, cnr, snrd_csf, snrd_gm, snrd_wm, snrd_total;
	TissueStats bg;
} AirMetrics;

#ifdef QC_HAVE_AIR
// `scratch` is the caller's nvox-float buffer (the RAS permutation keeps the count).
static int qc_air(const nifti_image *nt1, const char *ftmpl,
                  const TissueStats *scsf, const TissueStats *sgm, const TissueStats *swm,
                  float *scratch, AirMetrics *am) {
	nifti_image *ras = NULL, *head = NULL, *dist = NULL;
	mat44 v2t;
	memset(am, 0, sizeof(*am));
	am->qi_1 = am->fber = am->cnr = am->snrd_csf = am->snrd_gm = am->snrd_wm = am->snrd_total = NAN;
	if (nii_qc_air_masks_f32(nt1, ftmpl, &ras, &head, &dist, &v2t)) {
		printf("qc: air mask pipeline failed (RAS / registration to '%s' / head mask)\n", ftmpl);
		return 1;
	}
	int rc = 1;
	int nx = (int)ras->nx, ny = (int)ras->ny, nz = (int)ras->nz;
	size_t nxy = (size_t)nx * ny, nvox = ras->nvox;
	float *img = (float *)ras->data;
	const float *ds = (const float *)dist->data;
	uint8_t *hat = (uint8_t *)calloc(nvox, 1), *art = (uint8_t *)calloc(nvox, 1),
	        *tmp = (uint8_t *)calloc(nvox, 1), *rot = (uint8_t *)calloc(nvox, 1),
	        *hd = (uint8_t *)malloc(nvox);
	float *vals = (float *)malloc(nvox * sizeof(float));
	if (!hat || !art || !tmp || !rot || !hd || !vals) {
		printf("qc: out of memory allocating air masks\n");
		goto done;
	}
	for (size_t o = 0; o < nvox; o++) hd[o] = ((const float *)head->data)[o] > 0;
	for (int d = 0; d < QC_HEAD_DILATE; d++) {
		qc_morph6(hd, 1, tmp, nx, ny, nz, 0);
		memcpy(hd, tmp, nvox);
	}
	// nii_qc sets raw voxels <= 0 to -1; -ras only permutes, so the flag survives
	for (size_t o = 0; o < nvox; o++)
		if (img[o] < 0.0f) { rot[o] = 1; img[o] = 0.0f; }
	// 1. the hat: MRIQC's slabs, int() truncation included.
	mat44 t2v = nifti_mat44_inverse(v2t);
	int kInion = (int)(t2v.m[2][0] * 0 + t2v.m[2][1] * -120 + t2v.m[2][2] * -14 + t2v.m[2][3]);
	int kGlabella = (int)(t2v.m[2][0] * 0 + t2v.m[2][1] * 90 + t2v.m[2][2] * -14 + t2v.m[2][3]);
	long nHat = 0;
	float distMax = 0;
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++) {
			if (k < kInion || (j >= ny / 2 && k < kGlabella)) continue;
			size_t base = (size_t)j * nx + (size_t)k * nxy;
			for (int i = 0; i < nx; i++) {
				size_t o = base + i;
				if (hd[o]) continue;
				hat[o] = 1;
				nHat++;
				if (ds[o] > distMax) distMax = ds[o];
			}
		}
	if (nHat < QC_MIN_VOX) { // e.g. skull-stripped input: no usable background
		printf("qc: fewer than %d air voxels above the landmark plane; air metrics omitted\n", QC_MIN_VOX);
		rc = 0;
		goto done;
	}
	am->present = 1;
	// 2. artifacts, flagged against the MAD of the hat's positive voxels.
	int n = 0;
	for (size_t o = 0; o < nvox; o++)
		if (hat[o] && img[o] > 0) vals[n++] = img[o];
	TissueStats s0;
	qc_stats(vals, n, n, scratch, &s0);
	if (s0.mad > 0 && distMax > 0)
		for (size_t o = 0; o < nvox; o++)
			if (hat[o] && img[o] > 0 && ds[o] / distMax >= 0.1f && img[o] / s0.mad > QC_AIR_ZSCORE)
				art[o] = 1;
	qc_morph6(art, 1, tmp, nx, ny, nz, 1); // open: erode, then dilate
	qc_morph6(tmp, 1, art, nx, ny, nz, 0);
	long nArt = 0;
	for (size_t o = 0; o < nvox; o++)
		if (art[o] && hat[o]) { hat[o] = 0; nArt++; }
	am->qi_1 = (double)nArt / nHat;
	// 3. background statistics over the pruned air.
	n = qc_gather(img, hat, 1, nvox, vals);
	qc_stats(vals, n, n, scratch, &am->bg);
	if (n < 1) { rc = 0; goto done; }
	// 4. Dietrich SNR: median over the air MAD (stdv when the MAD is degenerate).
	double sigma = am->bg.mad > 1.0 ? am->bg.mad : am->bg.stdv;
	if (sigma > 1e-3) {
		const TissueStats *ts[3] = {scsf, sgm, swm};
		double *out[3] = {&am->snrd_csf, &am->snrd_gm, &am->snrd_wm};
		double sum = 0;
		int cnt = 0;
		for (int t = 0; t < 3; t++) {
			if (!qc_finite(ts[t]->median)) continue;
			*out[t] = QC_DIETRICH * ts[t]->median / sigma;
			sum += *out[t];
			cnt++;
		}
		if (cnt) am->snrd_total = sum / cnt;
	}
	// 5. FBER: median squared intensity inside the head over the same in all the air
	// outside the head (MRIQC: not the hat, but excluding the rotation mask).
	int nFg = 0, nBg = 0;
	for (size_t o = 0; o < nvox; o++) {
		if (hd[o]) scratch[nFg++] = img[o] * img[o];
		else if (!rot[o]) vals[nBg++] = img[o] * img[o];
	}
	if (nFg && nBg) {
		double fg = qc_pctl_inplace(scratch, nFg, 0.5);
		double bg = qc_pctl_inplace(vals, nBg, 0.5);
		am->fber = bg < 1e-3 ? -1.0 : fg / bg;
	}
	// 6. CNR with the air term (MRIQC's definition; sigma_air is <1 % of it in practice).
	if (swm->ok && sgm->ok && qc_finite(am->bg.stdv))
		am->cnr = fabs(swm->median - sgm->median) /
		          sqrt(am->bg.stdv * am->bg.stdv + sgm->stdv * sgm->stdv + swm->stdv * swm->stdv);
	rc = 0;
done:
	free(hat); free(art); free(tmp); free(rot); free(hd); free(vals);
	nifti_image_free(ras); nifti_image_free(head); nifti_image_free(dist);
	return rc;
}
#endif // QC_HAVE_AIR

// ---- output: one ordered value table feeds both the TSV and the JSON ----
typedef struct { char key[32]; double v; int is_count; } QcVal;
#define QC_MAX_VALS 64 // 53 with --air; the count is fixed at compile time

static void qc_push(QcVal *tab, int *n, const char *key, double v) {
	assert(*n < QC_MAX_VALS);
	snprintf(tab[*n].key, sizeof tab->key, "%s", key);
	tab[*n].v = v;
	tab[*n].is_count = 0;
	(*n)++;
}
static void qc_push_count(QcVal *tab, int *n, const char *key, long v) {
	qc_push(tab, n, key, (double)v);
	tab[*n - 1].is_count = 1;
}

// summary_<tissue>_* in MRIQC's order.
static void qc_push_tissue(QcVal *tab, int *n, const char *tn, const TissueStats *st, int soft) {
	const char *suf[7] = {"mean", "stdv", "median", "mad", "p05", "p95", "k"};
	double v[7] = {st->mean, st->stdv, st->median, st->mad, st->p05, st->p95, st->kurt};
	char key[32];
	for (int i = 0; i < 7; i++) {
		snprintf(key, sizeof key, "summary_%s_%s", tn, suf[i]);
		qc_push(tab, n, key, v[i]);
	}
	snprintf(key, sizeof key, "summary_%s_n", tn);
	if (soft) qc_push(tab, n, key, st->n);
	else qc_push_count(tab, n, key, (long)st->n);
}

static int qc_write_tsv(const char *fout, const QcVal *vals, int n) {
	FILE *f = fopen(fout, "w");
	if (!f) { printf("qc: cannot open output '%s'\n", fout); return 1; }
	for (int i = 0; i < n; i++) fprintf(f, "%s%s", i ? "\t" : "", vals[i].key);
	fputc('\n', f);
	for (int i = 0; i < n; i++) {
		if (i) fputc('\t', f);
		// counts print exactly: %.6g would round 1030101 to 1.0301e+06
		if (vals[i].is_count) fprintf(f, "%ld", (long)vals[i].v);
		else if (qc_finite(vals[i].v)) fprintf(f, "%.6g", vals[i].v);
		else fputs("nan", f);
	}
	fputc('\n', f);
	int bad = ferror(f);
	bad |= fclose(f);
	if (bad) printf("qc: failed while writing output '%s'\n", fout);
	return bad;
}

// Print a string as a JSON literal; the values here are file names from argv.
static void qc_json_str(FILE *f, const char *s) {
	fputc('"', f);
	for (; *s; s++) {
		if ((unsigned char)*s < 0x20) { fprintf(f, "\\u%04x", (unsigned char)*s); continue; }
		if (*s == '"' || *s == '\\') fputc('\\', f);
		fputc(*s, f);
	}
	fputc('"', f);
}

// MRIQC-style report: metrics flat at the top level (JSON has no NaN, so non-finite
// values are null), image geometry, and provenance. Mirrors MRIQC's <sub>_T1w.json.
static int qc_write_json(const char *fout, const QcVal *vals, int n, const nifti_image *nt1,
                         const int *csf, int ncsf, const int *wm, int nwm, int pve, const char *ftmpl,
                         const QcPrep *prep) {
	FILE *f = fopen(fout, "w");
	if (!f) { printf("qc: cannot open output '%s'\n", fout); return 1; }
	fputs("{\n", f);
	for (int i = 0; i < n; i++) {
		fprintf(f, "  \"%s\": ", vals[i].key);
		if (vals[i].is_count) fprintf(f, "%ld", (long)vals[i].v);
		else if (qc_finite(vals[i].v)) fprintf(f, "%.15g", vals[i].v);
		else fputs("null", f);
		fputs(",\n", f);
	}
	fprintf(f, "  \"size_x\": %lld,\n  \"size_y\": %lld,\n  \"size_z\": %lld,\n",
	        (long long)nt1->nx, (long long)nt1->ny, (long long)nt1->nz);
	fprintf(f, "  \"spacing_x\": %.15g,\n  \"spacing_y\": %.15g,\n  \"spacing_z\": %.15g,\n",
	        (double)nt1->dx, (double)nt1->dy, (double)nt1->dz);
	fputs("  \"provenance\": {\n    \"software\": \"niimath --qc\",\n", f);
	if (pve) fputs("    \"pve\": true", f);
	else {
		fputs("    \"csf_labels\": [", f);
		for (int i = 0; i < ncsf; i++) fprintf(f, "%s%d", i ? ", " : "", csf[i]);
		fputs("],\n    \"wm_labels\": [", f);
		for (int i = 0; i < nwm; i++) fprintf(f, "%s%d", i ? ", " : "", wm[i]);
		fputs("]", f);
	}
	if (prep) fprintf(f, ",\n    \"intensity\": {\"clip\": [%.15g, %.15g], \"n4\": %s, \"wm_median\": %.15g}",
	                  prep->lo, prep->hi, prep->n4 ? "true" : "false", prep->wm_median);
	else fputs(",\n    \"intensity\": \"raw\"", f);
	if (ftmpl) { fputs(",\n    \"air_template\": ", f); qc_json_str(f, ftmpl); }
	fputs("\n  }\n}\n", f);
	int bad = ferror(f);
	bad |= fclose(f);
	if (bad) printf("qc: failed while writing output '%s'\n", fout);
	return bad;
}

static void qc_usage(void) {
	printf("  usage: niimath --qc <t1> --seg <seg> --csf <i[,j..]> --wm <i[,j..]> [--erode 0|1] [--mask <brain>] [--air <template>] [--raw] [--out qc.tsv] [--json qc.json]\n");
	printf("     or: niimath --qc <t1> --pve <csf> <gm> <wm> [--mask <brain>] [--air <template>] [--raw] [--out qc.tsv] [--json qc.json]\n");
	printf("  --mask: brain(+CSF) mask; its voxels not GM or WM count as CSF (MRIQC's extra-cerebral CSF)\n");
	printf("  --raw: raw intensities instead of MRIQC's clip (0..255) + N4 + WM-median-1000 scaling\n");
}

// `other` must be a single 3D image on the T1's voxel grid.
static int qc_same_grid(const nifti_image *nt1, const nifti_image *other, const char *what) {
	int nvox3d = 0;
	if (nii_nvox3d_int(other, &nvox3d) || other->nvox != nvox3d) {
		printf("qc: %s must be a single 3D image of at most INT_MAX voxels\n", what);
		return 0;
	}
	if (other->nx != nt1->nx || other->ny != nt1->ny || other->nz != nt1->nz) {
		printf("qc: %s %lldx%lldx%lld does not match image %lldx%lldx%lld\n", what,
		       (long long)other->nx, (long long)other->ny, (long long)other->nz,
		       (long long)nt1->nx, (long long)nt1->ny, (long long)nt1->nz);
		return 0;
	}
	// No separate unit-code equality check: max_displacement_mm() normalises each
	// transform to mm via its own xyz_units, so a physically identical grid stored
	// with different but valid unit codes (e.g. an mm image and a 0.001 m image)
	// yields ~0 displacement and must be accepted, not rejected.
	float grid_mm = max_displacement_mm(nt1, other);
	if (!(grid_mm >= 0.0f && grid_mm <= 0.001f)) {
		printf("qc: %s spatial grid differs from image (maximum corner displacement %.6g mm)\n",
		       what, grid_mm);
		return 0;
	}
	return 1;
}

int nii_qc(int argc, char *argv[]) {
	const char *ft1 = NULL, *fseg = NULL, *csfstr = NULL, *wmstr = NULL;
	const char *fout = NULL, *fjson = NULL, *ftmpl = NULL, *fmask = NULL, *fpve[3] = {NULL, NULL, NULL};
	nifti_image *nt1 = NULL, *nseg = NULL, *nmask = NULL, *npve[3] = {NULL, NULL, NULL};
	float *wts = NULL, *tv = NULL, *tw = NULL;
	uint8_t *classes = NULL, *eroded = NULL, *rot = NULL;
	float *vals = NULL, *scratch = NULL;
	int rc = EXIT_FAILURE;
	int do_erode = 0, erode_set = 0, raw = 0;
	QcPrep prep = {NAN, NAN, NAN, 0};
	for (int i = 2; i < argc; i++) {
		const char *a = argv[i];
		if ((!strcmp(a, "-i") || !strcmp(a, "--in")) && i + 1 < argc) ft1 = argv[++i];
		else if ((!strcmp(a, "-seg") || !strcmp(a, "--seg")) && i + 1 < argc) fseg = argv[++i];
		else if ((!strcmp(a, "-csf") || !strcmp(a, "--csf")) && i + 1 < argc) csfstr = argv[++i];
		else if ((!strcmp(a, "-wm") || !strcmp(a, "--wm")) && i + 1 < argc) wmstr = argv[++i];
		else if ((!strcmp(a, "-erode") || !strcmp(a, "--erode")) && i + 1 < argc) {
			const char *v = argv[++i];
			if (!strcmp(v, "0")) do_erode = 0;
			else if (!strcmp(v, "1")) do_erode = 1;
			else { printf("qc: --erode must be 0 or 1 (got '%s')\n", v); return EXIT_FAILURE; }
			erode_set = 1;
		}
		else if ((!strcmp(a, "-pve") || !strcmp(a, "--pve")) && i + 3 < argc) {
			fpve[0] = argv[++i]; fpve[1] = argv[++i]; fpve[2] = argv[++i];
		}
		else if ((!strcmp(a, "-o") || !strcmp(a, "-out") || !strcmp(a, "--out")) && i + 1 < argc) fout = argv[++i];
		else if ((!strcmp(a, "-json") || !strcmp(a, "--json")) && i + 1 < argc) fjson = argv[++i];
		else if ((!strcmp(a, "-air") || !strcmp(a, "--air")) && i + 1 < argc) ftmpl = argv[++i];
		else if ((!strcmp(a, "-mask") || !strcmp(a, "--mask")) && i + 1 < argc) fmask = argv[++i];
		else if (!strcmp(a, "-raw") || !strcmp(a, "--raw")) raw = 1;
		else if (a[0] != '-' && !ft1) ft1 = a; // positional T1
		else {
			printf("qc: unsupported option '%s'\n", a);
			qc_usage();
			return EXIT_FAILURE;
		}
	}
	int pve = fpve[0] != NULL;
	if (pve && (fseg || csfstr || wmstr || erode_set)) {
		printf("qc: --pve replaces --seg/--csf/--wm/--erode\n");
		return EXIT_FAILURE;
	}
	if (!ft1 || (!pve && (!fseg || !csfstr || !wmstr))) {
		printf("qc: missing required arguments.\n");
		qc_usage();
		printf("  labels: 0 = non-brain (excluded); --csf/--wm list the CSF/WM label values; every other non-zero label is GM.\n");
		return EXIT_FAILURE;
	}
	if (!fout && !fjson) fout = "qc.tsv";
#ifndef QC_HAVE_AIR
	if (ftmpl) {
		printf("qc: --air needs a build with allineate and conform (registration and -ras)\n");
		return EXIT_FAILURE;
	}
#endif
	int csfset[64], wmset[64], ncsf = 0, nwm = 0;
	if (!pve) {
		ncsf = qc_parse_labels(csfstr, csfset, 64);
		nwm = qc_parse_labels(wmstr, wmset, 64);
		if (ncsf < 0 || nwm < 0) return EXIT_FAILURE;
	}
	for (int i = 0; i < ncsf; i++) {
		if (qc_in_set(csfset[i], wmset, nwm)) {
			printf("qc: label %d appears in both CSF and WM sets\n", csfset[i]);
			return EXIT_FAILURE;
		}
	}

	nt1 = qc_read_as(ft1, DT_FLOAT32, "float32");
	if (!nt1) goto cleanup;
	TissueStats scsf, swm, sgm;
	int nx = (int)nt1->nx, ny = (int)nt1->ny, nz = (int)nt1->nz;
	float mm = (float)xyz_units_to_mm(nt1->xyz_units); // N4's spacings are mm
	int nvox3d = 0;
	if (nii_nvox3d_int(nt1, &nvox3d) || nt1->nvox != nvox3d) {
		printf("qc: only single 3D images of at most INT_MAX voxels are supported\n");
		goto cleanup;
	}
	if (!raw && (nx < 3 || ny < 3 || nz < 3)) { // qc_clip's median filter reflects once at the borders
		printf("qc: MRIQC intensity processing needs at least 3 voxels per axis (got %dx%dx%d); use --raw\n", nx, ny, nz);
		goto cleanup;
	}
	size_t nvox = (size_t)nvox3d;
	float *t1 = (float *)nt1->data;
	rot = (uint8_t *)malloc(nvox);
	if (!rot) {
		printf("qc: out of memory allocating masks\n");
		goto cleanup;
	}
	for (size_t i = 0; i < nvox; i++) {
		if (!qc_finite(t1[i])) {
			printf("qc: T1 contains a non-finite value at voxel %zu\n", i);
			goto cleanup;
		}
		rot[i] = t1[i] <= 0.0f;
	}
	if (fmask) {
		nmask = qc_read_as(fmask, DT_FLOAT32, "float32");
		if (!nmask || !qc_same_grid(nt1, nmask, "brain mask")) goto cleanup;
	}
	const float *brain = nmask ? (const float *)nmask->data : NULL;
	if (pve) {
		size_t maxn = 1;
		for (int t = 0; t < 3; t++) {
			npve[t] = qc_read_as(fpve[t], DT_FLOAT32, "float32");
			if (!npve[t] || !qc_same_grid(nt1, npve[t], "partial-volume map")) goto cleanup;
			float *f = (float *)npve[t]->data;
			for (size_t i = 0; i < nvox; i++) {
				if (!qc_finite(f[i])) {
					printf("qc: partial-volume map '%s' has a non-finite value at voxel %zu\n", fpve[t], i);
					goto cleanup;
				}
				f[i] = f[i] < 0.0f ? 0.0f : f[i] > 1.0f ? 1.0f : f[i];
			}
		}
		if (brain) {
			float *fc = (float *)npve[0]->data;
			const float *fg = (const float *)npve[1]->data, *fw = (const float *)npve[2]->data;
			for (size_t i = 0; i < nvox; i++)
				if (brain[i] > 0.0f && fc[i] + fg[i] + fw[i] < 1.0f) fc[i] = 1.0f - fg[i] - fw[i];
		}
		for (int t = 0; t < 3; t++) {
			const float *f = (const float *)npve[t]->data;
			size_t n = 0;
			for (size_t i = 0; i < nvox; i++) n += f[i] > 0.0f;
			if (n > maxn) maxn = n;
		}
		vals = (float *)malloc(maxn * sizeof(float));
		wts = (float *)malloc(maxn * sizeof(float));
		tv = (float *)malloc(maxn * sizeof(float));
		tw = (float *)malloc(2 * maxn * sizeof(float));
		scratch = (float *)malloc(nvox * sizeof(float));
		if (!vals || !wts || !tv || !tw || !scratch) {
			printf("qc: out of memory allocating tissue statistics buffers\n");
			goto cleanup;
		}
		if (!raw) { // one map: 0 no tissue, 1 any fraction (N4's weight), 2 WM fraction >= 0.9 (harmonize)
			const float *fc = (const float *)npve[0]->data, *fg = (const float *)npve[1]->data,
			            *fwm = (const float *)npve[2]->data;
			uint8_t *map = (uint8_t *)malloc(nvox);
			if (!map) {
				printf("qc: out of memory allocating masks\n");
				goto cleanup;
			}
			for (size_t i = 0; i < nvox; i++) map[i] = fwm[i] >= 0.9f ? 2 : fc[i] + fg[i] + fwm[i] > 0.0f;
			int bad = qc_prep(t1, map, 2, nx, ny, nz, nt1->dx * mm, nt1->dy * mm, nt1->dz * mm, scratch, &prep);
			free(map);
			if (bad) goto cleanup;
		}
		qc_pve_tissue("CSF", t1, (float *)npve[0]->data, nvox, vals, wts, tv, tw, &scsf);
		qc_pve_tissue("GM", t1, (float *)npve[1]->data, nvox, vals, wts, tv, tw, &sgm);
		qc_pve_tissue("WM", t1, (float *)npve[2]->data, nvox, vals, wts, tv, tw, &swm);
		goto metrics;
	}
	// Keep labels in float64: converting an integer label map through float32
	// aliases distinct labels above 2^24 and defeats the integer validation below.
	nseg = qc_read_as(fseg, DT_FLOAT64, "float64");
	if (!nseg) goto cleanup;

	if (!qc_same_grid(nt1, nseg, "segmentation")) goto cleanup;
	const double *seg = (const double *)nseg->data;

	// Compact classification: 0 background, 1 CSF, 2 GM, 3 WM.
	classes = (uint8_t *)calloc(nvox, 1);
	if (!classes) {
		printf("qc: out of memory allocating tissue masks\n");
		goto cleanup;
	}
	long nraw_csf = 0, nraw_gm = 0, nraw_wm = 0;
	for (size_t i = 0; i < nvox; i++) {
		double sv = seg[i];
		if (!(sv >= INT_MIN && sv <= INT_MAX) || sv != trunc(sv)) {
			printf("qc: segmentation must contain finite integer labels (voxel %zu is %.17g)\n", i, sv);
			goto cleanup;
		}
		int lbl = (int)sv;
		if (lbl == 0) {
			if (brain && brain[i] > 0.0f) { classes[i] = 1; nraw_csf++; }
			continue;
		}
		if (qc_in_set(lbl, csfset, ncsf)) { classes[i] = 1; nraw_csf++; }
		else if (qc_in_set(lbl, wmset, nwm)) { classes[i] = 3; nraw_wm++; }
		else { classes[i] = 2; nraw_gm++; }
	}
	// The float64 segmentation is dead after classification; free it before the
	// peak allocations below so its nvox*8 bytes are not held alongside the
	// erosion/scratch/value buffers (the high-water mark matters for WASM).
	nifti_image_free(nseg); nseg = NULL; seg = NULL;
	eroded = (uint8_t *)calloc(nvox, 1);
	long max_tissue = nraw_csf;
	if (nraw_gm > max_tissue) max_tissue = nraw_gm;
	if (nraw_wm > max_tissue) max_tissue = nraw_wm;
	size_t stat_bytes = 0, scratch_bytes = 0;
	if (nii_mul_size((size_t)(max_tissue > 0 ? max_tissue : 1), sizeof(float), &stat_bytes)) {
		printf("qc: tissue statistics allocation size overflow\n");
		goto cleanup;
	}
	if (nii_mul_size(nvox, sizeof(float), &scratch_bytes)) {
		printf("qc: percentile scratch allocation size overflow\n");
		goto cleanup;
	}
	vals = (float *)malloc(stat_bytes);
	scratch = (float *)malloc(scratch_bytes);
	if (!eroded || !vals || !scratch) {
		printf("qc: out of memory allocating tissue statistics buffers\n");
		goto cleanup;
	}

	if (!raw && qc_prep(t1, classes, 3, nx, ny, nz, nt1->dx * mm, nt1->dy * mm, nt1->dz * mm, scratch, &prep)) goto cleanup;
	qc_tissue("CSF", t1, classes, 1, nraw_csf, eroded, nvox, nx, ny, nz,
	          do_erode, vals, scratch, &scsf);
	qc_tissue("WM", t1, classes, 3, nraw_wm, eroded, nvox, nx, ny, nz,
	          do_erode, vals, scratch, &swm);
	qc_tissue("GM", t1, classes, 2, nraw_gm, eroded, nvox, nx, ny, nz,
	          do_erode, vals, scratch, &sgm);

metrics:;
	// ---- metrics ----
	double cjv = NAN, cnr = NAN;
	double dmu = fabs(swm.median - sgm.median);
	if (swm.ok && sgm.ok && dmu > 1e-9) {
		cjv = (swm.mad + sgm.mad) / dmu;
		double cnr_denom = sqrt(swm.stdv * swm.stdv + sgm.stdv * sgm.stdv);
		if (cnr_denom > 0.0) cnr = dmu / cnr_denom; // s_bg = 0
	}
	double snr_csf = NAN, snr_wm = NAN, snr_gm = NAN, snr_total = NAN;
	if (scsf.ok && scsf.stdv > 0 && scsf.n > 1) snr_csf = scsf.median / (scsf.stdv * sqrt((double)scsf.n / (scsf.n - 1)));
	if (swm.ok && swm.stdv > 0 && swm.n > 1) snr_wm = swm.median / (swm.stdv * sqrt((double)swm.n / (swm.n - 1)));
	if (sgm.ok && sgm.stdv > 0 && sgm.n > 1) snr_gm = sgm.median / (sgm.stdv * sqrt((double)sgm.n / (sgm.n - 1)));
	if (qc_finite(snr_csf) && qc_finite(snr_wm) && qc_finite(snr_gm))
		snr_total = (snr_csf + snr_wm + snr_gm) / 3.0;

	// WM2MAX: median_WM / P99.95(whole image)
	double wm2max = NAN;
	{
		memcpy(scratch, t1, nvox * sizeof(float));
		double p9995 = qc_pctl_inplace(scratch, (int)nvox, QC_WM2MAX_PCTL);
		if (swm.ok && qc_finite(p9995) && fabs(p9995) > 1e-12)
			wm2max = swm.median / p9995;
	}

	// EFC over the voxels outside the rotation mask (raw > 0; MRIQC's framemask).
	double efc = NAN;
	{
		long N = 0;
		double sumsq = 0.0;
		for (size_t i = 0; i < nvox; i++) if (!rot[i]) { N++; sumsq += (double)t1[i] * t1[i]; }
		double bmax = sqrt(sumsq);
		if (N > 1 && bmax > 0.0) {
			double efc_max = (double)N * (1.0 / sqrt((double)N)) * log(1.0 / sqrt((double)N));
			double s = 0.0;
			for (size_t i = 0; i < nvox; i++)
				if (!rot[i]) s += (t1[i] / bmax) * log(((double)t1[i] + 1e-16) / bmax);
			if (efc_max != 0.0) efc = s / efc_max;
		}
	}

	// ICV fractions and absolute volumes use FULL (un-eroded) tissue counts.
	double icv_total = (double)(scsf.nraw + sgm.nraw + swm.nraw);
	double icvs_csf = icv_total > 0 ? scsf.nraw / icv_total : NAN;
	double icvs_gm = icv_total > 0 ? sgm.nraw / icv_total : NAN;
	double icvs_wm = icv_total > 0 ? swm.nraw / icv_total : NAN;
	double unit_to_mm = 1.0;
	if (nt1->xyz_units == NIFTI_UNITS_METER) unit_to_mm = 1000.0;
	else if (nt1->xyz_units == NIFTI_UNITS_MICRON) unit_to_mm = 0.001;
	else if (nt1->xyz_units == NIFTI_UNITS_UNKNOWN)
		printf("qc: spatial units are unspecified; assuming voxel dimensions are millimetres\n");
	else if (nt1->xyz_units != NIFTI_UNITS_MM) {
		printf("qc: unsupported spatial units code %d; absolute volumes reported as nan\n", nt1->xyz_units);
		unit_to_mm = NAN;
	}
	double vvol = fabs((double)nt1->dx * nt1->dy * nt1->dz) *
	              unit_to_mm * unit_to_mm * unit_to_mm;
	if (!(vvol > 0.0 && qc_finite(vvol))) vvol = NAN;
	double vol_csf = scsf.nraw * vvol, vol_gm = sgm.nraw * vvol, vol_wm = swm.nraw * vvol;

	AirMetrics air = {0};
#ifdef QC_HAVE_AIR
	if (ftmpl) {
		for (size_t i = 0; i < nvox; i++) if (rot[i]) t1[i] = -1.0f; // rotation-mask flag for qc_air
		if (qc_air(nt1, ftmpl, &scsf, &sgm, &swm, scratch, &air)) goto cleanup;
	}
#endif

	// ---- collect, then write ----
	QcVal tab[QC_MAX_VALS];
	int nvals = 0;
	qc_push(tab, &nvals, "cjv", cjv);
	qc_push(tab, &nvals, "cnr_noair", cnr);
	if (air.present) qc_push(tab, &nvals, "cnr", air.cnr);
	qc_push(tab, &nvals, "snr_csf", snr_csf); qc_push(tab, &nvals, "snr_wm", snr_wm);
	qc_push(tab, &nvals, "snr_gm", snr_gm); qc_push(tab, &nvals, "snr_total", snr_total);
	if (air.present) {
		qc_push(tab, &nvals, "snrd_csf", air.snrd_csf); qc_push(tab, &nvals, "snrd_wm", air.snrd_wm);
		qc_push(tab, &nvals, "snrd_gm", air.snrd_gm); qc_push(tab, &nvals, "snrd_total", air.snrd_total);
		qc_push(tab, &nvals, "fber", air.fber); qc_push(tab, &nvals, "qi_1", air.qi_1);
	}
	qc_push(tab, &nvals, "wm2max", wm2max); qc_push(tab, &nvals, "efc_brain", efc);
	qc_push(tab, &nvals, "icvs_csf", icvs_csf); qc_push(tab, &nvals, "icvs_gm", icvs_gm);
	qc_push(tab, &nvals, "icvs_wm", icvs_wm);
	qc_push(tab, &nvals, "vol_csf_mm3", vol_csf); qc_push(tab, &nvals, "vol_gm_mm3", vol_gm);
	qc_push(tab, &nvals, "vol_wm_mm3", vol_wm);
	qc_push_tissue(tab, &nvals, "csf", &scsf, pve);
	qc_push_tissue(tab, &nvals, "gm", &sgm, pve);
	qc_push_tissue(tab, &nvals, "wm", &swm, pve);
	if (air.present) qc_push_tissue(tab, &nvals, "bg", &air.bg, 0);

	if (fout && qc_write_tsv(fout, tab, nvals)) goto cleanup;
	if (fjson && qc_write_json(fjson, tab, nvals, nt1, csfset, ncsf, wmset, nwm, pve,
	                           air.present ? ftmpl : NULL, raw ? NULL : &prep)) goto cleanup;

	if (fout) printf("qc: wrote %s\n", fout);
	if (fjson) printf("qc: wrote %s\n", fjson);
	printf("qc: CJV=%.4g  CNR(noair)=%.4g  SNR(total)=%.4g  WM2MAX=%.4g  EFC=%.4g\n",
	       cjv, cnr, snr_total, wm2max, efc);
	if (air.present)
		printf("qc: SNRd(total)=%.4g  FBER=%.4g  QI1=%.4g  CNR=%.4g\n",
		       air.snrd_total, air.fber, air.qi_1, air.cnr);
	if (pve) printf("qc: summed fractions CSF=%.1f GM=%.1f WM=%.1f\n", scsf.nraw, sgm.nraw, swm.nraw);
	else printf("qc: voxels CSF=%.0f GM=%.0f WM=%.0f (erode=%d)\n", scsf.nraw, sgm.nraw, swm.nraw, do_erode);
	rc = EXIT_SUCCESS;

cleanup:
	free(classes);
	free(eroded);
	free(rot);
	free(vals);
	free(scratch);
	free(wts); free(tv); free(tw);
	for (int t = 0; t < 3; t++) if (npve[t]) nifti_image_free(npve[t]);
	if (nmask) nifti_image_free(nmask);
	if (nt1) nifti_image_free(nt1);
	if (nseg) nifti_image_free(nseg);
	return rc;
}

#endif // HAVE_QC

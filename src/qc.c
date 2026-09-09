// niimath --qc : MRIQC-style anatomical quality metrics from a T1 image + segmentation.
//
// Given a T1-weighted intensity image, a matching integer tissue segmentation,
// and the label values that denote CSF and white matter, this computes the
// subset of MRIQC anatomical Image Quality Metrics (IQMs) that do NOT require a
// background/air noise distribution (our backgrounds are masked to zero) and do
// NOT require soft partial-volume maps (we have a hard label map). The AGENTS.md
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
// the masked voxels (the binary special case). To emulate the boundary
// suppression MRIQC gets from soft pvms, tissue masks are eroded one voxel before
// the intensity statistics (toggle with --erode 0). ICV fractions and absolute
// volumes use the FULL (un-eroded) tissue extent.

#ifdef HAVE_QC

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core.h" // set_input_hdr, nifti_image_change_datatype
#ifdef HAVE_ALLINEATE
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
	long n;    // voxels used for the statistics (post-erosion when eroded)
	long nraw; // full (un-eroded) tissue voxel count, for ICV/volume
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
static void qc_stats(const float *vals, int n, long nraw, float *scratch, TissueStats *st) {
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

// One-voxel 6-connected erosion of one class in the compact class map.
static void qc_erode6(const uint8_t *classes, uint8_t tissue, uint8_t *out,
                      int nx, int ny, int nz) {
	size_t nxy = (size_t)nx * ny;
	for (int z = 0; z < nz; z++)
		for (int y = 0; y < ny; y++)
			for (int x = 0; x < nx; x++) {
				size_t idx = (size_t)z * nxy + (size_t)y * nx + x;
				if (classes[idx] != tissue) { out[idx] = 0; continue; }
				int keep = (x > 0 && classes[idx - 1] == tissue) &&
				           (x < nx - 1 && classes[idx + 1] == tissue) &&
				           (y > 0 && classes[idx - nx] == tissue) &&
				           (y < ny - 1 && classes[idx + nx] == tissue) &&
				           (z > 0 && classes[idx - nxy] == tissue) &&
				           (z < nz - 1 && classes[idx + nxy] == tissue);
				out[idx] = keep ? 1 : 0;
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

// Fill one tissue, falling back to its raw class if erosion leaves too few voxels.
static void qc_tissue(const char *name, const float *t1, const uint8_t *classes,
                      uint8_t tissue, long nraw, uint8_t *eroded, size_t nvox,
                      int nx, int ny, int nz, int do_erode, float *vals,
                      float *scratch, TissueStats *st) {
	const uint8_t *use = classes;
	uint8_t selected = tissue;
	int n = 0;
	if (do_erode) {
		qc_erode6(classes, tissue, eroded, nx, ny, nz);
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

// ---- air ("hat") metrics: MRIQC's ArtifactMask, on the RAS-canonical image ----
// The hat is the head-free air superior to MRIQC's landmark plane: its glabella
// [0,90,-14] and inion [0,-120,-14] share template z = -14, so the two axis-aligned
// slab fills MRIQC uses approximate that one plane, tested here directly (exact under
// pitch, roll and yaw). Artifacts are hat voxels brighter than QC_AIR_ZSCORE MADs of
// the full hat, outside the 10 % shell nearest the head, after a 6-connected opening;
// the background statistics are taken over air = hat - artifacts.
#define QC_LANDMARK_PLANE_Z (-14.0)
// 1 / sqrt(2 / (4 - pi)): Dietrich's Rayleigh correction (mriqc DIETRICH_FACTOR).
#define QC_DIETRICH 0.6551364
#define QC_AIR_ZSCORE 10.0

typedef struct {
	int present;
	double qi_1, fber, cnr, snrd_csf, snrd_gm, snrd_wm, snrd_total;
	TissueStats bg;
} AirMetrics;

#ifdef HAVE_ALLINEATE
static int qc_air(const nifti_image *nt1, const char *ftmpl, const char *ft1,
                  const TissueStats *scsf, const TissueStats *sgm, const TissueStats *swm,
                  AirMetrics *am) {
	nifti_image *ras = NULL, *head = NULL, *dist = NULL;
	mat44 v2t;
	memset(am, 0, sizeof(*am));
	am->qi_1 = am->fber = am->cnr = am->snrd_csf = am->snrd_gm = am->snrd_wm = am->snrd_total = NAN;
	if (nii_qc_air_masks_f32(nt1, ftmpl, ft1, &ras, &head, &dist, &v2t)) {
		printf("qc: air mask pipeline failed (RAS / registration to '%s' / head mask)\n", ftmpl);
		return 1;
	}
	int rc = 1;
	int nx = (int)ras->nx, ny = (int)ras->ny, nz = (int)ras->nz;
	size_t nxy = (size_t)nx * ny, nvox = ras->nvox;
	const float *img = (const float *)ras->data, *hd = (const float *)head->data,
	            *ds = (const float *)dist->data;
	uint8_t *hat = (uint8_t *)calloc(nvox, 1), *art = (uint8_t *)calloc(nvox, 1),
	        *tmp = (uint8_t *)calloc(nvox, 1);
	float *vals = (float *)malloc(nvox * sizeof(float)), *scratch = (float *)malloc(nvox * sizeof(float));
	if (!hat || !art || !tmp || !vals || !scratch) {
		printf("qc: out of memory allocating air masks\n");
		goto done;
	}
	// 1. the hat. Only the z row of voxel->template matters: a plane equation in (i,j,k).
	long nHat = 0;
	float distMax = 0;
	for (int k = 0; k < nz; k++)
		for (int j = 0; j < ny; j++) {
			double zjk = v2t.m[2][1] * j + v2t.m[2][2] * k + v2t.m[2][3];
			size_t base = (size_t)j * nx + (size_t)k * nxy;
			for (int i = 0; i < nx; i++) {
				if (v2t.m[2][0] * i + zjk < QC_LANDMARK_PLANE_Z) continue;
				size_t o = base + i;
				if (hd[o] > 0) continue;
				hat[o] = 1;
				nHat++;
				if (ds[o] > distMax) distMax = ds[o];
			}
		}
	if (nHat < 10) { // e.g. skull-stripped input: no usable background
		printf("qc: fewer than 10 air voxels above the landmark plane; air metrics omitted\n");
		rc = 0;
		goto done;
	}
	am->present = 1;
	// 2. artifacts, flagged against the MAD of the FULL hat.
	int n = qc_gather(img, hat, 1, nvox, vals);
	TissueStats s0;
	qc_stats(vals, n, n, scratch, &s0);
	if (s0.mad > 0 && distMax > 0)
		for (size_t o = 0; o < nvox; o++)
			if (hat[o] && img[o] > 0 && ds[o] / distMax >= 0.1f && img[o] / s0.mad > QC_AIR_ZSCORE)
				art[o] = 1;
	qc_erode6(art, 1, tmp, nx, ny, nz);
	memset(art, 0, nvox);
	for (int k = 1; k < nz - 1; k++)
		for (int j = 1; j < ny - 1; j++)
			for (int i = 1; i < nx - 1; i++) {
				size_t o = (size_t)k * nxy + (size_t)j * nx + i;
				if (!tmp[o]) continue;
				art[o] = art[o - 1] = art[o + 1] = art[o - nx] = art[o + nx] = art[o - nxy] = art[o + nxy] = 1;
			}
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
	// 5. FBER: median squared intensity inside the head over the same in the air.
	int nFg = 0;
	for (size_t o = 0; o < nvox; o++)
		if (hd[o] > 0) scratch[nFg++] = img[o] * img[o];
	if (nFg) {
		double fg = qc_pctl_inplace(scratch, nFg, 0.5);
		for (int i = 0; i < n; i++) vals[i] *= vals[i];
		double bg = qc_pctl_inplace(vals, n, 0.5);
		am->fber = bg < 1e-3 ? -1.0 : fg / bg;
	}
	// 6. CNR with the air term (MRIQC's definition; sigma_air is <1 % of it in practice).
	if (swm->ok && sgm->ok && qc_finite(am->bg.stdv))
		am->cnr = fabs(swm->median - sgm->median) /
		          sqrt(am->bg.stdv * am->bg.stdv + sgm->stdv * sgm->stdv + swm->stdv * swm->stdv);
	rc = 0;
done:
	free(hat); free(art); free(tmp); free(vals); free(scratch);
	nifti_image_free(ras); nifti_image_free(head); nifti_image_free(dist);
	return rc;
}
#endif // HAVE_ALLINEATE

// ---- output: one ordered value table feeds both the TSV and the JSON ----
typedef struct { const char *key; double v; int is_count; } QcVal;
#define QC_MAX_VALS 96

static void qc_push(QcVal *vals, int *n, const char *key, double v) {
	if (*n < QC_MAX_VALS) { vals[*n].key = key; vals[*n].v = v; vals[*n].is_count = 0; (*n)++; }
}
static void qc_push_count(QcVal *vals, int *n, const char *key, long v) {
	if (*n < QC_MAX_VALS) { vals[*n].key = key; vals[*n].v = (double)v; vals[*n].is_count = 1; (*n)++; }
}

// summary_<tissue>_* in MRIQC's order. The key strings must outlive the table, so
// they are formatted into caller-owned storage.
static void qc_push_tissue(QcVal *vals, int *n, const char *tn, const TissueStats *st, char keys[8][32]) {
	const char *suf[8] = {"mean", "stdv", "median", "mad", "p05", "p95", "k", "n"};
	double v[7] = {st->mean, st->stdv, st->median, st->mad, st->p05, st->p95, st->kurt};
	for (int i = 0; i < 8; i++) snprintf(keys[i], 32, "summary_%s_%s", tn, suf[i]);
	for (int i = 0; i < 7; i++) qc_push(vals, n, keys[i], v[i]);
	qc_push_count(vals, n, keys[7], st->n);
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
	int bad = ferror(f) | fclose(f);
	if (bad) printf("qc: failed while writing output '%s'\n", fout);
	return bad;
}

// Print a string as a JSON literal; the values here are file names from argv.
static void qc_json_str(FILE *f, const char *s) {
	fputc('"', f);
	for (; *s; s++) {
		if (*s == '"' || *s == '\\') fputc('\\', f);
		fputc(*s, f);
	}
	fputc('"', f);
}

// MRIQC-style report: metrics flat at the top level (JSON has no NaN, so non-finite
// values are null), image geometry, and provenance. Mirrors MRIQC's <sub>_T1w.json.
static int qc_write_json(const char *fout, const QcVal *vals, int n, const nifti_image *nt1,
                         const int *csf, int ncsf, const int *wm, int nwm, const char *ftmpl) {
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
	fputs("  \"provenance\": {\n    \"software\": \"niimath --qc\",\n    \"csf_labels\": [", f);
	for (int i = 0; i < ncsf; i++) fprintf(f, "%s%d", i ? ", " : "", csf[i]);
	fputs("],\n    \"wm_labels\": [", f);
	for (int i = 0; i < nwm; i++) fprintf(f, "%s%d", i ? ", " : "", wm[i]);
	fputs("]", f);
	if (ftmpl) { fputs(",\n    \"air_template\": ", f); qc_json_str(f, ftmpl); }
	fputs("\n  }\n}\n", f);
	int bad = ferror(f) | fclose(f);
	if (bad) printf("qc: failed while writing output '%s'\n", fout);
	return bad;
}

int nii_qc(int argc, char *argv[]) {
	const char *ft1 = NULL, *fseg = NULL, *csfstr = NULL, *wmstr = NULL;
	const char *fout = NULL, *fjson = NULL, *ftmpl = NULL;
	nifti_image *nt1 = NULL, *nseg = NULL;
	uint8_t *classes = NULL, *eroded = NULL;
	float *vals = NULL, *scratch = NULL;
	int rc = EXIT_FAILURE;
	int do_erode = 1;
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
		}
		else if ((!strcmp(a, "-o") || !strcmp(a, "-out") || !strcmp(a, "--out")) && i + 1 < argc) fout = argv[++i];
		else if ((!strcmp(a, "-json") || !strcmp(a, "--json")) && i + 1 < argc) fjson = argv[++i];
		else if ((!strcmp(a, "-air") || !strcmp(a, "--air")) && i + 1 < argc) ftmpl = argv[++i];
		else if (a[0] != '-' && !ft1) ft1 = a; // positional T1
		else {
			printf("qc: unsupported option '%s'\n", a);
			printf("  usage: niimath --qc <t1> --seg <seg> --csf <i[,j..]> --wm <i[,j..]> [--erode 0|1] [--air <template>] [--out qc.tsv] [--json qc.json]\n");
			return EXIT_FAILURE;
		}
	}
	if (!ft1 || !fseg || !csfstr || !wmstr) {
		printf("qc: missing required arguments.\n");
		printf("  usage: niimath --qc <t1> --seg <seg> --csf <i[,j..]> --wm <i[,j..]> [--erode 0|1] [--air <template>] [--out qc.tsv] [--json qc.json]\n");
		printf("  labels: 0 = non-brain (excluded); --csf/--wm list the CSF/WM label values; every other non-zero label is GM.\n");
		return EXIT_FAILURE;
	}
	if (!fout && !fjson) fout = "qc.tsv";
#ifndef HAVE_ALLINEATE
	if (ftmpl) {
		printf("qc: --air needs a build with allineate (the head-to-template registration)\n");
		return EXIT_FAILURE;
	}
#endif
	int csfset[64], wmset[64];
	int ncsf = qc_parse_labels(csfstr, csfset, 64);
	int nwm = qc_parse_labels(wmstr, wmset, 64);
	if (ncsf < 0 || nwm < 0) return EXIT_FAILURE;
	for (int i = 0; i < ncsf; i++) {
		if (qc_in_set(csfset[i], wmset, nwm)) {
			printf("qc: label %d appears in both CSF and WM sets\n", csfset[i]);
			return EXIT_FAILURE;
		}
	}

	nt1 = qc_read_as(ft1, DT_FLOAT32, "float32");
	if (!nt1) goto cleanup;
	// Keep labels in float64: converting an integer label map through float32
	// aliases distinct labels above 2^24 and defeats the integer validation below.
	nseg = qc_read_as(fseg, DT_FLOAT64, "float64");
	if (!nseg) goto cleanup;

	int nvox3d = 0, nvox3d_seg = 0;
	if (nii_nvox3d_int(nt1, &nvox3d) || nii_nvox3d_int(nseg, &nvox3d_seg)) {
		printf("qc: invalid or oversized image dimensions (QC requires at most INT_MAX voxels)\n");
		goto cleanup;
	}
	if (nt1->nvox != nvox3d || nseg->nvox != nvox3d_seg) {
		printf("qc: only single 3D images are supported (image nvox=%lld, seg nvox=%lld)\n",
		       (long long)nt1->nvox, (long long)nseg->nvox);
		goto cleanup;
	}
	if (nseg->nx != nt1->nx || nseg->ny != nt1->ny || nseg->nz != nt1->nz ||
	    nvox3d_seg != nvox3d) {
		printf("qc: segmentation %lldx%lldx%lld does not match image %lldx%lldx%lld\n",
		       (long long)nseg->nx, (long long)nseg->ny, (long long)nseg->nz,
		       (long long)nt1->nx, (long long)nt1->ny, (long long)nt1->nz);
		goto cleanup;
	}
	// No separate unit-code equality check: max_displacement_mm() normalises each
	// transform to mm via its own xyz_units, so a physically identical grid stored
	// with different but valid unit codes (e.g. an mm image and a 0.001 m image)
	// yields ~0 displacement and must be accepted, not rejected.
	float grid_mm = max_displacement_mm(nt1, nseg);
	if (!(grid_mm >= 0.0f && grid_mm <= 0.001f)) {
		printf("qc: segmentation spatial grid differs from image (maximum corner displacement %.6g mm)\n",
		       grid_mm);
		goto cleanup;
	}

	int nx = (int)nt1->nx, ny = (int)nt1->ny, nz = (int)nt1->nz;
	size_t nvox = (size_t)nvox3d;
	const float *t1 = (const float *)nt1->data;
	const double *seg = (const double *)nseg->data;

	// Compact classification: 0 background, 1 CSF, 2 GM, 3 WM.
	classes = (uint8_t *)calloc(nvox, 1);
	if (!classes) {
		printf("qc: out of memory allocating tissue masks\n");
		goto cleanup;
	}
	long nraw_csf = 0, nraw_gm = 0, nraw_wm = 0;
	for (size_t i = 0; i < nvox; i++) {
		if (!qc_finite(t1[i])) {
			printf("qc: T1 contains a non-finite value at voxel %zu\n", i);
			goto cleanup;
		}
		double sv = seg[i];
		if (!(sv >= INT_MIN && sv <= INT_MAX) || sv != trunc(sv)) {
			printf("qc: segmentation must contain finite integer labels (voxel %zu is %.17g)\n", i, sv);
			goto cleanup;
		}
		int lbl = (int)sv;
		if (lbl == 0) continue;
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

	TissueStats scsf, swm, sgm;
	qc_tissue("CSF", t1, classes, 1, nraw_csf, eroded, nvox, nx, ny, nz,
	          do_erode, vals, scratch, &scsf);
	qc_tissue("WM", t1, classes, 3, nraw_wm, eroded, nvox, nx, ny, nz,
	          do_erode, vals, scratch, &swm);
	qc_tissue("GM", t1, classes, 2, nraw_gm, eroded, nvox, nx, ny, nz,
	          do_erode, vals, scratch, &sgm);

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

	// EFC over non-zero (in-brain-signal) voxels: framemask = (t1 == 0). The metric
	// assumes non-negative intensities; a mean-centered/negative T1 reports nan.
	double efc = NAN;
	{
		long N = 0;
		double sumsq = 0.0;
		for (size_t i = 0; i < nvox; i++) if (t1[i] != 0.0f) { N++; sumsq += (double)t1[i] * t1[i]; }
		double bmax = sqrt(sumsq);
		if (N > 1 && bmax > 0.0) {
			double efc_max = (double)N * (1.0 / sqrt((double)N)) * log(1.0 / sqrt((double)N));
			double s = 0.0;
			for (size_t i = 0; i < nvox; i++)
				if (t1[i] != 0.0f) s += (t1[i] / bmax) * log(((double)t1[i] + 1e-16) / bmax);
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
#ifdef HAVE_ALLINEATE
	if (ftmpl && qc_air(nt1, ftmpl, ft1, &scsf, &sgm, &swm, &air)) goto cleanup;
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
	char keys[4][8][32];
	qc_push_tissue(tab, &nvals, "csf", &scsf, keys[0]);
	qc_push_tissue(tab, &nvals, "gm", &sgm, keys[1]);
	qc_push_tissue(tab, &nvals, "wm", &swm, keys[2]);
	if (air.present) qc_push_tissue(tab, &nvals, "bg", &air.bg, keys[3]);

	if (fout && qc_write_tsv(fout, tab, nvals)) goto cleanup;
	if (fjson && qc_write_json(fjson, tab, nvals, nt1, csfset, ncsf, wmset, nwm,
	                           air.present ? ftmpl : NULL)) goto cleanup;

	if (fout) printf("qc: wrote %s\n", fout);
	if (fjson) printf("qc: wrote %s\n", fjson);
	printf("qc: CJV=%.4g  CNR(noair)=%.4g  SNR(total)=%.4g  WM2MAX=%.4g  EFC=%.4g\n",
	       cjv, cnr, snr_total, wm2max, efc);
	if (air.present)
		printf("qc: SNRd(total)=%.4g  FBER=%.4g  QI1=%.4g  CNR=%.4g\n",
		       air.snrd_total, air.fber, air.qi_1, air.cnr);
	printf("qc: voxels CSF=%ld GM=%ld WM=%ld (erode=%d)\n", scsf.nraw, sgm.nraw, swm.nraw, do_erode);
	rc = EXIT_SUCCESS;

cleanup:
	free(classes);
	free(eroded);
	free(vals);
	free(scratch);
	if (nt1) nifti_image_free(nt1);
	if (nseg) nifti_image_free(nseg);
	return rc;
}

#endif // HAVE_QC

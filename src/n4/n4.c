/* niimath -n4: C port of ITK/ANTs N4 bias-field correction.
 *
 *  Derived from the Insight Toolkit (Copyright 1999-2019 Insight Software Consortium,
 *  Copyright 2020-present NumFOCUS) and ANTs (Copyright the ANTs contributors).
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file
 *  except in compliance with the License. You may obtain a copy of the License at
 *      https://www.apache.org/licenses/LICENSE-2.0.txt  (also LICENSE in this directory)
 *  Unless required by applicable law or agreed to in writing, software distributed under the
 *  License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 *  either express or implied. See the License for the specific language governing
 *  permissions and limitations under the License. See NOTICE in this directory.
 */
#include "n4.h"
#include "n4_priv.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef N4_DUMP
// Parity-bench builds (plan_n4.md): with N4_DUMP=<dir>, write stage arrays as raw little-endian
// float32 <dir>/<name>.f32 (x fastest), matching the ITK dumps.
void n4_dump(const char *name, const float *v, size_t n) {
	const char *dir = getenv("N4_DUMP");
	if (!dir) return;
	char path[4096];
	snprintf(path, sizeof path, "%s/%s.f32", dir, name);
	FILE *f = fopen(path, "wb");
	if (!f) return;
	fwrite(v, sizeof(float), n, f);
	fclose(f);
}
#endif

// itk::ShrinkImageFilter: floor(n/s) samples (at least 1) whose centre matches the input's, so the
// first sample is input index round_half_up((n-1)/2 - s*(m-1)/2).
static void n4_shrink_axis(int n, int s, int *m, int *off) {
	*m = n / s > 0 ? n / s : 1;
	int o = (int)floor((n - 1) / 2.0 - s * (*m - 1) / 2.0 + 0.5);
	*off = o > 0 ? o : 0;
}

void n4_default_params(N4Params *p) {
	*p = (N4Params){.shrink = 4, .levels = 4, .iterations = {50, 50, 50, 50}, .convergence = 0.0,
	                .control_points = {4, 4, 4}, .bins = 200, .fwhm = 0.15, .wiener = 0.01,
	                .spline_distance = 0.0, .rescale = 1};
}

// itk::N4BiasFieldCorrectionImageFilter::SharpenImage: Wiener-deconvolve the histogram of the
// included voxels by a Gaussian of the bias field and map each voxel to E[u|v]. Types follow ITK
// (RealType float; FFT and complex arithmetic double) so the result tracks it to rounding.
static int n4_sharpen(const float *u, const uint8_t *incl, size_t nvox, const N4Params *p, float *out) {
	const int nb = p->bins;
	float binMax = -FLT_MAX, binMin = FLT_MAX;
	for (size_t i = 0; i < nvox; i++)
		if (incl[i]) {
			// ITK's else-if: a voxel that raises the maximum is never tested against the minimum
			if (u[i] > binMax) binMax = u[i];
			else if (u[i] < binMin) binMin = u[i];
		}
	const float slope = n4_fdiv(binMax - binMin, (float)(nb - 1));
	if (!(slope > 0.0f)) { // flat histogram (ITK divides by zero here): nothing to sharpen
		for (size_t i = 0; i < nvox; i++) out[i] = incl[i] ? u[i] : 0.0f;
		return 0;
	}
	const float exponent = (float)(ceil(logf((float)nb) / log(2.0)) + 1);
	const int np = (int)(powf(2.0f, exponent) + 0.5f), off = (int)(0.5 * (np - nb));
	float *H = (float *)calloc(nb, sizeof(float)), *E = (float *)malloc(np * sizeof(float));
	double *buf = (double *)calloc(10 * (size_t)np, sizeof(double));
	if (!H || !E || !buf) {
		free(H); free(E); free(buf);
		return 1;
	}
	double *Vr = buf, *Vi = Vr + np, *Fr = Vi + np, *Fi = Fr + np, *Ur = Fi + np, *Ui = Ur + np,
	       *Nr = Ui + np, *Ni = Nr + np, *Dr = Ni + np, *Di = Dr + np;
	for (size_t i = 0; i < nvox; i++) // triangular Parzen window
		if (incl[i]) {
			float cidx = n4_fdiv(u[i] - binMin, slope);
			unsigned idx = (unsigned)(int)floorf(cidx); // a sub-minimum voxel wraps and is dropped, as in ITK
			float frac = cidx - (float)idx;
			if (frac == 0.0f) H[idx] = (float)(H[idx] + 1.0);
			else if (idx < (unsigned)nb - 1) {
				H[idx] = (float)(H[idx] + (1.0 - frac));
				H[idx + 1] += frac;
			}
		}
	for (int n = 0; n < nb; n++) Vr[n + off] = H[n];
	n4_fft(Vr, Vi, np, 0);
	// Gaussian of the bias field, in bins
	const float scaledFWHM = n4_fdiv((float)p->fwhm, slope);
	const float expFactor = (float)(4.0 * log(2.0) / (scaledFWHM * scaledFWHM));
	const float scaleFactor = (float)(2.0 * sqrt(log(2.0) / M_PI) / scaledFWHM);
	Fr[0] = scaleFactor;
	const int half = np / 2;
	for (int n = 1; n <= half; n++)
		Fr[n] = Fr[np - n] = scaleFactor * expf(-((float)n * (float)n) * expFactor);
	if (np % 2 == 0) Fr[half] = scaleFactor * exp(0.25 * -((float)np * (float)np) * expFactor);
	n4_fft(Fr, Fi, np, 0);
	// Wiener deconvolution: only Re(G) = Re(conj(F) / (|F|^2 + noise)) is used
	const double noise = (double)(float)p->wiener;
	for (int n = 0; n < np; n++) {
		double g = Fr[n] / (Fr[n] * Fr[n] + Fi[n] * Fi[n] + noise);
		Ur[n] = Vr[n] * g;
		Ui[n] = Vi[n] * g;
	}
	n4_fft(Ur, Ui, np, 1);
	for (int n = 0; n < np; n++) {
		Ur[n] = Ur[n] > 0.0 ? Ur[n] : 0.0;
		Ui[n] = 0.0;
		Nr[n] = (binMin + ((float)n - (float)off) * slope) * Ur[n];
		Dr[n] = Ur[n];
	}
	// E[u|v] = (G * (u U)) / (G * U), convolved by multiplying with the Gaussian's transform
	n4_fft(Nr, Ni, np, 0);
	n4_fft(Dr, Di, np, 0);
	for (int n = 0; n < np; n++) {
		double r = Nr[n] * Fr[n] - Ni[n] * Fi[n];
		Ni[n] = Nr[n] * Fi[n] + Ni[n] * Fr[n];
		Nr[n] = r;
		r = Dr[n] * Fr[n] - Di[n] * Fi[n];
		Di[n] = Dr[n] * Fi[n] + Di[n] * Fr[n];
		Dr[n] = r;
	}
	n4_fft(Nr, Ni, np, 1);
	n4_fft(Dr, Di, np, 1);
	for (int n = 0; n < nb; n++) E[n] = Dr[n + off] != 0.0 ? (float)(Nr[n + off] / Dr[n + off]) : 0.0f;
#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < nvox; i++) {
		if (!incl[i]) {
			out[i] = 0.0f;
			continue;
		}
		float cidx = n4_fdiv(u[i] - binMin, slope);
		unsigned idx = (unsigned)(int)floorf(cidx);
		out[i] = idx < (unsigned)nb - 1 ? E[idx] + (E[idx + 1] - E[idx]) * (cidx - (float)idx) : E[nb - 1];
	}
	free(H); free(E); free(buf);
	return 0;
}

// itk::N4BiasFieldCorrectionImageFilter::CalculateConvergenceMeasurement: coefficient of
// variation of exp(old - new) over the included voxels.
static double n4_convergence(const float *oldf, const float *newf, const uint8_t *incl, size_t n) {
	double mu = 0.0, sigma = 0.0, N = 0.0;
	for (size_t k = 0; k < n; k++) {
		if (!incl[k]) continue;
		double v = exp(oldf[k] - newf[k]);
		N += 1.0;
		if (N > 1.0) sigma += (v - mu) * (v - mu) * (N - 1.0) / N;
		mu = mu * (1.0 - 1.0 / N) + v / N;
	}
	return sqrt(sigma / (N - 1.0)) / mu;
}

// The filter itself: voxels with mask 0 take no part and are left unchanged (ANTs -x; also the
// zero padding of -bspline <mm>); n4_correct applies the -rescale rule.
static int n4_core(float *data, const uint8_t *mask, const float *confidence, int nx, int ny, int nz,
                   float dx, float dy, float dz, const N4Params *p, float *field_out) {
	const int s = p->shrink > 0 ? p->shrink : 1, n[3] = {nx, ny, nz};
	int m[3], off[3];
	for (int a = 0; a < 3; a++) n4_shrink_axis(n[a], s, &m[a], &off[a]);
	if (m[0] < 2 || m[1] < 2 || m[2] < 2) { // the B-spline domain needs two samples per axis
		printf("n4: image %dx%dx%d is too small for shrink %d (each axis needs at least %d voxels)\n", nx, ny, nz, s, 2 * s);
		return 1;
	}
	const size_t nxy = (size_t)nx * ny, nvox = nxy * nz, snvox = (size_t)m[0] * m[1] * m[2];
	const double sp[3] = {(double)dx * s, (double)dy * s, (double)dz * s};
	int nc[3] = {p->control_points[0], p->control_points[1], p->control_points[2]};
	double maxlatd = 1.0; // lattice grows 2*nc-3 per refinement
	for (int d = 0; d < 3; d++) {
		double c = nc[d];
		for (int L = 1; L < p->levels; L++) c = 2.0 * c - 3.0;
		maxlatd *= c;
	}
	double partial = 1.0; // the fit keeps 4 final-level control planes per shrunk slice
	for (int d = 0; d < 2; d++) {
		double c = nc[d];
		for (int L = 1; L < p->levels; L++) c = 2.0 * c - 3.0;
		partial *= c;
	}
	partial *= 4.0 * m[2];
	if (maxlatd > N4_MAX_LATTICE || partial > N4_MAX_PARTIAL) {
		printf("n4: the B-spline lattice is too fine for this image (%.0f control points, %.0f fit partials; limits %d, %d); use fewer levels, a coarser -bspline or a larger -shrink\n",
		       maxlatd, partial, N4_MAX_LATTICE, N4_MAX_PARTIAL);
		return 1;
	}
	const size_t maxlat = (size_t)maxlatd;
	float *logInput = (float *)malloc(snvox * sizeof(float)), *logUnc = (float *)malloc(snvox * sizeof(float)),
	      *resid = (float *)malloc(snvox * sizeof(float)), *field = (float *)calloc(snvox, sizeof(float)),
	      *newField = (float *)malloc(snvox * sizeof(float)), *lattice = (float *)calloc(maxlat, sizeof(float)),
	      *phi = (float *)malloc(maxlat * sizeof(float)),
	      *sconf = confidence ? (float *)malloc(snvox * sizeof(float)) : NULL, // B-spline point weights
	      *full = field_out ? field_out : (float *)malloc(nvox * sizeof(float));
	uint8_t *incl = (uint8_t *)malloc(snvox);
	int rc = 1;
	if (!logInput || !logUnc || !resid || !field || !newField || !lattice || !phi || !incl || !full ||
	    (confidence && !sconf)) {
		printf("n4: out of memory\n");
		goto done;
	}
	size_t k = 0, nIncluded = 0;
	for (int z = 0; z < m[2]; z++)
		for (int y = 0; y < m[1]; y++)
			for (int x = 0; x < m[0]; x++, k++) {
				size_t i = (size_t)(z * s + off[2]) * nxy + (size_t)(y * s + off[1]) * nx + (x * s + off[0]);
				logInput[k] = data[i];
				// ITK: included iff mask != 0 and confidence > 0; only included voxels > 0 are logged
				incl[k] = (!mask || mask[i]) && (!confidence || confidence[i] > 0.0f);
				if (sconf) sconf[k] = confidence[i];
			}
	n4_dump("shrunk_input", logInput, snvox);
#ifdef N4_DUMP
	{ // the shrunk weight and mask, as ITK's dumps hold them
		float *w = (float *)malloc(snvox * sizeof(float));
		for (k = 0; w && k < snvox; k++) w[k] = sconf ? sconf[k] : 1.0f;
		if (w) n4_dump("shrunk_weight", w, snvox);
		for (k = 0; w && k < snvox; k++) {
			int z = (int)(k / ((size_t)m[0] * m[1])), y = (int)(k / m[0] % m[1]), x = (int)(k % m[0]);
			w[k] = (!mask || mask[(size_t)(z * s + off[2]) * nxy + (size_t)(y * s + off[1]) * nx + (x * s + off[0])]) ? 1.0f : 0.0f;
		}
		if (w) n4_dump("shrunk_mask", w, snvox);
		free(w);
	}
#endif
	for (k = 0; k < snvox; k++)
		if (incl[k]) {
			nIncluded++;
			if (logInput[k] > 0.0f) logInput[k] = logf(logInput[k]);
		}
	n4_dump("log_input", logInput, snvox);
	if (nIncluded < 2) {
		printf("n4: fewer than 2 voxels in the mask/weight region after shrinking\n");
		goto done;
	}
	memcpy(logUnc, logInput, snvox * sizeof(float));
	for (int L = 0; L < p->levels; L++) {
		double conv = DBL_MAX;
		for (int it = 0; it < p->iterations[L] && conv > p->convergence; it++) {
			const int first = L == 0 && it == 0; // bench dumps only
			(void)first;
			if (n4_sharpen(logUnc, incl, snvox, p, resid)) goto done;
			if (first) n4_dump("sharpened_L0I1", resid, snvox);
			for (k = 0; k < snvox; k++) resid[k] = logUnc[k] - resid[k];
			if (first) n4_dump("residual_L0I1", resid, snvox);
			if (n4_bspline_fit(resid, incl, sconf, m, sp, nc, phi)) {
				printf("n4: B-spline fit failed (a voxel fell outside the parametric domain, or out of memory)\n");
				goto done;
			}
			const size_t nlat = (size_t)nc[0] * nc[1] * nc[2];
			for (size_t l = 0; l < nlat; l++) lattice[l] += phi[l]; // lattice starts at zero
			if (first) n4_dump("lattice_L0I1", lattice, nlat);
			if (n4_bspline_eval(lattice, nc, m, sp, newField)) goto done;
			if (first) n4_dump("field_L0I1", newField, snvox);
			conv = n4_convergence(field, newField, incl, snvox);
			float *t = field; field = newField; newField = t;
			for (k = 0; k < snvox; k++) logUnc[k] = logInput[k] - field[k];
		}
#ifdef N4_DUMP
		char name[64];
		snprintf(name, sizeof name, "logfield_L%d", L);
		n4_dump(name, field, snvox);
#endif
		if (L < p->levels - 1) { // double the lattice resolution for the next level
			n4_bspline_refine(lattice, nc, phi);
			for (int d = 0; d < 3; d++) nc[d] = 2 * nc[d] - 3;
			float *t = lattice; lattice = phi; phi = t;
			if (L == 0) n4_dump("lattice_L0_refined", lattice, (size_t)nc[0] * nc[1] * nc[2]);
		}
	}
	n4_dump("lattice_final", lattice, (size_t)nc[0] * nc[1] * nc[2]);
	if (n4_bspline_eval(lattice, nc, n, (const double[3]){dx, dy, dz}, full)) goto done;
	n4_dump("logfield_full", full, nvox);
#pragma omp parallel for schedule(static)
	for (size_t i = 0; i < nvox; i++)
		if (!mask || mask[i]) data[i] /= expf(full[i]);
	rc = 0;
done:
	free(logInput); free(logUnc); free(resid); free(field); free(newField); free(lattice); free(phi);
	free(sconf); free(incl);
	if (full != field_out) free(full);
	return rc;
}

// ANTs N4BiasFieldCorrection around the filter: -b <distance> padding and the -x mask output rules.
int n4_correct(float *data, const uint8_t *mask, const float *confidence, int nx, int ny, int nz,
               float dx, float dy, float dz, const N4Params *p, float *field_out) {
	const size_t nvox = (size_t)nx * ny * nz;
	const float vox[3] = {dx, dy, dz};
	for (int d = 0; d < 3; d++) // a degenerate spacing makes the B-spline coordinates NaN
		if (!(isfinite(vox[d]) && vox[d] >= 1e-6f && vox[d] <= 1e6f)) {
			printf("n4: voxel size %g on axis %d is not usable\n", vox[d], d + 1);
			return 1;
		}
	double omin = DBL_MAX, omax = -DBL_MAX; // input range inside the mask, for -rescale
	if (mask && p->rescale)
		for (size_t i = 0; i < nvox; i++)
			if (mask[i]) {
				if (data[i] < omin) omin = data[i];
				if (data[i] > omax) omax = data[i];
			}
	int rc;
	if (!(p->spline_distance > 0.0)) rc = n4_core(data, mask, confidence, nx, ny, nz, dx, dy, dz, p, field_out);
	else {
		// ANTs: pad each axis so its domain holds whole spans of the requested length (float math as
		// the CLI), lower half first; the padding is zero and outside the mask; crop back afterwards.
		const int n[3] = {nx, ny, nz};
		const float sp[3] = {dx, dy, dz};
		int lo[3], pn[3];
		N4Params q = *p;
		for (int d = 0; d < 3; d++) { // bound the request in double before ANTs' float arithmetic
			double sd = ceil((double)(n[d] - 1) * sp[d] / p->spline_distance);
			if (!(sd >= 1.0 && sd <= 4096.0 && (sd * p->spline_distance) / sp[d] < INT_MAX / 2 - n[d])) {
				printf("n4: -bspline %g gives an unusable number of spans on axis %d\n", p->spline_distance, d + 1);
				return 1;
			}
		}
		for (int d = 0; d < 3; d++) {
			float domain = (float)(n[d] - 1) * sp[d];
			unsigned spans = (unsigned)ceilf(domain / (float)p->spline_distance);
			unsigned long extra = (unsigned long)((spans * (float)p->spline_distance - domain) / sp[d] + 0.5f);
			lo[d] = (int)(0.5 * extra);
			pn[d] = n[d] + (int)extra;
			q.control_points[d] = (int)spans + N4_SPLINE_ORDER;
		}
		if ((double)pn[0] * pn[1] * pn[2] > INT_MAX) {
			printf("n4: -bspline %g would pad the image beyond %d voxels\n", p->spline_distance, INT_MAX);
			return 1;
		}
		const size_t pnvox = (size_t)pn[0] * pn[1] * pn[2];
		float *pd = (float *)calloc(pnvox, sizeof(float)), *pc = confidence ? (float *)calloc(pnvox, sizeof(float)) : NULL;
		float *pf = field_out ? (float *)malloc(pnvox * sizeof(float)) : NULL;
		uint8_t *pm = (uint8_t *)calloc(pnvox, 1);
		rc = !pd || !pm || (confidence && !pc) || (field_out && !pf);
		if (!rc) {
			for (int z = 0; z < nz; z++)
				for (int y = 0; y < ny; y++)
					for (int x = 0; x < nx; x++) {
						size_t i = ((size_t)z * ny + y) * nx + x;
						size_t j = ((size_t)(z + lo[2]) * pn[1] + (y + lo[1])) * pn[0] + (x + lo[0]);
						pd[j] = data[i];
						pm[j] = mask ? mask[i] : 1;
						if (pc) pc[j] = confidence[i];
					}
			rc = n4_core(pd, pm, pc, pn[0], pn[1], pn[2], dx, dy, dz, &q, pf);
			if (!rc)
				for (int z = 0; z < nz; z++)
					for (int y = 0; y < ny; y++)
						for (int x = 0; x < nx; x++) {
							size_t i = ((size_t)z * ny + y) * nx + x;
							size_t j = ((size_t)(z + lo[2]) * pn[1] + (y + lo[1])) * pn[0] + (x + lo[0]);
							data[i] = pd[j];
							if (pf) field_out[i] = pf[j];
						}
		} else
			printf("n4: out of memory\n");
		free(pd); free(pc); free(pf); free(pm);
	}
	if (!rc && mask && p->rescale) { // ANTs -r: map the corrected range inside the mask onto the input's
		double cmin = DBL_MAX, cmax = -DBL_MAX;
		for (size_t i = 0; i < nvox; i++)
			if (mask[i]) {
				if (data[i] < cmin) cmin = data[i];
				if (data[i] > cmax) cmax = data[i];
			}
		if (cmax > cmin) {
			double slope = (omax - omin) / (cmax - cmin);
			for (size_t i = 0; i < nvox; i++)
				if (mask[i]) data[i] = (float)(omax - slope * (cmax - data[i]));
		}
	}
	return rc;
}

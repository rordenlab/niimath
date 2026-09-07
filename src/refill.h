// refill.h - REFILL dynamic distortion correction for niimath (-refill-gefm, -refill-epifm,
//            -refill-centre, -refill-unwarp)
//
// C port of the MATLAB reference implementation of REFILL (Robinson et al., "Improved dynamic
// distortion correction for fMRI using single-echo EPI and a readout-reversed first image",
// Hum Brain Mapp 2023; https://github.com/simon-mri/REFILL-Dynamic-Distortion-Correction, MIT).
// Each operation reproduces one MATLAB stage on the stored voxel grid; the contract, the oracle
// and the per-stage tolerances are in test/refill_reference_manifest.md.  Read that before
// changing any constant here.
//
// Guarded by HAVE_REFILL; requires HAVE_ROMEO (the unwrapping is niimath's own ROMEO port).
// Ordinary FP (rides the fast-math source line); DT32 only; byte-identical across thread counts.

#ifndef REFILL_H
#define REFILL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"

typedef struct {
	double clip_lo, clip_hi;  // field-map limits in rad/s; outside -> missing (default -600, 2000)
	double s;                 // smoothn smoothness (default 2); 0 = no smoothing/filling
	double qthresh;           // -refill-epifm: ROMEO quality threshold for the mask (default 0.5)
	int echo1, echo2;         // -refill-gefm: 1-based echoes used (default 1, 3)
	const char *steps;        // directory for MATLAB-named intermediates, or NULL
	int ramp_fix;             // -refill-epifm: apply the readout ramp with the paper's sign (Eq. 7) instead of the reference code's
} refill_opts;

refill_opts refill_opts_default(void);

// -refill-gefm <phase> <mask> <te1_ms> <te2_ms>: working image = multi-echo FLASH magnitude.
// Output: 3D static field map in rad/s, masked, filled and smoothed (= steps/ge_fm_masked.nii).
int refill_gefm(nifti_image *nim, const char *phasefile, const char *maskfile,
	double te1_ms, double te2_ms, const refill_opts *o, gzModes gzMode);

// -refill-epifm <phase> <refill_phase> <te_ms>: working image = EPI magnitude time series.
// Output: 4D dynamic field maps in rad/s (= steps/epi_fm_masked.nii); side outputs
// <out>_quality (ROMEO combined quality map) and <out>_mask (ROMEO robustmask).
int refill_epifm(nifti_image *nim, const char *phasefile, const char *refillfile,
	double te_ms, const refill_opts *o, gzModes gzMode);

// -refill-centre <mask> <dwell_s> <y|y->: working image = 4D field maps.  Unwarps them with
// their own voxel-shift map, takes the in-mask median over all volumes, subtracts it.
int refill_centre(nifti_image *nim, const char *maskfile, double dwell, const char *dir,
	const refill_opts *o, gzModes gzMode);

// -refill-unwarp <fm> <dwell_s> <y|y->: working image = any 3D/4D image on the field map's
// grid; forward-resampled along the phase-encoding axis (aspire_unwarp).
int refill_unwarp(nifti_image *nim, const char *fmfile, double dwell, const char *dir,
	const refill_opts *o, gzModes gzMode);

#ifdef __cplusplus
}
#endif

#endif // REFILL_H

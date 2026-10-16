// confounds.h - fMRIPrep confounds table (-confounds), for fastprep.
//
// Ports the computations fMRIPrep runs in Python (nipype ComputeDVARS / ACompCor / TCompCor,
// fMRIPrep FSLMotionParams / FSLRMSDeviation / FramewiseDisplacement, niworkflows ExpandModel /
// SpikeRegressors) onto the motion-corrected series in boldref space, and writes the
// desc-confounds_timeseries.tsv columns.  Built with -moco (HAVE_MOCO).

#ifndef CONFOUNDS_H
#define CONFOUNDS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"
#include "resample.h"

typedef struct {
	double tr;            // RepetitionTime (s)
	int dummy;            // leading non-steady-state volumes
	const double *hmc;    // nt rows of -moco -estimate parameters (required)
	const char *ref;      // boldref (motion-parameter centre of gravity)
	const char *brain;    // boldref-space brain mask
	const char *acc[3];   // T1w-space CSF, WM, CSF+WM aCompCor maps (GM-excluded), or NULL
	const char *t1mask;   // T1w-space brain mask (crown); with acc
	const double *ftm;    // allineate fixed_to_moving (T1w world -> boldref world); with acc
	const resample_vsm *vsm;   // susceptibility shift applied with the motion correction, or NULL
	const nifti_image *bwarp;  // BOLD-world greedy field (see resample.h), or NULL
	const char *tsnr;     // side output: tSNR map on the boldref grid, or NULL
	const char *tsnr_json;     // with tsnr: summary (median in brain, model)
} confounds_opts;

// `nim`: the (slice-time corrected) raw series, float32, untouched.  Returns 0 on success.
int nii_confounds(nifti_image *nim, const confounds_opts *o, const char *tsv);

#ifdef __cplusplus
}
#endif

#endif // CONFOUNDS_H

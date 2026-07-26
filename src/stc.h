// stc.h - slice-time correction for 4D datasets (-stc)
//
// Clean-room implementation of the slice-time correction published in AFNI's 3dTshift
// documentation: the default "detrend -> interpolate -> retrend" Fourier method.  AFNI's
// 3dTshift.c, its shifting engine and its FFT are Medical College of Wisconsin copyrighted and
// GPL-2; they were NOT read, translated or paraphrased.  They served only as a black-box oracle.
// Every convention the published help does not fix was measured through the 3dTshift executable
// and is recorded in the moco_bench repository's test/stc_reference_manifest.md.
//
// Guarded by HAVE_STC.

#ifndef STC_H
#define STC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"

// -stc --slicetiming <list|@file> [-tzero <seconds>]: ordinary chain operation.
//
// Shifts every voxel time series of a scalar 4D image so that all slices share one temporal
// origin.  `nim` must be float32 and 4D with nt >= 5.  `times` holds exactly `ntimes` slice
// acquisition times in SECONDS, one per storage-axis-k slice in slice-index order; `ntimes` must
// equal nz.  When `have_tzero` is 0 the common time point defaults to the arithmetic mean of
// `times`; otherwise `tzero` (also seconds) is used and must lie within [min(times), max(times)].
//
// A slice whose fractional shift falls below the skip threshold is copied VERBATIM, which is the
// measured oracle behaviour and is what makes an all-equal timing list a bit-exact no-op.  Every
// other slice is corrected, and there a time series holding any non-finite sample is written out
// as all-NaN (a Fourier shift is global, so no sample of such a series is defined).  The two
// rules meet at the skip threshold: a non-finite sample inside a skipped slice survives as
// itself, unchanged, because that slice is never transformed.
//
// On success `nim->data` is replaced by the corrected series and `nim->toffset` is set to tzero
// (expressed in the header's own time units).  Every other header field, including the spatial
// transforms, is untouched.  On any failure `nim` is left exactly as it was and a diagnostic is
// printed.  Returns 0 on success, non-zero on failure.
int nii_stc(nifti_image *nim, const double *times, int ntimes, int have_tzero, double tzero);

#ifdef __cplusplus
}
#endif

#endif // STC_H

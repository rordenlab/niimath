// fmap.h - B0 fieldmap EPI distortion correction (-fugue)
//
// Clean-room implementation of the observable behaviour of FSL's `fugue` in its
// "--loadfmap + --dwell + --unwarpdir" unwarping mode.  FSL is licensed under the University of
// Oxford's non-commercial licence, which is incompatible with niimath's BSD-2-Clause.  FSL's
// sources for fugue/prelude (fugue.cc, unwarpfns.cc, fsl_prepare_fieldmap.tcl) were NOT read,
// grepped, translated or paraphrased -- the executables served ONLY as a black-box oracle, driven
// with synthetic inputs constructed to isolate one convention at a time.  The method itself is
// published: Jezzard & Balaban, Magn Reson Med 34:65-73 (1995), plus FSL's own documentation at
// https://fsl.fmrib.ox.ac.uk/fsl/docs/registration/fugue.html
//
// Every measured convention, the experiment behind it, and the two deliberate divergences are
// recorded in the fmap_bench repository's test/fmap_reference_manifest.md.  Read that before
// changing any constant here.
//
// Guarded by HAVE_FMAP.

#ifndef FMAP_H
#define FMAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include "core.h"

// -fugue <fieldmap> <dwell> <unwarpdir>: ordinary chain operation, DT32 only.
//
// Corrects susceptibility-induced geometric distortion in an EPI image by resampling it along a
// single phase-encoding axis.  `nim` must be float32, 3D or 4D, and share a voxel grid with the
// fieldmap (dimensions and world transform).  `fmapfile` is a 3D B0 fieldmap in rad/s -- the units
// fsl_prepare_fieldmap writes -- and is applied identically to every volume of a 4D input.
// `dwell` is the effective echo spacing in SECONDS (BIDS EffectiveEchoSpacing).  `unwarpdir` is
// one of x y z (or the BIDS spellings i j k) with an optional trailing '-'.
//
//   s(v)   = fmap(v) / (2*pi) * dwell * N_axis                     [voxels]
//   out(v) = linear_interp_1d(in, v + sigma * s(v) * e_axis)       sigma = +1, or -1 for "a-"
//
// The shift is sampled at the OUTPUT voxel -- a pure pull, so signal is not conserved under
// compression, which is the reference behaviour.  Samples falling outside the field of view
// contribute 0.  No intensity (Jacobian) correction is applied; that is opt-in in the reference
// and is not implemented here.
//
// Voxels whose |fieldmap| falls below FMAP_EPS carry no data, and the shift field is extrapolated
// over them along each line parallel to the unwarp axis before resampling.  Skipping that step
// does not merely lose accuracy at the edges: it changes a third of the image.
//
// The working image is transformed IN PLACE, one line at a time, so peak memory is the input plus
// one shift volume rather than two copies of the input.  Byte-identical across thread counts:
// every line is independent and no reduction is used.
//
// Fails closed, leaving `nim` untouched, on: a non-float32 or >4D working image; a fieldmap that
// is missing, unreadable, not 3D, oversized, or off-grid; a non-finite fieldmap voxel; a
// non-finite or non-positive dwell; or an unrecognised unwarpdir.  Returns 0 on success.
int fmap_unwarp(nifti_image *nim, const char *fmapfile, double dwell, const char *unwarpdir);

#ifdef __cplusplus
}
#endif

#endif // FMAP_H

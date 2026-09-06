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
// Every measured convention, the experiment behind it, and the three deliberate divergences are
// recorded in the fmap_bench repository's test/fmap_reference_manifest.md (which currently has no
// git remote -- treat its figures as internal measurements until it is published).  Read that before
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
// Fails closed on: a non-float32 or >4D working image; a fieldmap that is missing, unreadable,
// not 3D, oversized, or off-grid; a non-finite fieldmap voxel; a non-finite or non-positive
// dwell; or an unrecognised unwarpdir.  Returns 0 on success.
//
// Every failure path is all-or-nothing: rejections are detected before any voxel is touched, and
// the resampling loop is entered by the whole thread team or by none of it, so an allocation
// failure there leaves `nim` exactly as it arrived.  (That team-wide decision is load-bearing --
// branching on a per-thread allocation result would put the team on different worksharing
// regions, which hangs.  See the barrier in fmap.c.)
//
// Byte-identical across thread counts, and byte-stable run to run WITHIN a build.  NOT bit-stable
// ACROSS builds: the interior-gap fill is FMA-contracted under the tree's -ffast-math, so a
// different compiler can move roughly 1 voxel in 2400 by ~1 ULP of the shift.
int fmap_unwarp(nifti_image *nim, const char *fmapfile, double dwell, const char *unwarpdir);

#ifdef HAVE_ROMEO
// -fmapprep <brain_magnitude> <deltaTE_ms>: ordinary chain operation, DT32 only.
//
// Turns a two-echo phase-difference image into a B0 fieldmap in rad/s, the units `-fugue` and
// FSL's fugue consume.  The working image is the wrapped phase difference; `magfile` is the
// BRAIN-EXTRACTED magnitude belonging to it (bet output, or any image whose nonzero voxels are
// the region to keep), which supplies both the mask and the anatomical weighting for unwrapping.
// `delta_te_ms` is the echo time difference in MILLISECONDS.
//
//   mask      = (magnitude != 0)                       used verbatim; no erosion, no dilation
//   phase_rad = the working image rescaled so its observed range spans exactly 2*pi
//   unwrapped = ROMEO 3D spatial unwrapping over the mask, magnitude-weighted
//   out(v)    = unwrapped(v)/deltaTE_s - median_mask(...)   inside the mask, 0 outside
//
// There is deliberately NO regularisation: no median filter, no despiking, no smoothing and no
// erosion.  That is measured reference behaviour, not an omission -- a single-voxel spike passes
// through the reference intact.  The ONE post-processing step, 2*pi branch-outlier correction, is
// a divergence forced by using ROMEO rather than PRELUDE and is described in fmap.c.
//
// The median is the UPPER of the two central values for an even population (`sorted[n/2]`), not
// their average.  The distinction is not cosmetic: on a field whose mask splits into two equal
// populations the averaging median is wrong by half the field's range.
//
// Unwrapping is ROMEO, not a reimplementation of the reference's PRELUDE, so the two disagree at
// poorly conditioned voxels.  That is why the acceptance gate is on the final unwarped EPI.
//
// Fails closed on: a non-float32, non-3D or oversized working image; a missing, unreadable,
// non-3D, oversized or off-grid magnitude; an empty mask; a non-finite or constant phase image;
// a non-finite or non-positive delta_te_ms; or an unwrapping failure.  Returns 0 on success.
//
// Unlike fmap_unwarp this is NOT all-or-nothing: the phase is rescaled to radians IN PLACE
// before ROMEO is called, so a failure at or after the unwrap leaves `nim->data` overwritten.
// No output is written either way, which is safe under the op loop's free-without-saving
// contract -- but a caller outside that loop must not reuse `nim` after a failure.
// `debranch` enables the 2*pi branch-outlier correction described in fmap.c (default on; the CLI
// spells the opt-out `-no-debranch`).  It is applied in fmap.c AFTER romeo_unwrap_frame() returns,
// so romeo.c is untouched and neither --medic nor -romeo is affected by it -- MEDIC keeps a
// faithful ROMEO.  The opt-out exists so the raw ROMEO field can be recovered when tracing a
// divergence against a reference MEDIC implementation.
int fmap_prepare(nifti_image *nim, const char *magfile, double delta_te_ms, int debranch);
#endif

#ifdef __cplusplus
}
#endif

#endif // FMAP_H

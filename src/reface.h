/*----------------------------------------------------------------------------
 * reface.h — face-replacement (anonymization) compositing stage: public API.
 *
 * CLEAN-ROOM reimplementation of the AFNI afni_refacer2.csh "reface" recipe.
 * The AFNI script (P.A. Taylor / R.W. Cox, NIH public domain) chains
 * 3dAllineate / 3dcalc / 3dBlurInMask; the arithmetic recipe reproduced here
 * (brightness match `ifac`, sign-based composite, masked edge blend) is trivial
 * and independently written. The Gaussian edge blend uses niimath's OWN
 * BSD/public-domain float32 smoother (nifti_smooth_gauss_f32, via core32.h) as a
 * normalized convolution — it is NOT a port of AFNI's 3dBlurInMask. Unlike
 * qwarp.c (an attributed AFNI port) this file borrows no AFNI source.
 *
 * Leading-edge divergence (2026-07): standalone-only until back-ported to
 * niimath; add reface.c/.h to the parity set at that point (see AGENTS.md).
 *--------------------------------------------------------------------------*/
#ifndef REFACE_H
#define REFACE_H

#include "nifti_io.h"

/* Zero isolated nonzero voxels (AFNI 3dcalc -isola): any voxel whose full
   26-neighborhood is zero is itself set to zero. A single pass tested against a
   snapshot of the original nonzero pattern, so removing one voxel never affects
   another voxel's isolation test. `v` is the float voxel data of an nx*ny*nz
   volume, modified in place. Best-effort: a scratch-allocation failure leaves
   `v` unchanged (a denoiser, never a correctness requirement). */
void reface_isola(float *v, int nx, int ny, int nz);

/* Composite an anonymized ("refaced") image from the original subject and a
   shell that has ALREADY been back-projected onto the subject grid (same dims).
   Both must be plain DT_FLOAT32 single-volume images; neither is modified.
   Shell voxel semantics: >0 replace (scaled), ==0 preserve subject, <0 zero.

   Recipe (AFNI reface): ifac = mean(subject over shell>0 && subject!=0) /
   mean(shell over shell>0); out = shell>0 ? shell*ifac : shell<0 ? 0 :
   max(subject,0); then a normalized Gaussian blur (FWHM 2.666 mm) INSIDE the
   shell>0 mask, preserving voxels outside it.

   On success writes a newly allocated finite float32 image on the subject grid
   to *result (caller frees with nifti_image_free) and returns 0. On any error
   (bad/mismatched/non-float inputs, empty overlap, non-finite ifac, allocation
   or blur failure, non-finite output) returns nonzero and sets *result = NULL —
   no partial result is produced. */
int reface_apply(const nifti_image *subject, const nifti_image *shell_on_subject,
                 nifti_image **result);

#endif /* REFACE_H */

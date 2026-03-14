#ifndef ALLINEATE_H
#define ALLINEATE_H

/* Affine image registration using local Pearson correlation
   Adapted from AFNI's 3dAllineate by RW Cox (public domain)
   Uses lpc+ZZ cost function with twopass coarse-to-fine optimization */

#include "nifti_io.h"

/* Register source image to base image grid using affine (12 DOF) alignment.
   source: the moving image (will be modified in-place: data replaced, dims updated)
   base: the stationary/reference image
   Returns 0 on success, nonzero on error. */
int nii_allineate(nifti_image *source, nifti_image *base);

/* Deface: register template to input, warp mask to input space, zero masked voxels.
   input: the image to deface (modified in-place, stays in its own space)
   tmpl: template image (moving image for registration)
   mask: mask in template space (non-zero = keep, zero = remove face)
   cost_mode: 0 = lpa+ZZ (cross-modal), 1 = lpc+ZZ (structural-to-EPI), 2 = Hellinger (fast)
   Returns 0 on success, nonzero on error. */
int nii_deface(nifti_image *input, nifti_image *tmpl, nifti_image *mask, int cost_mode);

#endif /* ALLINEATE_H */

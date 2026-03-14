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

#endif /* ALLINEATE_H */

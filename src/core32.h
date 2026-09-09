#ifndef CORE32_H
#define CORE32_H

#include "nifti_io.h"

#ifdef __cplusplus
extern "C" {
#endif

int main32(int argc, char *argv[]);

#if defined(HAVE_QC) && defined(HAVE_ALLINEATE)
/* --qc --air helper: RAS T1, head mask, air distance field and the RAS-voxel -> template-mm
   affine, from a float32 3D image. Defined in coreFLT.c (DT32); see the comment there. */
int nii_qc_air_masks_f32(const nifti_image *t1, const char *template, const char *t1name,
                         nifti_image **ras, nifti_image **head, nifti_image **dist, mat44 *vox2tmpl);
#endif

int nifti_smooth_gauss_f32(float *data, int nx, int ny, int nz, int nvol,
							 float dx, float dy, float dz,
							 float sigma_x_mm, float sigma_y_mm, float sigma_z_mm,
							 float kernel_width);

#ifdef __cplusplus
}
#endif

#endif // CORE32_H

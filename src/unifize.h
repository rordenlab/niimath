#ifndef UNIFIZE_H
#define UNIFIZE_H

/* Bias field correction (intensity uniformization)
   Adapted from AFNI's 3dUnifize by RW Cox (public domain) */

int unifize_image(float *data, int nx, int ny, int nz, float dx, float dy, float dz);

#endif /* UNIFIZE_H */

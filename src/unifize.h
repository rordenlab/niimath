#ifndef UNIFIZE_H
#define UNIFIZE_H

/* Bias field correction (intensity uniformization)
   Adapted from AFNI's 3dUnifize by RW Cox (public domain) */

/* do_gm != 0 additionally applies AFNI's -GM global gray-matter scaling. */
int unifize_image(float *data, int nx, int ny, int nz, float dx, float dy, float dz, int do_gm);

#endif /* UNIFIZE_H */

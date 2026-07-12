#ifndef CORE32_H
#define CORE32_H

#ifdef __cplusplus
extern "C" {
#endif

int main32(int argc, char *argv[]);

int nifti_smooth_gauss_f32(float *data, int nx, int ny, int nz, int nvol,
							 float dx, float dy, float dz,
							 float sigma_x_mm, float sigma_y_mm, float sigma_z_mm,
							 float kernel_width);

#ifdef __cplusplus
}
#endif

#endif // CORE32_H

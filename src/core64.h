#ifndef CR_CORE64_H
#define CR_CORE64_H

#ifdef  __cplusplus
extern "C" {
#endif

int main64(int argc, char *argv[]);

int nifti_smooth_gauss_f64(double *data, int nx, int ny, int nz, int nvol,
							 double dx, double dy, double dz,
							 double sigma_x_mm, double sigma_y_mm, double sigma_z_mm,
							 double kernel_width);

#ifdef  __cplusplus
}
#endif

#endif // CR_CORE64_H

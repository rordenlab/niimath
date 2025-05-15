#ifndef EX_CONFORM_H
#define EX_CONFORM_H

#ifdef  __cplusplus
//extern "C" {
#endif

#include <stdbool.h>
#include <nifti2_io.h>

int toRAS(nifti_image *nim);
int conform(nifti_image *nim);
int reslice(nifti_image *nim, nifti_image *nim2, int isLinear);
int comply(nifti_image *nim, const int outDims[3], const float outPixDims[3], float f_high, int isLinear);

#ifdef  __cplusplus
//}
#endif

#endif // EX_CONFORM_H

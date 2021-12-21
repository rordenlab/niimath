#ifndef DG_TENSOR_DECOMP_H
#define DG_TENSOR_DECOMP_H

#ifdef  __cplusplus
extern "C" {
#endif

void EIG_tsfunc(
                          int npts, float ts[],
                          float * val, int isUpperTriangle );

#ifdef  __cplusplus
}
#endif

#endif //DG_TENSOR_DECOMP_H





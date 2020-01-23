#ifndef DG_TENSOR_DECOMP_H
#define DG_TENSOR_DECOMP_H

#ifdef  __cplusplus
extern "C" {
#endif

void EIG_tsfunc( double tzero, double tdelta ,
                          int npts, float ts[],
                          double ts_mean, double ts_slope,
                          void * ud, int nbriks, float * val, int isUpperTriangle );

#ifdef  __cplusplus
}
#endif

#endif //DG_TENSOR_DECOMP_H





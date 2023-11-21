#ifndef EX_FDR_H
#define EX_FDR_H

#ifdef  __cplusplus
//extern "C" {
#endif

#include <stdbool.h>

float fdr(float *pvals, double qval, size_t nvox);

#ifdef  __cplusplus
//}
#endif

#endif // EX_FDR_H

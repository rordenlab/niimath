#ifndef EX_CONFORM_H
#define EX_CONFORM_H

#ifdef  __cplusplus
//extern "C" {
#endif

#include <stdbool.h>
#include <nifti2_io.h>

int conform(nifti_image *nim);

#ifdef  __cplusplus
//}
#endif

#endif // EX_CONFORM_H

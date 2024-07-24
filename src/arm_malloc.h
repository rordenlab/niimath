#ifndef ARM_MALLOC_H
#define ARM_MALLOC_H

#ifdef  __cplusplus
extern "C" {
#endif

#ifdef __aarch64__
 #include <stdio.h>
 #include <stdlib.h>
 #include <stdint.h>

  #define _mm_malloc(size, alignment) malloc(size)
  #define _mm_free(ptr) free(ptr)
#endif

#ifdef  __cplusplus
}
#endif

#endif // ARM_MALLOC_H

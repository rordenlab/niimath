#ifndef AL_SIZE_GUARD_H
#define AL_SIZE_GUARD_H

#include <stdint.h>

/*
 * Voxel-count guard for buffers whose element type is float.
 *
 * INT_MAX is not a sufficient allocation bound on 32-bit targets:
 * counts above SIZE_MAX/sizeof(float) wrap a `count * sizeof(float)` malloc.
 * AL_SIZE_MAX is overrideable only so the native C-API harness can exercise
 * the wasm32/native-32 boundary without allocating a giant image.
 */
#ifndef AL_SIZE_MAX
#define AL_SIZE_MAX SIZE_MAX
#endif

static inline int al_float_nvox_fits(uint64_t nvox) {
    return nvox <= (uint64_t)(AL_SIZE_MAX / sizeof(float));
}

#endif /* AL_SIZE_GUARD_H */

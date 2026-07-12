# zlib-ng (https://github.com/zlib-ng/zlib-ng) — a modern, actively maintained zlib
# replacement with per-architecture optimized deflate/inflate/CRC (NEON + ARMv8 CRC32 on
# arm64; SSE2/SSSE3/PCLMUL/AVX on x86_64; correct on MSVC). Built in ZLIB_COMPAT mode so it
# installs a drop-in `zlib.h` + static `libz`/`zlibstatic` that find_package(ZLIB) resolves
# exactly like stock zlib — no niimath source change needed.
#
# Why it is the default over the Cloudflare fork: Cloudflare's `gcc.amd64` branch is x86/GCC
# inline-asm tuned, has no arm64 acceleration (it falls back to generic on Apple Silicon), and
# emits a CORRUPT deflate stream under MSVC (hence External-CLOUDFLARE-ZLIB.cmake substitutes
# madler zlib on Windows). zlib-ng is fast on x86_64 AND arm64 AND correct on MSVC, so it gives
# one consistent, fast zlib across every release target. Measured ~2x stock system zlib and
# ~1.3x the Cloudflare fork on arm64 for niimath's `.nii.gz` write path.

set(ZLIBNG_TAG 2.3.3) # pinned stable release

set(ZLIB_CMAKE_ARGS
    -Wno-dev
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX:PATH=${DEP_INSTALL_DIR}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DZLIB_COMPAT:BOOL=ON           # drop-in `zlib.h`/`gzopen`/`deflate` API + libz naming
    -DZLIB_ENABLE_TESTS:BOOL=OFF
    -DWITH_GTEST:BOOL=OFF
)

# Static MSVC runtime so niimath.exe stays self-contained (same rationale as the Cloudflare
# external: without CMP0091=NEW + MultiThreaded, a /MD static lib's __imp_* CRT imports fail to
# link into the /MT niimath.exe).
if (MSVC)
    list(APPEND ZLIB_CMAKE_ARGS
        -DCMAKE_POLICY_DEFAULT_CMP0091:STRING=NEW
        -DCMAKE_MSVC_RUNTIME_LIBRARY:STRING=MultiThreaded)
endif()

# Only add OS X architectures if explicitly defined (the universal build drives x86_64 and
# arm64 as separate single-arch SuperBuild passes and lipo-combines them).
if(CMAKE_OSX_ARCHITECTURES)
    list(APPEND ZLIB_CMAKE_ARGS -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES})
endif()

ExternalProject_Add(zlib
    GIT_REPOSITORY "${git_protocol}github.com/zlib-ng/zlib-ng.git"
    GIT_TAG "${ZLIBNG_TAG}"
    SOURCE_DIR zlib-ng
    BINARY_DIR zlib-ng-build
    CMAKE_ARGS ${ZLIB_CMAKE_ARGS}
)

set(ZLIB_ROOT ${DEP_INSTALL_DIR})

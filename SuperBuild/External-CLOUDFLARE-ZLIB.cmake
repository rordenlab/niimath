set(CLOUDFLARE_BRANCH gcc.amd64) # Cloudflare zlib branch

# Always disable shared libraries and enforce static MSVC runtime on Windows
set(ZLIB_CMAKE_ARGS
    -Wno-dev
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX:PATH=${DEP_INSTALL_DIR}
    -DBUILD_SHARED_LIBS:BOOL=OFF
)

# Only set MSVC static runtime if on Windows
if (MSVC)
    list(APPEND ZLIB_CMAKE_ARGS -DCMAKE_C_FLAGS_RELEASE:STRING=-MT)
endif()

# Only add OS X architectures if explicitly defined
if(CMAKE_OSX_ARCHITECTURES)
    list(APPEND ZLIB_CMAKE_ARGS -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES})
endif()

ExternalProject_Add(zlib
    GIT_REPOSITORY "${git_protocol}github.com/ningfei/zlib.git"
    GIT_TAG "${CLOUDFLARE_BRANCH}"
    SOURCE_DIR cloudflare-zlib
    BINARY_DIR cloudflare-zlib-build
    CMAKE_ARGS ${ZLIB_CMAKE_ARGS}
)

set(ZLIB_ROOT ${DEP_INSTALL_DIR})

set(CLOUDFLARE_BRANCH gcc.amd64) # Cloudflare zlib branch

# Always disable shared libraries and enforce static MSVC runtime on Windows
set(ZLIB_CMAKE_ARGS
    -Wno-dev
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX:PATH=${DEP_INSTALL_DIR}
    -DBUILD_SHARED_LIBS:BOOL=OFF
)

# Only set MSVC static runtime if on Windows. Use the CMAKE_MSVC_RUNTIME_LIBRARY
# abstraction (with CMP0091=NEW) rather than a raw -MT flag: newer CMake sets the
# runtime via that variable and its default (MultiThreadedDLL = /MD) otherwise
# overrides a bare -MT, producing a /MD static lib whose __imp_* CRT imports fail
# to link into the /MT niimath.exe (LNK2019). Static CRT here keeps niimath.exe
# self-contained.
if (MSVC)
    list(APPEND ZLIB_CMAKE_ARGS
        -DCMAKE_POLICY_DEFAULT_CMP0091:STRING=NEW
        -DCMAKE_MSVC_RUNTIME_LIBRARY:STRING=MultiThreaded)
endif()

# Only add OS X architectures if explicitly defined
if(CMAKE_OSX_ARCHITECTURES)
    list(APPEND ZLIB_CMAKE_ARGS -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES})
endif()

# The Cloudflare fork's `gcc.amd64` branch is optimized with GCC/AMD64 inline
# assembly. Under MSVC it compiles but its DEFLATE encoder emits a corrupt
# stream for non-trivial data (a defaced/`-add` volume was unreadable even by
# niimath's own inflate; simple near-binary data like `-edt` happened to
# survive). INFLATE is fine, so reading Cloudflare-written files still works.
# Use upstream madler zlib on MSVC; keep the fast Cloudflare fork on gcc/clang
# (macOS/Linux), where its deflate is correct.
if(MSVC)
    ExternalProject_Add(zlib
        GIT_REPOSITORY "${git_protocol}github.com/madler/zlib.git"
        GIT_TAG "v1.3.1"
        SOURCE_DIR madler-zlib
        BINARY_DIR madler-zlib-build
        CMAKE_ARGS ${ZLIB_CMAKE_ARGS}
    )
else()
    ExternalProject_Add(zlib
        GIT_REPOSITORY "${git_protocol}github.com/ningfei/zlib.git"
        GIT_TAG "${CLOUDFLARE_BRANCH}"
        SOURCE_DIR cloudflare-zlib
        BINARY_DIR cloudflare-zlib-build
        CMAKE_ARGS ${ZLIB_CMAKE_ARGS}
    )
endif()

set(ZLIB_ROOT ${DEP_INSTALL_DIR})

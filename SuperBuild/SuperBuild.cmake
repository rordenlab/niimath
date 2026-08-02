# Check if git exists
find_package(Git)
if(NOT GIT_FOUND)
    message(FATAL_ERROR "Cannot find Git. Git is required for Superbuild")
endif()

# Use git protocol or not
option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use https instead." OFF)
if(USE_GIT_PROTOCOL)
    set(git_protocol "ssh://git@")
else()
    set(git_protocol "https://")
endif()

# emulate fslmaths behavior, add pigz support
option(FSLSTYLE "FSL behavior, pigz support" ON)
if(FSLSTYLE)
   ADD_DEFINITIONS(-DFSLSTYLE)
   ADD_DEFINITIONS(-DPIGZ)
   ADD_DEFINITIONS(-DREJECT_COMPLEX)
endif()

if(NOT BUILD_FLAVOR)
    set(BUILD_FLAVOR "all" CACHE STRING
        "Choose the flavor of build, options are: all tiny nano." FORCE)
    set_property(CACHE BUILD_FLAVOR PROPERTY STRINGS "all;tiny;nano")
endif()

# Basic CMake build settings
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(USE_STATIC_RUNTIME "Use static runtime" ON)

# Ensure static MSVC runtime and static libraries
set(BUILD_SHARED_LIBS OFF)
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded")
endif()

option(ENABLE_ZSTD "Enable zstd (.nii.zst) compression support" ON)
option(BUILD_BMP "Build bitmap (PNG) output support" ON)
option(ENABLE_GPL "Enable optional copyleft module (-spm_coreg/-spm_deface); GPL-2-or-later binary" OFF)
# These mirror src/CMakeLists.txt options; declare and forward them here so the
# documented top-level build (cmake .. on the repo root) actually honours them
# instead of silently reporting an unused variable.
option(USE_OPENMP "Build with OpenMP support" ON)
option(ENABLE_QC "Enable anatomical QC metrics (--qc)" ON)
option(ENABLE_MEDIC "Enable MEDIC multi-echo distortion correction (--medic, -unwarp)" ON)
# -moco and -stc are ON for every target; their former Apple-Silicon-only gate was lifted once
# release_smoke.py (Windows/Ubuntu/macOS via AppVeyor) and js/tests/temporal.test.ts (Emscripten)
# began checking both numerically on every build. These are forwarded to the inner project as
# explicit -D cache entries, which its own option() default can never override, so any future
# per-target rule must be applied HERE as well as in src/CMakeLists.txt.
option(ENABLE_MOCO "Enable rigid-body motion correction (-moco)" ON)
option(ENABLE_STC "Enable slice-time correction (-stc)" ON)
# Forwarded as an explicit -D below: src/ own option() default can never override a
# SuperBuild cache entry, so a platform rule must be mirrored in BOTH files.
option(ENABLE_SKULLSTRIP "Enable AFNI-style surface skull stripping (-skullstrip)" OFF)
# Mirrored from src/CMakeLists.txt, per the rule in AGENTS.md known-issue 2: this value is
# forwarded as an explicit -D cache entry that src/'s own option() default can never override,
# so a platform rule enforced only there would let the top level configure cleanly and then
# hard-fail inside the inner project. Reject the unsupported request here too.
if(ENABLE_SKULLSTRIP AND (CMAKE_SIZEOF_VOID_P EQUAL 4 OR NOT BUILD_FLAVOR STREQUAL "all"))
  message(FATAL_ERROR "ENABLE_SKULLSTRIP=ON requires a 64-bit native BUILD_FLAVOR=all build.")
endif()
# Mirror the Emscripten guard too, not only the pointer-size/flavor one: wasm64 reports 8-byte
# pointers and would pass the test above, so an outer configure could succeed and the inner
# project then hard-fail.
if(ENABLE_SKULLSTRIP AND (EMSCRIPTEN OR CMAKE_SYSTEM_NAME STREQUAL "Emscripten"))
  message(FATAL_ERROR "ENABLE_SKULLSTRIP=ON is not supported for WebAssembly targets.")
endif()
option(ENABLE_ROMEO "Enable ROMEO phase unwrapping (-romeo)" ON)
option(ENABLE_ALLINEATE "Enable allineate affine registration" ON)
option(ENABLE_QWARP "Enable -qwarp nonlinear (deformable) registration" OFF)
if(ENABLE_QWARP AND NOT ENABLE_ALLINEATE)
  message(FATAL_ERROR "ENABLE_QWARP=ON requires ENABLE_ALLINEATE=ON (qwarp uses the allineate/NEWUOA optimizer).")
endif()
# OPENMP_XCODE is the AppleClang-only legacy alias (src/CMakeLists.txt gates OpenMP
# on `USE_OPENMP AND OPENMP_XCODE`). Its DEFAULT follows USE_OPENMP so a plain
# top-level build enables OpenMP on Apple (matching the docs), but it stays an
# explicit, consumed override: the macOS universal-release scripts pass
# `-DOPENMP_XCODE=OFF` to disable OpenMP for the cross-arch slices, where the
# single-arch Homebrew libomp cannot link into a universal binary. Declared after
# USE_OPENMP so its default can reference it; forwarded as ${OPENMP_XCODE} below.
option(OPENMP_XCODE "AppleClang OpenMP (defaults to USE_OPENMP; set OFF for universal builds)" ${USE_OPENMP})

include(ExternalProject)

set(DEPENDENCIES)

set(DEP_INSTALL_DIR ${CMAKE_BINARY_DIR})
# zlib-ng is the release baseline: fast on x86_64 AND arm64, and correct on MSVC (unlike the
# Cloudflare fork, which has no arm64 acceleration and needs a madler fallback on Windows).
set(ZLIB_IMPLEMENTATION "zlib-ng" CACHE STRING "Choose zlib implementation.")
set_property(CACHE ZLIB_IMPLEMENTATION PROPERTY STRINGS  "zlib-ng;Cloudflare;System;Custom")
if(${ZLIB_IMPLEMENTATION} STREQUAL "zlib-ng")
    message("-- Build with zlib-ng: ON")
    include(${CMAKE_SOURCE_DIR}/SuperBuild/External-ZLIBNG.cmake)
    list(APPEND DEPENDENCIES zlib)
    set(BUILD_ZLIBNG TRUE)
    message("--   Will build zlib-ng (ZLIB_COMPAT) from github")
elseif(${ZLIB_IMPLEMENTATION} STREQUAL "Cloudflare")
    message("-- Build with Cloudflare zlib: ON")
    include(${CMAKE_SOURCE_DIR}/SuperBuild/External-CLOUDFLARE-ZLIB.cmake)
    list(APPEND DEPENDENCIES zlib)
    set(BUILD_CLOUDFLARE-ZLIB TRUE)
    message("--   Will build Cloudflare zlib from github")
elseif(${ZLIB_IMPLEMENTATION} STREQUAL "Custom")
    set(ZLIB_ROOT ${ZLIB_ROOT} CACHE PATH "Specify custom zlib root directory.")
    if(NOT ZLIB_ROOT)
        message(FATAL_ERROR "ZLIB_ROOT needs to be set to locate custom zlib!")
    endif()
endif()

ExternalProject_Add(src
    DEPENDS ${DEPENDENCIES}
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/src
    BINARY_DIR src-build
    CMAKE_ARGS
        -Wno-dev
        --no-warn-unused-cli
        -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_VERBOSE_MAKEFILE:BOOL=${CMAKE_VERBOSE_MAKEFILE}
        -DBUILD_FLAVOR:STRING=${BUILD_FLAVOR}
        -DOPENMP_XCODE:BOOL=${OPENMP_XCODE}
        -DUSE_STATIC_RUNTIME:BOOL=${USE_STATIC_RUNTIME}
        -DZLIB_IMPLEMENTATION:STRING=${ZLIB_IMPLEMENTATION}
        -DZLIB_ROOT:PATH=${ZLIB_ROOT}
        -DENABLE_ZSTD:BOOL=${ENABLE_ZSTD}
        -DZSTD_ROOT:PATH=${ZSTD_ROOT}
        -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}
        -DENABLE_GPL:BOOL=${ENABLE_GPL}
        -DUSE_OPENMP:BOOL=${USE_OPENMP}
        -DENABLE_QC:BOOL=${ENABLE_QC}
        -DENABLE_MEDIC:BOOL=${ENABLE_MEDIC}
        -DENABLE_MOCO:BOOL=${ENABLE_MOCO}
        -DENABLE_STC:BOOL=${ENABLE_STC}
        -DENABLE_SKULLSTRIP:BOOL=${ENABLE_SKULLSTRIP}
        -DENABLE_ROMEO:BOOL=${ENABLE_ROMEO}
        -DCMAKE_INTERPROCEDURAL_OPTIMIZATION:BOOL=${CMAKE_INTERPROCEDURAL_OPTIMIZATION}
        -DENABLE_ALLINEATE:BOOL=${ENABLE_ALLINEATE}
        -DENABLE_QWARP:BOOL=${ENABLE_QWARP}
        -DBUILD_BMP:BOOL=${BUILD_BMP}
        # forward static runtime and static linking
        -DCMAKE_MSVC_RUNTIME_LIBRARY:STRING=${CMAKE_MSVC_RUNTIME_LIBRARY}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/bin/ DESTINATION ${SKBUILD_PROJECT_NAME}/bin
        USE_SOURCE_PERMISSIONS)


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

if(${CMAKE_C_COMPILER_ID} STREQUAL "AppleClang")
    option(OPENMP_XCODE "Build with OpenMP support" OFF)
endif()

include(ExternalProject)

set(DEPENDENCIES)

set(DEP_INSTALL_DIR ${CMAKE_BINARY_DIR})
set(ZLIB_IMPLEMENTATION "Cloudflare" CACHE STRING "Choose zlib implementation.")
set_property(CACHE ZLIB_IMPLEMENTATION PROPERTY STRINGS  "Cloudflare;Miniz;System;Custom")
if(${ZLIB_IMPLEMENTATION} STREQUAL "Cloudflare")
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
        # forward static runtime and static linking
        -DCMAKE_MSVC_RUNTIME_LIBRARY:STRING=${CMAKE_MSVC_RUNTIME_LIBRARY}
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/bin/ DESTINATION ${SKBUILD_PROJECT_NAME}/bin
        USE_SOURCE_PERMISSIONS)

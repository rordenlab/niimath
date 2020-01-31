cmake_minimum_required(VERSION 3.1.0)

project(console C)

# Option Choose whether to use static runtime
include(ucm.cmake)
option(USE_STATIC_RUNTIME "Use static runtime" ON)
if(USE_STATIC_RUNTIME)
    ucm_set_runtime(STATIC)
else()
    ucm_set_runtime(DYNAMIC)
endif()

# Basic CMake build settings
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS  "Debug;Release;RelWithDebInfo;MinSizeRel")
endif()


if(${CMAKE_C_COMPILER_ID} STREQUAL "AppleClang")
    # using AppleClang
    add_definitions(-fno-caret-diagnostics)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-dead_strip")
elseif(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
    # using GCC
    set(CMAKE_C_STANDARD 11)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    endif()
elseif(MSVC)
    # using Visual Studio C++
    add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:8388608")  # set "Stack Reserve Size" to 8MB (default value is 1MB)
endif()

# Compiler dependent flags
include (CheckCCompilerFlag)
if(UNIX)
    check_c_compiler_flag(-msse2 HAS_SSE2)
    if(HAS_SSE2)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2 -mfpmath=sse")
    endif()
endif()

set(PROGRAMS niimath)

set(NIIMATH_SRCS
    niimath.c
    core.c
    tensor.c
    core32.c
    core64.c
    niftilib/nifti2_io.c
    znzlib/znzlib.c
    )

include_directories(niftilib)
include_directories(znzlib)
add_executable(niimath ${NIIMATH_SRCS})

set(ZLIB_IMPLEMENTATION "Miniz" CACHE STRING "Choose zlib implementation.")
set_property(CACHE ZLIB_IMPLEMENTATION PROPERTY STRINGS  "Miniz;System;Custom")
if(NOT ${ZLIB_IMPLEMENTATION} STREQUAL "Miniz")
    if(NOT ${ZLIB_IMPLEMENTATION} STREQUAL "System")
        set(ZLIB_ROOT ${ZLIB_ROOT} CACHE PATH "Specify custom zlib root directory.")
        if(NOT ZLIB_ROOT)
            message(FATAL_ERROR "ZLIB_ROOT needs to be set to locate custom zlib!")
        endif()
    endif()
    find_package(ZLIB REQUIRED)
    add_definitions(-DHAVE_ZLIB)
    target_include_directories(niimath PRIVATE ${ZLIB_INCLUDE_DIRS})
    target_link_libraries(niimath ${ZLIB_LIBRARIES})
endif()

#link math library
target_link_libraries(niimath m)

install(TARGETS ${PROGRAMS} DESTINATION bin)

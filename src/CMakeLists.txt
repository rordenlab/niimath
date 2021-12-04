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

#emulate fslmaths behavior, add pigz support
option(FSLSTYLE "FSL behavior, pigz support" ON)
if(FSLSTYLE)
   ADD_DEFINITIONS(-DFSLSTYLE)
   ADD_DEFINITIONS(-DREJECT_COMPLEX)
   ADD_DEFINITIONS(-DPIGZ)
endif()

#e.g.  cmake -DOPENMP_XCODE=ON ..
option(OPENMP_XCODE "Build the static library" OFF)
if(${CMAKE_C_COMPILER_ID} STREQUAL "AppleClang")
    # using AppleClang
    add_definitions(-fno-caret-diagnostics)
    if (OPENMP_XCODE)
        set(CMAKE_C_FLAGS "-I/usr/local/include -Xclang -fopenmp ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/lib -lomp")
    endif()
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
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4018")   # '<': signed/unsigned mismatch
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4244")   # 'initializing': conversion from 'double' to 'int', possible loss of data
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4267")   # 'initializing': conversion from 'size_t' to 'int', possible loss of data
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /wd4305")   # 'argument': truncation from 'double' to 'float'
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:8388608")  # set "Stack Reserve Size" to 8MB (default value is 1MB)
endif()

# Compiler dependent flags
include (CheckCCompilerFlag)
if(UNIX)
    check_c_compiler_flag(-march=armv8-a+crc ARM_CRC)
    if(ARM_CRC)
        # wrong answer for Apple Silicon: check_cxx_compiler_flag(-msse2 HAS_SSE2)
    else()
	    check_c_compiler_flag(-msse2 HAS_SSE2)
	    if(HAS_SSE2)
	        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2 -mfpmath=sse")
	    endif()
    endif()
endif()

# support crc intrinsics for ARM CPUs
if(UNIX)
    check_c_compiler_flag(-march=armv8-a+crc HAS_CRC)
    if(HAS_CRC)
    	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=armv8-a+crc")
    endif()
endif()

set(PROGRAMS niimath)

set(NIIMATH_SRCS
    niimath.c
    core.c
    tensor.c
    bwlabel.c
    bw.c
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

if(NOT MSVC)
    # Link math library
    target_link_libraries(niimath m)
endif()

if(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
    if(USE_STATIC_RUNTIME)
        execute_process(COMMAND ${CMAKE_C_COMPILER} "--print-file-name=libgomp.a"
                        OUTPUT_VARIABLE OPENMP_LIBRARY OUTPUT_STRIP_TRAILING_WHITESPACE)
        target_link_libraries(niimath ${OPENMP_LIBRARY})
        if(MINGW)
            target_link_libraries(niimath "-Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive")
        else()
            target_link_libraries(niimath ${CMAKE_DL_LIBS})
        endif()
    endif()
endif()

install(TARGETS ${PROGRAMS} DESTINATION bin)

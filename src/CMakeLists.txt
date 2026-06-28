cmake_minimum_required(VERSION 3.5)

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
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")
endif()

#emulate fslmaths behavior, add pigz support
option(FSLSTYLE "FSL behavior, pigz support" ON)
if(FSLSTYLE)
   ADD_DEFINITIONS(-DFSLSTYLE)
   ADD_DEFINITIONS(-DPIGZ)
   ADD_DEFINITIONS(-DREJECT_COMPLEX)
endif()

if(NOT BUILD_FLAVOR)
    set(BUILD_FLAVOR "all" CACHE STRING
        "Choose the flavor of build, options are: all tiny nano." FORCE)
    set_property(CACHE BUILD_FLAVOR PROPERTY STRINGS  "all;tiny;nano")
endif()

if(${BUILD_FLAVOR} STREQUAL "all")
    ADD_DEFINITIONS(-DHAVE_64BITS)
    set(ADDITIONAL_SRCS core64.c)
    ADD_DEFINITIONS(-DHAVE_BUTTERWORTH)
    set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} bw.c)
    ADD_DEFINITIONS(-DHAVE_FORMATS)
    set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} base64.c)
    ADD_DEFINITIONS(-DHAVE_TENSOR)
    set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} tensor.c)
    ADD_DEFINITIONS(-DHAVE_DTIFIT)
    set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} dtifit.c)
    ADD_DEFINITIONS(-DNII2MESH)
    set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} bwlabel.c fdr.c meshify.c quadric.c radixsort.c)
    ADD_DEFINITIONS(-DHAVE_CONFORM)
    set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} conform.c)
    set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} unifize.c)
    option(ENABLE_ALLINEATE "Enable allineate affine registration" ON)
    if(ENABLE_ALLINEATE)
        ADD_DEFINITIONS(-DHAVE_ALLINEATE)
        set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} allineate.c powell_newuoa.c)
        # allineate benefits from fast-math (preserving NaN/Inf detection)
        if(NOT MSVC)
            set_source_files_properties(allineate.c powell_newuoa.c PROPERTIES
                COMPILE_FLAGS "-ffast-math -fno-finite-math-only")
        endif()
    endif()

    # Optional GPL spm_coreg module (-spm_coreg/-spm_deface). OFF so default builds stay
    # BSD-2; ON produces a GPL-2 combined-work binary. Sources from the niimath_gpl submodule.
    option(ENABLE_GPL "Enable optional GPL spm_coreg module (-spm_coreg/-spm_deface); GPL-2 binary" OFF)
    if(ENABLE_GPL)
        if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/GPL/spm_coreg.c")
            message(FATAL_ERROR "ENABLE_GPL set but src/GPL is empty — run: git submodule update --init")
        endif()
        ADD_DEFINITIONS(-DHAVE_GPL)
        include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GPL ${CMAKE_CURRENT_SOURCE_DIR})
        # GPL estimate only (returns the rigid transform); reslice uses BSD code.
        set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS}
            GPL/cost.c GPL/hist2.c GPL/loaduint8.c GPL/matrix.c GPL/powell.c
            GPL/smooth.c GPL/spm_coreg.c GPL/spmcoreg_niimath.c)
    endif()

    option(USE_CLASSIC_CUBES "Use classic or marching cubes" ON)
    if(USE_CLASSIC_CUBES)
        ADD_DEFINITIONS(-DUSE_CLASSIC_CUBES)
        set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} oldcubes.c)
    else()
        set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} MarchingCubes.c)
    endif()
elseif(${BUILD_FLAVOR} STREQUAL "tiny")
    ADD_DEFINITIONS(-DNII2MESH)
    set(ADDITIONAL_SRCS bwlabel.c fdr.c meshify.c oldcubes.c quadric.c radixsort.c)

    option(USE_CLASSIC_CUBES "Use classic or marching cubes" ON)
    if(USE_CLASSIC_CUBES)
        ADD_DEFINITIONS(-DUSE_CLASSIC_CUBES)
        set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} oldcubes.c)
    else()
        set(ADDITIONAL_SRCS ${ADDITIONAL_SRCS} MarchingCubes.c)
    endif()
endif()

if(${CMAKE_C_COMPILER_ID} STREQUAL "AppleClang")
    #e.g.  cmake -DOPENMP_XCODE=ON ..
    # Default ON when allineate is enabled (benefits from parallel coarse search)
    if(ENABLE_ALLINEATE)
        option(OPENMP_XCODE "Build with OpenMP support" ON)
    else()
        option(OPENMP_XCODE "Build with OpenMP support" OFF)
    endif()
    # using AppleClang
    add_definitions(-fno-caret-diagnostics)
    if (OPENMP_XCODE)
        execute_process(COMMAND brew --prefix OUTPUT_VARIABLE HOMEBREW_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT EXISTS "${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a")
            message(FATAL_ERROR "Install Homebrew and run 'brew install libomp' to enable OpenMP")
        else()
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xclang -fopenmp -I${HOMEBREW_PREFIX}/opt/libomp/include")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.a")
            message(STATUS "OpenMP support enabled")
        endif()
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
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=armv8-a+crc")
    else()
	    check_c_compiler_flag(-msse2 HAS_SSE2)
	    if(HAS_SSE2)
	        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2 -mfpmath=sse")
	    endif()
    endif()
endif()

set(PROGRAMS niimath)

# Option to build bitmap support (ON by default)
option(BUILD_BMP "Build bitmap output and include filter.c bmp.c" ON)

if(BUILD_BMP)
    message(STATUS "Bitmap support: ON (adding filter.c, bmp.c, spng.c, and -DHAVE_BMP)")
    add_definitions(-DHAVE_BMP)
    set(BMP_SRCS filter.c bmp.c spng.c)
else()
    message(STATUS "Bitmap support: OFF")
    set(BMP_SRCS "")  # empty so appending is safe
endif()

set(NIIMATH_SRCS
    niimath.c
    core.c
    core32.c
    nifti_io.c
    ${ADDITIONAL_SRCS}
    ${BMP_SRCS}
)
add_executable(niimath ${NIIMATH_SRCS})

# Default to System zlib for standalone `cmake -S src` builds: there is no bundled
# miniz, and spng.c (bitmap support, on by default) needs zlib's deflate/inflate, so
# the old "Miniz" default failed to link. The SuperBuild overrides this with its
# Cloudflare zlib via -DZLIB_IMPLEMENTATION, so release/CI/wheel builds are unaffected.
set(ZLIB_IMPLEMENTATION "System" CACHE STRING "Choose zlib implementation.")
set_property(CACHE ZLIB_IMPLEMENTATION PROPERTY STRINGS  "System;Custom")
# There is no bundled miniz; spng.c (bitmap) and nifti gz I/O need a real zlib, so
# reject the legacy "Miniz" value instead of silently producing an unlinkable build.
# (The SuperBuild forwards "Cloudflare" with a ZLIB_ROOT, which the Custom path below
# resolves.) zlib is always linked.
if(${ZLIB_IMPLEMENTATION} STREQUAL "Miniz")
    message(FATAL_ERROR "ZLIB_IMPLEMENTATION=Miniz is not supported (no bundled miniz). Use System (default) or Custom with ZLIB_ROOT.")
endif()
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

if(NOT MSVC)
    # Link math library
    target_link_libraries(niimath m)
endif()

# Zstd compression support (optional, enabled by default if found)
option(ENABLE_ZSTD "Enable zstd (.nii.zst) compression support" ON)
if(ENABLE_ZSTD)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(ZSTD QUIET libzstd)
    endif()
    if(NOT ZSTD_FOUND)
        find_library(ZSTD_LIBRARIES NAMES zstd)
        find_path(ZSTD_INCLUDE_DIRS NAMES zstd.h)
        if(ZSTD_LIBRARIES AND ZSTD_INCLUDE_DIRS)
            set(ZSTD_FOUND TRUE)
        endif()
    endif()
    if(ZSTD_FOUND)
        # pkg-config returns a bare name ("zstd") plus a separate library dir. Resolve a
        # full path so the correct library links from non-default prefixes (e.g.
        # /opt/homebrew, which is not on the default link path) and so cross-built or
        # universal release libraries are picked instead of a host-arch homebrew lib.
        set(ZSTD_LINK ${ZSTD_LIBRARIES})
        if(ZSTD_LIBRARY_DIRS)
            find_library(ZSTD_LIB_FULL NAMES zstd libzstd PATHS ${ZSTD_LIBRARY_DIRS} NO_DEFAULT_PATH)
            if(ZSTD_LIB_FULL)
                set(ZSTD_LINK ${ZSTD_LIB_FULL})
            endif()
        endif()
        add_definitions(-DHAVE_ZSTD)
        target_include_directories(niimath PRIVATE ${ZSTD_INCLUDE_DIRS})
        target_link_libraries(niimath ${ZSTD_LINK})
        message(STATUS "Zstd support: ON (${ZSTD_LINK})")
    else()
        message(STATUS "Zstd support: OFF (libzstd not found, install with: brew install zstd)")
    endif()
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

# For Python package, we need to install to the package directory
if(SKBUILD)
    # scikit-build sets SKBUILD environment variable
    # Install directly in the niimath package directory
    install(TARGETS ${PROGRAMS} DESTINATION niimath)
else()
    install(TARGETS ${PROGRAMS} DESTINATION bin)
endif()

cmake_minimum_required(VERSION 3.5)

if(COMMAND CMAKE_POLICY)
    CMAKE_POLICY(SET CMP0003 NEW)
endif()

if(NOT SKBUILD_PROJECT_NAME)
  set(SKBUILD_PROJECT_NAME niimath)
endif()
project(${SKBUILD_PROJECT_NAME})

# These must be set BEFORE any add_subdirectory or include() that compiles code
set(BUILD_SHARED_LIBS OFF)
if (MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded") # static MSVC CRT
endif()

# Optionally enable Link Time Optimization (LTO) if supported
option(USE_LTO "Enable Link Time Optimization if supported" ON)
if (USE_LTO)
  if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.9")
    if (POLICY CMP0069)
      cmake_policy(SET CMP0069 NEW)
    endif()
    include(CheckIPOSupported)
    check_ipo_supported(RESULT lto_supported OUTPUT lto_output)
    if (lto_supported)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
      message(STATUS "Enabling Link Time Optimization (LTO)")
    endif()
  endif()
endif()

include(${CMAKE_SOURCE_DIR}/SuperBuild/SuperBuild.cmake)

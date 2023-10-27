# Minimum cmake version:
# - 3.12 is needed for https://cmake.org/cmake/help/latest/command/target_link_libraries.html#id8
# - 3.12 is needed for https://cmake.org/cmake/help/latest/command/add_compile_definitions.html
cmake_minimum_required(VERSION 3.12)

project(fabm
  VERSION 2.0.0
  DESCRIPTION "Framework for Aquatic Biogeochemical Models"
  HOMEPAGE_URL https://fabm.net
  LANGUAGES Fortran
 )

if(FABM_FORCED_HOST)
  # FABM_HOST was specified by a higher level CMakeLists.txt
  unset(FABM_HOST CACHE)
  set(FABM_HOST ${FABM_FORCED_HOST})
else()
  # FABM_HOST is user-configurable.
  # Make sure FABM_HOST is set to a subdir of src/drivers
  file(GLOB HOSTNAMES RELATIVE "${CMAKE_CURRENT_LIST_DIR}/src/drivers" "${CMAKE_CURRENT_LIST_DIR}/src/drivers/*")
  set(FABM_HOST "gotm" CACHE STRING "Host that FABM should be compiled for")
  set_property(CACHE FABM_HOST PROPERTY STRINGS ${HOSTNAMES})
  message(STATUS "FABM host: ${FABM_HOST}")
endif()

option(FABM_EMBED_VERSION "Embed FABM version information" OFF)

# Use host-prescribed real kind if set.
if(FABM_FORCED_REAL_KIND)
  # FABM_REAL_KIND was specified by a higher level CMakeLists.txt
  unset(FABM_REAL_KIND CACHE)
  set(FABM_REAL_KIND ${FABM_FORCED_REAL_KIND})
else()
  # FABM_REAL_KIND is user-configurable [as advanced variable]
  set(FABM_REAL_KIND "selected_real_kind(13)" CACHE STRING "Fortran real kind for host API")
  mark_as_advanced(FABM_REAL_KIND)
endif()
set(FABM_REAL_KIND_EXTERNAL ${FABM_REAL_KIND})
set(FABM_REAL_KIND_INTERNAL ${FABM_REAL_KIND_EXTERNAL} CACHE STRING "Fortran real kind for biogeochemistry.")
mark_as_advanced(FABM_REAL_KIND_INTERNAL)

option(FABM_USE_IEEE_ARITHMETIC "Use IEEE_ARITHMETIC module" ON)
option(FABM_USE_DO_CONCURRENT "Use DO CONCURRENT statement (Fortran 2008)" ON)
option(FABM_USE_CONTIGUOUS "Use CONTIGUOUS attribute (Fortran 2008)" ON)
mark_as_advanced(FABM_USE_IEEE_ARITHMETIC FABM_USE_DO_CONCURRENT FABM_USE_CONTIGUOUS)

option(FABM_BUILD_C_API "Build C interfaces" OFF)
mark_as_advanced(FABM_BUILD_C_API)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  add_compile_options(-ffree-line-length-none)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(WIN32)
    add_compile_options("$<$<CONFIG:DEBUG>:/Od>")
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  add_compile_options(-N 1023)
endif()

# Customize preprocessor definitions
add_compile_definitions($<$<NOT:$<CONFIG:DEBUG>>:NDEBUG>)

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS OR FABM_BUILD_C_API)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Set default installation prefix.
if(WIN32)
  if(DEFINED ENV{LOCALAPPDATA})
    set(DEFAULT_PREFIX "$ENV{LOCALAPPDATA}/fabm/${FABM_HOST}")
  else()
    set(DEFAULT_PREFIX "$ENV{APPDATA}/fabm/${FABM_HOST}")
  endif()
else()
  set(DEFAULT_PREFIX "$ENV{HOME}/local/fabm/${FABM_HOST}")
endif()

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Use "include" rather than "add_subdirectory" to avoid inserting a separate
# "src" level in the current binary directory
include("${CMAKE_CURRENT_LIST_DIR}/src/CMakeLists.txt")


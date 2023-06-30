cmake_minimum_required(VERSION 3.0)

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create main project.
project(python_fabm Fortran)

# Specify position-independent code since we will build a shared library.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(WIN32)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
  else()
    # Do not warn about Windows-specific export directives
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable 7841")
  endif()
endif()

if(MSVC)
  add_compile_options("/libs:static")
endif()

set(FABM_BASE ${PROJECT_SOURCE_DIR}/../../.. CACHE PATH "Directory with FABM source code")

set(PYFABM_NAME "fabm_c" CACHE STRING "Base name of library")
set(PYFABM_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "Directory to create library in")
set(PYFABM_DIM_COUNT 0 CACHE STRING "Number of spatial dimensions")
add_definitions(-D_FABM_DIMENSION_COUNT_=${PYFABM_DIM_COUNT})

# Force FABM_HOST to python and include FABM CMakeLists.txt.
set(FABM_FORCED_HOST python)
set(FABM_BUILD_C_API ON CACHE BOOL "Build C interfaces" FORCE)
add_subdirectory(${FABM_BASE} fabm)

set_target_properties(fabm_c PROPERTIES
  OUTPUT_NAME "${PYFABM_NAME}"
  RUNTIME_OUTPUT_DIRECTORY "$<1:${PYFABM_DIR}>"
  LIBRARY_OUTPUT_DIRECTORY "$<1:${PYFABM_DIR}>"
)


cmake_minimum_required(VERSION 3.12)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(BUILD_SHARED_LIBS ON)

# Create main project.
project(python_fabm Fortran)

# Customize compiler flags
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
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

set(PYFABM_NAME "fabm" CACHE STRING "Base name of library")
set(PYFABM_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "Directory to create library in")
set(PYFABM_DEFINITIONS "" CACHE STRING "FABM host preprocessor definitions")
if(PYFABM_DEFINITIONS)
  add_compile_definitions(${PYFABM_DEFINITIONS})
endif()

# Force FABM_HOST to python and include FABM CMakeLists.txt.
set(FABM_FORCED_HOST python)
set(FABM_BUILD_C_API ON CACHE BOOL "Build C interfaces" FORCE)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../../.. fabm)

# Override cmake's default directory and name for the pyfabm dynamic library
# Note: the output directory is set with generator expressions to prevent
# multi-output generators (e.g. Visual Studio) from appending a per-configuration
# subdirectory. RUNTIME_OUTPUT_DIRECTORY applies on DLL platforms (Windows),
# LIBRARY_OUTPUT_DIRECTORY on other platforms.
# See also:
# * https://cmake.org/cmake/help/latest/prop_tgt/RUNTIME_OUTPUT_DIRECTORY.html
# * https://cmake.org/cmake/help/latest/prop_tgt/LIBRARY_OUTPUT_DIRECTORY.html
set_target_properties(fabm PROPERTIES
  OUTPUT_NAME "${PYFABM_NAME}"
  RUNTIME_OUTPUT_DIRECTORY "$<1:${PYFABM_DIR}>"
  LIBRARY_OUTPUT_DIRECTORY "$<1:${PYFABM_DIR}>"
)


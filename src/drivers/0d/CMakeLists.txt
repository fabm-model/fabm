cmake_minimum_required(VERSION 3.0)

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create main project.
project(fabm0d Fortran)

# Find GOTM directory.
# This directory will also be configurable in the CMake GUI or interactive mode.
find_path(GOTM_BASE
          NAMES src/gotm/gotm.F90
          HINTS ENV GOTMDIR
          PATHS ${PROJECT_SOURCE_DIR}/../../../../.. ${PROJECT_SOURCE_DIR}/../../..
          PATH_SUFFIXES GOTM/gotm-git gotm-git gotm
          DOC "Directory with GOTM source code (www.gotm.net)")

# Make NetCDF usage a configurable option.
option(FABM_USE_NetCDF "Enable NetCDF output" ON)

# Specify additional compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(WIN32)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
  endif()
endif()

set(GOTM_USE_NetCDF ${FABM_USE_NetCDF} CACHE BOOL "Enable output in NetCDF format" FORCE)
set(GOTM_USE_FABM OFF CACHE BOOL "use FABM" FORCE)
set(GOTM_USE_SEAGRASS OFF CACHE BOOL "use seagrass" FORCE)
set(GOTM_USE_STIM OFF CACHE BOOL "use STIM" FORCE)
mark_as_advanced(GOTM_USE_FABM GOTM_USE_NetCDF GOTM_USE_SEAGRASS GOTM_USE_STIM)
add_subdirectory(${GOTM_BASE} gotm EXCLUDE_FROM_ALL)

# Force FABM_HOST to 0d and include FABM CMakeLists.txt.
set(FABM_FORCED_HOST 0d)
add_subdirectory(${PROJECT_SOURCE_DIR}/../.. fabm EXCLUDE_FROM_ALL)

# Build 0d executable.
add_executable (fabm0d
                main.F90
                fabm0d.F90
                shared.F90
                output.F90
                fabm_0d.h
                )

# Link with FABM library, field/output manager and GOTM.
target_link_libraries(fabm0d PRIVATE fabm output_manager field_manager airsea input_manager util)

# Set module directory (output of *.mod)
set_property(TARGET fabm0d PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Set include directories with modules (*.mod)
target_include_directories(fabm0d PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/modules")

# Set include directories with header files (*.h) for GOTM - version.h, cppdefs.h, ode_solver_template.F90.
target_include_directories(fabm0d PRIVATE "${GOTM_BASE}/include")
target_include_directories(fabm0d PRIVATE "${GOTM_BASE}/src/util")

# Set include directories with header files (*.h) - our own and FABM's.
target_include_directories(fabm0d PRIVATE "${PROJECT_SOURCE_DIR}")
target_include_directories(fabm0d PRIVATE "${PROJECT_SOURCE_DIR}/../../../include")
target_include_directories(fabm0d PRIVATE $<TARGET_PROPERTY:yaml,INTERFACE_INCLUDE_DIRECTORIES>)

# Set default installation prefix.
if(WIN32)
  if(DEFINED ENV{LOCALAPPDATA})
    set(DEFAULT_PREFIX "$ENV{LOCALAPPDATA}/fabm/0d")
  else()
    set(DEFAULT_PREFIX "$ENV{APPDATA}/fabm/0d")
  endif()
else()
  set(DEFAULT_PREFIX "$ENV{HOME}/local/fabm/0d")
endif()

# Set cache variable with installation prefix.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # None set - set our default prefix.
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Directory to install FABM's 0d driver in" FORCE)
else()
  # User specified a prefix - just update the variable description.
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Directory to install FABM's 0d driver in")
endif()
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")

install(TARGETS fabm0d
  DESTINATION bin)

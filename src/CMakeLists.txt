cmake_minimum_required (VERSION 2.8)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(TARGET PROPERTY FOLDER library)

# Create main project.
project (fabm Fortran)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Default module directory (used by individual models)
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/model_modules)

# Add compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "-ffree-line-length-none")
endif()

# Make sure FABMHOST is set to a subdir of src/drivers
FILE(GLOB HOSTNAMES RELATIVE "${PROJECT_SOURCE_DIR}/drivers" "${PROJECT_SOURCE_DIR}/drivers/*")
set (FABMHOST "gotm" CACHE STRING "Host that FABM should be compiled for")
set_property(CACHE FABMHOST PROPERTY STRINGS ${HOSTNAMES})
message(STATUS "FABM host: ${FABMHOST}")

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS OR ${FABMHOST} STREQUAL "python")
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Set default installation prefix.
if(WIN32)
  set(DEFAULT_PREFIX "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/fabm/${FABMHOST}")
else()
  set(DEFAULT_PREFIX "$ENV{HOME}/local/fabm/${FABMHOST}")
endif()

# Re-initialize instalaltion prefix if needed.
mark_as_advanced(FABM_DEFAULT_PREFIX)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Installation prefix is still set to cmake default - overwrite it.
  set(FABM_INITIALIZE_PREFIX ON)
else()
  # Check whether cached actual and default prefixes match.
  # If so, we'll overwrite prefix with our current preferred (host-specific) default.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install FABM in")
  set(FABM_DEFAULT_PREFIX  ${DEFAULT_PREFIX} CACHE PATH "Default directory to install FABM in")
  if (CMAKE_INSTALL_PREFIX STREQUAL FABM_DEFAULT_PREFIX)
    set(FABM_INITIALIZE_PREFIX ON)
  endif()
endif()
if (FABM_INITIALIZE_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install FABM in" FORCE)
  set(FABM_DEFAULT_PREFIX  ${DEFAULT_PREFIX} CACHE PATH "Default directory to install FABM in" FORCE)
endif()
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")

# Global include directories
include_directories ("${PROJECT_SOURCE_DIR}/drivers/${FABMHOST}")
include_directories ("${PROJECT_SOURCE_DIR}/../include")
include_directories ("${CMAKE_BINARY_DIR}/modules")

add_library(fabm_base OBJECT
            fabm_driver.F90
            fabm_standard_variables.F90
            fabm_properties.F90
            fabm_types.F90
            fabm_particle.F90
            fabm_expressions.F90
           )

add_library(fabm
            $<TARGET_OBJECTS:fabm_base>
            $<TARGET_OBJECTS:aed>
            $<TARGET_OBJECTS:au>
            $<TARGET_OBJECTS:bb>
            $<TARGET_OBJECTS:examples>
            $<TARGET_OBJECTS:gotm>
            $<TARGET_OBJECTS:hzg>
            $<TARGET_OBJECTS:iow>
            $<TARGET_OBJECTS:klimacampus>
            $<TARGET_OBJECTS:metu>
            $<TARGET_OBJECTS:msi>
            $<TARGET_OBJECTS:niva>
            $<TARGET_OBJECTS:pml>
            fabm_library.F90
            fabm.F90
            config/fabm_config_types.F90
            config/fabm_yaml.F90
            config/fabm_config.F90
           )

# Store FABM *.mod in separate directory, so these files can later be installed.
set_property(TARGET fabm_base PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
set_property(TARGET fabm      PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Add directory with model-specific *.mod to include directories of main FABM library.
target_include_directories(fabm PRIVATE "${CMAKE_BINARY_DIR}/model_modules")

add_subdirectory(models)

if(EXISTS ${PROJECT_SOURCE_DIR}/drivers/${FABMHOST}/CMakeLists.txt)
  add_subdirectory (${PROJECT_SOURCE_DIR}/drivers/${FABMHOST})
endif()

if (NOT DEFINED FABM_NO_LIBRARY_INSTALL)
  install(TARGETS fabm
    DESTINATION lib)
  install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/\${BUILD_TYPE}/
    DESTINATION include)
  install(FILES ${PROJECT_SOURCE_DIR}/../include/fabm.h ${PROJECT_SOURCE_DIR}/drivers/${FABMHOST}/fabm_driver.h
    DESTINATION include)
endif()


cmake_minimum_required (VERSION 2.8)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(TARGET PROPERTY FOLDER library)

project (fabm Fortran )

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

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS OR ${FABMHOST} STREQUAL "python")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Set default installation prefix.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
if(WIN32)
  SET(CMAKE_INSTALL_PREFIX
    "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/fabm" CACHE PATH "Directory to install FABM in" FORCE
    )
else()
  SET(CMAKE_INSTALL_PREFIX
    "$ENV{HOME}/local/fabm" CACHE PATH "Directory to install FABM in" FORCE
    )
endif()
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

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

add_subdirectory (models)

if(EXISTS ${PROJECT_SOURCE_DIR}/drivers/${FABMHOST}/CMakeLists.txt)
add_subdirectory (${PROJECT_SOURCE_DIR}/drivers/${FABMHOST})
endif()

install(TARGETS fabm
  DESTINATION ${FABMHOST}/lib)
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/\${BUILD_TYPE}/
  DESTINATION ${FABMHOST}/include)
INSTALL(FILES ${PROJECT_SOURCE_DIR}/../include/fabm.h ${PROJECT_SOURCE_DIR}/drivers/${FABMHOST}/fabm_driver.h
  DESTINATION ${FABMHOST}/include)

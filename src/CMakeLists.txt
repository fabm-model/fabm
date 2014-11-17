cmake_minimum_required(VERSION 2.8.8)

# Create a list of all institutes to include.
# When adding new institutes to the source tree, please add them here as well.
set(FABM_INSTITUTES
    aed          # Aquatic Eco Dynamics, University of Western Australia
    au           # University of Aarhus
    bb           # Bolding & Burchard
    examples     # Examples supplied with FABM itself
    gotm         # Models ported from original GOTM/BIO library
    hzg          # Helmholtz-Zentrum Geesthacht
    iow          # Leibniz Institute for Baltic Sea Research
    klimacampus  # KlimaCampus Hamburg
    metu         # Middle East Technical University
    msi          # Marine Systems Institute, Tallinn University of Technology
    niva         # Norsk Institutt for Vannforskning
    pml          # Plymouth Marine Laboratory
   )

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create main project.
project(fabm Fortran)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(NOT WIN32)
    # Do not warn about Windows-specific export directives
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable 7841")
  else()
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 1023")
endif()

if(FABM_FORCED_HOST)
  set(FABM_HOST ${FABM_FORCED_HOST})
else()
  # Make sure FABM_HOST is set to a subdir of src/drivers
  FILE(GLOB HOSTNAMES RELATIVE "${PROJECT_SOURCE_DIR}/drivers" "${PROJECT_SOURCE_DIR}/drivers/*")
  set (FABM_HOST "gotm" CACHE STRING "Host that FABM should be compiled for")
  set_property(CACHE FABM_HOST PROPERTY STRINGS ${HOSTNAMES})
  message(STATUS "FABM host: ${FABM_HOST}")
endif()

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS OR ${FABM_HOST} STREQUAL "python")
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

# Global include directories
include_directories("${PROJECT_SOURCE_DIR}/drivers/${FABM_HOST}"
                    "${PROJECT_SOURCE_DIR}/../include"
                    "${CMAKE_CURRENT_BINARY_DIR}/modules"
                   )

add_library(fabm_base OBJECT
            fabm_driver.F90
            fabm_standard_variables.F90
            fabm_properties.F90
            fabm_types.F90
            fabm_particle.F90
            fabm_expressions.F90
            fabm_builtin_models.F90
            ../include/fabm.h
            ${PROJECT_SOURCE_DIR}/drivers/${FABM_HOST}/fabm_driver.h
           )

# Create a list of institute-specific object libraries.
foreach(institute ${FABM_INSTITUTES}) 
  add_subdirectory(models/${institute})
  list(APPEND MODEL_OBJECTS $<TARGET_OBJECTS:fabm_models_${institute}>)
endforeach(institute) 

add_library(fabm
            $<TARGET_OBJECTS:fabm_base>
            ${MODEL_OBJECTS}
            fabm_library.F90
            fabm_coupling.F90
            fabm.F90
            config/fabm_config_types.F90
            config/fabm_yaml.F90
            config/fabm_config.F90
           )

# Group config source files in a folder in IDEs that support this.
source_group(config REGULAR_EXPRESSION config/)

# Store FABM *.mod in separate directory, so these files can later be installed.
set_property(TARGET fabm_base fabm PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Add directories with model-specific *.mod to include directories of main FABM library.
foreach(institute ${FABM_INSTITUTES})
  set_property(TARGET fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/models/${institute}")
endforeach(institute)

# Run host-specific cmake scripts if present.
if(EXISTS ${PROJECT_SOURCE_DIR}/drivers/${FABM_HOST}/CMakeLists.txt)
  add_subdirectory (${PROJECT_SOURCE_DIR}/drivers/${FABM_HOST})
endif()

# Place projects in specific solution folders (Visual Studio only).
set_property(TARGET fabm_base fabm PROPERTY FOLDER fabm)
foreach(institute ${FABM_INSTITUTES})
  set_property(TARGET fabm_models_${institute} PROPERTY FOLDER fabm/models)
endforeach(institute)

if(NOT FABM_FORCED_HOST)

# Re-initialize installation prefix if needed.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR CMAKE_INSTALL_PREFIX STREQUAL PREVIOUS_DEFAULT_PREFIX)
  get_filename_component(DEFAULT_PREFIX "${DEFAULT_PREFIX}" ABSOLUTE)
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Directory to install FABM in" FORCE)
  set(PREVIOUS_DEFAULT_PREFIX "${DEFAULT_PREFIX}" CACHE INTERNAL "Default directory to install FABM in")
endif()
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")

# Install FABM library and include files, unless the host overruled this.
if (NOT DEFINED FABM_NO_LIBRARY_INSTALL)
  install(TARGETS fabm
    DESTINATION lib)
  install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/\${BUILD_TYPE}/
    DESTINATION include)
  install(FILES ${PROJECT_SOURCE_DIR}/../include/fabm.h ${PROJECT_SOURCE_DIR}/drivers/${FABM_HOST}/fabm_driver.h
    DESTINATION include)
endif()

# Begin testing section

add_custom_target(test_all COMMENT "building tests - nothing yet :-)")

# End of testing section

endif()

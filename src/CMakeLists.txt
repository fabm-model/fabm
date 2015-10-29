cmake_minimum_required(VERSION 2.8.8)

# List of contributing institutes that should be included in compilation.
# When adding new institutes to the source tree, please add them here as well.
# You can exclude institute directories from compilation by commenting them out.
set(DEFAULT_INSTITUTES
    au           # University of Aarhus, Denmark
    bb           # Bolding & Bruggeman - formerly Bolding & Burchard
    examples     # Examples supplied with FABM itself
    gotm         # Models ported from original GOTM/BIO library
    pclake       # The PCLake model - reference implementation
    pml          # Plymouth Marine Laboratory, United Kingdom
    ersem        # European Regional Seas Ecosytem Model
    niva         # Norsk Institutt for Vannforskning, Norway
)

# List of institutes that use manual edits to fabm_library.F90.
# This is DEPRECATED. Please use the new conventions described on
# this wiki and add your institute to the above FABM_INSTITUTES instead.
set(FABM_INSTITUTES_OLD
    aed          # Aquatic Eco Dynamics, University of Western Australia
    hzg          # Helmholtz-Zentrum Geesthacht, Germany
    iow          # Leibniz Institute for Baltic Sea Research, Germany
    klimacampus  # KlimaCampus Hamburg, Germany
    metu         # Middle East Technical University, Turkey
    msi          # Marine Systems Institute, Tallinn University of Technology, Estonia
   )

set(FABM_INSTITUTES ${DEFAULT_INSTITUTES} CACHE STRING "Institute directories to include during compilation.")

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create main project.
project(fabm Fortran)

# Create fabm_library.F90 with all active institutes,
# and set final paths to institute source directories (FABM_${INSTITUTE_UC}_DIR)
set(FABM_LIBRARY_USE "")
set(FABM_LIBRARY_ADD "")
set(ACTIVE_INSTITUTES ${FABM_INSTITUTES_OLD})
foreach(institute ${FABM_INSTITUTES})
  # First retrieve the effective source directory
  set(DEFAULT_INSTITUTE_BASE "${PROJECT_SOURCE_DIR}/models/${institute}")
  string(TOUPPER ${institute} INSTITUTE_UC)
  find_path(FABM_${INSTITUTE_UC}_BASE
    NAMES CMakeLists.txt
    HINTS "${DEFAULT_INSTITUTE_BASE}"
    DOC "Path to source code of ${institute} models"
    NO_DEFAULT_PATH
  )
  if(EXISTS "${DEFAULT_INSTITUTE_BASE}/CMakeLists.txt")
    mark_as_advanced(FABM_${INSTITUTE_UC}_BASE)
  endif()

  if(FABM_${INSTITUTE_UC}_BASE)
    # Either a custom institute directory has been specified, or the default institute directory exists.
    set(FABM_LIBRARY_USE "${FABM_LIBRARY_USE}
   use ${institute}_model_library")
    set(FABM_LIBRARY_ADD "${FABM_LIBRARY_ADD}
      call self%add(${institute}_model_factory,'${institute}')")
    set(ACTIVE_INSTITUTES ${ACTIVE_INSTITUTES} ${institute})
  else()
    #message(STATUS "${institute} will not be compiled because FABM_${INSTITUTE_UC}_BASE is not set. Please visit http://fabm.net/${institute} for more information.")
  endif()
endforeach(institute)
configure_file(${PROJECT_SOURCE_DIR}/fabm_library.F90.in
               ${CMAKE_CURRENT_BINARY_DIR}/fabm_library.F90)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  if(WIN32)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Od")
  endif()
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -N 1023")
endif()

if(FABM_FORCED_HOST)
  # FABM_HOST was specified by a higher level CMakeLists.txt
  set(FABM_HOST ${FABM_FORCED_HOST})
else()
  # FABM_HOST is user-configurable.
  # Make sure FABM_HOST is set to a subdir of src/drivers
  FILE(GLOB HOSTNAMES RELATIVE "${PROJECT_SOURCE_DIR}/drivers" "${PROJECT_SOURCE_DIR}/drivers/*")
  set (FABM_HOST "gotm" CACHE STRING "Host that FABM should be compiled for")
  set_property(CACHE FABM_HOST PROPERTY STRINGS ${HOSTNAMES})
  message(STATUS "FABM host: ${FABM_HOST}")
endif()

# Use host-prescribed real kind if set.
if(FABM_REAL_KIND)
  add_definitions(-D_FABM_REAL_KIND_=${FABM_REAL_KIND})
endif()

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
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
            ../include/fabm_private.h
            ${PROJECT_SOURCE_DIR}/drivers/${FABM_HOST}/fabm_driver.h
           )

# Create a list of institute-specific object libraries.
foreach(institute ${ACTIVE_INSTITUTES})
  string(TOUPPER ${institute} INSTITUTE_UC)
  if(NOT DEFINED FABM_${INSTITUTE_UC}_BASE)
    set(FABM_${INSTITUTE_UC}_BASE "${PROJECT_SOURCE_DIR}/models/${institute}")
  endif()
  add_subdirectory(${FABM_${INSTITUTE_UC}_BASE} ${CMAKE_CURRENT_BINARY_DIR}/models/${institute})
  set_property(TARGET fabm_models_${institute} PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/models/${institute}")
  set_property(TARGET fabm_models_${institute} APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/models/${institute}")
  list(APPEND MODEL_OBJECTS $<TARGET_OBJECTS:fabm_models_${institute}>)
endforeach(institute) 

add_library(fabm
            $<TARGET_OBJECTS:fabm_base>
            ${MODEL_OBJECTS}
            ${CMAKE_CURRENT_BINARY_DIR}/fabm_library.F90
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
foreach(institute ${ACTIVE_INSTITUTES})
  set_property(TARGET fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/models/${institute}")
endforeach(institute)

# Place projects in specific solution folders (Visual Studio only).
set_property(TARGET fabm_base fabm PROPERTY FOLDER fabm)
foreach(institute ${ACTIVE_INSTITUTES})
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

add_executable(test_yaml EXCLUDE_FROM_ALL
               config/test_yaml.F90
               config/fabm_config_types.F90
               config/fabm_yaml.F90
              )

# End of testing section

endif()

# List of contributing institutes that should be included in compilation.
# When adding new institutes to the source tree, please add them here as well.
# You can exclude institute directories from compilation by commenting them out.
set(DEFAULT_INSTITUTES
    akvaplan     # Akvaplan-niva, Norway
    au           # University of Aarhus, Denmark
    bb           # Bolding & Bruggeman - formerly Bolding & Burchard
    csiro        # Commonwealth Scientific and Industrial Research Organisation, Australia
    ersem        # European Regional Seas Ecosystem Model, https://doi.org/10.5194/gmd-9-1293-2016
    examples     # Examples supplied with FABM itself
    gotm         # Models ported from original GOTM/BIO library
    ihamocc      # isopycnic coordinate HAMburg Ocean Carbon Cycle model, https://doi.org/10.5194/gmd-13-2393-2020
    iow          # Leibniz Institute for Baltic Sea Research, Germany
    jrc          # EC - Joint Research Centre, Ispra, Italy
    mops         # Model of Oceanic Pelagic Stoichiometry, https://doi.org/10.5194/gmd-8-2929-2015
    msi          # Marine Systems Institute, Tallinn University of Technology, Estonia
    nersc        # Nansen Environmental and Remote Sensing Center
    niva         # Norwegian Institute for Water Research, Norway
    ogs          # National Institute of Oceanography and Applied Geophysics - OGS, Italy
    pclake       # The PCLake model - reference implementation
    pisces       # Pelagic Interactions Scheme for Carbon and Ecosystem Studies, https://doi.org/10.5194/gmd-8-2465-2015
    pml          # Plymouth Marine Laboratory, United Kingdom
    selma        # Simple EcoLogical Model for the Aquatic - PROGNOS
    su           # Swansea University, United Kingdom
    uhh          # University of Hamburg, Germany
)

set(OPTIONAL_INSTITUTES ersem mops pclake pisces ihamocc ogs)

if(CMAKE_CURRENT_LIST_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(FATAL_ERROR "Please point cmake to the top-level directory with FABM source code, not to its src subdirectory.")
endif()

# Ensure FABM_INSTITUTES tracks the up-to-date DEFAULT_INSTITUTES list unless the user manually customized FABM_INSTITUTES.
if(NOT FABM_DEFAULT_INSTITUTES OR NOT FABM_INSTITUTES)
  set(FABM_DEFAULT_INSTITUTES "${DEFAULT_INSTITUTES}" CACHE INTERNAL "Default institutes at the time FABM_INSTITUTES was initialized or updated automatically.")
endif()
set(FABM_INSTITUTES ${DEFAULT_INSTITUTES} CACHE STRING "Institute directories to include during compilation.")
if(FABM_INSTITUTES STREQUAL FABM_DEFAULT_INSTITUTES)
  # FABM_INSTITUTES matches its original default value; update both with the current DEFAULT_INSTITUTES list.
  set(FABM_INSTITUTES ${DEFAULT_INSTITUTES} CACHE STRING "Institute directories to include during compilation." FORCE)
  set(FABM_DEFAULT_INSTITUTES ${DEFAULT_INSTITUTES} CACHE INTERNAL "Default institutes at the time FABM_INSTITUTES was initialized or updated automatically.")
endif()

set(FABM_EXTRA_INSTITUTES "" CACHE STRING "Additional institute directories to include during compilation")

# Create fabm_library.F90 with all active institutes,
# and set final paths to institute source directories (FABM_${INSTITUTE_UC}_DIR)
set(FABM_LIBRARY_USE "")
set(FABM_LIBRARY_ADD "")
set(ACTIVE_INSTITUTES "")
set(ALL_INSTITUTES ${FABM_INSTITUTES} ${FABM_EXTRA_INSTITUTES})
list(REMOVE_DUPLICATES ALL_INSTITUTES)
foreach(institute IN LISTS ALL_INSTITUTES)
  # First retrieve the effective source directory
  string(TOUPPER ${institute} INSTITUTE_UC)
  find_path(FABM_${INSTITUTE_UC}_BASE
    NAMES CMakeLists.txt
    HINTS "${CMAKE_CURRENT_LIST_DIR}/../extern/${institute}" "${CMAKE_CURRENT_LIST_DIR}/models/${institute}"
    DOC "Path to source code of ${institute} models"
    NO_DEFAULT_PATH
    NO_CMAKE_FIND_ROOT_PATH
  )
  if(institute IN_LIST DEFAULT_INSTITUTES AND NOT institute IN_LIST OPTIONAL_INSTITUTES)
    mark_as_advanced(FABM_${INSTITUTE_UC}_BASE)
  endif()

  if(FABM_${INSTITUTE_UC}_BASE)
    # Either a custom institute directory has been specified, or the default institute directory exists.
    if (EXISTS "${FABM_${INSTITUTE_UC}_BASE}")
      # The institute directory exists
      string(APPEND FABM_LIBRARY_USE "
   use ${institute}_model_library")
      string(APPEND FABM_LIBRARY_ADD "
      call self%add(${institute}_model_factory, '${institute}')")
      list(APPEND ACTIVE_INSTITUTES ${institute})
    else()
      # The institute directory does not exist and it was explicitly specified (because different from default).
      message(FATAL_ERROR "Source directory ${FABM_${INSTITUTE_UC}_BASE} specified for institute ${institute} not found. Please ensure FABM_${INSTITUTE_UC}_BASE is set to a valid directory.")
    endif()
  elseif(institute IN_LIST OPTIONAL_INSTITUTES)
    #message(STATUS "${institute} will not be compiled because FABM_${INSTITUTE_UC}_BASE is not set. Please visit https://fabm.net/${institute} for more information.")
  else()
    message(FATAL_ERROR "Source directory for institute ${institute} not found. Please specify FABM_${INSTITUTE_UC}_BASE.")
  endif()
endforeach(institute)
configure_file(${CMAKE_CURRENT_LIST_DIR}/fabm_library.F90.in
               ${CMAKE_CURRENT_BINARY_DIR}/fabm_library.F90)

if(FABM_EMBED_VERSION)
  # Use git to get description of current FABM tag/commit id and branch
  set(GET_GIT_INFO_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules/GetGitInfo.cmake")
  add_custom_target(fabm_version
     ${CMAKE_COMMAND} -DINFILE=${CMAKE_CURRENT_LIST_DIR}/fabm_version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/fabm_version.F90 -P ${GET_GIT_INFO_SCRIPT}
     BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/fabm_version.F90
     WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
     COMMENT "Retrieving description of last FABM commit..."
     SOURCES ${CMAKE_CURRENT_LIST_DIR}/fabm_version.F90.in ${GET_GIT_INFO_SCRIPT}
     VERBATIM
  )
  set_property(TARGET fabm_version PROPERTY FOLDER fabm)
else()
  # Insert dummy current FABM tag/commit id and branch
  set(GIT_COMMIT_ID unknown)
  set(GIT_BRANCH_NAME unknown)
  configure_file(${CMAKE_CURRENT_LIST_DIR}/fabm_version.F90.in ${CMAKE_CURRENT_BINARY_DIR}/fabm_version.F90)
endif()

configure_file(${CMAKE_CURRENT_LIST_DIR}/fabm_parameters.F90.in ${CMAKE_CURRENT_BINARY_DIR}/fabm_parameters.F90)

add_library(fabm_base OBJECT
            ${CMAKE_CURRENT_BINARY_DIR}/fabm_parameters.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_driver.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_standard_variables.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_properties.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_types.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_particle.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_expressions.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/constant.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/reduction.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/sum.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/scale.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/source.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/time_filter.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/depth_integral.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/relaxation.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/depth_mapping.F90
            ${CMAKE_CURRENT_LIST_DIR}/builtin/model_library.F90
            ${CMAKE_CURRENT_LIST_DIR}/../include/fabm.h
            ${CMAKE_CURRENT_LIST_DIR}/../include/fabm_private.h
            ${CMAKE_CURRENT_LIST_DIR}/../include/fabm_version.h
            ${CMAKE_CURRENT_LIST_DIR}/drivers/${FABM_HOST}/fabm_driver.h
           )
source_group(builtin REGULAR_EXPRESSION /builtin/)
target_link_libraries(fabm_base PUBLIC yaml)
target_include_directories(fabm_base PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/drivers/${FABM_HOST}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../include>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS fabm_base EXPORT fabmConfig)

# Define preprocessor macros. These are set on the fabm_base target, but also
# at directory level (add_compile_definitions) for institutes that still use
# add_dependencies instead of target_link_libraries to link to fabm_base.
if(NOT FABM_USE_IEEE_ARITHMETIC)
  add_compile_definitions(_FABM_NO_IEEE_ARITHMETIC_)
  target_compile_definitions(fabm_base PUBLIC _FABM_NO_IEEE_ARITHMETIC_)
endif()
if(NOT FABM_USE_DO_CONCURRENT)
  add_compile_definitions(_FABM_NO_DO_CONCURRENT_)
  target_compile_definitions(fabm_base PUBLIC _FABM_NO_DO_CONCURRENT_)
endif()
if(NOT FABM_USE_CONTIGUOUS)
  add_compile_definitions(_FABM_NO_CONTIGUOUS_)
  target_compile_definitions(fabm_base PUBLIC _FABM_NO_CONTIGUOUS_)
endif()

add_library(fabm
            ${CMAKE_CURRENT_BINARY_DIR}/fabm_library.F90
            ${CMAKE_CURRENT_BINARY_DIR}/fabm_version.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_job.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_task_order.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_graph.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_coupling.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_schedule.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_config.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_debug.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_work.F90
            ${CMAKE_CURRENT_LIST_DIR}/fabm_v0_compatibility.F90
           )
target_link_libraries(fabm PUBLIC yaml fabm_base)

if(FABM_EMBED_VERSION)
  add_dependencies(fabm fabm_version)
endif()

# Store FABM *.mod in separate directory, so these files can later be installed.
set_property(TARGET fabm_base fabm PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

# Global include directories (deprecated; set only for institutes that link to
# fabm_base using add_dependencies instead of target_link_libraries)
include_directories($<TARGET_PROPERTY:fabm_base,INTERFACE_INCLUDE_DIRECTORIES>)

# Create a list of institute-specific object libraries.
set(CMAKE_FOLDER fabm/models)
foreach(institute IN LISTS ACTIVE_INSTITUTES)
  string(TOUPPER ${institute} INSTITUTE_UC)
  add_subdirectory(${FABM_${INSTITUTE_UC}_BASE} ${CMAKE_CURRENT_BINARY_DIR}/models/${institute})
  set_property(TARGET fabm_models_${institute} PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/models/${institute}")
  target_include_directories(fabm_models_${institute} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/models/${institute}>)
  target_link_libraries(fabm PRIVATE fabm_models_${institute})
  install(TARGETS fabm_models_${institute} EXPORT fabmConfig)
endforeach(institute)

if(FABM_BUILD_C_API)
  add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/c" c)
  source_group(c REGULAR_EXPRESSION /c/)
endif()

# Tests
add_executable(test_host EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/test/host.F90
               ${CMAKE_CURRENT_LIST_DIR}/test/models.F90
              )
target_link_libraries(test_host PRIVATE fabm)

# Place projects in specific solution folders (Visual Studio only).
set_property(TARGET fabm_base fabm test_host PROPERTY FOLDER fabm)

if(NOT FABM_FORCED_HOST)
  # Re-initialize installation prefix if needed.
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR CMAKE_INSTALL_PREFIX STREQUAL PREVIOUS_DEFAULT_PREFIX)
    get_filename_component(DEFAULT_PREFIX "${DEFAULT_PREFIX}" ABSOLUTE)
    set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Directory to install FABM in" FORCE)
    set(PREVIOUS_DEFAULT_PREFIX "${DEFAULT_PREFIX}" CACHE INTERNAL "Default directory to install FABM in")
  endif()
  message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
endif()

# Install FABM library and include files, unless the host overruled this.
if (NOT DEFINED FABM_NO_LIBRARY_INSTALL)
  # Note: using an implicit destination for install(TARGETS ...) below,
  # or TYPE instead of DESTINATION, requires cmake >= 3.14
  install(TARGETS fabm
    EXPORT fabmConfig
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/\${BUILD_TYPE}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(FILES
    ${CMAKE_CURRENT_LIST_DIR}/../include/fabm.h
    ${CMAKE_CURRENT_LIST_DIR}/../include/fabm_version.h
    ${CMAKE_CURRENT_LIST_DIR}/drivers/${FABM_HOST}/fabm_driver.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

export(EXPORT fabmConfig)
install(EXPORT fabmConfig DESTINATION cmake)

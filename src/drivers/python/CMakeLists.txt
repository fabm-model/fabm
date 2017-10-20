cmake_minimum_required(VERSION 2.8.8)

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

# Force FABM_HOST to python and include FABM CMakeLists.txt.
set(FABM_FORCED_HOST python)
add_subdirectory(${PROJECT_SOURCE_DIR}/../.. fabm)

# Add Python-FABM library
add_library(python_fabm SHARED
            python_fabm.F90
            variable.F90
            link_list.F90
            c_helper.F90
            helper.F90
           )

# Link in FABM itself
target_link_libraries(python_fabm fabm)

# Set module directory (output of *.mod)
set_property(TARGET python_fabm PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Set include directories with modules (*.mod) - our own and FABM's.
set_property(TARGET python_fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/modules")
set_property(TARGET python_fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/fabm/modules")

# Set include directories with header files (*.h) - our own and FABM's.
set_property(TARGET python_fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}")
set_property(TARGET python_fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/../../../include")

# Query Python itself for default installation prefix.
find_package(PythonInterp)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -m site --user-site
                OUTPUT_VARIABLE SITE_PACKAGES_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)
set(DEFAULT_PREFIX "${SITE_PACKAGES_DIR}")

# Set cache variable with installation prefix.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # None set - set our default prefix.
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Directory to install python_fabm in" FORCE)
else()
  # User specified a prefix - just update the variable description.
  set(CMAKE_INSTALL_PREFIX "${DEFAULT_PREFIX}" CACHE PATH "Directory to install python_fabm in")
endif()
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")

add_custom_command(OUTPUT wheel.stamp POST_BUILD COMMAND pip2 wheel ${PROJECT_SOURCE_DIR} DEPENDS python_fabm VERBATIM)

add_custom_target(wheel DEPENDS wheel.stamp)

#KB add_custom_target(wheel_upload COMMAND cp *.whl /tmp DEPENDS wheel) # I think it is called twine

install(TARGETS python_fabm
  DESTINATION pyfabm)
install(DIRECTORY pyfabm/
  DESTINATION pyfabm)


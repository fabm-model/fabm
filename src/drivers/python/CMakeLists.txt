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

# Force FABM_HOST to python and include FABM CMakeLists.txt.
set(FABM_FORCED_HOST python)
add_subdirectory(${PROJECT_SOURCE_DIR}/../.. fabm EXCLUDE_FROM_ALL)

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

# Find target Python interpreter.
find_package(PythonInterp)

# Build wheel
add_custom_command(OUTPUT wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pyfabm" pyfabm
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:python_fabm>" pyfabm/
                   COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --user --force-reinstall --upgrade pyfabm WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")

add_custom_target(wheel_upload COMMAND twine upload dist/*.whl DEPENDS wheel)


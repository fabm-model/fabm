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
set(FABM_BUILD_C_API ON CACHE BOOL "Build C interfaces" FORCE)
add_subdirectory(${PROJECT_SOURCE_DIR}/../.. fabm EXCLUDE_FROM_ALL)

# Add Python-FABM library
#add_library(python_fabm SHARED
#            pyfabm/__init__.py
#           )
#source_group(pyfabm FILES REGULAR_EXPRESSION ".*\\.py")

# Find target Python interpreter.
find_package(PythonInterp)

# Build wheel
add_custom_command(OUTPUT wheel.stamp
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pyfabm" pyfabm
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:fabm_c>" pyfabm/
                   COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp SOURCES pyfabm/__init__.py)

# Let pip install the wheel
install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist --user --force-reinstall --upgrade pyfabm WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")

add_custom_target(wheel_upload COMMAND twine upload dist/*.whl DEPENDS wheel)


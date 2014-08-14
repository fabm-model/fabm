# Find GOTM directory.
# This directory will also be configurable in the CMake GUI or interactive mode.
find_path(GOTMDIR
          NAMES src/gotm/gotm.F90
          HINTS ENV GOTMDIR
          PATH_SUFFIXES GOTM/gotm-git gotm-git
          DOC "Directory with GOTM source code (www.gotm.net)")

# generate GOTM's version.h
configure_file("${GOTMDIR}/include/version.h.in" "${CMAKE_BINARY_DIR}/gotm/version.h")

# Build 0d executable.
add_executable (fabm0d main.F90 fabm0d.F90 shared.F90 output.F90 
                ${GOTMDIR}/src/airsea/short_wave_radiation.F90
                ${GOTMDIR}/src/input/input.F90
                ${GOTMDIR}/src/util/eqstate.F90
                ${GOTMDIR}/src/util/gridinterpol.F90
                ${GOTMDIR}/src/util/time.F90
                )

# Group GOTM source files in a folder in IDEs that support this.
source_group(gotm REGULAR_EXPRESSION ${GOTMDIR})

# Give exectuable project its own folder in Visual Studio solutions.
set_property(TARGET fabm0d PROPERTY FOLDER executable)

# Link with FABM library.
target_link_libraries(fabm0d fabm)

# Set module directory (output of *.mod) and include directories (input pf *.h, *.mod).
set_property(TARGET fabm0d PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/0d_modules)
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${GOTMDIR}/include")
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${GOTMDIR}/src/util")
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/0d_modules")
set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/gotm")
if (WIN32)
  set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES "${GOTMDIR}/compilers/vs2010")
endif(WIN32)

# Make NetCDF use a configurable option.
option(USE_NETCDF "Enable NetCDF output" ON)
if(USE_NETCDF)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${GOTMDIR}/src/cmake/Modules/")
  find_package(NetCDF REQUIRED)
  target_link_libraries(fabm0d ${NetCDF_LIBRARIES})
  set_property(TARGET fabm0d APPEND PROPERTY INCLUDE_DIRECTORIES ${NetCDF_INCLUDE_DIRS})
  set_property(TARGET fabm0d APPEND PROPERTY COMPILE_DEFINITIONS NETCDF4)
  if (NetCDF_STATIC_MSVC_BUILD)
    set_property(TARGET fabm0d PROPERTY LINK_FLAGS_DEBUG "/NODEFAULTLIB:\"libcmtd\"")
    target_compile_options(fabm0d PRIVATE "/libs:static")
  endif()
endif(USE_NETCDF)

install(TARGETS fabm0d
  DESTINATION bin)

set(FABM_NO_LIBRARY_INSTALL ON PARENT_SCOPE)

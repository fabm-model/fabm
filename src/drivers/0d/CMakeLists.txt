# Find GOTM directory.
# This directory will also be configurable in the CMake GUI or interactive mode.
find_path(GOTMDIR
          NAMES src/gotm/gotm.F90
          HINTS ENV GOTMDIR
          PATH_SUFFIXES GOTM/gotm-git gotm-git
          DOC "Directory with GOTM source code (www.gotm.net)")

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
target_include_directories(fabm0d PRIVATE "${GOTMDIR}/include")
target_include_directories(fabm0d PRIVATE "${GOTMDIR}/src/util")
target_include_directories(fabm0d PRIVATE "${CMAKE_BINARY_DIR}/0d_modules")
if (WIN32)
  target_include_directories(fabm0d PRIVATE "${GOTMDIR}/compilers/vs2010")
endif(WIN32)

# Make NetCDF use a configurable option.
option(USE_NETCDF "Enable NetCDF output" ON)

if(USE_NETCDF)

# First try to locate nf-config.
find_program( NETCDF_CONFIG_EXECUTABLE
   NAMES nf-config
   HINTS ENV NETCDF_ROOT
   PATH_SUFFIXES bin Bin
   DOC "NetCDF config program. Used to detect NetCDF include directory and linker flags." )
mark_as_advanced(NETCDF_CONFIG_EXECUTABLE)

if(NETCDF_CONFIG_EXECUTABLE)

# Found nf-config - use it to retrieve include directory and linking flags.
# Mark NetCDF paths as advanced configuration options in CMake (hidden by default).
exec_program(${NETCDF_CONFIG_EXECUTABLE} ARGS --includedir OUTPUT_VARIABLE includedir)
exec_program(${NETCDF_CONFIG_EXECUTABLE} ARGS --flibs      OUTPUT_VARIABLE flibs)
set(NETCDF_INCLUDE_DIRS ${includedir} CACHE STRING "NetCDF include directories")
set(NETCDF_LIBRARIES ${flibs} CACHE STRING "NetCDF linking flags")
mark_as_advanced(NETCDF_INCLUDE_DIRS NETCDF_FLIBS)

elseif (WIN32)

# On Windows: use CMake to locte paths; default to NetCDF static library provided with GOTM.
find_library(NETCDF_LIBRARIES NAMES netcdfs
             HINTS ${GOTMDIR}/extras/netcdf/win32/3.6.3/lib
             DOC "NetCDF library")
find_path(NETCDF_INCLUDE_DIRS netcdf.mod
          HINTS ${GOTMDIR}/extras/netcdf/win32/3.6.3/include ENV NETCDFINC
          DOC "NetCDF include directory")

else()

# Use GOTM environment variables: NETCDFLIBNAME, NETCDFLIBDIR, NETCDFINC
if(DEFINED ENV{NETCDFLIBNAME})
  set(flibs $ENV{NETCDFLIBNAME})
else()
  set(flibs "-lnetcdf")
endif()
if(DEFINED ENV{NETCDFLIBDIR})
  set(flibs "${flibs} -L$ENV{NETCDFLIBDIR}")
endif()
set(NETCDF_LIBRARIES ${flibs} CACHE STRING "NetCDF linking flags")
find_path(NETCDF_INCLUDE_DIRS netcdf.mod
          HINTS ENV NETCDFINC
          DOC "NetCDF include directory")

endif()

target_include_directories(fabm0d PRIVATE ${NETCDF_INCLUDE_DIRS})
target_link_libraries(fabm0d ${NETCDF_LIBRARIES})
target_compile_definitions(fabm0d PRIVATE NETCDF4)

endif(USE_NETCDF)

install(TARGETS fabm0d
  DESTINATION bin)

set(FABM_NO_LIBRARY_INSTALL ON PARENT_SCOPE)

add_library(fabm_omp
            fabm_omp.F90
           )

# Link in FABM itself
target_link_libraries(fabm_omp fabm)

# Set module directory (where *.mod will be saved by the compiler)
set_property(TARGET fabm_omp PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories(fabm_omp PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>)

set_property(TARGET fabm_omp PROPERTY FOLDER fabm)

find_package(OpenMP)
target_link_libraries(fabm_omp OpenMP::OpenMP_Fortran)


add_library(fabm_c SHARED
            fabm_c.F90
            variable.F90
            parameter.F90
            link_list.F90
            c_helper.F90
            helper.F90
            integrate.F90
           )

# Link in FABM itself
target_link_libraries(fabm_c fabm)

# Set module directory (where *.mod will be saved by the compiler)
set_property(TARGET fabm_c PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories(fabm_c PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules> $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../include> $<INSTALL_INTERFACE:include>)

set_property(TARGET fabm_c PROPERTY FOLDER fabm)
cmake_minimum_required(VERSION 2.8.8)

# Use solution folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Create main project.
project(fabm_particle_driver Fortran)

# Force FABM_HOST to 0d and include FABM CMakeLists.txt.
set(FABM_FORCED_HOST particles)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(${PROJECT_SOURCE_DIR}/../.. fabm)

# Build 0d executable.
add_library(fabm_particle_driver SHARED
            particle_driver.F90
           )

# Link with FABM library.
target_link_libraries(fabm_particle_driver PRIVATE fabm_particles)

set_property(TARGET fabm_particle_driver PROPERTY FOLDER fabm_particles)

# Set module directory (output of *.mod)
set_property(TARGET fabm_particle_driver PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Set include directories with modules (*.mod) - our own and FABM's.
set_property(TARGET fabm_particle_driver APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/modules")
set_property(TARGET fabm_particle_driver APPEND PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/fabm/modules")

# Set include directories with header files (*.h) - our own and FABM's.
set_property(TARGET fabm_particle_driver APPEND PROPERTY INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}")
set_property(TARGET fabm_particle_driver APPEND PROPERTY INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/../../../include")

#set_property(TARGET fabm_particle_driver PROPERTY Fortran_VISIBILITY_PRESET hidden)


cmake_minimum_required(VERSION 2.8.8)

# Create main project.
project(${PROJECT_PREFIX}_yaml Fortran)

add_library(${PROJECT_PREFIX}_yaml OBJECT
            yaml_types.F90
            yaml.F90
           )

add_executable(test_${PROJECT_PREFIX}_yaml EXCLUDE_FROM_ALL
               $<TARGET_OBJECTS:${PROJECT_PREFIX}_yaml>
               test_yaml.F90
              )

cmake_minimum_required(VERSION 2.8.8)

# Create main project.
project(yaml Fortran)

add_library(yaml OBJECT
            yaml_types.F90
            yaml.F90
           )

add_executable(test_yaml EXCLUDE_FROM_ALL
               $<TARGET_OBJECTS:yaml>
               test_yaml.F90
              )

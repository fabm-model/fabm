cmake_minimum_required(VERSION 2.8.8)

if(POLICY CMP0048)
   cmake_policy(SET CMP0048 NEW)
endif()

if(TARGET yaml)
   return()
endif(TARGET yaml)

# Create main project.
project(yaml Fortran)

add_library(yaml OBJECT
            yaml_types.F90
            yaml.F90
           )
set(YAML_SOURCES $<TARGET_OBJECTS:yaml> PARENT_SCOPE)
set_property(TARGET yaml PROPERTY Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_include_directories(yaml
   PUBLIC
      $<INSTALL_INTERFACE:include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
)

add_executable(test_yaml EXCLUDE_FROM_ALL
               ${YAML_SOURCES}
               test_yaml.F90
              )

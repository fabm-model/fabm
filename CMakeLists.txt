cmake_minimum_required(VERSION 3.0)

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

if(NOT PROJECT_NAME STREQUAL fabm)
  project(fabm Fortran)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/src/CMakeLists.txt")


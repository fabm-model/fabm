cmake_minimum_required(VERSION 3.0)

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()

if(NOT PROJECT_NAME STREQUAL fabm)
  project(fabm VERSION 1.0.0 LANGUAGES Fortran)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/src/CMakeLists.txt")


# 3.14 for https://cmake.org/cmake/help/latest/module/FindPython3.html with NumPy
cmake_minimum_required(VERSION 3.14)

find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_custom_command(OUTPUT fabm_base.c fabm_base.h
  COMMAND cython fabm_base.pyx -o ${CMAKE_CURRENT_BINARY_DIR}/fabm_base.c -3
  MAIN_DEPENDENCY fabm_base.pyx
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sys;import os;print(sys.prefix+os.pathsep+sys.exec_prefix if sys.prefix!=sys.exec_prefix else sys.prefix)"
  OUTPUT_VARIABLE PYTHONHOME
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)

configure_file(python_parameters.F90.in python_parameters.F90)

enable_language(C)
add_library(fabm_models_python_c OBJECT
  ${CMAKE_CURRENT_BINARY_DIR}/fabm_base.c
  ${CMAKE_CURRENT_BINARY_DIR}/fabm_base.h
  embed.c
)
target_link_libraries(fabm_models_python_c Python3::Python Python3::NumPy)
target_include_directories(fabm_models_python_c PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(fabm_models_python_c PRIVATE -DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION)
set_property(TARGET fabm_models_python_c PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
install(TARGETS fabm_models_python_c EXPORT fabmConfig)

add_library(fabm_models_python OBJECT
  python_model_library.F90
  wrapped_python_model.F90
  python_parameters.F90.in
  ${CMAKE_CURRENT_BINARY_DIR}/python_parameters.F90
)
target_link_libraries(fabm_models_python PRIVATE fabm_base fabm_models_python_c)
target_sources(fabm PRIVATE $<TARGET_OBJECTS:fabm_models_python_c>)

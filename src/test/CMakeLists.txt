add_executable(test_host EXCLUDE_FROM_ALL
               ${CMAKE_CURRENT_LIST_DIR}/shared.F90
               ${CMAKE_CURRENT_LIST_DIR}/host.F90
               ${CMAKE_CURRENT_LIST_DIR}/models.F90
)
target_link_libraries(test_host PRIVATE fabm)

if(FABM_BUILD_OMP_API)
  add_executable(test_omp EXCLUDE_FROM_ALL
                ${CMAKE_CURRENT_LIST_DIR}/shared.F90
                ${CMAKE_CURRENT_LIST_DIR}/omp.F90
  )
  find_package(OpenMP)
  target_link_libraries(test_omp PRIVATE fabm fabm_omp OpenMP::OpenMP_Fortran)
endif()

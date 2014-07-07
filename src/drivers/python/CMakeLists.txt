add_library (fabm_python SHARED fabm_python.F90 helper.F90 )
target_link_libraries(fabm_python fabm)
install(TARGETS fabm_python
  DESTINATION ${FABMHOST}/python/pyfabm)
install(DIRECTORY pyfabm
  DESTINATION ${FABMHOST}/python)

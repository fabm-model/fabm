add_library (python_fabm SHARED python_fabm.F90 helper.F90 )
target_link_libraries(python_fabm fabm)
install(TARGETS python_fabm
  DESTINATION ${FABMHOST}/python/pyfabm)
install(DIRECTORY pyfabm
  DESTINATION ${FABMHOST}/python)

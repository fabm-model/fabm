add_library (python_fabm SHARED python_fabm.F90 helper.F90 )
target_link_libraries(python_fabm fabm)
set_property(TARGET python_fabm PROPERTY FOLDER python)

install(TARGETS python_fabm
  DESTINATION pyfabm)
install(DIRECTORY pyfabm/
  DESTINATION pyfabm)

set(FABM_NO_LIBRARY_INSTALL ON PARENT_SCOPE)
 

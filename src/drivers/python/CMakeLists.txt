add_library (python_fabm SHARED python_fabm.F90 helper.F90 )
target_link_libraries(python_fabm fabm)
set_property(TARGET python_fabm PROPERTY FOLDER python)

find_package(PythonInterp)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -m site --user-site
                OUTPUT_VARIABLE SITE_PACKAGES_DIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)
set(DEFAULT_PREFIX "${SITE_PACKAGES_DIR}" PARENT_SCOPE)

install(TARGETS python_fabm
  DESTINATION pyfabm)
install(DIRECTORY pyfabm/
  DESTINATION pyfabm)

set(FABM_NO_LIBRARY_INSTALL ON PARENT_SCOPE)


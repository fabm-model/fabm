file(COPY "${SOURCE_DIR}/pyfabm" DESTINATION .)
file(COPY "${PYFABM_LIB}" DESTINATION pyfabm)
configure_file("${SOURCE_DIR}/setup.py.in" setup.py)
execute_process(COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel)

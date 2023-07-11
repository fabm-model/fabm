target_sources(fabm PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/fabm_c.F90
            ${CMAKE_CURRENT_LIST_DIR}/variable.F90
            ${CMAKE_CURRENT_LIST_DIR}/parameter.F90
            ${CMAKE_CURRENT_LIST_DIR}/link_list.F90
            ${CMAKE_CURRENT_LIST_DIR}/c_helper.F90
            ${CMAKE_CURRENT_LIST_DIR}/helper.F90
            ${CMAKE_CURRENT_LIST_DIR}/integrate.F90
           )

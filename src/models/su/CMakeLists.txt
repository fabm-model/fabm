if(TARGET fabm_models_ersem)
   set(ERSEM_FILES ersem/ersem_model_library.F90 ersem/mesozooplankton.F90)
   add_definitions(-DSU_ERSEM)
endif()

add_library(fabm_models_su OBJECT
            su_model_library.F90
            npz.F90
            mixo.F90
            algae.F90
            bacteria.F90
            pools.F90
            light_cycle.F90
            ${ERSEM_FILES}
           )

add_dependencies(fabm_models_su fabm_base)
source_group(ersem REGULAR_EXPRESSION ersem/)

if(TARGET fabm_models_ersem)
   add_dependencies(fabm_models_su fabm_models_ersem)
   target_include_directories(fabm_models_su PRIVATE $<TARGET_PROPERTY:fabm_models_ersem,INTERFACE_INCLUDE_DIRECTORIES>)
endif()

#-----------------------------------------------------------------------
# Copyright (C) 2020 - Suzana Leles
#-----------------------------------------------------------------------

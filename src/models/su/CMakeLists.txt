add_library(fabm_models_su OBJECT
            su_model_library.F90
            npz.F90
            mixo.F90
            algae.F90
            bacteria.F90
            pools.F90
            light_cycle.F90
           )

target_link_libraries(fabm_models_su PRIVATE fabm_base)
source_group(ersem REGULAR_EXPRESSION ersem/)

if(TARGET fabm_models_ersem)
   target_link_libraries(fabm_models_su PRIVATE fabm_models_ersem)
   target_sources(fabm_models_su PRIVATE ersem/ersem_model_library.F90 ersem/mesozooplankton.F90)
   target_compile_definitions(fabm_models_su PRIVATE SU_ERSEM)
endif()

#-----------------------------------------------------------------------
# Copyright (C) 2020 - Suzana Leles
#-----------------------------------------------------------------------

add_library(fabm_models_selma OBJECT
            selma_model_library.F90
            selma.F90
            phytoplankton.F90
            zooplankton.F90
           )

target_link_libraries(fabm_models_selma PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2017 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

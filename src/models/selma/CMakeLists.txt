add_library(fabm_models_selma OBJECT
            selma_model_library.F90
            selma.F90
            phytoplankton.F90
            zooplankton.F90
           )

add_dependencies(fabm_models_selma fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2017 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

add_library(fabm_models_au OBJECT
            pclake/au_abiotic_sediment.F90
            pclake/au_abiotic_water.F90
            pclake/au_auxilary.F90
            pclake/au_dummy_state.F90
            pclake/au_foodweb_sediment.F90
            pclake/au_foodweb_water.F90
            pclake/au_macrophytes.F90
            pclake/au_pclake_model_library.F90
            pclake/au_phytoplankton_sediment.F90
            pclake/au_phytoplankton_water.F90
            pclake/au_utility.F90
            prey_predator/prey_predator.F90
           )

add_dependencies(fabm_models_au fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

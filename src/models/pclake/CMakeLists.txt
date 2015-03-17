add_library(fabm_models_pclake OBJECT
            utility.F90
            abiotic_sediment.F90
            abiotic_water.F90
            phytoplankton_sediment.F90
            phytoplankton_water.F90
            macrophytes.F90
            foodweb_sediment.F90
            foodweb_water.F90
            auxilary.F90
            model_library.F90
           )

add_dependencies(fabm_models_pclake fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

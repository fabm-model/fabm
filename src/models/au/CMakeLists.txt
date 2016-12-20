add_library(${PROJECT_PREFIX}_models_au OBJECT
            model_library.F90
            pclake/model_library.F90
            pclake/abiotic_sediment.F90
            pclake/abiotic_water.F90
            pclake/auxiliary.F90
            pclake/zoobenthos.F90
            pclake/zooplankton.F90
            pclake/fish.F90
            pclake/macrophytes.F90
            pclake/phytoplankton_sediment.F90
            pclake/phytoplankton_water.F90
            pclake/utility.F90
            prey_predator/prey_predator.F90
            prey_predator/lotka_volterra.F90
            prey_predator/jacob_monod.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_au ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

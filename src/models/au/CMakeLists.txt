add_library(fabm_models_au OBJECT
            model_library.F90
            prey_predator/prey_predator.F90
            prey_predator/lotka_volterra.F90
            prey_predator/jacob_monod.F90
           )

add_dependencies(fabm_models_au fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

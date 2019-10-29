add_library(fabm_models_gotm OBJECT
            gotm_model_library.F90
            ergom.F90
            fasham.F90
            npzd.F90
            light.F90
           )

add_dependencies(fabm_models_gotm fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

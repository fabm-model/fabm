add_library(fabm_models_gotm OBJECT
            gotm_model_library.F90
            ergom/ergom.F90
            fasham/fasham.F90
            npzd/npzd.F90
            light/light.F90
           )

add_dependencies(fabm_models_gotm fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

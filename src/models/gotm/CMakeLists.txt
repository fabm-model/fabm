add_library(${PROJECT_PREFIX}_models_gotm OBJECT
            gotm_model_library.F90
            ergom/ergom.F90
            fasham/fasham.F90
            npzd/npzd.F90
            light/light.F90
            photo/photo.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_gotm ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

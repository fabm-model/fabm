add_library(${PROJECT_PREFIX}_models_jrc OBJECT
            jrc_model_library.F90
            bsem/bsem.F90
            med_ergom/med_ergom.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_jrc ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2016 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

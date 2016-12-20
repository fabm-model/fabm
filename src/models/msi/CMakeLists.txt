add_library(${PROJECT_PREFIX}_models_msi OBJECT
            ergom1/ergom1.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_msi ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

add_library(${PROJECT_PREFIX}_models_metu OBJECT
            mnemiopsis/mnemiopsis.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_metu ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

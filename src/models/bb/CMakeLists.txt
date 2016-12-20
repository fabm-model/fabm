add_library(${PROJECT_PREFIX}_models_bb OBJECT
            bb_model_library.F90
            filter_feeder/filter_feeder.F90
            passive/passive.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_bb ${PROJECT_PREFIX}_base)


#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

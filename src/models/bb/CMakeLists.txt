add_library(fabm_models_bb OBJECT
            bb_model_library.F90
            filter_feeder/filter_feeder.F90
            lorenz63/lorenz63.F90
            passive/passive.F90
            gas/gas.F90
           )

target_link_libraries(fabm_models_bb PRIVATE fabm_base)


#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

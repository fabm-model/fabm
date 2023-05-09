add_library(fabm_models_metu OBJECT
            mnemiopsis/mnemiopsis.F90
           )

target_link_libraries(fabm_models_metu PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

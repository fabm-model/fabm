add_library(fabm_models_hzg OBJECT
            omexdia_p/omexdia_p.F90
           )

target_link_libraries(fabm_models_hzg PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

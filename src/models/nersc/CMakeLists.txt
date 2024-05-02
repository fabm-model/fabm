add_library(fabm_models_nersc OBJECT
            nersc_model_library.F90
            ecosmo/ecosmo_operational.F90
           )

target_link_libraries(fabm_models_nersc PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

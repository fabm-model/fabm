add_library(fabm_models_msi OBJECT
            msi_model_library.F90
            ergom1.F90
           )

target_link_libraries(fabm_models_msi PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

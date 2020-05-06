add_library(fabm_models_csiro OBJECT
            csiro_model_library.F90
            seagrass.F90
           )

add_dependencies(fabm_models_csiro fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

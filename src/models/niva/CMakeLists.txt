add_library(fabm_models_niva OBJECT
            niva_model_library.F90
            oxydep/oxydep.F90
           )

add_dependencies(fabm_models_niva fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

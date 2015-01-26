add_library(fabm_models_pml OBJECT
            pml_model_library.F90
            carbonate/carbonate.F90
            carbonate/co2_dyn.F90
           )

add_dependencies(fabm_models_pml fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

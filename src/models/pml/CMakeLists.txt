add_library(${PROJECT_PREFIX}_models_pml OBJECT
            pml_model_library.F90
            carbonate/carbonate.F90
            carbonate/co2_dyn.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_pml ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

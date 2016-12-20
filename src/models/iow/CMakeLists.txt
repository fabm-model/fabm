add_library(${PROJECT_PREFIX}_models_iow OBJECT
           iow_model_library.F90
           age/iow_age.F90
           ergom/ergom_cgt.F90
           spm/spm.F90
           tnt/tnt.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_iow ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

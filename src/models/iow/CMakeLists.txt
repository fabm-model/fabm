add_library(fabm_models_iow OBJECT
           iow_model_library.F90
           age/iow_age.F90
           ergom/ergom_cgt.F90
           spm/spm.F90
           )

target_link_libraries(fabm_models_iow PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

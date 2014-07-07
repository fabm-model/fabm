add_library(iow OBJECT
           iow_model_library.F90
           age/iow_age.F90
           ergom/ergom_cgt.F90
           spm/spm.F90
           )

add_dependencies(iow fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

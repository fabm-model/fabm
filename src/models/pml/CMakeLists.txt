add_library(pml OBJECT
            carbonate/carbonate.F90
            carbonate/co2_dyn.F90
           )

add_dependencies(pml fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

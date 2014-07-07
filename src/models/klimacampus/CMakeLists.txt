add_library(klimacampus OBJECT
            phy_feedback/phy_feedback.F90
           )

add_dependencies(klimacampus fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

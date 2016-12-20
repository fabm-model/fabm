add_library(${PROJECT_PREFIX}_models_klimacampus OBJECT
            phy_feedback/phy_feedback.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_klimacampus ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

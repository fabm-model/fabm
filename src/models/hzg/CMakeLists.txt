add_library(${PROJECT_PREFIX}_models_hzg OBJECT
            omexdia_p/omexdia_p.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_hzg ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

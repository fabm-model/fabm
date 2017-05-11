add_library(fabm_models_examples OBJECT
            examples_model_library.F90
            benthic_predator.F90
            duplicator.F90
            mean.F90
            npzd/examples_npzd_model_library.F90
            npzd/det.F90
            npzd/nut.F90
            npzd/phy.F90
            npzd/zoo.F90
            nonlocal.F90
            light_cycle.F90
           )

source_group(npzd REGULAR_EXPRESSION npzd/)

add_dependencies(fabm_models_examples fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

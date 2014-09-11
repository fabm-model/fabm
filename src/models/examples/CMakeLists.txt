add_library(fabm_models_examples OBJECT
            examples_model_library.F90
            benthic_predator/benthic_predator.F90
            duplicator/duplicator.F90
            mean/mean.F90
            npzd/examples_npzd.F90
            npzd/det.F90
            npzd/nut.F90
            npzd/phy.F90
            npzd/zoo.F90
            npzd_f2003/npzd_f2003.F90
           )

add_dependencies(fabm_models_examples fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

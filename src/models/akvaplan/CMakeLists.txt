add_library(fabm_models_akvaplan OBJECT
            akvaplan_model_library.F90
            tracer.F90
            plume_injection.F90
            tracer_sed.F90
            antiparasitic.F90
           )

add_dependencies(fabm_models_akvaplan fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2016 - Akvaplan-niva
#-----------------------------------------------------------------------

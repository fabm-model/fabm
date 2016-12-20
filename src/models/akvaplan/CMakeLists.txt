add_library(${PROJECT_PREFIX}_models_akvaplan OBJECT
            akvaplan_model_library.F90
            tracer.F90
            plume_injection.F90
           )

add_dependencies(${PROJECT_PREFIX}_models_akvaplan ${PROJECT_PREFIX}_base)

#-----------------------------------------------------------------------
# Copyright (C) 2016 - Akvaplan-niva
#-----------------------------------------------------------------------

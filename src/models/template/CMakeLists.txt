add_library(fabm_models_template OBJECT
            model_library.F90
            mymodel.F90
           )

target_link_libraries(fabm_models_template PRIVATE fabm_base)

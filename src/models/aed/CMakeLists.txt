add_library(fabm_models_aed OBJECT
            aed.h
            aed_carbon.F90
            aed_chlorophylla.F90
            aed_core.F90
            aed_iron.F90
            aed_models.F90
            aed_nitrogen.F90
            aed_organic_matter.F90
            aed_oxygen.F90
            aed_pathogens.F90
            aed_phosphorus.F90
            aed_phytoplankton.F90
            aed_phyto_utils.F90
            aed_sedflux.F90
            aed_silica.F90
            aed_sulfur.F90
            aed_totals.F90
            aed_tracer.F90
            aed_util.F90
            aed_zooplankton.F90
            aed_zoop_utils.F90
           )

add_dependencies(fabm_models_aed fabm_base)

# Make sure aed.h can be found.
set_property(TARGET fabm_models_aed APPEND PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

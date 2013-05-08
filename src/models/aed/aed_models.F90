!###############################################################################
!#                                                                             #
!# aed_models.F90                                                              #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created June 2012                                                           #
!#                                                                             #
!###############################################################################

#ifdef _FABM_F2003_

#include "aed.h"

MODULE aed_models
!-------------------------------------------------------------------------------
! A wrapper module for initialising all aed modules
!-------------------------------------------------------------------------------
   USE fabm_types

   USE aed_sedflux
   USE aed_chlorophylla
   USE aed_oxygen
   USE aed_silica
   USE aed_carbon
   USE aed_nitrogen
   USE aed_phosphorus
   USE aed_organic_matter
   USE aed_phytoplankton
   USE aed_pathogens
   USE aed_zooplankton
   USE aed_iron
   USE aed_sulfur
   USE aed_tracer
!  USE aed_geochemistry
!  USE aed_seddiagenesis
   USE aed_bacteria
   USE aed_viruses
   USE aed_totals

   IMPLICIT NONE

   PRIVATE   ! By default make everything private

   PUBLIC aed_create_model

CONTAINS
!===============================================================================


!###############################################################################
FUNCTION aed_create_model(namlst,modelname,instancename,parent) RESULT(model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CHARACTER(*),INTENT(in) :: modelname,instancename
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_model_info),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN
    NULLIFY(model)

    ! print *,'**** Initialising model ',name
    SELECT case (modelname)
       case ('aed_sedflux');        model => aed_sedflux_create(namlst,instancename,parent);
       case ('aed_chlorophylla');   model => aed_chla_create(namlst,instancename,parent);
       case ('aed_oxygen');         model => aed_oxygen_create(namlst,instancename,parent);
       case ('aed_silica');         model => aed_silica_create(namlst,instancename,parent);
       case ('aed_carbon');         model => aed_carbon_create(namlst,instancename,parent);
       case ('aed_nitrogen');       model => aed_nitrogen_create(namlst,instancename,parent);
       case ('aed_phosphorus');     model => aed_phosphorus_create(namlst,instancename,parent);
       case ('aed_organic_matter'); model => aed_organic_matter_create(namlst,instancename,parent);
       case ('aed_phytoplankton');  model => aed_phytoplankton_create(namlst,instancename,parent);
       case ('aed_zooplankton');    model => aed_zooplankton_create(namlst,instancename,parent);
       case ('aed_pathogens');      model => aed_pathogens_create(namlst,instancename,parent);
       case ('aed_iron');           model => aed_iron_create(namlst,instancename,parent);
       case ('aed_sulfur');         model => aed_sulfur_create(namlst,instancename,parent);
       case ('aed_tracer');         model => aed_tracer_create(namlst,instancename,parent);
!      case ('aed_geochemistry');   model => aed_geochemistry_create(namlst,instancename,parent);
!      case ('aed_seddiagenesis');  model => aed_seddiagenesis_create(namlst,instancename,parent);
       case ('aed_bacteria');       model => aed_bacteria_create(namlst,instancename,parent);
       case ('aed_viruses');        model => aed_viruses_create(namlst,instancename,parent);
       case ('aed_totals');         model => aed_totals_create(namlst,instancename,parent);
       case default
           print *,'*** Unknown module ',modelname
    END SELECT

!   IF (ASSOCIATED(model)) CALL initialize_model_info(model)

END FUNCTION aed_create_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE aed_models
#endif

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

#include "aed.h"

MODULE aed_models
!-------------------------------------------------------------------------------
! A wrapper module for initialising all aed modules
!-------------------------------------------------------------------------------
   USE aed_core

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
   USE aed_totals

   IMPLICIT NONE

   PRIVATE   ! By default make everything private

   TYPE,EXTENDS(type_base_model_factory) :: type_factory
      CONTAINS
      PROCEDURE :: create
   END TYPE

   TYPE (type_factory),SAVE,TARGET,PUBLIC :: aed_model_factory

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE create(self,name,model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_factory),INTENT(in) :: self
   CHARACTER(*),        INTENT(in) :: name
   CLASS (type_base_model),POINTER :: model
!-------------------------------------------------------------------------------
!BEGIN

    ! print *,'**** Initialising model ', modelname
    SELECT case (name)
       CASE ('aed_sedflux');        ALLOCATE(aed_type_sedflux::model)
       CASE ('aed_chlorophylla');   ALLOCATE(aed_type_chla::model)
       CASE ('aed_oxygen');         ALLOCATE(aed_type_oxygen::model)
       CASE ('aed_silica');         ALLOCATE(aed_type_silica::model)
       CASE ('aed_carbon');         ALLOCATE(aed_type_carbon::model)
       CASE ('aed_nitrogen');       ALLOCATE(aed_type_nitrogen::model)
       CASE ('aed_phosphorus');     ALLOCATE(aed_type_phosphorus::model)
       CASE ('aed_organic_matter'); ALLOCATE(aed_type_organic_matter::model)
       CASE ('aed_phytoplankton');  ALLOCATE(aed_type_phytoplankton::model)
       CASE ('aed_zooplankton');    ALLOCATE(aed_type_zooplankton::model)
       CASE ('aed_pathogens');      ALLOCATE(aed_type_pathogens::model)
       CASE ('aed_iron');           ALLOCATE(aed_type_iron::model)
       CASE ('aed_sulfur');         ALLOCATE(aed_type_sulfur::model)
       CASE ('aed_tracer');         ALLOCATE(aed_type_tracer::model)
       CASE ('aed_totals');         ALLOCATE(aed_type_totals::model)
       CASE DEFAULT ;               NULLIFY(model)
    END SELECT

END SUBROUTINE create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE aed_models

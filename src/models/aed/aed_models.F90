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

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: aed_model_factory

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

    SELECT case (name)
       case ('aed_sedflux');        allocate(type_aed_sedflux::model)
       case ('aed_chlorophylla');   allocate(type_aed_chla::model)
       case ('aed_oxygen');         allocate(type_aed_oxygen::model)
       case ('aed_silica');         allocate(type_aed_silica::model)
       case ('aed_carbon');         allocate(type_aed_carbon::model)
       case ('aed_nitrogen');       allocate(type_aed_nitrogen::model)
       case ('aed_phosphorus');     allocate(type_aed_phosphorus::model)
       case ('aed_organic_matter'); allocate(type_aed_organic_matter::model)
       case ('aed_phytoplankton');  allocate(type_aed_phytoplankton::model)
       case ('aed_zooplankton');    allocate(type_aed_zooplankton::model)
       case ('aed_pathogens');      allocate(type_aed_pathogens::model)
       case ('aed_iron');           allocate(type_aed_iron::model)
       case ('aed_sulfur');         allocate(type_aed_sulfur::model)
       case ('aed_tracer');         allocate(type_aed_tracer::model)
!      case ('aed_geochemistry');   allocate(type_aed_geochemistry::model)
!      case ('aed_seddiagenesis');  allocate(type_aed_seddiagenesis::model)
       case ('aed_bacteria');       allocate(type_aed_bacteria::model)
       case ('aed_viruses');        allocate(type_aed_viruses::model)
       case ('aed_totals');         allocate(type_aed_totals::model)
    END SELECT

END SUBROUTINE create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
END MODULE aed_models

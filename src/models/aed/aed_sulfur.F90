!###############################################################################
!#                                                                             #
!# aed_sulfur.F90                                                              #
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
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_sulfur
!-------------------------------------------------------------------------------
! aed_sulfur --- sulfur biogeochemical model
!
! The AED module sulfur contains equations that describe exchange of
! soluable reactive sulfur across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_type_sulfur
!
   TYPE,extends(type_base_model) :: aed_type_sulfur
!     Variable identifiers
      TYPE (type_state_variable_id)      :: id_so4
      TYPE (type_dependency_id)          :: id_temp
      TYPE (type_diagnostic_variable_id) :: id_sed_so4

!     Model parameters
      AED_REAL :: Fsed_dic,Ksed_dic,theta_sed_dic
      LOGICAL  :: use_oxy,use_dic

      CONTAINS     ! Model Methods
        PROCEDURE :: initialize               => aed_init_sulfur
        PROCEDURE :: do                       => aed_sulfur_do
        PROCEDURE :: do_ppdd                  => aed_sulfur_do_ppdd
        PROCEDURE :: do_benthos               => aed_sulfur_do_benthos
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_sulfur(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_sulfur),TARGET,INTENT(inout) :: self
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status
   INTEGER           :: num_sulfurs
   AED_REAL          :: decay(100)
   AED_REAL          :: settling(100)
   AED_REAL          :: Fsed(100)

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_sulfur/ num_sulfurs,decay,settling,Fsed
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_sulfur,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_sulfur'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   PRINT *,'AED_SULFUR : Note this module has not been completed. Stopping.'


!  ! Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature)

END SUBROUTINE aed_init_sulfur
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sulfur_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_sulfur model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_sulfur),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
!  AED_REAL           :: dic,diff_dic
!  AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
!  _GET_(self%id_dic,dic) ! sulfur

   ! Set temporal derivatives
!  diff_dic = 0.

!  _SET_ODE_(self%id_dic,diff_dic)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_sulfur_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sulfur_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of sulfur biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_sulfur),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL                   :: dic
   AED_REAL                   :: diff_dic
   AED_REAL, parameter        :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
!  _GET_(self%id_dic,dic) ! sulfur

   ! Set temporal derivatives
!  diff_dic = 0.

!  _SET_PP_(self%id_dic,self%id_dic,diff_dic)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_sulfur_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sulfur_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED sulfur.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_sulfur),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
!  AED_REAL :: dic,oxy

   ! Temporary variables
!  AED_REAL :: dic_flux

   ! Parameters
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
!  _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
!  _GET_(self%id_dic,dic) ! sulfur

!  IF (self%use_oxy) THEN
!     ! Sediment flux dependent on oxygen and temperature
!     _GET_(self%id_oxy,oxy)
!     dic_flux = self%Fsed_dic * self%Ksed_dic/(self%Ksed_dic+oxy) * (self%theta_sed_dic**(temp-20.0))
!  ELSE
!     ! Sediment flux dependent on temperature only.
!     dic_flux = self%Fsed_dic * (self%theta_sed_dic**(temp-20.0))
!  ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to FABM.
   !_SET_BOTTOM_FLUX_(self%id_dic,dic_flux/secs_pr_day)
   !_SET_ODE_SED_FLUX_(self%id_dic,dic_flux)
!  _SET_BOTTOM_EXCHANGE_(self%id_dic,dic_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_dic,-dic_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
!  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_dic,dic_flux)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

END SUBROUTINE aed_sulfur_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_sulfur

!###############################################################################
!#                                                                             #
!# aed_silica.F90                                                              #
!#                                                                             #
!# Created 24 August 2011                                                      #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!###############################################################################

#include "aed.h"

MODULE aed_silica
!-------------------------------------------------------------------------------
! Nitrogen module contains equations for nitrification and deitrification
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE   ! By default make everything private.
!
   PUBLIC aed_type_silica
!
   TYPE,extends(type_base_model) :: aed_type_silica
!     Variable identifiers
      TYPE (type_state_variable_id)      :: id_rsi,id_oxy
      TYPE (type_horizontal_dependency_id)  :: id_Fsed_rsi
      TYPE (type_dependency_id)          :: id_temp
      TYPE (type_horizontal_diagnostic_variable_id) :: id_sed_rsi

!     Model parameters
      AED_REAL :: Fsed_rsi,Ksed_rsi,theta_sed_rsi
      LOGICAL  :: use_oxy,use_rsi,use_sed_model

      CONTAINS      ! Model Methods
        PROCEDURE :: initialize               => aed_init_silica
        PROCEDURE :: do                       => aed_silica_do
        PROCEDURE :: do_ppdd                  => aed_silica_do_ppdd
        PROCEDURE :: do_benthos               => aed_silica_do_benthos
   END TYPE

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_silica(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: namlst
   CLASS (aed_type_silica),TARGET,INTENT(inout) :: self

!
!LOCALS
   INTEGER  :: status

   AED_REAL          :: rsi_initial=4.5
   AED_REAL          :: Fsed_rsi = 3.5
   AED_REAL          :: Ksed_rsi = 30.0
   AED_REAL          :: theta_sed_rsi = 1.0
   CHARACTER(len=64) :: silica_reactant_variable=''
   CHARACTER(len=64) :: Fsed_rsi_variable=''

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_silica/ rsi_initial,Fsed_rsi,Ksed_rsi,theta_sed_rsi,silica_reactant_variable, &
                         Fsed_rsi_variable
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_silica,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_silica'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%Fsed_rsi  = Fsed_rsi/secs_pr_day
   self%Ksed_rsi  = Ksed_rsi
   self%theta_sed_rsi = theta_sed_rsi
   self%use_oxy = silica_reactant_variable .NE. '' !This means oxygen module switched on

   ! Register state variables
   CALL self%register_state_variable(self%id_rsi,'rsi','mmol/m**3', 'silica',     &
                                    rsi_initial,minimum=zero_,no_river_dilution=.false.)

   ! Register external state variable dependencies
   IF (self%use_oxy) &
      CALL self%register_state_dependency(self%id_oxy,silica_reactant_variable)

   self%use_sed_model = Fsed_rsi_variable .NE. ''
   IF (self%use_sed_model) THEN
      CALL self%register_horizontal_dependency(self%id_Fsed_rsi,Fsed_rsi_variable)
      CALL self%request_coupling(self%id_Fsed_rsi,Fsed_rsi_variable)
   ENDIF

   ! Register diagnostic variables
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_rsi,'sed_rsi','mmol/m**2/d', &
                     'reactive silica')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature)
END SUBROUTINE aed_init_silica
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_silica_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_silica model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_silica),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!ARGUMENT
   !AED_REAL                   :: rsi,oxy,temp,tss !State variables
   AED_REAL, parameter        :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_


   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_silica_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_silica_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of silica biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_silica),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!ARGUMENT
   AED_REAL, parameter        :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_


   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_silica_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_silica_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED silica.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_silica),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: rsi,oxy

   ! Temporary variables
   AED_REAL :: rsi_flux, Fsed_rsi

   ! Parameters
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_(self%id_rsi,rsi) ! silica

   IF (self%use_sed_model) THEN
       _GET_HORIZONTAL_(self%id_Fsed_rsi,Fsed_rsi)
   ELSE
       Fsed_rsi = self%Fsed_rsi
   ENDIF

   IF (self%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
       _GET_(self%id_oxy,oxy)
       rsi_flux = Fsed_rsi * self%Ksed_rsi/(self%Ksed_rsi+oxy) * (self%theta_sed_rsi**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
       rsi_flux = Fsed_rsi * (self%theta_sed_rsi**(temp-20.0))
   ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to FABM.
   _SET_BOTTOM_EXCHANGE_(self%id_rsi,rsi_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_rsi,-rsi_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_rsi,rsi_flux)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

END SUBROUTINE aed_silica_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_silica

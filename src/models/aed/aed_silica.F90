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

#ifdef _FABM_F2003_

#include "aed.h"

MODULE aed_silica
!-------------------------------------------------------------------------------
! Nitrogen module contains equations for nitrification and deitrification
!-------------------------------------------------------------------------------
   USE fabm_types
   USE fabm_driver

   IMPLICIT NONE

   PRIVATE   ! By default make everything private.
!
   PUBLIC type_aed_silica, aed_silica_create
!
   TYPE,extends(type_base_model) :: type_aed_silica
!     Variable identifiers
      type (type_state_variable_id)      :: id_rsi,id_oxy
      type (type_bottom_state_variable_id)  :: id_Fsed_rsi
      type (type_dependency_id)          :: id_temp
      type (type_horizontal_diagnostic_variable_id) :: id_sed_rsi

!     Model parameters
      real(rk) :: Fsed_rsi,Ksed_rsi,theta_sed_rsi
      LOGICAL  :: use_oxy,use_rsi,use_sed_model

      CONTAINS      ! Model Methods
!       procedure :: initialize               => aed_silica_init
        procedure :: do                       => aed_silica_do
        procedure :: do_ppdd                  => aed_silica_do_ppdd
        procedure :: do_benthos               => aed_silica_do_benthos
        procedure :: get_conserved_quantities => aed_silica_get_conserved_quantities
   END TYPE

!===============================================================================
CONTAINS


!###############################################################################
FUNCTION aed_silica_create(namlst,name,parent) RESULT(self)
!-------------------------------------------------------------------------------
! Initialise the AED model
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: namlst
   CHARACTER(len=*),INTENT(in)              :: name
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_aed_silica),POINTER :: self

   real(rk)          :: rsi_initial=4.5
   real(rk)          :: Fsed_rsi = 3.5
   real(rk)          :: Ksed_rsi = 30.0
   real(rk)          :: theta_sed_rsi = 1.0
   CHARACTER(len=64) :: silica_reactant_variable=''
   CHARACTER(len=64) :: Fsed_rsi_variable=''

   real(rk),PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_silica/ rsi_initial,Fsed_rsi,Ksed_rsi,theta_sed_rsi,silica_reactant_variable, &
                         Fsed_rsi_variable
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL initialize_model_info(self,name,parent)

   ! Read the namelist
   read(namlst,nml=aed_silica,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%Fsed_rsi  = Fsed_rsi/secs_pr_day
   self%Ksed_rsi  = Ksed_rsi
   self%theta_sed_rsi = theta_sed_rsi

   ! Register state variables
   call self%register_state_variable(self%id_rsi,'rsi','mmol/m**3', 'silica',     &
                                    rsi_initial,minimum=0.0_rk,no_river_dilution=.false.)

   ! Register external state variable dependencies
   self%use_oxy = silica_reactant_variable .NE. '' !This means oxygen module switched on
   IF (self%use_oxy) &
      call self%register_state_dependency(self%id_oxy,silica_reactant_variable)

   self%use_sed_model = Fsed_rsi_variable .NE. ''
   IF (self%use_sed_model) &
      call self%register_bottom_state_dependency(self%id_Fsed_rsi,Fsed_rsi_variable)

   ! Register diagnostic variables
   call self%register_horizontal_diagnostic_variable(self%id_sed_rsi,'sed_rsi','mmol/m**2/d', &
                     'reactive silica',                                              &
                     time_treatment=time_treatment_step_integrated)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   RETURN

99 CALL fatal_error('aed_silica_init','Error reading namelist aed_silica')

END FUNCTION aed_silica_create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_silica_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_silica model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_silica),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!ARGUMENT
   !real(rk)                   :: rsi,oxy,temp,tss !State variables
   real(rk), parameter        :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_


   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_silica_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_silica_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of silica biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_silica),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!ARGUMENT
   real(rk), parameter        :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_


   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_silica_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_silica_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED silica.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_silica),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   ! Environment
   real(rk) :: temp

   ! State
   real(rk) :: rsi,oxy

   ! Temporary variables
   real(rk) :: rsi_flux, Fsed_rsi

   ! Parameters
   real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

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
   _FABM_HORIZONTAL_LOOP_END_

END SUBROUTINE aed_silica_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_silica_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities (currently only silica)
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_silica),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   real(rk) :: rsi
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_rsi,rsi) ! silica

   ! Total nutrient is simply the sum of all variables.
!   _SET_CONSERVED_QUANTITY_(self%id_totP,rsi)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_silica_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_silica
#endif

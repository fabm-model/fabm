!###############################################################################
!#                                                                             #
!# aed_tracer.F90                                                              #
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
MODULE aed_tracer
!-------------------------------------------------------------------------------
! aed_tracer --- tracer biogeochemical model
!
! The AED module tracer contains equations that describe exchange of
! soluable reactive tracer across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_type_tracer
!
   TYPE,extends(type_base_model) :: aed_type_tracer
!     Variable identifiers
      TYPE (type_state_variable_id),ALLOCATABLE :: id_ss(:)
      TYPE (type_state_variable_id)             :: id_retain
      TYPE (type_dependency_id)                 :: id_temp
      TYPE (type_horizontal_dependency_id)      :: id_taub
      LOGICAL                                   :: retention_time
      LOGICAL                                   :: resuspension

!     Model parameters
      AED_REAL,ALLOCATABLE :: decay(:), settling(:), Fsed(:)
      AED_REAL,ALLOCATABLE :: epsilon(:), tau_0(:), tau_r(:), Ke_ss(:)

      CONTAINS      ! Model Methods
        PROCEDURE :: initialize               => aed_init_tracer
        PROCEDURE :: do                       => aed_tracer_do
!       PROCEDURE :: do_ppdd                  => aed_tracer_do_ppdd
        PROCEDURE :: do_benthos               => aed_tracer_do_benthos
        PROCEDURE :: get_light_extinction     => aed_tracer_get_light_extinction
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_tracer(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_type_tracer),TARGET,INTENT(inout) :: self
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: num_tracers
   AED_REAL :: decay(100)
   AED_REAL :: settling(100)
   AED_REAL :: Fsed(100)
   AED_REAL :: epsilon(100)
   AED_REAL :: tau_0(100)
   AED_REAL :: tau_r(100)
   AED_REAL :: Ke_ss(100)
   AED_REAL :: trace_initial = zero_
   INTEGER  :: i
   LOGICAL  :: retention_time = .FALSE.
   LOGICAL  :: resuspension = .TRUE.
   CHARACTER(4) :: trac_name

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_tracer/ num_tracers,decay,settling,Fsed,resuspension,epsilon,tau_0,tau_r,Ke_ss,retention_time
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_tracer,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_tracer'

   ! Store parameter values in our own derived type

   IF ( num_tracers > 0 ) THEN
      ALLOCATE(self%id_ss(num_tracers))
      ALLOCATE(self%decay(num_tracers))    ; self%decay(1:num_tracers)    = decay(1:num_tracers)
      ALLOCATE(self%settling(num_tracers)) ; self%settling(1:num_tracers) = settling(1:num_tracers)
      ALLOCATE(self%Fsed(num_tracers))     ; self%Fsed(1:num_tracers)     = Fsed(1:num_tracers)

      ALLOCATE(self%epsilon(num_tracers))  ; self%epsilon(1:num_tracers)  = epsilon(1:num_tracers)
      ALLOCATE(self%tau_0(num_tracers))    ; self%tau_0(1:num_tracers)    = tau_0(1:num_tracers)
      ALLOCATE(self%tau_r(num_tracers))    ; self%tau_r(1:num_tracers)    = tau_r(1:num_tracers)
      ALLOCATE(self%Ke_ss(num_tracers))    ; self%Ke_ss(1:num_tracers)    = Ke_ss(1:num_tracers)

      trac_name = 'ss0'
      ! Register state variables
      DO i=1,num_tracers
         trac_name(3:3) = CHAR(ICHAR('0') + i)
                                             ! divide settling by secs_pr_day to convert m/d to m/s
         CALL self%register_state_variable(self%id_ss(i),TRIM(trac_name),'mmol/m**3','tracer', &
                                   trace_initial,minimum=zero_,no_river_dilution=.false.,vertical_movement=(settling(i)/secs_pr_day))
      ENDDO
      self%retention_time = retention_time
   ENDIF

   IF (retention_time) THEN
      CALL self%register_state_variable(self%id_retain, "ret",'sec','tracer', &
                                   trace_initial,minimum=zero_)
   ENDIF

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature)
   CALL self%register_dependency(self%id_taub,standard_variables%bottom_stress)

   self%resuspension = resuspension
END SUBROUTINE aed_init_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_tracer_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_tracer),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_

!
!-------------------------------------------------------------------------------
!BEGIN
   IF ( .NOT. self%retention_time ) RETURN

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   _SET_ODE_(self%id_retain, 1.0)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_tracer_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_tracer_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED tracer.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_tracer),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: ss, bottom_stress = 0.

   ! Temporary variables
   AED_REAL :: ss_flux, theta_sed_ss = 1.0, resus_flux = 0.
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN

   IF ( .NOT. ALLOCATED(self%id_ss) ) RETURN

   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature
   _GET_HORIZONTAL_(self%id_taub,bottom_stress)
   bottom_stress = MIN(bottom_stress, 100.0_rk)

    DO i=1,ubound(self%id_ss,1)
    ! Retrieve current (local) state variable values.
       _GET_(self%id_ss(i),ss)

      IF ( self%resuspension ) THEN
         IF (bottom_stress > self%tau_0(i)) THEN
            resus_flux = self%epsilon(i) * ( bottom_stress - self%tau_0(i)) / self%tau_r(i)
         ELSE
            resus_flux = 0.
         ENDIF
      ENDIF

      ! Sediment flux dependent on temperature only.
      ss_flux = self%Fsed(i) * (theta_sed_ss**(temp-20.0))

      ! Transfer sediment flux value to model.
      _SET_BOTTOM_EXCHANGE_(self%id_ss(i), ss_flux + resus_flux)

   ENDDO

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

END SUBROUTINE aed_tracer_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_tracer_get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_tracer),INTENT(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
!LOCALS
   AED_REAL :: ss
   INTEGER  :: ss_i
!
!-----------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   DO ss_i=1,ubound(self%id_ss,1)
      ! Retrieve current (local) state variable values.
      _GET_(self%id_ss(ss_i), ss)

      ! Self-shading with explicit contribution from background tracer concentration.
      _SET_EXTINCTION_(self%Ke_ss(ss_i)*ss)
   ENDDO

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_tracer_get_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_tracer

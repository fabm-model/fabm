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
      TYPE (type_dependency_id)                 :: id_temp

!     Model parameters
      AED_REAL,ALLOCATABLE :: decay(:),settling(:), Fsed(:)

      CONTAINS      ! Model Methods
        PROCEDURE :: initialize               => aed_init_tracer
        PROCEDURE :: do                       => aed_tracer_do
        PROCEDURE :: do_ppdd                  => aed_tracer_do_ppdd
        PROCEDURE :: do_benthos               => aed_tracer_do_benthos
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
   AED_REAL :: trace_initial = zero_
   INTEGER  :: i
   CHARACTER(4) :: trac_name

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_tracer/ num_tracers,decay,settling,Fsed
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_tracer,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_tracer'

   ! Store parameter values in our own derived type

   ALLOCATE(self%id_ss(num_tracers))
   ALLOCATE(self%decay(num_tracers))    ; self%decay(1:num_tracers)    = decay(1:num_tracers)
   ALLOCATE(self%settling(num_tracers)) ; self%settling(1:num_tracers) = settling(1:num_tracers)
   ALLOCATE(self%Fsed(num_tracers))     ; self%Fsed(1:num_tracers)     = Fsed(1:num_tracers)

   trac_name = 'ss0'
   ! Register state variables
   DO i=1,num_tracers
      trac_name(3:3) = CHAR(ICHAR('0') + i)
      CALL self%register_state_variable(self%id_ss(i),TRIM(trac_name),'mmol/m**3','tracer', &
                                   trace_initial,minimum=zero_,no_river_dilution=.false.)
   ENDDO

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature)
END SUBROUTINE aed_init_tracer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_tracer_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_tracer model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_tracer),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_tracer_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_tracer_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of tracer biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_tracer),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_tracer_do_ppdd
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
   AED_REAL :: ss

   ! Temporary variables
   AED_REAL :: ss_flux, theta_sed_ss = 1.0
   INTEGER  :: i

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    DO i=1,ubound(self%id_ss,1)
    ! Retrieve current (local) state variable values.
       _GET_(self%id_ss(i),ss)

      ! Sediment flux dependent on temperature only.
      ss_flux = self%Fsed(i) * (theta_sed_ss**(temp-20.0))

      ! Transfer sediment flux value to FABM.
      _SET_BOTTOM_EXCHANGE_(self%id_ss(i),ss_flux)

   ENDDO

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

END SUBROUTINE aed_tracer_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_tracer

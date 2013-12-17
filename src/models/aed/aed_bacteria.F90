!###############################################################################
!#                                                                             #
!# aed_bacteria.F90                                                            #
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

#ifdef _FABM_F2003_

#include "aed.h"

!
MODULE aed_bacteria
!-------------------------------------------------------------------------------
! aed_bacteria --- bacteria biogeochemical model
!
! The AED module bacteria contains equations that describe exchange of
! dynamics of hetertrophic bacteria
!-------------------------------------------------------------------------------
   USE fabm_types

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_bacteria

   type,extends(type_base_model) :: type_aed_bacteria
!     Variable identifiers
      type (type_state_variable_id)      :: id_bact
      type (type_dependency_id)          :: id_temp
    ! type (type_diagnostic_variable_id) :: id_sed_bact
      type (type_conserved_quantity_id)  :: id_totB

!     Model parameters
      real(rk) :: growth,mortality

      CONTAINS
!     Model Procedures
        procedure :: initialize               => aed_bacteria_init
        procedure :: do                       => aed_bacteria_do
        procedure :: do_ppdd                  => aed_bacteria_do_ppdd
        procedure :: do_benthos               => aed_bacteria_do_benthos
        procedure :: get_conserved_quantities => aed_bacteria_get_conserved_quantities
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
subroutine aed_bacteria_init(self,configunit)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_aed_bacteria),TARGET,INTENT(INOUT) :: self
   INTEGER,INTENT(in)                             :: configunit
!
!LOCALS

   real(rk)          :: growth
   real(rk)          :: mortality

   real(rk),PARAMETER :: secs_pr_day = 86400.
   real(rk),PARAMETER :: bact_initial = 0.
   NAMELIST /aed_bacteria/ growth,mortality

!-------------------------------------------------------------------------------
!BEGIN
   print *,"WARNING! aed_bacteria model is currently under development"

   ! Read the namelist
   read(configunit,nml=aed_bacteria,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%growth = growth
   self%mortality = mortality

   ! Register state variables
   call self%register_state_variable(self%id_bact,'bact','mmol/m**3','bacteria',     &
                                    bact_initial,minimum=0.0_rk,no_river_dilution=.false.)

   ! Register diagnostic variables
!  self%id_sed_bact = self%register_diagnostic_variable('sed_bact','mmol/m**2/d',   &
!                    'Filterable reactive bacteria',                                     &
!                    time_treatment=time_treatment_step_integrated, shape=shape_hz)

   ! Register conserved quantities
   call self%register_conserved_quantity(self%id_totB,'TB','mmol/m**3','Total bacteria')

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   RETURN

99 CALL self%fatal_error('aed_bacteria_init','Error reading namelist aed_bacteria')

end subroutine aed_bacteria_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bacteria_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_bacteria model
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_bacteria),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   real(rk)           :: bact,diff_bact
   real(rk),PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_bact,bact) ! bacteria

   ! Set temporal derivatives
   diff_bact = 0.

   PRINT *,' WARNING: Bacteria module is currently inactive'

   _SET_ODE_(self%id_bact,diff_bact)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_bacteria_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_bacteria_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of bacteria biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_bacteria),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   real(rk)                   :: bact
   real(rk)                   :: diff_bact
   real(rk), parameter        :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_bact,bact) ! bacteria

   ! Set temporal derivatives
   diff_bact = 0.

   _SET_PP_(self%id_bact,self%id_bact,diff_bact)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_bacteria_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_bacteria_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED bacteria.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_bacteria),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   ! Environment
   real(rk) :: temp

   ! State
   real(rk) :: bact

   ! Temporary variables
   real(rk) :: bact_flux

   ! Parameters
   real(rk),PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_(self%id_bact,bact) ! bacteria

!  _SET_BOTTOM_EXCHANGE_(self%id_bact,bact_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_bact,-bact_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
!  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_bact,bact_flux)

   ! Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_

END SUBROUTINE aed_bacteria_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_bacteria_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities (currently only bacteria)
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_bacteria),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   real(rk) :: bact
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_bact,bact) ! bacteria

   ! Total nutrient is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totB,bact)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_bacteria_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed_bacteria
#endif

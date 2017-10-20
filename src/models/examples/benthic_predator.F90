#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: examples_benthic_predator --- test for benthic interfaces in FABM
!
! !INTERFACE:
   module examples_benthic_predator
!
! !DESCRIPTION:
! This is a very simple model for a benthic predator, grazing according to a
! Monod/Michaelis-Menten functional response on a pelagic prey, and
! respiring/dying according to a linear loss term. Variables for the prey
! (e.g., phytoplankon or zooplankton) and the target for the losses (typically
! a detrital or mineral pool) must be provided by coupling to an external model.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC TYPES:
   type,extends(type_base_model),public :: type_examples_benthic_predator
!     Variable identifiers
      type (type_bottom_state_variable_id) :: id_biomass
      type (type_state_variable_id)        :: id_prey,id_waste
      type (type_bottom_state_variable_id) :: id_bottom_prey,id_bottom_waste

!     Model parameters: maximum grazing rate, half-saturation prey density, loss rate
      real(rk) :: g_max,K,h
      logical  :: interact_with_pelagic
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the benthic predator model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples\_benthic\_predator namelist is read and te variables
!  exported by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_examples_benthic_predator),intent(inout),target :: self
   integer,                               intent(in)           :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                  :: pred_initial=0.01
   real(rk)                  :: g_max = 1., K=1., h=0.05
   character(len=64)         :: waste_target_variable='',prey_source_variable=''

   real(rk), parameter :: days_per_sec = 1.0_rk/86400
   namelist /examples_benthic_predator/ waste_target_variable,prey_source_variable,pred_initial,g_max,K,h
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   if (configunit>0) read(configunit,nml=examples_benthic_predator,err=99,end=100)

   ! Retrieve final authorative parameter values (e.g., from fabm.yaml)
   ! NB: all rates must be provided in values per day, and are converted here to values per second.
   call self%get_parameter(self%g_max,'g_max','d-1','maximum grazing rate',default=self%g_max,scale_factor=days_per_sec)
   call self%get_parameter(self%h,    'h',    'd-1','mortality',           default=self%h,    scale_factor=days_per_sec)
   call self%get_parameter(self%interact_with_pelagic,'interact_with_pelagic','','interact directly with pelagic',default=.true.)
   if (self%interact_with_pelagic) then
      ! Predator gets prey from the pelagic and puts waste into a bottom-bound variable.
      call self%get_parameter(self%K,'K','mmol m-3','prey half-saturation',default=self%K)
      call self%register_state_dependency(self%id_prey, 'prey', 'mmol m-3','prey')
      call self%register_state_dependency(self%id_waste,'waste','mmol m-3','waste')
      if (prey_source_variable /='') call self%request_coupling(self%id_prey, prey_source_variable)
      if (waste_target_variable/='') call self%request_coupling(self%id_waste,waste_target_variable)
   else
      ! Predator gets prey from the bottom and puts waste into a bottom-bound variable.
      call self%get_parameter(self%K,'K','mmol m-2','prey half-saturation',default=self%K)
      call self%register_state_dependency(self%id_bottom_prey, 'prey', 'mmol m-2','prey')
      call self%register_state_dependency(self%id_bottom_waste,'waste','mmol m-2','waste')
   end if

   ! Register state variable for predator biomass - it is benthic because its identifier
   ! self%id_pred is of type type_bottom_state_variable_id.
   call self%register_state_variable(self%id_biomass,'bm','mmol/m**2','predator density',pred_initial,minimum=0.0_rk)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_biomass)

   return

99 call self%fatal_error('initialize','Error reading namelist examples_benthic_predator')
100 call self%fatal_error('initialize','Namelist examples_benthic_predator was not found')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
! This routine calculates the benthic sink and source terms, as well as
! (matching) bottom fluxes for pelagic variables. Both have units mmol/m**2/s.
!
! !INPUT PARAMETERS:
   class (type_examples_benthic_predator),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk) :: prey,biomass,g
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

      ! Retrieve current predator density (bottom-bound state variable)
      _GET_HORIZONTAL_(self%id_biomass,biomass)

      if (self%interact_with_pelagic) then
         ! Retrieve current prey density (pelagic state variable)
         _GET_(self%id_prey,prey)
      else
         ! Retrieve current prey density (bottom-bound state variable)
         _GET_HORIZONTAL_(self%id_bottom_prey,prey)
      end if

      ! Calculate grazing rate
      g = self%g_max*biomass*prey/(prey+self%K)

      ! Set local temporal derivatives of benthic variables
      _SET_BOTTOM_ODE_(self%id_biomass,g-self%h*biomass)

      if (self%interact_with_pelagic) then
         ! Set bottom fluxes of pelagic variables (these mirror local benthic derivatives)
         _SET_BOTTOM_EXCHANGE_(self%id_prey,-g)
         _SET_BOTTOM_EXCHANGE_(self%id_waste,self%h*biomass)
      else
         _SET_BOTTOM_ODE_(self%id_bottom_prey,-g)
         _SET_BOTTOM_ODE_(self%id_bottom_waste,self%h*biomass)
      end if

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------

   end module examples_benthic_predator

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

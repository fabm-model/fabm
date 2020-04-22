#include "fabm_driver.h"

module examples_benthic_predator
   ! This is a very simple model for a benthic predator, grazing according to a
   ! Monod/Michaelis-Menten functional response on a pelagic prey, and
   ! respiring/dying according to a linear loss term. Variables for the prey
   ! (e.g., phytoplankon or zooplankton) and the target for the losses (typically
   ! a detrital or mineral pool) must be provided by coupling to an external model.

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_examples_benthic_predator
      ! Variables
      type (type_bottom_state_variable_id) :: id_biomass
      type (type_state_variable_id)        :: id_prey, id_waste
      type (type_bottom_state_variable_id) :: id_bottom_prey, id_bottom_waste

      ! Parameters: maximum grazing rate, half-saturation prey density, loss rate
      real(rk) :: g_max, K, h
      logical  :: interact_with_pelagic
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_examples_benthic_predator), intent(inout), target :: self
      integer,                                intent(in)            :: configunit

      real(rk), parameter :: days_per_sec = 1.0_rk/86400

      ! Retrieve parameter values
      ! NB: all rates must be provided in values per day, and are converted here to values per second.
      call self%get_parameter(self%interact_with_pelagic, 'interact_with_pelagic', '', 'food and waste are pelagic (not benthic)', default=.true.)
      call self%get_parameter(self%g_max, 'g_max', 'd-1', 'maximum population growth rate', default=1._rk, scale_factor=days_per_sec)
      if (self%interact_with_pelagic) then
         ! Predator gets prey from the pelagic and puts waste into a pelagic variable.
         call self%get_parameter(self%K,'K', 'mmol m-3', 'prey half-saturation', default=1._rk)
         call self%register_state_dependency(self%id_prey, 'prey', 'mmol m-3', 'prey')
         call self%register_state_dependency(self%id_waste, 'waste', 'mmol m-3', 'waste')
      else
         ! Predator gets prey from the bottom and puts waste into a benthic variable.
         call self%get_parameter(self%K, 'K', 'mmol m-2', 'prey half-saturation', default=1._rk)
         call self%register_state_dependency(self%id_bottom_prey, 'prey', 'mmol m-2', 'prey')
         call self%register_state_dependency(self%id_bottom_waste, 'waste', 'mmol m-2', 'waste')
      end if
      call self%get_parameter(self%h, 'h', 'd-1', 'mortality', default=0.05_rk, scale_factor=days_per_sec)

      ! Register state variable for predator biomass - it is benthic because its identifier is of type type_bottom_state_variable_id
      call self%register_state_variable(self%id_biomass, 'bm', 'mmol m-2', 'predator density', 0.01_rk, minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_biomass)
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_examples_benthic_predator),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: prey, biomass, g

      ! Enter spatial loops over the horizontal domain (if any).
      _BOTTOM_LOOP_BEGIN_

         ! Retrieve current predator density (bottom-bound state variable)
         _GET_BOTTOM_(self%id_biomass,biomass)

         if (self%interact_with_pelagic) then
            ! Retrieve current prey density (pelagic state variable)
            _GET_(self%id_prey, prey)
         else
            ! Retrieve current prey density (bottom-bound state variable)
            _GET_BOTTOM_(self%id_bottom_prey, prey)
         end if

         ! Calculate grazing rate (mmol m-2 s-1)
         g = self%g_max * biomass * prey / (prey + self%K)

         ! Set local temporal derivatives of benthic variables
         _ADD_BOTTOM_SOURCE_(self%id_biomass, g - self%h * biomass)

         if (self%interact_with_pelagic) then
            ! Set bottom fluxes of pelagic variables (these mirror local benthic derivatives and have unit mmol m-2 s-1)
            _ADD_BOTTOM_FLUX_(self%id_prey, -g)
            _ADD_BOTTOM_FLUX_(self%id_waste, self%h * biomass)
         else
            _ADD_BOTTOM_SOURCE_(self%id_bottom_prey, -g)
            _ADD_BOTTOM_SOURCE_(self%id_bottom_waste, self%h * biomass)
         end if

      ! Leave spatial loops over the horizontal domain (if any).
      _BOTTOM_LOOP_END_

   end subroutine do_bottom

end module examples_benthic_predator

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module examples_sediment_pom
   ! This is a simple one compartment model for particulate organic matter
   ! in the sediment.

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_examples_sediment_pom
      ! Variables
      type (type_bottom_state_variable_id) :: id_c
      type (type_state_variable_id)        :: id_mintarget

      ! Parameters
      real(rk) :: r
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_examples_sediment_pom), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400

      ! Retrieve parameter values
      ! All rates must be provided in values per day in fabm.yaml.
      ! They are converted here to values per second (scale_factor argument) to ensure
      ! that the sources and sinks calculated from these parameters will be in per second as FABM expects.
      call self%get_parameter(self%r, 'r', 'd-1', 'remineralization rate', default=0.003_rk, scale_factor=d_per_s)

      ! Register state variables
      call self%register_state_variable(self%id_c, 'c','mmol m-2', 'density', initial_value=0.0_rk)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c)

      ! Register dependencies on external state variables
      call self%register_state_dependency(self%id_mintarget, 'mineralisation_target', 'mmol m-3', 'pelagic sink for remineralized matter')

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_examples_sediment_pom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c

      ! Enter spatial loops over the horizontal domain (if any).
      _BOTTOM_LOOP_BEGIN_

         ! Retrieve current predator density (bottom-bound state variable)
         _GET_BOTTOM_(self%id_c, c)

         ! Local source terms
         _ADD_BOTTOM_SOURCE_(self%id_c, -self%r * c)

         ! Bottom fluxes of pelagic variables
         _ADD_BOTTOM_FLUX_(self%id_mintarget, self%r * c)

      ! Leave spatial loops over the horizontal domain (if any).
      _BOTTOM_LOOP_END_

   end subroutine do_bottom

end module examples_sediment_pom

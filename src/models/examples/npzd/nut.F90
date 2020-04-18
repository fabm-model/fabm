#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - nutrient component
! This is a general nutrient (passive non-sinking, non-floating tracer), characterized by
! an initial concentration only.

module examples_npzd_nut
   use fabm_types

   implicit none

   private

   type, extends(type_base_model),public :: type_examples_npzd_nut
      ! Variable identifiers
      type (type_state_variable_id) :: id_n
   contains
      procedure :: initialize
   end type

contains

   subroutine initialize(self, configunit)
      class (type_examples_npzd_nut), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      ! Register state variables
      call self%register_state_variable(self%id_n, 'c', 'mmol m-3', 'concentration', 1.0_rk, minimum=0.0_rk, no_river_dilution=.true.)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_n)
   end subroutine initialize

end module examples_npzd_nut

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

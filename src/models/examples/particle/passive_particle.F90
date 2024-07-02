#include "fabm_driver.h"

module examples_particle_passive_particle

   use fabm_particle
   use fabm_types

   implicit none

   private

   type, extends(type_particle_model), public :: type_passive_particle
      type (type_state_variable_id) :: id_c
   contains
      procedure :: initialize
   end type

contains

   subroutine initialize(self, configunit)
      class (type_passive_particle), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      call self%register_state_variable(self%id_c, 'c', 'mmol m-3', 'concentration')
   end subroutine

end module

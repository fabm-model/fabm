#include "fabm_driver.h"

module fabm_builtin_product
   use fabm_types
   use fabm_particle

   implicit none

   private

   type, extends(type_particle_model), public :: type_product
      type (type_dependency_id)          :: id_term1
      type (type_dependency_id)          :: id_term2
      type (type_diagnostic_variable_id) :: id_result
      integer :: output = output_instantaneous
   contains
      procedure :: initialize     => product_initialize
      procedure :: after_coupling => product_after_coupling
      procedure :: do             => product_do
   end type type_product

contains

   subroutine product_initialize(self, configunit)
      class (type_product), intent(inout), target :: self
      integer,              intent(in)            :: configunit

      call self%register_dependency(self%id_term1, 'term1', '', 'term1')
      call self%register_dependency(self%id_term2, 'term2', '', 'term2')
      call self%register_diagnostic_variable(self%id_result, 'result', '', 'result', output=self%output)
   end subroutine product_initialize

   subroutine product_after_coupling(self)
      class (type_product), intent(inout) :: self

      if (associated(self%id_result%link%target, self%id_result%link%original)) &
         self%id_result%link%target%units = trim(self%id_term1%link%target%units) // ' ' // trim(self%id_term2%link%target%units)
   end subroutine

   subroutine product_do(self, _ARGUMENTS_DO_)
      class (type_product), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: term1, term2

      _LOOP_BEGIN_
         _GET_(self%id_term1, term1)
         _GET_(self%id_term2, term2)
         _SET_DIAGNOSTIC_(self%id_result, term1 * term2)
      _LOOP_END_
   end subroutine product_do

end module

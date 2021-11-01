#include "fabm_driver.h"

module fabm_builtin_reduction
   use fabm_types

   implicit none

   private

   public type_reduction_operator

   type, extends(type_base_model) :: type_reduction_operator
   contains
      procedure :: merge_components => reduction_operator_merge_components
   end type

contains

   subroutine reduction_operator_merge_components(self, log_unit)
      class (type_reduction_operator), intent(inout) :: self
      integer, optional,               intent(in)    :: log_unit
   end subroutine

end module

#include "fabm_driver.h"

module fabm_builtin_relaxation
   use fabm_types

   implicit none

   private

   public type_interior_relaxation

   type, extends(type_base_model) :: type_interior_relaxation
      type (type_state_variable_id) :: id_original
      type (type_dependency_id)     :: id_target
      type (type_dependency_id)     :: id_rate
      real(rk) :: rate
      logical :: rate_is_variable
   contains
      procedure :: initialize => interior_relaxation_initialize
      procedure :: do         => interior_relaxation_do
   end type

contains

   subroutine interior_relaxation_initialize(self,configunit)
      class (type_interior_relaxation),intent(inout),target :: self
      integer,                         intent(in)           :: configunit

      call self%register_state_dependency(self%id_original, 'original', '', 'variable that is to be relaxed')
      call self%register_dependency(self%id_target, 'target', '', 'target to relax towards')
      call self%get_parameter(self%rate_is_variable, 'rate_is_variable', '', 'use variable relaxation rate', default=.false.)
      if (self%rate_is_variable) then
         call self%register_dependency(self%id_rate, 'rate', 's-1', 'relaxation rate')
      else
         call self%get_parameter(self%rate, 'rate', 's-1', 'relaxation rate', minimum=0._rk)
      end if
   end subroutine interior_relaxation_initialize

   subroutine interior_relaxation_do(self, _ARGUMENTS_DO_)
      class (type_interior_relaxation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: original_value, target_value, relaxation_rate

      relaxation_rate = self%rate
      _LOOP_BEGIN_
         _GET_(self%id_original, original_value)
         _GET_(self%id_target, target_value)
         if (self%rate_is_variable) _GET_(self%id_rate, relaxation_rate)
         _ADD_SOURCE_(self%id_original, relaxation_rate * (target_value - original_value))
      _LOOP_END_
   end subroutine interior_relaxation_do

end module

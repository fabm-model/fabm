#include "fabm_driver.h"

module fabm_builtin_scale
   use fabm_types

   implicit none

   private

   public type_scaled_interior_variable, type_scaled_horizontal_variable

   type, extends(type_base_model) :: type_scaled_interior_variable
      type (type_dependency_id)          :: id_source
      type (type_diagnostic_variable_id) :: id_result
      real(rk)                           :: weight = 1.0_rk
      real(rk)                           :: offset = 0.0_rk
      logical                            :: include_background = .false.
   contains
      procedure :: initialize     => scaled_interior_variable_initialize
      procedure :: do             => scaled_interior_variable_do
      procedure :: after_coupling => scaled_interior_variable_after_coupling
   end type

   type, extends(type_base_model) :: type_scaled_horizontal_variable
      type (type_horizontal_dependency_id)          :: id_source
      type (type_horizontal_diagnostic_variable_id) :: id_result
      real(rk)                                      :: weight = 1.0_rk
      real(rk)                                      :: offset = 0.0_rk
      logical                                       :: include_background = .false.
   contains
      procedure :: initialize     => scaled_horizontal_variable_initialize
      procedure :: do_horizontal  => scaled_horizontal_variable_do_horizontal
      procedure :: after_coupling => scaled_horizontal_variable_after_coupling
   end type

contains

   subroutine scaled_interior_variable_initialize(self, configunit)
      class (type_scaled_interior_variable), intent(inout), target :: self
      integer,                               intent(in)            :: configunit
      call self%register_implemented_routines((/source_do/))
   end subroutine

   subroutine scaled_interior_variable_do(self, _ARGUMENTS_DO_)
      class (type_scaled_interior_variable), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value

      _CONCURRENT_LOOP_BEGIN_
         _GET_(self%id_source, value)
         _SET_DIAGNOSTIC_(self%id_result, self%offset + self%weight * value)
      _LOOP_END_
   end subroutine scaled_interior_variable_do

   subroutine scaled_interior_variable_after_coupling(self)
      class (type_scaled_interior_variable), intent(inout) :: self

      if (self%include_background) then
         self%offset = self%offset + self%weight * self%id_source%background
      else
         call self%id_result%link%target%background_values%set_value(self%weight * self%id_source%background)
      end if
   end subroutine scaled_interior_variable_after_coupling

   subroutine scaled_horizontal_variable_initialize(self, configunit)
      class (type_scaled_horizontal_variable), intent(inout), target :: self
      integer,                                 intent(in)            :: configunit
      call self%register_implemented_routines((/source_do_horizontal/))
   end subroutine

   subroutine scaled_horizontal_variable_do_horizontal(self, _ARGUMENTS_HORIZONTAL_)
      class (type_scaled_horizontal_variable), intent(in) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_

      real(rk) :: value

      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_source, value)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_result, self%offset + self%weight * value)
      _HORIZONTAL_LOOP_END_
   end subroutine scaled_horizontal_variable_do_horizontal

   subroutine scaled_horizontal_variable_after_coupling(self)
      class (type_scaled_horizontal_variable), intent(inout) :: self

      if (self%include_background) then
         self%offset = self%offset + self%weight * self%id_source%background
      else
         call self%id_result%link%target%background_values%set_value(self%weight * self%id_source%background)
      end if
   end subroutine scaled_horizontal_variable_after_coupling

end module

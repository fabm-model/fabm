#include "fabm_driver.h"

module fabm_builtin_scale
   use fabm_types
   use fabm_builtin_source, only: copy_fluxes, copy_horizontal_fluxes

   implicit none

   private

   public type_scaled_interior_variable, type_scaled_horizontal_variable

   type, extends(type_base_model) :: type_scaled_variable
      real(rk)                        :: weight = 1.0_rk
      real(rk)                        :: offset = 0.0_rk
      character(len=attribute_length) :: units         = ''
      real(rk)                        :: missing_value = -2.e20_rk
      logical                         :: act_as_state_variable = .false.
      logical                         :: include_background = .false.
      integer                         :: result_output = output_instantaneous
   end type

   type, extends(type_scaled_variable) :: type_scaled_interior_variable
      type (type_dependency_id)          :: id_source
      type (type_diagnostic_variable_id) :: id_result
   contains
      procedure :: initialize     => scaled_interior_variable_initialize
      procedure :: do             => scaled_interior_variable_do
      procedure :: after_coupling => scaled_interior_variable_after_coupling
   end type

   type, extends(type_scaled_variable) :: type_scaled_horizontal_variable
      type (type_horizontal_dependency_id)          :: id_source
      type (type_horizontal_diagnostic_variable_id) :: id_result
      integer                                       :: domain = domain_horizontal
   contains
      procedure :: initialize     => scaled_horizontal_variable_initialize
      procedure :: do_horizontal  => scaled_horizontal_variable_do_horizontal
      procedure :: after_coupling => scaled_horizontal_variable_after_coupling
   end type

contains

   subroutine scaled_interior_variable_initialize(self, configunit)
      class (type_scaled_interior_variable), intent(inout), target :: self
      integer,                               intent(in)            :: configunit

      call self%register_dependency(self%id_source, 'source', '', 'source variable')
      call self%register_diagnostic_variable(self%id_result, 'result', self%units, 'result', &
         missing_value=self%missing_value, output=self%result_output, act_as_state_variable=self%act_as_state_variable)
      if (self%act_as_state_variable) call copy_fluxes(self, self%id_result, self%id_source%link%target%name, scale_factor=1.0_rk / self%weight)
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

      call self%register_dependency(self%id_source, 'source', '', 'source variable')
      call self%register_diagnostic_variable(self%id_result, 'result', self%units, 'result', &
         missing_value=self%missing_value, output=self%result_output, act_as_state_variable= &
         self%act_as_state_variable, source=source_do_horizontal, domain=self%domain)
      if (self%act_as_state_variable) call copy_horizontal_fluxes(self, self%id_result, self%id_source%link%target%name, scale_factor=1.0_rk / self%weight)
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

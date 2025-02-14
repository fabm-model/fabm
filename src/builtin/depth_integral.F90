#include "fabm_driver.h"

module fabm_builtin_depth_integral
   use fabm_types

   implicit none

   private

   public type_depth_integral, type_bounded_depth_integral

   type, extends(type_base_model) :: type_depth_integral
      type (type_dependency_id)                     :: id_input
      type (type_dependency_id)                     :: id_thickness
      type (type_horizontal_diagnostic_variable_id) :: id_output
      logical                                       :: average = .false.
   contains
      procedure :: initialize     => depth_integral_initialize
      procedure :: do_column      => depth_integral_do_column
      procedure :: after_coupling => depth_integral_after_coupling
   end type

   type, extends(type_depth_integral) :: type_bounded_depth_integral
      real(rk) :: minimum_depth = 0.0_rk
      real(rk) :: maximum_depth = huge(1.0_rk)
   contains
      procedure :: initialize => bounded_depth_integral_initialize
      procedure :: do_column  => bounded_depth_integral_do_column
   end type

contains

   subroutine depth_integral_initialize(self, configunit)
      class (type_depth_integral), intent(inout), target :: self
      integer,                     intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_column/))
      call self%get_parameter(self%average, 'average', '', 'compute average by dividing integral by total height', default=self%average)
      call self%register_dependency(self%id_input, 'source', '', 'source')
      call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
      call self%register_diagnostic_variable(self%id_output, 'result', '', 'result', source=source_do_column)
   end subroutine depth_integral_initialize

   subroutine depth_integral_after_coupling(self)
      class (type_depth_integral), intent(inout) :: self

      if (associated(self%id_output%link%target, self%id_output%link%original)) then
         self%id_output%link%target%units = self%id_input%link%target%units
         if (.not. self%average) self%id_output%link%target%units = trim(self%id_output%link%target%units) // ' m'
      end if
   end subroutine depth_integral_after_coupling

   subroutine depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_depth_integral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: h, value, result, depth

      result = 0._rk
      depth = 0._rk
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_thickness, h)
         _GET_(self%id_input, value)
         depth = depth + h
         result = result + h * value
      _VERTICAL_LOOP_END_

      if (self%average .and. depth > 0._rk) result = result / depth
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output, result)
   end subroutine depth_integral_do_column

   subroutine bounded_depth_integral_initialize(self, configunit)
      class (type_bounded_depth_integral), intent(inout), target :: self
      integer,                             intent(in)            :: configunit

      call depth_integral_initialize(self, configunit)
      call self%get_parameter(self%minimum_depth, 'minimum_depth', 'm', 'minimum depth (distance from surface)', default=self%minimum_depth)
      call self%get_parameter(self%maximum_depth, 'maximum_depth', 'm', 'maximum depth (distance from surface)', default=self%maximum_depth)
   end subroutine bounded_depth_integral_initialize

   subroutine bounded_depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_bounded_depth_integral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: h, value, result, depth, layer_top, layer_bottom

      result = 0._rk
      depth = 0._rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_thickness, h)
         _GET_(self%id_input, value)
         layer_top = max(depth, self%minimum_depth)
         depth = depth + h
         layer_bottom = min(depth, self%maximum_depth)
         result = result + max(layer_bottom - layer_top, 0._rk) * value
      _DOWNWARD_LOOP_END_

      if (self%average .and. depth > self%minimum_depth) result = result / (min(self%maximum_depth, depth) - self%minimum_depth)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output, result)
   end subroutine bounded_depth_integral_do_column

end module

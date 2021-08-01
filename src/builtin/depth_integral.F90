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
      procedure :: do_column => bounded_depth_integral_do_column
   end type

contains

   subroutine depth_integral_initialize(self, configunit)
      class (type_depth_integral), intent(inout), target :: self
      integer,                     intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_column/))
      call self%register_dependency(self%id_input, 'source', '', 'source')
      call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
      call self%register_diagnostic_variable(self%id_output, 'result', '', 'result', source=source_do_column)
   end subroutine depth_integral_initialize

   subroutine depth_integral_after_coupling(self)
      class (type_depth_integral), intent(inout) :: self

      if (associated(self%id_output%link%target, self%id_output%link%original)) self%id_output%link%target%units = trim(self%id_input%link%target%units)//'*m'
   end subroutine depth_integral_after_coupling

   subroutine depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_depth_integral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: h, value, result, depth

      result = 0
      depth = 0
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_thickness, h)
         _GET_(self%id_input, value)
         depth = depth + h
         result = result + h * value
      _VERTICAL_LOOP_END_

      if (self%average) result = result / depth
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output, result)
   end subroutine depth_integral_do_column

   subroutine bounded_depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_bounded_depth_integral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: h, value, cum, depth
      logical :: started

      cum = 0
      depth = 0
      started = self%minimum_depth <= 0
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_thickness, h)
         _GET_(self%id_input, value)

         depth = depth + h
         if (.not.started) then
            ! Not yet at minimum depth before
            if (depth >= self%minimum_depth) then
               ! Now crossing minimum depth interface
               started = .true.
               h = depth-self%minimum_depth
            end if
         elseif (depth > self%maximum_depth) then
            ! Now crossing maximum depth interface; subtract part of layer height that is not included
            h = h - (depth - self%maximum_depth)
            _VERTICAL_LOOP_EXIT_
         end if
         cum = cum + h * value
      _VERTICAL_LOOP_END_

      if (.not. self%average) then
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output, cum)
      elseif (depth > self%minimum_depth) then
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output, cum / (min(self%maximum_depth, depth) - self%minimum_depth))
      endif
   end subroutine bounded_depth_integral_do_column

end module

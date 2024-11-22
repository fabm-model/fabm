#include "fabm_driver.h"

module fabm_builtin_constant
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_interior_constant
      type (type_diagnostic_variable_id) :: id_constant
   contains
      procedure :: initialize => interior_constant_initialize
   end type

   type, extends(type_base_model), public :: type_horizontal_constant
      type (type_horizontal_diagnostic_variable_id) :: id_constant
   contains
      procedure :: initialize => horizontal_constant_initialize
   end type

   type, extends(type_base_model), public :: type_surface_constant
      type (type_surface_diagnostic_variable_id) :: id_constant
   contains
      procedure :: initialize => surface_constant_initialize
   end type

   type, extends(type_base_model), public :: type_bottom_constant
      type (type_bottom_diagnostic_variable_id) :: id_constant
   contains
      procedure :: initialize => bottom_constant_initialize
   end type

   type, extends(type_base_model), public :: type_global_constant
   contains
      procedure :: initialize => global_constant_initialize
   end type

contains

   subroutine interior_constant_initialize(self, configunit)
      class (type_interior_constant), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk)                               :: value
      type (type_interior_standard_variable) :: standard_variable

      call self%register_implemented_routines()
      call self%get_parameter(standard_variable%name, 'standard_name', '', 'standard name', default='')
      call self%get_parameter(value, 'value', '', 'value')
      if (standard_variable%name /= '') then
         ! Note: for Cray 8.3.7, standard_variable needs to be declared separately.
         ! It cannot be constructed on the fly within the function call
         call self%register_diagnostic_variable(self%id_constant, 'data', '', 'data', missing_value=value, &
            output=output_none, standard_variable=standard_variable, source=source_constant)
      else
         call self%register_diagnostic_variable(self%id_constant, 'data', '', 'data', missing_value=value, &
            output=output_none, source=source_constant)
      end if
   end subroutine interior_constant_initialize

   subroutine horizontal_constant_initialize(self, configunit)
      class (type_horizontal_constant), intent(inout), target :: self
      integer,                          intent(in)            :: configunit

      real(rk)                                 :: value
      type (type_horizontal_standard_variable) :: standard_variable

      call self%register_implemented_routines()
      call self%get_parameter(standard_variable%name, 'standard_name', '', 'standard name', default='')
      call self%get_parameter(value, 'value', '', 'value')
      if (standard_variable%name /= '') then
         ! Note: for Cray 8.3.7, standard_variable needs to be declared separately.
         ! It cannot be constructed on the fly within the function call
         call self%register_diagnostic_variable(self%id_constant, 'data', '', 'data', missing_value=value, &
            output=output_none, standard_variable=standard_variable, source=source_constant)
      else
         call self%register_diagnostic_variable(self%id_constant, 'data', '', 'data', missing_value=value, &
            output=output_none, source=source_constant)
      end if
   end subroutine horizontal_constant_initialize

   subroutine surface_constant_initialize(self, configunit)
      class (type_surface_constant), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      real(rk)                              :: value
      type (type_surface_standard_variable) :: standard_variable

      call self%register_implemented_routines()
      call self%get_parameter(standard_variable%name, 'standard_name', '', 'standard name', default='')
      call self%get_parameter(value, 'value', '', 'value')
      if (standard_variable%name /= '') then
         ! Note: for Cray 8.3.7, standard_variable needs to be declared separately.
         ! It cannot be constructed on the fly within the function call
         call self%register_diagnostic_variable(self%id_constant, 'data', '', 'data', missing_value=value, &
            output=output_none, standard_variable=standard_variable, source=source_constant)
      else
         call self%register_diagnostic_variable(self%id_constant, 'data', '', 'data', missing_value=value, &
            output=output_none, source=source_constant)
      end if
   end subroutine surface_constant_initialize

   subroutine bottom_constant_initialize(self, configunit)
      class (type_bottom_constant), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      real(rk)                             :: value
      type (type_bottom_standard_variable) :: standard_variable

      call self%register_implemented_routines()
      call self%get_parameter(standard_variable%name, 'standard_name', '', 'standard name', default='')
      call self%get_parameter(value, 'value', '', 'value')
      if (standard_variable%name /= '') then
         ! Note: for Cray 8.3.7, standard_variable needs to be declared separately.
         ! It cannot be constructed on the fly within the function call
         call self%register_diagnostic_variable(self%id_constant, 'data', '', 'data', missing_value=value, &
            output=output_none, standard_variable=standard_variable, source=source_constant)
      else
         call self%register_diagnostic_variable(self%id_constant, 'data', '', 'data', missing_value=value, &
            output=output_none, source=source_constant)
      end if
   end subroutine bottom_constant_initialize

   subroutine global_constant_initialize(self, configunit)
      class (type_global_constant), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      real(rk)                             :: value
      type (type_global_standard_variable) :: standard_variable

      call self%register_implemented_routines()
      call self%get_parameter(standard_variable%name, 'standard_name', '', 'standard name', default='')
      call self%get_parameter(value, 'value', '', 'value')
      if (standard_variable%name /= '') then
         ! Note: for Cray 8.3.7, standard_variable needs to be declared separately.
         ! It cannot be constructed on the fly within the function call
         call self%add_scalar_variable('data', '', 'data', fill_value=value, &
            output=output_none, source=source_constant, standard_variable=standard_variable)
      else
         call self%add_scalar_variable('data', '', 'data', fill_value=value, &
            output=output_none, source=source_constant)
      end if
   end subroutine global_constant_initialize

end module

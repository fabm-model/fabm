#include "fabm_driver.h"

module fabm_builtin_time_filter
   use fabm_types
   use fabm_expressions

   implicit none

   private

   public type_interior_running_mean, type_surface_running_mean

   type, extends(type_base_model) :: type_interior_running_mean
      type (type_diagnostic_variable_id) :: id_meanout
      type (type_dependency_id)          :: id_source, id_meanin
   contains
      procedure :: initialize => interior_running_mean_initialize
      procedure :: do         => interior_running_mean_do
   end type

   type, extends(type_base_model) :: type_surface_running_mean
      type (type_surface_diagnostic_variable_id) :: id_meanout
      type (type_horizontal_dependency_id)       :: id_source, id_meanin
   contains
      procedure :: initialize => surface_running_mean_initialize
      procedure :: do_surface => surface_running_mean_do_surface
   end type

contains

   subroutine interior_running_mean_initialize(self, configunit)
      class (type_interior_running_mean), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      real(rk) :: window, missing_value
      integer :: n

      call self%get_parameter(window, 'window', 's', 'window size')
      call self%get_parameter(n, 'n', '', 'number of bins')
      call self%get_parameter(missing_value, 'missing_value', '', 'missing value to until the full window size has been covered', default=-2e20_rk)

      call self%register_dependency(self%id_source, 'source', '', 'variable for which to compute running mean')
      call self%register_dependency(self%id_meanin, temporal_mean(self%id_source, period=window, resolution=window / n, missing_value=missing_value))
      call self%register_diagnostic_variable(self%id_meanout, 'mean',  '', 'running mean')
   end subroutine

   subroutine interior_running_mean_do(self,_ARGUMENTS_DO_)
      class (type_interior_running_mean), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: mean

      _LOOP_BEGIN_
         _GET_(self%id_meanin, mean)
         _SET_DIAGNOSTIC_(self%id_meanout ,mean)
      _LOOP_END_
   end subroutine

   subroutine surface_running_mean_initialize(self, configunit)
      class (type_surface_running_mean), intent(inout), target :: self
      integer,                           intent(in)            :: configunit

      real(rk) :: window, missing_value
      integer :: n

      call self%get_parameter(window, 'window', 's', 'window size')
      call self%get_parameter(n, 'n', '', 'number of bins')
      call self%get_parameter(missing_value, 'missing_value', '', 'missing value to until the full window size has been covered', default=-2e20_rk)

      call self%register_dependency(self%id_source, 'source', '', 'variable for which to compute running mean')
      call self%register_dependency(self%id_meanin, temporal_mean(self%id_source, period=window, resolution=window / n))
      call self%register_diagnostic_variable(self%id_meanout, 'mean',  '', 'running mean')
   end subroutine

   subroutine surface_running_mean_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_surface_running_mean), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: mean

      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_meanin, mean)
         _SET_SURFACE_DIAGNOSTIC_(self%id_meanout ,mean)
      _SURFACE_LOOP_END_
   end subroutine

end module

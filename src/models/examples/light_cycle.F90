#include "fabm_driver.h"

module examples_light_cycle

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_examples_light_cycle
      type (type_global_dependency_id)           :: id_yearday
      type (type_surface_diagnostic_variable_id) :: id_swr_sf

      real(rk) :: light_period
      real(rk) :: offset
      real(rk) :: swr
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_examples_light_cycle), intent(inout), target :: self
      integer,                           intent(in)            :: configunit

      call self%register_dependency(self%id_yearday, standard_variables%number_of_days_since_start_of_the_year)
      call self%register_diagnostic_variable(self%id_swr_sf, 'swr_sf', 'W m-2', 'surface downward shortwave radiation', standard_variable=standard_variables%surface_downwelling_shortwave_flux)
      call self%get_parameter(self%swr,          'swr',          'W m-2', 'surface downward shortwave radiation during light phase', minimum=0.0_rk)
      call self%get_parameter(self%light_period, 'light_period', 'h',     'duration of light phase',                                 minimum=0.0_rk,   maximum=24._rk, scale_factor=1._rk/24._rk)
      call self%get_parameter(self%offset,       'offset',       'h',     'start of light period (relative to midnight)',            minimum=-24.0_rk, maximum=24._rk, scale_factor=1._rk/24._rk, default=self%light_period*12)
   end subroutine initialize

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_examples_light_cycle), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: yearday

      _SURFACE_LOOP_BEGIN_
         _GET_GLOBAL_(self%id_yearday,yearday)
         if (modulo(yearday-self%offset, 1.0_rk) < self%light_period) then
            _SET_SURFACE_DIAGNOSTIC_(self%id_swr_sf, self%swr)
         else
            _SET_SURFACE_DIAGNOSTIC_(self%id_swr_sf, 0.0_rk)
         end if
      _SURFACE_LOOP_END_
   end subroutine do_surface

end module

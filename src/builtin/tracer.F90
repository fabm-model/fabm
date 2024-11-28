#include "fabm_driver.h"

module fabm_builtin_tracer

   ! This model describe a single passive tracer. Optionally, a vertical velocity
   ! (sinking/floating) and light attenuation coefficient can be specified.
   ! The unit is mol/m^3 by default, but may be explicitly configured as well.

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_tracer
      ! Variables
      type (type_state_variable_id) :: id_c
   contains
      procedure initialize
   end type

contains

   subroutine initialize(self, configunit)
      class (type_tracer), intent(inout), target :: self
      integer,             intent(in)            :: configunit

      character(len=64)               :: units
      real(rk)                        :: vertical_velocity
      real(rk)                        :: specific_light_attenuation
      real(rk), parameter             :: days_per_second = 1.0_rk/86400.0_rk
      logical                         :: conserved
      character(len=attribute_length) :: standard_name

      call self%register_implemented_routines()

      ! Retrieve parameter values
      call self%get_parameter(units, 'units', default='mol m-3')
      call self%get_parameter(vertical_velocity, 'vertical_velocity', 'm d-1', 'vertical velocity (negative for settling, positive for rising)', default=0.0_rk, scale_factor=days_per_second)
      call self%get_parameter(specific_light_attenuation, 'specific_light_attenuation', 'm-1 ('//trim(units)//')-1', 'specific light attenuation', default=0.0_rk)
      call self%get_parameter(conserved, 'conserved', '', 'treat tracer as conserved quantity (activates budget tracking in host)', default=.false.)

      ! Register state variables
      call self%register_state_variable(self%id_c, 'c', units, 'concentration', 1.0_rk, &
         minimum=0.0_rk, vertical_movement=vertical_velocity, specific_light_extinction=specific_light_attenuation)

      if (conserved) then
         ! User has activated conservation checking.
         ! We create a "conserved quantity" for this tracer, which should induce the host model to compute and save the integral across the domain.
         standard_name = get_safe_name(trim(self%get_path()) // '_total')
         call self%add_to_aggregate_variable(type_universal_standard_variable(name=standard_name(2:), units=units, aggregate_variable=.true., conserved=.true.), self%id_c)
      end if
   end subroutine initialize

end module

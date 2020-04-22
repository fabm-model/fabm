#include "fabm_driver.h"

module bb_passive

   ! This model describe a single passive tracer. Optionally, a vertical velocity
   ! (sinking/floating), light attenuation coefficient and surface flux can be specified.
   ! The unit is mol/m\^3 by default, but may be explicitly configured as well.

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_bb_passive
      ! Variables
      type (type_state_variable_id) :: id_tracer

      ! Parameters
      real(rk) :: surface_flux
   contains
      procedure initialize
      procedure do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_bb_passive), intent(inout), target :: self
      integer,                 intent(in)            :: configunit

      real(rk)                  :: vertical_velocity
      real(rk)                  :: specific_light_attenuation
      character(len=64)         :: units
      real(rk), parameter       :: days_per_second = 1.0_rk/86400.0_rk
      logical                   :: conserved
      character(len=attribute_length) :: standard_name

      ! Retrieve parameter values
      call self%get_parameter(units, 'units', default='mol m-3')
      call self%get_parameter(vertical_velocity, 'vertical_velocity', 'm d-1', 'vertical velocity (negative for settling, positive for rising)', default=0.0_rk, scale_factor=days_per_second)
      call self%get_parameter(specific_light_attenuation, 'specific_light_attenuation', 'm-1 ('//trim(units)//')-1', 'specific light attenuation', default=0.0_rk)
      call self%get_parameter(self%surface_flux, 'surface_flux', trim(units)//' m d-1', 'surface flux (positive for into the water)', default=0.0_rk, scale_factor=days_per_second)
      call self%get_parameter(conserved, 'conserved', '', 'treat tracer as conserved quantity (activates budget tracking in host)', default=.false.)

      ! Register state variables
      call self%register_state_variable(self%id_tracer, 'c', units, 'concentration', 1.0_rk, &
         minimum=0.0_rk, vertical_movement=vertical_velocity, specific_light_extinction=specific_light_attenuation)

      if (conserved) then
         ! User has activated conservation checking.
         ! We create a "conserved quantity" for this tracer, which should induce the host model to compute and save the integral across the domain.
         standard_name = get_safe_name(trim(self%get_path()) // '_total')
         call self%add_to_aggregate_variable(type_universal_standard_variable(name=standard_name(2:), units=units, aggregate_variable=.true., conserved=.true.), self%id_tracer)
      end if
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_bb_passive), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      ! Enter spatial loops (if any)
      _SURFACE_LOOP_BEGIN_

         ! Transfer surface exchange value to FABM.
         _ADD_SURFACE_FLUX_(self%id_tracer, self%surface_flux)

      ! Leave spatial loops (if any)
      _SURFACE_LOOP_END_
   end subroutine do_surface

end module bb_passive

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------


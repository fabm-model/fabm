#include "fabm_driver.h"

module akvaplan_plume_injection

   ! A model for injection of a tracer in the form of a plume with specific density, originating from the bottom.
   !
   ! It assumes no ambient water is entrained in the plume, allowing it to rise up to the height where the ambient density drops below the [original] density of the plume.
   ! It there is no such height (the ambient density is everywhere greater than that of the plume), the plume will rise all the way to the water surface.
   ! Note that the assumption of zero entrainment is highly simplistic; more detailed models exist, such as http://dx.doi.org/10.1098/rspa.1956.0011
   !
   ! The model takes a depth-integrated (2D) tracer flux (quantity m-2 time-1), identifies the layer to which the plume will rise, and then computes the depth-explicit (3D)
   ! tracer source (quantity m-3 time-1), which is non-zero only in the layer of injection, where it equals the depth-integrated flux divided by layer height.
   ! To actually associate that source field with a tracer, create an "interior_source" model instance coupled to the target tracer and the source field.
   !
   ! Copyright (C) 2016 - Akvaplan-niva

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_plume_injection
      type (type_horizontal_dependency_id) :: id_flux_int
      type (type_dependency_id)            :: id_rho
      type (type_dependency_id)            :: id_h
      type (type_diagnostic_variable_id)   :: id_flux

      ! Parameters
      real(rk) :: rho
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_column
   end type

contains

   subroutine initialize(self,configunit)
      class (type_plume_injection),intent(inout),target :: self
      integer,                     intent(in)           :: configunit

      ! Obtain parameter values
      call self%get_parameter(self%rho, 'rho', 'kg m-3', 'plume density')

      ! Register dependencies
      call self%register_dependency(self%id_flux_int,'flux_int','quantity m-2 s-1','depth-integrated tracer flux')
      call self%register_dependency(self%id_rho,standard_variables%density)
      call self%register_dependency(self%id_h,standard_variables%cell_thickness)

      ! Register output (depth-explicit tracer source)
      call self%register_diagnostic_variable(self%id_flux,'flux','quantity m-3 s-1','depth-explicit tracer flux',source=source_do_column,prefill_value=0.0_rk)
   end subroutine initialize

   subroutine do_column(self,_ARGUMENTS_DO_COLUMN_)
      class (type_plume_injection),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: flux_int,rho,h

      ! Retrieve tracer flux at bottom (quantity m-2 s-1)
      _GET_HORIZONTAL_(self%id_flux_int,flux_int)

      ! If the tracer flux is zero, there is no need to search for the layer to which the plume will rise - just return.
      if (flux_int==0.0_rk) return

      ! Loop from bottom to surface to find the first layer with density < plume density.
      ! That's the layer we want to inject the plume contents into.
      _UPWARD_LOOP_BEGIN_
         _GET_(self%id_rho,rho)
         if (rho<self%rho) then
            ! Target layer found - convert the depth-integrated tracer flux into a local source term by dividing by layer height, then return.
            _GET_(self%id_h,h)
            _SET_DIAGNOSTIC_(self%id_flux,flux_int/h)
            return
         end if
      _UPWARD_LOOP_END_

      ! We have reached the water surface without finding a suitable layer [with density < plume density] to inject into.
      ! Inject the plume contents into the surface layer as last resort.
      _MOVE_TO_SURFACE_
      _GET_(self%id_h,h)
      _SET_DIAGNOSTIC_(self%id_flux,flux_int/h)
   end subroutine do_column

end module

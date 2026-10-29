#include "fabm_driver.h"

module bb_mixed_layer_depth

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_bb_mixed_layer_depth
      type (type_dependency_id)                  :: id_rho
      type (type_dependency_id)                  :: id_dz
      type (type_surface_diagnostic_variable_id) :: id_mld

      ! Parameters
      real(rk) :: drho
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_column
   end type type_bb_mixed_layer_depth

contains

   subroutine initialize(self, configunit)
      class (type_bb_mixed_layer_depth), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      call self%get_parameter(self%drho, 'drho', 'kg m-3', 'density difference between surface and base of the mixed layer', &
        default=0.125_rk, minimum=0.0_rk) 

      call self%register_diagnostic_variable(self%id_mld, 'mld', 'm', 'mixed layer depth', source=source_do_column)

      call self%register_dependency(self%id_rho, standard_variables%density)
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
   end subroutine

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_bb_mixed_layer_depth), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: rho, rho_sf, rho_prev, z, dz, dz_prev

      z = 0.0_rk
      dz_prev = 0.0_rk
      rho_sf = -1.0_rk
      rho_prev = -1.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer thickness (m)
         _GET_(self%id_rho,rho)   ! Density (kg m-3)

         if (rho_sf == -1.0_rk) rho_sf = rho
         if (rho >= rho_sf + self%drho) then
            z = z + (rho_sf + self%drho - rho_prev) / (rho - rho_prev) * 0.5_rk * (dz_prev + dz)
            _SET_SURFACE_DIAGNOSTIC_(self%id_mld,z)
            return
         end if
         z = z + 0.5_rk * (dz_prev + dz)
         rho_prev = rho
         dz_prev = dz
      _DOWNWARD_LOOP_END_
      z = z + 0.5_rk * dz_prev
      _SET_SURFACE_DIAGNOSTIC_(self%id_mld,z)
   end subroutine do_column

end module bb_mixed_layer_depth

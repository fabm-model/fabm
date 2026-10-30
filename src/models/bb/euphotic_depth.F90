#include "fabm_driver.h"

module bb_euphotic_depth

   ! Euphotic depth model based on PAR attenuation
   ! E.g. Kirk (1994; https://doi.org/10.1017/CBO9781139168212)

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_bb_euphotic_depth
      type (type_dependency_id)                  :: id_dz
      type (type_dependency_id)                  :: id_K_d
      type (type_surface_diagnostic_variable_id) :: id_z_eu

      ! Parameters
      real(rk) :: logf
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_column
   end type type_bb_euphotic_depth

contains

   subroutine initialize(self, configunit)
      class (type_bb_euphotic_depth), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk) :: f

      call self%get_parameter(f, 'f', '1', 'fraction of photosynthetically active radiation at bottom of euphotic layer', &
        default=0.01_rk, minimum=0.0_rk, maximum=1.0_rk) 
      self%logf = -log(f)

      call self%register_diagnostic_variable(self%id_z_eu, 'z_eu', 'm', 'euphotic depth', source=source_do_column)

      call self%register_dependency(self%id_K_d, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_dz,  standard_variables%cell_thickness)
   end subroutine

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_bb_euphotic_depth), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: K_d_int, z, dz, K_d

      K_d_int = 0.0_rk
      z = 0.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer thickness (m)
         _GET_(self%id_K_d,K_d)   ! PAR attenuation (m-1)

         K_d_int = K_d_int + K_d * dz
         z = z + dz
         if (K_d_int >= self%logf) then
            ! Target value for depth-integrated PAR attenuation is found
            ! in the current layer. Linearly interpolate attentuation to
            ! find the exact depth where the target is reached.
            z = z - (K_d_int - self%logf) / K_d
            _VERTICAL_LOOP_EXIT_
         end if
      _DOWNWARD_LOOP_END_
      _SET_SURFACE_DIAGNOSTIC_(self%id_z_eu,z)
   end subroutine do_column

end module bb_euphotic_depth

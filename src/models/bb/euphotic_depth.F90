#include "fabm_driver.h"

module bb_euphotic_depth

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_bb_euphotic_depth
      type (type_dependency_id)                  :: id_dz
      type (type_dependency_id)                  :: id_ext
      type (type_surface_diagnostic_variable_id) :: id_zeu

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

      call self%register_diagnostic_variable(self%id_zeu, 'zeu', 'm', 'euphotic depth', source=source_do_column)

      call self%register_dependency(self%id_ext,  standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_dz,   standard_variables%cell_thickness)
   end subroutine

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_bb_euphotic_depth), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: ext_int, z, dz, ext

      ext_int = 0.0_rk
      z = 0.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer thickness (m)
         _GET_(self%id_ext,ext)   ! PAR attenuation (m-1)

         ext_int = ext_int + ext * dz
         z = z + dz
         if (ext_int >= self%logf) then
            z = z - (ext_int - self%logf) / ext
            _VERTICAL_LOOP_EXIT_
         end if
      _DOWNWARD_LOOP_END_
      _SET_SURFACE_DIAGNOSTIC_(self%id_zeu,z)
   end subroutine do_column

end module bb_euphotic_depth

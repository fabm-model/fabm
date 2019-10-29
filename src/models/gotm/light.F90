#include "fabm_driver.h"

module gotm_light

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_gotm_light
      ! Identifiers for dependencies [model inputs]
      type (type_horizontal_dependency_id)          :: id_swr0 ! Surface shortwave radiation
      type (type_dependency_id)                     :: id_dz   ! Cell thickness
      type (type_dependency_id)                     :: id_ext  ! Attentuation coefficient for PAR

      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)            :: id_par  ! Photosynthetically active radiation
      type (type_diagnostic_variable_id)            :: id_swr  ! Shortwave radiation
      type (type_horizontal_diagnostic_variable_id) :: id_par0 ! Surface photosynthetically active radiation

      ! Parameters
      real(rk) :: a,g1,g2
   contains
!     Model procedures
      procedure :: initialize
      procedure :: get_light
   end type type_gotm_light

contains

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_gotm_light),intent(inout),target :: self
      integer,                intent(in)           :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call self%get_parameter(self%a, 'a', '-','non-visible fraction of shortwave radiation',default=0.58_rk) 
      call self%get_parameter(self%g1,'g1','m','e-folding depth of non-visible fraction',    default=0.35_rk)
      call self%get_parameter(self%g2,'g2','m','e-folding depth of visible fraction',        default=23.0_rk) 

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_swr,'swr','W/m^2','shortwave radiation', &
              standard_variable=standard_variables%downwelling_shortwave_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_par,'par','W/m^2','photosynthetically active radiation', &
              standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux,source=source_do_column)
      call self%register_diagnostic_variable(self%id_par0,'par0','W/m^2','surface photosynthetically active radiation', &
              standard_variable=standard_variables%surface_downwelling_photosynthetic_radiative_flux,source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_swr0,standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_ext, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_dz,  standard_variables%cell_thickness)
   end subroutine
   
   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
      class (type_gotm_light),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_

      real(rk) :: swr0,dz,swr,par,z,ext,bioext

      _GET_HORIZONTAL_(self%id_swr0,swr0)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par0,swr0*(1-self%a))
      z = 0
      bioext = 0
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_ext,ext)   ! PAR attenuation (m-1)

         ! Set depth to centre of layer
         z = z + dz/2
         bioext = bioext + ext*dz/2

         ! Calculate photosynthetically active radiation (PAR), shortwave radiation, and PAR attenuation.
         par = swr0*(1-self%A)*exp(-z/self%g2-bioext)
         swr = par+swr0*self%A*exp(-z/self%g1)

         ! Move to bottom of layer
         z = z + dz/2
         bioext = bioext + ext*dz/2

         _SET_DIAGNOSTIC_(self%id_swr,swr) ! Shortwave radiation at layer centre
         _SET_DIAGNOSTIC_(self%id_par,par) ! Photosynthetically active radiation at layer centre
      _VERTICAL_LOOP_END_

   end subroutine get_light

end module gotm_light

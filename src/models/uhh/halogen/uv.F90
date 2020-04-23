#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_uv - UV irradiance model
!
! !INTERFACE:
   module fabm_uhh_uv
!
! !DESCRIPTION:
!
! 
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_uv
!     Variable identifiers
      type (type_diagnostic_variable_id)    :: id_uv
      type (type_dependency_id)             :: id_pres
      type (type_horizontal_dependency_id)  :: id_I0
!     Model parameters
      real(rk) :: aw
      real(rk) :: I0_const
      real(rk) :: uv_fraction

     contains

      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   contains


   !> Initialise the UV model
   subroutine initialize(self,configunit)
   class (type_uhh_uv), intent(inout), target :: self
   integer,                     intent(in)            :: configunit

   real(rk)          :: I0_const=0.0_rk
   real(rk)          :: UV_fraction=0.0_rk
   real(rk)          :: attenuation_coefficient=0.0_rk
   character(len=64) :: I0_variable

   namelist /uhh_uv/ I0_const, UV_fraction, I0_variable, attenuation_coefficient

   I0_variable = ''
   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_uv,err=99)

   call self%get_parameter(self%I0_const, 'I0_const', default=I0_const) 

   call self%get_parameter(self%uv_fraction, 'UV_fraction', default=UV_fraction)
 
   call self%get_parameter(self%aw, 'attenuation_coefficient', default=attenuation_coefficient)
   
   ! Register dependencies
   call self%register_dependency(self%id_I0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_pres, standard_variables%pressure)

   ! Register UV variable
   call self%register_diagnostic_variable(self%id_uv,'rad','W/m2', &
         'ultra-violet irradiance', output=output_instantaneous)

   return
99 call self%fatal_error('fabm_uhh_uv','Error reading namelist uhh_uv')
   end subroutine initialize

   !> Right hand sides of halogen model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_uv), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
   real(rk) :: pres, I0, uv

   _LOOP_BEGIN_

   _GET_(self%id_pres, pres)
   _GET_HORIZONTAL_(self%id_I0, I0)

   ! calculate UV irradiance based on pressure (dbar)
   ! here assumed [dbar] appr. [m]
   uv = self%uv_fraction*I0 * exp(-self%aw*pres)

   ! Set diagnostic variables
   _SET_DIAGNOSTIC_(self%id_uv,uv)

   _LOOP_END_
   end subroutine do

end module fabm_uhh_uv

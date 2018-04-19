#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bb_passive --- passive tracer model
!
! !INTERFACE:
   module bb_passive
!
! !DESCRIPTION:
! This model describe a single passive tracer. Optionally, a vertical velocity
! (sinking/floating), light attenuation coefficient and surface flux can be specified.
! The unit is mol/m\^3 by default, but may be explicitly configured as well.
!
! !USES:
   use fabm_types

   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model), public :: type_bb_passive
!     Variable identifiers
      type (type_state_variable_id) :: id_tracer

!     Model parameters
      real(rk)                      :: surface_flux

      contains

      procedure initialize
      procedure get_surface_exchange
   end type
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the passive tracer model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the bb\_passive namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_bb_passive), intent(inout), target :: self
   integer,                 intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk)                  :: vertical_velocity
   real(rk)                  :: specific_light_attenuation
   character(len=64)         :: units
   real(rk), parameter       :: days_per_second = 1.0_rk/86400.0_rk
   logical                   :: conserved
   character(len=attribute_length) :: standard_name
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%get_parameter(units,'units',default='mol m-3')
   call self%get_parameter(vertical_velocity,'vertical_velocity','m d-1','vertical velocity (negative for settling, positive for rising)',default=0.0_rk,scale_factor=days_per_second)
   call self%get_parameter(specific_light_attenuation,'specific_light_attenuation','m-1 ('//trim(units)//')-1','specific light attenuation',default=0.0_rk)
   call self%get_parameter(self%surface_flux,'surface_flux',trim(units)//' m d-1','surface flux (positive for into the water)',default=0.0_rk,scale_factor=days_per_second)
   call self%get_parameter(conserved,'conserved','','whether this is a conserved quantity (activates budget tracking in host)',default=.false.)

   ! Register state variables
   call self%register_state_variable(self%id_tracer, &
                    'c',units,'concentration', &
                    1.0_rk,minimum=0.0_rk, &
                    vertical_movement=vertical_velocity, &
                    specific_light_extinction=specific_light_attenuation)

   if (conserved) then
      standard_name = get_safe_name(trim(self%get_path())//'_total')
      call self%add_to_aggregate_variable(type_bulk_standard_variable(name=standard_name(2:),units=units,aggregate_variable=.true.,conserved=.true.),self%id_tracer)
   end if

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Air-sea exchange for the passive tracer model
!
! !INTERFACE:
   subroutine get_surface_exchange(self,_ARGUMENTS_DO_SURFACE_)
!
! !INPUT PARAMETERS:
   class (type_bb_passive), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Transfer surface exchange value to FABM.
   _SET_SURFACE_EXCHANGE_(self%id_tracer,self%surface_flux)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine get_surface_exchange
!EOC

!-----------------------------------------------------------------------

   end module bb_passive

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------


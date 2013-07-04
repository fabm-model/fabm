#ifdef _FABM_F2003_

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_bb_passive --- passive tracer model
!
! !INTERFACE:
   module fabm_bb_passive
!
! !DESCRIPTION:
! This model describe a single passive tracer. Optionally, a vertical velocity
! (sinking/floating) and a light absorption coefficient can be specified. The
! unit is mol/m\^3 by default, but may be explicitly set in the namelist
! instead.
!
! !USES:
   use fabm_types
   use fabm_driver

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
   real(rk)                  :: initial_concentration
   real(rk)                  :: vertical_velocity
   real(rk)                  :: specific_light_absorption
   real(rk)                  :: surface_flux
   character(len=64)         :: unit
   real(rk), parameter       :: secs_pr_day = 86400.0_rk

   namelist /bb_passive/     initial_concentration,vertical_velocity, &
                             specific_light_absorption,surface_flux
!EOP
!-----------------------------------------------------------------------
!BOC
   initial_concentration     = 1.0_rk
   vertical_velocity         = 0.0_rk
   specific_light_absorption = 0.0_rk
   surface_flux              = 0.0_rk
   unit                      = 'mol/m**3'

   ! Read the namelist
   read(configunit,nml=bb_passive,err=99,end=100)

   self%surface_flux = surface_flux/secs_pr_day

   ! Register state variables
   call self%register_state_variable(self%id_tracer, &
                    'tracer',unit,'tracer', &
                    initial_concentration,minimum=0.0_rk, &
                    vertical_movement=vertical_velocity/secs_pr_day, &
                    specific_light_extinction=specific_light_absorption)

   return

99 call fatal_error('bb_passive_create','Error reading namelist bb_passive')

100 call fatal_error('bb_passive_create','Namelist bb_passive was not found.')

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

   end module fabm_bb_passive

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#endif

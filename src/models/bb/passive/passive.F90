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

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_bb_passive, bb_passive_init, bb_passive_get_surface_exchange
!
! !PUBLIC DERIVED TYPES:
   type type_bb_passive
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_ :: id_tracer
      REALTYPE                 :: surface_flux
   end type
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the passive tracer model
!
! !INTERFACE:
   subroutine bb_passive_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the bb\_passive namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   type (type_bb_passive),       intent(out)   :: self
   _CLASS_ (type_model_info),    intent(inout) :: modelinfo
   integer,                      intent(in)    :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                  :: initial_concentration     = _ONE_
   REALTYPE                  :: vertical_velocity         = _ZERO_
   REALTYPE                  :: specific_light_absorption = _ZERO_
   REALTYPE                  :: surface_flux              = _ZERO_
   character(len=64)         :: unit                      = 'mol/m**3'
   REALTYPE, parameter       :: secs_pr_day=86400.
   namelist /bb_passive/ initial_concentration,vertical_velocity,specific_light_absorption,surface_flux
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=bb_passive,err=99)

   ! Register state variables
   self%id_tracer = register_state_variable(modelinfo,'tracer',unit,'tracer',     &
                                            initial_concentration,minimum=_ZERO_, &
                                            vertical_movement=vertical_velocity/secs_pr_day, &
                                            specific_light_extinction=specific_light_absorption)

   self%surface_flux = surface_flux/secs_pr_day

   return

99 call fatal_error('bb_passive_init','Error reading namelist bb_passive')

   end subroutine bb_passive_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Air-sea exchange for the passive tracer model
!
! !INTERFACE:
   subroutine bb_passive_get_surface_exchange(self,_FABM_ARGS_GET_SURFACE_EXCHANGE_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_bb_passive), intent(in)    :: self
   _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_HZ_LOOP_BEGIN_

   ! Transfer surface exchange value to FABM.
   _SET_SURFACE_EXCHANGE_(self%id_tracer,self%surface_flux)

   ! Leave spatial loops (if any)
   _FABM_HZ_LOOP_END_

   end subroutine bb_passive_get_surface_exchange
!EOC

!-----------------------------------------------------------------------

   end module fabm_bb_passive

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

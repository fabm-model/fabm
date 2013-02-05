#include "fabm_driver.h"

#ifdef _FABM_F2003_

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
   public bb_passive_create, type_bb_passive
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model) :: type_bb_passive
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_ :: id_tracer
      REALTYPE                 :: surface_flux
      
      contains
      
      procedure get_surface_exchange
   end type type_bb_passive
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
   function bb_passive_create(configunit,name,parent) result(self)
!
! !DESCRIPTION:
!  Here, the bb\_passive namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   integer,                        intent(in)    :: configunit
   character(len=*),               intent(in)    :: name
   class (type_base_model),target, intent(inout) :: parent
   type (type_bb_passive),         pointer       :: self
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

   namelist /bb_passive/     initial_concentration,vertical_velocity, &
                             specific_light_absorption,surface_flux
!EOP
!-----------------------------------------------------------------------
!BOC
   allocate(self)
   call self%initialize(name,parent)

   ! Read the namelist
   read(configunit,nml=bb_passive,err=99,end=100)

   self%surface_flux = surface_flux/secs_pr_day

   ! Register state variables
   call self%register_state_variable(self%id_tracer, &
                    'tracer',unit,'tracer', &
                    initial_concentration,minimum=_ZERO_, &
                    vertical_movement=vertical_velocity/secs_pr_day, &
                    specific_light_extinction=specific_light_absorption)

   return

99 call fatal_error('bb_passive_create','Error reading namelist bb_passive')

100 call fatal_error('bb_passive_create','Namelist bb_passive was not found.')

   end function bb_passive_create
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Air-sea exchange for the passive tracer model
!
! !INTERFACE:
   subroutine get_surface_exchange(self,_FABM_ARGS_GET_SURFACE_EXCHANGE_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   class (type_bb_passive), intent(in)    :: self
   _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Transfer surface exchange value to FABM.
   _SET_SURFACE_EXCHANGE_(self%id_tracer,self%surface_flux)

   ! Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_

   end subroutine get_surface_exchange
!EOC

!-----------------------------------------------------------------------

   end module fabm_bb_passive

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#endif

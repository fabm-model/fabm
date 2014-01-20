#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_iow_age --- age tracer model
!
! !INTERFACE:
   module fabm_iow_age
!
! !DESCRIPTION:
! This model describes an age tracer. Mostly, you track the age of a river  discharge
!  or something specified by lateral boundary conditions. However, you can also 
!  track the age of water masses which had last contact with the surface / bottom.
!  You only have to set track_[surface|bottom]_age.
!
! !USES:
   use fabm_types
   
! !REVISION HISTORY:!
!  Original author(s): Ulf Gräwe (ulf.graewe@io-warnemuende.de)
!
   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model), public :: type_iow_age
!     Variable identifiers
      type (type_state_variable_id) :: id_age

!     Model parameters
      logical                      :: track_surface_age
      logical                      :: track_bottom_age

      contains

      procedure :: initialize
      procedure :: do
      procedure :: check_surface_state
      procedure :: check_bottom_state
      
   end type
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the age tracer model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the iow\_age namelist is read and the variables exported
!  by the model are registered within FABM.
!
! !INPUT PARAMETERS:
   class (type_iow_age), intent(inout), target :: self
   integer,              intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk)                  :: initial_age
   logical                   :: track_surface_age
   logical                   :: track_bottom_age
   character(len=64)         :: units

   namelist /iow_age/     initial_age,track_surface_age,track_bottom_age
!EOP
!-----------------------------------------------------------------------
!BOC
   initial_age       = 0.0_rk
   units             = 'days'
   track_surface_age = .false.
   track_bottom_age  = .false.

   ! Read the namelist
   if (configunit>0) read(configunit,nml=iow_age,err=99,end=100)

   self%track_surface_age = track_surface_age
   self%track_bottom_age  = track_bottom_age

   ! Register state variables
   call self%register_state_variable(self%id_age, &
                    'age',units,'age of water mass', &
                    initial_age,minimum=0.0_rk)

   return

99 call self%fatal_error('iow_age_create','Error reading namelist iow_age')

100 call self%fatal_error('iow_age_create','Namelist iow_age was not found.')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of age tracer
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_iow_age), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk), parameter       :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Set temporal derivatives
   _SET_ODE_(self%id_age,1.0_rk/secs_pr_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: track age of surface water
!
! !INTERFACE:

   subroutine check_surface_state(self,_ARGUMENTS_CHECK_SURFACE_STATE_)
      class (type_iow_age), intent(in) :: self
      _DECLARE_ARGUMENTS_CHECK_SURFACE_STATE_

   if ( self%track_surface_age ) then
      ! set the age in the surface cell to zero
     _HORIZONTAL_LOOP_BEGIN_
        _SET_(self%id_age,0.0_rk)
     _HORIZONTAL_LOOP_END_
     
   endif
     
   end subroutine check_surface_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: track age of bottom water
!
! !INTERFACE:

   subroutine check_bottom_state(self,_ARGUMENTS_CHECK_BOTTOM_STATE_)
      class (type_iow_age), intent(in) :: self
      _DECLARE_ARGUMENTS_CHECK_BOTTOM_STATE_

   if ( self%track_bottom_age ) then
      ! set the age in the bottom cell to zero
     _HORIZONTAL_LOOP_BEGIN_
        _SET_(self%id_age,0.0_rk)
     _HORIZONTAL_LOOP_END_
     
   endif
     
   end subroutine check_bottom_state

!EOC

!-----------------------------------------------------------------------

   end module fabm_iow_age

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------


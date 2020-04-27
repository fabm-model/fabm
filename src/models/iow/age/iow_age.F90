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
!  It is also possible to track the age of a external tracer. Simply specify : tracer_age_variable=???
!  For detail of the age concept, have a look into;
!     England, 1995, Journal of Physical Oceanography, 25, 2756-2777
!     http://www.elic.ucl.ac.be/repomodx/cart/
!     Delhez E.J.M., J.-M. Campin, A.C. Hirst and E. Deleersnijder, 1999, Toward a general
!        theory of the age in ocean modelling, Ocean Modelling, 1, 17-27
!     Deleersnijder E., J.-M. Campin and E.J.M. Delhez, 2001, The concept of age in marine
!        modelling: I. Theory and preliminary model results, Journal of Marine Systems, 28, 229-267
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
      type (type_state_variable_id)       :: id_age
      type (type_state_variable_id)       :: id_age_alpha
      type (type_dependency_id)           :: id_tracer
      type (type_diagnostic_variable_id)  :: id_tracer_age

!     Model parameters
      logical                      :: track_surface_age
      logical                      :: track_bottom_age
      logical                      :: external_tracer

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
   character(len=64)         :: tracer_age_variable=''

   namelist /iow_age/     initial_age,track_surface_age,track_bottom_age, &
                          tracer_age_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   initial_age       = 0.0_rk
   units             = 'days'
   track_surface_age = .false.
   track_bottom_age  = .false.

   ! Read the namelist
   if (configunit>0) read(configunit,nml=iow_age,err=99,end=100)

   call self%get_parameter(self%track_surface_age,'track_surface_age','','track surface age',default=track_surface_age)
   call self%get_parameter(self%track_bottom_age, 'track_bottom_age', '','track bottom age', default=track_bottom_age)
   call self%get_parameter(self%external_tracer, 'external_tracer', '','whether to track age of a tracer rather than of the water parcel', default=tracer_age_variable/='')

   ! Register state variables

   if ( self%external_tracer ) then
      if ( tracer_age_variable=='' ) tracer_age_variable = 'tracer'
      call self%register_dependency(self%id_tracer,trim(tracer_age_variable),'',trim(tracer_age_variable))
      ! we need an auxillery variable
      call self%register_state_variable(self%id_age_alpha, &
                    'age_alpha',' ','age_alpha of water mass', &
                    0.0_rk,minimum=0.0_rk)
      ! for an external tracer, age is derived quantity
      call self%register_diagnostic_variable(self%id_tracer_age, &
                    'age_of_'//trim(tracer_age_variable),units,'age of tracer')
   else
      call self%register_state_variable(self%id_age, &
                    'age_of_water',units,'age of water mass', &
                    initial_age,minimum=0.0_rk)
   endif

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
   real(rk)                  :: tracer
   real(rk)                  :: alpha
   real(rk)                  :: age
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   if ( self%external_tracer ) then
      _GET_(self%id_tracer,tracer)
      _GET_(self%id_age_alpha,alpha)
      _ADD_SOURCE_(self%id_age_alpha,tracer)
      age = 0.0_rk
      if ( tracer .gt. 0.0_rk ) age = alpha/tracer
      _SET_DIAGNOSTIC_(self%id_tracer_age,age/secs_pr_day)
   else
      _ADD_SOURCE_(self%id_age,1.0_rk/secs_pr_day)
   endif

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

!EOP
!-----------------------------------------------------------------------
!BOC
   if ( .not. self%external_tracer ) then
      if ( self%track_surface_age ) then
         ! set the age in the surface cell to zero
        _HORIZONTAL_LOOP_BEGIN_
           _SET_(self%id_age,0.0_rk)
        _HORIZONTAL_LOOP_END_

      endif
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
!EOP
!-----------------------------------------------------------------------
!BOC

   if ( .not. self%external_tracer ) then
      if ( self%track_bottom_age ) then
         ! set the age in the bottom cell to zero
        _HORIZONTAL_LOOP_BEGIN_
           _SET_(self%id_age,0.0_rk)
        _HORIZONTAL_LOOP_END_

      endif
   endif

   end subroutine check_bottom_state

!EOC

!-----------------------------------------------------------------------

   end module fabm_iow_age

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------


#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bb_filter_feeder --- implicit filter feeder at fixed position in pelagic
!
! !INTERFACE:
   module bb_filter_feeder
!
! !DESCRIPTION:
! This model describes the food intake by an implicit predator, placed at
! fixed position in the pelagic. The presence of the predator is prescribed in
! the form of a spatial field of clearance rates. Clearance removes the specified
! prey for the water (the prey field must be a state variable of another biogeochemical
! model running within FABM). The model keeps track of ingested prey by including
! a state variable for time-integrated consumed prey.
!
! Note: the clearance rate is typically calculated as (volume filtered) individual-1 s-1,
! When an external filtration rate is applied it is given as d-1 and during initialisation
! it is converted to s-1.
! multiplied by the concentration of predator (individuals volume-1).
! Thus, the final result has unit s-1 - it can be thought of as the fraction of the water
! volume that is filtered by predators per unit time.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model), public :: type_bb_filter_feeder
      type (type_state_variable_id) :: id_consumed_prey
      type (type_state_variable_id) :: id_prey
      type (type_dependency_id)     :: id_clearance_rate

      logical  :: use_external_clearance_rate
      real(rk) :: clearance_rate

   contains

      procedure :: initialize
      procedure :: do
   end type
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the benthic predator model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the bb\_filter\_feeder type is initialized.
!
! !INPUT PARAMETERS:
   class (type_bb_filter_feeder),intent(inout),target :: self
   integer,                      intent(in)           :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk), parameter       :: days_per_second = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%get_parameter(self%use_external_clearance_rate,'use_external_clearance_rate',' ','use external clearance rate',default=.false.)
   ! Get clearance rate as constant or extenral dependency.
   if (self%use_external_clearance_rate) then
      call self%register_dependency(self%id_clearance_rate,'clearance_rate','s-1','clearance rate')
   else
      call self%get_parameter(self%clearance_rate,'clearance_rate','d-1','constant clearance rate',default=0.0_rk,scale_factor=days_per_second)
   end if

   ! Register state variables
   call self%register_state_variable(self%id_consumed_prey,'consumed_prey','','consumed prey',0.0_rk,minimum=0.0_rk)
   call self%set_variable_property(self%id_consumed_prey,'disable_transport',.true.)

   ! Register link to external pelagic prey.
   call self%register_state_dependency(self%id_prey,'prey')
   return
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of filter_feeder model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
! This routine calculates the sink and source terms.
!
! !INPUT PARAMETERS:
   class (type_bb_filter_feeder), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                   :: prey,clearance_rate
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any).
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_prey,prey)                     ! prey density
   if (self%use_external_clearance_rate) then
      _GET_(self%id_clearance_rate,clearance_rate) ! prescribed clearance rate
   else
      clearance_rate = self%clearance_rate
   end if

   ! Set fluxes of pelagic variables.
   _SET_ODE_(self%id_prey,-prey*clearance_rate)
   _SET_ODE_(self%id_consumed_prey,prey*clearance_rate)

   ! Leave spatial loops (if any).
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------

   end module bb_filter_feeder

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

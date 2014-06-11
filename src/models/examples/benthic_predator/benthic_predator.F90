#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: examples_benthic_predator --- test for benthic interfaces in FABM
!
! !INTERFACE:
   module examples_benthic_predator
!
! !DESCRIPTION:
! This is a very simple model for a benthic predator, grazing according to a
! Monod/Michaelis-Menten functional response on a pelagic prey, and
! respiring/dying according to a linear loss term. Variables for the prey
! (e.g., phytoplankon or zooplankton) and the target for the losses (typically
! a detrital or mineral pool) must be provided by an external model, e.g.,
! gotm\_npzd.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC TYPES:
   type,extends(type_base_model),public :: type_examples_benthic_predator
!     Variable identifiers
      type (type_bottom_state_variable_id) :: id_pred
      type (type_state_variable_id)        :: id_prey,id_nut

!     Model parameters: maximum grazing rate, half-saturation prey density, loss rate
      real(rk) :: g_max,K,h
   contains
      procedure :: initialize
      procedure :: do_bottom
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
!  Here, the examples\_benthic\_predator namelist is read and te variables
!  exported by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_examples_benthic_predator),intent(inout),target :: self
   integer,                               intent(in)           :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                  :: pred_initial=0.01
   real(rk)                  :: g_max = 1., K=1., h=0.05
   character(len=64)         :: waste_target_variable='',prey_source_variable=''

   real(rk), parameter :: secs_pr_day = 86400.
   namelist /examples_benthic_predator/ waste_target_variable,prey_source_variable,pred_initial,g_max,K,h
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   if (configunit>0) read(configunit,nml=examples_benthic_predator,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%g_max = g_max/secs_pr_day
   self%h     = h/secs_pr_day
   self%K     = K

   ! Register state variables
   ! NOTE the benthic=.true. argument, which specifies the variable is benthic.
   call self%register_state_variable(self%id_pred,'pred','mmol/m**2','predator density', &
                                          pred_initial,minimum=0.0_rk)

   ! Register link to external pelagic prey and mineral pools.
   ! Prey will be used to feed upon, mineral pool to place waste products in.
   call self%register_state_dependency(self%id_prey,prey_source_variable)
   call self%register_state_dependency(self%id_nut,waste_target_variable)

   return

99 call self%fatal_error('examples_benthic_predator_init','Error reading namelist examples_benthic_predator')
100 call self%fatal_error('examples_benthic_predator_init','Namelist examples_benthic_predator was not found')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
! This routine calculates the benthic sink and source terms, as well as
! (matching) bottom fluxes for pelagic variables. Both have units mmol/m**2/s.
!
! !INPUT PARAMETERS:
   class (type_examples_benthic_predator),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                   :: prey,pred,g
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_prey,prey)     ! prey density - pelagic
   _GET_HORIZONTAL_(self%id_pred,pred) ! predator density - benthic

   ! Calculate grazing rate
   g = self%g_max*pred*prey/(prey+self%K)

   ! Set local temporal derivatives of benthic variables
   _SET_ODE_BEN_(self%id_pred,g-self%h*pred)

   ! Set bottom fluxes of pelagic variables (these mirror local benthic derivatives)
   _SET_BOTTOM_EXCHANGE_(self%id_prey,-g)
   _SET_BOTTOM_EXCHANGE_(self%id_nut,self%h*pred)

   ! Leave spatial loops over the horizontal domain (if any).
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

!-----------------------------------------------------------------------

   end module examples_benthic_predator

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

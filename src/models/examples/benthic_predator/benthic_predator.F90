!$Id: benthic_predator.F90 119 2010-12-27 14:23:18Z jornbr $
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_benthic_predator --- test for benthic interfaces in FABM
!
! !INTERFACE:
   module fabm_benthic_predator
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   use fabm_driver
   
   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_benthic_predator, benthic_predator_init, benthic_predator_do_benthos, benthic_predator_get_conserved_quantities
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type type_benthic_predator
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_prey,id_pred,id_nut
      _TYPE_CONSERVED_QUANTITY_ID_  :: id_totmass
      
!     Model parameters
      REALTYPE :: g_max,K,h
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the benthic_predator model
!
! !INTERFACE:
   subroutine benthic_predator_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the benthic_predator namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_benthic_predator),intent(out)   :: self
   type (type_model_info),      intent(inout) :: modelinfo
   integer,                     intent(in)    :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                  :: pred_initial=0.01
   REALTYPE                  :: g_max = 2., K=1., h=.1
   character(len=64)         :: nut_variable='',prey_variable=''

   REALTYPE, parameter :: secs_pr_day = 86400.
   namelist /benthic_predator/ nut_variable,prey_variable,pred_initial,g_max,K,h
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=benthic_predator,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%g_max = g_max/secs_pr_day
   self%h     = h/secs_pr_day
   self%K     = K
   
   ! Register state variables
   !self%id_nut  = register_state_variable(modelinfo,'nut','mmol/m**3', 'nutrient', &
   !                                 prey_initial,minimum=_ZERO_,no_river_dilution=.true.)
   !self%id_prey = register_state_variable(modelinfo,'prey','mmol/m**3','prey', &
   !                                 prey_initial,minimum=_ZERO_)
   self%id_pred = register_state_variable(modelinfo,'pred','mmol/m**3','benthic predator', &
                                    pred_initial,minimum=_ZERO_,benthic=.true.)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   self%id_nut  = register_state_dependency(modelinfo,nut_variable)
   self%id_prey = register_state_dependency(modelinfo,prey_variable)

   ! Register conserved quantities
   self%id_totmass = register_conserved_quantity(modelinfo,'mass','mmol/m**3','mass')

   return

99 call fatal_error('benthic_predator_init','Error reading namelist benthic_predator')
   
   end subroutine benthic_predator_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of benthic_predator model
!
! !INTERFACE:
   subroutine benthic_predator_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   type (type_benthic_predator),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: prey,pred,g
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_HZ_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_prey,prey) ! prey density
   _GET_STATE_BEN_(self%id_pred,pred) ! predator density
   
   ! Calculate grazing rate
   g = self%g_max*pred*prey/(prey+self%K)

   ! Set local temporal derivatives of benthic variables
   _SET_ODE_BEN_(self%id_pred,g-self%h*pred)
   
   ! Set bottom fluxes of pelagic variables (these mirror local benthic derivatives)
   _SET_BOTTOM_EXCHANGE_(self%id_prey,-g)
   _SET_BOTTOM_EXCHANGE_(self%id_nut,self%h*pred)
   
   ! Leave spatial loops (if any)
   _FABM_HZ_LOOP_END_

   end subroutine benthic_predator_do_benthos
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   pure subroutine benthic_predator_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   type (type_benthic_predator), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: nut,prey
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_nut, nut)
   _GET_STATE_(self%id_prey,prey)
   
   ! Total nutrient is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totmass,nut+prey)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine benthic_predator_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------

   end module fabm_benthic_predator

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

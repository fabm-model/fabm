#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
!  fabm_prey_pred --- a simple example of biological interaction between
!  two conponents.based on http://www.scholarpedia.org/article/Predator-prey_model
!  and adapted for FABM,implemented by  Fenjuan Hu(Joint PhD stududent of
!  Southern Denmark University&Aarhus university),instructed by Karsten Bolding July,2013.
!
! !INTERFACE:
   module au_prey_predator
!
! !DESCRIPTION: au_prey_predator --- a simple example of biological process, one is uptaken by another
!  the au_prey_predator model described the interactions between two biological communities, especially
!  in the case that one is uptaken by another. Here we implemented two model approches for the prey-
!  predator interaction:
!  Lotka_Volterra model(L_V): explains the smiple and ideal interaction between the prey and predator. In this
!  model, the population of the prey and the predator is only influenced by each other and its own growth rate
!  Jacob-Monod model(J_M): it is manily about the nutrients being teken up by the phytoplankton, bacteria and so on.
!  The use of this model is to describe a continuous-flow growth device, such as a chemostat,
!  where there is continuous removal of nutrient and feeders and a continuous supply of fresh nutirent.
!  USES:
   use fabm_types

   implicit none

   private
!  PUBLIC DERIVED TYPE
   type,extends(type_base_model),public :: type_au_prey_predator
!  Variable identifiers
!  id_prey:     the component which is uptaken by predator;
!  id_predator: the predator uptaking the prey
   type (type_state_variable_id)   :: id_prey,id_predator
!
!  Model parameters
!  V,K,Y are the parameters in Jacob-Monod model:V-the uptake velocity, d-1, for different nutrient have different range; K-saturation constant;Y-the yield of predator per prey taken up
!  b,p,d,r are the parameters in Lotka-Volterra model: b is the natural growth rate of prey in the absence of predation;p measures the impact of predation on prey;
!  d is the death (or emigration) rate of predator in the absence of interaction with prey; r is the efficiency of turning predated prey into pred.
   integer         :: model_type
   real(rk)        :: V,K,Y   ! Jacob-Monod
   real(rk)        :: b,p,r,d ! Lotka-Volterra
!
   contains
!
! Model procedures
   procedure :: initialize
   procedure :: do

   end type type_au_prey_predator

!  private data members
   real(rk),parameter :: secs_pr_day=86400.0_rk
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the prey and predator variables and parameters
!
! !INTERFACE:
!
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the au_prey_predator namelist is read and the variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_au_prey_predator), intent(inout), target :: self
   integer,                          intent(in)         :: configunit
!
! !LOCAL VARIABLES:
   integer  :: model_type=1
   real(rk) :: prey_initial=1.0_rk
   real(rk) :: pred_initial=1.0_rk
   real(rk) :: V=1.0_rk
   real(rk) :: K=1.0_rk
   real(rk) :: Y=0.5_rk
   real(rk) :: b=1.0_rk
   real(rk) :: p=0.05_rk
   real(rk) :: r=0.02_rk
   real(rk) :: d=0.5_rk
!
   namelist /au_prey_predator/ model_type,prey_initial,pred_initial, &
                               V,K,Y,b,p,r,d
!EOP
!-----------------------------------------------------------------------
!BOC
!
   if (configunit>0) read(configunit,nml=au_prey_predator,err=99,end=100)
!
!  Store parameter values in our own derived type
   self%model_type=model_type
   self%V=V/secs_pr_day
   self%K=K
   self%Y=Y
   self%b=b/secs_pr_day
   self%p=p/secs_pr_day
   self%r=r/secs_pr_day
   self%d=d/secs_pr_day

!  Register state variables
   call self%register_state_variable(self%id_prey,'prey','mmol/m**3','nutrient',     &
                                    prey_initial,minimum=0.0_rk,no_river_dilution=.FALSE.)
   call self%register_state_variable(self%id_predator,'predator','mmol/m**3','phytoplankton',     &
                                    pred_initial,minimum=0.0_rk,no_river_dilution=.FALSE.)

!  Register conserved quantities
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_prey)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_predator)

   return

99 call self%fatal_error('au_prey_predator','Error reading namelist au_prey_predator')

100 call self%fatal_error('au_prey_predator','Namelist au_prey_predator was not found')

   end subroutine initialize

!-----------------------------------------------------------------------
! !IROUTINE:the type bound precedure: do(),right hand sides of prey and predator model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_au_prey_predator),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
!  state variables
   real(rk)   :: prey,predator
!  g is changing of predator
!  f is changing of prey
   real(rk)   :: g, f
!EOP
!-----------------------------------------------------------------------
!BOC
! Enter spatial loops (if any)
   _LOOP_BEGIN_

!  Retrieve current state variables values
!  prey density - pelagic
   _GET_(self%id_prey,prey)
!  predator density - pelagic
   _GET_(self%id_predator,predator)

   if (self%model_type .eq. 1) then
!     Use Jacob-Monod model
!     Calculate change of predator
      g = self%V*predator*prey/(prey+self%K)
!     Calculate change of prey
      f = -g/self%Y
   else if (self%model_type .eq. 2) then
!     Use Lotka-Volterra model
!     Calculate change of predator
      g = (self%r*prey-self%d)*predator
!     Calculate change of prey
      f = (self%b-self%p*predator)*prey
   end if
!
!  Set temporal derivatives,  subtracted by the secs_pr_day
   _SET_ODE_(self%id_predator,g)
   _SET_ODE_(self%id_prey,f)
!  Export diagnostic variables
!  Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------

end module au_prey_predator

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

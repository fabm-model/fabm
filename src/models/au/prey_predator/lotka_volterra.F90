#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module au_pp_lotka_volterra
!
! !DESCRIPTION:
!  Lotka_Volterra model(L_V): explains the simple and ideal interaction
!  between the prey and predator. In this model, the population of the
!  prey and the predator is only influenced by each other and its own
!  growth rate.
!
!  USES:
   use fabm_types

   implicit none

   private

!  PUBLIC DERIVED TYPE
   type,extends(type_base_model),public :: type_au_pp_lotka_volterra
!     Variable identifiers
!     id_prey:     the component which is uptaken by predator;
!     id_predator: the predator uptaking the prey
      type (type_state_variable_id)   :: id_prey,id_predator
!
!     Model parameters
      real(rk)        :: b ! natural growth rate of prey
      real(rk)        :: p ! impact of predation
      real(rk)        :: r ! efficiency
      real(rk)        :: d ! death rate
!
      contains
!
! Model procedures
      procedure :: initialize
      procedure :: do

   end type type_au_pp_lotka_volterra

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
!
!
! !INPUT PARAMETERS:
   class (type_au_pp_lotka_volterra), intent(inout), target :: self
   integer,                           intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
!   Register model parameters
   call self%get_parameter(self%b,'b','d-1','growth rate of prey',   default=1.00_rk,scale_factor=d_per_s)
   call self%get_parameter(self%p,'p','d-1','impact of predation',   default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r,'r','d-1','growth efficiency rate',default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%d,'d','d-1','death rate',            default=0.50_rk,scale_factor=d_per_s)

!  Register state variables
   call self%register_state_variable(self%id_prey,'prey','mmol/m**3','nutrient',     &
                                    0._rk,minimum=0.0_rk,no_river_dilution=.FALSE.)
   call self%register_state_variable(self%id_predator,'predator','mmol/m**3','phytoplankton',     &
                                    0._rk,minimum=0.0_rk,no_river_dilution=.FALSE.)

!  Register conserved quantities
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_prey)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_predator)

   return
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
   class (type_au_pp_lotka_volterra),intent(in) :: self
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

!  Use Lotka-Volterra model
!  Calculate change of predator
   g = (self%r*prey-self%d)*predator
!  Calculate change of prey
   f = (self%b-self%p*predator)*prey
!
!  Set temporal derivatives,  subtracted by the secs_pr_day
   _ADD_SOURCE_(self%id_predator,g)
   _ADD_SOURCE_(self%id_prey,f)
!  Export diagnostic variables
!  Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------

end module au_pp_lotka_volterra

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

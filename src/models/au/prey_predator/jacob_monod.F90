#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module au_pp_jacob_monod
!
! !DESCRIPTION:
!  Jacob-Monod model(J_M): it is mainly about the nutrients being taken
!  up by the phytoplankton, bacteria and so on.
!  The use of this model is to describe a continuous-flow growth device,
!  such as a chemostat, where there is continuous removal of nutrient
!  and feeders and a continuous supply of fresh nutirent.
!
!  USES:
   use fabm_types

   implicit none

   private

!  PUBLIC DERIVED TYPE
   type,extends(type_base_model),public :: type_au_pp_jacob_monod
!     Variable identifiers
!     id_prey:     the component which is uptaken by predator;
!     id_predator: the predator uptaking the prey
      type (type_state_variable_id)   :: id_prey,id_predator
!
!     Model parameters
      real(rk)        :: V   ! uptake velocity
      real(rk)        :: K   ! saturation constant
      real(rk)        :: Y   ! predator yield
!
      contains
!
! Model procedures
      procedure :: initialize
      procedure :: do

   end type type_au_pp_jacob_monod

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
!  Here, parameter values are read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_au_pp_jacob_monod), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Register model parameters
   call self%get_parameter(self%V,'V','d-1',      'uptake velocity',default=1.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%K,'K','mmol/m**3','saturation constant', default=1.0_rk)
   call self%get_parameter(self%Y,'Y','mmol/m**3','predator yield',default=0.5_rk)

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
   class (type_au_pp_jacob_monod),intent(in) :: self
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

!  Use Jacob-Monod model
!  Calculate change of predator
   g = self%V*predator*prey/(prey+self%K)
!  Calculate change of prey
   f = -g/self%Y
!
!  Set temporal derivatives,  subtracted by the secs_pr_day
   _ADD_SOURCE_(self%id_predator,g)
   _ADD_SOURCE_(self%id_prey,f)
!  Export diagnostic variables
!  Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------

end module au_pp_jacob_monod

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

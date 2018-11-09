#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module bb_lorenz63
!
! !DESCRIPTION:
!
!  USES:
   use fabm_types
!
   implicit none
!
   private
!
!  PUBLIC DERIVED TYPE
   type,extends(type_base_model),public :: type_bb_lorenz63
!     Variable identifiers
      type (type_state_variable_id)   :: id_x, id_y, id_z
!
!     Model parameters
      real(rk)        :: sigma !
      real(rk)        :: rho   !
      real(rk)        :: beta  !
!
      contains
!
! Model procedures
      procedure :: initialize
      procedure :: do
!
   end type type_bb_lorenz63
!
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
   class (type_bb_lorenz63), intent(inout), target :: self
   integer,                  intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk)        :: beta
!EOP
!-----------------------------------------------------------------------
!BOC
!   Register model parameters
   beta = 8._rk/3._rk
   call self%get_parameter(self%sigma,'sigma','','',default=10.0_rk)
   call self%get_parameter(self%rho,  'rho'  ,'','',default=28._rk)
   call self%get_parameter(self%beta, 'beta' ,'','',default=beta)

!  Register state variables
   call self%register_state_variable(self%id_x,'x','','',0._rk,minimum=0.0_rk)
   call self%register_state_variable(self%id_y,'y','','',0._rk,minimum=0.0_rk)
   call self%register_state_variable(self%id_z,'z','','',0._rk,minimum=0.0_rk)

!  Register conserved quantities
!   call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_predator)

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
   class (type_bb_lorenz63),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
!  state variables
   real(rk)   :: x, y, z
!EOP
!-----------------------------------------------------------------------
!BOC
! Enter spatial loops (if any)
   _LOOP_BEGIN_

!  Retrieve current state variables values
   _GET_(self%id_x,x)
   _GET_(self%id_y,y)
   _GET_(self%id_z,z)

!  Set temporal derivatives
   _SET_ODE_(self%id_x, self%sigma*(y-x))
   _SET_ODE_(self%id_y, x*(self%rho-z)-y)
   _SET_ODE_(self%id_z, x*y-self%beta*z)

!  Export diagnostic variables
!  Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------

end module bb_lorenz63

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

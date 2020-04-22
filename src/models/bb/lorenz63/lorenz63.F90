#include "fabm_driver.h"

module bb_lorenz63

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_bb_lorenz63
      ! Variable identifiers
      type (type_state_variable_id) :: id_x, id_y, id_z

      ! Model parameters
      real(rk) :: sigma
      real(rk) :: rho
      real(rk) :: beta
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type type_bb_lorenz63

contains

   subroutine initialize(self, configunit)
      class (type_bb_lorenz63), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      ! Register model parameters
      call self%get_parameter(self%sigma, 'sigma', '-', 'sigma', default=10.0_rk)
      call self%get_parameter(self%rho,   'rho'  , '-', 'rho',   default=28._rk)
      call self%get_parameter(self%beta,  'beta' , '-', 'beta',  default=8._rk / 3._rk)

      ! Register state variables
      call self%register_state_variable(self%id_x, 'x', '-', 'x', 0._rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_y, 'y', '-', 'y', 0._rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_z, 'z', '-', 'z', 0._rk, minimum=0.0_rk)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_bb_lorenz63), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: x, y, z

      _LOOP_BEGIN_

         ! Retrieve current state variables values
         _GET_(self%id_x, x)
         _GET_(self%id_y, y)
         _GET_(self%id_z, z)

         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_x, self%sigma*(y-x))
         _ADD_SOURCE_(self%id_y, x*(self%rho-z)-y)
         _ADD_SOURCE_(self%id_z, x*y-self%beta*z)
      
      _LOOP_END_
   end subroutine do

end module bb_lorenz63

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

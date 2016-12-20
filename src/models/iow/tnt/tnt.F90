#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_iow_tnt
!
! !INTERFACE:
   module fabm_iow_tnt
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_iow_tnt
      type (type_state_variable_id)                   :: id_conc
      type (type_dependency_id)                       :: id_temp

      real(rk) :: k
      real(rk) :: E_a
      real(rk) :: T_ref
   contains
      procedure :: initialize
      procedure :: do
   end type type_iow_tnt
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise the sediment model
!
! !INTERFACE:
   subroutine initialize(self,configunit)

   integer,              intent(in)            :: configunit
   class (type_iow_tnt), intent(inout), target :: self

   real(rk) :: w
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   call self%get_parameter(self%k,     'k',     'd-1',      'decay rate (half-life = ln(2)/k)', default=0.0_rk, scale_factor=1.0_rk/86400)
   call self%get_parameter(self%E_a,   'E_a',   'J mol-1',  'activation energy for decay (0: no temperature dependence)', minimum=0.0_rk)
   call self%get_parameter(self%T_ref, 'T_ref', 'degree_C', 'temperature for which decay rate k is given')
   call self%get_parameter(w,          'w',     'd-1',      'vertical velocity (<0 for sinking)', default=0.0_rk, scale_factor=1.0_rk/86400)

   call self%register_dependency(self%id_temp, standard_variables%temperature)

   call self%register_state_variable(self%id_conc, 'conc', '', 'concentration', vertical_movement=w)
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of spm model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class(type_iow_tnt), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_

  real(rk), parameter :: R = 8.3144598_rk    ! universal gas constant (J mol-1 K-1)
  real(rk), parameter :: Kelvin = 273.15_rk  ! offset between degrees Celsius and Kelvin

  real(rk) :: T, c, k
!EOP
!-----------------------------------------------------------------------
!BOC
   _LOOP_BEGIN_
      _GET_(self%id_temp, T)
      _GET_(self%id_conc, c)
      k = self%k * exp(self%E_a/(R*(Kelvin+self%T_ref))-self%E_a/(R*(Kelvin+T)))
      _SET_ODE_(self%id_conc, -k*c)
   _LOOP_END_

   end subroutine do
!EOC

end module fabm_iow_tnt

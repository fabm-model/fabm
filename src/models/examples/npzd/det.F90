#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - detritus component
! This model features a single detritus variable, characterized by a rate of decay (rdn)
! and a sinking rate. Mineralized detritus feeds into a dissolved mineral pool that must
! be provided by an external model (e.g., examples_npzd_nut).

module examples_npzd_det
   use fabm_types

   implicit none

   private

   type, extends(type_base_model),public :: type_examples_npzd_det
      ! Variable identifiers
      type (type_state_variable_id)     :: id_d
      type (type_state_variable_id)     :: id_mintarget

      ! Model parameters
      real(rk) :: rdn
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
   end type

contains

   subroutine initialize(self, configunit)
      class (type_examples_npzd_det), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_d, kc

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(w_d,      'w_d', 'm d-1',     'vertical velocity (<0 for sinking)', default=-5.0_rk, scale_factor=d_per_s)
      call self%get_parameter(kc,       'kc',  'm2 mmol-1', 'specific light extinction',          default=0.03_rk)
      call self%get_parameter(self%rdn, 'rdn', 'd-1',       'remineralization rate',              default=0.003_rk, scale_factor=d_per_s)

      ! Register state variables
      call self%register_state_variable(self%id_d, 'c','mmol m-3',  'concentration', 4.5_rk, &
         minimum=0.0_rk, vertical_movement=w_d, specific_light_extinction=kc)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_d)

      ! Register dependencies on external state variables
      call self%register_state_dependency(self%id_mintarget, 'mineralisation_target', 'mmol m-3', 'sink for remineralized matter')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_examples_npzd_det), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: d

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_d, d) ! detritus

         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_d, -self%rdn*d)
         _ADD_SOURCE_(self%id_mintarget, self%rdn*d)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

   subroutine do_ppdd(self, _ARGUMENTS_DO_PPDD_)
      class (type_examples_npzd_det), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_PPDD_

      real(rk) :: d

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_d, d) ! detritus

         ! Assign destruction rates to different elements of the destruction matrix.
         ! By assigning with _SET_DD_SYM_ [as opposed to _SET_DD_], assignments to dd(i,j)
         ! are automatically assigned to pp(j,i) as well.
         _SET_DD_SYM_(self%id_d, self%id_mintarget, self%rdn*d)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do_ppdd

end module examples_npzd_det

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

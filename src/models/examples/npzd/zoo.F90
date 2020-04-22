#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - zooplankton component

module examples_npzd_zoo
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_examples_npzd_zoo
      ! Variable identifiers
      type (type_state_variable_id) :: id_z
      type (type_state_variable_id) :: id_exctarget, id_morttarget, id_grztarget

      ! Model parameters
      real(rk) :: z0, gmax, iv, rzn, rzd
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
   end type

contains

   subroutine initialize(self, configunit)
      class (type_examples_npzd_zoo), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%z0,   'z0',   'mmol m-3',  'background concentration',      default=0.0225_rk)
      call self%get_parameter(self%gmax, 'gmax', 'd-1',       'maximum specific grazing rate', default=0.5_rk,  scale_factor=d_per_s)
      call self%get_parameter(self%iv,   'iv',   'm3 mmol-1', 'Ivlev grazing constant',        default=1.1_rk)
      call self%get_parameter(self%rzn,  'rzn',  'd-1',       'excretion rate',                default=0.01_rk, scale_factor=d_per_s)
      call self%get_parameter(self%rzd,  'rzd',  'd-1',       'mortality',                     default=0.02_rk, scale_factor=d_per_s)

      ! Register state variables
      call self%register_state_variable(self%id_z, 'c', 'mmol m-3', 'concentration', 0.0_rk, minimum=0.0_rk)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_z)

      ! Register dependencies on external state variables.
      call self%register_state_dependency(self%id_grztarget,  'grazing_target',   'mmol m-3', 'prey source')
      call self%register_state_dependency(self%id_exctarget,  'excretion_target', 'mmol m-3', 'sink for excreted matter')
      call self%register_state_dependency(self%id_morttarget, 'mortality_target', 'mmol m-3', 'sink for dead matter')
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_examples_npzd_zoo), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: p, z, g

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_z, z)         ! zooplankton
         _GET_(self%id_grztarget, p) ! prey

         ! Grazing rate
         g = fpz(self%gmax, self%iv, p, z + self%z0)

         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_z, g - self%rzn*z - self%rzd*z)
         _ADD_SOURCE_(self%id_grztarget, -g)
         _ADD_SOURCE_(self%id_morttarget, self%rzd*z)
         _ADD_SOURCE_(self%id_exctarget, self%rzn*z)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

   subroutine do_ppdd(self, _ARGUMENTS_DO_PPDD_)
      class (type_examples_npzd_zoo), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_PPDD_

      real(rk) :: p, z

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_z,z)         ! zooplankton
         _GET_(self%id_grztarget,p) ! prey

         ! Assign destruction rates to different elements of the destruction matrix.
         ! By assigning with _SET_DD_SYM_ [as opposed to _SET_DD_], assignments to dd(i,j)
         ! are automatically assigned to pp(j,i) as well.
          _SET_DD_SYM_(self%id_grztarget,self%id_z,fpz(self%gmax,self%iv,p,z+self%z0))
          _SET_DD_SYM_(self%id_z,self%id_exctarget,self%rzn*z)
          _SET_DD_SYM_(self%id_z,self%id_morttarget,self%rzd*z)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do_ppdd

   ! Ivlev formulation for zooplankton grazing on phytoplankton
   elemental real(rk) function fpz(gmax, iv, p, z)
      real(rk), intent(in) :: gmax, iv, p, z

      fpz = gmax * (1.0_rk - exp(-iv * iv * p * p)) * z
   end function fpz

end module examples_npzd_zoo

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

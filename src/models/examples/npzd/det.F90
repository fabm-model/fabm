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
      type (type_state_variable_id) :: id_c
      type (type_state_variable_id) :: id_mintarget
      type (type_bottom_state_variable_id) :: id_sedtarget

      ! Model parameters
      real(rk) :: rdn, w_d
      logical :: sedimentation
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_examples_npzd_det), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: kc

      ! Store parameter values in our own derived type
      ! All rates must be provided in values per day in fabm.yaml.
      ! They are converted here to values per second (scale_factor argument) to ensure
      ! that the sources and sinks calculated from these parameters will be in per second as FABM expects.
      call self%get_parameter(self%w_d, 'w_d', 'm d-1',     'vertical velocity (<0 for sinking)', default=-5.0_rk, scale_factor=d_per_s)
      call self%get_parameter(kc,       'kc',  'm2 mmol-1', 'specific light extinction',          default=0.03_rk)
      call self%get_parameter(self%rdn, 'rdn', 'd-1',       'remineralization rate',              default=0.003_rk, scale_factor=d_per_s)
      call self%get_parameter(self%sedimentation, 'sedimentation', '', 'transfer into sediment pool at the bed', default=.false.)

      ! Register state variables
      call self%register_state_variable(self%id_c, 'c','mmol m-3', 'concentration', initial_value=4.5_rk, &
         minimum=0.0_rk, vertical_movement=self%w_d, specific_light_extinction=kc)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c)

      ! Register dependencies on external state variables
      call self%register_state_dependency(self%id_mintarget, 'mineralisation_target', 'mmol m-3', 'sink for remineralized matter')
      if (self%sedimentation) call self%register_state_dependency(self%id_sedtarget, 'sedimentation_target', 'mmol m-2', 'sink for sedimented matter')

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_examples_npzd_det), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: c

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_c, c) ! detritus

         ! Local source terms
         _ADD_SOURCE_(self%id_c, -self%rdn*c)
         _ADD_SOURCE_(self%id_mintarget, self%rdn*c)

      ! Leave spatial loops (if any)
      _LOOP_END_

   end subroutine do

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_examples_npzd_det), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c

      if (.not. self%sedimentation) return

      ! Enter spatial loops over the horizontal domain (if any).
      _BOTTOM_LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_c, c)

         ! Bottom fluxes of pelagic variables
         _ADD_BOTTOM_FLUX_(self%id_c, self%w_d * c)

         ! Local source terms of benthic variables
         _ADD_BOTTOM_SOURCE_(self%id_sedtarget, -self%w_d * c)

      ! Leave spatial loops over the horizontal domain (if any).
      _BOTTOM_LOOP_END_

   end subroutine do_bottom

end module examples_npzd_det

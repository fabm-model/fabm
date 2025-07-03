#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module examples_npzd_phy
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_examples_npzd_phy
      ! Variable identifiers
      type (type_state_variable_id)      :: id_c
      type (type_state_variable_id)      :: id_exctarget, id_morttarget, id_upttarget
      type (type_dependency_id)          :: id_par
      type (type_surface_dependency_id)  :: id_I_0
      type (type_diagnostic_variable_id) :: id_PPR, id_NPR, id_dPAR

      ! Model parameters
      real(rk) :: p0, z0, kc, i_min, rmax, gmax, iv, alpha, rpn, rpdu, rpdl
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_examples_npzd_phy), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      ! Store parameter values in our own derived type
      ! All rates must be provided in values per day in fabm.yaml.
      ! They are converted here to values per second (scale_factor argument) to ensure
      ! that the sources and sinks calculated from these parameters will be in per second as FABM expects.
      call self%get_parameter(self%p0,    'p0',    'mmol m-3',  'background concentration ',                default=0.0225_rk)
      call self%get_parameter(self%kc,    'kc',    'm2 mmol-1', 'specific light extinction',                default=0.03_rk)
      call self%get_parameter(self%i_min, 'i_min', 'W m-2',     'minimum light intensity in euphotic zone', default=25.0_rk)
      call self%get_parameter(self%rmax,  'rmax',  'd-1',       'maximum specific growth rate',             default=1.0_rk,  scale_factor=d_per_s)
      call self%get_parameter(self%alpha, 'alpha', 'mmol m-3',  'half-saturation nutrient concentration',   default=0.3_rk)
      call self%get_parameter(self%rpn,   'rpn',   'd-1',       'excretion rate',                           default=0.01_rk, scale_factor=d_per_s)
      call self%get_parameter(self%rpdu,  'rpdu',  'd-1',       'mortality in euphotic zone',               default=0.02_rk, scale_factor=d_per_s)
      call self%get_parameter(self%rpdl,  'rpdl',  'd-1',       'mortality below euphotic zone',            default=0.1_rk,  scale_factor=d_per_s)
      call self%get_parameter(w_p,        'w_p',   'm d-1',     'vertical velocity (<0 for sinking)',       default=-1.0_rk, scale_factor=d_per_s)

      ! Register state variables
      call self%register_state_variable(self%id_c, 'c', 'mmol m-3', 'concentration', &
         initial_value=0.0_rk, minimum=0.0_rk, vertical_movement=w_p)

      ! Register contribution of state to global aggregate variables.
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_c)

      ! Register dependencies on external state variables
      call self%register_state_dependency(self%id_upttarget,  'uptake_target',    'mmol m-3', 'nutrient source')
      call self%register_state_dependency(self%id_exctarget,  'excretion_target', 'mmol m-3', 'sink for excreted matter')
      call self%register_state_dependency(self%id_morttarget, 'mortality_target', 'mmol m-3', 'sink for dead matter')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_PPR,  'PPR', 'mmol m-3 d-1', 'gross primary production rate')
      call self%register_diagnostic_variable(self%id_NPR,  'NPR', 'mmol m-3 d-1', 'net community production rate')
      call self%register_diagnostic_variable(self%id_dPAR, 'PAR', 'W m-2',        'photosynthetically active radiation')

      ! Register environmental dependencies
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

      ! Contribute to light attentuation
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_c, scale_factor=self%kc)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%p0 * self%kc)
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_examples_npzd_phy), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: n, c, par, I_0
      real(rk)            :: iopt, rpd, primprod
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_c, c)         ! phytoplankton concentration
         _GET_(self%id_upttarget, n) ! nutrients

         ! Retrieve current environmental conditions.
         _GET_(self%id_par, par)         ! local photosynthetically active radiation
         _GET_SURFACE_(self%id_I_0, I_0) ! surface photosynthetically active radiation

         ! Light acclimation formulation based on surface light intensity.
         iopt = max(0.25*I_0, self%I_min)

         ! Loss rate of phytoplankton to detritus depends on local light intensity.
         if (par >= self%I_min) then
            rpd = self%rpdu
         else
            rpd = self%rpdl
         end if

         ! Phytoplankton growth limited by light and nutrient availability (multiplicative)
         primprod = self%rmax * par / iopt * exp(1.0_rk - par / iopt) * n / (self%alpha + n) * (c + self%p0)

         ! Set temporal derivatives
         _ADD_SOURCE_(self%id_c, primprod - self%rpn*c - rpd*c)
         _ADD_SOURCE_(self%id_upttarget, -primprod)
         _ADD_SOURCE_(self%id_morttarget, rpd*c)
         _ADD_SOURCE_(self%id_exctarget, self%rpn*c)

         ! Diagnostics
         _SET_DIAGNOSTIC_(self%id_dPAR, par)
         _SET_DIAGNOSTIC_(self%id_PPR, primprod*secs_pr_day)
         _SET_DIAGNOSTIC_(self%id_NPR, (primprod - self%rpn*c)*secs_pr_day)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

end module examples_npzd_phy

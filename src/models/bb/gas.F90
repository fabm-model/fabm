#include "fabm_driver.h"

module bb_gas

   ! OMIP protocol for inert gas as defined in https://doi.org/10.5194/gmd-10-2169-2017

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_bb_gas
      type (type_state_variable_id)              :: id_c
      type (type_dependency_id)                  :: id_temp, id_salt
      type (type_surface_dependency_id)          :: id_patm, id_wind, id_ice, id_x_A
      type (type_surface_diagnostic_variable_id) :: id_F
      real(rk) :: A, B, C, D, E
      real(rk) :: a1, a2, a3, a4, b1, b2, b3
   contains
      procedure initialize
      procedure do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_bb_gas), intent(inout), target :: self
      integer,             intent(in)            :: configunit

      call self%get_parameter(self%A, 'A', '', 'coefficient A for Schmidt number polynomial (Table 1 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%B, 'B', '', 'coefficient B for Schmidt number polynomial (Table 1 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%C, 'C', '', 'coefficient C for Schmidt number polynomial (Table 1 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%D, 'D', '', 'coefficient D for Schmidt number polynomial (Table 1 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%E, 'E', '', 'coefficient E for Schmidt number polynomial (Table 1 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%a1, 'a1', '', 'coefficient a1 for solubility function (Table 2 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%a2, 'a2', '', 'coefficient a2 for solubility function (Table 2 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%a3, 'a3', '', 'coefficient a3 for solubility function (Table 2 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%a4, 'a4', '', 'coefficient a4 for solubility function (Table 2 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%b1, 'b1', '', 'coefficient b1 for solubility function (Table 2 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%b2, 'b2', '', 'coefficient b2 for solubility function (Table 2 in https://doi.org/10.5194/gmd-10-2169-2017)')
      call self%get_parameter(self%b3, 'b3', '', 'coefficient b3 for solubility function (Table 2 in https://doi.org/10.5194/gmd-10-2169-2017)')

      write(*,*) "Schmidt number at 20 degrees Celsius:", Schmidt_number(self, 20.0_rk)

      call self%register_state_variable(self%id_c, 'c', 'mol m-3', 'concentration', 0.0_rk)

      call self%register_diagnostic_variable(self%id_F, 'F', 'mol m-2 s-1', 'air-to-sea flux')

      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
      call self%register_dependency(self%id_wind, standard_variables%wind_speed)
      call self%register_dependency(self%id_patm, standard_variables%surface_air_pressure)
      call self%register_dependency(self%id_ice, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_x_A, 'x_A', '1', 'mole fraction in dry air')
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_bb_gas), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: T_c, S, P, u, Fice, x_A, Csurf
      real(rk) :: Sc, Kw, cTK, log_phi_A0, Csat, F
      real(rk), parameter :: a = 0.251_rk * 0.01_rk / 3600.0_rk   ! Orr et al. p 2178, converting from cm hr-1 to m s-1
      real(rk), parameter :: P0 = 101325.0_rk   ! 1 atm in Pa

      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_temp, T_c)          ! surface temperature (degrees Celsius) - formally must be in-situ
         _GET_(self%id_salt, S)            ! surface practical salinity
         _GET_SURFACE_(self%id_patm, P)    ! surface air pressure (Pa)
         _GET_SURFACE_(self%id_wind, u)    ! wind speed (m s-1)
         _GET_SURFACE_(self%id_ice, Fice)  ! ice area fraction (1)
         _GET_SURFACE_(self%id_x_A, x_A)   ! mole fraction of target gas in dry air (1)
         _GET_(self%id_c, Csurf)           ! surface concentration of target gas

         ! Schmidt number
         Sc = Schmidt_number(self, T_c)

         ! Instantaneous gas transfer velocity (m s-1)
         Kw = a * sqrt(660.0_rk / Sc) * u * u * (1.0 - Fice)

         ! Natural logarithm of solubility function
         cTK = (T_c + 273.15_rk) * 0.01_rk   ! Absolute temperature (K) / 100
         log_phi_A0 = self%a1 + self%a2 / cTK + self%a3 * log(cTK) + self%a4 * cTK * cTK &
            + S * (self%b1 + cTK * (self%b2 + self%b3 * cTK))

         ! Saturation concentration (mol m-3)
         ! NB we multiply with 1000 to convert from mol L-1 to mol m-3
         Csat = P * exp(log_phi_A0) * x_A * (1000.0_rk / P0)

         ! Air -> sea flux (mol m-2 s-1)
         F = Kw * (Csat - Csurf)

         _ADD_SURFACE_FLUX_(self%id_c, F)
         _SET_SURFACE_DIAGNOSTIC_(self%id_F, F)
      _SURFACE_LOOP_END_
   end subroutine do_surface

   elemental function Schmidt_number(self, T_c) result(Sc)
      class (type_bb_gas), intent(in) :: self
      real(rk),            intent(in) :: T_c
      real(rk)                        :: Sc

      Sc = self%A + T_c * (self%B + T_c * (self%C + T_c * (self%D + T_c * self%E)))
   end function

end module bb_gas

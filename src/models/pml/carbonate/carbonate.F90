#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_pml_carbonate --- shell around carbonate chemistry model by
! Jerry Blackford (Plymouth Marine Laboratory), adapted for FABM by Jorn Bruggeman
!
! This code is maintained as example and to support old model setups.
! New simulations should use the updated version of this code distributed with
! ERSEM (http://ersem.com).
!
! !INTERFACE:
module fabm_pml_carbonate
!
! !DESCRIPTION:
! Carbonate chemistry model system model based on PML code.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_pml_carbonate
!     Variable identifiers
      type (type_state_variable_id)              :: id_dic, id_alk
      type (type_dependency_id)                  :: id_temp, id_salt, id_pres, id_dens
      type (type_surface_dependency_id)          :: id_wind, id_pco2_surf
      type (type_diagnostic_variable_id)         :: id_ph, id_pco2, id_CarbA, id_Bicarb, &
                                                    id_Carb, id_Om_cal, id_Om_arg, id_alk_diag
      type (type_surface_diagnostic_variable_id) :: id_co2_flux

!     Model parameters
      real(rk) :: TA_offset, TA_slope, pCO2a
      logical  :: alk_param
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Read parameter values and store settings in the model's derived type.
!
! !INPUT PARAMETERS:
   class (type_pml_carbonate), intent(inout),target :: self
   integer,                    intent(in )          :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day, and are converted here to values per second.
   call self%get_parameter(self%alk_param, 'alk_param', '', 'compute alkalinity as linear function of salinity', default=.true.)
   if (self%alk_param) then
      call self%get_parameter(self%TA_offset, 'alk_offset', 'mEq m-3', 'offset for alkalinity as linear function of salinity', default=520.1_rk)
      call self%get_parameter(self%TA_slope, 'alk_slope', 'mEq m-3', 'scale factor for alkalinity as linear function of salinity', default=51.24_rk)
   end if
   call self%get_parameter(self%pCO2a, 'pCO2a', 'ppm', 'mole fraction of atmospheric CO2', default=0.0_rk)

   ! First state variable: total dissolved inorganic carbon
   call self%register_state_variable(self%id_dic, 'dic', 'mmol m-3', 'total dissolved inorganic carbon', &
                                2185.0_rk, minimum=0.0_rk, no_precipitation_dilution=.false., no_river_dilution=.true., &
                                standard_variable=standard_variables%mole_concentration_of_dissolved_inorganic_carbon)

   if (self%alk_param) then
     ! Alkalinity is diagnosed from temperature and salinity. Register it as output variable.
     call self%register_diagnostic_variable(self%id_alk_diag, 'alk', 'mEq m-3', 'alkalinity')
   else
     ! Alkalinity is a state variable.
     call self%register_state_variable(self%id_alk, 'alk', 'mEq m-3', 'alkalinity', &
                                  2333.0_rk, minimum=0.0_rk, no_precipitation_dilution=.false., no_river_dilution=.true.)
   end if

   ! Register diagnostic variables.
   call self%register_diagnostic_variable(self%id_ph,       'pH',       '-',            'pH'                           )
   call self%register_diagnostic_variable(self%id_pco2,     'pCO2',     'ppm',          'CO2 partial pressure'         )
   call self%register_diagnostic_variable(self%id_CarbA,    'CarbA',    'mmol m-3',     'carbonic acid concentration'  )
   call self%register_diagnostic_variable(self%id_Bicarb,   'Bicarb',   'mmol m-3',     'bicarbonate ion concentration')
   call self%register_diagnostic_variable(self%id_Carb,     'Carb',     'mmol m-3',     'carbonate ion concentration'  )
   call self%register_diagnostic_variable(self%id_Om_cal,   'Om_cal',   '-',            'calcite saturation state'     )
   call self%register_diagnostic_variable(self%id_Om_arg,   'Om_arg',   '-',            'aragonite saturation state'   )
   call self%register_diagnostic_variable(self%id_co2_flux, 'CO2_flux', 'mmol m-2 s-1', 'surface CO2 flux'             )

   ! Register external dependencies.
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
   call self%register_dependency(self%id_pres, standard_variables%pressure)
   call self%register_dependency(self%id_dens, standard_variables%density)
   call self%register_dependency(self%id_wind, standard_variables%wind_speed)
   if (self%pCO2a==0.0_rk) call self%register_dependency(self%id_pco2_surf, standard_variables%mole_fraction_of_carbon_dioxide_in_air)

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of carbonate system model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!  Calculate carbonate system equilibrium from DIC and alkalinity, and
!  store values of all carbonate system parameters (e.g., pCO2, pH,
!  concentrations of inorganic carbon species, calcite and aragonite saturation
!  states.
!
! !INPUT PARAMETERS:
   class (type_pml_carbonate),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   ! Environment
   real(rk) :: temp, salt, pres, dens
   real(rk) :: dic, TA
   real(rk) :: PCO2WATER, pH, HENRY, ca, bc, cb, Om_cal, Om_arg
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Get environmental variables.
   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_(self%id_pres,pres)
   _GET_(self%id_dens,dens)

   ! Get current value for total dissolved inorganic carbon (our own state variable).
   _GET_(self%id_dic,dic)

   if (self%alk_param) then
      ! Linearly approximate alkalinity (uEq/kg) from salinity.
      TA = self%TA_offset + self%TA_slope*salt
   else
      ! Alkalinity (mEq/m**3) is a separate state variable.
      ! Divide by density/1000 to get alkalinity in uEq/kg, as expected by underlying carbonate system model.
      _GET_(self%id_alk,TA)
      TA = TA/dens*1.0e3_rk
   end if

   ! Calculate carbonate system equilibrium.
   call CO2DYN(dic/1.0e3_rk/dens, TA/1.0e6_rk, temp, salt, PCO2WATER, pH, HENRY, ca, bc, cb)

   ! Calculate calcite and aragonite calcification states.
   call CaCO3_Saturation (temp, salt, pres, cb, Om_cal, Om_arg)

   ! Store diagnostic variables.
   _SET_DIAGNOSTIC_(self%id_ph    ,ph)
   _SET_DIAGNOSTIC_(self%id_pco2  ,PCO2WATER*1.0e6_rk)        ! to ppm
   _SET_DIAGNOSTIC_(self%id_CarbA ,ca       *1.0e3_rk*dens)   ! from mol/kg to mmol/m**3
   _SET_DIAGNOSTIC_(self%id_Bicarb,bc       *1.0e3_rk*dens)   ! from mol/kg to mmol/m**3
   _SET_DIAGNOSTIC_(self%id_Carb  ,cb       *1.0e3_rk*dens)   ! from mol/kg to mmol/m**3
   _SET_DIAGNOSTIC_(self%id_Om_cal,Om_cal)
   _SET_DIAGNOSTIC_(self%id_Om_arg,Om_arg)
   if (self%alk_param) _SET_DIAGNOSTIC_(self%id_alk_diag,TA*dens*1.0e-3_rk)    ! from uEg/kg to mmol/m**3

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Air-sea exchange for the carbonate system model
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
! Calculate air -> sea CO2 flux.
!
! !INPUT PARAMETERS:
   class (type_pml_carbonate), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   ! Environment
   real(rk) :: temp, salt, wnd, dens

   ! State
   real(rk) :: dic, TA

   ! Temporary variables
   real(rk) :: PCO2WATER, pH, HENRY, ca, bc, cb, fl, pCO2a

   ! Parameters
   real(rk), parameter :: days_per_sec = 1.0_rk / 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _SURFACE_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_(self%id_dens,dens)
   _GET_SURFACE_(self%id_wind,wnd)
   if (self%pCO2a==0.0_rk) then
      _GET_SURFACE_(self%id_pco2_surf,pCO2a)
   else
      pCO2a = self%pCO2a
   end if

   _GET_(self%id_dic,dic)

   if (self%alk_param) then
      ! Linearly approximate alkalinity (uEq/kg) from salinity.
      TA = self%TA_offset + self%TA_slope*salt
   else
      ! Alkalinity (mEq/m**3) is a separate state variable
      ! Divide by density/1000 to get alkalinity in uEq/kg.
      _GET_(self%id_alk,TA)
      TA = TA/dens*1.0e3_rk
   end if

   ! Calculate carbonate system equilibrium to get pCO2 and Henry constant.
   call CO2DYN(dic/1.0e3_rk/dens, TA/1.0e6_rk, temp, salt, PCO2WATER, pH, HENRY, ca, bc, cb)

   ! Calculate air-sea exchange of CO2 (positive flux is from atmosphere to water)
   call Air_sea_exchange(temp, wnd, PCO2WATER*1.0e6_rk, pCO2a, Henry, dens/1.0e3_rk, fl)

   ! Transfer surface exchange value to FABM.
   _ADD_SURFACE_FLUX_(self%id_dic,fl * days_per_sec)

   ! Store surface flux as diagnostic variable.
   _SET_SURFACE_DIAGNOSTIC_(self%id_co2_flux,fl * days_per_sec)

   ! Leave spatial loops (if any)
   _SURFACE_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------

end module fabm_pml_carbonate

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

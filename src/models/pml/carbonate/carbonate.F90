#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_pml_carbonate --- shell around carbonate chemistry model by
! Jerry Blackford (Plymouth Marine Laboratory), adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
module fabm_pml_carbonate
!
! !DESCRIPTION:
! Carbonate chemistry model system model based on PML code.
!
! !USES:
   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_pml_carbonate, pml_carbonate_init, pml_carbonate_do, pml_carbonate_get_surface_exchange
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type type_pml_carbonate
!     Variable identifiers
      type (type_state_variable_id)      :: id_dic, id_alk
      type (type_dependency_id)          :: id_temp, id_salt, id_pres, id_dens
      type (type_horizontal_dependency_id)  :: id_wind, id_pco2_surf
      type (type_diagnostic_variable_id) :: id_ph, id_pco2, id_CarbA, id_Bicarb, &
                                     & id_Carb, id_Om_cal, id_Om_arg, id_alk_diag
      type (type_horizontal_diagnostic_variable_id) :: id_co2_flux

!     Model parameters
      real(rk) :: TA_offset, TA_slope, pCO2a
      logical  :: alk_param
   end type
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine pml_carbonate_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the bio namelist is read and
!  various variables (rates and settling velocities)
!  are transformed into SI units.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_pml_carbonate), intent(out)   :: self
   _CLASS_ (type_model_info), intent(inout) :: modelinfo
   integer,                   intent(in )   :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk) :: dic_initial, alk_initial
   real(rk)  :: alk_offset = 520.1, alk_slope = 51.24, pCO2a = _ZERO_
   logical :: alk_param = .true.
   namelist /pml_carbonate/ dic_initial, alk_initial, alk_param, alk_offset, alk_slope, pCO2a
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=pml_carbonate,err=99,end=100)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day, and are converted here to values per second.
   self%TA_offset = alk_offset
   self%TA_slope  = alk_slope
   self%alk_param = alk_param
   self%pCO2a     = pCO2a

   ! First state variable: total dissolved inorganic carbon
   call register_state_variable(modelinfo,self%id_dic,'dic','mmol/m**3','total dissolved inorganic carbon', &
                                    dic_initial,minimum=_ZERO_,no_precipitation_dilution=.false.,no_river_dilution=.true.)

   ! Alkalinity may be apoproximated from salinity and temperature, or feature as separate state variable.
   if (.not. alk_param) then
     call register_state_variable(modelinfo,self%id_alk,'alk','mEq/m**3','alkalinity', &
                                    alk_initial,minimum=_ZERO_,no_precipitation_dilution=.false.,no_river_dilution=.true.)
   else
     call register_diagnostic_variable(modelinfo,self%id_alk_diag,'alk', 'mEq/m**3','alkalinity', &
                                                     time_treatment=time_treatment_averaged)
   end if

   ! Register diagnostic variables.
   call register_diagnostic_variable(modelinfo,self%id_ph, 'pH',      '-',          'pH',                              &
                         time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_pco2,'pCO2',    'ppm',        'CO2 partial pressure',           &
                         time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_CarbA,'CarbA',   'mmol/m**3',  'carbonic acid concentration',   &
                         time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_Bicarb,'Bicarb',  'mmol/m**3',  'bicarbonate ion concentration',&
                         time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_Carb,'Carb',    'mmol/m**3',  'carbonate ion concentration',    &
                         time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_Om_cal,'Om_cal',  '-',          'calcite saturation state',     &
                         time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_Om_arg,'Om_arg',  '-',          'aragonite saturation state',   &
                         time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_co2_flux,'CO2_flux','mmol/m**2/s','surface CO2 flux',   &
                         time_treatment=time_treatment_averaged)

   ! Register external dependencies.
   call register_dependency(modelinfo,self%id_temp,varname_temp)
   call register_dependency(modelinfo,self%id_salt,varname_salt)
   call register_dependency(modelinfo,self%id_pres,varname_pres)
   call register_dependency(modelinfo,self%id_dens,varname_dens)
   call register_dependency(modelinfo,self%id_wind,varname_wind_sf)
   if (self%pCO2a.eq._ZERO_) call register_dependency(modelinfo,self%id_pco2_surf,'pco2_surf')

   return

99 call fatal_error('pml_carbonate_init','Error reading namelist pml_carbonate')
100 call fatal_error('pml_carbonate_init','Namelist pml_carbonate was not found')

   end subroutine pml_carbonate_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of carbonate system model
!
! !INTERFACE:
   subroutine pml_carbonate_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_pml_carbonate),    intent(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   ! Environment
   real(rk) :: temp, salt, pres, dens
   real(rk) :: dic, TA
   real(rk) :: PCO2WATER, pH, HENRY, ca, bc, cb,Om_cal,Om_arg
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

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
      ! Divide by density/1000 to get alkalinity in uEq/kg.
      _GET_(self%id_alk,TA)
      TA = TA/dens*1.0D3
   end if

   ! Calculate carbonate system equilibrium.
   call CO2DYN(dic/1.0D3/dens, TA/1.0D6, temp, salt, PCO2WATER, pH, HENRY, ca, bc, cb)

   ! Call carbonate saturation state subroutine to calculate calcite and aragonite calcification states.
   call CaCO3_Saturation (temp, salt, pres, cb, Om_cal, Om_arg)

   ! Store diagnostic variables.
   _SET_DIAG_(self%id_ph    ,ph)
   _SET_DIAG_(self%id_pco2  ,PCO2WATER*1.0D6)        ! to ppm
   _SET_DIAG_(self%id_CarbA ,ca       *1.0D3*dens)   ! from mol/kg to mmol/m**3
   _SET_DIAG_(self%id_Bicarb,bc       *1.0D3*dens)   ! from mol/kg to mmol/m**3
   _SET_DIAG_(self%id_Carb  ,cb       *1.0D3*dens)   ! from mol/kg to mmol/m**3
   _SET_DIAG_(self%id_Om_cal,Om_cal)
   _SET_DIAG_(self%id_Om_arg,Om_arg)
   if (self%alk_param) then
      _SET_DIAG_(self%id_alk_diag,TA*dens*1.0D-3)    ! from uEg/kg to mmol/m**3
   end if

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine pml_carbonate_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Air-sea exchange for the carbonate system model
!
! !INTERFACE:
   subroutine pml_carbonate_get_surface_exchange(self,_FABM_ARGS_GET_SURFACE_EXCHANGE_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_pml_carbonate), intent(in)    :: self
   _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   ! Environment
   real(rk) :: temp, salt, wnd, dens

   ! State
   real(rk) :: dic, TA

   ! Temporary variables
   real(rk) :: PCO2WATER, pH, HENRY, ca, bc, cb, fl, pCO2a

   ! Parameters
   real(rk), parameter :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)
   _GET_(self%id_salt,salt)
   _GET_(self%id_dens,dens)
   _GET_HORIZONTAL_(self%id_wind,wnd)
   if (self%pCO2a.eq._ZERO_) then
      _GET_HORIZONTAL_(self%id_pco2_surf,pCO2a)
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
      TA = TA/dens*1.0d3
   end if

   ! Calculate carbonate system equilibrium.
   call CO2DYN(dic/1.0D3/dens, TA/1.0D6, temp, salt, PCO2WATER, pH, HENRY, ca, bc, cb)

   ! Calculate air-sea exchange of CO2 (positive flux is from atmosphere to water)
   call Air_sea_exchange(temp, wnd, PCO2WATER*1.0D6, pCO2a, Henry, dens/1.0D3, fl)

   ! Transfer surface exchange value to FABM.
   _SET_SURFACE_EXCHANGE_(self%id_dic,fl/secs_pr_day)

   ! Also store surface flux as diagnostic variable.
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_co2_flux,fl/secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_

   end subroutine pml_carbonate_get_surface_exchange
!EOC

!-----------------------------------------------------------------------

end module fabm_pml_carbonate

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

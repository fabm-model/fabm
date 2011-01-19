!$Id: co2sys.F90 119 2010-12-27 14:23:18Z jornbr $
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_co2sys --- shell around CO2/carbonate system model by
! Jerry Blackford (Plymouth Marine Laboratory), adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
module fabm_co2sys
!
! !DESCRIPTION:
! CO2 system model based on PML code.
!
! !USES:
   use fabm_types
   use fabm_driver

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_co2sys, co2sys_init, co2sys_do, co2sys_get_surface_exchange
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type type_co2sys
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_dic, id_alk
      _TYPE_DEPENDENCY_ID_          :: id_temp, id_salt, id_pres, id_wind, id_dens
      _TYPE_DIAGNOSTIC_VARIABLE_ID_ :: id_ph, id_pco2, id_CarbA, id_Bicarb, &
                                     & id_Carb, id_Om_cal, id_Om_arg, id_co2_flux, id_alk_diag
                                     
!     Model parameters
      REALTYPE :: TA_offset, TA_slope, pCO2a
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
   subroutine co2sys_init(self,modelinfo,namlst)
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
   type (type_co2sys), intent(out)     :: self
   type (type_model_info),intent(inout) :: modelinfo
   integer,               intent(in )   :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE :: dic_initial, alk_initial
   REALTYPE  :: alk_offset = 520.1, alk_slope = 51.24, pCO2a = 390.
   logical :: alk_param = .true.
   namelist /co2sys/ dic_initial, alk_initial, alk_param, alk_offset, alk_slope, pCO2a
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=co2sys,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day, and are converted here to values per second.
   self%TA_offset = alk_offset
   self%TA_slope  = alk_slope
   self%pCO2a     = pCO2a
   self%alk_param  = alk_param
      
   ! First state variable: total dissolved inorganic carbon
   self%id_dic = register_state_variable(modelinfo,'dic','mol C/m**3','total dissolved inorganic carbon', &
                                    dic_initial,minimum=_ZERO_,no_precipitation_dilution=.true.,no_river_dilution=.true.)
                                    
   ! Alkalinity may be apoproximated from salinity and temperature, or feature as separate state variable.
   if (.not. alk_param) then
     self%id_alk = register_state_variable(modelinfo,'alk','mEq/m**3','alkalinity', &
                                    alk_initial,minimum=_ZERO_,no_precipitation_dilution=.true.,no_river_dilution=.true.)
   else
     self%id_alk_diag = register_diagnostic_variable(modelinfo, 'alk', 'mEq/m**3','alkalinity', &
                                                     time_treatment=time_treatment_averaged)
   end if
                                    
   ! Register diagnostic variables.
   self%id_ph       = register_diagnostic_variable(modelinfo, 'pH',      '-',          'pH',                           &
                         time_treatment=time_treatment_averaged)
   self%id_pco2     = register_diagnostic_variable(modelinfo, 'pCO2',    'ppm',        'CO2 partial pressure',         &
                         time_treatment=time_treatment_averaged)
   self%id_CarbA    = register_diagnostic_variable(modelinfo, 'CarbA',   'mmol/m**3',  'carbonic acid concentration',  &
                         time_treatment=time_treatment_averaged)
   self%id_Bicarb   = register_diagnostic_variable(modelinfo, 'Bicarb',  'mmol/m**3',  'bicarbonate ion concentration',&
                         time_treatment=time_treatment_averaged)
   self%id_Carb     = register_diagnostic_variable(modelinfo, 'Carb',    'mmol/m**3',  'carbonate ion concentration',  &
                         time_treatment=time_treatment_averaged)
   self%id_Om_cal   = register_diagnostic_variable(modelinfo, 'Om_cal',  '-',          'calcite saturation state',     &
                         time_treatment=time_treatment_averaged)
   self%id_Om_arg   = register_diagnostic_variable(modelinfo, 'Om_arg',  '-',          'aragonite saturation state',   &
                         time_treatment=time_treatment_averaged)
   self%id_co2_flux = register_diagnostic_variable(modelinfo, 'CO2_flux','mmol/m**2/s','surface CO2 flux',             &
                         time_treatment=time_treatment_averaged, shape=shape_hz)
   
   ! Register external dependencies.
   self%id_temp = register_dependency(modelinfo, varname_temp)
   self%id_salt = register_dependency(modelinfo, varname_salt)
   self%id_pres = register_dependency(modelinfo, varname_pres)
   self%id_dens = register_dependency(modelinfo, varname_dens)
   self%id_wind = register_dependency(modelinfo, varname_wind_sf,shape=shape_hz)

   return

99 call fatal_error('co2sys_init','Error reading namelist co2sys')
   
   end subroutine co2sys_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of carbonate system model
!
! !INTERFACE:
   subroutine co2sys_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_co2sys),    intent(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   ! Environment
   REALTYPE :: temp, salt, pres, dens
   REALTYPE :: dic, TA
   REALTYPE :: PCO2WATER, pH, HENRY, ca, bc, cb,Om_cal,Om_arg
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Get environmental variables.
   _GET_DEPENDENCY_(self%id_temp,temp)
   _GET_DEPENDENCY_(self%id_salt,salt)
   _GET_DEPENDENCY_(self%id_pres,pres)
   _GET_DEPENDENCY_(self%id_dens,dens)

   ! Get current value for total dissolved inorganic carbon (our own state variable).
   _GET_STATE_(self%id_dic,dic)

   if (self%alk_param) then
      ! Linearly approximate alkalinity (uEq/kg) from salinity.
      TA = self%TA_offset + self%TA_slope*salt
   else
      ! Alkalinity (mEq/m**3) is a separate state variable.
      ! Divide by density/1000 to get alkalinity in uEq/kg.
      _GET_STATE_(self%id_alk,TA)
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

   end subroutine co2sys_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Air-sea exchange for the carbonate system model
!
! !INTERFACE:
   subroutine co2sys_get_surface_exchange(self,_FABM_ARGS_GET_SURFACE_EXCHANGE_)
!
! !DESCRIPTION:
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_co2sys), intent(in)    :: self
   _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   ! Environment
   REALTYPE :: temp, salt, wnd, dens
   
   ! State
   REALTYPE :: dic, TA
   
   ! Temporary variables
   REALTYPE :: PCO2WATER, pH, HENRY, ca, bc, cb, fl

   ! Parameters
   REALTYPE, parameter :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_HZ_LOOP_BEGIN_

   _GET_DEPENDENCY_(self%id_temp,temp)
   _GET_DEPENDENCY_(self%id_salt,salt)
   _GET_DEPENDENCY_(self%id_dens,dens)
   _GET_DEPENDENCY_HZ_(self%id_wind,wnd)

   _GET_STATE_(self%id_dic,dic)

   if (self%alk_param) then
      ! Linearly approximate alkalinity (uEq/kg) from salinity.
      TA = self%TA_offset + self%TA_slope*salt
   else
      ! Alkalinity (mEq/m**3) is a separate state variable
      ! Divide by density/1000 to get alkalinity in uEq/kg.
      _GET_STATE_(self%id_alk,TA)
      TA = TA/dens*1.0d3
   end if

   ! Calculate carbonate system equilibrium.
   call CO2DYN(dic/1.0D3/dens, TA/1.0D6, temp, salt, PCO2WATER, pH, HENRY, ca, bc, cb)

   ! Calculate air-sea exchange of CO2 (positive flux is from atmosphere to water)
   call Air_sea_exchange(temp, wnd, PCO2WATER*1.0D6, self%pCO2a, Henry, dens/1.0D3, fl)
   
   ! Transfer surface exchange value to FABM.
   _SET_SURFACE_EXCHANGE_(self%id_dic,fl/secs_pr_day)

   ! Also store surface flux as diagnostic variable.
   _SET_DIAG_HZ_(self%id_co2_flux,fl/secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_HZ_LOOP_END_

   end subroutine co2sys_get_surface_exchange
!EOC

!-----------------------------------------------------------------------

end module fabm_co2sys

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

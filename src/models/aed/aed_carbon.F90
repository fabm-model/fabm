!###############################################################################
!#                                                                             #
!# aed_carbon.F90                                                              #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#ifdef _FABM_F2003_

#include "fabm_driver.h"
#include "aed.h"

!
MODULE aed_carbon
!-------------------------------------------------------------------------------
! aed_carbon --- carbon biogeochemical model
!
! The AED module carbon contains equations that describe exchange of
! soluable reactive carbon across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE fabm_types
   USE fabm_driver

   USE aed_util,  ONLY: aed_gas_piston_velocity

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_carbon, aed_carbon_create
!
   TYPE,extends(type_base_model) :: type_aed_carbon
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_dic, id_pH, id_ch4, id_oxy, id_Fsed_dic
      _TYPE_DEPENDENCY_ID_          :: id_temp, id_salt, id_wind
      _TYPE_DIAGNOSTIC_VARIABLE_ID_ :: id_sed_dic, id_ch4ox, id_atm_co2_exch
      _TYPE_CONSERVED_QUANTITY_ID_  :: id_totC

!     Model parameters
      REALTYPE :: Fsed_dic,Ksed_dic,theta_sed_dic
      REALTYPE :: Fsed_ch4,Ksed_ch4,theta_sed_ch4
      REALTYPE :: Rch4ox,Kch4ox,vTch4ox,atmco2,ionic
      LOGICAL  :: use_oxy,use_dic,use_sed_model
      LOGICAL  :: simDIC, simCH4

      CONTAINS  ! Model Parameters
!       procedure :: initialize               => aed_carbon_init
        procedure :: do                       => aed_carbon_do
        procedure :: do_ppdd                  => aed_carbon_do_ppdd
        procedure :: do_benthos               => aed_carbon_do_benthos
        procedure :: get_conserved_quantities => aed_carbon_get_conserved_quantities
        procedure :: get_surface_exchange     => aed_carbon_get_surface_exchange
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
FUNCTION aed_carbon_create(namlst,name,parent) RESULT(self)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: namlst
   CHARACTER(len=*),INTENT(in)              :: name
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_aed_carbon),POINTER :: self

   REALTYPE          :: pH_initial=7.5
   REALTYPE          :: ionic = 0.0
   REALTYPE          :: dic_initial=4.5
   REALTYPE          :: Fsed_dic = 3.5
   REALTYPE          :: Ksed_dic = 30.0
   REALTYPE          :: theta_sed_dic = 1.0
   CHARACTER(len=64) :: Fsed_dic_variable=''
   REALTYPE          :: ch4_initial=4.5
   REALTYPE          :: Fsed_ch4 = 3.5
   REALTYPE          :: Ksed_ch4 = 30.0
   REALTYPE          :: theta_sed_ch4 = 1.0
   CHARACTER(len=64) :: Fsed_ch4_variable=''
   REALTYPE          :: Rch4ox = 0.01
   REALTYPE          :: Kch4ox = 0.01
   REALTYPE          :: vTch4ox= 1.05
   REALTYPE          :: atmco2 = 367e-6
   CHARACTER(len=64) :: methane_reactant_variable=''


   REALTYPE,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_carbon/ dic_initial,pH_initial,ionic,Fsed_dic,Ksed_dic,theta_sed_dic,Fsed_dic_variable, &
                         ch4_initial,Fsed_ch4,Ksed_ch4,theta_sed_ch4,Fsed_ch4_variable, &
                         atmco2,Rch4ox,Kch4ox,vTch4ox,methane_reactant_variable

!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL self%initialize(name,parent)

   ! Read the namelist
   read(namlst,nml=aed_carbon,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%Fsed_dic      = Fsed_dic/secs_pr_day
   self%Ksed_dic      = Ksed_dic
   self%theta_sed_dic = theta_sed_dic
   self%ionic         = ionic
   self%Fsed_ch4      = Fsed_ch4/secs_pr_day
   self%Ksed_ch4      = Ksed_ch4
   self%theta_sed_ch4 = theta_sed_ch4
   self%Rch4ox        = Rch4ox/secs_pr_day
   self%Kch4ox        = Kch4ox
   self%vTch4ox       = vTch4ox
   self%atmco2        = atmco2
   self%simDIC        = .false.
   self%simCH4        = .false.


   ! Register state variables
   IF(dic_initial>MISVAL) THEN
     self%id_dic = self%register_state_variable('dic','mmol/m**3','dissolved inorganic carbon',     &
                                      dic_initial,minimum=_ZERO_,no_river_dilution=.false.)
     self%simDIC = .true.
     self%id_pH = self%register_state_variable('pH','-','pH',     &
                                      pH_initial,minimum=_ZERO_,no_river_dilution=.true.)
   END IF

   IF(ch4_initial>MISVAL) THEN
     self%id_ch4 = self%register_state_variable('ch4','mmol/m**3','methane',    &
                                    ch4_initial,minimum=_ZERO_,no_river_dilution=.false.)
     self%simCH4 = .true.
   END IF

   ! Register external state variable dependencies
   self%use_oxy = methane_reactant_variable .NE. '' !This means oxygen module switched on
   IF (self%use_oxy) THEN
     self%id_oxy = self%register_state_dependency(methane_reactant_variable)
   ENDIF

   self%use_sed_model = Fsed_dic_variable .NE. ''
   IF (self%use_sed_model) &
       self%id_Fsed_dic = self%register_state_dependency(Fsed_dic_variable,benthic=.true.)

   ! Register diagnostic variables
   self%id_ch4ox  = self%register_diagnostic_variable('ch4ox','/d',                    &
                                                           'methane oxidation rate',   &
                                        time_treatment=time_treatment_step_integrated)
   self%id_sed_dic = self%register_diagnostic_variable('sed_dic','mmol/m**2/d',        &
                                                      'Filterable reactive carbon',    &
                          time_treatment=time_treatment_step_integrated, shape=shape_hz)

   self%id_atm_co2_exch = self%register_diagnostic_variable('atm_co2_exch',            &
                             'mmol/m**2/d', 'CO2 exchange across atm/water interface', &
                          time_treatment=time_treatment_step_integrated, shape=shape_hz)

   ! Register conserved quantities
   self%id_totC = self%register_conserved_quantity('TP','mmol/m**3','Total carbon')

   ! Register environmental dependencies
   self%id_temp = self%register_dependency(varname_temp)
   self%id_salt = self%register_dependency(varname_salt)
   self%id_wind = self%register_dependency(varname_wind_sf,  shape=shape_hz)

   RETURN

99 CALL fatal_error('aed_carbon_init','Error reading namelist aed_carbon')

END FUNCTION aed_carbon_create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_carbon_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_carbon),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   REALTYPE           :: dic,ch4,oxy,temp
   REALTYPE           :: ch4oxidation

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   IF(self%simDIC .AND. self%simCH4) THEN
      ! Retrieve current (local) state variable values.
      _GET_STATE_(self%id_dic,dic) ! carbon
      _GET_STATE_(self%id_ch4,ch4) ! carbon

      ! Retrieve current dependent state variable values.
      IF (self%use_oxy) THEN ! & use_oxy
         _GET_STATE_(self%id_oxy,oxy) ! oxygen
      ELSE
         oxy = 0.0
      ENDIF

      ! Retrieve current environmental conditions.
      _GET_DEPENDENCY_(self%id_temp,temp)  ! temperature

      ! Define some intermediate quantities units mmol C/m3/day
      ch4oxidation = aed_carbon_fch4ox(self,oxy,temp)

      ! Set temporal derivatives
      _SET_ODE_(self%id_dic,ch4*ch4oxidation)
      _SET_ODE_(self%id_ch4,-ch4*ch4oxidation)

      ! If an externally maintained oxygen pool is present, take nitrification from it
      IF (self%use_oxy) then ! & use_oxy
         _SET_ODE_(self%id_oxy,-(32./12.)*ch4*ch4oxidation)
      ENDIF

      ! Export diagnostic variables
      _SET_DIAG_(self%id_ch4ox, ch4oxidation)
   ENDIF

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_carbon_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_carbon_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of carbon biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_carbon),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   REALTYPE                   :: dic
   REALTYPE                   :: diff_dic
   REALTYPE, parameter        :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   IF(self%simDIC .AND. self%simCH4) THEN
      ! Retrieve current (local) state variable values.
      _GET_STATE_(self%id_dic,dic) ! carbon

      ! Set temporal derivatives
      diff_dic = 0.

      _SET_PP_(self%id_dic,self%id_dic,diff_dic)
   END IF

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_carbon_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_carbon_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED carbon.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_carbon),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   ! Environment
   REALTYPE :: temp

   ! State
   REALTYPE :: dic,oxy

   ! Temporary variables
   REALTYPE :: dic_flux, Fsed_dic

   ! Parameters
   REALTYPE,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.self%simDIC) RETURN

   ! Enter spatial loops (if any)
   _FABM_HZ_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_DEPENDENCY_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_dic,dic) ! carbon

   IF ( self%use_sed_model ) THEN
      _GET_STATE_BEN_(self%id_Fsed_dic,Fsed_dic)
   ELSE
       Fsed_dic = self%Fsed_dic
   ENDIF

   IF (self%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      _GET_STATE_(self%id_oxy,oxy)
      dic_flux = Fsed_dic * self%Ksed_dic/(self%Ksed_dic+oxy) * (self%theta_sed_dic**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
      dic_flux = Fsed_dic * (self%theta_sed_dic**(temp-20.0))
   ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to FABM.
   !_SET_BOTTOM_FLUX_(self%id_dic,dic_flux/secs_pr_day)
   !_SET_ODE_SED_FLUX_(self%id_dic,dic_flux)
   _SET_BOTTOM_EXCHANGE_(self%id_dic,dic_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_dic,-dic_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
   _SET_DIAG_HZ_(self%id_sed_dic,dic_flux)

   ! Leave spatial loops (if any)
   _FABM_HZ_LOOP_END_

END SUBROUTINE aed_carbon_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_carbon_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities (currently only carbon)
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_carbon),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   REALTYPE :: dic,ch4
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   dic = _ZERO_
   ch4 = _ZERO_

   ! Retrieve current (local) state variable values.
   IF(self%simDIC) THEN
     _GET_STATE_(self%id_dic,dic) ! DIC
   END IF
   IF(self%simCH4) THEN
     _GET_STATE_(self%id_ch4,ch4) ! CH4
   END IF

   ! Total is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totC,dic+ch4)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_carbon_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
PURE REALTYPE FUNCTION aed_carbon_fch4ox(self,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for methane oxidation
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_carbon),INTENT(in) :: self
   REALTYPE,INTENT(in)                  :: oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (self%use_oxy) THEN
      aed_carbon_fch4ox = self%Rch4ox * oxy/(self%Kch4ox+oxy) * (self%vTch4ox**(temp-20.0))
   ELSE
      aed_carbon_fch4ox = self%Rch4ox * (self%vTch4ox**(temp-20.0))
   ENDIF

END FUNCTION aed_carbon_fch4ox
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE aed_carbon_get_surface_exchange(self,_FABM_ARGS_GET_SURFACE_EXCHANGE_)
!-------------------------------------------------------------------------------
! Air-sea exchange for the aed carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_carbon),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_
!
!LOCALS
   ! Environment
   REALTYPE :: temp, salt, wind

   ! State
   REALTYPE :: dic,ph

   ! Temporary variables
   REALTYPE :: pCO2,FCO2
   REALTYPE :: Ko, KCO2
   REALTYPE :: Tabs,windHt

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.self%simDIC) RETURN

   ! Enter spatial loops (if any)
   _FABM_HZ_LOOP_BEGIN_

   !Get dependent state variables from physical driver
   _GET_DEPENDENCY_(self%id_temp,temp)     ! Temperature (degrees Celsius)
   _GET_DEPENDENCY_(self%id_salt,salt)     ! Salinity (psu)
   _GET_DEPENDENCY_HZ_(self%id_wind,wind)  ! Wind speed at 10 m above surface (m/s)
   windHt = 10.

    ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_dic,dic) ! Concentration of carbon in surface layer
   _GET_STATE_(self%id_pH,ph) ! Concentration of carbon in surface layer

   kCO2 = aed_gas_piston_velocity(windHt,wind,temp,salt)

   ! Solubility, Ko (mol/L/atm)
   Tabs = temp + 273.15
   Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
          + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
   Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
   Ko = exp(Ko)

   ! pCO2 in surface water layer
   pCO2 = aed_carbon_co2(self,temp,dic,ph) / Ko

   ! FCO2 = kCO2 * Ko * (pCO2 - PCO2a)
   ! pCO2a = 367e-6 (Keeling & Wharf, 1999)

   ! mmol/m2/s = m/s * mmol/L/atm * atm
   FCO2 = kCO2 * Ko*1e6 * (pCO2 - self%atmco2)

   ! Transfer surface exchange value to FABM (mmmol/m2) converted by driver.
   _SET_SURFACE_EXCHANGE_(self%id_dic,FCO2)

   ! Also store oxygen flux across the atm/water interface as diagnostic variable (mmmol/m2).
   _SET_DIAG_HZ_(self%id_atm_co2_exch,FCO2)

   ! Leave spatial loops (if any)
   _FABM_HZ_LOOP_END_

END SUBROUTINE aed_carbon_get_surface_exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
PURE REALTYPE FUNCTION aed_carbon_co2(self,temp,dic,pH)
!-------------------------------------------------------------------------------
! CO2 concentration of DIC at fixed T
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_carbon),INTENT(in) :: self
   REALTYPE, INTENT(IN)                 :: dic, temp, pH
!
!LOCALS
   ! Temporary variables
   REALTYPE :: K_h, Kw, Ka1, Ka2, i_f
   REALTYPE :: H, CO2, HCO3, CO3, TA
!-------------------------------------------------------------------------------
!BEGIN

   ! Acidity constants temperature dependence

   ! pKh =  -0.000075324675x2 + 0.016279653680x + 1.110424242424
   ! pKa1 = 0.000142121212x2 - 0.012648181818x + 6.577539393939
   ! pKa2 =  0.000113679654x2 - 0.014687186147x + 10.625769696970
   ! pKw =   0.000201991342x2 - 0.043419653680x + 14.949709090909

   K_h = -0.000075324675*temp*temp + 0.016279653680*temp + 1.110424242424
   Ka1 = 0.000142121212*temp*temp - 0.012648181818*temp + 6.577539393939
   Ka2 = 0.000113679654*temp*temp - 0.014687186147*temp + 10.625769696970
   Kw = 0.000201991342*temp*temp - 0.043419653680*temp + 14.949709090909


   ! Ionic strength dependence

   ! 1st calculate function f
   i_f = (((SQRT(self%ionic)) / (1+SQRT(self%ionic))) -0.20*self%ionic) * &
                       (298.0/(temp+273.))**0.666667

   ! pKh = pKh(0) + bI
   ! b = 0.105 (Butler, 1982)
   K_h = K_h + 0.105*self%ionic

   ! pKw = pKw(0) - f
   Kw = Kw - i_f

   ! pKa1 = pKa1(0) - f - bI
   Ka1 = Ka1 - i_f - 0.105*self%ionic

   !pKa2 = pKa2(0) - 2f
   Ka2 = Ka2 + 2.0*i_f

   ! Convert from pK etc to Kh, Kw, Ka1, Ka2
   K_h  = 10.**(-K_h)
   Ka1 = 10.**(-Ka1)
   Ka2 = 10.**(-Ka2)
   Kw  = 10.**(-Kw)


   ! Calculate the speciation to know the molar mass of DIC                                                             !
   H    = 10.**(-pH)
   CO3  = (Ka1*Ka2)/(H*H + Ka1*H + Ka1*Ka2)
   HCO3 = (Ka1*H)/(H*H + Ka1*H + Ka1*Ka2)
   CO2  = (H*H)/(H*H + Ka1*H + Ka1*Ka2)


   ! and update speciation (mol C/L)
   CO3  = dic*CO3
   HCO3 = dic*HCO3
   CO2  = dic*CO2

   ! calculate TA for the previous timestep
   TA = dic * (Ka1*H + 2.0*Ka1*Ka2) / (H*H + Ka1*H + Ka1*Ka2)
   TA = TA + (Kw/H) - H

   aed_carbon_co2 = CO2

END FUNCTION aed_carbon_co2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





END MODULE aed_carbon
#endif

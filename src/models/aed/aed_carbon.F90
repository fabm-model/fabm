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

#include "aed.h"

!
MODULE aed_carbon
!-------------------------------------------------------------------------------
! aed_carbon --- carbon biogeochemical model
!
! The AED module carbon contains equations that describe exchange of
! soluable reactive carbon across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,  ONLY: aed_gas_piston_velocity

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_type_carbon
!
   TYPE,extends(type_base_model) :: aed_type_carbon
!     Variable identifiers
      TYPE (type_state_variable_id)      :: id_dic, id_pH, id_ch4, id_oxy
      TYPE (type_horizontal_dependency_id) :: id_Fsed_dic
      TYPE (type_dependency_id)          :: id_temp, id_salt
      TYPE (type_horizontal_dependency_id)  :: id_wind
      TYPE (type_diagnostic_variable_id) :: id_ch4ox
      TYPE (type_horizontal_diagnostic_variable_id) :: id_sed_dic
      TYPE (type_horizontal_diagnostic_variable_id) :: id_atm_co2_exch,id_atm_ch4_exch

!     Model parameters
      AED_REAL :: Fsed_dic,Ksed_dic,theta_sed_dic
      AED_REAL :: Fsed_ch4,Ksed_ch4,theta_sed_ch4
      AED_REAL :: Rch4ox,Kch4ox,vTch4ox,atmco2,ionic
      LOGICAL  :: use_oxy,use_dic,use_sed_model
      LOGICAL  :: simDIC, simCH4

      CONTAINS  ! Model Parameters
        PROCEDURE :: initialize               => aed_init_carbon
        PROCEDURE :: do                       => aed_carbon_do
        PROCEDURE :: do_ppdd                  => aed_carbon_do_ppdd
        PROCEDURE :: do_benthos               => aed_carbon_do_benthos
        PROCEDURE :: get_surface_exchange     => aed_carbon_get_surface_exchange
!       PROCEDURE :: check_state              => aed_carbon_check_state
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_carbon(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_carbon),TARGET,INTENT(inout) :: self
   INTEGER,INTENT(in) :: namlst
!
!LOCALS

   INTEGER  :: status

   AED_REAL          :: pH_initial=7.5
   AED_REAL          :: ionic = 0.0
   AED_REAL          :: dic_initial=4.5
   AED_REAL          :: Fsed_dic = 3.5
   AED_REAL          :: Ksed_dic = 30.0
   AED_REAL          :: theta_sed_dic = 1.0
   CHARACTER(len=64) :: Fsed_dic_variable=''
   AED_REAL          :: ch4_initial=4.5
   AED_REAL          :: Fsed_ch4 = 3.5
   AED_REAL          :: Ksed_ch4 = 30.0
   AED_REAL          :: theta_sed_ch4 = 1.0
   CHARACTER(len=64) :: Fsed_ch4_variable=''
   AED_REAL          :: Rch4ox = 0.01
   AED_REAL          :: Kch4ox = 0.01
   AED_REAL          :: vTch4ox= 1.05
!  AED_REAL          :: atmco2 = 367e-6
   AED_REAL          :: atmco2 = 367.
   CHARACTER(len=64) :: methane_reactant_variable=''


   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_carbon/ dic_initial,pH_initial,ionic,Fsed_dic,Ksed_dic,theta_sed_dic,Fsed_dic_variable, &
                         ch4_initial,Fsed_ch4,Ksed_ch4,theta_sed_ch4,Fsed_ch4_variable, &
                         atmco2,Rch4ox,Kch4ox,vTch4ox,methane_reactant_variable

!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_carbon,iostat=status)
   IF (status /= 0) THEN
      print *,'Error reading namelist aed_carbon'
      STOP
   ENDIF

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
   IF (dic_initial>MISVAL) THEN
     CALL self%register_state_variable(self%id_dic,'dic','mmol/m**3','dissolved inorganic carbon',     &
                                      dic_initial,minimum=zero_,no_river_dilution=.false.)
     self%simDIC = .true.
     CALL self%register_state_variable(self%id_pH,'pH','-','pH',     &
                                      pH_initial,minimum=zero_,no_river_dilution=.true.)
   ENDIF

   IF (ch4_initial>MISVAL) THEN
     CALL self%register_state_variable(self%id_ch4,'ch4','mmol/m**3','methane',    &
                                    ch4_initial,minimum=zero_,no_river_dilution=.false.)
     self%simCH4 = .true.
   ENDIF

   !# Register external state variable dependencies
   self%use_oxy = methane_reactant_variable .NE. '' !This means oxygen module switched on
   IF (self%use_oxy) THEN
     CALL self%register_state_dependency(self%id_oxy,methane_reactant_variable)
   ENDIF

   self%use_sed_model = Fsed_dic_variable .NE. ''
   IF (self%use_sed_model) THEN
      CALL self%register_horizontal_dependency(self%id_Fsed_dic,Fsed_dic_variable)
      CALL self%request_coupling(self%id_Fsed_dic,Fsed_dic_variable)
   ENDIF

   !# Register diagnostic variables
   CALL self%register_diagnostic_variable(self%id_ch4ox,'ch4ox','/d', 'methane oxidation rate')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_dic,'sed_dic','mmol/m**2/d',        &
                                                      'Filterable reactive carbon')

   CALL self%register_horizontal_diagnostic_variable(self%id_atm_co2_exch,'atm_co2_exch',            &
                             'mmol/m**2/d', 'CO2 exchange across atm/water interface')
   CALL self%register_horizontal_diagnostic_variable(self%id_atm_ch4_exch,'atm_ch4_exch',            &
                             'mmol/m**2/d', 'CH4 exchange across atm/water interface')

   !# Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature)
   CALL self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   CALL self%register_dependency(self%id_wind,standard_variables%wind_speed)
END SUBROUTINE aed_init_carbon
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_carbon_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_carbon),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   AED_REAL :: dic,ch4,oxy,temp
   AED_REAL :: ch4oxidation
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   IF(self%simDIC .AND. self%simCH4) THEN
      ! Retrieve current (local) state variable values.
      _GET_(self%id_dic,dic) ! carbon
      _GET_(self%id_ch4,ch4) ! carbon

      ! Retrieve current dependent state variable values.
      IF (self%use_oxy) THEN ! & use_oxy
         _GET_(self%id_oxy,oxy) ! oxygen
      ELSE
         oxy = 0.0
      ENDIF

      ! Retrieve current environmental conditions.
      _GET_(self%id_temp,temp)  ! temperature

      ! Define some intermediate quantities units mmol C/m3/day
      ch4oxidation = aed_carbon_fch4ox(self%use_oxy,self%Rch4ox,self%Kch4ox,self%vTch4ox,oxy,temp)

      ! Set temporal derivatives
      _SET_ODE_(self%id_dic,ch4*ch4oxidation)
      _SET_ODE_(self%id_ch4,-ch4*ch4oxidation)

      ! If an externally maintained oxygen pool is present, take nitrification from it
      IF (self%use_oxy) then ! & use_oxy
         _SET_ODE_(self%id_oxy,-(32./12.)*ch4*ch4oxidation)
      ENDIF

      ! Export diagnostic variables
      _SET_DIAGNOSTIC_(self%id_ch4ox, ch4oxidation)
   ENDIF

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_carbon_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_carbon_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of carbon biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_carbon),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL                   :: dic
   AED_REAL                   :: diff_dic
   AED_REAL, parameter        :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   IF(self%simDIC .AND. self%simCH4) THEN
      ! Retrieve current (local) state variable values.
      _GET_(self%id_dic,dic) ! carbon

      ! Set temporal derivatives
      diff_dic = 0.

      _SET_PP_(self%id_dic,self%id_dic,diff_dic)
   END IF

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_carbon_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_carbon_get_surface_exchange(self,_ARGUMENTS_DO_SURFACE_)
!-------------------------------------------------------------------------------
! Air-sea exchange for the aed carbon model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_carbon),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind

   ! State
   AED_REAL :: dic,ph,ch4

   ! Temporary variables
   AED_REAL :: pCO2,FCO2,FCH4
   AED_REAL :: Ko,kCH4,KCO2, CH4solub
   AED_REAL :: Tabs,windHt,atm
   AED_REAL :: A1,A2,A3,A4,B1,B2,B3,logC

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.self%simDIC .AND. .NOT.self%simCH4) RETURN

   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   !Get dependent state variables from physical driver
   _GET_(self%id_temp,temp)     ! Temperature (degrees Celsius)
   _GET_(self%id_salt,salt)     ! Salinity (psu)
   _GET_HORIZONTAL_(self%id_wind,wind)  ! Wind speed at 10 m above surface (m/s)
   windHt = 10.


   Tabs = temp + 273.15

   ! CO2 flux
   IF(self%simDIC) THEN

      ! Retrieve current (local) state variable values.
     _GET_(self%id_dic,dic) ! Concentration of carbon in surface layer
     _GET_(self%id_pH,ph)   ! Concentration of carbon in surface layer

     kCO2 = aed_gas_piston_velocity(windHt,wind,temp,salt,schmidt_model=2)

     ! Solubility, Ko (mol/L/atm)
     Ko = -58.0931+90.5069*(100.0/Tabs) + 22.294*log(Tabs/100.0) &
          + 0.027766*salt - 0.025888*salt*(Tabs/100.0)
     Ko = Ko + 0.0050578*salt*(Tabs/100.0)*(Tabs/100.0)
     Ko = exp(Ko)

     ! pCO2 in surface water layer
     pCO2 = aed_carbon_co2(self%ionic,temp,dic,ph) / Ko

     ! FCO2 = kCO2 * Ko * (pCO2 - PCO2a)
     ! pCO2a = 367e-6 (Keeling & Wharf, 1999)

     !------ Yanti correction (20/5/2013) ----------------------------------------
     ! pCO2 is actually in uatm (=ppm)
     ! mmol/m2/s = m/s * mmol/L/atm * atm
     FCO2 = kCO2 * Ko * (pCO2 - self%atmco2)

     ! FCO2 = - kCO2 * Ko*1e6 * ((pCO2 * 1e-6) - self%atmco2) ! dCO2/dt
     !----------------------------------------------------------------------------

     ! Transfer surface exchange value to FABM (mmmol/m2) converted by driver.
     _SET_SURFACE_EXCHANGE_(self%id_dic,-FCO2)

     ! Also store oxygen flux across the atm/water interface as diagnostic variable (mmmol/m2).
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_atm_co2_exch,FCO2)
   END IF

   ! CH4 flux
   IF(self%simCH4) THEN
     ! Algorithm from Arianto Santoso <abs11@students.waikato.ac.nz>

     ! Concentration of methane in surface layer
     _GET_(self%id_ch4,ch4)

     ! Piston velocity for CH4
     kCH4 = aed_gas_piston_velocity(windHt,wind,temp,salt,schmidt_model=4)

     ! Solubility, Ko (mol/L/atm)

     atm = 1.76 * 1e-6 !## current atmospheric CH4 data from NOAA (in ppm)

     A1 = -415.2807
     A2 = 596.8104
     A3 = 379.2599
     A4 = -62.0757
     B1 = -0.05916
     B2 = 0.032174
     B3 = -0.0048198

     logC = (log(atm)) + A1 + (A2 * (100./Tabs)) + (A3 * log (Tabs/100.)) + (A4 * (Tabs/100.)) + &
                 salt * (B1 + (B2  * (Tabs/100.)) + (B3 * (Tabs/100.)*(Tabs/100.)))

     CH4solub = exp(logC) * 1e-3


     ! mmol/m2/s = m/s * mmol/m3
     FCH4 = kCH4 *  (ch4 - CH4solub)

     !----------------------------------------------------------------------------

     ! Transfer surface exchange value to FABM (mmmol/m2) converted by driver.
     _SET_SURFACE_EXCHANGE_(self%id_ch4,-FCH4)

     ! Also store ch4 flux across the atm/water interface as diagnostic variable (mmmol/m2).
     _SET_HORIZONTAL_DIAGNOSTIC_(self%id_atm_ch4_exch,FCH4)
   END IF


   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_
END SUBROUTINE aed_carbon_get_surface_exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_carbon_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED carbon.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_carbon),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: dic,oxy

   ! Temporary variables
   AED_REAL :: dic_flux, Fsed_dic

   ! Parameters
   AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN

   IF(.NOT.self%simDIC) RETURN

   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_(self%id_dic,dic) ! carbon

   IF ( self%use_sed_model ) THEN
      _GET_HORIZONTAL_(self%id_Fsed_dic,Fsed_dic)
   ELSE
       Fsed_dic = self%Fsed_dic
   ENDIF

   IF (self%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      _GET_(self%id_oxy,oxy)
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
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_dic,dic_flux)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_
END SUBROUTINE aed_carbon_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_carbon_check_state(self,_ARGUMENTS_CHECK_STATE_)
!-------------------------------------------------------------------------------
! Update pH after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_carbon),INTENT(in) :: self
   _DECLARE_ARGUMENTS_CHECK_STATE_
!
!LOCALS
   ! State
   AED_REAL :: dic, pH

!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT.self%simDIC) RETURN

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

    ! Retrieve current (local) state variable values.
!  _GET_(self%id_dic,dic) ! Concentration of carbon in surface layer
!  _GET_(self%id_pH, pH) ! Concentration of carbon in surface layer

!print*,"new pH = ",pH
!  _SET_STATE_(self%id_pH, pH)

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_carbon_check_state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed_carbon_fch4ox(use_oxy,Rch4ox,Kch4ox,vTch4ox,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for methane oxidation
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rch4ox,Kch4ox,vTch4ox,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      aed_carbon_fch4ox = Rch4ox * oxy/(Kch4ox+oxy) * (vTch4ox**(temp-20.0))
   ELSE
      aed_carbon_fch4ox = Rch4ox * (vTch4ox**(temp-20.0))
   ENDIF

END FUNCTION aed_carbon_fch4ox
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION aed_carbon_co2(ionic, temp, dic, pH)
!-------------------------------------------------------------------------------
! CO2 concentration of DIC at fixed T
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL, INTENT(IN) :: ionic, dic, temp, pH
!
!LOCALS
   ! Temporary variables
   AED_REAL :: K_h, Kw, Ka1, Ka2, i_f
   AED_REAL :: H, CO2, HCO3, CO3, TA
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
   i_f = (((SQRT(ionic)) / (1+SQRT(ionic))) -0.20*ionic) * &
                       (298.0/(temp+273.))**0.666667

   ! pKh = pKh(0) + bI
   ! b = 0.105 (Butler, 1982)
   K_h = K_h + 0.105*ionic

   ! pKw = pKw(0) - f
   Kw = Kw - i_f

   ! pKa1 = pKa1(0) - f - bI
   Ka1 = Ka1 - i_f - 0.105*ionic

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

!-----------------------------------------------------------------------
! OXYDEP is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Jorn Bruggemann
!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_niva_oxydep --- OXYDEP biogeochemical model based upon
! Yakushev et al, 2013  with  modifications by Jorn Bruggeman
! and adapted for FABM by Jorn Bruggeman
!
! !INTERFACE:
   module fabm_niva_oxydep
!
! !DESCRIPTION:
!
! OXYDEP targests on the silmplest possible way of parameterization of the oxygen  (DO) fate in changeable redox conditions.
! It has a simplified ecosystem, and simulates production of DO due to photosynthesis and consumation of DO for biota respiraion,
! OM mineralization, nitrification, and oxidation of reduced specied of S, Mn, Fe, present in suboxic conditions.
! OXYDEP consists of 6 state variables ( in N-units):
! - PHY - all the phototrophic organisms (phytoplankton and bacteria).
!   PHY grows due to photosynthesis, loses inorganic matter
!   due to respiraion, and loses organic matter in dissolved (DOM) and particulate (POM)
!   forms due to metabolism and mortality. PHY growth is limited by irradiance, temperature and NUT availability.
! - HET - heterotrophs, can consume PHY and POP,  produce DOM and DON and respirate NUT.
! - NUT - represents oxydized forms of nutrients (i.e. NO3 and NO2 for N),
!   that doesn't need additional  oxygen for nitrification.
! - DOM - is dissolved organic matter. DOM  includes all kinds of labile dissolved organic matter
!   and reduced forms of inorganic nutrients (i.e. NH4 and Urea for N).
! - POM - is particular organic matter (less labile than DOM). Temperature affects DOM and POM mineralization.
! - OXY - is dissolved oxygen.
! For the details of  OXYDEP  implemented here (actually, a previous version withour HET), see (Yakushev et al, 2013)
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev, Jorn Bruggeman
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_oxydep
!     Variable identifiers
      type (type_state_variable_id)        :: id_oxy,id_phy,id_het,id_nut,id_pom,id_dom
      type (type_dependency_id)            :: id_par,id_temp, id_salt
      type (type_horizontal_dependency_id) :: id_windspeed
      type (type_diagnostic_variable_id)   :: id_MortHet,id_RespHet,id_GrazPhy,id_GrazPOM,id_GrowthPhy,id_MortPhy,id_ExcrPhy
      type (type_diagnostic_variable_id)   :: id_DOM_decay_ox,id_DOM_decay_denitr,id_POM_decay_ox,id_POM_decay_denitr
      type (type_diagnostic_variable_id)   :: id_LimT,id_LimP,id_LimN,id_LimLight, id_N_fixation
!     Model parameters
      !----Phy -----------!
       real(rk) :: Max_uptake, bm, cm, Knut, r_phy_nut, r_phy_pom, r_phy_dom, r_phy_om_anox ,ir_min, Iopt, O2_add_mor_phy
      !----Het -----------!
       real(rk) :: r_phy_het, Kphy, r_pop_het, Kpop, r_het_nut, r_het_pom, Uz, Hz
      !----DOM, POM   --- !
       real(rk) :: r_pom_dom, r_pom_nut_oxy, r_pom_nut_nut, r_dom_nut_oxy, r_dom_nut_nut, Tda, beta_da
      ! real(rk) :: r_om_nut_sul  = 0.005    ! Specific rate of OM anoxic decay   	[1/d]
      !----Oxy -----------!
       real(rk) :: O2_suboxic
      !  Lower boundary
       real(rk) :: Bu,Trel,b_ox,b_dom_ox,b_dom_anox,b_nut
      ! Upper boundary for oxygen flux calculations
       real(rk) :: pvel,a0,a1,a2
      !  Stochiometric coefficients
       real(rk) :: NtoB, OtoN, NtoN
      !---sinking-------------------------------!
       real(rk) :: Wphy, Whet, Wpom
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the OXYDEP namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_niva_oxydep), intent(inout), target :: self
   integer,                  intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk),parameter ::  d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
  ! PHY
   call self%get_parameter(self%Max_uptake,'Max_uptake','1/d', 'Maximum nutrient uptake rate',                               default=5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Knut,      'Knut',      'nd',  'Half-sat.const. for uptake of NUT by PHY for NUT/PHY ratio', default=0.1_rk)
   call self%get_parameter(self%bm,        'bm',      '1/gradC',     'Coefficient for uptake rate dependence on t',          default=0.12_rk)
   call self%get_parameter(self%cm,        'cm',      'nd',          'Coefficient for uptake rate dependence on t',          default=1.4_rk)
   call self%get_parameter(self%ir_min,    'ir_min',      'nd',      'bioshading parameter ',                                default=25._rk)
   call self%get_parameter(self%Iopt,      'Iopt',    'Watts/m**2/h',  'Optimal irradiance',                                 default=25.0_rk)
   call self%get_parameter(self%r_phy_nut,     'r_phy_nut',     '1/d', 'Specific Phy  respiration rate',                          default=0.04_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_pom,     'r_phy_pom',     '1/d', 'Specific Phy rate of mortality',                         default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_dom,     'r_phy_dom',     '1/d', 'Specific Phy rate of excretion',                         default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_om_anox, 'r_phy_om_anox', '1/d', 'Specific additional Phy mortality in suboxic/anoxic conditions', default=0.4_rk,scale_factor=d_per_s)
   call self%get_parameter(self%O2_add_mor_phy, 'O2_add_mor_phy', 'mmol/m3', 'Threshold O2 value for additional Phy mortality', default=20._rk)
  ! HET
   call self%get_parameter(self%r_phy_het, 'r_phy_het', '1/d', 'Max.spec. rate of grazing of Het on PHY',                    default=2.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Kphy,      'Kphy',      'nd',  'Half-sat.const.for grazing of Het on Phy for Phy/Het ratio ',default=0.1_rk)
   call self%get_parameter(self%r_pop_het, 'r_pop_het', '1/d', 'Max.spec. rate of grazing of Het on POM ',                   default=0.7_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Kpop,      'Kpop',      'nd',  'Half-sat.const.for grazing of Het on POM for POM/Het ratio ',default=2._rk)
   call self%get_parameter(self%r_het_nut, 'r_het_nut', '1/d', 'Specific Het respiration rate ',                             default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_het_pom, 'r_het_pom', '1/d', 'Specific Het mortality rate ',                               default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Uz,        'Uz',        'nd',  'Food absorbency for Het',                                    default=0.5_rk)
   call self%get_parameter(self%Hz,        'Hz',        'nd',  'Ratio betw. diss. and part. excretes of Het ',               default=0.5_rk)
  ! POM
   call self%get_parameter(self%r_pom_dom,     'r_pom_dom',     '1/d', 'Specific rate of POM decomposition ',                default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_pom_nut_oxy, 'r_pom_nut_oxy', '1/d', 'Specific rate of POM oxic decay  ',                  default=0.03_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_pom_nut_nut, 'r_pom_nut_nut', '1/d', 'Specific rate of POM denitrification  ',             default=0.006_rk,scale_factor=d_per_s)
  ! DOM
   call self%get_parameter(self%r_dom_nut_oxy, 'r_dom_nut_oxy', '1/d', 'Specific rate of DOM oxic decay  ',                  default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_dom_nut_nut, 'r_dom_nut_nut', '1/d', 'Specific rate of DOM denitrification  ',             default=0.002_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Tda,           'Tda',           'nd',  'Coefficient for dependence of mineralization on t ', default=13._rk)
   call self%get_parameter(self%beta_da,       'beta_da',       'nd',  'Coefficient for dependence of mineralization on t ', default=20._rk)
  !---sinking-------------------------------!
   call self%get_parameter(self%Wphy,     'Wphy',     'm/s', 'vertical velocity of PHY (<0 for sinking)',           default=-0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Whet,     'Whet',     'm/s', 'vertical velocity of het (<0 for sinking)',           default=-0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Wpom,     'Wpom',     'm/s', 'vertical velocity of POM (<0 for sinking)',           default=-1.0_rk,scale_factor=d_per_s)
  !  Lower boundary
   call self%get_parameter(self%Bu,        'Bu',         'nd',      'Burial coeficient for lower boundary',                  default=0.25_rk)
   call self%get_parameter(self%Trel,      'Trel',       's/m',     'Relaxation time for exchange with teh sediments',       default=1e5_rk)
   call self%get_parameter(self%O2_suboxic,    'O2_suboxic',     'mmol/m3', 'Limiting O2 value for oxic/suboxic switch',             default=40._rk)
   call self%get_parameter(self%b_ox,      'b_ox',       'mmol/m3', 'OXY in the sediment',                                   default=0._rk)
   call self%get_parameter(self%b_dom_ox,  'b_dom_ox',   'mmol/m3', 'OM in the sediment (oxic conditions)',                  default=2._rk)
   call self%get_parameter(self%b_dom_anox,'b_dom_anox', 'mmol/m3', 'OM in the sediment (anoxic conditions) ',               default=6._rk)
   call self%get_parameter(self%b_nut,     'b_nut',      'mmol/m3', 'NUT in the sediment',                                   default=0._rk)
  ! Upper boundary for oxygen flux calculations
   call self%get_parameter(self%pvel,       'pvel',       'm/s',     'wind speed',                                            default=5._rk)
   call self%get_parameter(self%a0,         'a0',         'mmol/m3', 'oxygen saturation parameter  ',                         default=31.25_rk)
   call self%get_parameter(self%a1,         'a1',         'nd',      'oxygen saturation parameter  ',                         default=14.603_rk)
   call self%get_parameter(self%a2,         'a2',         '1/degC',  'oxygen saturation parameter  ',                         default=0.4025_rk)
     !  Stochiometric coefficients
   call self%get_parameter(self%NtoB,     'NtoB',         'uM(N)/mgWW/m3', 'N[uM]/BIOMASS [mg/m3]',                         default=0.016_rk)
   call self%get_parameter(self%OtoN,     'OtoN',         'uM(O)/uM(N)',   'Redfield (138/16) to NO3',                      default=8.625_rk)
   call self%get_parameter(self%NtoN,     'NtoN',         'uM(N)/uM(N)',   'Richards denitrification (84.8/16.)',          default=5.3_rk)

   ! Register state variables
   call self%register_state_variable(self%id_oxy,'Oxy','mmol/m**3','OXY',  150.0_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_phy,'Phy','mmol/m**3','PHY',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Wphy)
   call self%register_state_variable(self%id_nut,'NUT','mmol/m**3','NUT',  1.0_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_pom,'POM','mmol/m**3','POM',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Wpom)
   call self%register_state_variable(self%id_dom,'DOM','mmol/m**3','DOM',  0.1_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_het,'Het','mmol/m**3','HET',  0.1_rk, minimum=0.0_rk, vertical_movement=self%Whet)
   ! Register the contribution of all state variables to total nitrogen
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_phy)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_nut)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_pom)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_dom)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_het)


   ! Register diagnostic variables
call self%register_diagnostic_variable(self%id_MortHet,'MortHet','mmol/m**3/d',  'MortHET,  Mortality of Het',           &
                    output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_RespHet,'RespHet','mmol/m**3/d',  'RespHET, Respiration rate of Het',     &
                output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazPhy,'GrazPhy','mmol/m**3/d',  'GrazPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazPOM,'GrazPOM','mmol/m**3/d',  'GrazPOM',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_MortPhy,'MortPhy','mmol/m**3/d',  'MortPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_ExcrPhy,'ExcrPhy','mmol/m**3/d',  'ExcrPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimN,'LimN','nd',  'LimN',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrowthPhy,'GrowthPhy','mmol/m**3/d', 'GrowthPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimT,'LimT','nd',  'LimT',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimLight,'LimLight','nd',  'LimLight',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_DOM_decay_ox,'DOM_decay_ox','mmol/m**3/d',  'DOM_decay_ox',             &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_DOM_decay_denitr,'DOM_decay_denitr','mmol/m**3/d',  'DOM_decay_denitr', &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_POM_decay_ox,'POM_decay_ox','mmol/m**3/d',  'POM_decay_ox',             &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_POM_decay_denitr,'POM_decay_denitr','mmol/m**3/d',  'POM_decay_denitr', &
            output=output_time_step_integrated)
!call self%register_diagnostic_variable(self%id_N_fixation,'N_fixation','mmol/m**3/d',  'N_fixation',           &
!            output=output_time_step_integrated)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_windspeed,standard_variables%wind_speed)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_oxydep),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk)                   :: oxy, nut, pom, dom, phy, het, t, iopt

   real(rk) :: doxy, dnut, ddom, dpom, dphy, dhet
 ! Rates of biogeochemical processes
 ! Phy
   real(rk) :: GrowthPhy               ! Nutrient uptake rate (1/d)
   real(rk) :: Iz                      ! Irradiance at certain depth
   real(rk) :: LimLight                ! Photosynthesis dependencs on irradiance
   real(rk) :: LimT                    ! Photosynthesis dependencs on temperature
   real(rk) :: LimN                    ! Photosynthesis dependencs on nutrient
   real(rk) :: MortPhy                 ! Mortality of Phy (1/d)
   real(rk) :: ExcrPhy                 ! Excretion of Phy (1/d)
   real(rk) :: RespPhy                 ! Respiration of Phy (1/d)
 ! Het
   real(rk) :: GrazPhy, GrazPOM, RespHet, MortHet
 ! POM
   real(rk) :: Autolys                 ! Autolysis of POM to DOM (1/d)
   real(rk) :: POM_decay_ox            ! oxic mineralization of POM and ammonification (1/d)
   real(rk) :: POM_decay_denitr        ! suboxic mineralization of POM (denitrification) (1/d)
 ! DOM
   real(rk) :: DOM_decay_ox            ! oxic mineralization of POM and ammonification (1/d)
   real(rk) :: DOM_decay_denitr        ! suboxic mineralization of POM (denitrification) (1/d)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_phy,phy)
   _GET_(self%id_het,het)
   _GET_(self%id_pom,pom)
   _GET_(self%id_nut,nut)
   _GET_(self%id_dom,dom)

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
   _GET_(self%id_temp,t)              ! temperature

!--------------------------------------------------------------
! Phy
!--------------------------------------------------------------
   ! Growth of Phy and uptake of NUT
       LimLight = Iz/self%Iopt*exp(1-Iz/self%Iopt)  !Dependence on Irradiance
       LimT     = exp(self%bm*t-self%cm)            !Dependence on Temperature
       LimN     = yy(self%Knut,nut/phy)             !Dependence on nutrients
    GrowthPhy = self%Max_uptake*LimLight*LimT*LimN*phy

   ! Respiraion of Phy and increase of NUT
     RespPhy=self%r_phy_nut*phy

   ! Methabolism of Phy and increase of DOM
     ExcrPhy=self%r_phy_dom*phy

   ! Mortality of Phy and increase of POM
     MortPhy = phy*(self%r_phy_pom &
   ! Additional  mortaliny in suboxic/anoxic conditions
              +(0.5-0.5*tanh(oxy-self%O2_add_mor_phy))*phy*self%r_phy_om_anox)

!--------------------------------------------------------------
! Het
!--------------------------------------------------------------
   ! Grazing of Het on Phy
    GrazPhy = het*self%r_phy_het*yy(self%Kphy,(max(0.0_rk,phy-0.01))/(het+0.0001))

   ! Grazing of Het on POM
    GrazPOM = self%r_pop_het*het*yy(self%Kpop,(max(0.0_rk,pom-0.01))/(het+0.0001))

   ! Respiraion of HET and increase of NUT
    RespHet = self%r_het_nut*het !*(0.5+0.5*tanh(oxy-20))

   ! Mortality of HET and increase of POM
    MortHet = het*(self%r_het_pom+(0.5-0.5*tanh(oxy-15))*0.3) 

!--------------------------------------------------------------
! POM
!--------------------------------------------------------------
   ! Decomposition of POM and increase of DOM
    Autolys = self%r_pom_dom*pom

   ! Oxic mineralization of POM and ammonification depend on T
    POM_decay_ox   = self%r_pom_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*pom
   ! Suboxic mineralization of OM (denitrification and anammox), depends on T,O2,NO3/NO2
    POM_decay_denitr = self%r_pom_nut_nut*(1.+self%beta_da*yy(self%tda,t)) &
                      * (0.5-0.5*tanh(self%O2_suboxic-oxy)) &
                      * (1-tanh(1.-nut))*pom

!--------------------------------------------------------------
! DOM
!--------------------------------------------------------------
! Oxic mineralization of DOM and ammonification depend on T
    DOM_decay_ox   = self%r_dom_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*dom
! Suboxic mineralization of OM (denitrification and anammox), depends on T,O2,NO3/NO2
    DOM_decay_denitr = self%r_dom_nut_nut*(1.+self%beta_da*yy(self%tda,t)) &
                           * (0.5-0.5*tanh(self%O2_suboxic-oxy)) &
                           * (1-tanh(10.-nut))*dom

   ! Now we can summarize processes and write state variables sink/sources:

!--------------------------------------------------------------
! OXY
!--------------------------------------------------------------
   ! Changes of OXY due to OM production and decay!
    doxy = -self%OtoN*(POM_decay_ox + DOM_decay_ox &
           + (POM_decay_denitr + DOM_decay_denitr) &
           - GrowthPhy*phy + RespPhy) &
   ! additional consumption of OXY due to oxidation of reduced froms of S,Mn,Fe etc.
   ! in suboxic conditions  equales consumption for NH4 oxidation (Yakushev et al, 2008 in Nauka Kubani)
           - (0.5-0.5*tanh(oxy-self%O2_add_mor_phy))*self%OtoN*DOM_decay_ox
!--------------------------------------------------------------
! NUT
!--------------------------------------------------------------
   dnut = -GrowthPhy + RespPhy + POM_decay_ox + POM_decay_denitr + DOM_decay_ox + DOM_decay_denitr &
          - self%NtoN*(POM_decay_denitr+DOM_decay_denitr) ! Changes of NUT (as NO3+NO2) due to denitrification!

   ddom = ExcrPhy + Autolys -DOM_decay_ox - DOM_decay_denitr + (GrazPhy+GrazPOM)*(1.-self%Uz)*self%Hz
   dpom = MortPhy-Autolys-POM_decay_ox - POM_decay_denitr + (GrazPhy+GrazPOM)*(1.-self%Uz)*(1.-self%Hz)-GrazPOM
   dphy = GrowthPhy-RespPhy-ExcrPhy-MortPhy-GrazPhy
   dhet = self%Uz*(GrazPhy+GrazPOM)-MortHet-RespHet
!--------------------------------------------------------------


!derivatives for FABM
   _SET_ODE_(self%id_oxy,doxy)
   _SET_ODE_(self%id_nut,dnut)
   _SET_ODE_(self%id_dom,ddom)
   _SET_ODE_(self%id_pom,dpom)
   _SET_ODE_(self%id_phy,dphy)
   _SET_ODE_(self%id_het,dhet)
   
   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   !if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic,self%dic_per_n*dn)
   !if (_AVAILABLE_(self%id_dic)) _SET_ODE_(self%id_dic,dnut*106._rk/16._rk)
   !if (_AVAILABLE_(self%id_alk)) _SET_ODE_(self%id_alk,ddom)
   
   ! Export diagnostic variables
_SET_DIAGNOSTIC_(self%id_MortHet,MortHet)
_SET_DIAGNOSTIC_(self%id_RespHet,RespHet)
_SET_DIAGNOSTIC_(self%id_GrazPhy,GrazPhy)
_SET_DIAGNOSTIC_(self%id_GrazPOM,GrazPOM)
_SET_DIAGNOSTIC_(self%id_MortPhy,MortPhy)
_SET_DIAGNOSTIC_(self%id_ExcrPhy,ExcrPhy)
_SET_DIAGNOSTIC_(self%id_LimN,LimN)
_SET_DIAGNOSTIC_(self%id_GrowthPhy,GrowthPhy)
_SET_DIAGNOSTIC_(self%id_LimT,LimT)
_SET_DIAGNOSTIC_(self%id_LimLight,LimLight)
_SET_DIAGNOSTIC_(self%id_DOM_decay_ox,DOM_decay_ox)
_SET_DIAGNOSTIC_(self%id_DOM_decay_denitr,DOM_decay_denitr)
_SET_DIAGNOSTIC_(self%id_POM_decay_ox,POM_decay_ox)
_SET_DIAGNOSTIC_(self%id_POM_decay_denitr,POM_decay_denitr)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
!  Oxygen air-water flux. Adopted from ERSEM

! !INPUT PARAMETERS:
   class (type_niva_oxydep),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

! !LOCAL VARIABLES:
   real(rk)                   :: O2, temp, salt, windspeed
   real(rk)                   :: Ox, Oa, TempT, Obe, Q_O2

   _HORIZONTAL_LOOP_BEGIN_
    _GET_(self%id_oxy,O2)
    _GET_(self%id_temp,temp)              ! temperature
    _GET_(self%id_salt,salt)              ! salinity
    _GET_HORIZONTAL_(self%id_windspeed,windspeed)

   Ox = 1953.4-128*temp+3.9918*temp*temp-0.050091*temp*temp*temp !(Wanninkoff, 1992)
     if (Ox>0) then
       Oa = 0.028*(windspeed**3.)*sqrt(400/Ox)   !
     else
       Oa = 0.
     endif
   ! Calculation of O2 saturation Obe according to UNESCO, 1986
   TempT = (temp+273.15)/100.
   Obe = exp(-173.4292+249.6339/TempT+143.3483*log(TempT)-21.8492*TempT+salt*(-0.033096+0.014259*TempT-0.0017*TempT*TempT)) !Osat
   Obe = Obe*1000./22.4  ! convert from ml/l into uM

!  Q_O2 = Oa*(Obe-O2)*0.24 ! 0.24 is to convert from [cm/h] to [m/day]
   Q_O2 = windspeed*(Obe-O2)/86400. !After (Burchard et al., 2005)

  _SET_SURFACE_EXCHANGE_(self%id_oxy,Q_O2)

_HORIZONTAL_LOOP_END_

   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
!

! !INPUT PARAMETERS:
   class (type_niva_oxydep),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: oxy,pom,dom,phy,het,nut

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_nut,nut)
   _GET_(self%id_phy,phy)
   _GET_(self%id_het,het)
   _GET_(self%id_pom,pom)
   _GET_(self%id_dom,dom)
   _GET_(self%id_oxy,oxy)

   ! BURYING into the sediments, mmol/m2/s (sinking rates "Wxxx" are in m/s and positive upward)
   _SET_BOTTOM_EXCHANGE_(self%id_pom,self%Bu*self%Wpom*pom)
   _SET_BOTTOM_EXCHANGE_(self%id_phy,self%Bu*self%Wphy*phy)
   _SET_BOTTOM_EXCHANGE_(self%id_het,self%Bu*self%Whet*het)

   ! we use here the relaxation condition with relaxation time Trel

   ! UPWARD fluxes of dissolved parameters
   !--- independent on redox conditions
   _SET_BOTTOM_EXCHANGE_(self%id_dom,-(dom-self%b_dom_ox)/self%Trel)

   !--- dependent on redox conditions, in suboxic and anoxic conditions upward flux of DOM increases, and OXY=0 in the pore water
   _SET_BOTTOM_EXCHANGE_(self%id_dom,-(1.-0.5*(1.-tanh(self%O2_suboxic-oxy)))*(dom-self%b_dom_anox)/self%Trel)
   _SET_BOTTOM_EXCHANGE_(self%id_oxy,-(1.-0.5*(1.-tanh(self%O2_suboxic-oxy)))*(oxy-0.)/self%Trel)

   ! DOWNWARD fluxes of dissolved parameters
   !--- dependent on redox conditions: in oxic conditions OXY is additionally consumed due to its flux from water to sediments
   _SET_BOTTOM_EXCHANGE_(self%id_oxy,-(1-tanh(self%O2_suboxic-oxy))*(oxy-min(self%b_ox,oxy))/self%Trel)

   ! in oxic conditions fluxes of NO3/NO2 for denitrification in the sediments
   _SET_BOTTOM_EXCHANGE_(self%id_nut,-(1-tanh(self%O2_suboxic-oxy))*(nut-min(self%b_nut,nut))/self%Trel)

   _HORIZONTAL_LOOP_END_

   end subroutine
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Saturation function squared
!
! !INTERFACE:
   real(rk) function yy(a,x)
!
! !DESCRIPTION:
! This is a squared Michaelis-Menten type of limiter:
! \begin{equation}\label{Y}
! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
! \end{equation}
!
! !IN2PUT PARAMETERS:
   real(rk), intent(in)                :: a,x
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   yy=x**2/(a**2+x**2)

   end function yy
!EOC

   end module fabm_niva_oxydep

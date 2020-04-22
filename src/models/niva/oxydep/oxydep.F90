#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:
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
! - ZOO - heterotrophs, can consume PHY and POP,  produce DOM and DON and respirate NUT.
! - NUT - represents oxydized forms of nutrients (i.e. NO3 and NO2 for N),
!   that doesn't need additional  oxygen for nitrification.
! - DOM - is dissolved organic matter. DOM  includes all kinds of labile dissolved organic matter
!   and reduced forms of inorganic nutrients (i.e. NH4 and Urea for N).
! - POM - is particular organic matter (less labile than DOM). Temperature affects DOM and POM mineralization.
! - OXY - is dissolved oxygen.
! For the details of  OXYDEP  implemented here (actually, a previous version withour ZOO), see (Yakushev et al, 2013)
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
      type (type_state_variable_id)        :: id_oxy,id_phy,id_zoo,id_nut,id_pom,id_dom
      type (type_state_variable_id)        :: id_dic,id_alk

      type (type_dependency_id)            :: id_par,id_temp, id_salt
      type (type_horizontal_dependency_id) :: id_I_0
      !type (type_diagnostic_variable_id)   :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR

!     Model parameters
      !----Phy -----------!
       real(rk) :: Max_uptake,bm,cm,Knut, r_phy_nut, r_phy_pom,r_phy_dom,r_phy_om_anox,i_min
          !real(rk) :: Io		     = 80.          ! Optimal Irradiance at the surface
          !real(rk) :: k_Erlov      = 0.10         ! Extinction coefficient
          ! real(rk) :: Iopt         = 25 !0.25     ! Optimal irradiance
          !real(rk) :: LatLight     = 50 !57.       ! Latitude of the region			[degree]
      !----Zoo -----------!
       real(rk) :: r_phy_zoo,Kphy,r_pop_zoo,Kpop,r_zoo_nut,r_zoo_pom,Uz,Hz
      !---Organic matter mineralization---- !
       real(rk) :: r_pom_dom,r_pom_nut_oxy,r_pom_nut_nut,r_dom_nut_oxy, r_dom_nut_nut, Tda, beta_da
      ! real(rk) :: r_om_nut_sul  = 0.005    ! Specific rate of OM anoxic decay   	[1/d]
      !  Lower boundary
       real(rk) :: Bu,Trel,O2LimC,b_ox,b_dom_ox,b_dom_anox,b_nut
      ! Upper boundary	 ! for oxygen flux calculations
       real(rk) :: pvel,a0,a1,a2
      !  Stochiometric coefficients
       real(rk) :: NtoB   = 0.016    ! N[uM]/BIOMASS [mg/m3] [uM(N)/mgWW/m3]
       real(rk) :: OtoN   = 8.625    ! Redfield' (138/16) to NO3 [uM(O)/uM(N)]
      ! real(rk) :: OtoN  = -6.625    ! Redfield' (106/16) [uM(O)/uM(N)]
       real(rk) :: NtoN   = 5.3      ! Richards' denitrification (84.8/16.) [uM(N)/uM(N)]
      !---sinking-------------------------------!
       real(rk) :: w_phy  = -1.5     ! PHY sinkng rate    [m/d]
       real(rk) :: w_zoo  = -0.1     ! ZOO sinkng rate    [m/d]
       real(rk) :: w_pom  = -1.5     ! POM sinkng rate    [m/d]
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
!
!
! !INPUT PARAMETERS:
   class (type_niva_oxydep), intent(inout), target :: self
   integer,                  intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
       real(rk),parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%w_phy,     'w_phy',     'm/s', 'vertical velocity of PHY (<0 for sinking)',default=-1.5_rk,scale_factor=d_per_s)
   call self%get_parameter(self%w_zoo,     'w_zoo',     'm/s', 'vertical velocity of ZOO (<0 for sinking)',  default=-0.1_rk,scale_factor=d_per_s)
   call self%get_parameter(self%w_pom,     'w_pom',     'm/s', 'vertical velocity of POM (<0 for sinking)',          default=-1.5_rk,scale_factor=d_per_s)
  ! PHY
   call self%get_parameter(self%Max_uptake,'Max_uptake','1/d', 'Maximum nutrient uptake rate',                       default=5.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Knut,      'Knut',      'nd',  'Half-sat.const. for uptake of NUT by PHY for NUT/PHY ratio', default=0.1_rk)
   call self%get_parameter(self%bm,        'bm',      '1/gradC',     'Coefficient for uptake rate dependence on t',          default=0.12_rk)
   call self%get_parameter(self%cm,        'cm',      'nd',          'Coefficient for uptake rate dependence on t',          default=1.4_rk)
   call self%get_parameter(self%i_min,      'i_min',      'nd',      'bioshading parameter ',                                default=25._rk)
   call self%get_parameter(self%r_phy_nut,     'r_phy_nut',     '1/d',    'Specific respiration rate',                       default=0.04_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_pom,     'r_phy_pom',     '1/d',    'Specific rate of mortality',                      default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_dom,     'r_phy_dom',     '1/d',    'Specific rate of excretion',                      default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_phy_om_anox, 'r_phy_om_anox', '1/d',    'Specific rate of mortality',                      default=0.4_rk,scale_factor=d_per_s)
  ! ZOO
   call self%get_parameter(self%r_phy_zoo, 'r_phy_zoo', '1/d', 'Max.spec. rate of grazing of ZOO on PHY',                    default=2.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Kphy,      'Kphy',      'nd',  'Half-sat.const.for grazing of Zoo on Phy for Phy/Zoo ratio ',default=0.1_rk)
   call self%get_parameter(self%r_pop_zoo, 'r_pop_zoo', '1/d', ' Max.spec. rate of grazing of ZOO on POM ',                  default=0.7_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Kpop,      'Kpop',      'nd',  'Half-sat.const.for grazing of ZOO on POM for POM/ZOO ratio ',default=2._rk)
   call self%get_parameter(self%r_zoo_nut, 'r_zoo_nut', '1/d', 'Specific ZOO respiration rate ',                             default=0.02_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_zoo_pom, 'r_zoo_pom', '1/d', 'Specific ZOO mortality rate ',                               default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Uz,        'Uz',        'nd',  'Food absorbency for Zoo',                                    default=0.5_rk)
   call self%get_parameter(self%Hz,        'Hz',        'nd',  'Ratio betw. diss. and part. excretes of Zoo ',               default=0.5_rk)
  ! POM
   call self%get_parameter(self%r_pom_dom,     'r_pom_dom',     '1/d', 'Specific rate of POM decomposition ',                default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_pom_nut_oxy, 'r_pom_nut_oxy', '1/d', 'Specific rate of POM oxic decay  ',                  default=0.03_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_pom_nut_nut, 'r_pom_nut_nut', '1/d', 'Specific rate of POM denitrification  ',             default=0.006_rk,scale_factor=d_per_s)
  ! DOM
   call self%get_parameter(self%r_dom_nut_oxy, 'r_dom_nut_oxy', '1/d', 'Specific rate of DOM oxic decay  ',                  default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%r_dom_nut_nut, 'r_dom_nut_nut', '1/d', 'Specific rate of DOM denitrification  ',             default=0.002_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Tda,           'Tda',           'nd',  'Coefficient for dependence on t ',                   default=13._rk)
   call self%get_parameter(self%beta_da,       'beta_da',       'nd',  'Coefficient for dependence on t ',                   default=20._rk)
  !  Lower boundary
   call self%get_parameter(self%Bu,        'Bu',         'nd',      'Burial coeficient for lower boundary',                  default=0.25_rk)
   call self%get_parameter(self%Trel,      'Trel',       's/m',     'Relaxation time for exchange with teh sediments',       default=1e5_rk)
   call self%get_parameter(self%O2LimC,    'O2LimC',     'mmol/m3', 'Limiting O2 value for oxic/suboxic switch',             default=40._rk)
   call self%get_parameter(self%b_ox,      'b_ox',       'mmol/m3', 'OXY in the sediment',                                   default=0._rk)
   call self%get_parameter(self%b_dom_ox,  'b_dom_ox',   'mmol/m3', 'OM in the sediment (oxic conditions)',                  default=2._rk)
   call self%get_parameter(self%b_dom_anox,'b_dom_anox', 'mmol/m3', 'OM in the sediment (anoxic conditions) ',               default=6._rk)
   call self%get_parameter(self%b_nut,     'b_nut',      'mmol/m3', 'NUT in the sediment',                                   default=0._rk)
     ! Upper boundary	 ! for oxygen flux calculations
   call self%get_parameter(self%pvel,       'pvel',       'm/s',     'wind speed',                                            default=5._rk)
   call self%get_parameter(self%a0,         'a0',         'mmol/m3', 'oxygen saturation parameter  ',                         default=31.25_rk)
   call self%get_parameter(self%a1,         'a1',         'nd',      'oxygen saturation parameter  ',                         default=14.603_rk)
   call self%get_parameter(self%a2,         'a2',         '1/degC',  'oxygen saturation parameter  ',                         default=0.4025_rk)
     !  Stochiometric coefficients
   call self%get_parameter(self%NtoB,     'NtoB',         'uM(N) / mgWW/m3', 'N[uM]/BIOMASS [mg/m3]',                         default=0.016_rk)
   call self%get_parameter(self%OtoN,     'OtoN',         'uM(O)/uM(N)',    ' Redfield (138/16) to NO3',                      default=8.625_rk)
   call self%get_parameter(self%NtoN,     'NtoN',         'uM(N)/uM(N)',    '  Richards denitrification (84.8/16.)',          default=5.3_rk)
      !---sinking-------------------------------!
   !call self%get_parameter(self%w_phy,     'w_phy',         'm/d',    ' PHY sinking rate ',         default=-1.5_rk , scale_factor=d_per_s )
   !call self%get_parameter(self%w_zoo,     'w_zoo',         'm/d',    ' ZOO sinking rate',          default=-0.1_rk, scale_factor=d_per_s )
   !call self%get_parameter(self%w_pom,     'w_pom',         'm/d',    ' POM sinking rate',          default=-1.5_rk, scale_factor=d_per_s )

   ! Register state variables
   call self%register_state_variable(self%id_oxy,'oxy','mmol/m**3','Dissolved Oxygen',        150.0_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_phy,'phy','mmol/m**3','Phytoplankton',             0.1_rk, minimum=0.0_rk, vertical_movement=self%w_phy)
   call self%register_state_variable(self%id_nut,'nut','mmol/m**3','Inorganic Nutrient',        1.0_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_pom,'pom','mmol/m**3','Particulate Organic Matter',0.1_rk, minimum=0.0_rk, vertical_movement=self%w_pom)
   call self%register_state_variable(self%id_dom,'dom','mmol/m**3','Dissolved Organic Matter',  0.1_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_zoo,'zoo','mmol/m**3','Zooplankton',               0.1_rk, minimum=0.0_rk, vertical_movement=self%w_zoo)

   ! Register the contribution of all state variables to total nitrogen
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_phy)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_nut)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_pom)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_dom)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_zoo)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   !call self%register_state_dependency(self%id_dic,'dic','mmol/m**3','total dissolved inorganic carbon',required=.false.)
   !call self%register_state_dependency(self%id_alk,'alk','mmol/m**3','total alkalinity',required=.false.)

   ! Register diagnostic variables
   !call self%register_diagnostic_variable(self%id_GPP,'GPP','mmol/m**3',  'gross primary production',           &
   !                  output=output_time_step_integrated)

   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_shortwave_flux)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)

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
   real(rk)                   :: oxy,nut,pom,dom,phy,zoo,t,iopt,I_0

   real(rk) :: doxy,dnut,ddom,dpom,dphy,dzoo
 ! Rates of biogeochemical processes
 ! Phy
   real(rk) :: Growthphy               ! Nutrient uptake rate			(1/d)
   real(rk) :: Iz                      ! Irradiance at certain depth
   real(rk) :: LimLight                ! Photosynthesis dependencs on irradiance
   real(rk) :: LimT                    ! Photosynthesis dependencs on temperature
   real(rk) :: LimNut                  ! Photosynthesis dependencs on nutrient
   real(rk) :: MortPhy                 ! Mortality of Phy (1/d)
   real(rk) :: ExcrPhy                 ! Excretion of Phy (1/d)
   real(rk) :: RespPhy                 ! Respiration of Phy (1/d)
 ! Zoo
   real(rk) :: GrazPhy,GrazPOM,RespZoo,MortZoo
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
   _GET_(self%id_zoo,zoo)
   _GET_(self%id_pom,pom)
   _GET_(self%id_nut,nut)
   _GET_(self%id_dom,dom)

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation
   _GET_(self%id_temp,t)              ! temperature

!--------------------------------------------------------------
! Phy
!--------------------------------------------------------------
! Growth of Phy and uptake of NUT
          iopt=max(0.25*I_0,self%I_min)
! change "julianday"!!!  *cos((LatLight+(23.5*sin(2*3.14*julianday/365.)))*3.14/180.)
          LimLight = Iz/Iopt*exp(1-Iz/Iopt)       ! Dependence on Irradiance
          LimT     = exp(self%bm*t-self%cm)             ! Dependence on Temperature
          LimNut   = yy(self%Knut,nut/phy) !Dependence on nutrients
     Growthphy = self%Max_uptake*LimLight*LimT*LimNut*phy
!--------------------------------------------------------------
! Respiraion of phy and increase of NUT
     RespPhy=self%r_phy_nut*phy
!--------------------------------------------------------------
! Methabolism of phy and increase of DOM
     ExcrPhy=self%r_phy_dom*phy
!--------------------------------------------------------------
! Mortality of phy and increase of POM
     MortPhy = phy*(self%r_phy_pom &
!            additional  mortaliny in suboxic/anoxic conditions
                +(0.5-0.5*tanh(oxy-15))**self%r_phy_om_anox)

!--------------------------------------------------------------
! Zoo
!--------------------------------------------------------------
!    grazing of Zoo on Phy
!    if (phy>1) then
         GrazPhy = zoo*self%r_phy_zoo*yy(self%Kphy,(max(0.0_rk,phy-0.01))/(zoo+0.0001))
       !else
       !    GrazPhy=0.
       !    endif
!--------------------------------------------------------------
!    grazing of Zoo on POM
     GrazPOM = self%r_pop_zoo*zoo*yy(self%Kpop,(max(0.0_rk,pom-0.01))/(zoo+0.0001))
!--------------------------------------------------------------
! Respiraion of ZOO and increase of NUT
     RespZoo = self%r_zoo_nut*zoo !*(0.5+0.5*tanh(oxy-20))
!--------------------------------------------------------------
! Mortality of ZOO and increase of POM
     MortZoo = zoo*(self%r_zoo_pom+(0.5-0.5*tanh(oxy-15))*0.3) !+ &     (0.5+0.4*tanh(H2S-10.))*0.45)

!--------------------------------------------------------------
! POM
!--------------------------------------------------------------
! Decomposition of POM and increase of DOM
      Autolys = self%r_pom_dom*pom
!--------------------------------------------------------------
! Oxic mineralization of POM and ammonification depend on T
   POM_decay_ox   = self%r_pom_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*pom
! Suboxic mineralization of OM (denitrification and anammox), depends on T,O2,NO3/NO2
   POM_decay_denitr = self%r_pom_nut_nut*(1.+self%beta_da*yy(self%tda,t)) &
                           * (0.5-0.5*tanh(self%O2LimC-oxy)) &
                           * (1-tanh(1.-nut))*pom
! Mineralization of OM, ammonification and growth of NUT

!--------------------------------------------------------------
! DOM
!--------------------------------------------------------------
! Oxic mineralization of DOM and ammonification depend on T
   DOM_decay_ox   = self%r_dom_nut_oxy*(1.+self%beta_da*yy(self%tda,t))*dom
! Suboxic mineralization of OM (denitrification and anammox), depends on T,O2,NO3/NO2
   DOM_decay_denitr = self%r_dom_nut_nut*(1.+self%beta_da*yy(self%tda,t)) &
                           * (0.5-0.5*tanh(self%O2LimC-oxy)) &
                           * (1-tanh(10.-nut))*dom
! Mineralization of OM, ammonification and growth of NUT

! Now we can summarize processes and write state variables sink/sources:
!--------------------------------------------------------------
! OXY
!--------------------------------------------------------------
! Changes of OXY due to OM production and decay!
   doxy = -self%OtoN*(POM_decay_ox + DOM_decay_ox +  &
           (POM_decay_denitr + DOM_decay_denitr) &
           - GrowthPhy*phy + RespPhy)
! additional consumption of OXY due to oxidation of reduced froms of S,Mn,Fe etc.
! in suboxia equales condumption for NH4 oxidation (Yakushev et al, 2008 in Nauka Kubani)
!                     + dd(dom,nut,ci)*(1.-0.5*(1.-tanh(O2LimC-oxy)))
   if (oxy.lt.30.) doxy = doxy-self%OtoN*DOM_decay_ox
!--------------------------------------------------------------
! NUT
!--------------------------------------------------------------
   dnut = -Growthphy + RespPhy + POM_decay_ox + POM_decay_denitr + DOM_decay_ox + DOM_decay_denitr &
          - self%NtoN*(POM_decay_denitr+DOM_decay_denitr) ! Changes of NUT (as NO3+NO2) due to denitrification!

!--------------------------------------------------------------
   ddom = ExcrPhy + Autolys -DOM_decay_ox - DOM_decay_denitr + (GrazPhy+GrazPOM)*(1.-self%Uz)*self%Hz
   dpom = MortPhy-Autolys-POM_decay_ox - POM_decay_denitr + (GrazPhy+GrazPOM)*(1.-self%Uz)*(1.-self%Hz)-GrazPOM
   dphy = GrowthPhy-RespPhy-ExcrPhy-MortPhy-GrazPhy
   dzoo = self%Uz*(GrazPhy+GrazPOM)-MortZoo-RespZoo
!--------------------------------------------------------------
   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   !if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic,self%dic_per_n*dn)

   ! Export diagnostic variables
   !_SET_DIAGNOSTIC_(self%id_dPAR,par)
!derivatives for FABM
   _ADD_SOURCE_(self%id_oxy,doxy)
   _ADD_SOURCE_(self%id_nut,dnut)
   _ADD_SOURCE_(self%id_dom,ddom)
   _ADD_SOURCE_(self%id_pom,dpom)
   _ADD_SOURCE_(self%id_phy,dphy)
   _ADD_SOURCE_(self%id_zoo,dzoo)

   !if (_AVAILABLE_(self%id_dic)) _ADD_SOURCE_(self%id_dic,dnut*106._rk/16._rk)
   !if (_AVAILABLE_(self%id_alk)) _ADD_SOURCE_(self%id_alk,ddom)

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
!
! !INPUT PARAMETERS:
   class (type_niva_oxydep),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

! !LOCAL VARIABLES:
   real(rk)                   :: O2, temp,salt
   real(rk)                   :: Ox,Oa,TempT,Obe, Q_O2

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_oxy,O2)
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_salt,salt)              ! salinity

!/*---------------------------------------------------O2 exchange with air */
   Ox = 1800.6-120.1*temp+3.7818*temp*temp &
       -0.047608*temp*temp*temp !Ox=Sc, Schmidt number
   if (Ox>0) then
    !    Oa = 0.028*7.6*7.6*7.6*sqrt(660/Ox)   ! Pvel for the Baltic Sea by Schneider
      Oa = 0.028*6.*6.*6.*sqrt(660/Ox)       ! Pvel for the Black Sea
   else
      Oa = 0.
   end if

      ! Calculation of O2 saturation Obe according to UNESCO, 1986
   TempT = (temp+273.15)/100.
   Obe = exp(-173.4292+249.6339/TempT+143.3483*log(TempT)-21.8492*TempT+salt*(-0.033096+0.014259*TempT-0.0017*TempT*TempT)) !Osat
   Obe = Obe*1000./22.4  ! - in uM

   Q_O2 = Oa*(Obe-O2)*0.24/86400. ! 0,24 is to transform from [cm/h] to [m/day]

   _ADD_SURFACE_FLUX_(self%id_oxy,Q_O2)

   _HORIZONTAL_LOOP_END_

   end subroutine
!-----------------------------------------------------------------------


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
   real(rk)                   :: t,oxy,pom,dom,phy,zoo,nut

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_nut,nut)
   _GET_(self%id_phy,phy)
   _GET_(self%id_zoo,zoo)
   _GET_(self%id_pom,pom)
   _GET_(self%id_dom,dom)
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_temp,t)              ! temperature

      !----------------------------------------
      !-burying into the sediments (sinking rates "w_xxx" are in m/s and positive upward)
   _ADD_BOTTOM_FLUX_(self%id_pom,self%Bu*self%w_pom*pom)   ! mmol/m2/s
   _ADD_BOTTOM_FLUX_(self%id_phy,self%Bu*self%w_phy*phy)   ! mmol/m2/s
   _ADD_BOTTOM_FLUX_(self%id_zoo,self%Bu*self%w_zoo*zoo)   ! mmol/m2/s
      !   cc(pon,1)=cc(pon,1)+w_pon/h(1)*cc(pon,1)*Bu from ROLM

      ! we use here the relaxation condition with relaxation time Trel
      !---------------------------------------- upward fluxes of dissolved parameters
      !---------------------------------------- independent on redox conditions
   _ADD_BOTTOM_FLUX_(self%id_dom,-(dom-self%b_dom_ox)/self%Trel)   ! mmol/m2/s

      !---------------------------------------- dependent on redox conditions
      !in suboxic and anoxic conditions upward flux of DOM increases
   _ADD_BOTTOM_FLUX_(self%id_dom,-(1.-0.5*(1.-tanh(self%O2LimC-oxy)))*(dom-self%b_dom_anox)/self%Trel)
   _ADD_BOTTOM_FLUX_(self%id_oxy,-(1.-0.5*(1.-tanh(self%O2LimC-oxy)))*(oxy-0.)/self%Trel)

!    cc(oxy,1) =cc(oxy,1) -(1-tanh(O2LimC-cc(oxy,1)))*(cc( oxy,1)-min(b_ox,cc(oxy,1))) &
!				/(Trel*timestep) !/86400)

!---------------------------------------- dependent on redox conditions downward fluxes

! if (h(1).gt.0.75) then ! that is NOT valid in the channel (<5 m)!
!! in oxic conditions OXY is additionally consumed due to its flux from water to sediments
   _ADD_BOTTOM_FLUX_(self%id_oxy,-(1-tanh(self%O2LimC-oxy))*(oxy-min(self%b_ox,oxy))/self%Trel)
!				/(Trel*timestep) !/86400)
!! in oxic conditions fluxes of NO3/NO2 for denitrification in the sediments
   _ADD_BOTTOM_FLUX_(self%id_nut,-(1-tanh(self%O2LimC-oxy))*(nut-min(self%b_nut,nut))/self%Trel)
!				/(Trel*timestep) !/86400)
! endif

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

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

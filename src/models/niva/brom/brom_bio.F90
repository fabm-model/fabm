#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_bio
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev, Elizaveta Protsenko, Jorn Bruggeman
!

! !PUBLIC DERIVED TYPES:
     type,extends(type_base_model),public :: type_niva_brom_bio
!     Variable identifiers
     type (type_state_variable_id)        :: id_Phy,id_Het
     type (type_state_variable_id)        :: id_O2,id_NO2,id_NO3,id_NH4,id_PO4
     type (type_state_variable_id)        :: id_PON,id_DON
     type (type_state_variable_id)        :: id_Baae,id_Baan,id_Bhae,id_Bhan
     type (type_state_variable_id)        :: id_DIC,id_H2S,id_Si,id_Sipart,id_Alk
     type (type_dependency_id)            :: id_temp,id_salt,id_par
     type (type_dependency_id)            :: id_Kp1,id_Kp2,id_Kp3,id_Knh4,id_KSi,id_Hplus      
     type (type_horizontal_dependency_id) :: id_windspeed
    !type (type_diagnostic_variable_id)   :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR
    type (type_diagnostic_variable_id)   :: id_MortHet,id_Grazing,id_RespHet,id_GrazBhae,id_GrazBhan,id_GrazBaae,id_GrazBaan
    type (type_diagnostic_variable_id)   :: id_GrazPhy,id_GrazPOP,id_GrazBact,id_MortPhy,id_ExcrPhy,id_LimNH4,id_LimN,id_GrowthPhy
    type (type_diagnostic_variable_id)   :: id_LimT,id_LimP,id_LimNO3,id_LimLight!,id_Iz
      !===========================================================================!
      !     Model parameters
      !-------------------------------------------------------------------------
      ! specific rates of biogeochemical processes
 
       !----Phy  ----------!
      real(rk) :: KNF,k_Erlov,kc,Io,Iopt,bm,cm,KFN,KFP,KFD,LatLight
       !----Het -----------!
      real(rk) :: KFZ,KFY,KPZ,KPP,KZN,KZP,Uz,Hz,limGrazBac       
       !---- O2--------! 
       ! Upper boundary	 ! for oxygen flux calculations
       real(rk) :: pvel          = 5.       ! wind speed  			      	     [m/s]
       real(rk) :: a0            = 31.25    !oxygen saturation Parameter         	      [uM]
       real(rk) :: a1            = 14.603   !oxygen saturation Parameter         	      [ -]
       real(rk) :: a2            = 0.4025   !oxygen saturation Parameter         	  [1/degC]
       !---- N---------!
      real(rk) ::  KPSI=1.46 !    % ! Strength of NH4 inhibition of NO3 uptake constant =1.46(Wroblev,G,01)
      real(rk) ::  KNO3=0.15 !    %! Half-sat.const.for uptake of NO3+NO2 =0.5 (G,01)
      real(rk) ::  KNH4=0.02 !    %! Half-sat.const.for uptake of NH4=0.2 (G,01)

       !---- P---------!
      real(rk) ::  KPO4=0.02      ! !Half-sat. constant for uptake of PO4 by Phy 
       !---- Sinking---!      
      real(rk) :: Wsed= 5. !1Rate of sinking of detritus (POP, PON)d-1 !!  Wdetr=1.5 (Savchuk, Wulff,1996),!Wdetr= 3.5; 20. (Gregoire,2000)
      real(rk) :: Wphy=0.1 !Rate of sinking of Phy m d-1 !! Wphy=1 Savchuk, Wulff,1996)
      real(rk) :: Whet=1.0 !Rate of sinking of Het m d-1 
 
       !---- Stoichiometric coefficients ----!
      real(rk) ::   Sp=0.001     !P[uM]/BIOMASS [mg/m3]
      real(rk) ::   SpZ=0.001    !P[uM]/BIOMASS Het [mg/m3]
      real(rk) ::   Sn=0.016     !N[uM]/BIOMASS [mg/m3]
      real(rk) ::   Ss=0.053     !S[uM]/BIOMASS [mg/m3] of bacteria, i.e. during ChemBaaeynthesis
      real(rk) ::   Sc=0.106     !C[uM]/BIOMASS [mg/m3] for bact should be *.5 due to large P content
      real(rk) ::   OkP=106.     !O[uM]/P[uM]
      real(rk) ::   NkP=16.      !N[uM]/P[uM]
      real(rk) ::   OkN=6.625    !O[uM]/N[uM]
      real(rk) ::   SkP=53.      !S[uM]/P[uM] 
      real(rk) ::   CkN=8.       !C[uM]/N[uM]  
      real(rk) ::   SikN=2.      !Si[uM]/N[uM] (Ivanenkov, 1979)      
!!===========================================================================!   
      
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
   class (type_niva_brom_bio), intent(inout), target :: self
   integer,                    intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
!-----------------------------------------------------------------------   
! Parameters, i.e. rate constants    
       !----Phy  ----------!
    !   call self%get_parameter(self%Wsed,'Wsed','m/d',' Rate of sinking of detritus (POP, PON)',default=5.0_rk)
   call self%get_parameter(self%KNF,     'KNF',     '1/d',         'Maximum specific growth rate',                    default=2.0_rk)  
   call self%get_parameter(self%k_Erlov, 'k_Erlov', '1/m',         'Extinction coefficient',                          default=0.05_rk)
   call self%get_parameter(self%Io,      'Io',      'Watts/m**2/h','Surface irradiance',                              default=80.0_rk)
   call self%get_parameter(self%Iopt,    'Iopt',    'Watts/m**2/h','Optimal irradiance',                              default=25.0_rk)
   call self%get_parameter(self%bm,      'bm',      '1/gradC',     'Coefficient for growth dependence on t',          default=0.12_rk)
   call self%get_parameter(self%cm,      'cm',      'nd',          'Coefficient for growth dependence on t',          default=1.4_rk)
   call self%get_parameter(self%kc,      'kc',      'm2/mmol N',   'Attenuation constant for the self shading effect',default=0.03_rk) !(Burchard 05)
   call self%get_parameter(self%KFN,     'KFN',     '1/d',         'Specific respiration rate',                       default=0.04_rk)
   call self%get_parameter(self%KFP,     'KFP',     '1/d',         'Specific rate of mortality',                      default=0.10_rk)
   call self%get_parameter(self%KFD,     'KFD',     '1/d',         'Specific rate of excretion',                      default=0.01_rk)
   call self%get_parameter(self%LatLight,'LatLight','degree',      'Latitude',                                        default=50.0_rk)
   
       !----Het  ----------!
    call self%get_parameter(self%KFZ,'KFZ','1/d','Max.spec. rate of grazing of Het on Phy',default=1.0_rk)
    call self%get_parameter(self%KFY,'KFY','nd','Half-sat.const.for grazing of Het on Phy for Phy/Het ratio',default=1.1_rk)
    call self%get_parameter(self%KPZ,'KPZ','mmol/m**3','Max.spec.rate of grazing of Het on POP',default=0.70_rk)
    call self%get_parameter(self%KPP,'KPP','nd','Half-sat.const.for grazing of Het on POP for POP/Het  ratio',default=0.2_rk)
    call self%get_parameter(self%KZN,'KZN','1/d','Specific respiration rate',default=0.02_rk)
    call self%get_parameter(self%KZP,'KZP','1/d','Maximum specific rate of mortality of Het',default=0.05_rk)
    call self%get_parameter(self%Uz,'Uz','nd','Food absorbency for Het',default=0.5_rk)
    call self%get_parameter(self%Hz,'Hz','nd','Ratio betw. diss. and part. excretes of Het ',default=0.5_rk)  
    call self%get_parameter(self%limGrazBac,'limGrazBac','mmol/m**3','Limiting parameter for bacteria grazing by Het ',default=2._rk)       
   !
      !!!!!!real(rk) :: Wsed= 5. !1Rate of sinking of detritus (POP, PON)d-1 !!  Wdetr=1.5 (Savchuk, Wulff,1996),!Wdetr= 3.5; 20. (Gregoire,2000)
      !!!!!!real(rk) :: Wphy=0.1 !Rate of sinking of Phy m d-1 !! Wphy=1 Savchuk, Wulff,1996)
      !!!!!!real(rk) :: Whet=1.0 !Rate of sinking of Het m d-1 
!-----------------------------------------------------------------------
! Register state variables
   call self%register_state_variable(self%id_Phy, 'Phy', 'mmol/m**3','Phy', minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_variable(self%id_Het, 'Het', 'mmol/m**3','Het',minimum=0.0_rk,vertical_movement=-self%Whet/86400._rk)
   call self%register_state_variable(self%id_NO3, 'NO3', 'mmol/m**3','NO3', minimum=0.0_rk)
   call self%register_state_variable(self%id_PO4, 'PO4', 'mmol/m**3','PO4', minimum=0.0_rk)
   call self%register_state_variable(self%id_NH4, 'NH4', 'mmol/m**3','NH4', minimum=0.0_rk)
   call self%register_state_variable(self%id_PON, 'PON', 'mmol/m**3','PON', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_DON, 'DON', 'mmol/m**3','DON', minimum=0.0_rk)
   call self%register_state_variable(self%id_O2,  'O2',  'mmol/m**3','O2', minimum=0.0_rk)
!-----------------------------------------------------------------------
   ! Register the contribution of all state variables to total nitrogen
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_bio)
!-----------------------------------------------------------------------
   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   !call self%register_state_dependency(self%id_dic,'dic','mmol/m**3','total dissolved inorganic carbon',required=.false.)
   call self%register_state_dependency(self%id_H2S, 'H2S', 'mmol/m**3','H2S')
   call self%register_state_dependency(self%id_NO2, 'NO2', 'mmol/m**3','NO2')
   call self%register_state_dependency(self%id_DIC, 'DIC', 'mmol/m**3','DIC')
   call self%register_state_dependency(self%id_Baae, 'Baae', 'mmol/m**3','aerobic autotrophic bacteria')
   call self%register_state_dependency(self%id_Bhae, 'Bhae', 'mmol/m**3','aerobic heterotrophic bacteria')
   call self%register_state_dependency(self%id_Baan, 'Baan', 'mmol/m**3','anaerobic aurotrophic bacteria')
   call self%register_state_dependency(self%id_Bhan, 'Bhan', 'mmol/m**3','anaerobic heterotrophic bacteria')
   call self%register_state_dependency(self%id_Si,   'Si', 'mmol/m**3', 'Si')
   call self%register_state_dependency(self%id_Sipart,   'Sipart', 'mmol/m**3', 'Si particulate')
   call self%register_state_dependency(self%id_Alk,  'Alk','mmol/m**3','Alk')
!-----------------------------------------------------------------------   
   ! Register diagnostic variables
   !call self%register_diagnostic_variable(self%id_GPP,'GPP','mmol/m**3',  'gross primary production',           &
   !                  output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_MortHet,'MortHet','mmol/m**3',  'Mortality of Het',           &
                    output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_Grazing,'Grazing','mmol/m**3',  'Grazing of Het',           &
                    output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_RespHet,'RespHet','mmol/m**3',  'Respiration rate of Het',           &
                output=output_time_step_integrated)   
call self%register_diagnostic_variable(self%id_GrazBhae,'GrazBhae','mmol/m**3',  'GrazBhae',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazBhan,'GrazBhan','mmol/m**3',  'GrazBhan',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazBaae,'GrazBaae','mmol/m**3',  'GrazBaae',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazBaan,'GrazBaan','mmol/m**3',  'GrazBaan',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazPhy,'GrazPhy','mmol/m**3',  'GrazPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazPOP,'GrazPOP','mmol/m**3',  'GrazPOP',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrazBact,'GrazBact','mmol/m**3',  'GrazBact',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_MortPhy,'MortPhy','mmol/m**3',  'MortPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_ExcrPhy,'ExcrPhy','mmol/m**3',  'ExcrPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimNH4,'LimNH4','mmol/m**3',  'LimNH4',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimN,'LimN','mmol/m**3',  'LimN',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_GrowthPhy,'GrowthPhy','mmol/m**3',  'GrowthPhy',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimT,'LimT','mmol/m**3',  'LimT',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimP,'LimP','mmol/m**3',  'LimP',           &
            output=output_time_step_integrated)
call self%register_diagnostic_variable(self%id_LimNO3,'LimNO3','mmol/m**3',  'LimNO3',           &
            output=output_time_step_integrated)

call self%register_diagnostic_variable(self%id_LimLight,'LimLight','mmol/m**3',  'LimLight',           &
            output=output_time_step_integrated)
!call self%register_diagnostic_variable(self%id_Iz,'Iz','mmol/m**3',  'Iz',           &
!            output=output_time_step_integrated)

!-----------------------------------------------------------------------
   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_windspeed,standard_variables%wind_speed)
!-----------------------------------------------------------------------   
   call self%register_dependency(self%id_Hplus, 'Hplus', 'mmol/m**3','H+ hydrogen')
   ! Register dependencies on equilibrium constants for phosphoric acid, ammonia.
   call self%register_dependency(self%id_Kp1,  'Kp1',  '-', '[H+][H2PO4-]/[H3PO4]')
   call self%register_dependency(self%id_Kp2,  'Kp2',  '-', '[H][HPO4]/[H2PO4]')
   call self%register_dependency(self%id_Kp3,  'Kp3',  '-', '[H][PO4]/[HPO4]')
   call self%register_dependency(self%id_Knh4, 'Knh4', '-', '[H+][NH3]/[NH4]')
   call self%register_dependency(self%id_KSi,  'KSi','-','[H+][H3SiO4-]/[Si(OH)4]')

   ! Specify that are rates computed in this module are per day (default: per second)
   self%dt = 86400.

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
   class (type_niva_brom_bio),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk) :: NH4,NO2,NO3,PO4,Phy,Het,H2S,O2,Baae,Baan,Bhae,Bhan,PON,Si,Sipart,Alk,Hplus
   real(rk) :: Knh4,Kp1,Kp2,Kp3,Kh2s1,KSi
   real(rk) :: temp,salt,Iz
   real(rk) :: LimLight,LimT,LimP,LimNO3,LimNH4,LimN,GrowthPhy,MortPhy,ExcrPhy,dAlk
   real(rk) :: GrazPhy,GrazPOP,GrazBaae,GrazBaan,GrazBhae,GrazBhan,GrazBact,Grazing,RespHet,MortHet
   !real(rk) :: KNF,k_Erlov,Io,Iopt,bm,cm,kc,KFN,KFP,KFD,LatLight
   integer(4) :: rer
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_NH4,NH4)
   _GET_(self%id_NO2,NO2)
   _GET_(self%id_NO3,NO3)
   _GET_(self%id_PO4,PO4)
   _GET_(self%id_Phy,Phy)
   _GET_(self%id_Het,Het)
   _GET_(self%id_H2S,H2S)
   _GET_(self%id_O2,O2)
   _GET_(self%id_Baae,Baae)
   _GET_(self%id_Baan,Baan)
   _GET_(self%id_Bhae,Bhae)
   _GET_(self%id_Bhan,Bhan)
   _GET_(self%id_PON,PON)
   _GET_(self%id_Si,Si)
   _GET_(self%id_Sipart,Sipart)
   _GET_(self%id_Alk,Alk)
   _GET_(self%id_Hplus,Hplus)
   
   ! Get equilibrium constants
   _GET_(self%id_Kp1,  Kp1)
   _GET_(self%id_Kp2,  Kp2)
   _GET_(self%id_Kp3,  Kp3)
   _GET_(self%id_Knh4, Knh4)
   _GET_(self%id_KSi, KSi)
   
   ! Retrieve current environmental conditions.
   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! temperature
!    _GET_(self%KNF,self%KNF)  
   
!%%
!%!---------------------------------------------------------------------------
!%!========P=H=Y==============================================================
!%!---------------------------------------------------------------------------
 
     LimLight    =Iz/self%Iopt*exp(1-Iz/self%Iopt)        ! Influence of the Irradiance on photosynthesis

     LimT        = exp(self%bm*temp-self%cm)  ! 0.2+0.22*(exp(bm*t(ci))-cm)/(cm+0.28*exp(bm*t(ci))) !exp(bm*t(ci)-cm)  ! Influence of Temperature on photosynthesis

!%! dependence of photosynthesis on P
!%     LimP          =th(cc(po4,ci),s_pho_po4,_ZERO_,_ONE_)* &
!%                    ((cc(po4, ci))/cc(phy, ci))/(cc(po4, ci)/cc(phy, ci)+KPO4/25.) 
     LimP   = & !%0.0001 + & %Cc(k,i_PO4)*Cc(k,i_PO4)/(KPO4*KPO4+Cc(k,i_PO4)*Cc(k,i_PO4)) * &
           yy(self%KPO4*self%NkP,PO4/Phy)
       
!%            (Cc(k,i_PO4)/Phy)/(Cc(k,i_PO4)/Phy+KPO4/1.*NkP) 

!%! dependence of photosynthesis on cc(nO3, ci),cc(nO2, ci):
!%     LimNO3        =th((cc(no3, ci)+cc(no2, ci)),s_pho_noX,_ZERO_,_ONE_)* &
!%                    (cc(no3, ci)+cc(no2, ci))/cc(phy, ci)/((cc(no3, ci)+cc(no2,ci))/cc(phy, ci)+KNO3/200.)
     LimNO3 = & !%0.0001 + & th((cc(no3, ci)+cc(no2, ci)),s_pho_noX,_ZERO_,_ONE_)* &
           yy(self%KNO3,(NO3+NO2)/Phy)
!% (Cc(k,i_NO3)+Cc(k,i_NO2))/Phy/((Cc(k,i_NO3)+Cc(k,i_NO2))/Phy+KNO3/1.)

!%! dependence of photosynthesis on cc(nh4, ci):
!%     LimNH4        =th(cc(nh4, ci),s_pho_nh4,_ZERO_,_ONE_)*&
!%                    ((cc(nh4, ci))/cc(phy, ci))/(cc(nh4, ci)/cc(phy, ci)+KNH4/200.)
     LimNH4 = & !%0.0001 + & %th(cc(nh4, ci),s_pho_nh4,_ZERO_,_ONE_)*&
                    yy(self%KNH4,NH4/Phy)
!%                    Cc(k,i_NH4)/Phy/(Cc(k,i_NH4)/Phy+KNH4/1.)

     LimN          =LimNO3+LimNH4  ! Influence of N on photosynthesis

 !    LimSi   = yy(1.,Cc(k,i_Si)/Phy)
!%!----------------------------------------     
     GrowthPhy     =self%KNF*LimLight*LimT*min(LimP,LimN)*Phy
!%!----------------------------------------
!%     MortPhy  =(0.8* (0.5*(1.-tanh(4-H2S))) +KFP)*Phy % ! Rate of mortality of phy
!  % ! Rate of mortality of phy
     MortPhy=(self%KFP+(0.5-0.5*tanh(O2-60.))*0.45+ &
            (0.5-0.5*tanh(O2-20))*0.45)*Phy
    
     ExcrPhy=self%KFD*Phy ! Excretion of phy


!%!---------------------------------------------------------------------------
!%!========Z=O=O==============================================================
!%!---------------------------------------------------------------------------
!%!---------------------------------------- Grazing of Het on phy, detritus and bacteria
     GrazPhy = self%KFZ*Het*yy(self%KFY,Phy/(Het+0.0001)) 

     GrazPOP = self%KPZ*Het*yy(self%KPP,PON/(Het+0.0001)) 

     GrazBaae  =1.0*self%KPZ*Het*yy(self%limGrazBac,Baae/(Het+0.0001)) 
     GrazBaan =0.5*self%KPZ*Het*yy(self%limGrazBac,Baan/(Het+0.0001)) 
     GrazBhae  =1.0*self%KPZ*Het*yy(self%limGrazBac, Bhae/(Het+0.0001)) 
     GrazBhan =1.3*self%KPZ*Het*yy(self%limGrazBac, Bhan/(Het+0.0001)) 
     
     GrazBact =GrazBaae+GrazBaan+GrazBhae+GrazBhan
!%!----------------------------------------
     Grazing = GrazPhy+GrazPOP+GrazBact !%+GrazBact
!%!----------------------------------------
    RespHet =self%KZN*Het*(0.5+0.5*tanh(O2-20))
        MortHet=(0.25+(0.5-0.5*tanh(O2-20))*0.3+ &
            (0.5+0.4*tanh(H2S-10.))*0.45)*Het

        
                    !+((Dc_OM_total-ChemBaae-ChemBaan)/self%NkP &
                    !+(fe_rd-Fe_ox-Fe_ox2+4.*DcDM_Fe+4.*DcPM_Fe)/2.7-(mn_ox2+mn_rd2-mn_ox-mn_rd)/0.67) 
!%!----------------------------------------
         dAlk =   & !  the ‘nutrient-H+-compensation principle’. Formulated by Wolf-Gladrow et al., 2007 :
!        " (i) an increase of alkalinity by 1 mole when nitrate or nitrite is the N source, 
                    + 1.*GrowthPhy*(LimNO3/LimN)  & ! decrease of H+ to compensate NO3 consumption
!        (ii) to a !decrease of alkalinity by 1 mole when ammonia is used"
                    -1.*GrowthPhy*(LimNH4/LimN)    ! &
                    !+ RespHet*Knh4/(Knh4+Hplus) &           ! +/- NH3t &
                    !+ ((-GrowthPhy+RespHet)/self%NkP)&                                        ! +/- PO4t 
                    !    *((Kp1*Kp2-Hplus*Hplus)*Hplus+2e0*Kp1*Kp2*Kp3)/(((Hplus+Kp1) &       ! +/- PO4t 
                    !    *Hplus+Kp1*Kp2)*Hplus+Kp1*Kp2*Kp3)  &                                ! +/- PO4t 
                    !+((-GrowthPhy+ExcrPhy)*self%SikN)*KSi/(KSi+Hplus)       !  ! +/- Si (tot)  
               !     + (GrowthPhy*(LimNH4/LimN)-RespHet)*(-1.))
                   ! (-GrowthPhy+RespHet)/self%NkP)&
                   ! !    *((Kp1*Kp2-Hplus*Hplus)*Hplus+2e0*Kp1*Kp2*Kp3)/(((Hplus+Kp1) & ! +/- PO4t 
                    !    *Hplus+Kp1*Kp2)*Hplus+Kp1*Kp2*Kp3)  &  ! +/- PO4t         
        
!%!----------------------------------------
! components of temporal derivarives calculated in this module:
          _SET_ODE_(self%id_Phy,(GrowthPhy-MortPhy-ExcrPhy-GrazPhy))
          _SET_ODE_(self%id_Het,(self%Uz*Grazing-MortHet-RespHet))
          _SET_ODE_(self%id_O2,(GrowthPhy-RespHet)*self%OkN)
          _SET_ODE_(self%id_DON,+ExcrPhy+Grazing*(1.-self%Uz)*self%Hz)
          _SET_ODE_(self%id_PON,MortPhy+MortHet+Grazing*(1.-self%Uz)*(1.-self%Hz)-GrazPOP) 
          _SET_ODE_(self%id_NH4,+RespHet-GrowthPhy*(LimNH4/LimN)) 
          _SET_ODE_(self%id_NO2,-GrowthPhy*(LimNO3/LimN)*(NO2/(0.00001+NO2+NO3))) 
          _SET_ODE_(self%id_NO3,-GrowthPhy*(LimNO3/LimN)*((NO3+0.00001)/(0.00001+NO2+NO3))) 
          _SET_ODE_(self%id_DIC,(-GrowthPhy+RespHet)*self%CkN)
          _SET_ODE_(self%id_PO4,(-GrowthPhy+RespHet)/self%NkP) 
          _SET_ODE_(self%id_Si,(-GrowthPhy+ExcrPhy)*self%SikN)    
          _SET_ODE_(self%id_Sipart,(MortPhy+GrazPhy)*self%SikN)  
          _SET_ODE_(self%id_Alk,dAlk)          
          _SET_ODE_(self%id_Baae,-GrazBaae)
          _SET_ODE_(self%id_Baan,-GrazBaan)
          _SET_ODE_(self%id_Bhae,-GrazBhae)
          _SET_ODE_(self%id_Bhan,-GrazBhan)   
!%!----------------------------------------
          
_SET_DIAGNOSTIC_(self%id_MortHet,MortHet)
_SET_DIAGNOSTIC_(self%id_Grazing,Grazing)
_SET_DIAGNOSTIC_(self%id_RespHet,RespHet)
_SET_DIAGNOSTIC_(self%id_GrazBhae,GrazBhae)
_SET_DIAGNOSTIC_(self%id_GrazBhan,GrazBhan)
_SET_DIAGNOSTIC_(self%id_GrazBaae,GrazBaae)
_SET_DIAGNOSTIC_(self%id_GrazBaan,GrazBaan)
_SET_DIAGNOSTIC_(self%id_GrazPhy,GrazPhy)
_SET_DIAGNOSTIC_(self%id_GrazPOP,GrazPOP)
_SET_DIAGNOSTIC_(self%id_GrazBact,GrazBact)
_SET_DIAGNOSTIC_(self%id_MortPhy,MortPhy)
_SET_DIAGNOSTIC_(self%id_ExcrPhy,ExcrPhy)
_SET_DIAGNOSTIC_(self%id_LimNH4,LimNH4)
_SET_DIAGNOSTIC_(self%id_LimN,LimN)
_SET_DIAGNOSTIC_(self%id_GrowthPhy,GrowthPhy)

_SET_DIAGNOSTIC_(self%id_LimT,LimT)
_SET_DIAGNOSTIC_(self%id_LimP,LimP)
_SET_DIAGNOSTIC_(self%id_LimNO3,LimNO3)

_SET_DIAGNOSTIC_(self%id_LimLight,LimLight)
!_SET_DIAGNOSTIC_(self%id_Iz,Iz)

!if ( L__ == 51) then
!    continue
!end if

_LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: 
!!
!! !INTERFACE:
!   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!!
!! !DESCRIPTION:
!! 
!! !INPUT PARAMETERS:
!   class (type_niva_brom_bio),intent(in) :: self
!   _DECLARE_ARGUMENTS_DO_SURFACE_
!!
!! !LOCAL VARIABLES:
! 
!   _HORIZONTAL_LOOP_BEGIN_
!   _HORIZONTAL_LOOP_END_
!   
!   end subroutine 

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
   class (type_niva_brom_bio),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk)                   :: O2, temp, salt, windHt, windspeed 
   real(rk)                   :: Ox, Oa, TempT, Obe, Q_O2
   integer(2) rer 
   _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_O2,O2)
      _GET_(self%id_temp,temp)              ! temperature
      _GET_(self%id_salt,salt)              ! salinity
      _GET_HORIZONTAL_(self%id_windspeed,windspeed)
      
!/*---------------------------------------------------O2 exchange with air */
      windHt=5.
      Ox = 1953.4-128*temp+3.9918*temp*temp-0.050091*temp*temp*temp !(Wanninkoff, 1992)
      if (Ox>0) then 
    !    Oa = 0.028*7.6*7.6*7.6*sqrt(660/Ox)   ! Pvel for the Baltic Sea by Schneider
	        Oa = 0.028*(windspeed**3.)*sqrt(400/Ox)   ! 
        else
            Oa = 0.
      endif
! Calculation of O2 saturation Obe according to UNESCO, 1986
  TempT = (temp+273.15)/100.
  Obe = exp(-173.4292+249.6339/TempT+143.3483*log(TempT)-21.8492*TempT+salt*(-0.033096+0.014259*TempT-0.0017*TempT*TempT)) !Osat
  Obe = Obe*1000./22.4  ! - in uM

  Q_O2 = Oa*(Obe-O2)*0.24/windHt! 0.24 is to convert from [cm/h] to [m/day]  

 _SET_SURFACE_EXCHANGE_(self%id_O2,Q_O2)

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
   class (type_niva_brom_bio),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   
   _HORIZONTAL_LOOP_BEGIN_
   _HORIZONTAL_LOOP_END_
   
   end subroutine 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!---------------------------------------------------------------------- 
! SUBROUTINE: Saturation function squared 
!-----------------------------------------------------------------------
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
! !INPUT PARAMETERS:
!  REALTYPE, intent(in) 
!  real :: a
   real(rk)  :: a,x
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

   end module fabm_niva_brom_bio

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

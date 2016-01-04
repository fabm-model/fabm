#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:

   module fabm_niva_brom_redox
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
   type,extends(type_base_model),public :: type_niva_brom_redox
!     Variable identifiers
    type (type_state_variable_id)        :: id_Mn2,id_Mn3,id_Mn4,id_MnCO3
    type (type_state_variable_id)        :: id_H2S
    type (type_state_variable_id)        :: id_Fe2,id_Fe3,id_FeS,id_MnS,id_NO2
    type (type_state_variable_id)        :: id_S0,id_S2O3,id_SO4
    type (type_state_variable_id)        :: id_Baae,id_Bhae,id_Baan,id_Bhan
    type (type_state_variable_id)        :: id_O2,id_NO3,id_NH4
    type (type_state_variable_id)        :: id_PON,id_DON,id_PO4,id_Si,id_Sipart
    type (type_state_variable_id)        :: id_DIC,id_Alk,id_CaCO3,id_FeS2

    type (type_dependency_id)            :: id_temp,id_salt
    type (type_dependency_id)            :: id_Kp1,id_Kp2,id_Kp3,id_Knh4,id_Kh2s1,id_Hplus,id_KSi
    type (type_dependency_id)            :: id_Ca,id_CO3,id_Om_Ca,id_Om_Ar !,!id_CaCO3_diss      
    !type (type_dependency_id)            :: id_par,id_temp,id_h
    !type (type_horizontal_dependency_id) :: id_I_0

    type (type_diagnostic_variable_id)  :: id_DcPM_O2,id_DcDM_O2,id_DcPM_NOX, id_DcDM_NOX,id_DcPM_SO4, id_DcDM_SO4
    type (type_diagnostic_variable_id)  ::  id_DcDM_Fe,id_DcPM_Fe,id_DcDM_Mn,id_DcPM_Mn
    type (type_diagnostic_variable_id)  :: id_ChemBaae, id_ChemBaan,  id_HetBhan,id_HetBhae,id_MortBaan,id_MortBhan,id_MortBaae,id_MortBhae
    type (type_diagnostic_variable_id)  :: id_fe_ox,id_autolis,id_anammox,id_fe_ox2 !id_mns_form
    type (type_diagnostic_variable_id)  ::  id_s23_rd_DM,id_s23_rd_PM,id_s23_no3,id_s0_no3,id_s4_rd_PM,id_s4_rd_DM,id_s0_ox,id_s23_ox
    type (type_diagnostic_variable_id)  :: id_Nitrif1, id_Nitrif2, id_fe_rd,id_Denitr1_PM,id_Disprop,id_hs_ox,id_sulfido,id_s4_rd,id_s23_rd
    type (type_diagnostic_variable_id)  :: id_Denitr1_DM,id_Denitr2_PM,id_Denitr2_DM,id_Denitr1, id_Denitr2
    type (type_diagnostic_variable_id)  :: id_mn_ox2,id_mn_ox,id_mn_rd,id_mn_rd2,id_mns_diss,id_mns_prec,id_mnco3_diss,id_mnco3_prec

!     Model parameters
      real(rk) :: Wsed= 5. !1Rate of sinking of detritus (POP, PON)d-1 !!  Wdetr=1.5 (Savchuk, Wulff,1996),!Wdetr= 3.5; 20. (Gregoire,2000)
      real(rk) :: Wbact=0.4 !Rate of sinking of bacteria (Bhae,Baae,Bhan,Baan) d-1
      real(rk) :: Wm= 7. !7. !15. !Rate of accelerated sinking of particles with settled Mn hydroxides d-1
      !===========================================================================!
      !-------------------------------------------------------------------------
      ! specific rates of biogeochemical processes
       !---- Mn---------!
      real(rk) :: K_mn_ox=0.1! 10. !0.01 !0.3 !0.75 !1.5 !10.; % Specific rate of oxidation of Mn2 to Mn3  with O2 (1/day).
      real(rk) :: K_mn_ox2=0.2 !10.! 0.01 !0.2 !-2 !10; %  Specific rate of oxidation of Mn3 to Mn4  with O2 (1/day)
      real(rk) :: K_mn_rd=0.5 !1. !0.5 !0.5 !0.4 ! 5;  %  Specific rate of reduction of Mn4 to Mn3  with H2S (1/day)
      real(rk) :: K_mn_rd2=1. !0.25 !0.25 !0.2 ! -2 !10; %  Specific rate of reduction of Mn3 to Mn2  with H2S (1/day)
      real(rk) :: K_mns= 0.02 !2000. ! %  Conditional equilibrium constant for MnS from Mn2 with H2S (M)
      real(rk) :: K_mns_diss=0.0001 !0.0001   !10; %  Specific rate of dissolution of MnS to Mn2 and H2S (1/day)   
      real(rk) :: K_mns_form=0.0002 !0.00002 ! %  Specific rate of formation of MnS from Mn2 with H2S (1/day)            
      real(rk) :: K_mnco3= 15. !10. !2.e-2 ! Conditional equilibrium constant %  1.8e-11 (M) (Internet)      1 uM2 for Mn2+CO3->MnCO3 (Meysman,2003)
      real(rk) :: K_mnco3_diss= 7.e-4 !1.e-6 ! Specific rate of dissolution of MnCO3 (1/day)=6.8e-4   !2.5 X 10-1 yr-1(vanCap06) ! 1x10-4 yr-1 (Hunter et al, 98) !10^(-2)-10^(3) 1/yr( Van Cappellen-Wang-96 )
      real(rk) :: K_mnco3_form= 3.e-4 !2.e-8 ! Specific rate of formation of MnCO3 (1/day)=2.7e-7  !1. X 10-4 yr-1(vanCap06) ! 1x10-4 yr-1 (Hunter et al, 98)
      real(rk) :: K_mnco3_ox=0.0027  ! Specific rate of oxidation of MnCO3 with O2 (1/day)=0.0027  ( 1x10^(-6) M/yr ( Wang-Van Cappellen-96).
      real(rk) :: K_DON_Mn=0.001 !0.0001 !-0.0006   ! %  Specific rate of oxidation of DON with Mn4 (1/day)
      real(rk) :: K_PON_Mn=0.001 !0.0001 !-0.0003 ! %  Specific rate of oxidation of PON with Mn4 (1/day)
      real(rk) :: s_mnox_mn2=0.01 !threshold of Mn2 oxidation (uM Mn) (Yakushev,2007)
      real(rk) :: s_mnox_mn3=0.01 !threshold of Mn3 oxidation (uM Mn) (Yakushev,2007)
      real(rk) :: s_mnrd_mn4=0.01 !threshold of Mn4 reduciton (uM Mn) (Yakushev,2007)
      real(rk) :: s_mnrd_mn3=0.01 !threshold of Mn3 reduciton (uM Mn) (Yakushev,2007)
      !---- Fe---------!
      real(rk) :: K_fe_ox=0.5 !1.      ! Specific rate of oxidation of Fe2 to Fe3  with O2 (1/day)  *=4. (Konovalov,05)
      real(rk) :: K_fe_ox2=0.001 !0.1    ! Specific rate of oxidation of Fe2 to Fe3  with MnO2 (1/day) *=0.74 (Konovalov,05); 3x10^6 1/(M yr) is estimated in Van Cappellen-Wang-96 
      real(rk) :: K_fe_rd=0.5 !0.5 !-2    ! Specific rate of reduction of Fe3 to Fe2  with H2S (1/day) *=0.05 (Konovalov,05)
      real(rk) :: K_FeS=2510.    ! FeS equilibrium constant (Solubility Product Constant) (uM)=2510  ( 2.51x10-6 mol cm-3, Bektursuniva,11)  
      real(rk) :: K_FeS_form=5.e-4 ! Specific rate of precipitation of FeS from Fe2 with H2S (1/day)=1.e-5 (4x10-3 1/yr, Bektursunova,11)
      real(rk) :: K_FeS_diss=1.e-6 ! Specific rate of dissollution of FeS to Fe2 and H2S  (1/day)=3.e-6 (1x10-3 1/yr, Bektursunova,11)
      real(rk) :: K_FeS_ox=0.001 ! Specific rate of oxidation of FeS with O2 (1/day)=0.001(3x10^5 1/(M yr),Van Cappellen Wang,96) 
      real(rk) :: K_DON_Fe=0.00005 !-0.0003 ! %  Specific rate of oxidation of DON with Fe3 (1/day)
      real(rk) :: K_PON_Fe=0.00001 !-0.0001 ! %  Specific rate of oxidation of PON with Fe3 (1/day)
      real(rk) :: K_FeS2_form=0.000001 ! specific rate of FeS2 formation by FeS oxidation by H2S (1/day)=0.000009 (10^(-4) L/mol/s (Rickard-97)
      real(rk) :: K_FeS2_ox=0.00044 !specific rate of pyrite oxidation by O2  (1/uM/d)=4.38x10^(-4) 1/micromolar/day (Wijsman et al -2002). 
      real(rk) :: s_feox_fe2=0.001 !threshold of Fe2 reduciton
      real(rk) :: s_ferd_fe3=0.01  !threshold of Fe3 reduciton  (uM Fe)
      !---- S ---------!
      real(rk) :: K_hs_ox=0.5    ! Specific rate of oxidation of H2S to S0  with O2 (1/day)
      real(rk) :: K_s0_ox=0.02 !-0.2    ! Specific rate of oxidation of S0 with O2
      real(rk) :: K_s23_ox=0.01   ! Specific rate of oxidation of S2O3 with O2
      real(rk) :: K_s4_rd=0.000005 !-0.00001!  Specific rate of OM sulfate reduction with sulfate
      real(rk) :: K_s23_rd=0.001 !-0.001 ! Specific rate of OM sulfate reduction with thiosulfate
      real(rk) :: K_dispro=0.001 !-0.01  ! Specific rate of S0 dispropotionation
      real(rk) :: K_s0_no3=0.9    ! Specific rate of oxidation of S0 with NO3
      real(rk) :: K_s23_no3=0.01   ! Specific rate of oxidation of S2O3 with NO3
      real(rk) :: k_mnrdHS=1.     !half sat. of Mn reduction (uM S)
      real(rk) :: k_ferdHS =1.     !half sat. of Fe reduction (uM S)      
       !---- N---------!
      real(rk) ::  K_DON_ox=0.01 !4 ! %  Specific rate of oxidation of DON with O2 (1/day) = 0.002(S,W,96)0.1-1(W,K,91)
      real(rk) ::  K_PON_ox=0.002 !2 ! %  Specific rate of oxidation of PON with O2 (1/day) =0.002 (S,W,96)
      real(rk) ::  Tda=13.       ! Temperature control coefficient for OM decay
      real(rk) ::  beta_da=20.   ! Temperature control coefficient for OM decay
      real(rk) ::  K_omox_o2=1. !  % !half sat. of o2 for OM mineralization (uM) 
      real(rk) ::  K_PON_DON=0.1 ! %  Specific rate of autolis of PON to DON (1/day)
      real(rk) ::  KN42=0.01   !    %! Spec.rate of 1st st. of nitrification =0.01(Sawchuk,96)0.1(Gregoire,01)
      real(rk) ::  KN23=0.1   !      %! Spec.rate of 2d st. of nitrification
      real(rk) ::  KN32=0.20  !    %! Spec.rate of 1 stage of denitrif =0.16(Y,98),0.5(S&W,96),0.015(Gregoire,01)
      real(rk) ::  KN24=0.25  !    %! Spec.rate of 2 stage of denitrif =0.22 (Y,98)
      real(rk) ::  k_omno_no3=0.001 !1 !  %!half sat. of no3 for OM denitr. (uM N)
      real(rk) ::  k_omno_no2=0.001 !2 ! %!half sat. of no2 for OM denitr. (uM N)
      real(rk) ::  KT=0.8          !   %! Spec.rate of thiodenitrification
      real(rk) ::  k_annamox=0.8 !0.8   ! Spec.rate of anammox
       !---- O2--------! 
      real(rk) ::  O2s_nf=4.488  ! half saturation for nitrification
      real(rk) ::  O2s_dn=10.      ! half saturation for denitrification
      real(rk) ::  s_omox_o2=0.01 !threshold of o2 for OM mineralization
      real(rk) ::  s_omno_o2=25.  !threshold of o2 for OM denitrification
      real(rk) ::  s_omso_o2=25.  !threshold of o2 for OM sulfate reduction
      real(rk) ::  s_omso_no=5.   !threshold of noX for OM sulfate reduction
      real(rk) ::  k_mnoxO2=2.    !half sat. of Mn oxidation (uM O) (Yakushev,2007)
       !---- C--------! 
      real(rk) ::  k_CaCO3_diss = 3.    !  CaCO3 dissollution rate constant (wide ranges are given in (Luff et al., 2001))
      real(rk) ::  k_CaCO3_prec = 0.0002 ! 0.0001 !  CaCO3 precipitation rate constant (wide ranges are given in (Luff et al., 2001))
       !---- Si-------! 
      real(rk) ::  k_Sipart_diss = 0.10    !  Si dissollution rate constant (1/day)=0.01 (Popova,11)
       !---- Bacteria-! 
      real(rk) ::  k_Baae_gro = 0.019    !  Baae maximum specific growth rate (1/day)
      real(rk) ::  k_Baae_mrt = 0.005    !  Baae specific rate of mortality (1/day)
      real(rk) ::  k_Baae_mrt_h2s = 0.899    !  Baae increased specific rate of mortality due to H2S (1/day)
      real(rk) ::  limBaae=2.   ! Limiting parameter for nutrient consumprion by Baae 

      real(rk) ::  k_Bhae_gro = 0.5    !  Bhae maximum specific growth rate (1/day)
      real(rk) ::  k_Bhae_mrt = 0.02    !  Bhae specific rate of mortality (1/day)
      real(rk) ::  k_Bhae_mrt_h2s = 0.799    !  Bhae increased specific rate of mortality due to H2S (1/day)      
      real(rk) ::  limBhae=5.   ! Limiting parameter for OM consumprion by Bhae 

      real(rk) ::  k_Baan_gro = 0.012 !0.060 !B= 0.012 !0.011 !0.017    !  Baan maximum specific growth rate (1/day)
      real(rk) ::  k_Baan_mrt = 0.012 !0.012 !B= 0.012    !  Baan specific rate of mortality (1/day)
      
      real(rk) ::  limBaan=2    ! Limiting parameter for nutrient consumprion by Baan 

      real(rk) ::  k_Bhan_gro = 0.1    !  Bhan maximum specific growth rate (1/day)
      real(rk) ::  k_Bhan_mrt = 0.007    !  Bhan specific rate of mortality (1/day)
      real(rk) ::  k_Bhan_mrt_o2 = 0.899    !  Bhan increased specific rate of mortality due to O2 (1/day)
      real(rk) ::  limBhan=2.   ! Limiting parameter for OM consumprion by Bhan 

      !===========================================================================!
      !---- Stoichiometric coefficients ----!
      real(rk) ::   Sp=0.001     !P[uM]/BIOMASS [mg/m3]
      real(rk) ::   SpZ=0.001    !P[uM]/BIOMASS ZOO [mg/m3]
      real(rk) ::   Sn=0.016     !N[uM]/BIOMASS [mg/m3]
      real(rk) ::   Ss=0.053     !S[uM]/BIOMASS [mg/m3] of bacteria, i.e. during ChemBaaeynthesis
      real(rk) ::   Sc=0.106     !C[uM]/BIOMASS [mg/m3] for bact should be *.5 due to large P content
      real(rk) ::   OkP=106.     !O[uM]/P[uM]
      real(rk) ::   NkP=16.      !N[uM]/P[uM]
      real(rk) ::   OkN=6.625    !O[uM]/N[uM]
      real(rk) ::   SkP=53.      !S[uM]/P[uM] 
      real(rk) ::   CkN=8.       !C[uM]/N[uM] . PROVER' !!!!!!
      real(rk) ::   SikN=2.      !Si[uM]/N[uM] (Ivanenkov, 1979)
      real(rk) ::   FekN=26.5    !Fe[uM]/N[uM] (Boudreau, 1996) 424/16
      real(rk) ::   MnkN=13.25   !Mn[uM]/N[uM] (Boudreau, 1996) 212/16
      real(rk) ::   f=0.66       ! conversion factor relating solid and dissolved species concentrations [-]

   contains
      procedure :: initialize
      procedure :: do
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
   class (type_niva_brom_redox), intent(inout), target :: self
   integer,                      intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   !call self%get_parameter(self%w_bio,'w_bio','m/s','vertical velocity of phytoplankton (<0 for sinking)',default=-1.5_rk)

   ! Register state variables
   call self%register_state_variable(self%id_Mn2, 'Mn2', 'mmol/m**3','Mn(II)', minimum=0.0_rk)
   call self%register_state_variable(self%id_Mn3, 'Mn3', 'mmol/m**3','Mn(III)',minimum=0.0_rk)
!   call self%register_state_variable(self%id_Mn4, 'Mn4', 'mmol/m**3','manganese IV', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   call self%register_state_variable(self%id_Mn4, 'Mn4', 'mmol/m**3','Mn(IV)', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)   
   call self%register_state_variable(self%id_H2S, 'H2S', 'mmol/m**3','H2S', minimum=0.0_rk)
   call self%register_state_variable(self%id_MnS, 'MnS', 'mmol/m**3','MnS', minimum=0.0_rk)
   call self%register_state_variable(self%id_MnCO3, 'MnCO3', 'mmol/m**3','MnCO3', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)

   call self%register_state_variable(self%id_Fe2, 'Fe2', 'mmol/m**3','Fe(II)', minimum=0.0_rk)
   call self%register_state_variable(self%id_Fe3, 'Fe3', 'mmol/m**3','Fe(III)', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   call self%register_state_variable(self%id_FeS, 'FeS', 'mmol/m**3','FeS', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   !!   call self%register_state_variable(self%id_FeS, 'FeS', 'mmol/m**3','iron sulfide', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)

   call self%register_state_variable(self%id_NO2, 'NO2', 'mmol/m**3','NO2', minimum=0.0_rk)
   call self%register_state_variable(self%id_S0 , 'S0',  'mmol/m**3','S0', minimum=0.0_rk)
   call self%register_state_variable(self%id_S2O3,'S2O3','mmol/m**3','S2O3', minimum=0.0_rk)
   call self%register_state_variable(self%id_SO4, 'SO4', 'mmol/m**3','SO4', minimum=0.0_rk)
   call self%register_state_variable(self%id_Si, 'Si', 'mmol/m**3','Si', minimum=0.0_rk)
   call self%register_state_variable(self%id_Sipart, 'Sipart', 'mmol/m**3','Si Particulate', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)

   call self%register_state_variable(self%id_Baae, 'Baae', 'mmol/m**3','Aerobic Autotrophic Bacteria', minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
   call self%register_state_variable(self%id_Bhae, 'Bhae', 'mmol/m**3','Aerobic Heterotrophic Bacteria', minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
   call self%register_state_variable(self%id_Baan, 'Baan', 'mmol/m**3','Anaerobic Autotrophic Bacteria', minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)
   call self%register_state_variable(self%id_Bhan, 'Bhan', 'mmol/m**3','Anaerobic Heterotrophic Bacteria', minimum=0.0_rk,vertical_movement=-self%Wbact/86400._rk)

   call self%register_state_variable(self%id_CaCO3, 'CaCO3', 'mmol/m**3','CaCO3', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_FeS2, 'FeS2', 'mmol/m**3','FeS2', minimum=0.0_rk, vertical_movement=-self%Wsed/86400._rk)
   ! Register the contribution of all state variables to total nitrogen
   !call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_bio)

   call self%register_dependency(self%id_Hplus, 'Hplus', 'mmol/m**3','H+')
   call self%register_dependency(self%id_Om_Ca,'Om_Ca','-','Omega CaCO3-Calcite')
   call self%register_dependency(self%id_Om_Ar,'Om_Ar','-','Omega CaCO3-Aragonite')
   call self%register_dependency(self%id_CO3,'CO3','mmol/m**3','CO3--')
   call self%register_dependency(self%id_Ca,'Ca','mmol/m**3','Ca++')
  !!   call self%register_dependency(self%id_Kh2s1,'Kh2s1','-', 'H2S <--> H+ + HS-')

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   call self%register_state_dependency(self%id_DIC,'DIC','mmol/m**3','total dissolved inorganic carbon',required=.false.)
   call self%register_state_dependency(self%id_Alk,'Alk','mmol/m**3','total alkalinity',required=.false.)
   call self%register_state_dependency(self%id_po4,'PO4','mmol/m**3','phosphate',required=.false.)

   call self%register_state_dependency(self%id_O2, 'O2', 'mmol/m**3','dissolved oxygen')
   call self%register_state_dependency(self%id_NH4,'NH4','mmol/m**3','ammonium')
   call self%register_state_dependency(self%id_NO3,'NO3','mmol/m**3','nitrate')
 !   call self%register_state_dependency(self%id_NO2,'NO2','mmol/m**3','nitrite')

   call self%register_state_dependency(self%id_PON,'PON','mmol/m**3','particulate organic nitrogen')
   call self%register_state_dependency(self%id_DON,'DON','mmol/m**3','dissolved organic nitrogen')

       !!!! Register diagnostic variables

   !call self%register_diagnostic_variable(self%id_GPP,'GPP','mmol/m**3',  'gross primary production',           &
   !                  output=output_time_step_integrated)

    call self%register_diagnostic_variable(self%id_DcPM_O2,'DcPM_O2','mmol/m**3',  'POM with O2 mineralization',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_DcDM_O2,'DcDM_O2','mmol/m**3',  'DOM with O2 mineralization',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_DcPM_NOX,'DcPM_NOX','mmol/m**3',  'POM denitrification (1+2 stage)',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_DcDM_NOX,'DcDM_NOX','mmol/m**3',  'DOM denitrification (1+2 stage)',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_DcPM_SO4,'DcPM_SO4','mmol/m**3',  'POM sulfatereduction (1+2 stage)',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_DcDM_SO4,'DcDM_SO4','mmol/m**3',  'DOM sulfatereduction (1+2 stage)',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_DcPM_Mn,'DcPM_Mn','mmol/m**3',  'POM with Mn(IV) mineralization ',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_DcDM_Mn,'DcDM_Mn','mmol/m**3',  'DOM with Mn(IV) mineralization',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_DcPM_Fe,'DcPM_Fe','mmol/m**3',  'POM with Fe(III) mineralization',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_DcDM_Fe,'DcDM_Fe','mmol/m**3',  'DOM with Fe(III) mineralization  ',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_ChemBaae,'ChemBaae','mmol/m**3',  'Growth of  Aerobic Autotrophic Bacteria',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_HetBhae,'HetBhae','mmol/m**3',  'Growth of  Aerobic Heterotrophic Bacteria',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_autolis,'autolis','mmol/m**3',  'Autolis',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_mn_ox,'mn_ox','mmol/m**3',  'Mn(II) with O2 oxidation ',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_mn_ox2,'mn_ox2','mmol/m**3',  'Mn(III) with O2 oxidation',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_mns_diss,'mns_diss','mmol/m**3',  'MnS dissolution',           &
                output=output_time_step_integrated)   
    call self%register_diagnostic_variable(self%id_mns_prec,'mns_prec','mmol/m**3',  'MnS formation',           &
                output=output_time_step_integrated)    
    call self%register_diagnostic_variable(self%id_mnco3_diss,'mnco3_diss','mmol/m**3',  'MnCO3 dissolusion',           &
                output=output_time_step_integrated)   
    call self%register_diagnostic_variable(self%id_mnco3_prec,'mnco3_prec','mmol/m**3',  'MnCO3 formation ',           &
                output=output_time_step_integrated)    
    call self%register_diagnostic_variable(self%id_mn_rd,'mn_rd','mmol/m**3',  'Mn(IV) with H2S reduction',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_mn_rd2,'mn_rd2','mmol/m**3',  'Mn(III) with H2S reduction',           &
                output=output_time_step_integrated)    
    call self%register_diagnostic_variable(self%id_HetBhan,'HetBhan','mmol/m**3',  'Growth of Anaerobic Heterotrophic Bacteria ',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_ChemBaan,'ChemBaan','mmol/m**3',  'Growth of Anaerobic Autotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_MortBhan,'MortBhan','mmol/m**3',  'Mortality of Anaerobic Heterotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_MortBaan,'MortBaan','mmol/m**3',  'Mortality of Anaerobic Autotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_MortBaae,'MortBaae','mmol/m**3',  'Mortality of Aerobic Autotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_MortBhae,'MortBhae','mmol/m**3',  'Mortality of Aerobic Heterotrophic Bacteria ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Nitrif1,'Nitrif1','mmol/m**3',  'Nitrification 1 stage',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_Nitrif2,'Nitrif2','mmol/m**3',  'Nitrification 2 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_anammox ,'anammox ','mmol/m**3',  'Anammox',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr1_PM ,'Denitr1_PM ','mmol/m**3',  'POM denitrification 1 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr1_DM ,'Denitr1_DM ','mmol/m**3',  'DOM denitrification 1 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr2_PM ,'Denitr2_PM ','mmol/m**3',  'POM denitrification 2 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr2_DM ,'Denitr2_DM ','mmol/m**3',  'DOM denitrification 2 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr1 ,'Denitr1 ','mmol/m**3',  '(POM+DOM) denitrification 1 stage',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Denitr2 ,'Denitr2 ','mmol/m**3',  '(POM+DOM) denitrification 2 stage',           &
                output=output_time_step_integrated)     
    call self%register_diagnostic_variable(self% id_fe_ox,'fe_ox','mmol/m**3',  'Fe(II) with O2 oxidation ',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self% id_fe_ox2,'fe_ox2','mmol/m**3',  'Fe(II)  with Mn(IV) oxidation',           &
                output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_fe_rd,'fe_rd','mmol/m**3',  'Fe (III) with H2S reduction',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_s23_rd_PM,'s23_rd_PM','mmol/m**3',  'POM sulfatereduction 2d stage',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_s23_rd_DM,'s23_rd_DM','mmol/m**3',  'DOM sulfatereduction 2d stage',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_s23_no3,'s23_no3','mmol/m**3',  ' S2O3 with NO3 oxidation',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_s0_no3,'s0_no3','mmol/m**3',  'S0 with NO3 oxidation',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_s23_ox,'s23_ox','mmol/m**3',  'S2O3  with O2oxidation',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_s0_ox,'s0_ox','mmol/m**3',  'S0 with O2 oxidation',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_s4_rd,'s4_rd','mmol/m**3',  '(POM+DOM) sulfatereduction 1st stage',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_s23_rd,'s23_rd','mmol/m**3',  '(POM+DOM) sulfatereduction 2d stage ',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_Disprop,'Disprop','mmol/m**3',  'S0 disproportionation',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_hs_ox,'hs_ox','mmol/m**3',  'H2S with O2 oxidation',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_sulfido,'sulfido','mmol/m**3',  'H2S with NO3 oxidation',           &
                output=output_time_step_integrated) 
    call self%register_diagnostic_variable(self%id_s4_rd_PM,'s4_rd_PM','mmol/m**3',  'POM sulfatereduction 1st stage',           &
                output=output_time_step_integrated)  
    call self%register_diagnostic_variable(self%id_s4_rd_DM,'s4_rd_DM','mmol/m**3',  'DOM sulfatereduction 1st stage',           &
                output=output_time_step_integrated) 

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)

   ! Register dependencies on equilibrium constants for phosphoric acid, ammonia, hydrogen sulfide.
   call self%register_dependency(self%id_Kp1,  'Kp1',  '-', '[H+][H2PO4-]/[H3PO4]')
   call self%register_dependency(self%id_Kp2,  'Kp2',  '-', '[H][HPO4]/[H2PO4]')
   call self%register_dependency(self%id_Kp3,  'Kp3',  '-', '[H][PO4]/[HPO4]')
   call self%register_dependency(self%id_Knh4, 'Knh4', '-', '[H+][NH3]/[NH4]')
   call self%register_dependency(self%id_Kh2s1,'Kh2s1','-', '[H+][HS-]/[H2S]')
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
   class (type_niva_brom_redox),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk)                   :: O2,PON,DON,Si,Sipart
   real(rk)                   :: Mn2,Mn3,Mn4,H2S,MnS,MnCO3
   real(rk)                   :: Fe2,Fe3,FeS,FeS2
   real(rk)                   :: NO2,NO3,NH4, PO4
   real(rk)                   :: S0,S2O3,SO4
   real(rk)                   :: Baae,Bhae,Baan,Bhan
   real(rk)                   :: DIC,Alk,Hplus,CaCO3,Ca,CO3,Om_Ca,Om_Ar
   real(rk)                   :: temp,salt

   real(rk) :: autolis,DcDM_O2,DcPM_O2, DcPM_NOX,DcDM_NOX
   real(rk) :: mn_ox,mn_ox2,mn_rd,mn_rd2,Om_MnS,mns_prec,mns_diss,DcDM_Mn,DcPM_Mn,Om_MnCO3,mnco3_prec,mnco3_diss,mnco3_ox
   real(rk) :: fe_ox,fe_rd,fe_ox2,Om_FeS,fes_prec,fes_diss,fes_ox,DcDM_Fe,DcPM_Fe,fes2_form,fes2_ox
   real(rk) :: Nitrif1,Nitrif2,anammox,Denitr1_PM,Denitr1_DM,Denitr2_PM,Denitr2_DM,Denitr1,Denitr2
   real(rk) :: Disprop,hs_ox,s0_ox,s0_no3,s23_ox,s23_no3,sulfido,s4_rd_PM,s4_rd_DM,s23_rd_PM,s23_rd_DM,s4_rd,s23_rd,DcPM_SO4,DcDM_SO4
   real(rk) :: ChemBaae,MortBaae,MortBhae,ChemBaan,MortBaan,MortBhan,HetBhan,HetBhae
   real(rk) :: Knh4,Kp1,Kp2,Kp3,Kh2s1,KSi
   real(rk) :: dAlk, Dc_OM_total,CaCO3_prec,CaCO3_diss

   
   integer  iter, oldcolor
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_O2,O2)

   _GET_(self%id_PON,PON)
   _GET_(self%id_DON,DON)

   _GET_(self%id_Mn2,Mn2)
   _GET_(self%id_Mn3,Mn3)
   _GET_(self%id_Mn4,Mn4)
   _GET_(self%id_MnS,MnS)
   _GET_(self%id_MnCO3,MnCO3)   
   _GET_(self%id_H2S,H2S)

   _GET_(self%id_Fe2,Fe2)
   _GET_(self%id_Fe3,Fe3)
   _GET_(self%id_FeS,FeS)
   _GET_(self%id_FeS2,FeS2)

   _GET_(self%id_PO4,PO4)

   _GET_(self%id_NH4,NH4)
   _GET_(self%id_NO2,NO2)
   _GET_(self%id_NO3,NO3)

   _GET_(self%id_S0,  S0)
   _GET_(self%id_S2O3,S2O3)
   _GET_(self%id_SO4, SO4)

   _GET_(self%id_Baae,Baae)
   _GET_(self%id_Bhae,Bhae)
   _GET_(self%id_Baan,Baan)
   _GET_(self%id_Bhan,Bhan)

   _GET_(self%id_DIC,DIC)
   _GET_(self%id_Alk,Alk)
   _GET_(self%id_Hplus,Hplus)
   _GET_(self%id_CaCO3,CaCO3)  
   _GET_(self%id_CO3,CO3)  
   _GET_(self%id_Om_Ca,Om_Ca)
   _GET_(self%id_Om_Ar,Om_Ar)   
   _GET_(self%id_Ca,Ca)  

   ! Get equilibrium constants
   _GET_(self%id_Kp1,  Kp1)
   _GET_(self%id_Kp2,  Kp2)
   _GET_(self%id_Kp3,  Kp3)
   _GET_(self%id_Kh2s1,Kh2s1)
   _GET_(self%id_Knh4, Knh4)
   _GET_(self%id_KSi, KSi)
   _GET_(self%id_Si,Si)
   _GET_(self%id_Sipart,Sipart)   

   ! Retrieve current environmental conditions.
   !_GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_salt,salt)              ! salinity   

!%%
!%!-------------------------------------------------------------------------
!%!========PON & DON========================================================
!%!-------------------------------------------------------------------------   
     autolis=self%K_PON_DON*PON
!%(CH2O)106(NH3)16H3PO4 + 106 O2 -> 106 CO2 + 106 H2O + 16 NH3 + H3PO4 : 
     DcDM_O2=self%K_DON_ox*DON*O2/(O2+self%K_omox_o2) &
     *(1.+self%beta_da*yy(self%tda,temp))
     !+++! *exp(0.15*t) !(Savchuk, Wulff,1996)       
!%(CH2O)106(NH3)16H3PO4 + 106 O2 -> 106 CO2 + 106 H2O + 16 NH3 + H3PO4 :
     DcPM_O2=self%K_PON_ox*PON*O2/(O2+self%K_omox_o2) &
     *(1.+self%beta_da*yy(self%tda,temp))
!%!-------------------------------------------------------------------------
!%!========Mn===============================================================
!%!-------------------------------------------------------------------------
!% Mn2 oxidation: 4Mn2+ + O2 + 4H+ -> 4Mn3+ + 2H2O (Canfield, 2005):
     mn_ox =max(0._rk, 0.5*(1.+tanh(Mn2-self%s_mnox_mn2))*self%K_mn_ox*Mn2*o2/(o2+self%k_mnoxO2) )
!% Mn3 oxidation: 2Mn3+ + 0.5O2 + 3H20 -> 2MnO2 + 6H+ (Tebo, 1997)	:
     mn_ox2=max(0._rk, 0.5*(1.+tanh(Mn3-self%s_mnox_mn3))*self%K_mn_ox2*Mn3*o2/(o2+self%k_mnoxO2) )
!% Mn4 reduction: 2MnO2 + 7H+ + HS- -> 2Mn3+ + 4H2O + S0 :  
     mn_rd =max(0._rk, 0.5*(1.+tanh(Mn4-self%s_mnrd_mn4))*self%K_mn_rd*Mn4*h2s/(h2s+self%k_mnrdHS) )
!% Mn3 reduction: 2Mn3+ + HS- -> 2Mn2+ + S0 + H+ :        
     mn_rd2=max(0._rk, 0.5*(1.+tanh(Mn3-self%s_mnrd_mn3))*self%K_mn_rd2*Mn3*h2s/(h2s+self%k_mnrdHS) )
!!!!!!!% solid MnS formation Mn2+ + H2S -> MnS + 2H+ :
!!!!!!!     mns_prec= self%K_mns_form*H2S*Mn2
!!!!!!!% MnS dissolution: MnS + 2O2 -> Mn2+ + SO42- :
!!!!!!!     mns_ox= self%K_mns_ox*o2*MnS

!MnS formation/dissollution (dSED) Mn2+  +  HS-  =  MnS(s)  +  H+
          Om_MnS=H2S*1.e-6*Mn2*1.e-6/(self%K_MnS*Hplus) 
!% Mn2+  +  HS-  ->  MnS(s)  +  H+ 
         mns_prec=self%K_mns_form*max(0._rk,(Om_MnS-1._rk))
!%   MnS(s)  +  H+ -> Mn2+  +  HS-  
         mns_diss=self%K_MnS_diss*MnS*max(0._rk,(1._rk-Om_MnS))  
!% MnS + 2O2 ->   Mn2+ + SO42-  ?

!  MnCO3 precipitation/dissolution
 !   Om_MnCO3=Mn2*HCO3/(self%K_MnCO3*CO2)
        Om_MnCO3=Mn2*CO3/(self%K_MnCO3)
!% ! Mn2+ + CO3-- <-> MnCO3 (vanCappelen,96): 
        mnco3_prec=self%K_mnco3_form*max(0._rk,(Om_MnCO3-1._rk))
!% 
        mnco3_diss=self%K_MnCO3_diss*MnCO3*max(0._rk,(1._rk-Om_MnCO3))
 !% 2 MnCO3(s)  +  O2  +  2 H2O  =  2 MnO2(s)  +  2 HCO3-  +  2 H+ (Morgan,05)       
        mnco3_ox=self%K_mnco3_ox*MnCO3*O2

!% (CH2O)106(NH3)16(H3PO4) +212MnO2 + 318CO2 +106H2O ->  424HCO3- + 212Mn2+ +16NH3 +H3PO4  (Boudreau, 1996)   ! in N units
     DcDM_Mn=max(0._rk,self%K_DON_Mn*DON &
                        *Mn4/(Mn4 +0.5) &
                        *(1.-0.5*(1+tanh(o2-self%O2s_dn)))  )
!% (CH2O)106(NH3)16(H3PO4) +212MnO2 + 318CO2 +106H2O ->  424HCO3- + 212Mn2+ +16NH3 +H3PO4  (Boudreau, 1996)
     DcPM_Mn=max(0._rk,self%K_PON_Mn*PON &
                        *Mn4/(Mn4 +0.5) &
                        *(1.-0.5*(1+tanh(o2-self%O2s_dn)))   )    

!%!-------------------------------------------------------------------------
!%!========Fe===============================================================
!%!-------------------------------------------------------------------------
!% Fe2 oxidation1:4Fe2+ + O2 +10H2O-> 4Fe(OH)3 +8H+    (vanCappelen,96):
    fe_ox= 0.5*(1.+tanh(Fe2-self%s_feox_fe2))*self%K_fe_ox*o2*Fe2   
!% Fe2 oxidation2: Fe2+ + MnO2 + 4H+ -> Fe3+ + Mn2+ + 2H2O (vanCappelen,96):
    fe_ox2=0.5*(1.+tanh(Fe2-self%s_feox_fe2))*self%K_fe_ox2*Mn4*Fe2 
!% Fe3 reduction:  2Fe(OH)3+HS-+5H+ ->  2Fe2++S0+6H2O  
    fe_rd=0.5*(1.+tanh(Fe3-self%s_ferd_fe3))*self%K_fe_rd*Fe3 *h2s/(h2s+self%k_ferdHS) !*Kh2s1/(Kh2s1+Hplus)?

!FeS formation/dissollution (Bektursunova, 11)
          Om_FeS=H2S*Fe2/(self%K_FeS*Hplus*1000000.) 
!% FeS formation Fe2+ + HS- -> FeS + H+ (Bektursunova, 11)
    fes_prec=self%K_FeS_form*max(0._rk,(Om_FeS-1._rk))
!% FeS dissollution FeS + H+ -> Fe2+ + HS (Bektursunova, 11)
    fes_diss=self%K_FeS_diss*FeS*max(0._rk,(1._rk-Om_FeS))     
!% FeS oxidation: FeS + 2.25O2 +2.5H2O -> Fe (OH)3 + 2H+ +SO42- (Soetaert,07) or FeS + 2O2 -> Fe2+ + SO42-(Bektursunova, 11) :
    fes_ox=self%K_FeS_ox*FeS*O2

! Pyrite formation by FeS oxidation by H2S 
! FeS + H2S -> FeS2 + H2     (Rickard,97)
      fes2_form=self%K_FeS2_form*H2S*FeS   
 !Pyrite oxidation by O2
 !  FeS2 + 3.5 O2 + H2O = Fe2+ + 2SO42- + 2H+  (Wijsman,02). 
      fes2_ox=self%K_FeS2_ox*FeS2*o2

!% (CH2O)106(NH3)16H3PO4 + 424Fe(OH)3 + 742CO2 -> 848HCO3-+ 424Fe2+ +318 H2O +16NH3 +H3PO4  (Boudreau, 1996)  Fe units
     DcDM_Fe=self%K_DON_Fe*DON &
                        *Fe3 /(Fe3 +self%k_omno_no3) &
                        *(1.-0.5*(1+tanh(o2-self%O2s_dn)))  
!% (CH2O)106(NH3)16H3PO4 + 424Fe(OH)3 + 742CO2 -> 848HCO3-+ 424Fe2+ +318 H2O +16NH3 +H3PO4  (Boudreau, 1996)
     DcPM_Fe=self%K_PON_Fe*PON &
                        *Fe3 /(Fe3 +self%k_omno_no3) &
                        *(1.-0.5*(1+tanh(o2-self%O2s_dn)))     
!%!-------------------------------------------------------------------------
!%!========N================================================================
!%!-------------------------------------------------------------------------
!%Nitrification 1st stage: NH4+ + 1.5 O2 -> NO2- + 2H+ + H2O (Canfield,05)
     Nitrif1 = self%KN42*NH4*o2*0.5*(1.+tanh(o2-self%O2s_nf))
!%Nitrification 2d stage: NO2- + 0.5 O2 -> NO3- (Canfield,05):
     Nitrif2 = self%KN23*NO2*o2*0.5*(1.+tanh(o2-self%O2s_nf))
!%!---------------------------------------- in suboxic conditions:
!%! Anammox : NO2- + NH4+ -> N2 + 2H2O (Canfield,05):
     anammox = self%k_annamox*NO2*NH4 &    
                      *(1.-0.5*(1+tanh(o2-self%O2s_dn)))
!%! OM denitrification (Richards, 1965): (CH2O)106(NH3)16H3PO4 + 84.8HNO3 = 106CO2 + 42.4N2 + 148.4H2O + 16NH3 + H3PO4  
!%! POM denitrification (1st stage) (Anderson,82): 1/2CH2O + NO3- -> NO2- + 1/2H2O + 1/2CO2 : 
     Denitr1_PM    = self%KN32*PON &
                      *(1.-0.5*(1+tanh(o2-self%O2s_dn)))  &
                      *NO3/(NO3 +self%k_omno_no3)
!%!---------------------------------------- 
!%! DOM denitrification (1st stage): 1/2CH2O + NO3- -> NO2- + 1/2H2O + 1/2CO2 (Anderson,82):
     Denitr1_DM    = self%KN32*DON &
                      *(1.-0.5*(1.+tanh(o2-self%O2s_dn)))  &
                      *NO3/(NO3+self%k_omno_no3)
!%!---------------------------------------- 
!%! POM denitrification (2d stage): 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 (Anderson,82):
     Denitr2_PM    = self%KN24*PON &
                      *(1.-0.5*(1.+tanh(o2-self%O2s_dn)))  &
                      *NO2/(NO2+self%k_omno_no2)
!%!---------------------------------------- 
!%! DOM denitrification (2d stage): 3/4CH2O + H+ + NO2- ->? 1/2N2 + 5/4H2O + 3/4CO2 (Anderson,82):
     Denitr2_DM    = self%KN24*DON &
                      *(1.-0.5*(1+tanh(o2-self%O2s_dn)))  &
                      *NO2/(NO2+self%k_omno_no2)
     Denitr1       =Denitr1_PM + Denitr1_DM
     Denitr2       =Denitr2_PM + Denitr2_DM
   
!% ! from Anderson  et al. and Richards: Denitrification POM (1+2 stage)
     DcPM_NOX      =16./212.*Denitr1_PM + 16./141.3*Denitr2_PM
     ! from Anderson  et al. and Richards: Denitrification DOM (1+2 stage)
     DcDM_NOX      =16./212.*Denitr1_DM + 16./141.3*Denitr2_DM
!%%
!%!-------------------------------------------------------------------------
!%!========S================================================================
!%!-------------------------------------------------------------------------
!%! S0 disproportionation: 4S0 + 3H2O -> 2H2S + S2O32= + 2H+ :
     Disprop =self%K_dispro*S0
!%!----------------------------------------
!%! HS oxidation with O2: 2H2S + O2 -> 2S0 + 2H2O :
     hs_ox=self%K_hs_ox*o2*H2S
!%!----------------------------------------
!%! S0 oxidation with O2: 2S0 + O2 + H2O -> S2O32= + 2H+ :
     s0_ox=self%K_s0_ox*o2*S0 
!%!----------------------------------------
!%! S0 oxidation with NO3: 4S0 + 3NO3- + 7H2O -> 4SO4= + 3NH4+ + 2H+  :
     s0_no3=self%K_s0_no3*NO3*S0 
!%!----------------------------------------
!%! S2O3 oxidation with O2: S2O32= + 2O2 + 2OH- -> 2SO42= + H2O :
     s23_ox=self%K_s23_ox*o2*S2O3 
!%!----------------------------------------
!%! S2O3 oxidation with NO3: S2O3= + NO3- + 2H2O --> 2SO4= + NH4+ :
     s23_no3=self%K_s23_no3*NO3*S2O3 
!%!----------------------------------------
!%! Thiodenitrification: 5H2S + 8NO3- + 2OH+ -> 5SO4 + 4N2 + 6H2O (Volkov, 1984): 
     sulfido       =self%KT*H2S*NO3
!%!---------------------------------------- in anoxic conditions:
!%! OM sulfatereduction (Boudreau, 1996)
!% (CH2O)106(NH3)16H3PO4 + 53SO42- = 106HCO3- + 16NH3 + H3PO4 + 53H2S  

!%! POM sulfatereduction (1st stage):
     s4_rd_PM      = (1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
                     *(1.-0.5*(1.+tanh(NO3-self%s_omso_no)))  &
                      *self%K_s4_rd*SO4*PON
!%! DOM sulfatereduction (1st stage):
     s4_rd_DM      =(1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
                     *(1.-0.5*(1.+tanh(NO3-self%s_omso_no)))  &
                      *self%K_s4_rd*SO4*DON

            if (o2.gt.10.) then
                        s4_rd_PM=0.  
                        s4_rd_DM=0. 
            endif

!%! POM sulfatereduction (2d stage): 
     s23_rd_PM     =(1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
                     *(1.-0.5*(1.+tanh(NO3-self%s_omso_no)))  &
                      *self%K_s23_rd*S2O3*PON
!%! DOM sulfatereduction (2d stage):
     s23_rd_DM     =(1.-0.5*(1.+tanh(o2-self%s_omso_o2)))  &
                     *(1.-0.5*(1.+tanh(NO3-self%s_omso_no)))  &
                      *self%K_s23_rd*S2O3*DON

           if (o2.gt.10.) then
                        s23_rd_PM=0.  
                        s23_rd_DM=0. 
           endif

     s4_rd         =s4_rd_PM  + s4_rd_DM  ! in S units
     s23_rd        =s23_rd_PM + s23_rd_DM

     DcPM_SO4      =16./53.*(s4_rd_PM+s23_rd_PM)  ! in N units
     DcDM_SO4      =16./53.*(s4_rd_DM+s23_rd_DM)
!%!---------------------------------------------------------------------------
!%!========B=a=c=t(ci)============================================================
!%!---------------------------------------------------------------------------

!%!---------------------------------------- OXIC CONDITIONS
!%!---------------------------------------- aerobic autotrophs
       
     ChemBaae =(Nitrif1+Nitrif2+mn_ox+fe_ox+s23_ox+s0_ox +anammox) &
              *self%k_Baae_gro*Baae*min(yy(self%limBaae,NH4/(Baae+0.0001)),yy(self%limBaae,PO4/(Baae+0.0001)))

     MortBaae=(self%k_Baae_mrt+ self%k_Baae_mrt_h2s*(0.5*(1.-tanh(1.-H2S))))*Baae*Baae    
!%!---------------------------------------- aerobic heterotroph

     HetBhae  =(DcPM_O2+DcDM_O2) &
             *self%k_Bhae_gro*Bhae*yy(self%limBhae,DON/(Bhae+0.0001))
                
     MortBhae =(self%k_Bhae_mrt+ self%k_Bhae_mrt_h2s*(0.5*(1.-tanh(1.-H2S))))*Bhae           


!%!---------------------------------------- ANOXIC CONDITIONS
!%!------------------------------------- anaerobic autotrophs
     
     ChemBaan =(mn_rd+mn_rd2+fe_rd+hs_ox +sulfido) &
              *self%k_Baan_gro*Baan*min(yy(self%limBaan,NH4/(Baan+0.0001)),yy(self%limBaan,PO4/(Baan+0.0001)))

     MortBaan = self%k_Baan_mrt*Baan*Baan

!%!---------------------------------------- anaerobic heterotroph
     HetBhan  =(DcPM_NOX+DcDM_NOX+DcPM_SO4+DcDM_SO4+DcDM_Mn+DcPM_Mn+DcPM_Fe+DcDM_Fe) &
                    *self%k_Bhan_gro*Bhan*yy(self%limBhan,DON/(Bhan+0.0001))
                
     MortBhan =(self%k_Bhan_mrt+ self%k_Bhan_mrt_o2*(0.5+0.5*(tanh(1.-O2))))*Bhan       
       
!%!---------------------------------------- Summariazed OM mineralization
     Dc_OM_total=DcDM_O2+DcPM_O2+DcPM_NOX+DcDM_NOX+DcDM_Mn+DcPM_Mn+DcDM_Fe+DcPM_Fe+DcDM_SO4+DcPM_SO4
     
!%!---------------------------------------- CaCO3 precipitation/dissolution     
      
     CaCO3_prec =self%k_CaCO3_prec*max(0._rk,(Om_Ar-1._rk))              ! (Luff et al., 2001)
     CaCO3_diss =CaCO3*self%k_CaCO3_diss*(max(0._rk,(1._rk-Om_Ar)))**4.5 ! (Luff et al., 2001)

    
!%!---------------------------------------- 
          
          _SET_ODE_(self%id_Mn2,-mn_ox+mn_rd2-mns_prec+mns_diss-mnco3_prec+mnco3_diss+0.5*fe_ox2+(DcDM_Mn+DcPM_Mn)*self%MnkN)
          _SET_ODE_(self%id_Mn3, mn_ox-mn_ox2+mn_rd-mn_rd2)
          _SET_ODE_(self%id_Mn4, mn_ox2-mn_rd-0.5*fe_ox2+mnco3_ox-(DcDM_Mn+DcPM_Mn)*self%MnkN)
          _SET_ODE_(self%id_MnS, mns_prec-mns_diss)
          _SET_ODE_(self%id_MnCO3, mnco3_prec-mnco3_diss-mnco3_ox)

          _SET_ODE_(self%id_Fe2,-fe_ox-fe_ox2+fe_rd-fes_prec+fes_diss+(DcDM_Fe+DcPM_Fe)*4.*self%FekN+feS2_ox)
          _SET_ODE_(self%id_Fe3,fe_ox+fe_ox2-fe_rd+fes_ox-(DcDM_Fe+DcPM_Fe)*4.*self%FekN)
          _SET_ODE_(self%id_FeS,fes_prec-fes_diss-fes_ox-feS2_form)
          _SET_ODE_(self%id_FeS2,feS2_form-feS2_ox)

          _SET_ODE_(self%id_H2S,-0.5*mn_rd-0.5*mn_rd2-0.5*fe_rd -hs_ox-fes_prec-mns_prec+mns_diss +0.5*Disprop-sulfido+s23_rd+fes_diss-feS2_form)
          _SET_ODE_(self%id_S0, hs_ox+0.5*mn_rd+0.5*mn_rd2+0.5*fe_rd -s0_ox -Disprop-s0_no3)
          _SET_ODE_(self%id_S2O3,0.5*s0_ox-s23_ox+0.25*Disprop +0.5*s4_rd-0.5*s23_rd-s23_no3)
          _SET_ODE_(self%id_SO4,sulfido-s4_rd+0.5*s23_ox+s0_no3+2.*s23_no3+fes_ox+2.*feS2_ox)

          _SET_ODE_(self%id_O2,(-DcDM_O2-DcPM_O2)*self%OkN -0.25*mn_ox-0.25*mn_ox2-0.25*fe_ox -0.5*hs_ox-0.5*s0_ox-0.5*s23_ox -1.5*Nitrif1-0.5*Nitrif2-2.25*fes_ox-3.5*feS2_ox-0.5*mnco3_ox )
          _SET_ODE_(self%id_DON,(autolis-DcDM_O2-DcDM_NOX-DcDM_SO4 -DcDM_Mn -DcDM_Fe -HetBhae-HetBhan))
          _SET_ODE_(self%id_PON,(-autolis-DcPM_O2-DcPM_NOX-DcPM_SO4 -DcPM_Mn -DcPM_Fe+MortBaae+MortBaan+MortBhae+MortBhan))
          _SET_ODE_(self%id_NH4,(Dc_OM_total-Nitrif1-anammox +0.75*s0_no3+s23_no3-ChemBaae-ChemBaan)) 
          _SET_ODE_(self%id_NO2,(Nitrif1-Nitrif2+Denitr1-Denitr2-anammox))
          _SET_ODE_(self%id_NO3,(Nitrif2-Denitr1-1.6*sulfido-0.75*s0_no3-s23_no3) )
  !        if (_AVAILABLE_(self%id_dic)) then
          _SET_ODE_(self%id_DIC,(Dc_OM_total-ChemBaae-ChemBaan)*self%CkN-CaCO3_prec+CaCO3_diss-mnco3_prec+mnco3_diss+mnco3_ox)
          _SET_ODE_(self%id_Si,self%k_Sipart_diss*Sipart)
          _SET_ODE_(self%id_Sipart,-self%k_Sipart_diss*Sipart)  !        end if
   !       if (_AVAILABLE_(self%id_PO4)) then
          _SET_ODE_(self%id_PO4,((Dc_OM_total-ChemBaae-ChemBaan)/self%NkP+(fe_rd-fe_ox-fe_ox2+4.*DcDM_Fe+4.*DcPM_Fe)/2.7-(mn_ox2+mn_rd2-mn_ox-mn_rd)/0.67)) 
  !        end if
          _SET_ODE_(self%id_Baae,ChemBaae-MortBaae)
          _SET_ODE_(self%id_Baan,ChemBaan-MortBaan)
          _SET_ODE_(self%id_Bhae,HetBhae-MortBhae)
          _SET_ODE_(self%id_Bhan,HetBhan-MortBhan)
          _SET_ODE_(self%id_CaCO3,CaCO3_prec-CaCO3_diss)

  !Here CaCO3=F(omega. etc.   CaCO3_prec -CaCO3_diss     
  ! CaCO3_prec = 0.5*(1.-tanh(Om_Ca-1.))*k_caco3_prec*co3*Ca
  ! CaCO3_diss = 0.5*(1.+tanh(Om_Ca-1.))*k_caco3_diss*CaCO3* (1.-Om_Ca)     !Jourabachi et al., 2005  
  ! k_caco3_diss=7.3x10-6...1.2x10-5 1/s when Ca2+ concentration varied from 345 to 850 mg/l
  ! k_caco3_prec =1.5x103 to 2.5x103 m3/(mol s) depending on the Ca2+ concentration

  !        _SET_ODE_(self%id_Si,0.) !(-MortPhy-ExcrPhy- &
  !!            -GrowthPhy+GrazPhy)/0.3 

  !! Alkalinity changes due to redox reactions:

         dAlk = (  & ! 0.01*dCc(k,i_DON) & ! part of OM charged (EPS etc.)
                    !+( Dc_OM_total-ChemBaae-ChemBaan) *Knh4/(Knh4+Hplus)  &  ! +/- NH3(tot) &
                    !
                    !+(-anammox  & !NO2- + NH4+ -> N2 + 2H2O                       ! +/- NH3(tot)
                    !  +0.75*s0_no3 & ! 4S0 + 3NO3- + 7H2O -> 4SO4-- + 3NH4+ + 2H+ ! +/- NH3(tot)
                    !  +s23_no3 &  ! S2O3-- + NO3- + 2H2O -> 2SO4-- + NH4+         ! +/- NH3(tot) 
                    ! ) *Knh4/(Knh4+Hplus) &                                       ! +/- NH3(tot) 
                    - 1.* Nitrif1 & !*(1.-Knh4/(Knh4+Hplus))  &    ! NH4+ + 1.5 O2 -> NO2- + 2H+ + H2O ! W-Gladrow,07   ! +/- HNO2t 
!                    - 2.* Nitrif2 & ! NO2- + 0.5 O2 -> NO3  !                                                          ! +/- HNO2t 
!: 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 or! 5 CH2O + 4 H+ + 4 NO3- -> 2 N2 + 5 CO2 + 7H2O ! W-Gladrow,07                    
!                    + 0.*(Denitr1_PM +Denitr1_DM)  &! &    ! +/- HNO2t 
                    + 1.*(Denitr2_PM +Denitr2_DM)  &! &    ! +/- HNO2t 
                    !+((Dc_OM_total-ChemBaae-ChemBaan)/self%NkP &                                            ! +/- PO4t 
                    !    +(fe_rd-fe_ox-fe_ox2+4.*DcDM_Fe+4.*DcPM_Fe)/2.7-(mn_ox2+mn_rd2-mn_ox-mn_rd)/0.67) & ! +/- PO4t 
                    !    *((Kp1*Kp2-Hplus*Hplus)*Hplus+2e0*Kp1*Kp2*Kp3)/(((Hplus+Kp1) &                      ! +/- PO4t 
                    !    *Hplus+Kp1*Kp2)*Hplus+Kp1*Kp2*Kp3)  &                                               ! +/- PO4t 
                    !+(-0.5*mn_rd-0.5*mn_rd2-0.5*fe_rd -hs_ox-fes_prec+fes_diss-mns_prec  &           ! +/- H2S (tot)  
                    !    +0.5*Disprop -sulfido +s23_rd-feS2_form)*Kh2s1/(Kh2s1+Hplus)  &              ! +/- H2S (tot)  
                    + 2.*(s4_rd + s23_rd) &      ! (CH2O)106(NH3)16H3PO4 + 53SO42- = 106HCO3- + 16NH3 + H3PO4 + 53H2S (Boudreau, 1996)                    
                    !+(self%k_Sipart_diss*Sipart)*KSi/(KSi+Hplus) &      !  ! +/- Si (tot)  

                    +     mn_ox   &  ! 4Mn2+ + O2 + 4H+ -> 4Mn3+ + 2H2O
                    - 3.*mn_ox2  &   ! 2Mn3+ + 3H2O  + 0.5 O2 -> 2MnO2 + 6H+
                    + 3.*mn_rd   &   ! 2MnO2 + 7H+ + HS- -> 2Mn3+ + 4H2O + S0
                    - 1.*mn_rd2  &   ! 2Mn3+ + HS- -> 2Mn2+ + S0 + H+
                    - 2.*mns_prec  & ! Mn2+ + H2S <-> MnS + 2H+                                           !!!!
                    + 2.*mns_diss  &  ! 
                    - 2.*mnco3_prec & ! Mn2+ + CO3-- <-> MnCO3 (vanCappelen,96):  
                    + 2.*mnco3_diss &
                    + 26.5*(DcDM_Mn +DcPM_Mn ) & ! & !DcDM_Mn is in N-units, i.e. 424/16
!!% (CH2O)106(NH3)16(H3PO4) +212MnO2 + 318CO2 +106H2O ->  424HCO3- + 212Mn2+ +16NH3 +H3PO4  (Boudreau, 1996)
                    - 2.*fe_ox  &    ! 4Fe2+ + O2 +10H2O-> 4Fe(OH)3 +8H+    (vanCappelen,96): 
                    - 1.*fe_ox2  &   ! 2Fe2+ + MnO2 +4H2O -> 2Fe(OH)3 + Mn2+ +2H+  (vanCappelen,96):                                          !!!!
                    + 2.*fe_rd  &   ! 2Fe(OH)3 + HS- + 5H+ ->  2Fe2++S0+6H2O  (here and below d(Alk_H2S) is excluded, as give before)                                        !!!!
                    - 1.*fes_prec  & ! Fe2+ + H2S  <-> FeS + H+ (Bektursunova, 11)
                    + 1.*fes_diss  & ! 
                    - 2.*fes_ox &  ! FeS + 2.25O2 +H2O -> 0.5Fe2O3 + 2H+ +SO42- (Soetaert,07)
                                   ! FeS + H2S -> FeS2 + H2     (Rickard, 1997,Meysmann,2003,Soetaert,07)
                    - 2.*fes2_ox & !  FeS2 + 3.5 O2 + H2O -> Fe2+ + 2 SO42- + 2 H+ (Wijsman,02)
                    + 53.*(DcDM_Fe +DcPM_Fe) & ! !DcDM_Fe is in N-units, i.e. 848/16
!% (CH2O)106(NH3)16H3PO4 + 424Fe(OH)3 + 742CO2 -> 848HCO3-+ 424Fe2+ +318 H2O +16NH3 +H3PO4  (Boudreau, 1996)        
                    - 0.5*Disprop  & ! 4S0 + 3H2O -> 2H2S + S2O3-- + 2H+                                         !!!!
                    - 1.*(-s0_ox)  & ! 2S0 + O2 + H2O -> S2O3-- + 2H+ 
                    - 0.5*s0_no3  &  ! 4S0 + 3NO3- + 7H2O -> 4SO4-- + 3NH4+ + 2H+ 
                    - 1.*s23_ox  &   ! S2O3-- + 2O2 + 2OH- -> 2SO4-- + H2O
                    - 0.4*sulfido &  ! 5H2S + 8NO3- + 2OH+ -> 5SO4-- + 4N2 + 6H2O (Volkov, 1984) 
                    - 2.*CaCO3_prec &   ! Ca2+ + CO32- -> CaCO3                                        !!!!
                    + 2.*CaCO3_diss &   ! CaCO3 -> Ca2+ + CO32- 
!                    + 0.15*(autolis-DcDM_O2-DcDM_NOX-DcDM_SO4 -DcDM_Mn -DcDM_Fe -HetBhae-HetBhan)*7. & ! part of OM charged (EPS etc.)
                    )
    _SET_ODE_(self%id_Alk,dAlk)

    _SET_DIAGNOSTIC_(self%id_DcPM_O2,DcPM_O2)  
    _SET_DIAGNOSTIC_(self%id_DcDM_O2,DcDM_O2)        
    _SET_DIAGNOSTIC_(self%id_DcPM_NOX,DcPM_NOX)
    _SET_DIAGNOSTIC_(self%id_DcDM_NOX,DcDM_NOX)
    _SET_DIAGNOSTIC_(self%id_DcDM_SO4,DcPM_SO4)
    _SET_DIAGNOSTIC_(self%id_DcPM_SO4,DcPM_SO4)
    _SET_DIAGNOSTIC_(self%id_DcPM_Mn,DcPM_Mn)
    _SET_DIAGNOSTIC_(self%id_DcDM_Mn,DcDM_Mn)
    _SET_DIAGNOSTIC_(self%id_DcDM_Fe,DcPM_Fe)
    _SET_DIAGNOSTIC_(self%id_DcPM_Fe,DcPM_Fe)
    _SET_DIAGNOSTIC_(self%id_ChemBaae,ChemBaae)
    _SET_DIAGNOSTIC_(self%id_HetBhae,HetBhae)
    _SET_DIAGNOSTIC_(self%id_ChemBaan,ChemBaan)
    _SET_DIAGNOSTIC_(self%id_HetBhan,HetBhan)
    _SET_DIAGNOSTIC_(self%id_MortBaan,MortBaan)
    _SET_DIAGNOSTIC_(self%id_MortBhan,MortBhan)  
    _SET_DIAGNOSTIC_(self%id_MortBaae,MortBaae)
    _SET_DIAGNOSTIC_(self%id_MortBhae,MortBhae)
    _SET_DIAGNOSTIC_(self%id_autolis,autolis)
    _SET_DIAGNOSTIC_(self%id_anammox ,anammox)
    _SET_DIAGNOSTIC_(self%id_Denitr1_PM ,Denitr1_PM)
    _SET_DIAGNOSTIC_(self%id_Denitr1_DM ,Denitr1_DM)
    _SET_DIAGNOSTIC_(self%id_Denitr2_PM ,Denitr2_PM)
    _SET_DIAGNOSTIC_(self%id_Denitr2_DM ,Denitr2_DM)
    _SET_DIAGNOSTIC_(self%id_Denitr1 ,Denitr1)
    _SET_DIAGNOSTIC_(self%id_Denitr2 ,Denitr2)
    _SET_DIAGNOSTIC_(self%id_fe_ox,fe_ox)
    _SET_DIAGNOSTIC_(self%id_fe_rd,fe_rd)
    _SET_DIAGNOSTIC_(self%id_s23_no3,s23_no3)
    _SET_DIAGNOSTIC_(self%id_s0_no3,s0_no3)
    _SET_DIAGNOSTIC_(self%id_s23_ox,s23_ox)
    _SET_DIAGNOSTIC_(self%id_s0_ox,s0_ox)
    _SET_DIAGNOSTIC_(self%id_s23_rd_PM,s23_rd_PM)
    _SET_DIAGNOSTIC_(self%id_s23_rd_DM,s23_rd_DM)
    _SET_DIAGNOSTIC_(self%id_s4_rd_PM,s4_rd_PM)
    _SET_DIAGNOSTIC_(self%id_s4_rd_DM,s4_rd_DM)
    _SET_DIAGNOSTIC_(self%id_s4_rd,s4_rd)
    _SET_DIAGNOSTIC_(self%id_s23_rd,s23_rd)
    _SET_DIAGNOSTIC_(self%id_Disprop,Disprop)
    _SET_DIAGNOSTIC_(self%id_hs_ox,hs_ox)
    _SET_DIAGNOSTIC_(self%id_sulfido,sulfido)
    _SET_DIAGNOSTIC_(self%id_fe_ox2,fe_ox2)
    _SET_DIAGNOSTIC_(self%id_mn_ox,mn_ox)
    _SET_DIAGNOSTIC_(self%id_mn_ox2,mn_ox2)
    _SET_DIAGNOSTIC_(self%id_mn_rd,mn_rd)
    _SET_DIAGNOSTIC_(self%id_mn_rd2,mn_rd2) 
    _SET_DIAGNOSTIC_(self%id_Nitrif1,Nitrif1)
    _SET_DIAGNOSTIC_(self%id_Nitrif2,Nitrif2)
    _SET_DIAGNOSTIC_(self%id_mns_diss,mns_diss)
    !_SET_DIAGNOSTIC_(self%id_FeS_ox,FeS_ox)
    !_SET_DIAGNOSTIC_(self%id_fes_prec,fes_prec) 
    !!_SET_DIAGNOSTIC_(self%id_mns_ox,mns_ox)
    _SET_DIAGNOSTIC_(self%id_mns_prec,mns_prec)
    !_SET_DIAGNOSTIC_(self%id_CaCO3_diss,CaCO3_diss) 
    _SET_DIAGNOSTIC_(self%id_mnco3_diss,mnco3_diss)    
    _SET_DIAGNOSTIC_(self%id_mnco3_prec,mnco3_prec)    

    _LOOP_END_

   end subroutine do
!EOC

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

   end module fabm_niva_brom_redox

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

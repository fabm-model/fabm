!$Id$
#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_rolm --- Yakushev biogeochemical model ROLM \label{sec:bio-rolm}
!
! !INTERFACE:
   module bio_rolm
!
! !DESCRIPTION:
!  The goal of this biogeochemical O-N-S-P-Mn-Fe RedOxLayer Model (ROLM)
!  (Yakushev et al., 2006) is to study the processes of these chemical 
!  elements cycling in the water column with suboxic and anoxic conditions.
!  The biogeochemical processes of organic matter (OM) formation 
!  (photosynthesis and chemosynthesis), OM decay (aerobic decomposition, 
!  denitrification, sulphate reduction), nitrification, anammox, hydrogen 
!  sulfide oxidation, the reduction and oxidation of manganese and iron 
!  species and the transformation of phosphorus species were parameterized.
!  Model consists of the following stat variables: Dissolved oxygen (O2), 
!  hydrogen sulfide (H2S), elemental sulfur (S0), thiosulfate (S2O3), 
!  sulfate (SO4), ammonia (NH4), nitrite (NO2), nitrate (NO3), particulate
!  organic nitrogen (PON), dissolved organic nitrogen (DON), phosphate 
!  (PO4), particulate organic phosphorus (POP), dissolved organic 
!  phosphorus (DOP), bivalent manganese (MnII), trivalent manganese 
!  (MnIII), quadrivalent manganese (MnIV), bivalent iron (FeII), trivalent
!  iron (FeIII), phytoplankton (Phy), zooplankton (Zoo), aerobic 
!  heterotrophic bacteria (Bhe), aerobic autotrophic bacteria (Bae), 
!  anaerobic heterotrophic bacteria (Bha), anaerobic autotrophic bacteria
!  (Baa). The concentrations are in uM (for O,N,P,S species) and 
!  mg*WetWeight*m-3 (for biological parameters). A detailed description 
!  of the model is available at
!  http://www.io-warnemuende.de/documents/mebe68\_2007-yakushev.pdf.
!  This implemented version of ROLM model includes slight modifications 
!  and was included into GOTM by E.Yakushev, O.Podymov and I.Kuznetsov
!  in September 2007.
!
!  when WRITEFINISH is defined, results of modeling are written into 
!  'finish' file (defined in bio\_rolm.nml) every end of a modeled year.
!  If you want to initialize a new calculation with this file, define
!  it as 'chem\_init' in bio\_rolm.nml or replace original initialization 
!  file with it.
#define WRITEFINISH

! !USES:
!  default: all is private.
   use bio_var
   use time,         only: calendar_date, julian_day, julianday, secondsofday
   use observations, only: read_profiles
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public init_bio_rolm, init_var_rolm,       &
          surface_fluxes_rolm, light_rolm, do_bio_rolm, end_bio_rolm
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): E.Yakushev, O.Podymov and I.Kuznetsov
!
!  Moscow, 1992 ! RV Knorr 15/03-2003 ! IOW 05/02-2007
!
! !LOCAL VARIABLES:
   logical, public   :: writeFinish
   logical           :: fluff
   REALTYPE, public  :: kc
   REALTYPE          :: i_min
   REALTYPE, public  :: iv
   REALTYPE          :: a0
   REALTYPE          :: a1
   REALTYPE          :: a2
   REALTYPE          :: aa
   REALTYPE          :: g2
   REALTYPE          :: pvel
   REALTYPE          :: sfl_po
   REALTYPE          :: sfl_am
   REALTYPE          :: sfl_ni
   character(len=72) :: chem_init ! initialization file with start parameters
   character(len=72) :: finish    ! file with annual results of model calculation
!----Phy  ----------!
   REALTYPE          :: KNF       ! Maximum specific growth rate
   REALTYPE          :: k_Erlov   ! Extinction coefficient
   REALTYPE          :: Io        ! Optimal irradiance 
   REALTYPE          :: Iopt      ! Optimal irradiance 
   REALTYPE          :: bm        ! Coefficient for dependence on t
   REALTYPE          :: cm        ! Coefficient for dependence on t
   REALTYPE          :: KFN       ! Specific respiration rate
   REALTYPE          :: KFP       ! Specific rate of mortality
   REALTYPE          :: KFD       ! Specific rate of excretion
!----Zoo -----------!
   REALTYPE          :: KFZ       ! Max.spec. rate of grazing of Zoo on Phy 
   REALTYPE          :: KF        ! Half-sat.const. for grazing of Zoo on Phy for Phy/Zoo ratio
   REALTYPE          :: FP        !
   REALTYPE          :: Phypor    !
   REALTYPE          :: KPZ       ! Max.spec.rate of grazing of Zoo on POP
   REALTYPE          :: KPP       ! Half-sat.const. for the graz of Zoo on POP for POP/Zoo  ratio
   REALTYPE          :: POPpor    !
   REALTYPE          :: KZN       ! Specific respiration rate
   REALTYPE          :: KZP       ! Maximum specific rate of mortality of Zoo
   REALTYPE          :: Uz        ! Food absorbency for Zoo
   REALTYPE          :: Hz        ! Ratio betw. diss. and part. excretes of Zoo
!---- OM --------!
   REALTYPE          :: KPD       ! Specific rate of decomposition of POM to DOM
!---- P ---------!
   REALTYPE          :: KPO4      ! Half-sat. constant for uptake of PO4 by Phy
!---- N ---------!
   REALTYPE          :: KPSI      ! Strength of NH4 inhibition of NO3 uptake constant = 1.46  (Wroblev,G,01)
   REALTYPE          :: KNO3      ! Half-sat.const.for uptake of NO3+NO2 =0.5 (G,01)
   REALTYPE          :: KNH4      ! Half-sat.const.for uptake of NH4=0.2 (G,01)
   REALTYPE          :: KND4      ! Spec.rate of DON decomp. 0.002(S,W,96)0.1-1(W,K,91)
   REALTYPE          :: KNP4      ! Spec.rate of PON decomp. =0.002 (S,W,96)
   REALTYPE          :: KN42      ! Spec.rate of 1st st. of nitrification =0.01(S&W,96)0.1(G,01)
   REALTYPE          :: KN23      ! Spec.rate of 2d st. of nitrification
   REALTYPE          :: KN32      ! Spec.rate of 1 stage of denitrif =0.16(Y,98),0.5(S&W,96),0.015(G,01)
   REALTYPE          :: KN24      ! Spec.rate of 1 stage of denitrif =0.22 (Y,98)
   REALTYPE          :: KT        ! Spec.rate of thiodenitrification
   REALTYPE          :: k_annamox !Spec.rate of anammox
!---- S ---------!
   REALTYPE          :: K_hs_ox   ! Specific rate of oxidation of H2S with O2
   REALTYPE          :: K_s0_ox   ! Specific rate of oxidation of S0 with O2
   REALTYPE          :: K_s23_ox  ! Specific rate of oxidation of S2O3 with O2
   REALTYPE          :: K_s4_rd   ! Specific rate of OM sulfate reduction with sulfate
   REALTYPE          :: K_s23_rd  ! Specific rate of OM sulfate reduction with thiosulfate
!---- Mn --------!
   REALTYPE          :: K_mn_ox   ! MnII oxidation with O2 constant
   REALTYPE          :: K_mn_ox2  ! MnIII oxidation with O2 constant
   REALTYPE          :: K_mn_rd   ! MnIV reduction with Sulfide constant
   REALTYPE          :: K_mn_rd2  ! MnIII reduction with Sulfide constant
!---- Fe --------!
   REALTYPE          :: K_fe_ox   ! Fe oxidation with O2 constant (K,05)
   REALTYPE          :: K_fe_nox  ! Fe oxidation with NO3 constant (K,05)
   REALTYPE          :: K_fe_mnox ! Fe oxidation with MnIV constant (K,05)
   REALTYPE          :: K_fe_rd   ! FeIII reduction by sulfide (K,05)
!  Stochiometric coefficients
   REALTYPE          :: Sp        ! P[uM]/BIOMASS [mg/m3]
   REALTYPE          :: SpZ       ! P[uM]/BIOMASS ZOO [mg/m3]
   REALTYPE          :: Sn        ! N[uM]/BIOMASS [mg/m3]
   REALTYPE          :: Ss        ! S[uM]/BIOMASS [mg/m3] of bacteria, i.e. during chemosynthesis
   REALTYPE          :: Sc        ! C[uM]/BIOMASS [mg/m3] for bact should be *.5 due to large P content
   REALTYPE          :: OkP       ! O[uM]/P[uM]
   REALTYPE          :: NkP       ! N[uM]/P[uM]
   REALTYPE          :: OkN       ! O[uM]/N[uM]
   REALTYPE          :: SkP       ! S[uM]/P[uM]
!  References: Yakushev, 1998; Savchuk, Wulff, 1996: Ward, Kilpatrik, 1991; 
!              Gregoire et al.2001, Konovalov et al., 2006;

   REALTYPE          :: Bu        ! Burial coeficient for lower boundary
   REALTYPE          :: Trel      ! Relaxation time for lower boundary
   REALTYPE          :: O2LimC    ! Limiting oxygen concentration for the lower boundary

   REALTYPE          :: LatLight  ! Latitude of modeled place (used for luminance calculation)

   REALTYPE          :: phy_initial,zoo_initial,bae_initial,bhe_initial,baa_initial,bha_initial,&
                        po4_initial,dop_initial,pop_initial,no3_initial,no2_initial,nh4_initial,&
                        pon_initial,don_initial,o2_initial,so4_initial,s2o3_initial,s0_initial,&
                        h2s_initial,mn4_initial,mn3_initial,mn2_initial,fe3_initial,fe2_initial

!  Sinking velosities [m d-1]:
   REALTYPE          :: w_phy,w_zoo,w_bae,w_bhe,w_baa,w_bha,w_mn4,w_fe3,w_pon,w_pop,w_s0

   REALTYPE          :: CNF,CNFI,CNFT,CPN,CDN,WF,WZ,WD,&
                        GrowthPhy,CZooPhy,CZooPOP,GrazBact,CZoo,LimNO3,LimNH4,LimP,&
                        LimN,LimT,LimLight,MortPhy,MortZoo,Autolis,AutolisP,AutolisN,Nitrif1,Nitrif2,&
                        Denitr1,Denitr2, PhosPOP,PhosDOP, ExcrPhy,Chemos,Hetero,Iz,&
                        hs_ox,s0_ox,s23_ox,s23_rd,s4_rd,mn_ox,mn_rd, mn_ox2,mn_rd2,&
                        Destr_OM, Sulfido, Sulfido2,mn_nox,SigmaN,snf,rnfp,&
                        Nfixation, anammox, Denitr1_PM, Denitr1_DM, Denitr2_PM,Denitr2_DM, &
                        s23_rd_PM,s23_rd_DM,s4_rd_PM,s4_rd_DM,ChemosA,HeteroA,&
                        MortBhet,MortBaut,MortBhetA,MortBautA,Disprop,fe_ox,fe_rd,fe_mnox,fe_nox,fe_n2ox,&
                        DcPM_O2,DcPM_NO3,DcPM_SO4,DcDM_O2,DcDM_NO3,DcDM_SO4,Mn_rd_DM,Mn_rd_PM,&
                        GrazPhy, GrazPOP, GrazBaut, GrazBhet, GrazBautA, GrazBheta, Grazing, AmmonPON,&
                        AmmonDON, NO3mi, NO2mi, O2nf, Bheta

!  parameters for "soft switches"
   REALTYPE          :: s_pho_po4
   REALTYPE          :: s_pho_noX
   REALTYPE          :: s_pho_nh4
   REALTYPE          :: s_anm_o2
   REALTYPE          :: s_nf1_O2
   REALTYPE          :: s_nf2_O2
   REALTYPE          :: s_omox_o2  ! threshold of o2 for OM mineralization
   REALTYPE          :: s_omno_o2  ! threshold of o2 for OM denitrification
   REALTYPE          :: s_omso_o2  ! threshold of o2 for OM sulfred.
   REALTYPE          :: s_omso_no  ! threshold of noX for OM sulfred.
   REALTYPE          :: s_zomr_hs  ! threshold of Zoo mortality
   REALTYPE          :: s_zobr_o2  ! threshold of Zoo breathing
   REALTYPE          :: s_mnox_mn2 ! threshold of Mn2 oxidation
   REALTYPE          :: s_mnox_mn3 ! threshold of Mn3 oxidation
   REALTYPE          :: s_mnrd_mn4 ! threshold of Mn4 reduction
   REALTYPE          :: s_mnrd_mn3 ! threshold of Mn3 reduction
   REALTYPE          :: s_feox_fe2 ! threshold of Fe2 reduction
   REALTYPE          :: s_bhe_on   ! threshold of TON for Bhe grouth
   REALTYPE          :: s_bbe_o2   ! threshold of O2 for Bhe & Bae mort.
   REALTYPE          :: s_bba_hs   ! threshold of H2S for Baa mortality
   REALTYPE          :: k_omox_o2  ! half sat. of o2 for OM miner.
   REALTYPE          :: k_omno_no3 ! half sat. of no3 for OM denitr.
   REALTYPE          :: k_omno_no2 ! half sat. of no2 for OM denitr.
   REALTYPE          :: k_mnoxO2   ! half sat. of Mn oxidation
   REALTYPE          :: k_mnrdHS   ! half sat. of Mn reduction
   REALTYPE          :: O2s_nf     ! half saturation for nitrification
   REALTYPE          :: s_bac_new  ! threshold for Bact re-apperence 
   REALTYPE          :: s_po4_srp  ! threshold for PO4 scaveging 

!  parameters for low boundary conditions
   REALTYPE          :: b_nh4, b_po4, b_h2s, b_mn2, b_fe2
   integer           :: out_unit
   REALTYPE, allocatable      :: ppi(:)
   integer, public, parameter :: o2=1, no3=2, no2=3, nh4=4, so4=5, s2o3=6, s0=7, h2s=8, mn4=9,&
                                 mn2=10, phy=11, zoo=12, mn3=13, po4=14, dop=15, pop=16,&
                                 don=17, pon=18, bae=19, bhe=20, baa=21, bha=22,&
                                 fe3=23, fe2=24, temp=25, salt=26


!------------------------------ namelists ------------------------

 namelist /bio_rolm_factors_nml/ &
   fluff, kc, i_min, iv,&
!----Phy ----------!
   KNF, k_Erlov, Io, Iopt, bm, cm, KFN, KFP, KFD,&
!----Zoo-----------!
   KFZ, KF, FP, Phypor, KPZ, KPP, POPpor, KZN, KZP, Uz, Hz,&
!----OM ---------!
   KPD,&
!----P---------!
   KPO4,&
!----N---------!
   KPSI, KNO3, KNH4, KND4, KNP4, KN42, KN23, KN32, KN24, KT, k_annamox, &
!----S---------!
   K_hs_ox, K_s0_ox, K_s23_ox, K_s4_rd, K_s23_rd,&
!----Mn---------!
   K_mn_ox, K_mn_ox2, K_mn_rd, K_mn_rd2,&
!----Fe---------!
   K_fe_ox, K_fe_nox, K_fe_mnox, K_fe_rd,& 
!Stoichiometric coefficients
   Sp, SpZ, Sn, Ss, Sc, OkP, NkP, OkN, SkP,&
! Burial coeff. and time of relaxation
   Bu, Trel,&
! Limit for oxygen concentration on the lower boundary
   O2LimC,&
! Latitude of modeled place
   LatLight

 namelist /bio_rolm_nml/ &
   numc,&
   w_phy,w_zoo,w_bae,w_bhe,w_baa,w_bha,w_mn4,w_fe3,w_pon,w_pop,w_s0,&
   surface_flux_method,a0,a1,a2,g2,aa,pvel,sfl_po,sfl_am,sfl_ni, chem_init, finish

 namelist /bio_rolm_switches_nml/ &
   s_pho_po4, s_pho_noX, s_pho_nh4, s_anm_o2, s_nf1_O2, s_nf2_O2, s_omox_o2, &
   s_omno_o2, s_omso_o2, s_omso_no, s_zomr_hs, s_zobr_o2, s_mnox_mn2, s_mnox_mn3, &
   s_mnrd_mn4, s_mnrd_mn3, s_feox_fe2, s_bhe_on, s_bbe_o2, s_bba_hs, mortality, &
   k_omox_o2, k_omno_no3, k_omno_no2, k_mnoxO2, k_mnrdHS, O2s_nf, s_bac_new, s_po4_srp

 namelist /bio_rolm_lowbound_nml/ &
   b_nh4, b_po4, b_h2s, b_mn2, b_fe2

!EOP
!-----------------------------------------------------------------------

 contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the bio module
!
! !INTERFACE:
   subroutine init_bio_rolm(namlst,fname,unit)
!
! !DESCRIPTION:
!  Here the bio namelist {\tt bio\_rolm.nml} is read and
!  various variables (rates and settling velocities)
!  are transformed into SI units.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,          intent(in)   :: namlst
   character(len=*), intent(in)   :: fname
   integer,          intent(in)   :: unit
!
! !REVISION HISTORY:
!  Original author(s): E.Yakushev, O.Podymov & I.Kuznetsov
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC
   LEVEL2 'init_bio_rolm'

   numc=24
   open(namlst,file=fname,action='read',status='old',err=98)
   read(namlst,nml=bio_rolm_nml,err=99)
   read(namlst,nml=bio_rolm_factors_nml,err=99)
   read(namlst,nml=bio_rolm_switches_nml,err=99)
   read(namlst,nml=bio_rolm_lowbound_nml,err=99)
   close(namlst)
   n_surface_fluxes = 4
#ifdef WRITEFINISH 
   writeFinish = .false.
#endif

   LEVEL3 'namelist "', fname, '" read'


!  Conversion from day to second
   w_phy = w_phy / secs_pr_day
   w_zoo = w_zoo / secs_pr_day
   w_bae = w_bae / secs_pr_day
   w_bhe = w_bhe / secs_pr_day
   w_baa = w_baa / secs_pr_day
   w_bha = w_bha / secs_pr_day
   w_mn4 = w_mn4 / secs_pr_day
   w_fe3 = w_fe3 / secs_pr_day
   w_pon = w_pon / secs_pr_day
   w_pop = w_pop / secs_pr_day
   w_s0  = w_s0 / secs_pr_day
   pvel  = pvel / secs_pr_day

!  initialize variable descriptions

   call bio_alloc_info


   var_names(1) = 'o2'
   var_units(1) = 'mmol o/m**3'
   var_long(1)  = 'dissolved_oxygen'

   var_names(2) = 'no3'
   var_units(2) = 'mmol n/m**3'
   var_long(2)  = 'nitrate'

   var_names(3) = 'no2'
   var_units(3) = 'mmol n/m**3'
   var_long(3)  = 'nitrite'

   var_names(4) = 'nh4'
   var_units(4) = 'mmol n/m**3'
   var_long(4)  = 'ammonia'
   var_names(5) = 'so4'
   var_units(5) = 'mmol s/m**3'
   var_long(5)  = 'sulfate'

   var_names(6) = 's2o3'
   var_units(6) = 'mmol s/m**3'
   var_long(6)  = 'thiosulfate'

   var_names(7) = 's0'
   var_units(7) = 'mmol s/m**3'
   var_long(7)  = 'elemental_sulfur'

   var_names(8) = 'h2s'
   var_units(8) = 'mmol s/m**3'
   var_long(8)  = 'hydrogen_sulfide'

   var_names(9) = 'mn4'
   var_units(9) = 'mmol mn/m**3'
   var_long(9)  = 'quadrivalent_manganese'

   var_names(10) = 'mn2'
   var_units(10) = 'mmol mn/m**3'
   var_long(10)  = 'bivalent_manganese'

   var_names(11) = 'phy'
   var_units(11) = 'mgWW m**-3'
   var_long(11)  = 'phytoplankton'

   var_names(12) = 'zoo'
   var_units(12) = 'mgWW m**-3'
   var_long(12)  = 'zooplankton'

   var_names(13) = 'mn3'
   var_units(13) = 'mmol mn/m**3'
   var_long(13)  = 'trivalent_manganese'

   var_names(14) = 'po4'
   var_units(14) = 'mmol p/m**3'
   var_long(14)  = 'phosphate'

   var_names(15) = 'dop'
   var_units(15) = 'mmol p/m**3'
   var_long(15)  = 'dissolved_organic_phosphorus'

   var_names(16) = 'pop'
   var_units(16) = 'mmol p/m**3'
   var_long(16)  = 'particulate_organic_phosphorus'

   var_names(17) = 'don'
   var_units(17) = 'mmol n/m**3'
   var_long(17)  = 'dissolved_organic_nitrogen'

   var_names(18) = 'pon'
   var_units(18) = 'mmol n/m**3'
   var_long(18)  = 'particulate_organic_nitrogen'

   var_names(19) = 'bae'
   var_units(19) = 'mgWW m**-3'
   var_long(19)  = 'aerobic_autotrophic_bacteria'

   var_names(20) = 'bhe'
   var_units(20) = 'mgWW m**-3'
   var_long(20)  = 'aerobic_heterotrophic_bacteria'

   var_names(21) = 'baa'
   var_units(21) = 'mgWW m**-3'
   var_long(21)  = 'anaerobic_autotrophic_bacteria'

   var_names(22) = 'bha'
   var_units(22) = 'mgWW m**-3'
   var_long(22)  = 'anaerobic_heterotrophic_bacteria'

   var_names(23) = 'fe3'
   var_units(23) = 'mmol fe/m**3'
   var_long(23)  = 'trivalent_iron'

   var_names(24) = 'fe2'
   var_units(24) = 'mmol fe/m**3'
   var_long(24)  = 'bivalent_iron'



   out_unit=unit

   LEVEL3 'module initialized'

   return

98 LEVEL2 'I could not open bio_rolm.nml'
   LEVEL2 'If thats not what you want you have to supply bio_rolm.nml'
   LEVEL2 'See the bio example on www.gotm.net for a working bio_rolm.nml'
   return
99 FATAL 'I could not read bio_rolm.nml'
   stop 'init_bio_rolm'

   end subroutine init_bio_rolm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the concentration variables
!
! !INTERFACE:
   subroutine init_var_rolm
!
! !DESCRIPTION:
!  Here, the the initial conditions are set and the settling velocities are
!  transferred to all vertical levels. All concentrations are declared
!  as non-negative variables. The initial arrays are read from a file "chem\_init.dat".
!  Furthermore, the primary production {\tt ppi} is allocated.
!
! !USES:
   IMPLICIT NONE

!
! !REVISION HISTORY:
!  Original author(s): E.Yakushev, O.Podymov & I.Kuznetsov

! !LOCAL VARIABLES:
  REALTYPE               :: z(0:nlev), cc1(0:nlev, numc)
  integer                :: i,rc,iret,j
  integer                :: yy, mm, dd, hh, min, ss
!EOP
!-----------------------------------------------------------------------
!BOC
! ------ Reading of the BGC data from an array ---------------
   open ( UNIT   = 198 &
         ,FILE   = chem_init &
         ,STATUS = 'OLD' &
         ,ACTION = 'READ' &
         ,IOSTAT = iret &
        )
   if(iret /= 0 ) then
      write(*,*) 'Error Opening File ',chem_init,' ',iret; STOP
   end if

   do i=1, nlev
     z(i) = -nlev+i-0.5 ! 1-meter depth step
   enddo

   call read_profiles(198, nlev, numc, yy, mm, dd, hh, min, ss, z, cc1, i, iret)

   do j=1, nlev
     do i=1, numc
       cc(i,  j) = cc1(j ,i)
     enddo
   enddo

  close (UNIT = 198, IOSTAT = iret)

  if(iret /= 0 ) then
    write(*,*) 'Error Closing File (chem_init.dat)',iret; STOP
  end if

! ------ new initial -------------

   do i=0,nlev
     ws(phy,i) = w_phy
     ws(zoo,i) = w_zoo
     ws(bae,i) = w_bae
     ws(bhe,i) = w_bhe
     ws(baa,i) = w_baa
     ws(bha,i) = w_bha

     ws(po4,i) = _ZERO_
     ws(dop,i) = _ZERO_
     ws(pop,i) = _ZERO_
     ws(no3,i) = _ZERO_
     ws(no2,i) = _ZERO_
     ws(nh4,i) = _ZERO_

     ws(pon,i) = w_pon
     ws(pop,i) = w_pop
     ws(o2,i)  = _ZERO_
     ws(so4,i) = _ZERO_
     ws(s2o3,i)= _ZERO_
     ws(s0,i)  = w_s0

     ws(h2s,i) = _ZERO_
     ws(mn4,i) = w_mn4
     ws(mn3,i) = _ZERO_
     ws(mn2,i) = _ZERO_
     ws(fe3,i) = w_fe3
     ws(fe2,i) = _ZERO_
   end do

   sfl = _ZERO_

   posconc(phy) = 1
   posconc(zoo) = 1
   posconc(bae) = 1
   posconc(bhe) = 1
   posconc(baa) = 1
   posconc(bha) = 1

   posconc(po4) = 1
   posconc(dop) = 1
   posconc(pop) = 1
   posconc(no3) = 1
   posconc(no2) = 1
   posconc(nh4) = 1

   posconc(pon) = 1
   posconc(don) = 1
   posconc(o2)  = 1
   posconc(so4) = 1
   posconc(s2o3)= 1
   posconc(s0)  = 1

   posconc(h2s) = 1
   posconc(mn4) = 1
   posconc(mn3) = 1
   posconc(mn2) = 1
   posconc(fe3) = 1
   posconc(fe2) = 1

   allocate(ppi(0:nlev),stat=rc)
   if (rc /= 0) stop 'init_var_rolm(): Error allocating ppi)'
   select case (surface_flux_method)
      case (-1)! absolutely nothing
      case (0) ! constant
         sfl(po4)= -sfl_po / secs_pr_day
         sfl(nh4)= -sfl_am / secs_pr_day
         sfl(no2)= -sfl_ni / secs_pr_day
      case (2) ! from file via sfl_read

      case default
   end select

   LEVEL3 'ROLM variables initialised .'

   return
   end subroutine init_var_rolm
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Step function
!
! !INTERFACE:
   REALTYPE function th(x,w,min,max)
!
! !DESCRIPTION:
!  Instead of the
!  heavyside switches used by e.g.\ \cite{Neumannetal2002}, we apply here 
!  a smoothed {\it tangens hyperbolicus} transition 
!  with prescribed width $x_w$:
!  \begin{equation}\label{theta}
!  \theta (x,x_w,y_{\min},y_{\max})= y_{\min}+(y_{\max}-y_{\min})
!  \frac12\left(1-\tanh \left(\frac{x}{x_w}   \right)      \right).
!  \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: x,w,min,max
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (w .gt. 1.e-10) then
      th=min+(max-min)*0.5*(1.+tanh(x/w))
   else
      if (x .gt. _ZERO_) then
         th=_ONE_
      else
         th=_ZERO_
      end if
   end if
   return
   end function th
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Saturation function squared
!
! !INTERFACE:
   REALTYPE function yy(a,x)
!
! !DESCRIPTION:
! This is a squared Michaelis-Menten type of limiter:
! \begin{equation}\label{Y}
! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   REALTYPE, intent(in)                :: a,x
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   yy=x**2/(a**2+x**2)
   return
   end function yy
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for the ROLM model
!
! !INTERFACE:
   subroutine surface_fluxes_rolm(nlev,t)
!
! !DESCRIPTION:
! Here, those surface fluxes which have been read from a file are transformed
! to SI units, and the surface oxygen flux is calculated by means of the 
! following formula:
! \begin{equation}\label{o2flux}
! F^s_9 = p_{vel} \left(O_{sat}-c_9 \right)
! \end{equation}
! with
! \begin{equation}\label{osat}
! O_{sat}= a_0\left(a_1-a_2T  \right).
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                              :: nlev
   integer                              :: i
   integer                              :: ci
   REALTYPE, intent(in)                 :: t
!
! !REVISION HISTORY:
!  Original author(s): E.Yakushev, O.Podymov & I.Kuznetsov
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
!BOC

!  Surface fluxes are constant for now.
!  In the near future we plan make them depend on seasonal variability.

   select case (surface_flux_method)
      case (-1)! absolutely nothing
      case (0) ! constant
         sfl(po4)= -sfl_po / secs_pr_day
         sfl(nh4)= -sfl_am / secs_pr_day
         sfl(no2)= -sfl_ni / secs_pr_day

!     case (1) ! from file via sfl_read
!        sfl(po4)= sfl_po / secs_pr_day
!        sfl(nh4)= sfl_am / secs_pr_day
!        sfl(no2)= sfl_ni / secs_pr_day

!     case (2) ! from file via sfl_read
!         sfl(ni) = -sfl_read(1)/secs_pr_day
!         sfl(am) = -sfl_read(2)/secs_pr_day
!         sfl(po) = -sfl_read(3)/secs_pr_day 
      case (3) ! sfl array filled externally - for 3D models
      case default
   end select

! surface oxygen flux
       sfl(o2) = pvel*(a0*(a1-a2*t)-cc(o2,nlev))
   return
   end subroutine surface_fluxes_rolm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Light properties for the ROLM model
!
! !INTERFACE:
   subroutine light_rolm(nlev, bioshade_feedback)
!
! !DESCRIPTION:
! Here, the photosynthetically available radiation is calculated
! by simply assuming that the short wave part of the total
! radiation is available for photosynthesis. The user should make
! sure that this is consistent with the light class given in the
! {\tt extinct} namelist of the {\tt obs.nml} file.
! The self-shading effect is also calculated in this subroutine,
! which may be used to consider the effect of bio-turbidity also
! in the temperature equation (if {\tt bioshade\_feedback} is set
! to true in {\tt bio.nml}).
! For details, see section \ref{sec:do-bio}.
!
! !USES:
   use bio_var, only: bioshade_
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: nlev
   logical, intent(in)                 :: bioshade_feedback
!
! !REVISION HISTORY:
!  Original author(s): E.Yakushev, O.Podymov & I.Kuznetsov
!
! !LOCAL VARIABLES:
   integer                   :: i
   REALTYPE                  :: zz,add
!EOP
!-----------------------------------------------------------------------
!BOC
   zz = _ZERO_
   add = _ZERO_
   do i=nlev,1,-1
      zz=zz+0.5*h(i)
      par(i)=rad(nlev)*(1.-aa)*exp(-zz/g2)*exp(-kc*add)
      zz=zz+0.5*h(i)
      if (bioshade_feedback) bioshade_(i)=exp(-kc*add)
   end do

   return
   end subroutine light_rolm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of the ROLM geobiochemical model\label{sec:bio-rolm-details}
!
! !INTERFACE:
   subroutine do_bio_rolm(first,numc,nlev,cc,pp,dd)
!
! !DESCRIPTION:
! The right hand sides of the \cite{Neumannetal2002} biogeochemical model are 
! coded in this soubroutine.
! First of all, based on (\ref{theta}) and (\ref{Y}), 
! we construct limiters for chemical
! reactions which depend on the availability of oxygen ($c_9$) and nitrate
! ($c_7$) and have to add up to unity:
!
! !USES:
   use time
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer                             :: numc,nlev,iret
   REALTYPE                            :: cc(1:numc,0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
   logical                             :: first
!
! !OUTPUT PARAMETERS:
   REALTYPE, intent(out)               :: pp(1:numc,1:numc,0:nlev)
   REALTYPE, intent(out)               :: dd(1:numc,1:numc,0:nlev)
!
! !REVISION HISTORY:
!  Original author(s): E.Yakushev, O.Podymov & I.Kuznetsov
!
! !LOCAL VARIABLES:
   REALTYPE, save                      :: iopt
   REALTYPE                            :: rat(0:nlev,0:nlev)
   REALTYPE                            :: psum,llda,llan,llsa,r1,r2,r3
   REALTYPE                            :: wo=30.,wn=0.1,dot2=0.2
   REALTYPE                            :: thopnp,thomnp,thomnm,thsum
   integer                             :: i,j,ci
   REALTYPE                            :: depthCi(nlev)
   integer                             :: y,m,d,julY,deny

!EOP
!----------------------------------------
!BOC
   if (first) then
      first =.false.
      iopt  =max(0.25*I_0,I_min)
      do ci=1,nlev
        ppi(ci)=par(ci)/iopt*exp(1.-par(ci)/iopt)
      end do
   end if

   pp = _ZERO_
   dd = _ZERO_
   rat=1.         ! fixed (in time  space) ratio between sink and source

   depthCi(nlev)=0.5
   do i=nlev-1,1,-1
     depthCi(i)=depthCi(i+1)+1
   enddo

!  writing down end-of-year results for a new initialization file
#ifdef WRITEFINISH
   call calendar_date(julianday,y,m,d)
   m=1
   d=1
   call julian_day(y,m,d,julY)
   deny=julianday-julY
   if (deny==364) writeFinish=.true.
   if ((deny==0) .and. (writeFinish)) then
     LEVEL1 'writing finish.dat'
     open ( UNIT = 198 &
             , FILE = finish &
             , STATUS = 'REPLACE' &
             , ACTION = 'WRITE' &
             , IOSTAT = iret &
           )
     if(iret /= 0 ) then
       write(*,*) 'Error Opening File ',finish,' ',iret; STOP
     end if
!    yyyy/mm/dd 00:00:01 nlev 1 - the header of the 'finish' file - beginning of a year
     write (198,'(I4.4,"/",I2.2,"/",I2.2," ",I2.2,":",I2.2,":",I2.2," ",I4," ",I1)') y,m,d,0,0,1,nlev,1
     do i=1,nlev
       write(198,'(F6.1,24(F13.5,1x))') -nlev+i-0.5, cc(o2,i),cc(no3,i),cc(no2,i),cc(nh4,i),&
                                  cc(so4,i),cc(s2o3,i),cc(s0,i),cc(h2s,i),&
                                  cc(mn4,i),cc(mn2,i),cc(phy,i),cc(zoo,i),cc(mn3,i),&
                                  cc(po4,i),cc(dop,i),cc(pop,i),cc(don,i),cc(pon,i),&
                                  cc(bae,i),cc(bhe,i),cc(baa,i),cc(bha,i),&
                                  cc(fe3,i),cc(fe2,i)
     enddo

     close (UNIT = 198, IOSTAT = iret)
     if(iret /= 0 ) then
       write(*,*) 'Error Closing File ',finish,' ',iret; STOP
     end if
     writeFinish =.false.
   endif
!---------------------------------------- End of writing new initial
#endif

!---------------------------------------- Begin of main loop
   do ci=1,nlev
!---------------------------------------------------------------------------
! Calcualtion of Rci - biogeochemical production/destruction rates
! Coefficients and rates of processes

! dependence of particulate MnIV sinking rate on MnIV concentration:
     ws(mn4,ci)=w_mn4 * cc(mn4,ci)/(cc(mn4,ci)+1.)

!---------------------------------------------------------------------------
!========P=H=Y==============================================================
!---------------------------------------------------------------------------
     Iz            =Io*exp(-k_Erlov*depthCi(ci)) ! Irradiance changing with depth
     LimLight      =Iz/Iopt*exp(1-Iz/Iopt)       ! Influence of the Irradiance on photosynthesis
     LimT          =exp(bm*t(ci)-cm)             ! Influence of Temperature on photosynthesis

! dependence of photosynthesis on P
     LimP          =th(cc(po4,ci),s_pho_po4,_ZERO_,_ONE_)*&
                    ((cc(po4, ci))/cc(phy, ci))/(cc(po4, ci)/cc(phy, ci)+KPO4/25.) 

! dependence of photosynthesis on cc(nO3, ci),cc(nO2, ci):
     LimNO3        =th((cc(no3, ci)+cc(no2, ci)),s_pho_noX,_ZERO_,_ONE_)*&
                    (cc(no3, ci)+cc(no2, ci))/cc(phy, ci)/((cc(no3, ci)+cc(no2,ci))/cc(phy, ci)+KNO3/200.)

! dependence of photosynthesis on cc(nh4, ci):
     LimNH4        =th(cc(nh4, ci),s_pho_nh4,_ZERO_,_ONE_)*&
                    ((cc(nh4, ci))/cc(phy, ci))/(cc(nh4, ci)/cc(phy, ci)+KNH4/200.)
!----------------------------------------
     LimN          =LimNO3+LimNH4 ! Influence of N on photosynthesis

     GrowthPhy     =KNF*LimLight*LimT*dmin1(LimP,LimN)*cc(phy, ci) &
                    *cos((LatLight+(23.5*sin(2*3.14*deny/365.)))*3.14/180.)

     pp(phy,po4,ci)=GrowthPhy*(1.-KFN)
     pp(o2,phy,ci) =OkP*Sp*pp(phy,po4,ci)
     dd(po4,phy,ci)=Sp*pp(phy,po4,ci)
     dd(nh4,phy,ci)=Sn*pp(phy,po4,ci)*(LimNH4/LimN)
     dd(no2,phy,ci)=Sn*pp(phy,po4,ci)*(LimNO3/LimN)*(cc(no2, ci)/(cc(no2, ci)+cc(no3, ci)))
     dd(no3,phy,ci)=Sn*pp(phy,po4,ci)*(LimNO3/LimN)*(cc(no3, ci)/(cc(no2, ci)+cc(no3, ci)))

     !addition of cc(o2, ci) due to cc(no3, ci)->cc(nh4, ci) convertion
     pp(o2,no3,ci) =2.*Sn*pp(phy,po4,ci)*(cc(no3, ci)/(cc(no2, ci)+cc(no3, ci)))*(LimNO3/LimN)

     !addition of cc(o2, ci) due to cc(no2, ci)->cc(nh4, ci) convertion
     pp(o2,no2,ci) =0.5*Sn*pp(phy,po4,ci)*(cc(no2, ci)/(cc(no2, ci)+cc(no3, ci)))*(LimNO3/LimN)

     MortPhy       =KFP*cc(phy, ci) ! Rate of mortality of phy

     dd(phy,pop,ci)=MortPhy
     pp(pop,phy,ci)=Sp*dd(phy,pop,ci)
     pp(pon,phy,ci)=Sn*dd(phy,pop,ci)

     ExcrPhy=KFD*cc(phy, ci) ! Excretion of phy

     dd(phy,dop,ci)=ExcrPhy
     pp(dop,phy,ci)=Sp*dd(phy,dop,ci)
     pp(don,phy,ci)=Sn*dd(phy,dop,ci)

!---------------------------------------------------------------------------
!========Z=O=O==============================================================
!---------------------------------------------------------------------------
!---------------------------------------- Grazing of zoo on phy, detritus and bacteria
     GrazPhy       =KFZ*cc(zoo, ci)*(cc(phy, ci)/cc(zoo, ci))/(cc(phy, ci)/cc(zoo, ci)+KF/15.)

     dd(phy,zoo,ci)=GrazPhy
     pp(zoo,phy,ci)=Uz*dd(phy,zoo,ci)

     GrazPOP       =KPZ*cc(zoo, ci)*(cc(pop, ci)/cc(zoo, ci))/(cc(pop, ci)/cc(zoo, ci)+KPP/0.002)

     dd(pop,zoo,ci)=GrazPOP
     pp(zoo,pop,ci)=Uz*dd(pop,zoo,ci)

     GrazBaut      =KPZ*cc(zoo, ci)*(cc(bae, ci)/cc(zoo, ci))/(cc(bae, ci)/cc(zoo, ci)+2.5)

     dd(bae,zoo,ci)=GrazBaut
     pp(zoo,bae,ci)=Uz*dd(bae,zoo,ci)

     GrazBhet      =1.7*KPZ*cc(zoo, ci)*(cc(bhe, ci)/cc(zoo, ci))/(cc(bhe, ci)/cc(zoo, ci)+1.5)
     dd(bhe,zoo,ci)=GrazBhet
     pp(zoo,bhe,ci)=Uz*dd(bhe,zoo,ci)

     GrazBautA     =1.3*KPZ*cc(zoo, ci)*(cc(baa, ci)/cc(zoo, ci))/(cc(baa, ci)/cc(zoo, ci)+2.)
     dd(baa,zoo,ci)=GrazBautA
     pp(zoo,baa,ci)=Uz*dd(baa,zoo,ci)

     GrazBhetA     =KPZ*cc(zoo, ci)*(cc(bha, ci)/cc(zoo, ci))/(cc(bha, ci)/cc(zoo, ci)+2.)
     dd(bha,zoo,ci)=GrazBhetA
     pp(zoo,bha,ci)=Uz*dd(bha,zoo,ci)

!----------------------------------------
     GrazBact      =dd(bae,zoo,ci)+dd(bhe,zoo,ci)+dd(baa,zoo,ci)+dd(bha,zoo,ci)
     Grazing       =dd(phy,zoo,ci)+dd(pop,zoo,ci)+GrazBact !Grazing = GrazPhy+GrazPOP+GrazBact
!----------------------------------------

     MortZoo       =(0.8*th(cc(h2s,ci),s_zomr_hs,_ZERO_,_ONE_)+KZP)*cc(zoo, ci)*cc(zoo, ci)

     dd(zoo,pop,ci)=MortZoo
     pp(pop,zoo,ci)=Sp*(MortZoo+Grazing*(1.-Uz)*(1.-Hz))
     pp(pon,zoo,ci)=Sn*(MortZoo+Grazing*(1.-Uz)*(1.-Hz))

     pp(dop,zoo,ci)=Sp*Grazing*(1.-Uz)*Hz
     pp(don,zoo,ci)=Sn*Grazing*(1.-Uz)*Hz

     dd(zoo,po4,ci)=KZN*cc(zoo, ci)
     pp(po4,zoo,ci)=Sp*KZN*cc(zoo, ci)
     pp(nh4,zoo,ci)=Sn*KZN*cc(zoo, ci)

!    dd(o2,zoo,ci) =th(cc(o2,ci),s_zobr_o2,_ZERO_,_ONE_)*OkP*Sp*KZN*cc(zoo, ci)
     dd(o2,zoo,ci) =OkP*Sp*KZN*cc(zoo, ci)
     if (cc(o2, ci).lt.0.1) &
      dd(o2,zoo,ci)=0.

!---------------------------------------------------------------------------
!========O=M=-=Destruction==================================================
!---------------------------------------------------------------------------
     AmmonPON      =0.
     AmmonDON      =0.
     PhosPOP       =0.
     PhosDOP       =0.
     Denitr1       =0.
     Denitr2       =0.
     s4_rd         =0.
     s23_rd        =0.
     Destr_OM      =0.
     DcPM_O2       =0.
     DcDM_O2       =0.
     DcPM_NO3      =0.
     DcDM_NO3      =0.
     DcPM_SO4      =0.
     DcDM_SO4      =0.

!---------------------------------------- In oxic conditions:
! POM oxic mineralization: (CH2O)106(NH3)16H3PO4 + 106O2 -> 106CO2 + 16NH3 + H3PO4 + 106H2O :
     DcPM_O2       =th(cc(o2, ci),s_omox_o2,_ZERO_,_ONE_)*exp(0.15*t(ci)) &
                    *KNP4*cc(pon, ci)*(cc(o2, ci))/(cc(o2, ci)+k_omox_o2) 

! DOM oxic mineralization: (CH2O)106(NH3)16H3PO4 + 106O2 -> 106CO2 + 16NH3 + H3PO4 + 106H2O :
     DcDM_O2       =th(cc(o2, ci),s_omox_o2,_ZERO_,_ONE_)*exp(0.15*t(ci)) &
                    *KND4*cc(don, ci)*(cc(o2, ci))/(cc(o2, ci)+k_omox_o2)

     Destr_OM      =OkN*(DcPM_O2+DcDM_O2)
     dd(o2,pon,ci) =OkN*DcPM_O2
     dd(o2,don,ci) =OkN*DcDM_O2


!---------------------------------------- in suboxic conditions:
! OM denitrification (Richards, 1965): (CH2O)106(NH3)16H3PO4 + 84.8HNO3 = 106CO2 + 42.4N2 + 148.4H2O + 16NH3 + H3PO4  
! POM denitrification (1st stage) (Anderson et al., 1982): 1/2CH2O + NO3 -> NO2- + 1/2H2O + 1/2CO2 : 

     Denitr1_PM    =th(cc(o2, ci),s_omno_o2,_ONE_,_ZERO_) &
                     *KN32*(1.-cc(o2, ci)/(25.*(26.-cc(o2, ci)))) &
                     *cc(no3, ci)/(cc(no3, ci) +k_omno_no3)*cc(pon, ci)

!----------------------------------------
! DOM denitrification (1st stage): 1/2CH2O + NO3 -> NO2- + 1/2H2O + 1/2CO2 : 

     Denitr1_DM    =th(cc(o2, ci),s_omno_o2,_ONE_,_ZERO_) &
                    *KN32*(1.-cc(o2, ci)/(25.*(26.-cc(o2, ci)))) &
                    *cc(no3, ci)/(cc(no3, ci)+k_omno_no3)*cc(don, ci)

!----------------------------------------
! POM denitrification (2d stage): 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2

    Denitr2_PM     =th(cc(o2, ci),s_omno_o2,_ONE_,_ZERO_) &
                    *KN24*(1.-cc(o2, ci)/(25.*(26.-cc(o2, ci)))) &
                    *cc(no2, ci)/(cc(no2, ci)+k_omno_no2)*cc(pon, ci)

!----------------------------------------
! DOM denitrification (2d stage): 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2

     Denitr2_DM    =th(cc(o2, ci),s_omno_o2,_ONE_,_ZERO_) &
                    *KN24*(1.-cc(o2, ci)/(25.*(26.-cc(o2, ci)))) &
                    *cc(no2, ci)/(cc(no2, ci)+k_omno_no2)*cc(don, ci)


     Denitr1       =Denitr1_PM + Denitr1_DM
     Denitr2       =Denitr2_PM + Denitr2_DM


     pp(no2,no3,ci)=Denitr1
     dd(no3,no2,ci)=pp(no2,no3,ci)

! from Anderson et al.:
     DcPM_NO3      =16./212.*Denitr1_PM + 16./141.3*Denitr2_PM
     DcDM_NO3      =16./212.*Denitr1_DM + 16./141.3*Denitr2_DM


!---------------------------------------- in suboxic conditions:
! OM denitrification (Richards, 1965): (CH2O)106(NH3)16H3PO4 + 84.8HNO3 = 106CO2 + 42.4N2 + 148.4H2O + 16NH3 + H3PO4  
! POM denitrification (1st stage) (Anderson et al., 1982): 1/2CH2O + NO3 -> NO2- + 1/2H2O + 1/2CO2 : 

     Denitr1_PM    =th(cc(o2, ci),s_omno_o2,_ONE_,_ZERO_) &
                    *KN32*(1.-cc(o2, ci)/(25.*(26.-cc(o2, ci)))) &
                    *cc(no3, ci)/(cc(no3, ci) +k_omno_no3)*cc(pon, ci)

!---------------------------------------- 
! DOM denitrification (1st stage): 1/2CH2O + NO3 -> NO2- + 1/2H2O + 1/2CO2 :

     Denitr1_DM    =th(cc(o2, ci),s_omno_o2,_ONE_,_ZERO_) &
                    *KN32*(1.-cc(o2, ci)/(25.*(26.-cc(o2, ci)))) &
                    *cc(no3, ci)/(cc(no3, ci)+k_omno_no3)*cc(don, ci)

!---------------------------------------- 
! POM denitrification (2d stage): 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2

     Denitr2_PM    =th(cc(o2, ci),s_omno_o2,_ONE_,_ZERO_) &
                    *KN24*(1.-cc(o2, ci)/(25.*(26.-cc(o2, ci)))) &
                    *cc(no2, ci)/(cc(no2, ci)+k_omno_no2)*cc(pon, ci)

!---------------------------------------- 
! DOM denitrification (2d stage): 3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2

     Denitr2_DM    =th(cc(o2, ci),s_omno_o2,_ONE_,_ZERO_) &
                    *KN24*(1.-cc(o2, ci)/(25.*(26.-cc(o2, ci)))) &
                    *cc(no2, ci)/(cc(no2, ci)+k_omno_no2)*cc(don, ci)


     Denitr1       =Denitr1_PM + Denitr1_DM
     Denitr2       =Denitr2_PM + Denitr2_DM

     pp(no2,no3,ci)=Denitr1
     dd(no3,no2,ci)=pp(no2,no3,ci)

! from Anderson  et al. and Richards:
     DcPM_NO3      =16./212.*Denitr1_PM + 16./141.3*Denitr2_PM
     DcDM_NO3      =16./212.*Denitr1_DM + 16./141.3*Denitr2_DM

!---------------------------------------- in anoxic conditions:
! OM sulfatereduction (Richards, 1965):  (CH2O)106(NH3)16H3PO4 + 53SO42- = 106CO2 + 106H2O + 16NH3 + H3PO4 + 53S2-

! POM sulfatereduction (1st stage):
     s4_rd_PM      =th(cc(o2, ci),s_omso_o2,_ONE_,_ZERO_) &
                    *th((cc(no3, ci)+cc(no2, ci)),s_omso_no,_ONE_,_ZERO_) &
                    *K_s4_rd*5000.*cc(pon, ci)

! DOM sulfatereduction (1st stage):
     s4_rd_DM      =th(cc(o2, ci),s_omso_o2,_ONE_,_ZERO_) &
                    *th((cc(no3, ci)+cc(no2, ci)),s_omso_no,_ONE_,_ZERO_) &
                    *K_s4_rd*5000.*cc(don, ci) 

! POM sulfatereduction (2d stage): 
     s23_rd_PM     =th(cc(o2, ci),s_omso_o2,_ONE_,_ZERO_) &
                    *th((cc(no3, ci)+cc(no2, ci)),s_omso_no,_ONE_,_ZERO_) &
                    *K_s23_rd*cc(pon, ci)*cc(s2o3, ci) 

! DOM sulfatereduction (2d stage):
     s23_rd_DM     =th(cc(o2, ci),s_omso_o2,_ONE_,_ZERO_) &
                    *th((cc(no3, ci)+cc(no2, ci)),s_omso_no,_ONE_,_ZERO_) &
                    *K_s23_rd*cc(don, ci)*cc(s2o3, ci) 

     s4_rd         =s4_rd_PM  + s4_rd_DM ! in S units
     s23_rd        =s23_rd_PM + s23_rd_DM

     pp(s2o3,so4,ci)=s4_rd
     dd(so4,s2o3,ci)=pp(s2o3,so4,ci)

     pp(h2s,s2o3,ci)=s23_rd
     dd(s2o3,h2s,ci)=pp(h2s,s2o3,ci)

     DcPM_SO4      =16./53.*(s4_rd_PM+s23_rd_PM) ! in N units
     DcDM_SO4      =16./53.*(s4_rd_DM+s23_rd_DM)

! total OM mineralization ("N" units)

     AmmonPON      =DcPM_O2+DcPM_NO3+DcPM_SO4
     AmmonDON      =DcDM_O2+DcDM_NO3+DcDM_SO4

     pp(nh4,don,ci)=AmmonDON
     dd(don,nh4,ci)=pp(nh4,don,ci)

     pp(nh4,pon,ci)=AmmonPON
     dd(pon,nh4,ci)=pp(nh4,pon,ci)

! total OM mineralization ("P" units)

     PhosPOP       =AmmonPON/16.
 
     pp(po4,pop,ci)=PhosPOP
     dd(pop,po4,ci)=pp(po4,pop,ci)

     PhosDOP     =AmmonDON/16.

     pp(po4,dop,ci)=PhosDOP
     dd(dop,po4,ci)=pp(po4,dop,ci)

!---------------------------------------------------------------------------
!========N==================================================================
!---------------------------------------------------------------------------
!Autolysis: POM -> DOM :
     AutolisN      =KPD*cc(pon, ci)
     pp(don,pon,ci)=AutolisN
     dd(pon,don,ci)=pp(don,pon,ci)

!---------------------------------------- 
! Nitrogen fixation  described as appearence of NH4 available for 
! phytoplankton: N2 -> NH4 :

     SigmaN        =cc(no3, ci)+cc(no2, ci)+cc(nh4, ci)
     snf           =1./(1.+(SigmaN/cc(po4, ci)*16.)**4.)
     rnfp          =cc(po4, ci)/(cc(po4, ci)+0.3)
     Nfixation     =10.*snf*rnfp*cc(phy, ci)*KNF*LimLight*LimT*Sn

     pp(nh4,nh4,ci)=Nfixation

!---------------------------------------- 
! Anammox : NO2 + NH4 -> N2 + 2H2O :

     anammox       =th(cc(o2, ci),s_anm_o2,_ONE_,_ZERO_)*cc(no2, ci)*cc(nh4, ci)*k_annamox

     dd(nh4,nh4,ci)=anammox
     dd(no2,no2,ci)=anammox+Denitr2

!----------------------------------------
! Nitrification 1st stage: NH4+ + 1.5 O2 -> NO2- + 2H+ + H2O

     Nitrif1       =th(cc(o2, ci),s_nf1_O2,_ZERO_,_ONE_)*KN42*cc(nh4, ci)*cc(o2, ci) /(O2s_nf+cc(o2, ci))

     pp(no2,nh4,ci)=Nitrif1
     dd(nh4,no2,ci)=pp(no2,nh4,ci)
     dd(o2,nh4,ci) =1.5*Nitrif1

!----------------------------------------
! Nitrification 2d stage: NO2- + 0.5 O2 -> NO3- :

     Nitrif2       =th(cc(o2, ci),s_nf2_O2,_ZERO_,_ONE_)*KN23*cc(no2, ci)*cc(o2, ci)/(O2s_nf+cc(o2, ci))

     pp(no3,no2,ci)=Nitrif2
     dd(no2,no3,ci)=pp(no3,no2,ci)
     dd(o2,no2,ci) =0.5*Nitrif2

!---------------------------------------------------------------------------
!========S==================================================================
!--------------------------------------------------------------------------

! S0 disproportionation: 4S0 + 3H2O -> 2H2S + S2O32-+ 2H+ :

     Disprop       =0.01*cc(s0, ci)

     dd(s0,s0,ci)  =Disprop
     pp(h2S,s0,ci) =0.5*dd(s0,s0,ci)
     pp(s2o3,s0,ci) =0.5*Disprop

!----------------------------------------
! HS oxidation with O2: 2H2S + O2 -> 2S0 + 2H2O :

     hs_ox         =K_hs_ox*cc(h2s, ci)*cc(o2, ci) 

     pp(s0,h2s,ci) =hs_ox
     dd(h2s,s0,ci) =pp(s0,h2s,ci)
     dd(o2,h2s,ci) =0.5*pp(s0,h2s,ci)

!----------------------------------------
! S0 oxidation with O2: 2S0 + O2 + H2O -> S2O32- + 2H+ :

     s0_ox         =K_s0_ox*cc(s0, ci)*cc(o2, ci) 

     dd(s0,s2o3,ci)=s0_ox
     dd(o2,s0,ci)  =dd(s0,s2o3,ci)
     pp(s2o3,s0,ci)=s0_ox+0.5*Disprop

!----------------------------------------
! S2O3 oxidation with O2: S2O32- + 2O2 + 2OH- -> 2SO42- + H2O :
 
     s23_ox        =K_s23_ox*cc(s2o3, ci)*cc(o2, ci)  

     pp(so4,s2o3,ci)=s23_ox
     dd(s2o3,so4,ci)=pp(so4,s2o3,ci)
     dd(o2,s2o3,ci) =2.*pp(so4,s2o3,ci)

!----------------------------------------
! Thiodenitrification: 3H2S + 4NO3- + 6OH- -> 3SO42- + 2N2 + 6H2O

     sulfido       =KT*cc(h2s, ci)*cc(no3, ci)

     pp(so4,h2s,ci)=sulfido
     dd(h2s,so4,ci)=pp(so4,h2s,ci)
     dd(no3,no3,ci)=1.25*sulfido


!---------------------------------------------------------------------------
!========Mn==================================================================
!---------------------------------------------------------------------------
 
! Mn2 oxidation: 4Mn2+ + O2 + 4H+ -> 4Mn3+ + 2H2O :

     mn_ox         =th(cc(mn2, ci),s_mnox_mn2,_ZERO_,_ONE_)&
                    *K_mn_ox*cc(mn2, ci) *cc(o2, ci)/(cc(o2, ci)+k_mnoxO2)

     pp(mn3,mn2,ci)=mn_ox
     dd(mn2,mn3,ci)=pp(mn3,mn2,ci)
     dd(o2,mn2,ci) =0.5*pp(mn3,mn2,ci)

!----------------------------------------------------
! Mn3 oxidation: 4Mn3+ + O2 +6OH- -> 4MnO2+6H2O	:

     mn_ox2        =th(cc(mn3, ci),s_mnox_mn3,_ZERO_,_ONE_)&
                    *K_mn_ox2*cc(mn3, ci) *cc(o2, ci)/(cc(o2, ci)+k_mnoxO2)

     pp(mn4,mn3,ci)=mn_ox2
     dd(mn3,mn4,ci)=pp(mn4,mn3,ci)
     dd(o2,mn3,ci) =0.5*pp(mn4,mn3,ci)

!----------------------------------------------------
! Mn4 reduction: 2MnO2 + 7H+ + HS- -> 2Mn3+ + 4H2O + S0 :

     mn_rd         =th(cc(mn4, ci),s_mnrd_mn4,_ZERO_,_ONE_)&
                    *K_mn_rd*cc(mn4, ci)*cc(h2s, ci)/(cc(h2s, ci)+k_mnrdHS)

     pp(mn3,mn4,ci)=mn_rd
     dd(mn4,mn3,ci)=pp(mn3,mn4,ci)
     dd(h2s,mn4,ci)=0.5*pp(mn3,mn4,ci)
     pp(s0,mn3,ci) =0.5*pp(mn3,mn4,ci)

!----------------------------------------------------
! Mn3 reduction: 2Mn3+ + HS- -> 2Mn2+ + S0 + H+ :

     mn_rd2        =th(cc(mn3, ci),s_mnrd_mn3,_ZERO_,_ONE_)&
                    *K_mn_rd2*cc(mn3, ci)*cc(h2s, ci)/(cc(h2s, ci)+k_mnrdHS)

     pp(mn2,mn3,ci)=mn_rd2
     dd(mn3,mn2,ci)=pp(mn2,mn3,ci)
     dd(h2s,mn3,ci)=0.5*pp(mn2,mn3,ci)
     pp(s0,mn2,ci) =0.5*pp(mn2,mn3,ci)
!---------------------------------------------------------------------------
!========Fe==================================================================
!---------------------------------------------------------------------------

! Fe2 oxidation by O2: 4Fe2+ + O2 + 2H2O -> 4Fe3+ + 4OH :
     fe_ox         =th(cc(fe2, ci),s_feox_fe2,_ZERO_,_ONE_)*K_fe_ox*cc(o2, ci)*cc(fe2, ci) 

     pp(fe3,fe2,ci)=fe_ox
     dd(fe2,fe3,ci)=pp(fe3,fe2,ci)
     dd(o2,fe2,ci) =pp(fe3,fe2,ci)

!----------------------------------------------------
! Fe2 oxidation by Mn4: 2Fe2++MnO2+2H2O ®2FeOOH+Mn2++2H+ :
     fe_mnox       =th(cc(fe2, ci),s_feox_fe2,_ZERO_,_ONE_)*K_fe_mnox*cc(mn4, ci)*cc(fe2, ci) 

     pp(fe3,mn2,ci)=fe_mnox
     dd(fe2,mn4,ci)=pp(fe3,mn2,ci)
     dd(mn4,fe2,ci)=0.5*pp(fe3,mn2,ci)
     pp(mn2,fe3,ci)=0.5*pp(fe3,mn2,ci)

!----------------------------------------------------
! Fe3 reduction: 2FeOOH + H2S -> 2Fe2+ + S0 + 4OH- :
     fe_rd         =K_fe_rd*cc(fe3, ci)*cc(h2s, ci) 

     pp(fe2,fe3,ci)=fe_rd
     dd(fe3,fe2,ci)=pp(fe2,fe3,ci)
     dd(h2s,fe3,ci)=0.5*pp(fe2,fe3,ci)
     pp(s0,fe3,ci)=dd(h2s,fe3,ci)

!---------------------------------------------------------------------------
!========B=a=c=t(ci)============================================================
!---------------------------------------------------------------------------
!---------------------------------------- OXIC CONDITIONS
!---------------------------------------- aerobic autotrophs

     Chemos        =(Nitrif1+Nitrif2+mn_ox+fe_ox+s23_ox+s0_ox +anammox) &
                     *10.*(cc(bae, ci)+0.000001*th(cc(bae, ci),s_bac_new,_ONE_,_ZERO_)) &
		     *cc(nh4, ci)/(cc(nh4, ci)+0.05)

     pp(bae,po4,ci)=Chemos
     dd(po4,bae,ci)=Sp*pp(bae,po4,ci)
     dd(nh4,bae,ci)=Sn*pp(bae,po4,ci)


     MortBaut      =(0.1+ 0.899*th(cc(o2, ci),s_bbe_o2,_ONE_,_ZERO_)) &
                    *cc(bae, ci)*cc(bae, ci)

     dd(bae,pop,ci)=MortBaut
     pp(pop,bae,ci)=Sp*0.3*dd(bae,pop,ci)
     pp(dop,bae,ci)=Sp*0.7*dd(bae,pop,ci)
     pp(pon,bae,ci)=Sn*0.3*dd(bae,pop,ci)
     pp(don,bae,ci)=Sn*0.7*dd(bae,pop,ci)

!---------------------------------------- aerobic heterotroph

     Hetero        = th((cc(don, ci)+cc(pon, ci)),s_bhe_on,_ZERO_,_ONE_)*(DcPM_O2+DcDM_O2) &
 		    *(cc(bhe, ci)+0.000001*th(cc(bhe, ci),s_bac_new,_ONE_,_ZERO_)) &
                    *(cc(don, ci))/(cc(don, ci)+2.5)

     pp(bhe,dop,ci)=Hetero
     dd(dop,bhe,ci)=Sp*pp(bhe,dop,ci)
     dd(don,bhe,ci)=Sn*pp(bhe,dop,ci)

     MortBhet      =cc(bhe, ci)*cc(bhe, ci)*(0.001) !+ 0.899*th(cc(o2, ci),s_bbe_o2,_ONE_,_ZERO_)) &

     dd(bhe,pop,ci)=MortBhet
     pp(pop,bhe,ci)=Sp*0.3*dd(bhe,pop,ci)
     pp(dop,bhe,ci)=Sp*0.7*dd(bhe,pop,ci)
     pp(pon,bhe,ci)=Sn*0.3*dd(bhe,pop,ci)
     pp(don,bhe,ci)=Sn*0.7*dd(bhe,pop,ci)

!---------------------------------------- ANOXIC CONDITIONS
!------------------------------------- anaerobic autotrophs
     ChemosA       =(1.*mn_rd+mn_rd2+fe_rd+hs_ox +sulfido+sulfido2) &
                     *3.5*cc(baa, ci) *cc(nh4, ci)/(cc(nh4, ci)+3.5)

     pp(baa,po4,ci)=ChemosA
     dd(po4,baa,ci)=Sp*pp(baa,po4,ci)
     dd(nh4,baa,ci)=Sn*pp(baa,po4,ci)

     MortBautA      =(0.05+ 0.899*th(cc(h2s, ci),s_bba_hs,_ZERO_,_ONE_)) &
                    *cc(baa, ci)*cc(baa, ci)

     dd(baa,pop,ci)=MortBautA
     pp(pop,baa,ci)=Sp*0.3*dd(baa,pop,ci)
     pp(dop,baa,ci)=Sp*0.7*dd(baa,pop,ci)
     pp(pon,baa,ci)=Sn*0.3*dd(baa,pop,ci)
     pp(don,baa,ci)=Sn*0.7*dd(baa,pop,ci)

!---------------------------------------- anaerobic heterotroph
     HeteroA       =(0.5*(DcPM_NO3+DcDM_NO3) + DcPM_SO4+DcDM_SO4) &
 		    *(cc(bha, ci)+0.000001*th(cc(bha, ci),s_bac_new,_ONE_,_ZERO_)) &
                    *60.*cc(don, ci)/(cc(don, ci)+0.1)
  
     pp(bha,dop,ci)=HeteroA
     dd(dop,bha,ci)=Sp*pp(bha,dop,ci)
     dd(don,bha,ci)=Sn*pp(bha,dop,ci)

     MortBhetA     =0. !00001*cc(bha, ci)*cc(bha, ci) 

     dd(bha,pop,ci)=MortBhetA
     pp(pop,bha,ci)=Sp*0.3*dd(bha,pop,ci)
     pp(dop,bha,ci)=Sp*0.7*dd(bha,pop,ci)
     pp(pon,bha,ci)=Sn*0.3*dd(bha,pop,ci)
     pp(don,bha,ci)=Sn*0.7*dd(bha,pop,ci)

!---------------------------------------------------------------------------
!========P===================================================================
!---------------------------------------------------------------------------

     AutolisP      =KPD*cc(pop, ci)
     pp(dop,pop,ci)=AutolisP
     dd(pop,dop,ci)=pp(dop,pop,ci)
     pp(po4,po4,ci)=th(cc(po4, ci),s_po4_srp,_ZERO_,_ONE_)*fe_rd/2.7+(mn_ox2 + mn_rd2)/0.67
     dd(po4,po4,ci)=th(cc(po4, ci),s_po4_srp,_ZERO_,_ONE_)*(fe_ox+fe_mnox)/2.7+(mn_ox+mn_rd)/0.67

     do i=1,numc
       do j=1,numc
         pp(i,j,ci)=pp(i,j,ci)/86400. !from days^-1 to sec^-1
         dd(i,j,ci)=dd(i,j,ci)/86400. !from days^-1 to sec^-1
       enddo
     enddo

   end do
!---------------------------------------- end of main loop

!---------------------------------------- Boundary conditions, low boundary
! we use here the relaxation conditiond with relaxation time Trel
!---------------------------------------- burial into the sediments
   dd(pon,pon,1)=(-1)*w_pon*cc(pon,1)/h(1)*Bu/Trel
   dd(pop,pop,1)=(-1)*w_pop*cc(pop,1)/h(1)*Bu/Trel
   dd(s0,s0,1)  =(-1)*w_s0*cc(s0,1)/h(1)*Bu/Trel
   dd(mn4,mn4,1)=(-1)*w_mn4*cc(mn4,1)/h(1)*Bu/Trel

!---------------------------------------- upward fluxes of dissolved parameters
!---------------------------------------- independent on redox conditions
   dd(nh4,nh4,1)=(+1)*(cc(nh4,1)-b_nh4)/Trel
   dd(po4,po4,1)=(+1)*(cc(po4,1)-b_po4)/Trel


!  y=1-0.5(1-tanh(1-x^2)), dd(h2s, h2s,1) = y * (cc(h2s,1) - SedConcH2S) / Trel. 
!  where y is a dependence of H2S on O2 concentration in the near-bottom layer
!  SedConcH2S is concentration of H2S in sediment,
!  x is oxygen concentration

!---------------------------------------- upward fluxes in oxic conditions
   dd(h2s,h2s,1)=(1-0.5*(1-tanh(O2LimC-cc(o2,1))))*(cc(h2s,1)-b_h2s)/Trel
   dd(mn2,mn2,1)=(1-0.5*(1-tanh(O2LimC-cc(o2,1))))*(cc(mn2,1)-b_mn2)/Trel
   dd(fe2,fe2,1)=(1-0.5*(1-tanh(O2LimC-cc(o2,1))))*(cc(fe2,1)-b_fe2)/Trel

!---------------------------------------- downward fluxes in oxic conditions
   dd(o2,o2,1)  =(+1)*0.05*(1-tanh(O2LimC-cc(o2,1)))*(cc(o2,1)-2.)/Trel

   return

   end subroutine do_bio_rolm
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finish the bio calculations
!
! !INTERFACE:
   subroutine end_bio_rolm
!
! !DESCRIPTION:
!  Nothing done yet --- supplied for completeness.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): E.Yakushev, O.Podymov & I.Kuznetsov
!
!EOP
!-----------------------------------------------------------------------
!BOC

   return
   end subroutine end_bio_rolm
!EOC

!-----------------------------------------------------------------------
   end module bio_rolm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

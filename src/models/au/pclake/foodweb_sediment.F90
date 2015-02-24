#include "fabm_driver.h"
!BOP
!!
! !INTERFACE:
   module au_pclake_foodweb_sediment
!
! !DESCRIPTION:
!-----------------------------------------------------------------------
! Module description
!-----------------------------------------------------------------------
!  The web_sed module describes the state variables regarding zoobenthos, thus
!  sDBen,sPBen,sNBen.local processes include consumption,migration, assimilation
!  respiration(only for sDBen), excretion(only for sNBen and sPBen) and mortality.
!  This module also discribes the processes which influence the state variables registered in
!  other modules, including:(aPhytS stands for all groups of settled phytoplankton)
!  Settled phytoplankton grazed by zoobenthos: aDPhytS==>sDBent,aNPhytS==>sNBent,aPPhytS==>sPBent
!  Detritus in the sediment grazed by zoobenthos: sDDetS==>sDBent,sNDetS==>sNBent,sPDetS==>sPBent
!  Detritus morted by zoobenthos: sDDetS<==sDBent,sNDetS<==sNBent,sPDetS<==sPBent
!  Nutrients excreted by zoobenthos: sNH4S<==sNBent,sPO4S<==sPBent,sNO3W<==0.0
!  This module also provide important diagnostic variable will be used in other modules, including:
!  Sediment detritus change, tDWebDetS, used by module:auxilary
!  environmental_correction_of_fish_correction,tDEnvFiAd,used by module: foodweb_water
!  food_limitation_function_of_adult_fish,aDSatFiAd,used by module: foodweb_water
!  adult fish egestion to detritus and nutrients in the water column(through sediment top, due to 
!  adult fish is predating zoobenthos on the bottom): sDFiAd==>sDDetW,sNFiAd==>sNDetW&sNH4W,sPFiAd==>sPDetW&sPO4W
! !USES:
   use fabm_types
   use fabm_expressions
   use au_pclake_utility, ONLY:uFunTmBio
   implicit none
!  default: all is private.
   private
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_au_pclake_foodweb_sediment
!     local state variable identifers
!     id_sDBent,zoobenthos concentration in dry-weight, gDW/m**2
!     id_sPBent,zoobenthos concentration in nitrogen element, gN/m**2
!     id_sNBent,zoobenthos concentration in phosphorus element, gP/m**2
      type (type_bottom_state_variable_id) :: id_sDBent,id_sPBent,id_sNBent
!     diagnostic variables for dependencies(without output)
      type (type_horizontal_diagnostic_variable_id)       :: id_tDWebDetS
      type (type_horizontal_diagnostic_variable_id)       :: id_tDEnvFiAd ,id_aDSatFiAd
!     state dependencies identifers
      type (type_bottom_state_variable_id)            :: id_DfoodDiatS,id_DfoodGrenS,id_DfoodBlueS,id_DDetpoolS
      type (type_bottom_state_variable_id)            :: id_NfoodDiatS,id_NfoodGrenS,id_NfoodBlueS,id_NDetpoolS
      type (type_bottom_state_variable_id)            :: id_PfoodDiatS,id_PfoodGrenS,id_PfoodBlueS,id_PDetpoolS
      type (type_bottom_state_variable_id)            :: id_NH4poolS,id_NO3poolS,id_PO4poolS,id_SiDetpoolS
      type (type_state_variable_id)                   :: id_DAdFish,id_NAdFish,id_PAdFish,id_DJvFish
      type (type_state_variable_id)                   :: id_NH4poolW,id_PO4poolW,id_DDetpoolW,id_NDetpoolW,id_PDetpoolW
!     environmental dependencies
      type (type_dependency_id)                       :: id_uTm
      type ( type_horizontal_dependency_id)           :: id_aCovVeg
      type (type_horizontal_dependency_id)     :: id_sDepthW
!     Model parameters
      real(rk)           :: cDBentIn,kMigrBent,cDCarrBent,kDAssBent,hDFoodBent,fDAssBent,fDissEgesBent
      real(rk)           :: kDRespBent,kMortBent,fDissMortBent,cTmOptBent,cSigTmBent,cPDBentRef
      real(rk)           :: cNDBentRef,cSiDDiat,fDAssFiAd,cPDFishRef,cNDFishRef,fDissEgesFish
      real(rk)           :: cSigTmFish,cTmOptFish
      real(rk)           :: cRelVegFish,hDBentFiAd,kMortFiAd,kDRespFiAd,kDAssFiAd,cDCarrFish
!     nutrient ratios parameter
      real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
      real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax

   contains
!     Module procedures
      procedure :: initialize
      procedure :: do_bottom
      end type type_au_pclake_foodweb_sediment
!  private data memebers(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero=0.000000000000000000000000000000001_rk
!   Lowest state variable value for foodweb
   real(rk),parameter :: WebZero=0.0001_rk
!EOP
!-----------------------------------------------------------------------
   contains
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the passive tracer model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
! !INPUT PARAMETERS:
   class (type_au_pclake_foodweb_sediment), intent(inout),      target :: self
   integer,                          intent(in)            :: configunit
!  LOCAL VARIABLES:
   real(rk)              :: sDBent_initial
   real(rk)              :: sPBent_initial
   real(rk)              :: sNBent_initial
   real(rk)              :: cDBentIn
   real(rk)              :: kMigrBent
   real(rk)              :: cDCarrBent
   real(rk)              :: kDAssBent
   real(rk)              :: hDFoodBent
   real(rk)              :: fDAssBent
   real(rk)              :: fDissEgesBent
   real(rk)              :: kDRespBent
   real(rk)              :: kMortBent
   real(rk)              :: fDissMortBent
   real(rk)              :: cTmOptBent
   real(rk)              :: cSigTmBent
   real(rk)              :: cPDBentRef
   real(rk)              :: cNDBentRef
   real(rk)              :: cSiDDiat
   real(rk)              :: fDAssFiAd
   real(rk)              :: cPDFishRef
   real(rk)              :: cNDFishRef
   real(rk)              :: fDissEgesFish
   real(rk)              :: kDAssFiAd
   real(rk)              :: cSigTmFish
   real(rk)              :: cTmOptFish
   real(rk)              :: cRelVegFish
   real(rk)              ::  hDBentFiAd
   real(rk)              :: kMortFiAd
   real(rk)              :: kDRespFiAd
   real(rk)              :: cDCarrFish
   real(rk)              ::cNDBlueMax
   real(rk)              ::cNDBlueMin
   real(rk)              ::cNDDiatMax
   real(rk)              ::cNDDiatMin
   real(rk)              ::cNDGrenMax
   real(rk)              ::cNDGrenMin
   real(rk)              ::cPDBlueMax
   real(rk)              ::cPDBlueMin
   real(rk)              ::cPDDiatMax
   real(rk)              ::cPDDiatMin
   real(rk)              ::cPDGrenMax
   real(rk)              ::cPDGrenMin
   character(len=64)     :: diatom_as_food_DW
   character(len=64)     :: green_as_food_DW
   character(len=64)     :: blue_as_food_DW
   character(len=64)     :: detritus_DW_pool_sediment
   character(len=64)     :: diatom_as_food_N
   character(len=64)     :: green_as_food_N
   character(len=64)     :: blue_as_food_N
   character(len=64)     :: detritus_N_pool_sediment
   character(len=64)     :: diatom_as_food_P
   character(len=64)     :: green_as_food_P
   character(len=64)     :: blue_as_food_P
   character(len=64)     :: detritus_P_pool_sediment
   character(len=64)     :: detritus_Si_pool_sediment
   character(len=64)     :: NH4_pool_sediment
   character(len=64)     :: NO3_pool_sediment
   character(len=64)     :: PO4_pool_sediment
   character(len=64)     :: adult_fish_biomass
   character(len=64)     :: adult_fish_nitrogen
   character(len=64)     :: adult_fish_phosphrus
   character(len=64)     :: young_fish_biomass
   character(len=64)     :: NH4_pool_water
   character(len=64)     :: PO4_pool_water
   character(len=64)     :: DDet_pool_water
   character(len=64)     :: NDet_pool_water
   character(len=64)     :: PDet_pool_water
   character(len=64)     :: vegetation_coverage
!  create namelist
   namelist /au_pclake_foodweb_sediment/ sDBent_initial,sPBent_initial,sNBent_initial, &
                           & cDBentIn,kMigrBent,cDCarrBent,kDAssBent,hDFoodBent,fDAssBent,fDissEgesBent, &
                           & kDRespBent,kMortBent,fDissMortBent,cTmOptBent,cSigTmBent,cPDBentRef,cNDBentRef,cSiDDiat,&
                           & fDAssFiAd,cPDFishRef,cNDFishRef,fDissEgesFish,cSigTmFish,cTmOptFish,&
                           & cRelVegFish,hDBentFiAd,kMortFiAd,kDRespFiAd,kDAssFiAd,cDCarrFish, &
                           & cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin, &
                           & cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax, &
                           & diatom_as_food_DW,green_as_food_DW,blue_as_food_DW,diatom_as_food_N,green_as_food_N,&
                           & blue_as_food_N,diatom_as_food_P,green_as_food_P,blue_as_food_P,&
                           & detritus_DW_pool_sediment,detritus_P_pool_sediment,detritus_N_pool_sediment,detritus_Si_pool_sediment,&
                           & NH4_pool_sediment,NO3_pool_sediment,PO4_pool_sediment, &
                           & adult_fish_biomass,adult_fish_nitrogen,adult_fish_phosphrus,young_fish_biomass, &
                           & NH4_pool_water,PO4_pool_water,DDet_pool_water,NDet_pool_water, &
                           & PDet_pool_water,vegetation_coverage
!EOP
!-----------------------------------------------------------------------
!BOC
!  initialize the parameters
   sDBent_initial=1.0_rk
   sPBent_initial=0.1_rk
   sNBent_initial=0.01_rk
   cDBentIn=0.01_rk
   kMigrBent=0.001_rk
   cDCarrBent=10.0_rk
   kDAssBent=0.1_rk
   hDFoodBent=200.0_rk
   fDAssBent=0.3_rk
   fDissEgesBent=0.25_rk
   kDRespBent=0.005_rk
   kMortBent=0.005_rk
   fDissMortBent=0.1_rk
   cTmOptBent=25.0_rk
   cSigTmBent=16.0_rk
   cPDBentRef=0.01_rk
   cNDBentRef=0.07_rk
   cSiDDiat=0.15_rk
   fDAssFiAd=0.4_rk
   cPDFishRef=0.022_rk
   cNDFishRef=0.1_rk
   fDissEgesFish=0.25_rk
   cTmOptFish=25.0_rk
   cSigTmFish=10.0_rk
   cRelVegFish=0.009_rk
   hDBentFiAd=0.023_rk
   kMortFiAd=0.00027_rk
   kDRespFiAd=0.004_rk
   kDAssFiAd=0.06_rk
   cDCarrFish=0.14_rk
   cNDBlueMax=0.15_rk
   cNDBlueMin=0.03_rk
   cNDDiatMax=0.05_rk
   cNDDiatMin=0.01_rk
   cNDGrenMax=0.1_rk
   cNDGrenMin=0.02_rk
   cPDBlueMax=0.025_rk
   cPDBlueMin=0.0025_rk
   cPDDiatMax=0.005_rk
   cPDDiatMin=0.0005_rk
   cPDGrenMax=0.015_rk
   cPDGrenMin=0.0015_rk
   diatom_as_food_DW=''
   green_as_food_DW=''
   blue_as_food_DW=''
   detritus_DW_pool_sediment=''
   diatom_as_food_N=''
   green_as_food_N=''
   blue_as_food_N=''
   detritus_N_pool_sediment=''
   diatom_as_food_P=''
   green_as_food_P=''
   blue_as_food_P=''
   detritus_P_pool_sediment=''
   detritus_Si_pool_sediment=''
   NH4_pool_sediment=''
   NO3_pool_sediment=''
   PO4_pool_sediment=''
   adult_fish_biomass=''
   adult_fish_nitrogen=''
   adult_fish_phosphrus=''
   young_fish_biomass=''
   NH4_pool_water=''
   PO4_pool_water=''
   DDet_pool_water=''
   NDet_pool_water=''
   PDet_pool_water=''
   vegetation_coverage=''
!EOP                             
!-----------------------------------------------------------------------
!BOC                             
!  Read parameters namelist
   if (configunit>0) read(configunit,nml=au_pclake_foodweb_sediment,err=99,end=100)
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cDBentIn,'cDBentIn',default=cDBentIn)
   call self%get_parameter(self%kMigrBent,'kMigrBent',default=kMigrBent,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrBent,'cDCarrBent',default=cDCarrBent)
   call self%get_parameter(self%kDAssBent,'kDAssBent',default=kDAssBent,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDFoodBent,'hDFoodBent',default=hDFoodBent)
   call self%get_parameter(self%fDAssBent,'fDAssBent',default=fDAssBent)
   call self%get_parameter(self%fDissEgesBent,'fDissEgesBent',default=fDissEgesBent)
   call self%get_parameter(self%kDRespBent,'kDRespBent',default=kDRespBent,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortBent,'kMortBent',default=kMortBent,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortBent,'fDissMortBent',default=fDissMortBent)
   call self%get_parameter(self%cTmOptBent,'cTmOptBent',default=cTmOptBent)
   call self%get_parameter(self%cSigTmBent,'cSigTmBent',default=cSigTmBent)
   call self%get_parameter(self%cPDBentRef,'cPDBentRef',default=cPDBentRef)
   call self%get_parameter(self%cNDBentRef,'cNDBentRef',default=cNDBentRef)
   call self%get_parameter(self%cSiDDiat,'cSiDDiat',default=cSiDDiat)
   call self%get_parameter(self%fDAssFiAd,'fDAssFiAd',default=fDAssFiAd)
   call self%get_parameter(self%cPDFishRef,'cPDFishRef',default=cPDFishRef)
   call self%get_parameter(self%cNDFishRef,'cNDFishRef',default=cNDFishRef)
   call self%get_parameter(self%fDissEgesFish,'fDissEgesFish',default=fDissEgesFish)
   call self%get_parameter(self%cTmOptFish,'cTmOptFish',default=cTmOptFish)
   call self%get_parameter(self%cSigTmFish,'cSigTmFish',default=cSigTmFish)
   call self%get_parameter(self%cRelVegFish,'cRelVegFish',default=cRelVegFish)
   call self%get_parameter(self%kDAssFiAd,'kDAssFiAd',default=kDAssFiAd,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDBentFiAd,'hDBentFiAd',default=hDBentFiAd)
   call self%get_parameter(self%kDRespFiAd,'kDRespFiAd',default=kDRespFiAd,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortFiAd,'kMortFiAd',default=kMortFiAd,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrFish,'cDCarrFish',default=cDCarrFish)
   call self%get_parameter(self%cNDDiatMin,'cNDDiatMin',default=cNDDiatMin)
   call self%get_parameter(self%cPDDiatMin,'cPDDiatMin',default=cPDDiatMin)
   call self%get_parameter(self%cNDGrenMin,'cNDGrenMin',default=cNDGrenMin)
   call self%get_parameter(self%cPDGrenMin,'cPDGrenMin',default=cPDGrenMin)
   call self%get_parameter(self%cNDBlueMin,'cNDBlueMin',default=cNDBlueMin)
   call self%get_parameter(self%cPDBlueMin,'cPDBlueMin',default=cPDBlueMin)
   call self%get_parameter(self%cNDBlueMax,'cNDBlueMax',default=cNDBlueMax)
   call self%get_parameter(self%cNDDiatMax,'cNDDiatMax',default=cNDDiatMax)
   call self%get_parameter(self%cNDGrenMax,'cNDGrenMax',default=cNDGrenMax)
   call self%get_parameter(self%cPDBlueMax,'cPDBlueMax',default=cPDBlueMax)
   call self%get_parameter(self%cPDDiatMax,'cPDDiatMax',default=cPDDiatMax)
   call self%get_parameter(self%cPDGrenMax,'cPDGrenMax',default=cPDGrenMax)
   
   
   
!  Register local state variable
   call self%register_state_variable(self%id_sDBent,'sDBent','g/m**2','zoobenthos_DW',     &
                                    sDBent_initial,minimum=WebZero)
   call self%register_state_variable(self%id_sPBent,'sPBent','g/m**2','zoobenthos_P',     &
                                    sPBent_initial,minimum=WebZero)
   call self%register_state_variable(self%id_sNBent,'sNBent','g/m**2','zoobenthos_N',     &
                                    sNBent_initial,minimum=WebZero)
!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_tDWebDetS,'tDWebDetS','g/m**2/s','tDWebDetS',           &
                   output=output_none)
   call self%register_diagnostic_variable(self%id_tDEnvFiAd,'tDEnvFiAd','g/m**2','tDEnvFiAd',           &
                   output=output_none)
   call self%register_diagnostic_variable(self%id_aDSatFiAd,'aDSatFiAd','g/m**2','aDSatFiAd',           &
                   output=output_none)
   
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNBent)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPBent)
!  regirster state variables dependencies (2 steps) (the external link to sediment is unsolved)
!  step1, Register dependencies on external state variables
   call self%register_state_dependency(self%id_DfoodDiatS, 'diatom_as_food_DW','g/m**2','diatom_as_food_DW')
   call self%register_state_dependency(self%id_DfoodGrenS, 'green_as_food_DW','g/m**2','green_as_food_DW')
   call self%register_state_dependency(self%id_DfoodBlueS, 'blue_as_food_DW','g/m**2','blue_as_food_DW')
   call self%register_state_dependency(self%id_NfoodDiatS, 'diatom_as_food_N','g/m**2','diatom_as_food_N')
   call self%register_state_dependency(self%id_NfoodGrenS, 'green_as_food_N','g/m**2','green_as_food_N')
   call self%register_state_dependency(self%id_NfoodBlueS, 'blue_as_food_N','g/m**2','blue_as_food_N')
   call self%register_state_dependency(self%id_PfoodDiatS, 'diatom_as_food_P','g/m**2','diatom_as_food_P')
   call self%register_state_dependency(self%id_PfoodGrenS, 'green_as_food_P','g/m**2','green_as_food_P')
   call self%register_state_dependency(self%id_PfoodBlueS, 'blue_as_food_P','g/m**2','blue_as_food_P')
   call self%register_state_dependency(self%id_DDetpoolS, 'detritus_DW_pool_sediment','g/m**2','detritus_DW_pool_sediment')
   call self%register_state_dependency(self%id_PDetpoolS, 'detritus_P_pool_sediment','g/m**2','detritus_P_pool_sediment')
   call self%register_state_dependency(self%id_NDetpoolS, 'detritus_N_pool_sediment','g/m**2','detritus_N_pool_sediment')
   call self%register_state_dependency(self%id_SiDetpoolS, 'detritus_Si_pool_sediment','g/m**2','detritus_Si_pool_sediment')
   call self%register_state_dependency(self%id_NH4poolS, 'NH4_pool_sediment','g/m**2','NH4_pool_sediment')
   call self%register_state_dependency(self%id_NO3poolS, 'NO3_pool_sediment','g/m**2','NO3_pool_sediment')
   call self%register_state_dependency(self%id_PO4poolS, 'PO4_pool_sediment','g/m**2','PO4_pool_sediment')
   call self%register_state_dependency(self%id_DAdFish, 'adult_fish_biomass','g/m**3','adult_fish_biomass')
   call self%register_state_dependency(self%id_NAdFish, 'adult_fish_nitrogen','g/m**3','adult_fish_nitrogen')
   call self%register_state_dependency(self%id_PAdFish, 'adult_fish_phosphrus','g/m**3','adult_fish_phosphrus')
   call self%register_state_dependency(self%id_NH4poolW, 'NH4_pool_water','g/m**3','NH4_pool_water')
   call self%register_state_dependency(self%id_PO4poolW, 'PO4_pool_water','g/m**3','PO4_pool_water')
   call self%register_state_dependency(self%id_DDetpoolW, 'DDet_pool_water','g/m**3','DDet_pool_water')
   call self%register_state_dependency(self%id_NDetpoolW, 'NDet_pool_water','g/m**3','NDet_pool_water')
   call self%register_state_dependency(self%id_PDetpoolW, 'PDet_pool_water','g/m**3','PDet_pool_water')
   call self%register_state_dependency(self%id_DJvFish, 'young_fish_biomass','g/m**3','young_fish_biomass')
!  step 2, Automatically couple dependencies if target variables have been specified.
   if (diatom_as_food_DW/='')     call self%request_coupling(self%id_DfoodDiatS,diatom_as_food_DW)
   if (green_as_food_DW/='')      call self%request_coupling(self%id_DfoodGrenS,green_as_food_DW)
   if (blue_as_food_DW/='')       call self%request_coupling(self%id_DfoodBlueS,blue_as_food_DW)
   if (diatom_as_food_N/='')      call self%request_coupling(self%id_NfoodDiatS,diatom_as_food_N)
   if (green_as_food_N/='')       call self%request_coupling(self%id_NfoodGrenS,green_as_food_N)
   if (blue_as_food_N/='')        call self%request_coupling(self%id_NfoodBlueS,blue_as_food_N)
   if (diatom_as_food_P/='')      call self%request_coupling(self%id_PfoodDiatS,diatom_as_food_P)
   if (green_as_food_P/='')       call self%request_coupling(self%id_PfoodGrenS,green_as_food_P)
   if (blue_as_food_P/='')        call self%request_coupling(self%id_PfoodBlueS,blue_as_food_P)
   if (detritus_DW_pool_sediment/='')   call self%request_coupling(self%id_DDetpoolS,detritus_DW_pool_sediment)
   if (detritus_P_pool_sediment/='')    call self%request_coupling(self%id_PDetpoolS,detritus_P_pool_sediment)
   if (detritus_N_pool_sediment/='')    call self%request_coupling(self%id_NDetpoolS,detritus_N_pool_sediment)
   if (detritus_Si_pool_sediment/='')   call self%request_coupling(self%id_SiDetpoolS,detritus_Si_pool_sediment)
   if (NH4_pool_sediment/='')  call self%request_coupling(self%id_NH4poolS,NH4_pool_sediment)
   if (NO3_pool_sediment/='')  call self%request_coupling(self%id_NO3poolS,NO3_pool_sediment)
   if (PO4_pool_sediment/='')  call self%request_coupling(self%id_PO4poolS,PO4_pool_sediment)
   if (adult_fish_biomass/='')    call self%request_coupling(self%id_DAdFish,adult_fish_biomass)
   if (adult_fish_nitrogen/='')   call self%request_coupling(self%id_NAdFish,adult_fish_nitrogen)
   if (adult_fish_phosphrus/='')  call self%request_coupling(self%id_PAdFish,adult_fish_phosphrus)
   if (NH4_pool_water/='')   call self%request_coupling(self%id_NH4poolW,NH4_pool_water)
   if (PO4_pool_water/='')   call self%request_coupling(self%id_PO4poolW,PO4_pool_water)
   if (DDet_pool_water/='')  call self%request_coupling(self%id_DDetpoolW,DDet_pool_water)
   if (NDet_pool_water/='')  call self%request_coupling(self%id_NDetpoolW,NDet_pool_water)
   if (PDet_pool_water/='')  call self%request_coupling(self%id_PDetpoolW,PDet_pool_water)
   if (young_fish_biomass/='')    call self%request_coupling(self%id_DJvFish,young_fish_biomass)
!  register diagnostic dependencies, 2 steps
!  step1, Register dependencies on external diagnostic variables
   call self%register_dependency(self%id_aCovVeg, 'vegetation_coverage','--','vegetation_coverage')
!  step 2, Automatically couple dependencies if target variables have been specified.
   if (vegetation_coverage/='') call self%request_coupling(self%id_aCovVeg,vegetation_coverage)
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,standard_variables%temperature)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
   return

99 call self%fatal_error('au_pclake_foodweb_sediment_init','Error reading namelist au_pclake_foodweb_sediment')

100 call self%fatal_error('au_pclake_foodweb_sediment_init','Namelist au_pclake_foodweb_sediment was not found.')

   end subroutine initialize
!EOC
!-----------------------------------------------------------------------
!BOP
!
 !IROUTINE: 
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_au_pclake_foodweb_sediment), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
! !LOCAL VARIABLES:
!  state variables value carriers
   real(rk)                   :: sDBent,sPBent,sNBent
!  environmental dependencies carriers
   real(rk)                   :: uTm,sDepthW
!  external links variable carriers
   real(rk)         :: sDDiatS,sDGrenS,sDBlueS,sDDetS
   real(rk)         :: sNDiatS,sNGrenS,sNBlueS,sNDetS
   real(rk)         :: sPDiatS,sPGrenS,sPBlueS,sPDetS
   real(rk)         :: sDFiAd,sDFiJv
!  variables for status auxilaries
   real(rk)          :: aDPhytS,aPPhytS,aNPhytS
   real(rk)          :: rPDBent,rNDBent
   real(rk)          :: rPDBlueS,rPDGrenS,rPDDiatS,rPDDetS,rPDFoodBent
   real(rk)          :: rNDDetS,rNDDiatS,rNDGrenS,rNDBlueS,rNDFoodBent
!  variables for temperature functions
   real(rk)          :: uFunTmBent,uFunTmFish
!  variables for DW fluxes
   real(rk)          :: tDWebBent,tDMigrBent,tDAssBent
   real(rk)          :: tDRespBent,tDMortBent,aDFoodBent
   real(rk)          :: aDSatBent,tDEnvBent,ukDIncrBent,tDConsBent
   real(rk)          :: tDConsBlueBent,tDConsGrenBent,tDConsDiatBent
!  variables for P fluxes
   real(rk)          :: tPWebBent,tPMigrBent,tPAssBent
   real(rk)          :: tPExcrBent,tPMortBent,aPFoodBent,afPAssBent
   real(rk)          :: tPConsBent,tPConsDetBent,tPConsPhytBent
   real(rk)          :: tPConsDiatBent,tPConsGrenBent,tPConsBlueBent
   real(rk)          :: tDConsDetBent
!  variables for N fluxes
   real(rk)          :: tNWebBent,tNMigrBent,tNAssBent
   real(rk)          :: tNExcrBent,tNMortBent,aNFoodBent,afNAssBent
   real(rk)          :: tNConsDiatBent,tNConsGrenBent,tNConsBlueBent
   real(rk)          :: tNConsPhytBent,tNConsDetBent,tNConsBent
!  variables for exchange of NH4S
   real(rk)          :: tNWebNH4S,tNEgesBentNH4,tNEgesBent,tNMortBentNH4
!  variables for exchange of NO3S
   real(rk)          :: tNWebNO3S
!  variables for exchange of PO4S
   real(rk)          :: tPWebPO4S,tPEgesBentPO4,tPEgesBent,tPMortBentPO4
!  variables for exchange for detritus
   real(rk)          :: tDWebDetS,tDEgesBent,tNWebDetS,tNEgesBentDet
   real(rk)          :: tNMortBentDet,tPWebDetS,tPEgesBentDet,tPMortBentDet
   real(rk)          :: tSiWebDetS,tSiConsDiatBent
!  variables for exchange for diatom
   real(rk)          :: tDWebDiatS,tNWebDiatS,tPWebDiatS
!  variables for exchange for green algae
   real(rk)          :: tDWebGrenS,tNWebGrenS,tPWebGrenS
!  variables for exchange for green algae
   real(rk)          :: tDWebBlueS,tNWebBlueS,tPWebBlueS
!  adult fish assimilation
   real(rk)          :: tDEnvFiAd,tDAssFiAd,ukDIncrFiAd,aFunVegFish,aDSatFiAd
   real(rk)          :: aCovVeg
!  variables related to adult fish assimilation and consumption
   real(rk)          :: tDConsFiAd,tNConsFiAd,tPConsFiAd
   real(rk)          :: afNAssFiAd,tNAssFiAd,afPAssFiAd,tPAssFiAd
   real(rk)          :: tDEgesFiAd,tNEgesFiAd,tPEgesFiAd
   real(rk)          :: tNEgesFiAdNH4,tPEgesFiAdPO4,tNEgesFiAdDet,tPEgesFiAdDet


!  Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_
!-----------------------------------------------------------------------
!  Retrieve current (local) state variable values.
!-----------------------------------------------------------------------
   _GET_HORIZONTAL_(self%id_sDBent,sDBent)
   _GET_HORIZONTAL_(self%id_sPBent,sPBent)
   _GET_HORIZONTAL_(self%id_sNBent,sNBent)
!-----------------------------------------------------------------------
!  Retrieve dependencis value
!-----------------------------------------------------------------------
!  Retrieve state dependencie value
   _GET_HORIZONTAL_(self%id_DDetpoolS,sDDetS)
   _GET_HORIZONTAL_(self%id_NDetpoolS,sNDetS)
   _GET_HORIZONTAL_(self%id_PDetpoolS,sPDetS)
   _GET_HORIZONTAL_(self%id_DfoodDiatS,sDDiatS)
   _GET_HORIZONTAL_(self%id_DfoodGrenS,sDGrenS)
   _GET_HORIZONTAL_(self%id_DfoodBlueS,sDBlueS)
   _GET_HORIZONTAL_(self%id_NfoodDiatS,sNDiatS)
   _GET_HORIZONTAL_(self%id_NfoodGrenS,sNGrenS)
   _GET_HORIZONTAL_(self%id_NfoodBlueS,sNBlueS)
   _GET_HORIZONTAL_(self%id_PfoodDiatS,sPDiatS)
   _GET_HORIZONTAL_(self%id_PfoodGrenS,sPGrenS)
   _GET_HORIZONTAL_(self%id_PfoodBlueS,sPBlueS)

   _GET_(self%id_DAdFish,sDFiAd)
   _GET_(self%id_DJvFish,sDFiJv)

   
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
! !retrieve diagnostic denpendency
   _GET_HORIZONTAL_(self%id_aCovVeg,aCovVeg)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
   sDFiAd=sDFiAd*sDepthW
   sDFiJv=sDFiJv*sDepthW
!----------------------------------------------------------------------
!  Current local nutrients ratios in zoobenthos(check the curent state)
!----------------------------------------------------------------------
   rPDBent=sPBent/(sDBent+NearZero)
   rNDBent=sNBent/(sDBent+NearZero)
   rPDDetS=sPDetS/(sDDetS+NearZero)
   rNDDetS=sNDetS/(sDDetS+NearZero)
   rPDBlueS=sPBlueS/(sDBlueS+NearZero)
   rPDGrenS=sPGrenS/(sDGrenS+NearZero)
   rPDDiatS=sPDiatS/(sDDiatS+NearZero)
   rNDDiatS=sNBlueS/(sDBlueS+NearZero)
   rNDGrenS=sNGrenS/(sDGrenS+NearZero)
   rNDBlueS=sNDiatS/(sDDiatS+NearZero)
!    check for phosphrus nutrient ratios
   if ( rPDDiatS .GT. self%cPDDiatMax)  then
       rPDDiatS=self%cPDDiatMax
   elseif (rPDDiatS .LT. self%cPDDiatMin)  then
       rPDDiatS = self%cPDDiatMin
   else
       rPDDiatS=rPDDiatS
   endif
   
   if ( rPDBlueS .GT. self%cPDBlueMax)  then
       rPDBlueS=self%cPDBlueMax
   elseif (rPDBlueS .LT. self%cPDBlueMin)  then
       rPDBlueS = self%cPDBlueMin
   else
       rPDBlueS=rPDBlueS
   endif

   if ( rPDGrenS .GT. self%cPDGrenMax)  then
       rPDGrenS=self%cPDGrenMax
   elseif (rPDGrenS .LT. self%cPDGrenMin)  then
       rPDGrenS = self%cPDGrenMin
   else
       rPDGrenS=rPDGrenS
   endif
!   check for nitrogen nutrient ratios
   if ( rNDBlueS .GT. self%cNDBlueMax)  then
       rNDBlueS=self%cNDBlueMax
   elseif (rNDBlueS .LT. self%cNDBlueMin)  then
       rNDBlueS = self%cNDBlueMin
   else
       rNDBlueS =rNDBlueS
   endif
   
   if ( rNDDiatS .GT. self%cNDDiatMax)  then
       rNDDiatS=self%cNDDiatMax
   elseif (rNDDiatS .LT. self%cNDDiatMin)  then
       rNDDiatS = self%cNDDiatMin
   else
       rNDDiatS=rNDDiatS
   endif
   
   
   if ( rNDGrenS .GT. self%cNDGrenMax)  then
       rNDGrenS=self%cNDGrenMax
   elseif (rNDGrenS .LT. self%cNDGrenMin)  then
       rNDGrenS = self%cNDGrenMin
   else
       rNDGrenS =rNDGrenS
   endif

!  auxilaries for phytoplankton
   aDPhytS=sDDiatS+sDGrenS+sDBlueS
   aPPhytS=sPDiatS+sPGrenS+sPBlueS
   aNPhytS=sNDiatS+sNGrenS+sNBlueS
!-----------------------------------------------------------------------
!  temperature function
!-----------------------------------------------------------------------
!  temp._function_of_zoobenthos
   uFunTmFish = uFunTmBio(uTm,self%cSigTmFish,self%cTmOptFish)
   uFunTmBent = uFunTmBio(uTm,self%cSigTmBent,self%cTmOptBent)
!---------------------------------------------------------------------------
!  zoobenthos migration
!---------------------------------------------------------------------------
!  migration_flux
   tDMigrBent = self%kMigrBent *(self%cDBentIn - sDBent)
!  net_migration_flux
   tPMigrBent = self%kMigrBent *(self%cPDBentRef*self%cDBentIn - sPBent)
!  Net_migration_flux
   tNMigrBent = self%kMigrBent *(self%cNDBentRef*self%cDBentIn - sNBent)
!---------------------------------------------------------------------------
!  zoobenthos assimilation,DW
!---------------------------------------------------------------------------
!  food_for_zoobenthos
   aDFoodBent = sDDetS + aDPhytS
!  food_limitation_function_of_zoobenthos
   aDSatBent = aDFoodBent /(self%hDFoodBent + aDFoodBent)
!  intrinsic_net_increase_rate_of_zoobenthos
   ukDIncrBent = (self%kDAssBent - self%kDRespBent) * uFunTmBent - self%kMortBent
!  environmental_correction_of_zoobenthos
   tDEnvBent = max(0.0_rk,ukDIncrBent / self%cDCarrBent * sDBent*sDBent)
!  assimilation_of_zoobenthos
   tDAssBent = aDSatBent *(self%kDAssBent * uFunTmBent * sDBent - tDEnvBent)
!  consumption_of_zoobenthos
   tDConsBent = tDAssBent / self%fDAssBent
!  detritus_consumption_by_zoobenthos
   tDConsDetBent = sDDetS / aDFoodBent * tDConsBent
!  diatoms_consumption_by_zoobenthos
   tDConsDiatBent = sDDiatS / aDFoodBent * tDConsBent
!  greens_consumption_by_zoobenthos
   tDConsGrenBent = sDGrenS / aDFoodBent * tDConsBent
!  blue-greens_consumption_by_zoobenthos
   tDConsBlueBent = sDBlueS / aDFoodBent * tDConsBent
!---------------------------------------------------------------------------
!  zoobenthos assimilation,P
!---------------------------------------------------------------------------
!  food_for_zoobenthos
   aPFoodBent = sPDetS + aPPhytS
!  average_P/D_ratio_of_zoobenthos_food
   rPDFoodBent = aPFoodBent /(aDFoodBent+NearZero)
!  detrital_P_consumption_by_zoobenthos
   tPConsDetBent = rPDDetS * tDConsDetBent
!  diatom_P_consumption_by_zoobenthos
   tPConsDiatBent = rPDDiatS * tDConsDiatBent
!  greens_P_consumption_by_zoobenthos
   tPConsGrenBent = rPDGrenS * tDConsGrenBent
!  blue-greens_P_consumption_by_zoobenthos
   tPConsBlueBent = rPDBlueS * tDConsBlueBent
!  phytoplankton_P_consumption_by_zoobenthos
   tPConsPhytBent = tPConsDiatBent + tPConsGrenBent + tPConsBlueBent
!  total_P_consumption_of_zoobenthos
   tPConsBent = tPConsDetBent + tPConsPhytBent
!  P_assim._efficiency_of_zoobenthos
   afPAssBent = min(1.0_rk,self%cPDBentRef / rPDFoodBent * self%fDAssBent)
!  P_assimilation_of_zoobenthos
   tPAssBent = afPAssBent * tPConsBent
!---------------------------------------------------------------------------
!  zoobenthos assimilation,N
!---------------------------------------------------------------------------
!  food_for_zoobenthos
   aNFoodBent = sNDetS + aNPhytS
!  average_N/D_ratio_of_zoobenthos_food
   rNDFoodBent = aNFoodBent /(aDFoodBent+NearZero)
!  detrital_N_consumption_by_zoobenthos
   tNConsDetBent = rNDDetS * tDConsDetBent
!  diatom_N_consumption_by_zoobenthos
   tNConsDiatBent = rNDDiatS * tDConsDiatBent
!  greens_N_consumption_by_zoobenthos
   tNConsGrenBent = rNDGrenS * tDConsGrenBent
!  blue-greens_N_consumption_by_zoobenthos
   tNConsBlueBent = rNDBlueS * tDConsBlueBent
!  phytoplankton_N_consumption_by_zoobenthos
   tNConsPhytBent = tNConsDiatBent + tNConsGrenBent + tNConsBlueBent
!  total_N_consumption_of_zoobenthos
   tNConsBent = tNConsDetBent + tNConsPhytBent
!  N_assim._efficiency_of_zoobenthos
   afNAssBent = min(1.0_rk,self%cNDBentRef / rNDFoodBent * self%fDAssBent)
!  N_assimilation_of_zoobenthos
   tNAssBent = afNAssBent * tNConsBent
!-----------------------------------------------------------------------
!  zoobenthos respiration and excretion, DW,P and N
!-----------------------------------------------------------------------
!  respiration_of_zoobenthos
   tDRespBent = (self%cPDBentRef / rPDBent) * self%kDRespBent * uFunTmBent * sDBent
!  P_excretion_of_zoobenthos
   tPExcrBent = (rPDBent / self%cPDBentRef) * self%kDRespBent * uFunTmBent * sPBent
!  N_excretion_of_zoobenthos
   tNExcrBent = (rNDBent / self%cNDBentRef) * self%kDRespBent * uFunTmBent * sNBent
!-----------------------------------------------------------------------
!  zoobenthos mortality
!-----------------------------------------------------------------------
!  zoobenthos_mortality_incl._environmental_correction
   tDMortBent = self%kMortBent*sDBent +(1.0_rk - aDSatBent) * tDEnvBent
!  mortality_of_zoobenthos
   tPMortBent = rPDBent * tDMortBent
!  mortality_of_zoobenthos
   tNMortBent = rNDBent * tDMortBent
!-----------------------------------------------------------------------
!  zoobenthos egestion
!-----------------------------------------------------------------------
!  egestion_of_zoobenthos_DW
   tDEgesBent = tDConsBent - tDAssBent
!  egestion_of_zoobenthos_N
   tNEgesBent = tNConsBent - tNAssBent
!-----------------------------------------------------------------------
!  adult fish assimilation_DW
!-----------------------------------------------------------------------
!  vegetation_dependence_of_fish_feeding
   aFunVegFish = max(0.0_rk,1.0_rk - self%cRelVegFish * aCovVeg)
!   for first time step check out
!   aFunVegFish = max(0.0_rk,1.0_rk - self%cRelVegFish * 0.2_rk)
!  food_limitation_function_of_adult_fish
   aDSatFiAd = (aFunVegFish * sDBent) *(aFunVegFish * sDBent) /(self%hDBentFiAd * &
   &self%hDBentFiAd + (aFunVegFish * sDBent) *(aFunVegFish * sDBent))
!  intrinsic_net_increase_rate_of_fish
   ukDIncrFiAd = (self%kDAssFiAd - self%kDRespFiAd) * uFunTmFish - self%kMortFiAd
!  environmental_correction_of_fish,in concentration
   tDEnvFiAd = max(0.0_rk,ukDIncrFiAd /(self%cDCarrFish - sDFiJv) * sDFiAd*sDFiAd)
!  assimilation_of_fish
   tDAssFiAd = aDSatFiAd *(self%kDAssFiAd * uFunTmFish * sDFiAd - tDEnvFiAd)
!-----------------------------------------------------------------------
!  adult fish assimilation_P
!-----------------------------------------------------------------------
!  zoobenthos_consumption_of_fish
   tDConsFiAd = tDAssFiAd / self%fDAssFiAd
!  (zoobenthos)_P_consumption_by_FiAd
   tPConsFiAd = rPDBent * tDConsFiAd
!  P_assim._efficiency_of_FiAd
   afPAssFiAd = min(1.0_rk,self%cPDFishRef / rPDBent * self%fDAssFiAd)
!  P_assimilation_of_FiAd
   tPAssFiAd = afPAssFiAd * tPConsFiAd
!-----------------------------------------------------------------------
!  adult fish assimilation_N
!-----------------------------------------------------------------------
!   (zoobenthos)_N_consumption_by_FiAd
    tNConsFiAd = rNDBent * tDConsFiAd
!   N_assim._efficiency_of_FiAd
    afNAssFiAd = min(1.0_rk,self%cNDFishRef / rNDBent * self%fDAssFiAd)
!   N_assimilation_of_FiAd
    tNAssFiAd = afNAssFiAd * tNConsFiAd
!-----------------------------------------------------------------------
!  external state variables change due to adult fish(egestion, from sediment top)
!-----------------------------------------------------------------------
!  egestion_of_fish,adult fish
   tDEgesFiAd = tDConsFiAd - tDAssFiAd
!  egestion_of_FiAd
   tPEgesFiAd = tPConsFiAd - tPAssFiAd
!  egestion_of_FiAd
   tNEgesFiAd = tNConsFiAd - tNAssFiAd
!  NH4_egestion_of_adult_fish
   tNEgesFiAdNH4 = self%fDissEgesFish * tNEgesFiAd
!  SRP_egestion_of_adult_fish
   tPEgesFiAdPO4 = self%fDissEgesFish * tPEgesFiAd
!  detrital_N_egestion_of_adult_fish
   tNEgesFiAdDet = tNEgesFiAd - tNEgesFiAdNH4
!  detrital_P_egestion_of_adult_fish
   tPEgesFiAdDet = tPEgesFiAd - tPEgesFiAdPO4
!-----------------------------------------------------------------------
!  total flux of web change to state variables
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_Zoobenthos
   tDWebBent = tDMigrBent + tDAssBent - tDRespBent - tDMortBent - tDConsFiAd
!  total_foodweb_flux_of_P_in_Zoobenthos
   tPWebBent = tPMigrBent + tPAssBent - tPExcrBent - tPMortBent - tPConsFiAd
!  total_foodweb_flux_of_N_in_Zoobenthos
   tNWebBent = tNMigrBent + tNAssBent - tNExcrBent - tNMortBent - tNConsFiAd
!=======================================================================
!  foodweb part relating to other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in sediment
!-----------------------------------------------------------------------
!  part_of_died_zoobenthos_N_becoming_ammonium-N
   tNMortBentNH4 = self%fDissMortBent*tNMortBent
!  NH4_egestion_of_zoobenthos
   tNEgesBentNH4 = self%fDissEgesBent * tNEgesBent
!  total_foodweb_flux_of_N_in_Pore_water_ammonium_in_lake_sediment
   tNWebNH4S = tNExcrBent + tNEgesBentNH4 + tNMortBentNH4
!-----------------------------------------------------------------------
!  Update NO3 in sediment
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_N_in_Pore_water_nitrate_in_lake_sediment
   tNWebNO3S = 0.0_rk
!-----------------------------------------------------------------------
!  Update PO4 in sediment
!-----------------------------------------------------------------------
!  part_of_died_zoobenthos_P_becoming_dissolved_P
   tPMortBentPO4 = self%fDissMortBent * tPMortBent
!  egestion_of_zoobenthos
   tPEgesBent = tPConsBent - tPAssBent
!  SRP_egestion_of_zoobenthos
   tPEgesBentPO4 = self%fDissEgesBent * tPEgesBent
!  total_foodweb_flux_of_P_in_Pore_water_P_in_lake_sediment
   tPWebPO4S = tPExcrBent + tPEgesBentPO4 + tPMortBentPO4
!-----------------------------------------------------------------------
!  Update detritus in sediment(DW,N,P)
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_Sediment_detritus_in_lake
   tDWebDetS = - tDConsDetBent + tDEgesBent + tDMortBent
!  part_of_died_zoobenthos_N_becoming_detrital_N
   tNMortBentDet = (1.0_rk-self%fDissMortBent)*tNMortBent
!  detrital_N_egestion_of_zoobenthos
   tNEgesBentDet = (1.0_rk - self%fDissEgesBent) * tNEgesBent
!  total_foodweb_flux_of_N_in_Sediment_N_in_lake_sediment
   tNWebDetS = - tNConsDetBent + tNEgesBentDet + tNMortBentDet
!  part_of_died_zoobenthos_P_becoming_detrital_P
   tPMortBentDet = (1.0_rk-self%fDissMortBent)*tPMortBent
!  detrital_P_egestion_of_zoobenthos
   tPEgesBentDet = (1.0_rk - self%fDissEgesBent) * tPEgesBent
!  total_foodweb_flux_of_P_in_Sediment_P_in_lake
   tPWebDetS = - tPConsDetBent + tPEgesBentDet + tPMortBentDet
!  diatom_consumption_by_zoobenthos
   tSiConsDiatBent = self%cSiDDiat * tDConsDiatBent
!  total_foodweb_flux_of_silica_in_sediment_detritus
   tSiWebDetS = tSiConsDiatBent
!-----------------------------------------------------------------------
!  Update diatom in sediment(DW,N,P)
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_sediment_diatoms_in_lake
   tDWebDiatS = - tDConsDiatBent
!  total_foodweb_flux_of_N_in_sediment_diatoms_in_lake
   tNWebDiatS = - tNConsDiatBent
!  total_foodweb_flux_of_P_in_sediment_diatoms_in_lake
   tPWebDiatS = - tPConsDiatBent 
!-----------------------------------------------------------------------
!  Update green algae in sediment(DW,N,P)
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_sediment_greens_in_lake
   tDWebGrenS = - tDConsGrenBent
!  total_foodweb_flux_of_N_in_sediment_greens_in_lake
   tNWebGrenS = - tNConsGrenBent
!  total_foodweb_flux_of_P_in_sediment_greens_in_lake
   tPWebGrenS = - tPConsGrenBent
!-----------------------------------------------------------------------
!  Update blue algae in sediment(DW,N,P)
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_sediment_blue-greens_in_lake
   tDWebBlueS = - tDConsBlueBent 
!  total_foodweb_flux_of_N_in_sediment_blue-greens_in_lake
   tNWebBlueS = - tNConsBlueBent 
!  total_foodweb_flux_of_P_in_sediment_blue-greens_in_lake
   tPWebBlueS = - tPConsBlueBent
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
!  update state variables
   _SET_ODE_BEN_(self%id_sDBent,tDWebBent)
   _SET_ODE_BEN_(self%id_sPBent,tPWebBent)
   _SET_ODE_BEN_(self%id_sNBent,tNWebBent)
!-----------------------------------------------------------------------
!  Update external links
!-----------------------------------------------------------------------
!  update abiotic variables in sediment
   _SET_ODE_BEN_(self%id_NH4poolS,tNWebNH4S)
   _SET_ODE_BEN_(self%id_NO3poolS,tNWebNO3S)
   _SET_ODE_BEN_(self%id_PO4poolS,tPWebPO4S)
   _SET_ODE_BEN_(self%id_DDetpoolS,tDWebDetS)
   _SET_ODE_BEN_(self%id_NDetpoolS,tNWebDetS)
   _SET_ODE_BEN_(self%id_PDetpoolS,tPWebDetS)
   _SET_ODE_BEN_(self%id_SiDetpoolS,tSiWebDetS)
!  update phytoplanktons in sediment
   _SET_ODE_BEN_(self%id_DfoodDiatS,tDWebDiatS)
   _SET_ODE_BEN_(self%id_NfoodDiatS,tNWebDiatS)
   _SET_ODE_BEN_(self%id_PfoodDiatS,tPWebDiatS)
   _SET_ODE_BEN_(self%id_DfoodGrenS,tDWebGrenS)
   _SET_ODE_BEN_(self%id_NfoodGrenS,tNWebGrenS)
   _SET_ODE_BEN_(self%id_PfoodGrenS,tPWebGrenS)
   _SET_ODE_BEN_(self%id_DfoodBlueS,tDWebBlueS)
   _SET_ODE_BEN_(self%id_NfoodBlueS,tNWebBlueS)
   _SET_ODE_BEN_(self%id_PfoodBlueS,tPWebBlueS)
!  update abiotic water section variable
   _SET_BOTTOM_EXCHANGE_(self%id_DAdFish,tDAssFiAd)
   _SET_BOTTOM_EXCHANGE_(self%id_NAdFish,tNAssFiAd)
   _SET_BOTTOM_EXCHANGE_(self%id_PAdFish,tPAssFiAd)
   _SET_BOTTOM_EXCHANGE_(self%id_NH4poolW,tNEgesFiAdNH4)
   _SET_BOTTOM_EXCHANGE_(self%id_PO4poolW,tPEgesFiAdPO4)
   _SET_BOTTOM_EXCHANGE_(self%id_DDetpoolW,tDEgesFiAd)
   _SET_BOTTOM_EXCHANGE_(self%id_NDetpoolW,tNEgesFiAdDet)
   _SET_BOTTOM_EXCHANGE_(self%id_PDetpoolW,tPEgesFiAdDet)
!-----------------------------------------------------------------------
!  output diagnostic variables for external links
!-----------------------------------------------------------------------
!  Export diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDWebDetS,tDWebDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDEnvFiAd,tDEnvFiAd)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aDSatFiAd,aDSatFiAd)
   


      _FABM_HORIZONTAL_LOOP_END_
   end subroutine do_bottom
! Spatial loop end
!
!EOC
!-----------------------------------------------------------------------
   end module au_pclake_foodweb_sediment
!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

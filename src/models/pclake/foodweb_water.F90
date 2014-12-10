#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_foodweb_water
!
! !DESCRIPTION:
!  The pclake_foodweb_water module describes the state variables regarding zooplanktons
!  white fish(juvenile and adult)and piscivorous fish. 
!  Therefore, the state variables and their related its local processes are:
!  sDZoo,sPZoo,sNZoo,(zooplankton in dry-weight, nitrogen element and phosphorus element respectively)
!  units: gDW/m**3, gP/m**3,gP/m**3
!  local processes:assimilation,respiration(only for sDZoo),excretion(only for sNZoo,sPZoo),
!  mortality and consumption by juvenile fish.
!  sDFiJv,sPFiJv,sNFiJv,(young fish in dry-weight, nitrogen element and phosphorus element respectively)
!  units: gDW/m**3, gP/m**3,gP/m**3
!  local processes: migration,reproduction(+,part of adult fish became young fish),assimilation(predation
!  on zooplankton),respiration(only for sDFiJv),exretion(only for sPFiJv,sNFiJv),mortality,consumption
!  by piscivorous fish and aging(-,part of young fish become adult fish)
!  sDFiAd,sPFiAd,sNFiAd,(adult fish in dry-weight, nitrogen element and phosphorus element respectively)
!  units:gDW/m**3, gP/m**3,gP/m**3
!  local processes: migration,reproduction(-,part of adult fish became young fish),respiration(only for
!  sDFiAd),exretion(only for sPFiAd and sNFiAd),mortality, consumption by piscivorous fish ,harvest
!  and aging(+,part of young fish become adult fish).(!! Notice the assimilation of adult fish is in the 
!   benthic module, where adult fish predating the zoobenthos.)
!  sDPisc,gDW/m**3)(piscivorous fish will have fixed N/D and P/D ratio,so the other two will be diagnostic variables)
!  local processes: migration,assimilation,respiration,mortality and harvest
!  This module also discribes the processes which influence the state variables registered in
!  other modules, including: (aPhytW stands for all groups of phytoplankton)
!  Zooplankton grazing on phytoplankton in water column: aDPhytW==>aDZoo,aNPhytW==>aNZoo,aPPhytW==>aPZoo
!  Zooplankton grazing on detritus in the water column: sDDetW==>aDZoo,sNDetW==>aNZoo,sPDetW==>aPZoo
!  Detritus morted by zooplankton,white fish and piscivorous fish, and egested by all fish:
!  sDDetW<==sDZoo&sDFiJv&sDFiAd&sDPisc,sNDetW<==sNZoo&sNFiJv&sNFiAd&aNPisc,sPDetW<==sPZoo&sPFiJv&sPFiAd&aPPisc
!  nutrients excreted by zooplankton, all fish and dissolved part of mortality and egestion:
!   sNH4W<==sNZoo&sNFiJv&sNFiAd&aNPisc,sPO4<==sPZoo&sPFiJv&sPFiAd&aPPisc,sNO3W<==0.0
! !USES:
   use fabm_types
   use pclake_utility, ONLY:uFunTmBio

   implicit none

!  default: all is private.
      private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_foodweb_water
!     local state variable identifers
!     id_sDZoo,zooplankton concentration in dry-weight, gDW/m**3
!     id_sPZoo,zooplankton concentration in nitrogen element, gN/m**3
!     id_sNZoo,zooplankton concentration in phosphorus element, gP/m**3
!     id_sDFiJv,juvenile white fish concentration in dry-weight, gDW/m**3
!     id_sPFiJv,juvenile white fish concentration in nitrogen element, gN/m**3
!     id_sNFiJv,juvenile white fish concentration in phosphorus element, gP/m**3
!     id_sDFiAd,adult white fish concentration in dry-weight, gDW/m**3
!     id_sPFiAd,adult white fish concentration in nitrogen element, gN/m**3
!     id_sNFiAd,adult white fish concentration in phosphorus element, gP/m**3
!     id_sDPisc,piscivorous fish concentration in dry-weight, gDW/m**3
      type (type_state_variable_id)            :: id_sDZoo,id_sPZoo,id_sNZoo
      type (type_state_variable_id)            :: id_sDFiJv,id_sPFiJv,id_sNFiJv
      type (type_state_variable_id)            :: id_sDFiAd,id_sPFiAd,id_sNFiAd,id_sDPisc
!     diagnostic variables for local output
!     id_aNPisc, piscivorous fish concentration in nitrogen element, gN/m**3
!     id_aPPisc, piscivorous fish concentration in phosphorus element, gP/m**3
      type (type_diagnostic_variable_id)       :: id_aNPisc,id_aPPisc
!     diagnostic variables for dependencies(without output)
!     state dependencies identifers
      type (type_state_variable_id)            :: id_DfoodDiat,id_DfoodGren,id_DfoodBlue,id_DDetpoolW
      type (type_state_variable_id)            :: id_NfoodDiat,id_NfoodGren,id_NfoodBlue,id_NDetpoolW
      type (type_state_variable_id)            :: id_PfoodDiat,id_PfoodGren,id_PfoodBlue,id_PDetpoolW,id_SiDetpoolW
      type (type_state_variable_id)            :: id_NH4poolW,id_NO3poolW,id_PO4poolW
!     environmental dependencies
      type (type_dependency_id)                :: id_uTm ,id_dz
      type (type_horizontal_dependency_id)     :: id_sDepthW
      type (type_global_dependency_id)         :: id_Day
!     diagnostic dependencies
      type (type_horizontal_dependency_id)                :: id_aDSubVeg,id_tDEnvFiAd,id_aDSatFiAd
!!    Model parameters
!     parameter for temperature
      real(rk)      :: cSigTmZoo,cTmOptZoo
!     parameters for zooplankton
      real(rk)      :: cDCarrZoo,kMortZoo,kDRespZoo,cPrefDiat
      real(rk)      :: cPrefGren,cPrefBlue,cPrefDet,hFilt
      real(rk)      :: fDAssZoo,cFiltMax
      real(rk)      :: cPDZooRef,cNDZooRef
!     parameters for fish
      real(rk)      :: kMigrFish,cDFiJvIn,cDFiAdIn
      real(rk)      :: cDPiscIn,kMigrPisc,fDBone
      real(rk)      :: fPBone,cDCarrFish,fDissEgesFish,fDissMortFish
      real(rk)      :: cTmOptFish,cSigTmFish,cDayReprFish,fReprFish
      real(rk)      :: fAgeFish,kDAssFiJv,hDZooFiJv,fDAssFiJv
      real(rk)      :: kDRespFiJv,kMortFiJv
      real(rk)      :: kDRespFiAd,kMortFiAd,cDCarrPiscMax,cDCarrPiscMin
      real(rk)      :: cDCarrPiscBare,cDPhraMinPisc,cCovVegMin
      real(rk)      :: cRelPhraPisc,cRelVegPisc,kDAssPisc,hDVegPisc
      real(rk)      :: hDFishPisc,fDAssPisc,fDissEgesPisc,kDRespPisc
      real(rk)      :: kMortPisc,fDissMortPisc,cTmOptPisc,cSigTmPisc
      real(rk)      :: cPDFishRef,cNDFishRef,cPDPisc,cNDPisc
      real(rk)      :: fDissEgesZoo,fDissMortZoo,cSiDDiat
!     nutrient ratios parameter
      real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
      real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax
!  harvest paramter: ,kHarvFishWin,kHarvFishSum,kHarvPiscSum,kHarvPiscWin
      
   contains
   
!     Model procedures
      procedure :: initialize
      procedure :: do
      !procedure :: do_bottom
      procedure :: get_light_extinction

   end type type_pclake_foodweb_water
   
!  private data memebers(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk)           :: Pi=3.14159265358979_rk
!   Lowest state variable value for foodweb
   real(rk),parameter :: WebZero=0.0001_rk
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_pclake_foodweb_water), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!  LOCAL VARIABLES:
   real(rk)              :: sDZoo_initial
   real(rk)              :: sPZoo_initial
   real(rk)              :: sNZoo_initial
   real(rk)              :: sDFiJv_initial
   real(rk)              :: sPFiJv_initial
   real(rk)              :: sNFiJv_initial
   real(rk)              :: sDFiAd_initial
   real(rk)              :: sPFiAd_initial
   real(rk)              :: sNFiAd_initial
   real(rk)              :: sDPisc_initial
   real(rk)              :: cSigTmZoo
   real(rk)              :: cTmOptZoo
   real(rk)              :: cDCarrZoo
   real(rk)              :: kMortZoo
   real(rk)              :: kDRespZoo
   real(rk)              :: cPrefDiat
   real(rk)              :: cPrefGren
   real(rk)              :: cPrefBlue
   real(rk)              :: cPrefDet
   real(rk)              :: hFilt
   real(rk)              :: fDAssZoo
   real(rk)              :: cFiltMax
   real(rk)              :: cPDZooRef
   real(rk)              :: cNDZooRef
   real(rk)              :: kMigrFish
   real(rk)              :: cDFiJvIn
   real(rk)              :: cDFiAdIn
   real(rk)              :: cDPiscIn
   real(rk)              :: kMigrPisc
   real(rk)              :: fDBone
   real(rk)              :: fPBone
   real(rk)              :: cDCarrFish
   real(rk)              :: fDissEgesFish
   real(rk)              :: fDissMortFish
   real(rk)              :: cTmOptFish
   real(rk)              :: cSigTmFish
   real(rk)              :: cDayReprFish
   real(rk)              :: fReprFish
   real(rk)              :: fAgeFish
   real(rk)              :: kDAssFiJv
   real(rk)              :: hDZooFiJv
   real(rk)              :: fDAssFiJv
   real(rk)              :: kDRespFiJv
   real(rk)              :: kMortFiJv
   real(rk)              :: kDRespFiAd
   real(rk)              :: kMortFiAd
   real(rk)              :: cDCarrPiscMax
   real(rk)              :: cDCarrPiscMin
   real(rk)              :: cDCarrPiscBare
   real(rk)              :: cDPhraMinPisc
   real(rk)              :: cCovVegMin
   real(rk)              :: cRelPhraPisc
   real(rk)              :: cRelVegPisc
   real(rk)              :: kDAssPisc
   real(rk)              :: hDVegPisc
   real(rk)              :: hDFishPisc
   real(rk)              :: fDAssPisc
   real(rk)              :: fDissEgesPisc
   real(rk)              :: kDRespPisc
   real(rk)              :: kMortPisc
   real(rk)              :: fDissMortPisc
   real(rk)              :: cTmOptPisc
   real(rk)              :: cSigTmPisc
   real(rk)              :: cPDFishRef
   real(rk)              :: cNDFishRef
   real(rk)              :: cPDPisc
   real(rk)              :: cNDPisc
   real(rk)              :: fDissEgesZoo
   real(rk)              :: fDissMortZoo
   real(rk)              :: cSiDDiat
   real(rk)              :: cNDBlueMax
   real(rk)              :: cNDBlueMin
   real(rk)              :: cNDDiatMax
   real(rk)              :: cNDDiatMin
   real(rk)              :: cNDGrenMax
   real(rk)              :: cNDGrenMin
   real(rk)              :: cPDBlueMax
   real(rk)              :: cPDBlueMin
   real(rk)              :: cPDDiatMax
   real(rk)              :: cPDDiatMin
   real(rk)              :: cPDGrenMax
   real(rk)              :: cPDGrenMin
   character(len=64)     :: diatom_as_food_DW
   character(len=64)     :: green_as_food_DW
   character(len=64)     :: blue_as_food_DW
   character(len=64)     :: detritus_DW_pool_water
   character(len=64)     :: diatom_as_food_N
   character(len=64)     :: green_as_food_N
   character(len=64)     :: blue_as_food_N
   character(len=64)     :: detritus_N_pool_water
   character(len=64)     :: diatom_as_food_P
   character(len=64)     :: green_as_food_P
   character(len=64)     :: blue_as_food_P
   character(len=64)     :: detritus_P_pool_water
   character(len=64)     :: detritus_Si_pool_water
   character(len=64)     :: NH4_pool_water
   character(len=64)     :: NO3_pool_water
   character(len=64)     :: PO4_pool_water
   character(len=64)     :: submerged_vegetation
   character(len=64)     :: env_correction_adfish
   character(len=64)     :: food_limit_function_adfish
!  create namelist
   namelist /pclake_foodweb_water/ sDZoo_initial,sPZoo_initial,sNZoo_initial,sDFiJv_initial,sPFiJv_initial,sNFiJv_initial,& 
                              & sDFiAd_initial,sPFiAd_initial,sNFiAd_initial,sDPisc_initial,cSigTmZoo,cTmOptZoo,&
                              & cDCarrZoo,kMortZoo,kDRespZoo,cPrefDiat,cPrefGren,cPrefBlue,cPrefDet,hFilt,fDAssZoo,cFiltMax, &
                              & cPDZooRef,cNDZooRef,kMigrFish,cDFiJvIn,cDFiAdIn,cDPiscIn,kMigrPisc, & 
                              & fDBone,fPBone,cDCarrFish,fDissEgesFish,fDissMortFish,cTmOptFish, & 
                              & cSigTmFish,cDayReprFish,fReprFish,fAgeFish,kDAssFiJv,hDZooFiJv,fDAssFiJv,&
                              & kDRespFiJv,kMortFiJv,kDRespFiAd,kMortFiAd,cDCarrPiscMax, &
                              & cDCarrPiscMin,cDCarrPiscBare,cDPhraMinPisc,cCovVegMin,cRelPhraPisc,cRelVegPisc,kDAssPisc,&
                              & hDVegPisc,hDFishPisc,fDAssPisc,fDissEgesPisc,kDRespPisc,kMortPisc,fDissMortPisc,cTmOptPisc,cSigTmPisc,&
                              & cPDFishRef,cNDFishRef,cPDPisc,cNDPisc,fDissEgesZoo,fDissMortZoo,cSiDDiat, &
                              & cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin, &
                              & cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax, &
                              & diatom_as_food_DW,green_as_food_DW,blue_as_food_DW,detritus_DW_pool_water,diatom_as_food_N,green_as_food_N,&
                              & blue_as_food_N,detritus_N_pool_water,diatom_as_food_P,green_as_food_P,blue_as_food_P,detritus_P_pool_water,detritus_Si_pool_water,&
                              & NH4_pool_water,NO3_pool_water,PO4_pool_water,submerged_vegetation,env_correction_adfish,food_limit_function_adfish
!EOP
!-----------------------------------------------------------------------
!BOC
!  initialize the parameters
   sDZoo_initial=0.05_rk
   sPZoo_initial=0.0005_rk
   sNZoo_initial=0.0035_rk
   sDFiJv_initial=1.0_rk
   sPFiJv_initial=1.0_rk
   sNFiJv_initial=1.0_rk
   sDFiAd_initial=1.0_rk
   sPFiAd_initial=1.0_rk
   sNFiAd_initial=1.0_rk
   sDPisc_initial=1.0_rk
   cSigTmZoo=13_rk 
   cTmOptZoo=25_rk
   cDCarrZoo=25_rk
   kMortZoo=0.04_rk
   kDRespZoo=0.15_rk
   cPrefDiat=0.75_rk
   cPrefGren=0.75_rk
   cPrefBlue=0.125_rk
   cPrefDet=0.25_rk
   hFilt=1.0_rk
   fDAssZoo =0.35_rk
   cFiltMax =4.5_rk
   cPDZooRef =0.01_rk
   cNDZooRef=0.07_rk
   kMigrFish=0.001_rk
   cDFiJvIn=0.005_rk
   cDFiAdIn=0.005_rk
   cDPiscIn=0.001_rk
   kMigrPisc=0.001_rk
   fDBone=0.35_rk
   fPBone=0.5_rk
   cDCarrFish=15.0_rk
   fDissEgesFish=0.25_rk
   fDissMortFish=0.1_rk
   cTmOptFish=25.0_rk
   cSigTmFish=10.0_rk
   cDayReprFish=120.0_rk
   fReprFish=0.02_rk
   fAgeFish=0.5_rk
   kDAssFiJv=0.12_rk
   hDZooFiJv=1.25_rk
   fDAssFiJv=0.4_rk
   kDRespFiJv=0.01_rk
   kMortFiJv=0.00137_rk
   kDRespFiAd=0.004_rk
   kMortFiAd=0.00027_rk
   cDCarrPiscMax=1.2_rk
   cDCarrPiscMin=0.1_rk
   cDCarrPiscBare=0.1_rk
   cDPhraMinPisc=50.0_rk
   cCovVegMin=40.0_rk
   cRelPhraPisc=0.075_rk
   cRelVegPisc=0.03_rk
   kDAssPisc=0.025_rk
   hDVegPisc=5.0_rk
   hDFishPisc=1.0_rk
   fDAssPisc=0.4_rk
   fDissEgesPisc=0.25_rk
   kDRespPisc=0.005_rk
   kMortPisc=0.00027_rk
   fDissMortPisc=0.1_rk
   cTmOptPisc=25.0_rk
   cSigTmPisc=10.0_rk
   cPDFishRef=0.022_rk
   cNDFishRef=0.1_rk
   cPDPisc=0.022_rk
   cNDPisc=0.1_rk
   fDissEgesZoo =0.25_rk
   fDissMortZoo = 0.1_rk
   cSiDDiat=0.15_rk
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
   detritus_DW_pool_water=''
   diatom_as_food_N=''
   green_as_food_N=''
   blue_as_food_N=''
   detritus_N_pool_water=''
   diatom_as_food_P=''
   green_as_food_P=''
   blue_as_food_P=''
   detritus_P_pool_water=''
   detritus_Si_pool_water=''
   NH4_pool_water=''
   NO3_pool_water=''
   PO4_pool_water=''
   submerged_vegetation=''
   env_correction_adfish=''
   food_limit_function_adfish=''
!EOP                             
!-----------------------------------------------------------------------
!BOC                             
!  Read parameters namelist
   if (configunit>0) read(configunit,nml=pclake_foodweb_water,err=99,end=100)
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cSigTmZoo,'cSigTmZoo',default=cSigTmZoo)
   call self%get_parameter(self%cTmOptZoo,'cTmOptZoo',default=cTmOptZoo)
   call self%get_parameter(self%cDCarrZoo,'cDCarrZoo',default=cDCarrZoo)
   call self%get_parameter(self%kMortZoo,'kMortZoo',default=kMortZoo,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespZoo,'kDRespZoo',default=kDRespZoo,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cPrefDiat,'cPrefDiat',default=cPrefDiat)
   call self%get_parameter(self%cPrefGren,'cPrefGren',default=cPrefGren)
   call self%get_parameter(self%cPrefBlue,'cPrefBlue',default=cPrefBlue)
   call self%get_parameter(self%cPrefDet,'cPrefDet',default=cPrefDet)
   call self%get_parameter(self%hFilt,'hFilt',default=hFilt)
   call self%get_parameter(self%fDAssZoo,'fDAssZoo',default=fDAssZoo)
   call self%get_parameter(self%cFiltMax,'cFiltMax',default=cFiltMax,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cPDZooRef,'cPDZooRef',default=cPDZooRef)
   call self%get_parameter(self%cNDZooRef,'cNDZooRef',default=cNDZooRef)
   call self%get_parameter(self%kMigrFish,'kMigrFish',default=kMigrFish,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDFiJvIn,'cDFiJvIn',default=cDFiJvIn)
   call self%get_parameter(self%cDFiAdIn,'cDFiAdIn',default=cDFiAdIn)
   call self%get_parameter(self%cDPiscIn,'cDPiscIn',default=cDPiscIn)
   call self%get_parameter(self%kMigrPisc,'kMigrPisc',default=kMigrPisc,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDBone,'fDBone',default=fDBone)
   call self%get_parameter(self%fPBone,'fPBone',default=fPBone)
   call self%get_parameter(self%cDCarrFish,'cDCarrFish',default=cDCarrFish)
   call self%get_parameter(self%fDissEgesFish,'fDissEgesFish',default=fDissEgesFish)
   call self%get_parameter(self%fDissMortFish,'fDissMortFish',default=fDissMortFish)
   call self%get_parameter(self%cTmOptFish,'cTmOptFish',default=cTmOptFish)
   call self%get_parameter(self%cSigTmFish,'cSigTmFish',default=cSigTmFish)
   call self%get_parameter(self%cDayReprFish,'cDayReprFish',default=cDayReprFish)
   call self%get_parameter(self%fReprFish,'fReprFish',default=fReprFish)
   call self%get_parameter(self%fAgeFish,'fAgeFish',default=fAgeFish)
   call self%get_parameter(self%kDAssFiJv,'kDAssFiJv',default=kDAssFiJv,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDZooFiJv,'hDZooFiJv',default=hDZooFiJv)
   call self%get_parameter(self%fDAssFiJv,'fDAssFiJv',default=fDAssFiJv)
   call self%get_parameter(self%kDRespFiJv,'kDRespFiJv',default=kDRespFiJv,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortFiJv,'kMortFiJv',default=kMortFiJv,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespFiAd,'kDRespFiAd',default=kDRespFiAd,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortFiAd,'kMortFiAd',default=kMortFiAd,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrPiscMax,'cDCarrPiscMax',default=cDCarrPiscMax)
   call self%get_parameter(self%cDCarrPiscMin,'cDCarrPiscMin',default=cDCarrPiscMin)
   call self%get_parameter(self%cDCarrPiscBare,'cDCarrPiscBare',default=cDCarrPiscBare)
   call self%get_parameter(self%cDPhraMinPisc,'cDPhraMinPisc',default=cDPhraMinPisc)
   call self%get_parameter(self%cCovVegMin,'cCovVegMin',default=cCovVegMin)
   call self%get_parameter(self%cRelPhraPisc,'cRelPhraPisc',default=cRelPhraPisc)
   call self%get_parameter(self%cRelVegPisc,'cRelVegPisc',default=cRelVegPisc)
   call self%get_parameter(self%kDAssPisc,'kDAssPisc',default=kDAssPisc,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDVegPisc,'hDVegPisc',default=hDVegPisc)
   call self%get_parameter(self%hDFishPisc,'hDFishPisc',default=hDFishPisc)
   call self%get_parameter(self%fDAssPisc,'fDAssPisc',default=fDAssPisc)
   call self%get_parameter(self%fDissEgesPisc,'fDissEgesPisc',default=fDissEgesPisc)
   call self%get_parameter(self%kDRespPisc,'kDRespPisc',default=kDRespPisc,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortPisc,'kMortPisc',default=kMortPisc,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortPisc,'fDissMortPisc',default=fDissMortPisc)
   call self%get_parameter(self%cTmOptPisc,'cTmOptPisc',default=cTmOptPisc)
   call self%get_parameter(self%cSigTmPisc,'cSigTmPisc',default=cSigTmPisc)
   call self%get_parameter(self%cPDFishRef,'cPDFishRef',default=cPDFishRef)
   call self%get_parameter(self%cNDFishRef,'cNDFishRef',default=cNDFishRef)
   call self%get_parameter(self%cPDPisc,'cPDPisc',default=cPDPisc)
   call self%get_parameter(self%cNDPisc,'cNDPisc',default=cNDPisc)
   call self%get_parameter(self%fDissEgesZoo,'fDissEgesZoo',default=fDissEgesZoo)
   call self%get_parameter(self%fDissMortZoo,'fDissMortZoo',default=fDissMortZoo)
   call self%get_parameter(self%cSiDDiat,'cSiDDiat',default=cSiDDiat)
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
   call self%register_state_variable(self%id_sDZoo,'sDZoo','g/m**3','zooplankton DW in water',     &
                                    sDZoo_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPZoo,'sPZoo','g/m**3','zooplankton P in water',     &
                                    sPZoo_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNZoo,'sNZoo','g/m**3','zooplankton N in water',     &
                                    sNZoo_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sDFiJv,'sDFiJv','g/m**3','juvenile fish DW in water',     &
                                    sDZoo_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPFiJv,'sPFiJv','g/m**3','juvenile fish P in water',     &
                                    sPFiJv_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNFiJv,'sNFiJv','g/m**3','juvenile fish N in water',     &
                                    sNFiJv_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sDFiAd,'sDFiAd','g/m**3','adult fish DW in water',     &
                                    sDFiAd_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPFiAd,'sPFiAd','g/m**3','adult fish P in water',     &
                                    sPFiAd_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNFiAd,'sNFiAd','g/m**3','adult fish N in water',     &
                                    sNFiAd_initial,minimum=WebZero,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sDPisc,'sDPisc','g/m**3','piscivorous fish DW in water',     &
                                    sDPisc_initial,minimum=WebZero,no_river_dilution=.TRUE.)
!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_aNPisc,'aNPisc','g/m**3','aNPisc',           &
                  output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_aPPisc,'aPPisc','g/m**3','aPPisc',           &
                  output=output_time_step_averaged)

!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNZoo)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNFiJv)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNFiAd)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_aNPisc)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPZoo)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPFiJv)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPFiAd)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_aPPisc)
!  regirster state variables dependencies (2 steps) (the external link to sediment is unsolved)
!  step1, Register dependencies on external state variables
   call self%register_state_dependency(self%id_DfoodDiat, 'diatom_as_food_DW','g/m**3','diatom_as_food_DW')
   call self%register_state_dependency(self%id_DfoodGren, 'green_as_food_DW','g/m**3','green_as_food_DW')
   call self%register_state_dependency(self%id_DfoodBlue, 'blue_as_food_DW','g/m**3','blue_as_food_DW')
   call self%register_state_dependency(self%id_NfoodDiat, 'diatom_as_food_N','g/m**3','diatom_as_food_N')
   call self%register_state_dependency(self%id_NfoodGren, 'green_as_food_N','g/m**3','green_as_food_N')
   call self%register_state_dependency(self%id_NfoodBlue, 'blue_as_food_N','g/m**3','blue_as_food_N')
   call self%register_state_dependency(self%id_PfoodDiat, 'diatom_as_food_P','g/m**3','diatom_as_food_P')
   call self%register_state_dependency(self%id_PfoodGren, 'green_as_food_P','g/m**3','green_as_food_P')
   call self%register_state_dependency(self%id_PfoodBlue, 'blue_as_food_P','g/m**3','blue_as_food_P')
   call self%register_state_dependency(self%id_DDetpoolW, 'detritus_DW_pool_water','g/m**3','detritus_DW_pool_water')
   call self%register_state_dependency(self%id_NDetpoolW, 'detritus_N_pool_water','g/m**3','detritus_N_pool_water')
   call self%register_state_dependency(self%id_PDetpoolW, 'detritus_P_pool_water','g/m**3','detritus_P_pool_water')
   call self%register_state_dependency(self%id_SiDetpoolW, 'detritus_Si_pool_water','g/m**3','detritus_Si_pool_water')
   call self%register_state_dependency(self%id_NH4poolW, 'NH4_pool_water','g/m**3','NH4_pool_water')
   call self%register_state_dependency(self%id_NO3poolW, 'NO3_pool_water','g/m**3','NO3_pool_water')
   call self%register_state_dependency(self%id_PO4poolW, 'PO4_pool_water','g/m**3','PO4_pool_water')
!  step 2, Automatically couple dependencies if target variables have been specified.
   if (diatom_as_food_DW/='') call self%request_coupling(self%id_DfoodDiat,diatom_as_food_DW)
   if (green_as_food_DW/='') call self%request_coupling(self%id_DfoodGren,green_as_food_DW)
   if (blue_as_food_DW/='') call self%request_coupling(self%id_DfoodBlue,blue_as_food_DW)
   if (diatom_as_food_N/='') call self%request_coupling(self%id_NfoodDiat,diatom_as_food_N)
   if (green_as_food_N/='') call self%request_coupling(self%id_NfoodGren,green_as_food_N)
   if (blue_as_food_N/='') call self%request_coupling(self%id_NfoodBlue,blue_as_food_N)
   if (diatom_as_food_P/='') call self%request_coupling(self%id_PfoodDiat,diatom_as_food_P)
   if (green_as_food_P/='') call self%request_coupling(self%id_PfoodGren,green_as_food_P)
   if (blue_as_food_P/='') call self%request_coupling(self%id_PfoodBlue,blue_as_food_P)
   if (detritus_DW_pool_water/='') call self%request_coupling(self%id_DDetpoolW,detritus_DW_pool_water)
   if (detritus_N_pool_water/='') call self%request_coupling(self%id_NDetpoolW,detritus_N_pool_water)
   if (detritus_P_pool_water/='') call self%request_coupling(self%id_PDetpoolW,detritus_P_pool_water)
   if (detritus_Si_pool_water/='') call self%request_coupling(self%id_SiDetpoolW,detritus_Si_pool_water)
   if (NH4_pool_water/='') call self%request_coupling(self%id_NH4poolW,NH4_pool_water)
   if (NO3_pool_water/='') call self%request_coupling(self%id_NO3poolW,NO3_pool_water)
   if (PO4_pool_water/='') call self%request_coupling(self%id_PO4poolW,PO4_pool_water)
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,standard_variables%temperature)
   call self%register_dependency(self%id_Day,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_dz,standard_variables%cell_thickness)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
!  register diagnostic dependencies
!  step1, Register dependencies on external diagnostic variables
   call self%register_dependency(self%id_tDEnvFiAd, 'env_correction_adfish','--','env_correction_adfish')
   call self%register_dependency(self%id_aDSubVeg, 'submerged_vegetation','g/m**2','submerged_vegetation')
   call self%register_dependency(self%id_aDSatFiAd, 'food_limit_function_adfish','--','food_limit_function_adfish')
   
!  step 2, Automatically couple dependencies if target variables have been specified.
   if (env_correction_adfish/='') call self%request_coupling(self%id_tDEnvFiAd,env_correction_adfish)
   if (submerged_vegetation/='') call self%request_coupling(self%id_aDSubVeg,submerged_vegetation)
   if (food_limit_function_adfish/='') call self%request_coupling(self%id_aDSatFiAd,food_limit_function_adfish)

   return


99  call self%fatal_error('pclake_foodweb_water_init','Error reading namelist pclake_foodweb_water')

100 call self%fatal_error('pclake_foodweb_water_init','Namelist pclake_foodweb_water was not found.')


   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP

! !IROUTINE: 
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
!  INPUT PARAMETERS:
   class (type_pclake_foodweb_water), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_
!  LOCAL VARIABLES:
!  Carriers for environment dependencies
   real(rk)     :: uTm,Day,sDepthW,dz
!  carriers for local state variables
   real(rk)      :: sDZoo,sNZoo,sPZoo
   real(rk)      :: sDFiJv,sPFiJv,sNFiJv
   real(rk)      :: sDFiAd,sPFiAd,sNFiAd,sDPisc
!  carriers for exteral link state variables
   real(rk)      :: sDDiatW,sDGrenW,sDBlueW,sDDetW
   real(rk)      :: sNDiatW,sNGrenW,sNBlueW,sNDetW
   real(rk)      :: sPDiatW,sPGrenW,sPBlueW,sPDetW
!  carriers for external link diagnostic variables
   real(rk)      :: aDSubVeg,tDEnvFiAd,aDSatFiAd
!  nutrient ratios variables
   real(rk)      :: rNDDiatW,rNDGrenW,rNDBlueW,rNDDetW
   real(rk)      :: rPDDiatW,rPDGrenW,rPDBlueW,rPDDetW
   real(rk)      :: rPDZoo,rNDZoo
   real(rk)      :: rPDFiJv,rNDFiJv,rPDFiAd,rNDFiAd
!  status auxilaries
   real(rk)      :: aDFish,aPFish,aNFish
!  variables for temperature function
   real(rk)      :: uFunTmZoo,uFunTmFish,uFunTmPisc
!  variables of Zooplankton DW fluxes
   real(rk)      :: aDSatZoo,ukDAssTmZoo,wDEnvZoo,ukDIncrZoo
   real(rk)      :: ukDRespTmZoo,oDFoodZoo
   real(rk)      :: oDOMW,oDPhytW
   real(rk)      :: aCorDRespZoo
!  variables for N assimilation
   real(rk)      :: afNAssZoo,wNConsDetZoo,wDConsDetZoo,wDConsZoo
   real(rk)      :: wNConsZoo,wNConsPhytZoo,wNConsDiatZoo, wNConsGrenZoo,wNConsBlueZoo
   real(rk)      :: wDConsDiatZoo,wDConsGrenZoo,wDConsBlueZoo
   real(rk)      :: rNDFoodZoo,oNFoodZoo
!  variables for P assimilation
   real(rk)      :: rPDFoodZoo,oPFoodZoo,afPAssZoo,wPConsZoo
   real(rk)      :: wPConsPhytZoo,wPConsDetZoo,wPConsDiatZoo
   real(rk)      :: wPConsGrenZoo,wPConsBlueZoo
!  variables for Zooplankton dw  fluxes
   real(rk)      :: wDWebZoo,wDAssZoo,wDRespZoo,wDMortZoo
!  variables for Zooplankton N fluxes
   real(rk)      :: wNWebZoo,wNAssZoo,wNExcrZoo,wNMortZoo
!  variables for Zooplankton P fluxes
   real(rk)      ::wPWebZoo,wPAssZoo,wPExcrZoo,wPMortZoo
!  PCLake_Osis, /m^2
   real(rk)      :: tDReprFish,tDAgeFish
   real(rk)      :: tPReprFish,tPAgeFish
   real(rk)      :: tNReprFish,tNAgeFish
!  variables for young fish flux_DW
   real(rk)     :: wDWebFiJv
!  PCLake_Osis, /m^2
   real(rk)     :: tDWebFiJv,tDMigrFiJv
   real(rk)     :: tDAssFiJv,tDRespFiJv,tDMortFiJv,tDConsFiJvPisc
   real(rk)     :: aDSatFiJv,tDEnvFiJv,ukDIncrFiJv,tDConsFiJv
!  variables for young fish flux_P
   real(rk)     :: wPWebFiJv
!  PCLake_Osis, /m^2
   real(rk)     :: tPWebFiJv,tPMigrFiJv,tPAssFiJv
   real(rk)     :: tPExcrFiJv,tPMortFiJv,tPConsFiJvPisc
   real(rk)     :: afPAssFiJv,tPConsFiJv,afNAssFiJv,tNConsFiJv
!  variables for young fish flux_N
   real(rk)     :: wNWebFiJv
!  PCLake_Osis, /m^2
   real(rk)     :: tNWebFiJv,tNMigrFiJv,tNAssFiJv
   real(rk)     :: tNExcrFiJv,tNMortFiJv,tNConsFiJvPisc
!  variables for adult fish flux_DW
   real(rk)     :: wDWebFiAd
!  PCLake_Osis, /m^2
   real(rk)     :: tDWebFiAd,tDMigrFiAd,tDRespFiAd,tDMortFiAd
   real(rk)     :: tDConsFiAdPisc
!  variables for adult fish flux_P
   real(rk)     :: wPWebFiAd
!  PCLake_Osis, /m^2
   real(rk)     :: tPWebFiAd,tPMigrFiAd,tPExcrFiAd,tPMortFiAd
   real(rk)     :: tPConsFiAdPisc 
!  assimilation
!  variables for adult fish flux_N
   real(rk)     :: wNWebFiAd
!   ,wNMigrFiAd,wNExcrFiAd,wNMortFiAd
!   real(rk)     :: wNConsFiAdPisc
!  PCLake_Osis, /m^2
   real(rk)     :: tNWebFiAd,tNMigrFiAd,tNExcrFiAd,tNMortFiAd
   real(rk)     :: tNConsFiAdPisc
   real(rk)      :: wDWebPisc 
!  PCLake_Osis, /m^2
   real(rk)     :: tDConsPisc,tDAssPisc,aDSatPisc,aFunVegPisc
   real(rk)     :: tDEnvPisc,akDIncrPisc,aDCarrPisc
   real(rk)     :: tDMigrPisc,tDRespPisc,tDMortPisc
   real(rk)     :: tDWebPisc
!  PCLake_Osis, /m^2
   real(rk)     :: aPPisc,tPConsPisc,rPDFoodPisc,afPAssPisc,tPAssPisc
   real(rk)     :: tPEgesPisc,tPExcrPisc,tPMortPisc,tPMigrPisc
!  PCLake_Osis, /m^2
   real(rk)     :: aNPisc,tNConsPisc,rNDFoodPisc,afNAssPisc,tNAssPisc
   real(rk)     :: tNEgesPisc,tNExcrPisc,tNMortPisc,tNMigrPisc
!  PCLake_Osis, /m^2
   real(rk)     :: wNWebNH4W,wNEgesZooNH4,wNEgesZoo,wNMortZooNH4,tNEgesFiJvNH4
   real(rk)     :: tNEgesFiJv,tNMortFishNH4,tNMortFishBot
   real(rk)     :: tNMortFish,tNEgesPiscNH4,tNMortPiscNH4,tNMortPiscBot
!  variables for exchange of NH4
   real(rk)     :: wNWebNO3W
!  PCLake_Osis, /m^2
   real(rk)     :: wPWebPO4W,wPEgesZooPO4,wPEgesZoo,wPMortZooPO4,tPEgesFiJvPO4
   real(rk)     :: tPEgesFiJv,tPMortFish,tPMortFishBot
   real(rk)     :: tPMortFishPO4,tPEgesPiscPO4,tPMortPiscPO4,tPMortPiscBot
!  PCLake_Osis, /m^2
   real(rk)     :: wDWebDetW,wDEgesZoo,tDEgesFiJv, tDMortFishDet
   real(rk)     :: tDMortFish,tDMortFishBot,tDEgesPisc,tDMortPiscDet,tDMortPiscBot
!  PCLake_Osis, /m^2
   real(rk)     :: wNWebDetW,wNEgesZooDet,wNMortZooDet,tNEgesFiJvDet,tNMortFishDet
   real(rk)     :: tNEgesPiscDet,tNMortPiscDet
!  variables for exchange of detritus P
!  PCLake_Osis, /m^2
   real(rk)     :: wPWebDetW,wPEgesZooDet,wPMortZooDet,tPEgesFiJvDet,tPMortFishDet
   real(rk)     :: tPEgesPiscDet,tPMortPiscDet
!  variables for exchange of detritus Si
   real(rk)     :: wSiWebDetW,wSiConsDiatZoo
!  variables for exchange of diatoms
   real(rk)     :: wDWebDiatW,wNWebDiatW,wPWebDiatW
!  variables for exchange of green algae
   real(rk)     :: wDWebGrenW,wNWebGrenW,wPWebGrenW
!  variables for exchange of green algae
   real(rk)     :: wDWebBlueW,wNWebBlueW,wPWebBlueW
!  adult fish assimilation
   real(rk)     :: ukDIncrFiAd

   
   
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_(self%id_sDZoo,sDZoo)
   _GET_(self%id_sNZoo,sNZoo)
   _GET_(self%id_sPZoo,sPZoo)
   _GET_(self%id_sDFiJv,sDFiJv)
   _GET_(self%id_sPFiJv,sPFiJv)
   _GET_(self%id_sNFiJv,sNFiJv)
   _GET_(self%id_sDFiAd,sDFiAd)
   _GET_(self%id_sPFiAd,sPFiAd)
   _GET_(self%id_sNFiAd,sNFiAd)
   _GET_(self%id_sDPisc,sDPisc)
!-----------------------------------------------------------------------
!  Retrieve dependencis value
!-----------------------------------------------------------------------
!  Retrieve state dependencie value
   _GET_(self%id_DfoodDiat,sDDiatW)
   _GET_(self%id_DfoodGren,sDGrenW)
   _GET_(self%id_DfoodBlue,sDBlueW)
   _GET_(self%id_DDetpoolW,sDDetW)
   _GET_(self%id_NfoodDiat,sNDiatW)
   _GET_(self%id_NfoodGren,sNGrenW)
   _GET_(self%id_NfoodBlue,sNBlueW)
   _GET_(self%id_NDetpoolW,sNDetW)
   _GET_(self%id_PfoodDiat,sPDiatW)
   _GET_(self%id_PfoodGren,sPGrenW)
   _GET_(self%id_PfoodBlue,sPBlueW)
   _GET_(self%id_PDetpoolW,sPDetW)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_GLOBAL_(self%id_Day,Day)
!    _GET_(self%id_dz,dz)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
! !retrieve diagnostic denpendency
   _GET_HORIZONTAL_(self%id_aDSubVeg,aDSubVeg)
   _GET_HORIZONTAL_(self%id_tDEnvFiAd,tDEnvFiAd)
   _GET_HORIZONTAL_(self%id_aDSatFiAd,aDSatFiAd)
!  convert FISH concentration to areal units
   sDFiJv=sDFiJv*sDepthW
   sPFiJv=sPFiJv*sDepthW
   sNFiJv=sNFiJv*sDepthW
   sDFiAd=sDFiAd*sDepthW
   sPFiAd=sPFiAd*sDepthW
   sNFiAd=sNFiAd*sDepthW
   sDPisc=sDPisc*sDepthW
!-------------------------------------------------------------------------
!  The orders for the processes. We try to orgnize the order
!  from zooplankton, then young fish, then adult fish, at last piscivorous
!  fish. And each group finish all their process before going to the next
!  groups. But there a special sections on zooplankton comsumption by 
!  young fish, due to the variables dependent on each other.Young fish 
!  starts with assimilation of DW after whole zooplankton group is finished
!  to provide wDAssFiJv for young fish's predation of zooplankton(DW,N,P).
!  The later process(predation) provide variable wDConsFiJv for fish
!  fish assimilation of N,P
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!  Current local nutrients ratios in zooplankton(check the curent state)
!-------------------------------------------------------------------------
!  P/C_ratio of food
   rPDDiatW=sPDiatW/(sDDiatw+NearZero)
   rPDGrenW=sPGrenW/(sDGrenW+NearZero)
   rPDBlueW=sPBlueW/(sDBlueW+NearZero)
   rPDDetW=sPDetW/(sDDetW+NearZero)
!  N/C_ratio of food
   rNDDiatW=sNDiatW/(sDDiatw+NearZero)
   rNDGrenW=sNGrenW/(sDGrenW+NearZero)
   rNDBlueW=sNBlueW/(sDBlueW+NearZero)
   rNDDetW=sNDetW/(sDDetW+NearZero)
   if ( rPDDiatW .GT. self%cPDDiatMax)  then
       rPDDiatW=self%cPDDiatMax
   elseif (rPDDiatW .LT. self%cPDDiatMin)  then
       rPDDiatW = self%cPDDiatMin
   else
       rPDDiatW=rPDDiatW
   endif
   
   if ( rPDBlueW .GT. self%cPDBlueMax)  then
       rPDBlueW=self%cPDBlueMax
   elseif (rPDBlueW .LT. self%cPDBlueMin)  then
       rPDBlueW = self%cPDBlueMin
   else
       rPDBlueW=rPDBlueW
   endif

   if ( rPDGrenW .GT. self%cPDGrenMax)  then
       rPDGrenW=self%cPDGrenMax
   elseif (rPDGrenW .LT. self%cPDGrenMin)  then
       rPDGrenW = self%cPDGrenMin
   else
       rPDGrenW=rPDGrenW
   endif
!   check for nitrogen nutrient ratios before calculating droop equations
   if ( rNDBlueW .GT. self%cNDBlueMax)  then
       rNDBlueW=self%cNDBlueMax
   elseif (rNDBlueW .LT. self%cNDBlueMin)  then
       rNDBlueW = self%cNDBlueMin
   else
       rNDBlueW =rNDBlueW
   endif
   
   if ( rNDDiatW .GT. self%cNDDiatMax)  then
       rNDDiatW=self%cNDDiatMax
   elseif (rNDDiatW .LT. self%cNDDiatMin)  then
       rNDDiatW = self%cNDDiatMin
   else
       rNDDiatW=rNDDiatW
   endif
   
   
   if ( rNDGrenW .GT. self%cNDGrenMax)  then
       rNDGrenW=self%cNDGrenMax
   elseif (rNDGrenW .LT. self%cNDGrenMin)  then
       rNDGrenW = self%cNDGrenMin
   else
       rNDGrenW =rNDGrenW
   endif
!  P/D_ratio_herb.zooplankton
   rPDZoo = sPZoo /(sDZoo+NearZero)
!  N/C_ratio_herb.zooplankton
   rNDZoo = sNZoo/(sDZoo+NearZero)
!  P/D_ratio_of_young_fish
   rPDFiJv = sPFiJv /(sDFiJv+NearZero)
!  P/D_ratio_of_adult_fish
   rPDFiAd = sPFiAd /(sDFiAd+NearZero)
!  N/D_ratio_of_young_fish
   rNDFiJv = sNFiJv /(sDFiJv+NearZero)
!  N/D_ratio_of_adult_fish
   rNDFiAd = sNFiAd /(sDFiAd+NearZero)
!-----------------------------------------------------------------------
!  status auxilaries---auxilaries for describing the current status, 
!  usually derivitives of state variables
!-----------------------------------------------------------------------
!  total_fish_biomass
   aDFish = sDFiJv + sDFiAd
!  total_fish_biomass
   aPFish = sPFiJv + sPFiAd
!  total_fish_biomass
   aNFish = sNFiJv + sNFiAd
!-----------------------------------------------------------------------
!  temperature function
!-----------------------------------------------------------------------
!  temp._function_of_zooplankton
   uFunTmZoo =uFunTmBio(uTm,self%cSigTmZoo,self%cTmOptZoo)
!  temp._function_of_fish
   uFunTmFish = uFunTmBio(uTm,self%cSigTmFish,self%cTmOptFish)
!  temp._function_of_Pisc
   uFunTmPisc = uFunTmBio(uTm,self%cSigTmPisc,self%cTmOptPisc)
!-----------------------------------------------------------------------
!  zooplankton assimilation DW
!-----------------------------------------------------------------------
!  organic_seston
   oDPhytW= sDDiatW+sDGrenW+sDBlueW
   oDOMW = sDDetW + oDPhytW
!  food_for_zooplankton
   oDFoodZoo = self%cPrefDiat * sDDiatW + self%cPrefGren * sDGrenW + self%cPrefBlue * sDBlueW + self%cPrefDet * sDDetW
!  max._assimilation_rate_of_zooplankton,temp._corrected
   ukDAssTmZoo = self%fDAssZoo * self%cFiltMax * uFunTmZoo * self%hFilt
!  food_saturation_function_of_zooplankton
   aDSatZoo = oDFoodZoo /(self%hFilt + oDOMW)
!  respiration_constant_of_zooplankton
   ukDRespTmZoo = self%kDRespZoo * uFunTmZoo
!  intrinsic_rate_of_increase_of_zooplankton
   ukDIncrZoo = ukDAssTmZoo - ukDRespTmZoo - self%kMortZoo
!  environmental_correction_of_zooplankton
   wDEnvZoo = max(0.0_rk,ukDIncrZoo / self%cDCarrZoo * sDZoo*sDZoo)
!  assimilation_of_zooplankton
   wDAssZoo = aDSatZoo *(ukDAssTmZoo * sDZoo - wDEnvZoo)
!-----------------------------------------------------------------------
!  zooplankton assimilation N
!-----------------------------------------------------------------------
!  Zooplankton_food
   oNFoodZoo = self%cPrefDiat*sNDiatW + self%cPrefGren*sNGrenW + &
               & self%cPrefBlue*sNBlueW + self%cPrefDet*sNDetW
!  consumption_of_zooplankton
   wDConsZoo = wDAssZoo / self%fDAssZoo
!  N/C_ratio_of_zooplankton_food
   rNDFoodZoo = oNFoodZoo /(oDFoodZoo+NearZero)
!  DW_diatoms_consumption_by_zooplankton
   wDConsDiatZoo = self%cPrefDiat*sDDiatW / oDFoodZoo * wDConsZoo
!  DW_greens_consumption_by_zooplankton
   wDConsGrenZoo = self%cPrefGren*sDGrenW / oDFoodZoo * wDConsZoo
!  DW_blue-greens_consumption_by_zooplankton
   wDConsBlueZoo = self%cPrefBlue*sDBlueW / oDFoodZoo * wDConsZoo
!  N_diatom_consumption_by_zoopl.
   wNConsDiatZoo = rNDDiatW*wDConsDiatZoo
!  N_green_consumption_by_zoopl.
   wNConsGrenZoo = rNDGrenW*wDConsGrenZoo
!  N_bluegreen_consumption_by_zoopl.
   wNConsBlueZoo = rNDBlueW*wDConsBlueZoo
!  total_N_phytoplankton_consumption_by_zoopl.
   wNConsPhytZoo = wNConsDiatZoo + wNConsGrenZoo + wNConsBlueZoo
!  DW_detritus_consumption_by_zooplankton
   wDConsDetZoo = self%cPrefDet*sDDetW / oDFoodZoo * wDConsZoo
!  consumption_of_detrital_N
   wNConsDetZoo = rNDDetW*wDConsDetZoo
!  total_N_consumption
   wNConsZoo = wNConsPhytZoo + wNConsDetZoo
!  N_assimilation_efficiency_of_herbivores
   afNAssZoo = min(1.0_rk,self%cNDZooRef / rNDFoodZoo * self%fDAssZoo)
!  assimilation_by_herbivores
   wNAssZoo = afNAssZoo*wNConsZoo
!-----------------------------------------------------------------------
!  zooplankton assimilation P
!-----------------------------------------------------------------------
!  Zooplankton_food
   oPFoodZoo = self%cPrefDiat*sPDiatW + self%cPrefGren*sPGrenW &
               & + self%cPrefBlue*sPBlueW + self%cPrefDet*sPDetW
!  P/D_ratio_of_zooplankton_food
   rPDFoodZoo = oPFoodZoo /(oDFoodZoo+NearZero)
!  P_diatom_consumption_by_zoopl.
   wPConsDiatZoo = rPDDiatW * wDConsDiatZoo
!  P_green_consumption_by_zoopl.
   wPConsGrenZoo = rPDGrenW * wDConsGrenZoo
!  P_bluegreen_consumption_by_zoopl.
   wPConsBlueZoo = rPDBlueW * wDConsBlueZoo
!  total_P_phytoplankton_consumption_by_zoopl.
   wPConsPhytZoo = wPConsDiatZoo + wPConsGrenZoo + wPConsBlueZoo
!  consumption_of_detrital_P
   wPConsDetZoo = rPDDetW * wDConsDetZoo
!  total_P_consumption
   wPConsZoo = wPConsPhytZoo + wPConsDetZoo
!  P_assimilation_efficiency_of_herbivores
   afPAssZoo = min(1.0_rk,self%cPDZooRef / rPDFoodZoo * self%fDAssZoo)
!  assimilation_by_herbivores
   wPAssZoo = afPAssZoo * wPConsZoo
!-----------------------------------------------------------------------
!  zooplankton respiration and excretion
!-----------------------------------------------------------------------
!  corr._factor_of_zoopl._respiration_for_P_and_N_content
   aCorDRespZoo = max(self%cPDZooRef / rPDZoo,self%cNDZooRef / rNDZoo)
!  zoopl._respiration_DW
   wDRespZoo = aCorDRespZoo * self%kDRespZoo * uFunTmZoo * sDZoo
!  N_excretion
   wNExcrZoo = rNDZoo / self%cNDZooRef * self%kDRespZoo * uFunTmZoo*sNZoo
!  P_excretion
   wPExcrZoo = rPDZoo / self%cPDZooRef * self%kDRespZoo * uFunTmZoo*sPZoo
!-----------------------------------------------------------------------
!  zooplankton mortality
!-----------------------------------------------------------------------
!  zoopl._mortality,incl._environmental_correction_DW
   wDMortZoo = self%kMortZoo * sDZoo +(1.0_rk - aDSatZoo) * wDEnvZoo
!  zoopl._mortality_N
   wNMortZoo = rNDZoo*wDMortZoo
!  zoopl._mortality_P
   wPMortZoo = rPDZoo * wDMortZoo
!-----------------------------------------------------------------------
!  zooplankton egestion
!-----------------------------------------------------------------------
!  egestion_of_zooplankton
   wDEgesZoo = wDConsZoo - wDAssZoo
!  N_egestion
   wNEgesZoo = wNConsZoo - wNAssZoo
!  P_egestion
    wPEgesZoo = wPConsZoo - wPAssZoo
!-----------------------------------------------------------------------
!  young fish assimilation_DW
!-----------------------------------------------------------------------
!  food_limitation_function_of_young_fish
!   aDSatFiJv = sDZoo *sDZoo /(self%hDZooFiJv * self%hDZooFiJv &
!   &+sDZoo * sDZoo )
!  this need to be modified and re-calibrate of the parameter when transfer to 1d
   aDSatFiJv = (sDZoo * sDepthW) *(sDZoo * sDepthW) /(self%hDZooFiJv * &
   &self%hDZooFiJv + (sDZoo * sDepthW) *(sDZoo * sDepthW))
!  intrinsic_net_increase_rate_of_fish
   ukDIncrFiJv = (self%kDAssFiJv - self%kDRespFiJv) * uFunTmFish - self%kMortFiJv
!  environmental_correction_of_fish
   tDEnvFiJv = max(0.0,ukDIncrFiJv /(self%cDCarrFish - sDFiAd) * sDFiJv*sDFiJv)
!  assimilation_of_fish
!  PCLake_Osis: 
   tDAssFiJv = aDSatFiJv *(self%kDAssFiJv * uFunTmFish * sDFiJv - tDEnvFiJv) 
!-----------------------------------------------------------------------
!  zooplankton predated by fish
!-----------------------------------------------------------------------
!  zooplankton_consumption_of_fish
!  PCLake_osis:sDFiJv,in g/m^2
    tDConsFiJv = tDAssFiJv / self%fDAssFiJv
!  (zooplankton)_P_consumption_by_FiJv
!  PCLake_osis:sPFiJv,in g/m^2 
   tPConsFiJv = rPDZoo * tDConsFiJv
!  (zooplankton)_N_consumption_by_FiJv
!  PCLake_osis:sNFiJv,in g/m^2
   tNConsFiJv = rNDZoo * tDConsFiJv
!-----------------------------------------------------------------------
!  young fish assimilation_P
!-----------------------------------------------------------------------
!  P_assim._efficiency_of_FiJv
   afPAssFiJv = min(1.0_rk,self%cPDFishRef / rPDZoo * self%fDAssFiJv)
!  P_assimilation_of_FiJv
!  PCLake_osis:sPFiJv,in g/m^2
   tPAssFiJv = afPAssFiJv * tPConsFiJv
!-----------------------------------------------------------------------
!  young fish assimilation_N
!-----------------------------------------------------------------------
!  N_assim._efficiency_of_FiJv
   afNAssFiJv = min(1.0_rk,self%cNDFishRef / rNDZoo * self%fDAssFiJv)
!  N_assimilation_of_FiJv
!  PCLake_osis:sNFiJv,in g/m^2
   tNAssFiJv = afNAssFiJv * tNConsFiJv
!-----------------------------------------------------------------------
!  young fish migration
!-----------------------------------------------------------------------
!  migration_flux of young fish, DW
!  PCLake_osis:sDFiJv,in g/m^2
   tDMigrFiJv = self%kMigrFish *(self%cDFiJvIn - sDFiJv)
!  net_migration_flux of young fish,P
!  PCLake_osis:sPFiJv,in g/m^2
   tPMigrFiJv = self%kMigrFish *(self%cPDFishRef * self%cDFiJvIn - sPFiJv)
!  net_migration_flux of young fish,N
!  PCLake_osis:sNFiJv,in g/m^2
   tNMigrFiJv = self%kMigrFish *(self%cNDFishRef * self%cDFiJvIn - sNFiJv)
!-----------------------------------------------------------------------
!  adult fish migration
!-----------------------------------------------------------------------
!  migration_flux of adult fish,DW
!  PCLake_osis:sDFiAd,in g/m^2
   tDMigrFiAd = self%kMigrFish *(self%cDFiAdIn - sDFiAd)  
!  net_migration_flux of adult fish, P
!  PCLake_osis:sPFiAd,in g/m^2
   tPMigrFiAd = self%kMigrFish *(self%cPDFishRef * self%cDFiAdIn - sPFiAd)
!  net_migration_flux of adult fish, N
!  PCLake_osis:sPFiAd,in g/m^2
   tNMigrFiAd = self%kMigrFish *(self%cNDFishRef * self%cDFiAdIn - sNFiAd)
!-----------------------------------------------------------------------
!  fish reproduction
!-----------------------------------------------------------------------
!  Reproduction_flux_DW
   if (Day >= self%cDayReprFish .and. Day < self%cDayReprFish + 1.0_rk) then
!  PCLake_osis:sDFiAd,in g/m^2
   tDReprFish = self%fReprFish * sDFiAd/secs_pr_day
!   wDReprFish = self%fReprFish * sDFiAd/secs_pr_day  ! sDFiAd,in g/m^3
   else
   tDReprFish =0.0_rk
!   wDReprFish = 0.0_rk 
   endif
!  Reproduction_flux_P
   tPReprFish = rPDFiAd * tDReprFish
!    wPReprFish = rPDFiAd * wDReprFish
!  Reproduction_flux_N
   tNReprFish = rNDFiAd * tDReprFish
!   wNReprFish = rNDFiAd * wDReprFish
!-----------------------------------------------------------------------
!  fish aging
!-----------------------------------------------------------------------
!  Ageing_DW
   if (Day >=  364.0_rk .AND. Day <= 365.0_rk) then
!  PCLake_osis:sDFiAd,in g/m^2
   tDAgeFish = self%fAgeFish * sDFiJv/secs_pr_day
!   wDAgeFish = self%fAgeFish * sDFiJv/secs_pr_day
   else
!   wDAgeFish = 0.0_rk
    tDAgeFish = 0.0_rk
   endif
!  Ageing_P
   tPAgeFish = rPDFiJv * tDAgeFish
!   wPAgeFish = rPDFiJv * wDAgeFish
!  Ageing_N
   tNAgeFish = rNDFiJv * tDAgeFish
!   wNAgeFish = rNDFiJv * wDAgeFish
!-----------------------------------------------------------------------
!  young fish respiration and excretion
!-----------------------------------------------------------------------
!  respiration_of_fish_DW
!  PCLake_osis:sDFiAd,in g/m^2
   tDRespFiJv = (self%cPDFishRef / rPDFiJv) * self%kDRespFiJv * uFunTmFish * sDFiJv
!  P_excretion_of_FiJv
!  PCLake_osis:sPFiAd,in g/m^2
   tPExcrFiJv = (rPDFiJv / self%cPDFishRef) * self%kDRespFiJv * uFunTmFish * sPFiJv
!  N_excretion_of_FiJv
!  PCLake_osis:sNFiAd,in g/m^2
   tNExcrFiJv = (rNDFiJv / self%cNDFishRef) * self%kDRespFiJv * uFunTmFish * sNFiJv
!-----------------------------------------------------------------------
!  adult fish respiration and excretion
!-----------------------------------------------------------------------
!  respiration_of_fish
!  PCLake_osis:sDFiAd,in g/m^2
   tDRespFiAd = (self%cPDFishRef / rPDFiAd) * self%kDRespFiAd * uFunTmFish * sDFiAd
!  P_excretion_of_FiAd
!  PCLake_osis:sPFiAd,in g/m^2
   tPExcrFiAd = (rPDFiAd / self%cPDFishRef) * self%kDRespFiAd * uFunTmFish * sPFiAd
!  N_excretion_of_FiAd
!  PCLake_osis:sNFiAd,in g/m^2
   tNExcrFiAd = (rNDFiAd / self%cNDFishRef) * self%kDRespFiAd * uFunTmFish * sNFiAd
!-----------------------------------------------------------------------
!  young fish mortality
!-----------------------------------------------------------------------
!  fish_mortality_incl._environmental_correction
!  PCLake_osis:sDFiAd,in g/m^2
   tDMortFiJv = self%kMortFiJv * sDFiJv +(1.0_rk - aDSatFiJv) * tDEnvFiJv
!  mortality_of_FiJv_P
!  PCLake_osis:sPFiAd,in g/m^2
   tPMortFiJv = rPDFiJv * tDMortFiJv
!  mortality_of_FiJv_N
!  PCLake_osis:sNFiAd,in g/m^2
   tNMortFiJv = rNDFiJv * tDMortFiJv
!-----------------------------------------------------------------------
!  adult fish mortality
!-----------------------------------------------------------------------
!  fish_mortality_incl._environmental_correction
!  PCLake_osis:sDFiAd,in g/m^2
   tDMortFiAd = self%kMortFiAd * sDFiAd +(1.0_rk - aDSatFiAd) * tDEnvFiAd
!  mortality_of_FiAd
!  PCLake_osis:sPFiAd,in g/m^2
   tPMortFiAd = rPDFiAd * tDMortFiAd
!   wPMortFiAd = rPDFiAd * wDMortFiAd
!  mortality_of_FiAd
!  PCLake_osis:sNFiAd,in g/m^2
   tNMortFiAd = rNDFiAd * tDMortFiAd
!-----------------------------------------------------------------------
!  fish egestion
!-----------------------------------------------------------------------
!  egestion_of_fish,young fish
!  PCLake_osis:sDFiAd,in g/m^2
   tDEgesFiJv = tDConsFiJv - tDAssFiJv
!  egestion_of_FiJv
!  PCLake_osis:sPFiAd,in g/m^2
   tNEgesFiJv = tNConsFiJv - tNAssFiJv
!  egestion_of_FiJv
!  PCLake_osis:sNFiAd,in g/m^2
   tPEgesFiJv = tPConsFiJv - tPAssFiJv
!---------------------------------------------------------------------------
!  Piscivorous fish assimilation( this whole area is calibrated in /m^2)
!---------------------------------------------------------------------------
!  vegetation_dependence_of_Pisc_growth_rate
   aFunVegPisc = aDSubVeg /(self%hDVegPisc + aDSubVeg + NearZero)
!  food_limitation_function_of_Pisc
   aDSatPisc = aDFish*aDFish /(self%hDFishPisc*self%hDFishPisc + aDFish*aDFish)
!  intrinsic_net_increase_rate_of_Pisc
   akDIncrPisc = (self%kDAssPisc * aFunVegPisc - self%kDRespPisc) * uFunTmPisc - self%kMortPisc
!  Carrying_capacity_of_Pisc_for_lake_without_marsh_zone
   aDCarrPisc = max(self%cDCarrPiscMin,min(self%cDCarrPiscMax,self%cDCarrPiscBare))
!  environmental_correction_of_Pisc
!  PCLake_osis:sDPisc,in g/m^2
   tDEnvPisc = max(0.0_rk,akDIncrPisc / aDCarrPisc * sDPisc*sDPisc)
!  assimilation_of_Pisc
!  PCLake_osis:sDPisc,in g/m^2
   tDAssPisc = aDSatPisc *(self%kDAssPisc * aFunVegPisc * uFunTmPisc * sDPisc - tDEnvPisc)
!-----------------------------------------------------------------------
!  Piscivorous fish consumption
!-----------------------------------------------------------------------
!  consumption_of_Pisc
!  PCLake_osis:sDPisc,in g/m^2
   tDConsPisc = tDAssPisc / self%fDAssPisc
!   wDConsPisc = wDAssPisc / self%fDAssPisc
!-----------------------------------------------------------------------
!  young fish predated by piscivirious fish
!-----------------------------------------------------------------------
!  young_fish_consumption_by_Pisc_DW
!  PCLake_osis:sDPisc,aDFish,in g/m^2
   tDConsFiJvPisc = sDFiJv / aDFish * tDConsPisc
!  young_fish_consumption_by_Pisc
   tPConsFiJvPisc = rPDFiJv * tDConsFiJvPisc
!  young_fish_consumption_by_Pisc
   tNConsFiJvPisc = rNDFiJv * tDConsFiJvPisc
!-----------------------------------------------------------------------
!  adult fish predated by piscivirious fish
!-----------------------------------------------------------------------
!  adult_fish_consumption_by_Pisc
!  PCLake_osis:sDPisc,aDFish,in g/m^2
   tDConsFiAdPisc = tDConsPisc - tDConsFiJvPisc
!  adult_fish_consumption_by_Pisc
   tPConsFiAdPisc = rPDFiAd * tDConsFiAdPisc
!  adult_fish_consumption_by_Pisc
   tNConsFiAdPisc = rNDFiAd * tDConsFiAdPisc
!-----------------------------------------------------------------------
!  piscivirious fish migration
!-----------------------------------------------------------------------
!  migration_flux
!  PCLake_osis:sDPisc,in g/m^2
   tDMigrPisc = self%kMigrPisc *(self%cDPiscIn - sDPisc)
!-----------------------------------------------------------------------
!  piscivirious fish respiration
!-----------------------------------------------------------------------
!  respiration_of_Pisc
!  PCLake_osis:sDPisc,in g/m^2
   tDRespPisc = self%kDRespPisc * uFunTmPisc * sDPisc
!-----------------------------------------------------------------------
!  piscivirious fish mortality
!-----------------------------------------------------------------------
!  mortality_of_Pisc(incl._environmental_correction)
!  PCLake_osis:sDPisc,in g/m^2
   tDMortPisc = self%kMortPisc * sDPisc +(1.0_rk - aDSatPisc) * tDEnvPisc
!---------------------------------------------------------------------------
!  piscivirious fish N process
!---------------------------------------------------------------------------
!  _Piscivorous_fish
    aPPisc = self%cPDPisc * sDPisc
!  total_P_consumption_by_Pisc
   tPConsPisc = tPConsFiJvPisc + tPConsFiAdPisc
!  average_P/D_ratio_of_Pisc_food
   rPDFoodPisc = tPConsPisc / tDConsPisc
!  P_assim._efficiency_of_Pisc
    afPAssPisc = min(1.0_rk,self%cPDPisc / rPDFoodPisc * self%fDAssPisc)
!  P_assimilation_of_Pisc
   tPAssPisc = afPAssPisc * tPConsPisc
!  respiration_of_Pisc
   tPExcrPisc = self%cPDPisc * tDRespPisc
!  mortality_of_Pisc
   tPMortPisc = self%cPDPisc * tDMortPisc
!  net_migration_flux
   tPMigrPisc = self%kMigrPisc *(self%cPDPisc * self%cDPiscIn - aPPisc)
!-----------------------------------------------------------------------
!  piscivirious fish N process
!-----------------------------------------------------------------------
!  Piscivorous_fish
    aNPisc = self%cNDPisc * sDPisc
!  total_N_consumption_by_Pisc
   tNConsPisc = tNConsFiJvPisc + tNConsFiAdPisc
!  average_N/D_ratio_of_Pisc_food
   rNDFoodPisc = tNConsPisc / tDConsPisc
!  N_assim._efficiency_of_Pisc
    afNAssPisc = min(1.0_rk,self%cNDPisc / rNDFoodPisc * self%fDAssPisc)
!  N_assimilation_of_Pisc
   tNAssPisc = afNAssPisc * tNConsPisc
!  egestion_of_Pisc
   tNEgesPisc = tNConsPisc - tNAssPisc
!  respiration_of_Pisc
   tNExcrPisc = self%cNDPisc * tDRespPisc
!  mortality_of_Pisc
   tNMortPisc = self%cNDPisc * tDMortPisc
!  net_migration_flux
   tNMigrPisc = self%kMigrPisc *(self%cNDPisc * self%cDPiscIn - aNPisc)
!-----------------------------------------------------------------------
!  piscivirious fish egestion
!-----------------------------------------------------------------------
!  egestion_of_Pisc
   tDEgesPisc = tDConsPisc - tDAssPisc
!  egestion_of_Pisc
   tNEgesPisc = tNConsPisc - tNAssPisc
!  egestion_of_Pisc
   tPEgesPisc = tPConsPisc - tPAssPisc
!-----------------------------------------------------------------------
!  total flux of web change to state variables
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_Herbivorous_zooplankton
   wDWebZoo = wDAssZoo - wDRespZoo - wDMortZoo - tDConsFiJv/sDepthW
!  total_foodweb_flux_of_N_in_Herbivorous_zooplankton
   wNWebZoo = wNAssZoo - wNExcrZoo - wNMortZoo  - tNConsFiJv / sDepthW
!  total_foodweb_flux_of_P_in_Herbivorous_zooplankton
   wPWebZoo = wPAssZoo - wPExcrZoo - wPMortZoo - tPConsFiJv / sDepthW
!  total_foodweb_flux_of_DW_in_Young_fish
   tDWebFiJv = tDMigrFiJv + tDReprFish + tDAssFiJv - tDRespFiJv - tDMortFiJv - tDConsFiJvPisc - tDAgeFish
!  temperal solution, vertial averaged
   wDWebFiJv=tDWebFiJv/sDepthW
!  total_foodweb_flux_of_P_in_Young_fish
   tPWebFiJv = tPMigrFiJv + tPReprFish  + tPAssFiJv - tPExcrFiJv - tPMortFiJv - tPConsFiJvPisc - tPAgeFish
!  temperal solution, vertial averaged
   wPWebFiJv = tPWebFiJv/sDepthW
!  total_foodweb_flux_of_N_in_Young_fish
   tNWebFiJv = tNMigrFiJv + tNReprFish + tNAssFiJv - tNExcrFiJv - tNMortFiJv - tNConsFiJvPisc - tNAgeFish
!  temperal solution, vertial averaged
   wNWebFiJv= tNWebFiJv/ sDepthW
!  total_foodweb_flux_of_DW_in_Adult_fish
   tDWebFiAd = tDMigrFiAd - tDRespFiAd - tDMortFiAd - tDReprFish - tDConsFiAdPisc + tDAgeFish
!  temperal solution, vertial averaged
   wDWebFiAd= tDWebFiAd/ sDepthW
!  total_foodweb_flux_of_P_in_Adult_fish
   tPWebFiAd = tPMigrFiAd  - tPExcrFiAd - tPMortFiAd - tPReprFish - tPConsFiAdPisc + tPAgeFish 
!  temperal solution, vertial averaged
   wPWebFiAd= tPWebFiAd/ sDepthW
!  total_foodweb_flux_of_N_in_Adult_fish
   tNWebFiAd = tNMigrFiAd - tNExcrFiAd - tNMortFiAd - tNReprFish - tNConsFiAdPisc + tNAgeFish
!  temperal solution, vertial averaged
   wNWebFiAd= tNWebFiAd/ sDepthW
!  total_foodweb_flux_of_DW_in_predatory_fish
   tDWebPisc = tDMigrPisc + tDAssPisc - tDRespPisc - tDMortPisc
   wDWebPisc=tDWebPisc/sDepthW
!=======================================================================
!  foodweb part relating to other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in water
!-----------------------------------------------------------------------
!  soluble_N_egestion
   wNEgesZooNH4 = self%fDissEgesZoo*wNEgesZoo
!  soluble_N_mortality
   wNMortZooNH4 = self%fDissMortZoo*wNMortZoo
!-----------------------------------------------------------------------
!  for fish it has t-, unit in /m^2
!-----------------------------------------------------------------------
!  NH4_egestion_of_young_fish
   tNEgesFiJvNH4 = self%fDissEgesFish * tNEgesFiJv
!  total fish mortality, N
   tNMortFish = tNMortFiJv + tNMortFiAd
!  part_of_died_fish_N_fixed_in_bones_AND_scales
   tNMortFishBot = self%fDBone * tNMortFish
!  part_of_died_fish_N_becoming_dissolved_N
   tNMortFishNH4 = self%fDissMortFish *(tNMortFish - tNMortFishBot)
!  SRN_egestion_of_Pisc
   tNEgesPiscNH4 = self%fDissEgesPisc * tNEgesPisc
!  part_of_died_Pisc_N_fixed_in_bones_AND_scales
   tNMortPiscBot = self%fDBone * tNMortPisc
!  part_of_died_fish_N_becoming_dissolved_N
   tNMortPiscNH4 = self%fDissMortPisc *(tNMortPisc - tNMortPiscBot)
!  total_foodweb_flux_of_N_in_ammonium_in_water_in_lake_water
   wNWebNH4W = wNExcrZoo + wNEgesZooNH4 + wNMortZooNH4 +(tNExcrFiJv + tNExcrFiAd + &
   &tNEgesFiJvNH4 + tNMortFishNH4 + tNExcrPisc + tNEgesPiscNH4 + tNMortPiscNH4)/sDepthW 
!-----------------------------------------------------------------------
!  Update NO3 in water   (no NO3????)
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_N_in_nitrate_in_water_in_lake_water
   wNWebNO3W = 0.0_rk
!-----------------------------------------------------------------------
!  Update PO4 in water
!-----------------------------------------------------------------------
!  soluble_P_egestion
   wPEgesZooPO4 = self%fDissEgesZoo*wPEgesZoo
!  soluble_P_mortality
   wPMortZooPO4 = self%fDissMortZoo * wPMortZoo
!-----------------------------------------------------------------------
!  for fish it has t-, unit in /m^2
!-----------------------------------------------------------------------
!  SRP_egestion_of_young_fish
   tPEgesFiJvPO4 = self%fDissEgesFish * tPEgesFiJv
!  total fish mortality
   tPMortFish = tPMortFiJv + tPMortFiAd
!  part_of_died_fish_P_fixed_in_bones_AND_scales
   tPMortFishBot = self%fPBone * tPMortFish
!  part_of_died_fish_P_becoming_dissolved_P
   tPMortFishPO4 = self%fDissMortFish *(tPMortFish - tPMortFishBot)
!  SRP_egestion_of_Pisc
   tPEgesPiscPO4 = self%fDissEgesPisc * tPEgesPisc
!  part_of_died_Pisc_P_fixed_in_bones_AND_scales
   tPMortPiscBot = self%fPBone * tPMortPisc
!  part_of_died_fish_P_becoming_dissolved_P
   tPMortPiscPO4 = self%fDissMortPisc *(tPMortPisc - tPMortPiscBot)
!  total_foodweb_flux_of_P_in_SRP_in_water_in_lake_water
   wPWebPO4W = wPExcrZoo + wPEgesZooPO4 + wPMortZooPO4 +(tPExcrFiJv + tPExcrFiAd + &
   &tPEgesFiJvPO4 + tPMortFishPO4 + tPExcrPisc + tPEgesPiscPO4 + tPMortPiscPO4)/sDepthW
!-----------------------------------------------------------------------
!  Update detrital DW in water
!  for fish it has t-, unit in /m^2
!-----------------------------------------------------------------------
!  bent._fish_mortality
   tDMortFish = tDMortFiJv + tDMortFiAd
!  part_of_died_fish_DW_fixed_in_bones_and_scales
   tDMortFishBot = self%fDBone * tDMortFish
!  part_of_died_fish_DW_becoming_detritus
   tDMortFishDet = tDMortFish - tDMortFishBot
!  part_of_died_fish_DW_fixed_in_bones_AND_scales
   tDMortPiscBot = self%fDBone * tDMortPisc
!  part_of_died_Pisc_DW_becoming_detritus
   tDMortPiscDet = tDMortPisc - tDMortPiscBot
!  total_foodweb_flux_of_DW_in_Detritus_in_lake_water
   wDWebDetW = - wDConsDetZoo + wDEgesZoo + wDMortZoo + (tDEgesFiJv + &
   &tDMortFishDet + tDEgesPisc + tDMortPiscDet)/sDepthW
!-----------------------------------------------------------------------
!  Update detrital N in water
!-----------------------------------------------------------------------
!  part_of_died_Pisc_N_becoming_detrital_N
   tNMortPiscDet = tNMortPisc - tNMortPiscBot - tNMortPiscNH4
!  detrital_N_egestion_of_Pisc
   tNEgesPiscDet = tNEgesPisc - tNEgesPiscNH4
!  part_of_died_fish_NW_becoming_detritus
   tNMortFishDet = tNMortFish - tNMortFishBot - tNMortFishNH4
!  detrital_N_egestion_of_young_fish
   tNEgesFiJvDet = tNEgesFiJv - tNEgesFiJvNH4
!  detrital_N_mortality
   wNMortZooDet = wNMortZoo - wNMortZooNH4
!  detrital_N_egestion
   wNEgesZooDet = wNEgesZoo - wNEgesZooNH4
!  total_foodweb_flux_of_N_in_Detritus_in_lake_water
   wNWebDetW = - wNConsDetZoo + wNEgesZooDet + wNMortZooDet +(tNEgesFiJvDet + &
   &tNMortFishDet + tNEgesPiscDet + tNMortPiscDet)/sDepthW
!-----------------------------------------------------------------------
!  Update detrital P in water
!-----------------------------------------------------------------------
!  part_of_died_Pisc_P_becoming_detrital_P
   tPMortPiscDet = tPMortPisc - tPMortPiscBot - tPMortPiscPO4
!  detrital_P_egestion_of_Pisc
   tPEgesPiscDet = tPEgesPisc - tPEgesPiscPO4
!  part_of_died_fish_PW_becoming_detritus
   tPMortFishDet = tPMortFish - tPMortFishBot - tPMortFishPO4
!  detrital_P_egestion_of_young_fish
   tPEgesFiJvDet = tPEgesFiJv - tPEgesFiJvPO4
!  detrital_P_mortality
   wPMortZooDet = wPMortZoo - wPMortZooPO4
!  detrital_P_egestion
   wPEgesZooDet = wPEgesZoo - wPEgesZooPO4
!  total_foodweb_flux_of_P_in_Detritus_in_lake_water
   wPWebDetW = - wPConsDetZoo + wPEgesZooDet + wPMortZooDet + (tPEgesFiJvDet + &
   &tPMortFishDet + tPEgesPiscDet + tPMortPiscDet)/sDepthW
!-----------------------------------------------------------------------
!  Update detrital Si in water
!-----------------------------------------------------------------------
!  consumption_of_diatoms
   wSiConsDiatZoo = self%cSiDDiat * wDConsDiatZoo
!  total_foodweb_flux_of_silica_in_lake_water_detritus
   wSiWebDetW = wSiConsDiatZoo
!-----------------------------------------------------------------------
!  Update diatom state variables
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_Diatoms_in_lake_water
   wDWebDiatW = - wDConsDiatZoo
!  total_foodweb_flux_of_N_in_Diatoms_in_lake_water
   wNWebDiatW = - wNConsDiatZoo
!  total_foodweb_flux_of_P_in_Diatoms_in_lake_water
   wPWebDiatW = - wPConsDiatZoo 
!-----------------------------------------------------------------------
!  Update green algae state variables
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_Greens_in_lake_water
   wDWebGrenW = - wDConsGrenZoo
!  total_foodweb_flux_of_N_in_Greens_in_lake_water
   wNWebGrenW = - wNConsGrenZoo
!  total_foodweb_flux_of_P_in_Greens_in_lake_water
   wPWebGrenW = - wPConsGrenZoo
!-----------------------------------------------------------------------
!  Update blue algae state variables
!-----------------------------------------------------------------------
!  total_foodweb_flux_of_DW_in_Blue-greens_in_lake_water
   wDWebBlueW = - wDConsBlueZoo
!  total_foodweb_flux_of_N_in_Blue-greens_in_lake_water
   wNWebBlueW = - wNConsBlueZoo
!  total_foodweb_flux_of_P_in_Blue-greens_in_lake_water
   wPWebBlueW = - wPConsBlueZoo
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
   _SET_ODE_(self%id_sDZoo,wDWebZoo)
   _SET_ODE_(self%id_sNZoo,wNWebZoo)
   _SET_ODE_(self%id_sPZoo,wPWebZoo)
   _SET_ODE_(self%id_sDFiJv,wDWebFiJv)
   _SET_ODE_(self%id_sPFiJv,wPWebFiJv)
   _SET_ODE_(self%id_sNFiJv,wNWebFiJv)
   _SET_ODE_(self%id_sDFiAd,wDWebFiAd)
   _SET_ODE_(self%id_sPFiAd,wPWebFiAd)
   _SET_ODE_(self%id_sNFiAd,wNWebFiAd)
   _SET_ODE_(self%id_sDPisc,wDWebPisc)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
!  update abiotic variables in water
   _SET_ODE_(self%id_NH4poolW,wNWebNH4W)
   _SET_ODE_(self%id_NO3poolW,wNWebNO3W)
   _SET_ODE_(self%id_PO4poolW,wPWebPO4W)
   _SET_ODE_(self%id_DDetpoolW,wDWebDetW)
   _SET_ODE_(self%id_NDetpoolW,wNWebDetW)
   _SET_ODE_(self%id_PDetpoolW,wPWebDetW)
   _SET_ODE_(self%id_SiDetpoolW,wSiWebDetW)
!  update phytoplankton in water
   _SET_ODE_(self%id_DfoodDiat,wDWebDiatW)
   _SET_ODE_(self%id_NfoodDiat,wNWebDiatW)
   _SET_ODE_(self%id_PfoodDiat,wPWebDiatW)
   _SET_ODE_(self%id_DfoodGren,wDWebGrenW)
   _SET_ODE_(self%id_NfoodGren,wNWebGrenW)
   _SET_ODE_(self%id_PfoodGren,wPWebGrenW)
   _SET_ODE_(self%id_DfoodBlue,wDWebBlueW)
   _SET_ODE_(self%id_NfoodBlue,wNWebBlueW)
   _SET_ODE_(self%id_PfoodBlue,wPWebBlueW)
!-----------------------------------------------------------------------
!  output diagnostic variables for external links
!-----------------------------------------------------------------------
!  Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_aNPisc,aNPisc)
   _SET_DIAGNOSTIC_(self%id_aPPisc,aPPisc)

   
   _LOOP_END_
! Spatial loop end
!
!EOP
!-----------------------------------------------------------------------

   end subroutine do
!BOP
!IROUTINE: Get the light extinction coefficient due to biogeochemical variables
! !DESCRIPTION:
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_pclake_foodweb_water), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:


!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   !! Retrieve current (local) state variable values.

   !
   !! Self-shading with explicit contribution from background
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------


   end module pclake_foodweb_water

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

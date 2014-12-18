#include "fabm_driver.h"
!
!BOP
!
! !INTERFACE:
   module au_pclake_auxilary
!
! !DESCRIPTION:
!  au_pclake_auxilary is created for the purpose of computing resuspension and sedimentation,
!  sediment burial processes. 
!  No local state variable is registed here.
!  resuspension and sedimentation involve: sDIMW<==>sDIMS,sD/N/PDetW<==>sD/N/PDetS,sPAIMW<==>sPAIMS
!  aD/N/PPhytW<==>aD/N/PPhytS,sDDiatW<==>sDDiatS,sNH4S==>sNH4W,sNO3S==>sNO3W,sPO4S==>sPO4W
!  Burial process involve: sDIMS==>,sD/N/PDetS==>,sPAIMS==>
!  feh: Sep.8
!  Diatom Si sedimentation and resuspension can't ben handled here, since SiDiat is not state variable both 
!  in water column and sediment. Something could be further considered.
! !USES:
   use fabm_types
   use fabm_expressions
   use au_pclake_utility, ONLY: uFunTmAbio,uFunTmBio

   implicit none

!  default: all is private.
      private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_au_pclake_auxilary
!     local state variable identifers

!     diagnostic variables for local output(should be discussed and added later???)
!     state dependencies identifers
!     SW:Sediment to Water
!     dependencies to pclake_abiotic_water state variables
      type (type_state_variable_id)            :: id_SWNH4,id_SWNO3,id_SWPO4,id_SWPAIM,id_SWO2,id_SWSiO2
      type (type_state_variable_id)            :: id_SWDIM,id_SWDDet,id_SWNDet,id_SWPDet,id_SWSiDet
!     dependencies to phy_wat state variables
      type (type_state_variable_id)            :: id_SWDBlue,id_SWNBlue,id_SWPBlue
      type (type_state_variable_id)            :: id_SWDDiat,id_SWNDiat,id_SWPDiat  !,id_SWSiDiat
      type (type_state_variable_id)            :: id_SWDGren,id_SWNGren,id_SWPGren
!     WS:Water to Sediment
!     dependencies to pclake_abiotic_sediment state variables
      type (type_bottom_state_variable_id)     :: id_WSPO4,id_WSPAIM,id_WSNH4,id_WSNO3
      type (type_bottom_state_variable_id)     :: id_WSDIM,id_WSDDet,id_WSNDet,id_WSPDet,id_WSSiDet
      type (type_bottom_state_variable_id)     :: id_WSDHum,id_WSNHum,id_WSPHum
!     dependencies to phy_sed state variables
      type (type_bottom_state_variable_id)     :: id_WSDBlue,id_WSNBlue,id_WSPBlue
      type (type_bottom_state_variable_id)     :: id_WSDDiat,id_WSNDiat,id_WSPDiat  !,id_WSSiDiat
      type (type_bottom_state_variable_id)     :: id_WSDGren,id_WSNGren,id_WSPGren
!     dependencies to foodweb_wat state variables(zooplankton only created for transportation purpose)
      type (type_state_variable_id)            :: id_TurbFish,id_DTranZoo,id_NTranZoo,id_PTranZoo
!     dependencies to vegetation state variables
      type (type_bottom_state_variable_id)     :: id_DragVeg
!     environmental dependencies
      type (type_global_dependency_id)         :: id_Day
      type (type_dependency_id)                :: id_uTm !,id_dz
      type (type_horizontal_dependency_id)     :: id_sDepthW,id_shear
      type (type_horizontal_diagnostic_variable_id)       :: id_tDBurIM,id_shearstress
      type (type_horizontal_diagnostic_variable_id)       :: id_aFunDimSusp,id_aFunTauSet,id_tDResusDead
!     diagnostic dependencies,due to burial process
      type ( type_horizontal_dependency_id)             :: id_tDAbioHumS
      type ( type_horizontal_dependency_id)             :: id_tDAbioDetS,id_tDPrimDetS,id_tDWebDetS,id_tDBedDetS
!!    Model parameters
!     logic variables for whether linking dependencies
!     diagnostic dependencies,due to resuspension
      logical      :: fish_module,vegetation_module
      real(rk)                   :: cDepthS
!     sediment properties parameters
      real(rk)                   :: fLutum,fLutumRef,bPorS
!     resuspension process parameters
      real(rk)                   :: kVegResus,kTurbFish
      real(rk)                   :: cSuspRef,cSuspMin,cSuspMax,cSuspSlope
      real(rk)                   :: hDepthSusp,cFetchRef,cFetch
      real(rk)                   :: cResusPhytExp,kResusPhytMax  !,cSiDDiat
!     sedimentation parameters
      real(rk)                   :: cThetaSet,cVSetIM,cVSetDet
      real(rk)                   :: cVSetDiat,cVSetGren,cVSetBlue
!     Burial process parameters
      real(rk)                   :: cRhoIM,cRhoOM,fDOrgSoil
      real(rk)                   :: cPO4Ground,cNH4Ground,cNO3Ground
!     parameter for fish temperature function
      real(rk)                   :: cSigTmFish,cTmOptFish
!!    variables need to be removed after compiled to 1d to 3d physical drivers
!     parameters for loadings
      real(rk)                   :: cLoadIM, cDLoadDet,cPLoadDet,cNLoadDet, cLoadPO4
      real(rk)                   :: cLoadPAIM,cLoadNH4,cLoadNO3,uQIn
!     bank erosion par
      real(rk)                   :: cDErosTot,fSedErosIM
!     nutrient ratios parameter
      real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
      real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax
!     parameters relating to new resuspension method
      real(rk)   :: crt_shear,ref_shear,alpha,eta,cVSetMain
!     variable for selecting different resuspension methods
      integer    :: resusp_meth
      
   contains
   
!     Model procedures
      procedure :: initialize
      procedure :: do_bottom
      procedure :: do
   end type type_au_pclake_auxilary
!  private data memebers(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk),parameter :: Pi=3.14159265358979_rk
   real(rk), parameter :: mmPerm=1000.0_rk
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
! !   DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_au_pclake_auxilary), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!   LOCAL VARIABLES:
   real(rk)                  :: cDepthS
   real(rk)                  :: cThetaSet
   real(rk)                  :: kVegResus
   real(rk)                  :: kTurbFish
   real(rk)                  :: cSuspRef
   real(rk)                  :: cSuspMin
   real(rk)                  :: cSuspMax
   real(rk)                  :: cSuspSlope
   real(rk)                  :: hDepthSusp
   real(rk)                  :: cFetchRef
   real(rk)                  :: cFetch
   real(rk)                  :: fLutum
   real(rk)                  :: fLutumRef
   real(rk)                  :: bPorS
   real(rk)                  :: cResusPhytExp
   real(rk)                  :: kResusPhytMax
   real(rk)                  :: cVSetIM
   real(rk)                  :: cVSetDet
   real(rk)                  :: cVSetDiat
   real(rk)                  :: cVSetGren
   real(rk)                  :: cVSetBlue
   real(rk)                  :: cRhoIM
   real(rk)                  :: cRhoOM
   real(rk)                  :: fDOrgSoil
   real(rk)                  :: cPO4Ground
   real(rk)                  :: cNH4Ground
   real(rk)                  :: cNO3Ground
   real(rk)                  :: cSigTmFish
   real(rk)                  :: cTmOptFish
   real(rk)                  ::cNDBlueMax
   real(rk)                  ::cNDBlueMin
   real(rk)                  ::cNDDiatMax
   real(rk)                  ::cNDDiatMin
   real(rk)                  ::cNDGrenMax
   real(rk)                  ::cNDGrenMin
   real(rk)                  ::cPDBlueMax
   real(rk)                  ::cPDBlueMin
   real(rk)                  ::cPDDiatMax
   real(rk)                  ::cPDDiatMin
   real(rk)                  ::cPDGrenMax
   real(rk)                  ::cPDGrenMin
   character(len=64)         :: Amonia_pool_in_water
   character(len=64)         :: Nitrates_pool_in_water
   character(len=64)         :: Phosphate_pool_in_water
   character(len=64)         :: Absorbed_phosphrus_in_water
   character(len=64)         :: Oxygen_pool_in_water
   character(len=64)         :: Inorg_pool_in_water
   character(len=64)         :: Detritus_DW_in_water
   character(len=64)         :: Detritus_N_in_water
   character(len=64)         :: Detritus_P_in_water
   character(len=64)         :: Detritus_Si_in_water
   character(len=64)         :: Amonia_pool_in_sediment
   character(len=64)         :: Nitrates_pool_in_sediment
   character(len=64)         :: Phosphate_pool_in_sediment
   character(len=64)         :: Absorbed_phosphrus_in_sediment
   character(len=64)         :: Inorg_pool_in_sediment
   character(len=64)         :: Detritus_DW_in_sediment
   character(len=64)         :: Detritus_N_in_sediment
   character(len=64)         :: Detritus_P_in_sediment
   character(len=64)         :: Detritus_Si_in_sediment
   character(len=64)         :: Diatom_DW_in_water
   character(len=64)         :: Green_DW_in_water
   character(len=64)         :: Blue_DW_in_water
   character(len=64)         :: Diatom_N_in_water
   character(len=64)         :: Green_N_in_water
   character(len=64)         :: Blue_N_in_water
   character(len=64)         :: Diatom_P_in_water
   character(len=64)         :: Green_P_in_water
   character(len=64)         :: Blue_P_in_water
   character(len=64)         :: Diatom_DW_in_sediment
   character(len=64)         :: Green_DW_in_sediment
   character(len=64)         :: Blue_DW_in_sediment
   character(len=64)         :: Diatom_N_in_sediment
   character(len=64)         :: Green_N_in_sediment
   character(len=64)         :: Blue_N_in_sediment
   character(len=64)         :: Diatom_P_in_sediment
   character(len=64)         :: Green_P_in_sediment
   character(len=64)         :: Blue_P_in_sediment
   character(len=64)         :: vegetation_DW
   character(len=64)         :: adult_fish_DW
   character(len=64)         :: Detritus_abiotic_update
   character(len=64)         :: Detritus_from_algae
   character(len=64)         :: Detritus_from_vegetation
   character(len=64)         :: Detritus_from_foodweb
!  Humus variables
   character(len=64)         :: Humus_DW_in_sediment
   character(len=64)         :: Humus_N_in_sediment
   character(len=64)         :: Humus_P_in_sediment
   character(len=64)         :: Humus_abiotic_update
!  zooplanktons, only for transportation
   character(len=64)         :: Zooplankton_DW
   character(len=64)         :: Zooplankton_N
   character(len=64)         :: Zooplankton_P
   character(len=64)         :: SiO2_pool_water
   
!  loading parameters
   real(rk)                  :: cLoadPO4
   real(rk)                  :: cLoadNO3
   real(rk)                  :: uQIn
   real(rk)                  :: cDErosTot
   real(rk)                  :: fSedErosIM
!  parameters for shear stress related resuspension and sedimentation
   real(rk)                  :: crt_shear
   real(rk)                  :: ref_shear
   real(rk)                  :: alpha
   real(rk)                  :: eta
   real(rk)                  :: cVSetMain
   integer                  :: resusp_meth

!  create namelist
   namelist /au_pclake_auxilary/ cDepthS,cThetaSet,kVegResus,kTurbFish,cSuspRef,cSuspMin, &
                           & cSuspMax,cSuspSlope,hDepthSusp,cFetchRef,cFetch,fLutum,fLutumRef, &
                           & bPorS,cResusPhytExp,kResusPhytMax,cVSetIM,cVSetDet,cVSetDiat,  &
                           & cVSetGren,cVSetBlue,cRhoIM,cRhoOM,fDOrgSoil,&
                           & cPO4Ground,cNH4Ground,cNO3Ground,cSigTmFish,cTmOptFish,&  !cSiDDiat,
                           & cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin, &
                           & cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax, &
                           & Amonia_pool_in_water,Nitrates_pool_in_water,Phosphate_pool_in_water, &
                           & Absorbed_phosphrus_in_water,Oxygen_pool_in_water,Inorg_pool_in_water,Detritus_DW_in_water,&
                           & Detritus_N_in_water,Detritus_P_in_water,Detritus_Si_in_water,Amonia_pool_in_sediment, &
                           & Nitrates_pool_in_sediment,Phosphate_pool_in_sediment,Absorbed_phosphrus_in_sediment,Inorg_pool_in_sediment,&
                           & Detritus_DW_in_sediment,Detritus_N_in_sediment,Detritus_P_in_sediment,Detritus_Si_in_sediment,&
                           & Diatom_DW_in_water,Green_DW_in_water,Blue_DW_in_water,Diatom_N_in_water,Green_N_in_water,Blue_N_in_water, &
                           & Diatom_P_in_water,Green_P_in_water,Blue_P_in_water,Diatom_DW_in_sediment,Green_DW_in_sediment,&
                           & Blue_DW_in_sediment,Diatom_N_in_sediment,Green_N_in_sediment,Blue_N_in_sediment, &  !Diatom_Si_in_sediment,Diatom_Si_in_water,
                           & Diatom_P_in_sediment,Green_P_in_sediment,Blue_P_in_sediment,vegetation_DW,adult_fish_DW, &
                           & Detritus_abiotic_update,Detritus_from_algae,Detritus_from_vegetation,Detritus_from_foodweb, &
!  humus linking strings
                           & Humus_DW_in_sediment,Humus_N_in_sediment,Humus_P_in_sediment, Humus_abiotic_update,&
!  zooplankton linking strings
                           & Zooplankton_DW,Zooplankton_N,Zooplankton_P,SiO2_pool_water, &
!  paramters for loadings
                           & cLoadPO4,cLoadNO3,uQIn,cDErosTot,fSedErosIM, &
                          & crt_shear,ref_shear,alpha,eta,cVSetMain,resusp_meth

!EOP
!-----------------------------------------------------------------------
!BOC
!  initialize the parameters
   cDepthS = 0.1_rk
   cThetaSet=1.01_rk
   kVegResus = 0.01_rk
   kTurbFish = 1.0_rk
   cSuspRef = 0.5_rk
   cSuspMin = 6.1_rk
   cSuspMax = 25.2_rk
   cSuspSlope = 2.1_rk
   hDepthSusp = 2.0_rk
   cFetchRef = 1000.0_rk
   cFetch = 1000.0_rk
   fLutum = 0.1_rk
   fLutumRef = 0.2_rk
   bPorS = 0.85_rk
   kResusPhytMax=0.25_rk
   cResusPhytExp=-0.379
   cVSetIM=1.0_rk
   cVSetDet=0.25_rk
   cVSetDiat= 0.5_rk
   cVSetGren= 0.2_rk
   cVSetBlue= 0.06_rk
   cRhoIM=2500000.0_rk
   cRhoOM=1400000.0_rk
   fDOrgSoil=0.1_rk
   cPO4Ground=0.1_rk
   cNH4Ground=1.0_rk
   cNO3Ground=0.1_rk
   cTmOptFish=25.0_rk
   cSigTmFish=10.0_rk
!   cSiDDiat=0.15_rk
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
   Amonia_pool_in_water= ''
   Nitrates_pool_in_water=''
   Phosphate_pool_in_water= ''
   Absorbed_phosphrus_in_water=''
   Oxygen_pool_in_water=''
   Inorg_pool_in_water=''
   Detritus_DW_in_water=''
   Detritus_N_in_water=''
   Detritus_P_in_water=''
   Detritus_Si_in_water=''
   Amonia_pool_in_sediment= ''
   Nitrates_pool_in_sediment=''
   Phosphate_pool_in_sediment= ''
   Absorbed_phosphrus_in_sediment=''
   Inorg_pool_in_sediment=''
   Detritus_DW_in_sediment=''
   Detritus_N_in_sediment=''
   Detritus_P_in_sediment=''
   Detritus_Si_in_sediment=''
   Diatom_DW_in_water=''
   Green_DW_in_water=''
   Blue_DW_in_water=''
   Diatom_N_in_water=''
   Green_N_in_water=''
   Blue_N_in_water=''
   Diatom_P_in_water=''
   Green_P_in_water=''
   Blue_P_in_water=''
!   Diatom_Si_in_water=''
   Diatom_DW_in_sediment=''
   Green_DW_in_sediment=''
   Blue_DW_in_sediment=''
   Diatom_N_in_sediment=''
   Green_N_in_sediment=''
   Blue_N_in_sediment=''
   Diatom_P_in_sediment=''
   Green_P_in_sediment=''
   Blue_P_in_sediment=''
!   Diatom_Si_in_sediment=''
   vegetation_DW=''
   adult_fish_DW=''
   Detritus_abiotic_update=''
   Detritus_from_algae=''
   Detritus_from_vegetation=''
   Detritus_from_foodweb=''
!  for humus
   Humus_DW_in_sediment=''
   Humus_N_in_sediment=''
   Humus_P_in_sediment=''
   Humus_abiotic_update=''
!  For zooplankton transport
   Zooplankton_DW=''
   Zooplankton_N=''
   Zooplankton_P=''
   SiO2_pool_water=''
!  for loadings
   cLoadPO4=0.0_rk
   cLoadNO3=0.0_rk
   uQIn=0.0_rk
!  bank erosion 
   cDErosTot=0.005_rk
   fSedErosIM=0.95_rk
!  parameters for shear stress related resuspension
   crt_shear=0.1_rk
   ref_shear=1.0_rk
   alpha=1.0_rk
   eta=1.0_rk
   cVSetMain= 0.5_rk
   resusp_meth=2
!  Read parameters namelist
   if (configunit>0) read(configunit,nml=au_pclake_auxilary,err=99,end=100)
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cDepthS,'cDepthS',default=cDepthS)
   call self%get_parameter(self%cThetaSet,'cThetaSet',default=cThetaSet)
   call self%get_parameter(self%kVegResus,'kVegResus',default=kVegResus)
   call self%get_parameter(self%kTurbFish,'kTurbFish',default=kTurbFish)  ! don't convert here, due to PCLake imperical function.
   call self%get_parameter(self%cSuspRef,'cSuspRef',default=cSuspRef)
   call self%get_parameter(self%cSuspMin,'cSuspMin',default=cSuspMin)
   call self%get_parameter(self%cSuspMax,'cSuspMax',default=cSuspMax)
   call self%get_parameter(self%cSuspSlope,'cSuspSlope',default=cSuspSlope)
   call self%get_parameter(self%hDepthSusp,'hDepthSusp',default=hDepthSusp)
   call self%get_parameter(self%cFetchRef,'cFetchRef',default=cFetchRef)
   call self%get_parameter(self%cFetch,'cFetch',default=cFetch)
   call self%get_parameter(self%fLutum ,'fLutum ',default=fLutum)
   call self%get_parameter(self%fLutumRef,'fLutumRef',default=fLutumRef)
   call self%get_parameter(self%bPorS,'bPorS',default=bPorS)
!  this parameter shouldn't convert to secs, due to the equation 
   call self%get_parameter(self%cResusPhytExp,'cResusPhytExp',default=cResusPhytExp) !,scale_factor=1.0_rk/secs_pr_day) 
   call self%get_parameter(self%kResusPhytMax,'kResusPhytMax',default=kResusPhytMax) !,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVSetIM,'cVSetIM',default=cVSetIM) !,scale_factor=1.0_rk/secs_pr_day scale it after equation
   call self%get_parameter(self%cVSetDet,'cVSetDet',default=cVSetDet)! ,scale_factor=1.0_rk/secs_pr_day)  scale it after equation
   call self%get_parameter(self%cVSetDiat,'cVSetDiat',default=cVSetDiat) !,scale_factor=1.0_rk/secs_pr_day)scale it after equation
   call self%get_parameter(self%cVSetGren,'cVSetGren',default=cVSetGren) !,scale_factor=1.0_rk/secs_pr_day)scale it after equation
   call self%get_parameter(self%cVSetBlue ,'cVSetBlue ',default=cVSetBlue) ! ,scale_factor=1.0_rk/secs_pr_day)scale it after equation
   call self%get_parameter(self%cRhoIM,'cRhoIM',default=cRhoIM)
   call self%get_parameter(self%cRhoOM,'cRhoOM',default=cRhoOM)
   call self%get_parameter(self%fDOrgSoil,'fDOrgSoil',default=fDOrgSoil)
   call self%get_parameter(self%cPO4Ground,'cPO4Ground',default=cPO4Ground)
   call self%get_parameter(self%cNH4Ground,'cNH4Ground',default=cNH4Ground)
   call self%get_parameter(self%cNO3Ground,'cNO3Ground',default=cNO3Ground)
   call self%get_parameter(self%cTmOptFish,'cTmOptFish',default=cTmOptFish)
   call self%get_parameter(self%cSigTmFish,'cSigTmFish',default=cSigTmFish)
!   call self%get_parameter(self%cSiDDiat,'cSiDDiat',default=cSiDDiat)
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
!  for loadings, in day unit
! scale_factor=1.0_rk/secs_pr_day only kept for loading and transporting
   call self%get_parameter(self%cLoadPO4,'cLoadPO4',default=cLoadPO4,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cLoadNO3,'cLoadNO3',default=cLoadNO3,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%uQIn,'uQIn',default=uQIn,scale_factor=1.0_rk/secs_pr_day)
!  for erosion
   call self%get_parameter(self%fSedErosIM,'fSedErosIM',default=fSedErosIM)
   call self%get_parameter(self%cDErosTot,'cDErosTot',default=cDErosTot,scale_factor=1.0_rk/secs_pr_day)
!  parameters for shear stress related resuspension
   call self%get_parameter(self%crt_shear,'crt_shear',default=crt_shear)
   call self%get_parameter(self%ref_shear,'ref_shear',default=ref_shear)
   call self%get_parameter(self%alpha,'alpha',default=alpha,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%eta,'eta',default=eta)
   call self%get_parameter(self%cVSetMain,'cVSetMain',default=cVSetMain,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%resusp_meth,'resusp_meth',default=resusp_meth)
!  regirster state variables dependencies
!  step1, Register dependencies on external state variables
!  Register dependencies to abiotic water module
   call self%register_state_dependency(self%id_SWNH4, 'Amonia_pool_in_water','g/m**3','Amonia_pool_in_water')
   call self%register_state_dependency(self%id_SWNO3, 'Nitrates_pool_in_water','g/m**3','Nitrates_pool_in_water')
   call self%register_state_dependency(self%id_SWPO4, 'Phosphate_pool_in_water','g/m**3','Phosphate_pool_in_water')
   call self%register_state_dependency(self%id_SWPAIM, 'Absorbed_phosphrus_in_water','g/m**3','Absorbed_phosphrus_in_water')
   call self%register_state_dependency(self%id_SWO2, 'Oxygen_pool_in_water','g/m**3','Oxygen_pool_in_water')
   call self%register_state_dependency(self%id_SWDIM, 'Inorg_pool_in_water','g/m**3','Inorg_pool_in_water')
   call self%register_state_dependency(self%id_SWDDet, 'Detritus_DW_in_water','g/m**3','Detritus_DW_in_water')
   call self%register_state_dependency(self%id_SWNDet, 'Detritus_N_in_water','g/m**3','Detritus_N_in_water')
   call self%register_state_dependency(self%id_SWPDet, 'Detritus_P_in_water','g/m**3','Detritus_P_in_water')
   call self%register_state_dependency(self%id_SWSiDet, 'Detritus_Si_in_water','g/m**3','Detritus_Si_in_water')
   call self%register_state_dependency(self%id_SWSiO2, 'SiO2_pool_water','g/m**3','SiO2_pool_water')
   
!  Register dependencies to abiotic sediment module
   call self%register_state_dependency(self%id_WSNH4, 'Amonia_pool_in_sediment','g/m**2','Amonia_pool_in_sediment')
   call self%register_state_dependency(self%id_WSNO3, 'Nitrates_pool_in_sediment','g/m**2','Nitrates_pool_in_sediment')
   call self%register_state_dependency(self%id_WSPO4, 'Phosphate_pool_in_sediment','g/m**2','Phosphate_pool_in_sediment')
   call self%register_state_dependency(self%id_WSPAIM, 'Absorbed_phosphrus_in_sediment','g/m**2','Absorbed_phosphrus_in_sediment')
   call self%register_state_dependency(self%id_WSDIM, 'Inorg_pool_in_sediment','g/m**2','Inorg_pool_in_sediment')
   call self%register_state_dependency(self%id_WSDDet, 'Detritus_DW_in_sediment','g/m**2','Detritus_DW_in_sediment')
   call self%register_state_dependency(self%id_WSNDet, 'Detritus_N_in_sediment','g/m**2','Detritus_N_in_sediment')
   call self%register_state_dependency(self%id_WSPDet, 'Detritus_P_in_sediment','g/m**2','Detritus_P_in_sediment')
   call self%register_state_dependency(self%id_WSSiDet, 'Detritus_Si_in_sediment','g/m**2','Detritus_Si_in_sediment')
   call self%register_state_dependency(self%id_WSDHum, 'Humus_DW_in_sediment','g/m**2','Humus_DW_in_sediment')
   call self%register_state_dependency(self%id_WSNHum, 'Humus_N_in_sediment','g/m**2','Humus_N_in_sediment')
   call self%register_state_dependency(self%id_WSPHum, 'Humus_P_in_sediment','g/m**2','Humus_P_in_sediment')
!  Register dependencies to phytoplankton in water column
   call self%register_state_dependency(self%id_SWDDiat, 'Diatom_DW_in_water','g/m**3','Diatom_DW_in_water')
   call self%register_state_dependency(self%id_SWDGren, 'Green_DW_in_water','g/m**3','Green_DW_in_water')
   call self%register_state_dependency(self%id_SWDBlue, 'Blue_DW_in_water','g/m**3','Blue_DW_in_water')
   call self%register_state_dependency(self%id_SWNDiat, 'Diatom_N_in_water','g/m**3','Diatom_N_in_water')
   call self%register_state_dependency(self%id_SWNGren, 'Green_N_in_water','g/m**3','Green_N_in_water')
   call self%register_state_dependency(self%id_SWNBlue, 'Blue_N_in_water','g/m**3','Blue_N_in_water')
   call self%register_state_dependency(self%id_SWPDiat, 'Diatom_P_in_water','g/m**3','Diatom_P_in_water')
   call self%register_state_dependency(self%id_SWPGren, 'Green_P_in_water','g/m**3','Green_P_in_water')
   call self%register_state_dependency(self%id_SWPBlue, 'Blue_P_in_water','g/m**3','Blue_P_in_water')
!   call self%register_state_dependency(self%id_SWSiDiat, 'Diatom_Si_in_water','g/m**3','Diatom_Si_in_water')
!  Register dependencies to phytoplankton in sediment
   call self%register_state_dependency(self%id_WSDDiat, 'Diatom_DW_in_sediment','g/m**2','Diatom_DW_in_sediment')
   call self%register_state_dependency(self%id_WSDGren, 'Green_DW_in_sediment','g/m**2','Green_DW_in_sediment')
   call self%register_state_dependency(self%id_WSDBlue, 'Blue_DW_in_sediment','g/m**2','Blue_DW_in_sediment')
   call self%register_state_dependency(self%id_WSNDiat, 'Diatom_N_in_sediment','g/m**2','Diatom_N_in_sediment')
   call self%register_state_dependency(self%id_WSNGren, 'Green_N_in_sediment','g/m**2','Green_N_in_sediment')
   call self%register_state_dependency(self%id_WSNBlue, 'Blue_N_in_sediment','g/m**2','Blue_N_in_sediment')
   call self%register_state_dependency(self%id_WSPDiat, 'Diatom_P_in_sediment','g/m**2','Diatom_P_in_sediment')
   call self%register_state_dependency(self%id_WSPGren, 'Green_P_in_sediment','g/m**2','Green_P_in_sediment')
   call self%register_state_dependency(self%id_WSPBlue, 'Blue_P_in_sediment','g/m**2','Blue_P_in_sediment')
!   call self%register_state_dependency(self%id_WSSiDiat, 'Diatom_Si_in_sediment','g/m**2','Diatom_Si_in_sediment')
!  register vegetation and fish for resuspension dependency
   call self%register_state_dependency(self%id_DragVeg, 'vegetation_DW','g/m**2','vegetation_DW')
   call self%register_state_dependency(self%id_TurbFish, 'adult_fish_DW','g/m**3','adult_fish_DW')
!  register zooplankton for transport purpose
   call self%register_state_dependency(self%id_DTranZoo, 'Zooplankton_DW','g/m**3','Zooplankton_DW')
   call self%register_state_dependency(self%id_PTranZoo, 'Zooplankton_P','g/m**3','Zooplankton_P')
   call self%register_state_dependency(self%id_NTranZoo, 'Zooplankton_N','g/m**3','Zooplankton_N')
!  step 2, Automatically couple dependencies if target variables have been specified.
!  Register dependencies to abiotic water module
   if (Amonia_pool_in_water/='') call self%request_coupling(self%id_SWNH4,  Amonia_pool_in_water)
   if (Nitrates_pool_in_water/='') call self%request_coupling(self%id_SWNO3,Nitrates_pool_in_water)
   if (Phosphate_pool_in_water/='') call self%request_coupling(self%id_SWPO4,Phosphate_pool_in_water)
   if (Absorbed_phosphrus_in_water/='') call self%request_coupling(self%id_SWPAIM,Absorbed_phosphrus_in_water)
   if (Oxygen_pool_in_water/='') call self%request_coupling(self%id_SWO2,Oxygen_pool_in_water)
   if (Inorg_pool_in_water/='') call self%request_coupling(self%id_SWDIM,Inorg_pool_in_water)
   if (Detritus_DW_in_water/='') call self%request_coupling(self%id_SWDDet,Detritus_DW_in_water)
   if (Detritus_N_in_water/='') call self%request_coupling(self%id_SWNDet,Detritus_N_in_water)
   if (Detritus_P_in_water/='') call self%request_coupling(self%id_SWPDet,Detritus_P_in_water)
   if (Detritus_Si_in_water/='') call self%request_coupling(self%id_SWSiDet,Detritus_Si_in_water)
   if (SiO2_pool_water/='') call self%request_coupling(self%id_SWSiO2,SiO2_pool_water)
!  Register dependencies to abiotic sediment module
   if (Amonia_pool_in_sediment/='') call self%request_coupling(self%id_WSNH4,  Amonia_pool_in_sediment)
   if (Nitrates_pool_in_sediment/='') call self%request_coupling(self%id_WSNO3,Nitrates_pool_in_sediment)
   if (Phosphate_pool_in_sediment/='') call self%request_coupling(self%id_WSPO4,Phosphate_pool_in_sediment)
   if (Absorbed_phosphrus_in_sediment/='') call self%request_coupling(self%id_WSPAIM,Absorbed_phosphrus_in_sediment)
   if (Inorg_pool_in_sediment/='') call self%request_coupling(self%id_WSDIM,Inorg_pool_in_sediment)
   if (Detritus_DW_in_sediment/='') call self%request_coupling(self%id_WSDDet,Detritus_DW_in_sediment)
   if (Detritus_N_in_sediment/='') call self%request_coupling(self%id_WSNDet,Detritus_N_in_sediment)
   if (Detritus_P_in_sediment/='') call self%request_coupling(self%id_WSPDet,Detritus_P_in_sediment)
   if (Detritus_Si_in_sediment/='') call self%request_coupling(self%id_WSSiDet,Detritus_Si_in_sediment)
   if (Humus_DW_in_sediment/='') call self%request_coupling(self%id_WSDHum,Humus_DW_in_sediment)
   if (Humus_N_in_sediment/='') call self%request_coupling(self%id_WSNHum,Humus_N_in_sediment)
   if (Humus_P_in_sediment/='') call self%request_coupling(self%id_WSPHum,Humus_P_in_sediment)
!  Register dependencies to phytoplankton in water column
   if (Diatom_DW_in_water/='') call self%request_coupling(self%id_SWDDiat,Diatom_DW_in_water)
   if (Green_DW_in_water/='') call self%request_coupling(self%id_SWDGren,Green_DW_in_water)
   if (Blue_DW_in_water/='') call self%request_coupling(self%id_SWDBlue,Blue_DW_in_water)
   if (Diatom_N_in_water/='') call self%request_coupling(self%id_SWNDiat,Diatom_N_in_water)
   if (Green_N_in_water/='') call self%request_coupling(self%id_SWNGren,Green_N_in_water)
   if (Blue_N_in_water/='') call self%request_coupling(self%id_SWNBlue,Blue_N_in_water)
   if (Diatom_P_in_water/='') call self%request_coupling(self%id_SWPDiat,Diatom_P_in_water)
   if (Green_P_in_water/='') call self%request_coupling(self%id_SWPGren,Green_P_in_water)
   if (Blue_P_in_water/='') call self%request_coupling(self%id_SWPBlue,Blue_P_in_water)
!   if (Diatom_Si_in_water/='') call self%request_coupling(self%id_SWSiDiat,Diatom_Si_in_water)
!  Register dependencies to phytoplankton in sediment
   if (Diatom_DW_in_sediment/='') call self%request_coupling(self%id_WSDDiat,Diatom_DW_in_sediment)
   if (Green_DW_in_sediment/='') call self%request_coupling(self%id_WSDGren,Green_DW_in_sediment)
   if (Blue_DW_in_sediment/='') call self%request_coupling(self%id_WSDBlue,Blue_DW_in_sediment)
   if (Diatom_N_in_sediment/='') call self%request_coupling(self%id_WSNDiat,Diatom_N_in_sediment)
   if (Green_N_in_sediment/='') call self%request_coupling(self%id_WSNGren,Green_N_in_sediment)
   if (Blue_N_in_sediment/='') call self%request_coupling(self%id_WSNBlue,Blue_N_in_sediment)
   if (Diatom_P_in_sediment/='') call self%request_coupling(self%id_WSPDiat,Diatom_P_in_sediment)
   if (Green_P_in_sediment/='') call self%request_coupling(self%id_WSPGren,Green_P_in_sediment)
   if (Blue_P_in_sediment/='') call self%request_coupling(self%id_WSPBlue,Blue_P_in_sediment)
!   if (Diatom_Si_in_sediment/='') call self%request_coupling(self%id_WSSiDiat,Diatom_Si_in_sediment)
!  register vegetation and fish for resuspension dependency
   if (vegetation_DW/='') call self%request_coupling(self%id_DragVeg,vegetation_DW)
   if (adult_fish_DW/='') call self%request_coupling(self%id_TurbFish,adult_fish_DW)
!  register zooplankton for transport purpose
   if (Zooplankton_DW/='') call self%request_coupling(self%id_DTranZoo,Zooplankton_DW)
   if (Zooplankton_P/='') call self%request_coupling(self%id_PTranZoo,Zooplankton_P)
   if (Zooplankton_N/='') call self%request_coupling(self%id_NTranZoo,Zooplankton_N)
   


!  register environmental dependencies
   call self%register_dependency(self%id_uTm,standard_variables%temperature)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
!   call self%register_dependency(self%id_dz,standard_variables%cell_thickness)
   call self%register_dependency(self%id_Day,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_shear,standard_variables%bottom_stress)
!  register diagnostic dependencies, 2 steps
!  step1, Register dependencies on external diagnostic variables
   call self%register_dependency(self%id_tDAbioDetS, 'Detritus_abiotic_update','g/m**2/s','Detritus_abiotic_update')
   call self%register_dependency(self%id_tDAbioHumS, 'Humus_abiotic_update','g/m**2/s','Humus_abiotic_update')
   call self%register_dependency(self%id_tDPrimDetS, 'Detritus_from_algae','--','Detritus_from_algae')
   call self%register_dependency(self%id_tDWebDetS, 'Detritus_from_foodweb','--','Detritus_from_foodweb')
   call self%register_dependency(self%id_tDBedDetS, 'Detritus_from_vegetation','--','Detritus_from_vegetation') 
     
!  step 2, Automatically couple dependencies if target variables have been specified.
   if (Detritus_abiotic_update/='') call self%request_coupling(self%id_tDAbioDetS,Detritus_abiotic_update)
   if (Humus_abiotic_update/='') call self%request_coupling(self%id_tDAbioHumS,Humus_abiotic_update)
   if (Detritus_from_algae/='') call self%request_coupling(self%id_tDPrimDetS,Detritus_from_algae)
   if (Detritus_from_foodweb/='') call self%request_coupling(self%id_tDWebDetS,Detritus_from_foodweb)
   if (Detritus_from_vegetation/='') call self%request_coupling(self%id_tDBedDetS,Detritus_from_vegetation)
!  
   call self%register_diagnostic_variable(self%id_tDBurIM,'tDBurIM','g/m**2/s','tDBurIM',           &
                 output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_shearstress,'shearstress','N/m**2','shearstress',           &
                 output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_aFunDimSusp,'aFunDimSusp','--','aFunDimSusp',           &
                 output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_tDResusDead,'tDResusDead','--','tDResusDead',           &
                 output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_aFunTauSet,'aFunTauSet','--','aFunTauSet',           &
                 output=output_time_step_averaged)
   


   return

99  call self%fatal_error('au_pclake_auxilary_init','Error reading namelist au_pclake_auxilary')

100 call self%fatal_error('au_pclake_auxilary_init','Namelist au_pclake_auxilary was not found.')

   end subroutine initialize
!
!EOC

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: 
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_au_pclake_auxilary), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
! !LOCAL VARIABLES:
! !carriers for state dependencies in different modules
!  in abiotic water column module
   real(rk)                   :: sDIMW,sDDetW,sNDetW,sPDetW,sSiDetW
   real(rk)                   :: sPO4W,sPAIMW,sNH4W,sNO3W
!  in phytoplankton water column module
   real(rk)                   :: sDDiatW,sDGrenW,sDBlueW
   real(rk)                   :: sNDiatW,sNGrenW,sNBlueW
   real(rk)                   :: sPDiatW,sPGrenW,sPBlueW
!  in abiotic water column module
   real(rk)                   :: sDIMS,sDDetS,sNDetS,sPDetS,sSiDetS
   real(rk)                   :: sPO4S,sPAIMS,sNH4S,sNO3S
   real(rk)                   :: sDHumS,sNHumS,sPHumS
!  in phytoplankton  sediment module
   real(rk)                   :: sDDiatS,sDGrenS,sDBlueS
   real(rk)                   :: sNDiatS,sNGrenS,sNBlueS
   real(rk)                   :: sPDiatS,sPGrenS,sPBlueS
!  in vegetation module
   real(rk)                   :: sDVeg
!  in foodweb water column module
   real(rk)                   :: sDFiAd
!  carriers for environmental dependencies
   real(rk)                   :: uTm,sDepthW !,dz  ! This is the depth for empirical suspension function, should be bottom_depth
!  carriers for diagnostic dependencies
   real(rk)                   :: tDAbioHumS
   real(rk)                   :: tDAbioDetS,tDPrimDetS,tDWebDetS,tDBedDetS
!  variables for nutrient rations for detritus in sediment
   real(rk)                   :: rPDDetS,rNDDetS,rSiDDetS
   real(rk)                   :: rPDHumS,rNDHumS
!  variables for nutrient rations for phytoplankton in water
   real(rk)                   :: rNDDiatW,rNDGrenW,rNDBlueW
   real(rk)                   :: rPDDiatW,rPDGrenW,rPDBlueW
!  variables for nutrient rationsfor phytoplankton in sediment
   real(rk)                   :: rNDDiatS,rNDGrenS,rNDBlueS
   real(rk)                   :: rPDDiatS,rPDGrenS,rPDBlueS
!  temperature related variables
   real(rk)                   :: uFunTmSet,uFunTmFish
!  variables related to resuspension(in the order of apperance)
   real(rk)                   :: aFunVegResus,tDTurbFish
   real(rk)                   :: aFunDimSusp,tDResusTauDead,tDResusBareDead
   real(rk)                   :: tDResusDead,tDResusIM,tDResusDet,tPResusDet
   real(rk)                   :: tNResusDet,tSiResusDet,tPResusPO4,tPResusAIM
   real(rk)                   :: tNResusNO3,tNResusNH4
!  variables for phytoplankton resuspension
   real(rk)                   :: akResusPhytRef
   real(rk)                   :: tDResusDiat,tDResusGren,tDResusBlue
   real(rk)                   :: tNResusDiat,tNResusGren,tNResusBlue
   real(rk)                   :: tPResusDiat,tPResusGren,tPResusBlue  !,tSiResusDiat
!  variables related to sedimentation(in the order of apperance)
   real(rk)                   :: aFunTauSet
   real(rk)                   :: uCorVSetIM,tDSetIM,tPSetAIM
   real(rk)                   :: uCorVSetDet,tDSetDet,tPSetDet,tNSetDet,tSiSetDet
!  variables for phytoplankton sedimentation
   real(rk)                   :: uCorVSetDiat,uCorVSetGren,uCorVSetBlue
   real(rk)                   :: tDSetDiat,tDSetGren,tDSetBlue
   real(rk)                   :: tNSetDiat,tNSetGren,tNSetBlue
   real(rk)                   :: tPSetDiat,tPSetGren,tPSetBlue !,tSiSetDiat
!  Variables related to burial process(ub the order of appearance)
   real(rk)                   :: tDIMS,tDDetS,vDeltaS
   real(rk)                   :: tDBurIM,tDBurDet,tPBurDet,tPBurAIM,tPBurPO4
   real(rk)                   :: tNBurDet,tNBurNH4,tNBurNO3,tSiBurDet
!  Humus variables
   real(rk)                   :: tDHumS,tDBurHum,tDBurOM,tNBurHum,tPBurHum
!  erosion variables
   real(rk)                   :: uDErosIM,uDErosIMW,uDErosIMS
   real(rk)                   :: uDErosOM,uPErosOM,uNErosOM
!  variables of new resuspension method
   real(rk)                   :: shear

!EOP
!-----------------------------------------------------------------------
!BOC
!------------------------------------------------------------------------
!  Spatial loop
   _FABM_HORIZONTAL_LOOP_BEGIN_
!-----------------------------------------------------------------------
!  Retrieve dependencis value
!-----------------------------------------------------------------------
! !Retrieve state dependencie value
!  from abiotic water column
   _GET_(self%id_SWNH4,sNH4W)
   _GET_(self%id_SWNO3,sNO3W)
   _GET_(self%id_SWPO4,sPO4W)
   _GET_(self%id_SWPAIM,sPAIMW)
   _GET_(self%id_SWDIM,sDIMW)
   _GET_(self%id_SWDDet,sDDetW)
   _GET_(self%id_SWNDet,sNDetW)
   _GET_(self%id_SWPDet,sPDetW)
   _GET_(self%id_SWSiDet,sSiDetW)
!  from phytoplankton in water column
   _GET_(self%id_SWDDiat,sDDiatW)
   _GET_(self%id_SWDGren,sDGrenW)
   _GET_(self%id_SWDBlue,sDBlueW)
   _GET_(self%id_SWNDiat,sNDiatW)
   _GET_(self%id_SWNGren,sNGrenW)
   _GET_(self%id_SWNBlue,sNBlueW)
   _GET_(self%id_SWPDiat,sPDiatW)
   _GET_(self%id_SWPGren,sPGrenW)
   _GET_(self%id_SWPBlue,sPBlueW)
!   _GET_(self%id_SWSiDiat,sSiDiatW)
!  from abiotic sediment
   _GET_HORIZONTAL_(self%id_WSNH4,sNH4S)
   _GET_HORIZONTAL_(self%id_WSNO3,sNO3S)
   _GET_HORIZONTAL_(self%id_WSPO4,sPO4S)
   _GET_HORIZONTAL_(self%id_WSPAIM,sPAIMS)
   _GET_HORIZONTAL_(self%id_WSDIM,sDIMS)
   _GET_HORIZONTAL_(self%id_WSDDet,sDDetS)
   _GET_HORIZONTAL_(self%id_WSNDet,sNDetS)
   _GET_HORIZONTAL_(self%id_WSPDet,sPDetS)
   _GET_HORIZONTAL_(self%id_WSSiDet,sSiDetS)
   _GET_HORIZONTAL_(self%id_WSDHum,sDHumS)
   _GET_HORIZONTAL_(self%id_WSNHum,sNHumS)
   _GET_HORIZONTAL_(self%id_WSPHum,sPHumS)
!  from phytoplankton in sediment
   _GET_HORIZONTAL_(self%id_WSDDiat,sDDiatS)
   _GET_HORIZONTAL_(self%id_WSDGren,sDGrenS)
   _GET_HORIZONTAL_(self%id_WSDBlue,sDBlueS)
   _GET_HORIZONTAL_(self%id_WSNDiat,sNDiatS)
   _GET_HORIZONTAL_(self%id_WSNGren,sNGrenS)
   _GET_HORIZONTAL_(self%id_WSNBlue,sNBlueS)
   _GET_HORIZONTAL_(self%id_WSPDiat,sPDiatS)
   _GET_HORIZONTAL_(self%id_WSPGren,sPGrenS)
   _GET_HORIZONTAL_(self%id_WSPBlue,sPBlueS)
!  vegatation influence on vegetation
   _GET_HORIZONTAL_(self%id_DragVeg,sDVeg)
!  fish predation influence on resuspension
   _GET_(self%id_TurbFish,sDFiAd)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_HORIZONTAL_(self%id_shear,shear)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
!  fish biomass converted to g/m^2
   sDFiAd=sDFiAd*sDepthW
!  retrieve diagnostic denpendency
   _GET_HORIZONTAL_(self%id_tDAbioDetS,tDAbioDetS)
   _GET_HORIZONTAL_(self%id_tDPrimDetS,tDPrimDetS)
   _GET_HORIZONTAL_(self%id_tDWebDetS,tDWebDetS)
   _GET_HORIZONTAL_(self%id_tDBedDetS,tDBedDetS)
   _GET_HORIZONTAL_(self%id_tDAbioHumS,tDAbioHumS)
!-----------------------------------------------------------------------
!  Current nutrients ratios(check the curent state)
!-----------------------------------------------------------------------
   rPDDetS=sPDetS/(sDDetS+NearZero)
   rNDDetS=sNDetS/(sDDetS+NearZero)
   rSiDDetS=sSiDetS/(sDDetS+NearZero)
   rPDHumS=sPHumS/(sDHumS+NearZero)
   rNDHumS=sNHumS/(sDHumS+NearZero)
! !external source status
!  for phytoplankton in water
   rPDDiatW = sPDiatW /(sDDiatW+NearZero)
   rNDDiatW = sNDiatW /(sDDiatW+NearZero)
   rPDGrenW = sPGrenW /(sDGrenW+NearZero)
   rNDGrenW = sNGrenW /(sDGrenW+NearZero)
   rPDBlueW = sPBlueW /(sDBlueW+NearZero)
   rNDBlueW = sNBlueW /(sDBlueW+NearZero)
!  for phytoplankton in sediment
   rPDDiatS = sPDiatS /(sDDiatS+NearZero)
   rNDDiatS = sNDiatS /(sDDiatS+NearZero)
   rPDGrenS = sPGrenS /(sDGrenS+NearZero)
   rNDGrenS = sNGrenS /(sDGrenS+NearZero)
   rPDBlueS = sPBlueS /(sDBlueS+NearZero)
   rNDBlueS = sNBlueS /(sDBlueS+NearZero)

!-----------------------------------------------------------------------
!  Temperature functions for sediment abiotic process
!-----------------------------------------------------------------------
!  temperature_correction_of_sedimentation
   uFunTmSet= uFunTmAbio(uTm,self%cThetaSet)
   uFunTmFish= uFunTmBio(uTm,self%cSigTmFish,self%cTmOptFish)
!-----------------------------------------------------------------------
!  Process related to other modules
!-----------------------------------------------------------------------
!  vegetation_dependence_of_resuspension
   aFunVegResus=max(0.0_rk,1.0_rk-self%kVegResus*sDVeg)
!------------------------------------------------------------------------------------------------------------
!  resuspension and sedimentation(PCLake method)
!------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------
!   The resuspended matter in the water column calculation
!------------------------------------------------------------------------------------------------------------
!  bioturbation_by_fish
      tDTurbFish=self%kTurbFish*uFunTmFish*sDFiAd
!  calculate resuspension rate, two methods
   select case(self%resusp_meth)
      case(1)
!         Empirical_suspended_matter_function_(logistic_fit_to_data), in day
          if (uTm >= 0.1_rk) then
              aFunDimSusp=self%cSuspRef*((self%cSuspMin+self%cSuspMax/(1.0_rk+exp(self%cSuspSlope*&
              & (sDepthW-self%hDepthSusp))))*((((self%cFetch +NearZero)/ self%cFetchRef) )** (0.5_rk)))
          else
              aFunDimSusp=0.0_rk
          endif
      case(2)
         if(shear <= self%crt_shear)then
            aFunDimSusp=0.0_rk
         else
            aFunDimSusp=(self%alpha/self%cVSetMain)*((shear-self%crt_shear)/self%ref_shear)**self%eta
         endif
   end select
   tDResusTauDead=min(aFunDimSusp, ((aFunDimSusp +NearZero )**(0.5_rk))) &
   &*((self%fLutum/ self%fLutumRef )** (0.5_rk))*self%bPorS
!  resuspension_due_to_shear_stress_AND_fish, in day
   tDResusBareDead=tDResusTauDead+tDTurbFish
!  resuspension,_corrected_for_vegetation_effect, in secs
!   tDResusDead=tDResusBareDead*aFunVegResus
   tDResusDead=tDResusBareDead*aFunVegResus/secs_pr_day
!------------------------------------------------------------------------------------------------------------
!  Different matter resuspension based on the suspended matter concentration in the water column
!------------------------------------------------------------------------------------------------------------
!  The inorganic matter resuspension
   tDResusIM=self%fLutum*sDIMS/(self%fLutum*sDIMS+sDDetS)*tDResusDead
!  detrital_resuspension_DW
   tDResusDet=sDDetS/(self%fLutum*sDIMS+sDDetS)*tDResusDead
!  detrital_resuspension_P
   tPResusDet=rPDDetS*tDResusDet
!  detrital_resuspension_N
   tNResusDet=rNDDetS*tDResusDet
!  detrital_resuspension_SI
   tSiResusDet=rSiDDetS*tDResusDet
!  resuspension_nutrient_P
   tPResusPO4=sPO4S/sDDetS*tDResusDet
!  resuspension_absorbed_PAIM
   tPResusAIM=sPAIMS/sDIMS*tDResusIM
!  resuspension_nutrient_NO3
   tNResusNO3=sNO3S/sDDetS*tDResusDet
!  resuspension_nutrient_NH4
   tNResusNH4=sNH4S/sDDetS*tDResusDet

!  convert the seconds rate to daily rate, for the equation purpose
   tDResusDead=tDResusDead*secs_pr_day
!  phytoplankton_resuspension_rate_constant, in day
   akResusPhytRef = self%kResusPhytMax * (1.0_rk - exp(self%cResusPhytExp * tDResusDead))
!  convert to secs
   akResusPhytRef=akResusPhytRef/secs_pr_day
!!  Algae group resuspension
!  reuspension of Diatom,DW
   tDResusDiat=akResusPhytRef*sDDiatS
!  reuspension of Green algae,DW
   tDResusGren=akResusPhytRef*sDGrenS
!  reuspension of Blue algae,DW
   tDResusBlue=akResusPhytRef*sDBlueS
!  reuspension of Diatom,N
   tNResusDiat = rNDDiatS * tDResusDiat
!  reuspension of Green algae,N
   tNResusGren = rNDGrenS * tDResusGren
!  reuspension of Blue algae,N
   tNResusBlue = rNDBlueS * tDResusBlue
!  reuspension of Diatom,P
   tPResusDiat = rPDDiatS * tDResusDiat
!  reuspension of Green algae,P
   tPResusGren = rPDGrenS * tDResusGren
!  reuspension of Blue algae,P
   tPResusBlue = rPDBlueS * tDResusBlue
!  Diatoms_sedimentation
!   tSiResusDiat = self%cSiDDiat * tDResusDiat
!-----------------------------------------------------------------------
!  The sedimentation calculation, based on resuspension
!-----------------------------------------------------------------------
!  correction_factor_for_settling_rate_(<=_1),basic settling rate, in day
   aFunTauSet=min(1.0_rk,1.0_rk/((aFunDimSusp +NearZero )**(0.5_rk)))
!-----------------------------------------------------------------------
!  Different matter sedimentation based on the basic settling rate
!-----------------------------------------------------------------------
!  sedimentation_velocity_of_IM, in day
   uCorVSetIM=aFunTauSet*((self%fLutumRef/self%fLutum)**(0.5))*uFunTmSet*self%cVSetIM
!  convert to seconds
   uCorVSetIM=uCorVSetIM/secs_pr_day
!  sedimentation_IM
   tDSetIM=uCorVSetIM*sDIMW
!  sedimentation_PAIM
   tPSetAIM=sPAIMW/(sDIMW +NearZero)*tDSetIM
!  sedimentation_velocity_of_detritus, in day
   uCorVSetDet=self%cVSetDet*aFunTauSet*uFunTmSet
!  convert to seconds
   uCorVSetDet=uCorVSetDet/secs_pr_day
!  sedimentation_flux_of_detritus
   tDSetDet=uCorVSetDet*sDDetW
!  sedimentation_detrital_P
   tPSetDet=uCorVSetDet*sPDetW
!  sedimentation_detrital_N
   tNSetDet=uCorVSetDet*sNDetW
!  sedimentation_detrital_Si
   tSiSetDet=uCorVSetDet*sSiDetW
!  corrected_sedimentation_velocity_of_Algae, in day
   uCorVSetDiat = self%cVSetDiat * aFunTauSet * uFunTmSet
!  convert to seconds
   uCorVSetDiat=uCorVSetDiat/secs_pr_day
!  sedimentation_flux_of_Diatom
   tDSetDiat = uCorVSetDiat * sDDiatW
!  corrected_sedimentation_velocity_of_Algae,in day
   uCorVSetGren = self%cVSetGren * aFunTauSet * uFunTmSet
!  convert to seconds
   uCorVSetGren=uCorVSetGren/secs_pr_day
!  sedimentation_flux_of_Algae
   tDSetGren = uCorVSetGren * sDGrenW
!  corrected_sedimentation_velocity_of_Algae, in day
   uCorVSetBlue = self%cVSetBlue * aFunTauSet * uFunTmSet
!  convert to seconds
   uCorVSetBlue=uCorVSetBlue/secs_pr_day
!  sedimentation_flux_of_Algae
   tDSetBlue = uCorVSetBlue * sDBlueW
!  sedimentation
   tNSetDiat = rNDDiatW * tDSetDiat
!  sedimentation
   tNSetGren = rNDGrenW * tDSetGren
!  sedimentation
   tNSetBlue = rNDBlueW * tDSetBlue
!  sedimentation
   tPSetDiat = rPDDiatW * tDSetDiat
!  sedimentation
   tPSetGren = rPDGrenW * tDSetGren
!  sedimentation
   tPSetBlue = rPDBlueW * tDSetBlue
!  Diatoms_sedimentation
!   tSiSetDiat = self%cSiDDiat * tDSetDiat
!-----------------------------------------------------------------------
!  Burial of sediment,contains erosion process
!-----------------------------------------------------------------------
!  IM_input_from_banks
   uDErosIM = (1.0 - self%fDOrgSoil) * self%cDErosTot
!  IM_input_to_sediment_from_banks
   uDErosIMS = self%fSedErosIM * uDErosIM
!  IM_input_to_water_column_from_banks
   uDErosIMW = uDErosIM - uDErosIMS
!  organic_matter_input_from_banks
   uDErosOM = self%fDOrgSoil * self%cDErosTot
!  organic_P_input_from_banks
   uPErosOM = 0.001_rk * uDErosOM  ! cPDSoilOM=0.001
!  organic_N_input_from_banks
   uNErosOM = 0.01_rk * uDErosOM  ! cNDSoilOM=0.01
!  increase_in_inorganic_matter_in_sediment
   ! original form looks like
   !tDIMS = tDAbioIMS
   ! tDAbioIMS = uDErosIMS + tDSetIM - tDResusIM
! due to the change for fabm, then:
   tDIMS= uDErosIMS + tDSetIM - tDResusIM
!  increase_in_sediment_humus_in_lake
!   tDAbioHumS = uDErosOM + fRefrDetS * tDMinDetS - tDMinHumS
!   and uDErosOM is calculated here
    tDHumS = uDErosOM+tDAbioHumS
! increase_in_sediment_detritus_in_lake
   ! Original form looks like:
   !tDDetS = tDAbioDetS   ! + tDPrimDetS + tDWebDetS + tDBedDetS
   !tDAbioDetS = tDSetDet - tDResusDet - tDMinDetS
   tDDetS= tDSetDet - tDResusDet+tDAbioDetS+ tDPrimDetS + tDWebDetS + tDBedDetS
!  turnover_depth_in_lake
   vDeltaS = (tDIMS / self%cRhoIM + (tDHumS + tDDetS) / self%cRhoOM)/(1.0_rk - self%bPorS)
!  burial_flux_of_DW_in_inorganic_matter_in_lake
   if (vDeltaS >= 0.0_rk) then
   tDBurIM = ((tDHumS + tDDetS) +(self%cRhoOM / self%cRhoIM) * tDIMS) / ((sDHumS + sDDetS) /sDIMS &
   & + self%cRhoOM / self%cRhoIM)
   else
   tDBurIM = ( (tDHumS + tDDetS) +(self%cRhoOM / self%cRhoIM) * tDIMS) / (self%fDOrgSoil &
   &/(1.0_rk - self%fDOrgSoil) + self%cRhoOM / self%cRhoIM) 
   endif
   
!  burial_flux_of_DW_in_organic_matter_in_lake
   if (vDeltaS >= 0.0) then
      tDBurOM = (sDHumS + sDDetS) / sDIMS * tDBurIM 
   else
      tDBurOM = self%fDOrgSoil /(1.0 - self%fDOrgSoil) * tDBurIM 
   endif
   
!  burial_flux_of_DW_in_detritus_in_lake
   if (vDeltaS >= 0.0) then
      tDBurDet = sDDetS /(sDHumS + sDDetS) * tDBurOM 
   else
      tDBurDet = 0.0 
   endif

!  burial_flux_of_P_in_detritus_in_lake
   if (vDeltaS >= 0.0_rk) then
   tPBurDet = rPDDetS * tDBurDet
   else
   tPBurDet = 0.0_rk 
   endif
  
!  burial_flux_of_P_absorbed_onto_inorganic_matter_in_lake
   if (vDeltaS >= 0.0_rk) then
   tPBurAIM = sPAIMS / sDIMS * tDBurIM 
   else
   tPBurAIM = 0.0_rk
   endif
  
!  burial_flux_of_dissolved_P_in_lake
   if (vDeltaS >= 0.0_rk) then
   tPBurPO4 = sPO4S *(vDeltaS / self%cDepthS) 
   else
   tPBurPO4 = self%cPO4Ground *(self%bPorS * vDeltaS) 
   endif
  
!  burial_flux_of_N_in_detritus_in_lake
   if (vDeltaS >= 0.0_rk) then
   tNBurDet =rNDDetS * tDBurDet 
   else
   tNBurDet = 0.0_rk
   endif
  
   
!  burial_flux_of_dissolved_NH4_in_lake
   if (vDeltaS >= 0.0_rk) then
   tNBurNH4 = sNH4S *(vDeltaS /self%cDepthS) 
   else
   tNBurNH4 = self%cNH4Ground *(self%bPorS * vDeltaS) 
   endif
  
!  burial_flux_of_dissolved_NO3_in_lake
   if (vDeltaS >= 0.0_rk) then
   tNBurNO3 = sNO3S *(vDeltaS / self%cDepthS) 
   else
   tNBurNO3 = self%cNO3Ground *(self%bPorS * vDeltaS) 
   endif
  
!  burial_flux_of_Si_in_detritus_in_lake
   if (vDeltaS >= 0.0_rk) then
   tSiBurDet = rSiDDetS * tDBurDet 
   else
   tSiBurDet = 0.0_rk
   endif
! Humus burial fluxes
!  burial_flux_of_DW_in_humus_in_lake
   if (vDeltaS >= 0.0) then
     tDBurHum = tDBurOM - tDBurDet 
   else
     tDBurHum = tDBurOM 
   endif

!  burial_flux_of_P_in_humus_in_lake
   if (vDeltaS >= 0.0) then
      tPBurHum = rPDHumS * tDBurHum 
   else
      tPBurHum = 0.001_rk * tDBurHum   ! cPDSoilOM=0.001
   endif
   
!  burial_flux_of_N_in_humus_in_lake
   if (vDeltaS >= 0.0) then
      tNBurHum = rNDHumS * tDBurHum 
   else
      tNBurHum = 0.01_rk * tDBurHum   !cNDSoilOM =0.01
   endif
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
!  update inorganic and organic matters in water column
   _SET_BOTTOM_EXCHANGE_(self%id_SWDIM,uDErosIMW+tDResusIM-tDSetIM)
   _SET_BOTTOM_EXCHANGE_(self%id_SWDDet,tDResusDet-tDSetDet)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNDet,tNResusDet-tNSetDet)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPDet,tPResusDet-tPSetDet)
   _SET_BOTTOM_EXCHANGE_(self%id_SWSiDet,tSiResusDet-tSiSetDet)
!  update dissoved nutrients in water column
   _SET_BOTTOM_EXCHANGE_(self%id_SWNH4,tNResusNH4)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNO3,tNResusNO3)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPO4,tPResusPO4)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPAIM,tPResusAIM-tPSetAIM)
!  update phytoplankton in water column
   _SET_BOTTOM_EXCHANGE_(self%id_SWDDiat,tDResusDiat-tDSetDiat)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNDiat,tNResusDiat-tNSetDiat)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPDiat,tPResusDiat-tPSetDiat)
!   _SET_BOTTOM_EXCHANGE_(self%id_SWSiDiat,tSiResusDiat-tSiSetDiat)
   _SET_BOTTOM_EXCHANGE_(self%id_SWDGren,tDResusGren-tDSetGren)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNGren,tNResusGren-tNSetGren)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPGren,tPResusGren-tPSetGren)
   _SET_BOTTOM_EXCHANGE_(self%id_SWDBlue,tDResusBlue-tDSetBlue)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNBlue,tNResusBlue-tNSetBlue)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPBlue,tPResusBlue-tPSetBlue)
!  update abiotic variables in sediment
   _SET_ODE_BEN_(self%id_WSDIM,uDErosIMS+tDSetIM-tDResusIM-tDBurIM)
   _SET_ODE_BEN_(self%id_WSDDet,tDSetDet-tDResusDet-tDBurDet)
   _SET_ODE_BEN_(self%id_WSPDet,tPSetDet-tPResusDet-tPBurDet)
   _SET_ODE_BEN_(self%id_WSNDet,tNSetDet-tNResusDet-tNBurDet)
   _SET_ODE_BEN_(self%id_WSSiDet,tSiSetDet-tSiResusDet-tSiBurDet)
   _SET_ODE_BEN_(self%id_WSPO4,-tPResusPO4-tPBurPO4)
   _SET_ODE_BEN_(self%id_WSPAIM,tPSetAIM-tPResusAIM-tPBurAIM)
   _SET_ODE_BEN_(self%id_WSNH4,-tNResusNH4-tNBurNH4)
   _SET_ODE_BEN_(self%id_WSNO3,-tNResusNO3-tNBurNO3)
   _SET_ODE_BEN_(self%id_WSDHum,uDErosOM-tDBurHum)
   _SET_ODE_BEN_(self%id_WSPHum,uPErosOM-tPBurHum)
   _SET_ODE_BEN_(self%id_WSNHum,uNErosOM-tNBurHum)
!  update settled phytoplankton
   _SET_ODE_BEN_(self%id_WSDDiat,tDSetDiat-tDResusDiat)
   _SET_ODE_BEN_(self%id_WSNDiat,tNSetDiat-tNResusDiat)
   _SET_ODE_BEN_(self%id_WSPDiat,tPSetDiat-tPResusDiat)
   _SET_ODE_BEN_(self%id_WSDGren,tDSetGren-tDResusGren)
   _SET_ODE_BEN_(self%id_WSNGren,tNSetGren-tNResusGren)
   _SET_ODE_BEN_(self%id_WSPGren,tPSetGren-tPResusGren)
   _SET_ODE_BEN_(self%id_WSDBlue,tDSetBlue-tDResusBlue)
   _SET_ODE_BEN_(self%id_WSNBlue,tNSetBlue-tNResusBlue)
   _SET_ODE_BEN_(self%id_WSPBlue,tPSetBlue-tPResusBlue)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBurIM,tDBurIM)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_shearstress,shear)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aFunDimSusp,aFunDimSusp)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDResusDead,tDResusDead)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aFunTauSet,aFunTauSet)
   
   

   _FABM_HORIZONTAL_LOOP_END_
! Spatial loop end
   end subroutine do_bottom
!EOC
!-----------------------------------------------------------------------
!BOP

! !IROUTINE: !feh temperal solution for loading and dilution of NH4, NO3, 
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_au_pclake_auxilary), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_
! !LOCAL VARIABLES:
!  carriers for state variable values
   real(rk)    :: sDIMW,sDDetW,sPDetW,sNDetW,sSiDetW,sSiO2W
   real(rk)    :: sNH4W,sNO3W,sPO4W,sPAIMW,sO2W
   real(rk)    :: sDDiatW,sNDiatW,sPDiatW
   real(rk)    :: sDGrenW,sNGrenW,sPGrenW
   real(rk)    :: sDBlueW,sNBlueW,sPBlueW
   real(rk)    :: sDZoo,sPZoo,sNZoo
!  environmental dependency carrier
   real(rk)   :: sDepthW,Day
!  loading variables
   real(rk)    :: uPLoadPO4,uNLoadNO3,uDLoadIM,uSiLoadSiO2
!  dilution variables
   real(rk)    :: ukDil,uQEv,uQDil,ukDilWat
   real(rk)    :: wDDilIM,wDDilDet,wPDilDet,wNDilDet,wSiDilDet
   real(rk)    :: wNDilNH4,wNDilNO3,wPDilPO4,wPDilAIM
!  transport flux variables
   real(rk)    :: wDTranIMW,wDTranDetW,wPTranDetW,wNTranDetW,wSiTranDetW
   real(rk)    :: wNTranNH4W,wNTranNO3W,wPTranPO4W,wPTranAIMW
   real(rk)    :: wDTranDiat,wNTranDiat,wPTranDiat
   real(rk)    :: wDTranGren,wNTranGren,wPTranGren
   real(rk)    :: wDTranBlue,wNTranBlue,wPTranBlue
   real(rk)    :: wDDilDiat,wNDilDiat,wPDilDiat
   real(rk)    :: wDDilGren,wNDilGren,wPDilGren
   real(rk)    :: wDDilBlue,wNDilBlue,wPDilBlue
   real(rk)    :: wO2TranW,wO2Inflow ,wO2Outfl
   real(rk)    :: wDTranZoo,wNTranZoo,wPTranZoo
   real(rk)    :: wSiDilSiO2,wSiTranSiO2
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _LOOP_BEGIN_
! !Retrieve state dependencie value
!  from abiotic water column
   _GET_(self%id_SWNH4,sNH4W)
   _GET_(self%id_SWNO3,sNO3W)
   _GET_(self%id_SWPO4,sPO4W)
   _GET_(self%id_SWO2,sO2W)
   _GET_(self%id_SWSiO2,sSiO2W)
   _GET_(self%id_SWPAIM,sPAIMW)
   _GET_(self%id_SWDIM,sDIMW)
   _GET_(self%id_SWDDet,sDDetW)
   _GET_(self%id_SWNDet,sNDetW)
   _GET_(self%id_SWPDet,sPDetW)
   _GET_(self%id_SWSiDet,sSiDetW)
!  from phytoplankton in water column
   _GET_(self%id_SWDDiat,sDDiatW)
   _GET_(self%id_SWDGren,sDGrenW)
   _GET_(self%id_SWDBlue,sDBlueW)
   _GET_(self%id_SWNDiat,sNDiatW)
   _GET_(self%id_SWNGren,sNGrenW)
   _GET_(self%id_SWNBlue,sNBlueW)
   _GET_(self%id_SWPDiat,sPDiatW)
   _GET_(self%id_SWPGren,sPGrenW)
   _GET_(self%id_SWPBlue,sPBlueW)
!  from zooplankton
   _GET_(self%id_DTranZoo,sDZoo)
   _GET_(self%id_PTranZoo,sPZoo)
   _GET_(self%id_NTranZoo,sNZoo)
!  retrieve environmental dependencies
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
   _GET_GLOBAL_(self%id_Day,Day)
!  P_load_PO4
   uPLoadPO4=self%cLoadPO4
!  N_load_NO3
   uNLoadNO3=self%cLoadNO3
   uQEv = 0.0_rk
!  dilution_rate_of_substances
   uQDil=self%uQIn-uQEv 
!  currently ignore algal loadings.
!  dilution_rate_of_substances
   ukDil = uQDil / mmPerm/sDepthW
!  loading_of_DW_of_inorg_matter
   uDLoadIM = 5.0_rk * self%uQIn / mmPerm  ! cDIMIn=5
!  dilution_of_DW_IM
   wDDilIM = ukDil * sDIMW
!  dilllution_of_detritus
   wDDilDet = ukDil * sDDetW
!  dilution_of_SRP
   wPDilPO4 = ukDil * sPO4W
!  dilution_of_detritus
   wPDilDet = ukDil*sPDetW
!  dilution_of_IM-ads._P
   wPDilAIM = ukDil * sPAIMW
!  dilution_of_ammonium
   wNDilNH4 = ukDil * sNH4W
!  dilution_of_nitrate
   wNDilNO3 = ukDil * sNO3W
!  dilution_of_detritus
   wNDilDet = ukDil * sNDetW
!  dilution_of_det_Si
   wSiDilDet= ukDil * sSiDetW
!  dilution_of_Diat
   wDDilDiat = ukDil * sDDiatW
!  dilution_of_Diat
   wPDilDiat = ukDil * sPDiatW
!  dilution_of_Diat
   wNDilDiat = ukDil * sNDiatW
!  dilution_of_Gren
   wDDilGren = ukDil * sDGrenW
!  dilution_of_Gren
   wPDilGren = ukDil * sPGrenW
!  dilution_of_Gren
   wNDilGren = ukDil * sNGrenW
!  dilution_of_Blue
   wDDilBlue = ukDil * sDBlueW
!  dilution_of_Blue
   wPDilBlue = ukDil * sPBlueW
!  dilution_of_Blue
   wNDilBlue = ukDil * sNBlueW
!  transport_flux_DW_in_IM
   wDTranIMW = uDLoadIM/sDepthW- wDDilIM
!  transport_flux_DW_in_detritus
   wDTranDetW =  - wDDilDet
!  transport_flux_of_P_in_PO4
   wPTranPO4W = uPLoadPO4/sDepthW  - wPDilPO4
!  transport_flux_of_P_in_AIM
   wPTranAIMW =   - wPDilAIM
!  transport_flux_of_P_in_detritus
   wPTranDetW =  - wPDilDet 
!  transport_flux_of_N_in_NH4
   wNTranNH4W =  - wNDilNH4 
!  transport_flux_of_N_in_NO3
   wNTranNO3W = uNLoadNO3/sDepthW - wNDilNO3
!  transport_flux_of_N_in_detritus
   wNTranDetW = - wNDilDet
!  transport_flux_of_Si_in_detritus
   wSiTranDetW = - wSiDilDet
!  transport_flux_of_D_in_Diat
   wDTranDiat=-wDDilDiat
!  transport_flux_of_N_in_Diat
   wNTranDiat=-wNDilDiat
!  transport_flux_of_P_in_Diat
   wPTranDiat=-wPDilDiat
!  transport_flux_of_D_in_Gren
   wDTranGren=-wDDilGren
!  transport_flux_of_N_in_Gren
   wNTranGren=-wNDilGren
!  transport_flux_of_P_in_Gren
   wPTranGren=-wPDilGren
!  transport_flux_of_D_in_Blue
   wDTranBlue=-wDDilBlue
!  transport_flux_of_N_in_Blue
   wNTranBlue=-wNDilBlue
!  transport_flux_of_P_in_Blue
   wPTranBlue=-wPDilBlue
!  dilution_rate_of_water
   ukDilWat = self%uQIn / mmPerm / sDepthW
!  oxygen_inflow
   wO2Inflow = ukDilWat * 5.0_rk
!  oxygen_outflow
   wO2Outfl = ukDil * sO2W
!  transport_flux_O2
   wO2TranW = wO2Inflow - wO2Outfl
!  net_migration_flux_of_D_in_Zoo
   wDTranZoo =( ukDilWat * 0.1_rk - ukDil*sDZoo)  ! cDZooIn= 0.1
!  net_migration_flux_of_P_in_ZOO
   wPTranZoo =(ukDilWat *0.01_rk*0.1_rk  - ukDil*sPZoo)  ! cPDZooRef=0.01
!  net_migration_flux_of_N_in_Zoo
   wNTranZoo =(ukDilWat * 0.1_rk*0.07_rk - ukDil * sNZoo) !cNDZooRef=0.07
!  total_transport_flux_of_Si_in_SiO2
   uSiLoadSiO2 = 3.0_rk * self%uQIn / mmPerm   !cSiO2In=3.0
!  Dilution_of_Si_in_SiO2
   wSiDilSiO2 = ukDil * sSiO2W
!  transport_flux_of_Si_in_SIO2
   wSiTranSiO2 = uSiLoadSiO2 / sDepthW - wSiDilSiO2
!  update transported state variables
   _SET_ODE_(self%id_SWDIM,wDTranIMW)
   _SET_ODE_(self%id_SWDDet,wDTranDetW)
   _SET_ODE_(self%id_SWNDet,wNTranDetW)
   _SET_ODE_(self%id_SWPDet,wPTranDetW)
   _SET_ODE_(self%id_SWSiDet,wSiTranDetW)
   _SET_ODE_(self%id_SWNH4,wNTranNH4W)
   _SET_ODE_(self%id_SWNO3,wNTranNO3W)
   _SET_ODE_(self%id_SWPO4,wPTranPO4W)
   _SET_ODE_(self%id_SWPAIM,wPTranAIMW)
   _SET_ODE_(self%id_SWO2,wO2TranW)
   _SET_ODE_(self%id_SWSiO2,wSiTranSiO2)
!  update phytoplankton in water column
   _SET_ODE_(self%id_SWDDiat,wDTranDiat)
   _SET_ODE_(self%id_SWNDiat,wNTranDiat)
   _SET_ODE_(self%id_SWPDiat,wPTranDiat)
   _SET_ODE_(self%id_SWDGren,wDTranGren)
   _SET_ODE_(self%id_SWNGren,wNTranGren)
   _SET_ODE_(self%id_SWPGren,wPTranGren)
   _SET_ODE_(self%id_SWDBlue,wDTranBlue)
   _SET_ODE_(self%id_SWNBlue,wNTranBlue)
   _SET_ODE_(self%id_SWPBlue,wPTranBlue)
!  update zooplankton group
   _SET_ODE_(self%id_DTranZoo,wDTranZoo)
   _SET_ODE_(self%id_PTranZoo,wPTranZoo)
   _SET_ODE_(self%id_NTranZoo,wNTranZoo)

   

   _LOOP_END_
!-----------------------------------------------------------------------
!  Spatial loop end
!-----------------------------------------------------------------------
   end subroutine do

!EOC


   end module au_pclake_auxilary

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

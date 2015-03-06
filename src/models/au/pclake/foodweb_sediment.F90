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



!EOP                             
!-----------------------------------------------------------------------
!BOC                             
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cDBentIn,     'cDBentIn',     'gDW m-2',  'external zoobenthos density',                                default=0.01_rk)
   call self%get_parameter(self%kMigrBent,    'kMigrBent',    'd-1',      'zoobenthos migration rate',                                  default=0.001_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrBent,   'cDCarrBent',   'gDW m-2',  'carrying capacity of zoobenthos',                            default=10.0_rk)
   call self%get_parameter(self%kDAssBent,    'kDAssBent',    'd-1',      'maximum assimilation rate',                                  default=0.1_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDFoodBent,   'hDFoodBent',   'g m-2',    'half-saturating food for zoobenthos',                        default=200.0_rk)
   call self%get_parameter(self%fDAssBent,    'fDAssBent',    '[-]',      'C ass. efficiency of zoobenthos',                            default=0.3_rk)
   call self%get_parameter(self%fDissEgesBent,'fDissEgesBent','[-]',      'soluble nutrient fraction of by zoobenthos egested food',    default=0.25_rk)
   call self%get_parameter(self%kDRespBent,   'kDRespBent',   'd-1',      'maint. respiration constant of zoobenthos',                  default=0.005_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortBent,    'kMortBent',    'd-1',      'mortality constant of zoobenthos',                           default=0.005_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortBent,'fDissMortBent','[-]',      'soluble P fraction of died zoobenthos P',                    default=0.1_rk)
   call self%get_parameter(self%cTmOptBent,   'cTmOptBent',   '°C',       'optimum temp. of zoobenthos',                                default=25.0_rk)
   call self%get_parameter(self%cSigTmBent,   'cSigTmBent',   '°C',       'temperature constant of zoobenthos(sigma in Gaussian curve)',default=16.0_rk)
   call self%get_parameter(self%cPDBentRef,   'cPDBentRef',   'mgP/mgDW', 'reference P/C ratio of zoobenthos',                          default=0.01_rk)
   call self%get_parameter(self%cNDBentRef,   'cNDBentRef',   'mgN/mgDW', 'reference N/C ratio of zoobenthos',                          default=0.07_rk)
   call self%get_parameter(self%cSiDDiat,     'cSiDDiat',     'mgSi/mgDW','Si/DW ratio of daitoms',                                     default=0.15_rk)
   call self%get_parameter(self%fDAssFiAd,    'fDAssFiAd',    '[-]',      'C assimilation efficiency of adult fish',                    default=0.4_rk)
   call self%get_parameter(self%cPDFishRef,   'cPDFishRef',   'mgP/mgDW', 'reference P/C ratio of Fish',                                default=0.022_rk)
   call self%get_parameter(self%cNDFishRef,   'cNDFishRef',   'mgN/mgDW', 'reference N/C ratio of Fish',                                default=0.1_rk)
   call self%get_parameter(self%fDissEgesFish,'fDissEgesFish','[-]',      'soluble nutrient fraction of by fish egested food',          default=0.25_rk)
   call self%get_parameter(self%cTmOptFish,   'cTmOptFish',   '°C',       'optimum temp. of fish',                                      default=25.0_rk)
   call self%get_parameter(self%cSigTmFish,   'cSigTmFish',   '°C',       'temperature constant of fish(sigma in Gaussian curve)',      default=10.0_rk)
   call self%get_parameter(self%cRelVegFish,  'cRelVegFish',  '[-]',      'decrease of fish feeding per vegetation cover(max. 0.01)', default=0.009_rk)
   call self%get_parameter(self%kDAssFiAd,    'kDAssFiAd',    'd-1',      'maximum assimilation rate of adult fish',                    default=0.06_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDBentFiAd,   'hDBentFiAd',   'g m-2',    'half-saturating zoobenthos biomass for adult fish predation',default=2.5_rk)
   call self%get_parameter(self%kDRespFiAd,   'kDRespFiAd',   'd-1',      'maintenance respiration constant of adult fish',             default=0.004_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortFiAd,    'kMortFiAd',    'd-1',      'specific mortality of adult fish',                           default=0.00027_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrFish,   'cDCarrFish',   'gDW m-2',  'carrying capacity of fish',                                  default=15.0_rk)
   call self%get_parameter(self%cNDDiatMin,   'cNDDiatMin',   'mgN/mgDW', 'minimum N/day ratio Diatoms',                                default=0.01_rk)
   call self%get_parameter(self%cPDDiatMin,   'cPDDiatMin',   'mgP/mgDW', 'minimum P/day ratio Diatoms',                                default=0.0005_rk)
   call self%get_parameter(self%cNDGrenMin,   'cNDGrenMin',   'mgN/mgDW', 'minimum N/day ratio greens',                                 default=0.02_rk)
   call self%get_parameter(self%cPDGrenMin,   'cPDGrenMin',   'mgP/mgDW', 'minimum P/day ratio greens',                                 default=0.0015_rk)
   call self%get_parameter(self%cNDBlueMin,   'cNDBlueMin',   'mgN/mgDW', 'minimum N/day ratio Bluegreens',                             default=0.03_rk)
   call self%get_parameter(self%cPDBlueMin,   'cPDBlueMin',   'mgP/mgDW', 'minimum P/day ratio Bluegreens',                             default=0.0025_rk)
   call self%get_parameter(self%cNDBlueMax,   'cNDBlueMax',   'mgN/mgDW', 'max. N/day ratio Bluegreens',                                default=0.15_rk)
   call self%get_parameter(self%cNDDiatMax,   'cNDDiatMax',   'mgN/mgDW', 'max. N/day ratio Diatoms',                                   default=0.005_rk)
   call self%get_parameter(self%cNDGrenMax,   'cNDGrenMax',   'mgN/mgDW', 'max. N/day ratio greens',                                    default=0.1_rk)
   call self%get_parameter(self%cPDBlueMax,   'cPDBlueMax',   'mgP/mgDW', 'max. P/day ratio blue-greens',                               default=0.025_rk)
   call self%get_parameter(self%cPDDiatMax,   'cPDDiatMax',   'mgP/mgDW', 'max. P/day ratio Diatoms',                                   default=0.05_rk)
   call self%get_parameter(self%cPDGrenMax,   'cPDGrenMax',   'mgP/mgDW', 'max. P/day ratio greens',                                    default=0.015_rk)

!  Register local state variable
   call self%register_state_variable(self%id_sDBent,'sDBent','g m-2','zoobenthos_DW',     &
                                    initial_value=1.0_rk,minimum=WebZero)
   call self%register_state_variable(self%id_sPBent,'sPBent','g m-2','zoobenthos_P',     &
                                    initial_value=0.1_rk,minimum=WebZero)
   call self%register_state_variable(self%id_sNBent,'sNBent','g m-2','zoobenthos_N',     &
                                    initial_value=0.01_rk,minimum=WebZero)
!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_tDWebDetS,'tDWebDetS','g m-2 s-1','tDWebDetS', output=output_none)
   call self%register_diagnostic_variable(self%id_tDEnvFiAd,'tDEnvFiAd','g m-2',    'tDEnvFiAd', output=output_none)
   call self%register_diagnostic_variable(self%id_aDSatFiAd,'aDSatFiAd','g m-2',    'aDSatFiAd', output=output_none)
   
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNBent)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPBent)
!  regirster state variables dependencies 
   call self%register_state_dependency(self%id_DfoodDiatS, 'diatom_as_food_DW',        'g m-2', 'diatom_as_food_DW')
   call self%register_state_dependency(self%id_DfoodGrenS, 'green_as_food_DW',         'g m-2', 'green_as_food_DW')
   call self%register_state_dependency(self%id_DfoodBlueS, 'blue_as_food_DW',          'g m-2', 'blue_as_food_DW')
   call self%register_state_dependency(self%id_NfoodDiatS, 'diatom_as_food_N',         'g m-2', 'diatom_as_food_N')
   call self%register_state_dependency(self%id_NfoodGrenS, 'green_as_food_N',          'g m-2', 'green_as_food_N')
   call self%register_state_dependency(self%id_NfoodBlueS, 'blue_as_food_N',           'g m-2', 'blue_as_food_N')
   call self%register_state_dependency(self%id_PfoodDiatS, 'diatom_as_food_P',         'g m-2', 'diatom_as_food_P')
   call self%register_state_dependency(self%id_PfoodGrenS, 'green_as_food_P',          'g m-2', 'green_as_food_P')
   call self%register_state_dependency(self%id_PfoodBlueS, 'blue_as_food_P',           'g m-2', 'blue_as_food_P')
   call self%register_state_dependency(self%id_DDetpoolS,  'detritus_DW_pool_sediment','g m-2', 'detritus_DW_pool_sediment')
   call self%register_state_dependency(self%id_PDetpoolS,  'detritus_P_pool_sediment', 'g m-2', 'detritus_P_pool_sediment')
   call self%register_state_dependency(self%id_NDetpoolS,  'detritus_N_pool_sediment', 'g m-2', 'detritus_N_pool_sediment')
   call self%register_state_dependency(self%id_SiDetpoolS, 'detritus_Si_pool_sediment','g m-2', 'detritus_Si_pool_sediment')
   call self%register_state_dependency(self%id_NH4poolS,   'NH4_pool_sediment',        'g m-2', 'NH4_pool_sediment')
   call self%register_state_dependency(self%id_NO3poolS,   'NO3_pool_sediment',        'g m-2', 'NO3_pool_sediment')
   call self%register_state_dependency(self%id_PO4poolS,   'PO4_pool_sediment',        'g m-2', 'PO4_pool_sediment')
   call self%register_state_dependency(self%id_DAdFish,    'adult_fish_biomass',       'g m-3', 'adult_fish_biomass')
   call self%register_state_dependency(self%id_NAdFish,    'adult_fish_nitrogen',      'g m-3', 'adult_fish_nitrogen')
   call self%register_state_dependency(self%id_PAdFish,    'adult_fish_phosphrus',     'g m-3', 'adult_fish_phosphrus')
   call self%register_state_dependency(self%id_NH4poolW,   'NH4_pool_water',           'g m-3', 'NH4_pool_water')
   call self%register_state_dependency(self%id_PO4poolW,   'PO4_pool_water',           'g m-3', 'PO4_pool_water')
   call self%register_state_dependency(self%id_DDetpoolW, 'DDet_pool_water',           'g m-3', 'DDet_pool_water')
   call self%register_state_dependency(self%id_NDetpoolW, 'NDet_pool_water',           'g m-3', 'NDet_pool_water')
   call self%register_state_dependency(self%id_PDetpoolW, 'PDet_pool_water',           'g m-3', 'PDet_pool_water')
   call self%register_state_dependency(self%id_DJvFish,   'young_fish_biomass',        'g m-3', 'young_fish_biomass')
!  register diagnostic dependencies
   call self%register_dependency(self%id_aCovVeg, 'vegetation_coverage','[-]','vegetation_coverage')
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,    standard_variables%temperature)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
   
   return

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

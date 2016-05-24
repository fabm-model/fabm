#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module au_pclake_fish
!
! !DESCRIPTION:
!  The au_pclake_fish module describes the state variables regarding zooplankton-
!  vorous fish(sDFiJv,sNFiJv,sPFiAd),benthivorous fish(sDFiAd,sNFiAd,sPFiAd)and
!  piscivorous fish(sDPisc).
!  Therefore, the state variables and their related its local processes are:
!  sDFiJv,sPFiJv,sNFiJv,(zooplanktivorous fish in dry-weight, nitrogen element and phosphor-
!  -us element respectively)
!  units: gDW/m**3, gP/m**3,gP/m**3
!  local processes: migration,reproduction(+,part of benthivorous fish became young fi-
!  -sh),assimilation(predation on zooplankton),respiration(only for sDFiJv),exr-
!  -etion(only for sPFiJv,sNFiJv),mortality,consumption by piscivorous fish and
!  aging(-,part of zooplanktivorous fish become benthivorous fish)
!  sDFiAd,sPFiAd,sNFiAd,(benthivorous fish in dry-weight, nitrogen element and phosphorus element respectively)
!  units:gDW/m**3, gP/m**3,gP/m**3
!  local processes: migration,reproduction(-,part of benthivorous fish became young fi-
!  -sh),respiration(only for sDFiAd),exretion(only for sPFiAd and sNFiAd),morta-
!  -lity, consumption by piscivorous fish ,harvest and aging(+,part of young fi-
!  -sh become benthivorous fish).(!! Notice the assimilation of benthivorous fish is in the
!   benthic module, where benthivorous fish predating the zoobenthos.)
!  sDPisc,gDW/m**3)(piscivorous fish will have fixed N/D and P/D ratio,so the
!  other two will be diagnostic variables)
!  local processes: migration,assimilation,respiration,mortality and harvest
!  This module also describes the processes which influence the state variables
!  registered in other modules, including:
!  Detritus morted white fish and piscivorous fish, and egested by all fish:
!  sDDetW<==sDFiJv&sDFiAd&sDPisc,sNDetW<==sNFiJv&sNFiAd&aNPisc,sPDetW<==sPFiJv&sPFiAd&aPPisc
!  nutrients excreted by fish and dissolved part of mortality and egestion:
!   sNH4W<==sNFiJv&sNFiAd&aNPisc,sPO4<==sPFiJv&sPFiAd&aPPisc
! !USES:
   use fabm_types
   use au_pclake_utility, ONLY:uFunTmBio

   implicit none

!  default: all is private.
      private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_au_pclake_fish
!     local state variable identifiers
!     id_sDFiJv,zooplanktivorous fish concentration in dry-weight, gDW/m**3
!     id_sPFiJv,zooplanktivorous fish concentration in nitrogen element, gN/m**3
!     id_sNFiJv,zooplanktivorous fish concentration in phosphorus element, gP/m**3
!     id_sDFiAd,benthivoros fish concentration in dry-weight, gDW/m**3
!     id_sPFiAd,benthivoros fish concentration in nitrogen element, gN/m**3
!     id_sNFiAd,benthivoros fish concentration in phosphorus element, gP/m**3
!     id_sDPisc,piscivorous fish concentration in dry-weight, gDW/m**3
      type (type_state_variable_id)            :: id_sDFiJv,id_sPFiJv,id_sNFiJv
      type (type_state_variable_id)            :: id_sDFiAd,id_sPFiAd,id_sNFiAd,id_sDPisc
!     Fish manipulation,changed biomass in fish
      type (type_state_variable_id)            :: id_ChangedFiAd,id_ChangedFiJv,id_ChangedPisc
!     diagnostic variables for local output
!     id_aNPisc, piscivorous fish concentration in nitrogen element, gN/m**3
!     id_aPPisc, piscivorous fish concentration in phosphorus element, gP/m**3
      type (type_diagnostic_variable_id)       :: id_aNPisc,id_aPPisc
!    diagnostic variables for modular fluxes
      type (type_diagnostic_variable_id)       :: id_wDFiJv,id_wPFiJv,id_wNFiJv
      type (type_diagnostic_variable_id)       :: id_wDFiAd ,id_wPFiAd,id_wNFiAd
      type (type_diagnostic_variable_id)       :: id_wDPisc ,id_wPFishPO4W,id_wNFishNH4W
      type (type_diagnostic_variable_id)       :: id_wDFishDetW,id_wNFishDetW,id_wPFishDetW
      type (type_diagnostic_variable_id)       :: id_wDFishZoo,id_wNFishZoo,id_wPFishZoo
!     state dependencies identifiers
      type (type_state_variable_id)            :: id_DDetpoolW,id_PDetpoolW,id_NDetpoolW
      type (type_state_variable_id)            :: id_NH4poolW,id_PO4poolW
      type (type_state_variable_id)            :: id_DFoodZoo,id_NFoodZoo,id_PFoodZoo
!     environmental dependencies
      type (type_dependency_id)                :: id_uTm ,id_dz
      type (type_horizontal_dependency_id)     :: id_sDepthW
      type (type_global_dependency_id)         :: id_Day
!     diagnostic dependencies
      type (type_horizontal_dependency_id)     :: id_aDSubVeg,id_tDEnvFiAd,id_aDSatFiAd
!     Fish manipulation, external fish manipulation rate
      type (type_horizontal_dependency_id)     :: id_ManFiAd,id_ManFiJv,id_ManPisc
!!    Model parameters
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
!     Fish manipulation parameters, switch for turned on/off fish manipulation
      logical    :: Manipulate_FiAd, Manipulate_FiJv, Manipulate_Pisc
   contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      !procedure :: do_bottom
      procedure :: get_light_extinction

   end type type_au_pclake_fish

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk)           :: Pi=3.14159265358979_rk
!   Lowest state variable value for fish
   real(rk),parameter :: FishZero=0.00001_rk
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
   class (type_au_pclake_fish), intent(inout), target :: self
   integer,                     intent(in)            :: configunit

!EOP
!-----------------------------------------------------------------------
!BOC

!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%kMigrFish,       'kMigrFish',      'd-1',       'fish migration rate',                                                            default=0.001_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDFiJvIn,        'cDFiJvIn',       'gDW m-2',   'external fish density',                                                          default=0.005_rk)
   call self%get_parameter(self%cDFiAdIn,        'cDFiAdIn',       'gDW m-2',   'external fish density',                                                          default=0.005_rk)
   call self%get_parameter(self%cDPiscIn,        'cDPiscIn',       'gDW m-2',   'external Pi_sc. density',                                                        default=0.001_rk)
   call self%get_parameter(self%kMigrPisc,       'kMigrPisc',      'd-1',       'Pi_sc. migration rate',                                                          default=0.001_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDBone,          'fDBone',         '[-]',       'fraction of fish C fixed in bones and scales',                                   default=0.35_rk)
   call self%get_parameter(self%fPBone,          'fPBone',         '[-]',       'fraction of fish P fixed in bones and scales',                                   default=0.5_rk)
   call self%get_parameter(self%cDCarrFish,      'cDCarrFish',     'gDW m-2',   'carrying capacity of fish',                                                      default=15.0_rk)
   call self%get_parameter(self%fDissEgesFish,   'fDissEgesFish',  '[-]',       'soluble nutrient fraction of by fish egested food',                              default=0.25_rk)
   call self%get_parameter(self%fDissMortFish,   'fDissMortFish',  '[-]',       'soluble nutrient fraction of died fish(excl. bones and scales',                  default=0.1_rk)
   call self%get_parameter(self%cTmOptFish,      'cTmOptFish',     '�C',        'optimum temp. of fish',                                                          default=25.0_rk)
   call self%get_parameter(self%cSigTmFish,      'cSigTmFish',     '�C',        'temperature constant of fish(sigma in Gaussian curve)',                          default=10.0_rk)
   call self%get_parameter(self%cDayReprFish,    'cDayReprFish',   '[-]',       'reproduction date of fish ',                                                     default=120.0_rk)
   call self%get_parameter(self%fReprFish,       'fReprFish',      '[-]',       'yearly reproduction fraction of benthivorous fish, daily rate',                  default=0.02_rk)
   call self%get_parameter(self%fAgeFish,        'fAgeFish',       '[-]',       'yearly ageing fraction of zooplanktivorous fish,.daily rate',                    default=0.5_rk)
   call self%get_parameter(self%kDAssFiJv,       'kDAssFiJv',      'd-1',       'maximum assimilation rate of zooplanktivorous fish',                             default=0.12_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDZooFiJv,       'hDZooFiJv',      'g m-2',     'half-saturating zooplankton biomass for zooplanktivorous fish predation',        default=1.25_rk)
   call self%get_parameter(self%fDAssFiJv,       'fDAssFiJv',      '[-]',       'C assimilation efficiency of zooplanktivorous fish',                             default=0.4_rk)
   call self%get_parameter(self%kDRespFiJv,      'kDRespFiJv',     'd-1',       'maintenance respiration constant of zooplanktivorous fish',                      default=0.01_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortFiJv,       'kMortFiJv',      'd-1',       'specific mortality of zooplanktivorous fish',                                    default=0.00137_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespFiAd,      'kDRespFiAd',     'd-1',       'maintenance respiration constant of benthivorous fish',                          default=0.004_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortFiAd,       'kMortFiAd',      'd-1',       'specific mortality of benthivorous fish',                                        default=0.00027_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrPiscMax,   'cDCarrPiscMax',  'gDW m-2',   'maximum carrying capacity of  Pi_sc',                                            default=1.2_rk)
   call self%get_parameter(self%cDCarrPiscMin,   'cDCarrPiscMin',  'gDW m-2',   'minimum carrying capacity of  Pi_sc',                                            default=0.1_rk)
   call self%get_parameter(self%cDCarrPiscBare,  'cDCarrPiscBare', 'gDW m-2',   'carrying capacity of  Pi sc for lake without marsh zone',                        default=0.1_rk)
   call self%get_parameter(self%cDPhraMinPisc,   'cDPhraMinPisc',  'gDW m-2',   'min. reed biomass for  Pi_sc',                                                   default=50.0_rk)
   call self%get_parameter(self%cCovVegMin,      'cCovVegMin',     '%',         'min. subm.veg. coverage for  Pi_sc',                                             default=40.0_rk)
   call self%get_parameter(self%cRelPhraPisc,    'cRelPhraPisc',   'gDW m-2',   'rel.  Pi_sc density per reed if subm.veg. absent',                              default=0.075_rk)
   call self%get_parameter(self%cRelVegPisc,     'cRelVegPisc',    'gDW m-2',   'extra rel.  Pi_sc density perreed if  aCovVeg  >  cCovVegMin',                  default=0.03_rk)
   call self%get_parameter(self%kDAssPisc,       'kDAssPisc',      'd-1',       'maximum assimilation rate',                                                      default=0.025_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDVegPisc,       'hDVegPisc',      'g m-2',     'half-sat. vegetation biomass for  Pi_sc growth',                                 default=5.0_rk)
   call self%get_parameter(self%hDFishPisc,      'hDFishPisc',     'g m-2',     'half-saturating DFish for  Pi_sc predation',                                     default=1.0_rk)
   call self%get_parameter(self%fDAssPisc,       'fDAssPisc',      '[-]',       'C ass. efficiency of  Pi_sc',                                                    default=0.4_rk)
   call self%get_parameter(self%fDissEgesPisc,   'fDissEgesPisc',  '[-]',       'soluble P fraction of by fish egested food',                                     default=0.25_rk)
   call self%get_parameter(self%kDRespPisc,      'kDRespPisc',     'd-1',       'maint. respiration constant of  Pi_sc',                                          default=0.005_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortPisc,       'kMortPisc',      'd-1',       'specific mortality of  Pi_sc',                                                   default=0.00027_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortPisc,   'fDissMortPisc',  '[-]',       'soluble nutrient fraction of died  Pi_sc(excl. bones and scales',                default=0.1_rk)
   call self%get_parameter(self%cTmOptPisc,      'cTmOptPisc',     '�C',        'optimum temp. of  Pi_sc',                                                        default=25.0_rk)
   call self%get_parameter(self%cSigTmPisc,      'cSigTmPisc',     '�C',        'temperature constant of  Pi_sc(sigma in Gaussian curve)',                        default=10.0_rk)
   call self%get_parameter(self%cPDFishRef,      'cPDFishRef',     'mgP/mgDW',  'reference P/C ratio of Fish',                                                    default=0.022_rk)
   call self%get_parameter(self%cNDFishRef,      'cNDFishRef',     'mgN/mgDW',  'reference N/C ratio of Fish',                                                    default=0.1_rk)
   call self%get_parameter(self%cPDPisc,         'cPDPisc',        'mgP/mgDW',  'reference P/C ratio of  Pi_sc',                                                  default=0.022_rk)
   call self%get_parameter(self%cNDPisc,         'cNDPisc',        'mgN/mgDW',  'reference N/C ratio of  Pi_sc ',                                                 default=0.1_rk)
!  Fish manipulation, register of switches
   call self%get_parameter(self%Manipulate_FiAd, 'Manipulate_FiAd', ' ',        'Turn on benthivorous fish manipulation',                                         default=.false.)
   call self%get_parameter(self%Manipulate_FiJv, 'Manipulate_FiJv', ' ',        'Turn on zooplanktivorous manipulation',                                          default=.false.)
   call self%get_parameter(self%Manipulate_Pisc,'Manipulate_Pisc',  ' ',        'Turn on piscivorous fish manipulation',                                          default=.false.)
!  Register local state variable
!  zooplanktivorous fish, transportation is turned off
   call self%register_state_variable(self%id_sDFiJv,'sDFiJv','g m-3','zooplanktivorous fish biomass',     &
                                    initial_value= 0.5_rk,minimum=NearZero,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sDFiJv,'disable_transport',.true.)
   call self%register_state_variable(self%id_sPFiJv,'sPFiJv','g m-3','zooplanktivorous fish phosphorus content',     &
                                    initial_value=0.011_rk,minimum=NearZero,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sPFiJv,'disable_transport',.true.)
   call self%register_state_variable(self%id_sNFiJv,'sNFiJv','g m-3','zooplanktivorous fish nitrogen content',     &
                                    initial_value=0.05_rk,minimum=NearZero,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_spFiJv,'disable_transport',.true.)
!  benthivoros fish, transportation turned off
   call self%register_state_variable(self%id_sDFiAd,'sDFiAd','g m-3','benthivorous fish biomass',     &
                                    initial_value=2.0_rk,minimum=NearZero,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sDFiAd,'disable_transport',.true.)
   call self%register_state_variable(self%id_sPFiAd,'sPFiAd','g m-3','benthivorous fish phosphorus content',     &
                                    initial_value=0.044_rk,minimum=NearZero,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sPFiAd,'disable_transport',.true.)
   call self%register_state_variable(self%id_sNFiAd,'sNFiAd','g m-3','benthivorous fish nitrogen content',     &
                                    initial_value=0.2_rk,minimum=NearZero,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sNFiAd,'disable_transport',.true.)
!  piscivorous fish
   call self%register_state_variable(self%id_sDPisc,'sDPisc','g m-3','piscivorous fish biomass', &
                                    initial_value=0.01_rk,minimum=NearZero,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sDPisc,'disable_transport',.true.)
!  Fish manipulation, if manipulation turned on, then registered state variable of fish biomass change
!  as well as register external fish manipulation rate
!  If benthivorous fish manipulation tured on
   If(self%Manipulate_FiAd) then
     call self%register_state_variable(self%id_ChangedFiAd, 'ChangedFiAd',          '',    'changed benthivorous fish biomass', 0.0_rk)
     call self%register_dependency(self%id_ManFiAd,         'manipulate_rate_FiAd', 's-1', 'benthivorous fish manipulate rate')
   endif
!  If zooplanktivorous turned on
   If(self%Manipulate_FiJv) then
     call self%register_state_variable(self%id_ChangedFiJv,  'ChangedFiJv',          '',    'changed zooplanktivorous biomass', 0.0_rk)
     call self%register_dependency(self%id_ManFiJv,          'manipulate_rate_FiJv', 's-1', 'zooplanktivorous manipulate rate')
   endif
!  If piscivorous fish turned on
   If(self%Manipulate_Pisc) then
     call self%register_state_variable(self%id_ChangedPisc,  'ChangedPisc',          '',    'changed piscivorous fish biomass', 0.0_rk)
     call self%register_dependency(self%id_ManPisc,          'manipulate_rate_Pisc', 's-1', 'piscivorous fish manipulate rate')
   endif
!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_aNPisc,    'aNPisc',      'g m-3',     'Piscivorous fish nitrogen content',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aPPisc,    'aPPisc',      'g m-3',     'Piscivorous fish phosphorus content', output=output_instantaneous)
!  register diagnostic variables for modular fluxes
   call self%register_diagnostic_variable(self%id_wDFiJv,     'wDFiJv',     'g m-3 s-1', 'fish_DFiJv_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPFiJv,     'wPFiJv',     'g m-3 s-1', 'fish_PFiJv_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNFiJv,     'wNFiJv',     'g m-3 s-1', 'fish_NFiJv_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDFiAd,     'wDFiAd',     'g m-3 s-1', 'fish_DFiAd_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPFiAd,     'wPFiAd',     'g m-3 s-1', 'fish_PFiAd_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNFiAd,     'wNFiAd',     'g m-3 s-1', 'fish_NFiAd_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDPisc,     'wDPisc',     'g m-3 s-1', 'fish_DPisc_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNFishNH4W, 'wNFishNH4W', 'g m-3 s-1', 'fish_NH4W_change',                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPFishPO4W, 'wPFishPO4W', 'g m-3 s-1', 'fish_PO4W_change',                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDFishDetW, 'wDFishDetW', 'g m-3 s-1', 'fish_DDetW_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNFishDetW, 'wNFishDetW', 'g m-3 s-1', 'fish_NDetW_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPFishDetW, 'wPFishDetW', 'g m-3 s-1', 'fish_PDetW_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDFishZoo,  'wDFishZoo',  'g m-3 s-1', 'fish_DZoo_change',                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNFishZoo,  'wNFishZoo',  'g m-3 s-1', 'fish_NZoo_change',                     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPFishZoo,  'wPFishZoo',  'g m-3 s-1', 'fish_PZoo_change',                     output=output_instantaneous)
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNFiJv)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNFiAd)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_aNPisc)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPFiJv)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPFiAd)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_aPPisc)
!  register state variables dependencies
   call self%register_state_dependency(self%id_DDetpoolW, 'detritus_DW_pool_water','g m-3', 'detritus_DW_pool_water')
   call self%register_state_dependency(self%id_NDetpoolW, 'detritus_N_pool_water', 'g m-3', 'detritus_N_pool_water')
   call self%register_state_dependency(self%id_PDetpoolW, 'detritus_P_pool_water', 'g m-3', 'detritus_P_pool_water')
   call self%register_state_dependency(self%id_NH4poolW,  'NH4_pool_water',        'g m-3', 'NH4_pool_water')
   call self%register_state_dependency(self%id_PO4poolW,  'PO4_pool_water',        'g m-3', 'PO4_pool_water')
   call self%register_state_dependency(self%id_DFoodZoo,  'Zooplankton_D_Food',    'g m-3', 'Zooplankton_D_Food')
   call self%register_state_dependency(self%id_NFoodZoo,  'Zooplankton_N_Food',    'g m-3', 'Zooplankton_N_Food')
   call self%register_state_dependency(self%id_PFoodZoo,  'Zooplankton_P_Food',    'g m-3', 'Zooplankton_P_Food')
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,    standard_variables%temperature)
   call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_dz,     standard_variables%cell_thickness)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
!  register diagnostic dependencies
   call self%register_dependency(self%id_tDEnvFiAd, 'env_correction_adfish',     '[-]',  'env_correction_adfish')
   call self%register_dependency(self%id_aDSubVeg,  'submerged_vegetation',      'g m-2','submerged_vegetation')
   call self%register_dependency(self%id_aDSatFiAd, 'food_limit_function_adfish','[-]',  'food_limit_function_adfish')


   return


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
   class (type_au_pclake_fish), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_
!  LOCAL VARIABLES:
!  Carriers for environment dependencies
   real(rk)     :: uTm,Day,dz,sDepthW
!  carriers for local state variables
   real(rk)      :: sDFiJv,sPFiJv,sNFiJv
   real(rk)      :: sDFiAd,sPFiAd,sNFiAd,sDPisc
!  carriers for exteral link state variables
  real(rk)      :: sDZoo,sNZoo,sPZoo
!  carriers for external link diagnostic variables
   real(rk)      :: aDSubVeg,tDEnvFiAd,aDSatFiAd
!  nutrient ratios variables
   real(rk)      :: rPDFiJv,rNDFiJv,rPDFiAd,rNDFiAd
   real(rk)      :: rPDZoo,rNDZoo
!  status auxiliaries
   real(rk)      :: aDFish,aPFish,aNFish
!  variables for temperature function
   real(rk)      :: uFunTmFish,uFunTmPisc
!  variables for fish(include JV and AD)
!  PCLake_Osis, /m^2
   real(rk)      :: tDReprFish,tDAgeFish
   real(rk)      :: tPReprFish,tPAgeFish
   real(rk)      :: tNReprFish,tNAgeFish
!  variables for zooplanktivorous fish flux_DW
   real(rk)     :: wDFiJv
!  PCLake_Osis, /m^2
   real(rk)     :: tDFiJv,tDMigrFiJv
   real(rk)     :: tDAssFiJv,tDRespFiJv,tDMortFiJv,tDConsFiJvPisc
   real(rk)     :: aDSatFiJv,tDEnvFiJv,ukDIncrFiJv,tDConsFiJv
!  variables for zooplanktivorous fish flux_P
   real(rk)     :: wPFiJv
!  PCLake_Osis, /m^2
   real(rk)     :: tPFiJv,tPMigrFiJv,tPAssFiJv
   real(rk)     :: tPExcrFiJv,tPMortFiJv,tPConsFiJvPisc
   real(rk)     :: afPAssFiJv,tPConsFiJv,afNAssFiJv,tNConsFiJv
!  variables for zooplanktivorous fish flux_N
   real(rk)     :: wNFiJv
!  ,wNMigrFiJv,wNAssFiJv
!   real(rk)     :: wNExcrFiJv,wNMortFiJv,wNConsFiJvPisc
!  PCLake_Osis, /m^2
   real(rk)     :: tNFiJv,tNMigrFiJv,tNAssFiJv
   real(rk)     :: tNExcrFiJv,tNMortFiJv,tNConsFiJvPisc
!  variables for benthivorous fish flux_DW
   real(rk)     :: wDFiAd
!  PCLake_Osis, /m^2
   real(rk)     :: tDFiAd,tDMigrFiAd,tDRespFiAd,tDMortFiAd
   real(rk)     :: tDConsFiAdPisc
   real(rk)     :: wPFiAd
!  PCLake_Osis, /m^2
   real(rk)     :: tPFiAd,tPMigrFiAd,tPExcrFiAd,tPMortFiAd
   real(rk)     :: tPConsFiAdPisc
!  assimilation
!  variables for benthivorous fish flux_N
   real(rk)     :: wNFiAd
!  PCLake_Osis, /m^2
   real(rk)     :: tNFiAd,tNMigrFiAd,tNExcrFiAd,tNMortFiAd
   real(rk)     :: tNConsFiAdPisc
!  variables for piscivorous fish ,DW process
   real(rk)      :: wDPisc
!  PCLake_Osis, /m^2
   real(rk)     :: tDConsPisc,tDAssPisc,aDSatPisc,aFunVegPisc
   real(rk)     :: tDEnvPisc,akDIncrPisc,aDCarrPisc
   real(rk)     :: tDMigrPisc,tDRespPisc,tDMortPisc
   real(rk)     :: tDPisc
!  variables for piscivorous fish P process
!  PCLake_Osis, /m^2
   real(rk)     :: aPPisc,tPConsPisc,rPDFoodPisc,afPAssPisc,tPAssPisc
   real(rk)     :: tPEgesPisc,tPExcrPisc,tPMortPisc,tPMigrPisc
!  variables for piscivorous fish N process
!  PCLake_Osis, /m^2
   real(rk)     :: aNPisc,tNConsPisc,rNDFoodPisc,afNAssPisc,tNAssPisc
   real(rk)     :: tNEgesPisc,tNExcrPisc,tNMortPisc,tNMigrPisc
!  variables for exchange of NH4
!  PCLake_Osis, /m^2
   real(rk)     :: wNFishNH4W,wNEgesZooNH4,wNEgesZoo,wNMortZooNH4,tNEgesFiJvNH4
   real(rk)     :: tNEgesFiJv,tNMortFishNH4,tNMortFishBot
   real(rk)     :: tNMortFish,tNEgesPiscNH4,tNMortPiscNH4,tNMortPiscBot
!  variables for exchange of PO4
!  PCLake_Osis, /m^2
   real(rk)     :: wPFishPO4W,tPEgesFiJvPO4,tPFishPO4W
   real(rk)     :: tPEgesFiJv,tPMortFish,tPMortFishBot
   real(rk)     :: tPMortFishPO4,tPEgesPiscPO4,tPMortPiscPO4,tPMortPiscBot
!  variables for exchange of Detritus DW
!  PCLake_Osis, /m^2
   real(rk)     :: wDFishDetW,tDEgesFiJv, tDMortFishDet
   real(rk)     :: tDMortFish,tDMortFishBot,tDEgesPisc,tDMortPiscDet,tDMortPiscBot
!  variables for exchange of Detritus N
!  PCLake_Osis, /m^2
   real(rk)     :: wNFishDetW,tNEgesFiJvDet,tNMortFishDet
   real(rk)     :: tNEgesPiscDet,tNMortPiscDet
!  variables for exchange of detritus P
!  PCLake_Osis, /m^2
   real(rk)     :: wPFishDetW,tPEgesFiJvDet,tPMortFishDet
   real(rk)     :: tPEgesPiscDet,tPMortPiscDet
!  variables for exchange of detritus Si
   real(rk)     :: wSiConsDiatZoo
!  variables for exchange of diatoms
   real(rk)     :: wDFishDiatW,wNFishDiatW,wPFishDiatW
!  variables for exchange of green algae
   real(rk)     :: wDFishGrenW,wNFishGrenW,wPFishGrenW
!  variables for exchange of green algae
   real(rk)     :: wDFishBlueW,wNFishBlueW,wPFishBlueW
!  benthivorous fish assimilation
   real(rk)     :: ukDIncrFiAd
!  Fish manipulation, manipulate rate variable
   real(rk)     :: rManFiAd,rManFiJv,rManPisc
   real(rk)     :: tDManFiAd,tDManFiJv,tDManPisc
   real(rk)     :: tNManFiAd,tNManFiJv
   real(rk)     :: tPManFiAd,tPManFiJv
   real(rk)     :: ChangedFiAd,ChangedFiJv,ChangedPisc
   integer, save :: n=0

!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_(self%id_sDFiJv,sDFiJv)
   _GET_(self%id_sPFiJv,sPFiJv)
   _GET_(self%id_sNFiJv,sNFiJv)
   _GET_(self%id_sDFiAd,sDFiAd)
   _GET_(self%id_sPFiAd,sPFiAd)
   _GET_(self%id_sNFiAd,sNFiAd)
   _GET_(self%id_sDPisc,sDPisc)
!-----------------------------------------------------------------------
!  Retrieve dependencies  value
!-----------------------------------------------------------------------
   _GET_(self%id_DFoodZoo,sDZoo)
   _GET_(self%id_NFoodZoo,sNZoo)
   _GET_(self%id_PFoodZoo,sPZoo)

!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_GLOBAL_(self%id_Day,Day)
   _GET_(self%id_dz,dz)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
! !retrieve diagnostic dependency
   _GET_HORIZONTAL_(self%id_aDSubVeg,aDSubVeg)
   _GET_HORIZONTAL_(self%id_tDEnvFiAd,tDEnvFiAd)
   _GET_HORIZONTAL_(self%id_aDSatFiAd,aDSatFiAd)
!  Fish manipulation, if turned on, get the dependencies
!  and state variables
!  If benthivorous fish manipulation tured on
  If(self%Manipulate_FiAd) then
    _GET_HORIZONTAL_(self%id_ManFiAd,rManFiAd)
    _GET_(self%id_ChangedFiAd,ChangedFiAd)
  endif
!  If zooplanktivorous manipulation tured on
  If(self%Manipulate_FiJv) then
     _GET_HORIZONTAL_(self%id_ManFiJv,rManFiJv)
     _GET_(self%id_ChangedFiJv,ChangedFiJv)
  endif
!  If piscivorous fish manipulation tured on
  If(self%Manipulate_Pisc) then
     _GET_HORIZONTAL_(self%id_ManPisc,rManPisc)
     _GET_(self%id_ChangedPisc,ChangedPisc)
  endif

!  convert fish concentration to areal units
   sDFiJv=sDFiJv*sDepthW
   sPFiJv=sPFiJv*sDepthW
   sNFiJv=sNFiJv*sDepthW
   sDFiAd=sDFiAd*sDepthW
   sPFiAd=sPFiAd*sDepthW
   sNFiAd=sNFiAd*sDepthW
   sDPisc=sDPisc*sDepthW
!-------------------------------------------------------------------------
!  The orders for the processes. We try to orgnize the order
!  from zooplanktivorous fish, to benthivorous fish, and at last piscivorous
!  fish. And each group finish all their process before going to the next
!  groups. But there a special sections on zooplankton comsumption by
!  zooplanktivorous fish, due to the variables dependent on each other.It
!  starts with assimilation of DW to provide wDAssFiJv for its predation of zoo-
!  -plankton(DW,N,P).The later process(predation) provide variable wDConsFiJv
!  for fish assimilation of N,P
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!  Current local nutrients ratios in zooplankton(check the current state)
!-------------------------------------------------------------------------
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
!  status auxiliaries---auxiliaries for describing the current status,
!  usually derivatives of state variables
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
!  temp._function_of_fish
   uFunTmFish = uFunTmBio(uTm,self%cSigTmFish,self%cTmOptFish)
!  temp._function_of_Pisc
   uFunTmPisc = uFunTmBio(uTm,self%cSigTmPisc,self%cTmOptPisc)
!-----------------------------------------------------------------------
!  zooplanktivorous fish assimilation_DW
!-----------------------------------------------------------------------
!  food_limitation_function_of_young_fish
   aDSatFiJv = (sDZoo * sDepthW) *(sDZoo * sDepthW) /(self%hDZooFiJv * &
   &self%hDZooFiJv + (sDZoo * sDepthW) *(sDZoo * sDepthW))
!  intrinsic_net_increase_rate_of_fish
   ukDIncrFiJv = (self%kDAssFiJv - self%kDRespFiJv) * uFunTmFish - self%kMortFiJv
!  environmental_correction_of_fish
   tDEnvFiJv = max(0.0_rk,ukDIncrFiJv /(self%cDCarrFish - sDFiAd) * sDFiJv*sDFiJv)
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
!  zooplanktivorous fish assimilation_P
!-----------------------------------------------------------------------
!  P_assim._efficiency_of_FiJv
   afPAssFiJv = min(1.0_rk,self%cPDFishRef / rPDZoo * self%fDAssFiJv)
!  P_assimilation_of_FiJv
!  PCLake_osis:sPFiJv,in g/m^2
   tPAssFiJv = afPAssFiJv * tPConsFiJv
!-----------------------------------------------------------------------
!  zooplanktivorous fish assimilation_N
!-----------------------------------------------------------------------
!  N_assim._efficiency_of_FiJv
   afNAssFiJv = min(1.0_rk,self%cNDFishRef / rNDZoo * self%fDAssFiJv)
!  N_assimilation_of_FiJv
!  PCLake_osis:sNFiJv,in g/m^2
   tNAssFiJv = afNAssFiJv * tNConsFiJv
!-----------------------------------------------------------------------
!  zooplanktivorous fish migration
!-----------------------------------------------------------------------
!  migration_flux of zooplanktivorous fish, DW
!  PCLake_osis:sDFiJv,in g/m^2
   tDMigrFiJv = self%kMigrFish *(self%cDFiJvIn - sDFiJv)
!  net_migration_flux of zooplanktivorous fish,P
!  PCLake_osis:sPFiJv,in g/m^2
   tPMigrFiJv = self%kMigrFish *(self%cPDFishRef * self%cDFiJvIn - sPFiJv)
!  net_migration_flux of zooplanktivorous fish,N
!  PCLake_osis:sNFiJv,in g/m^2
   tNMigrFiJv = self%kMigrFish *(self%cNDFishRef * self%cDFiJvIn - sNFiJv)
!-----------------------------------------------------------------------
!  benthivorous fish migration
!-----------------------------------------------------------------------
!  migration_flux of benthivorous fish,DW
!  PCLake_osis:sDFiAd,in g/m^2
   tDMigrFiAd = self%kMigrFish *(self%cDFiAdIn - sDFiAd)
!  net_migration_flux of benthivorous fish, P
!  PCLake_osis:sPFiAd,in g/m^2
   tPMigrFiAd = self%kMigrFish *(self%cPDFishRef * self%cDFiAdIn - sPFiAd)
!  net_migration_flux of benthivorous fish, N
!  PCLake_osis:sPFiAd,in g/m^2
   tNMigrFiAd = self%kMigrFish *(self%cNDFishRef * self%cDFiAdIn - sNFiAd)
!-----------------------------------------------------------------------
!  fish reproduction
!-----------------------------------------------------------------------
!  Reproduction_flux_DW
   if (Day >= self%cDayReprFish .and. Day < self%cDayReprFish + 1.0_rk) then
!  PCLake_osis:sDFiAd,in g/m^2
   tDReprFish = self%fReprFish * sDFiAd/secs_pr_day
   else
   tDReprFish =0.0_rk
   endif
!  Reproduction_flux_P
   tPReprFish = rPDFiAd * tDReprFish
!  Reproduction_flux_N
   tNReprFish = rNDFiAd * tDReprFish
!-----------------------------------------------------------------------
!  fish aging
!-----------------------------------------------------------------------
!  Ageing_DW
   if (Day >=  365.0_rk .AND. Day <= 365.0_rk) then
!  PCLake_osis:sDFiAd,in g/m^2
   tDAgeFish = self%fAgeFish * sDFiJv/secs_pr_day
   else
    tDAgeFish = 0.0_rk
   endif
!  Ageing_P
   tPAgeFish = rPDFiJv * tDAgeFish
!  Ageing_N
   tNAgeFish = rNDFiJv * tDAgeFish
!-----------------------------------------------------------------------
!  zooplanktivorous fish respiration and excretion
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
!  benthivorous fish respiration and excretion
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
!  zooplanktivorous fish mortality
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
!  benthivorous fish mortality
!-----------------------------------------------------------------------
!  fish_mortality_incl._environmental_correction
!  PCLake_osis:sDFiAd,in g/m^2
   tDMortFiAd = self%kMortFiAd * sDFiAd +(1.0_rk - aDSatFiAd) * tDEnvFiAd
!  mortality_of_FiAd
!  PCLake_osis:sPFiAd,in g/m^2
   tPMortFiAd = rPDFiAd * tDMortFiAd
!  mortality_of_FiAd
!  PCLake_osis:sNFiAd,in g/m^2
   tNMortFiAd = rNDFiAd * tDMortFiAd
!-----------------------------------------------------------------------
!  fish egestion
!-----------------------------------------------------------------------
!  egestion_of_fish,zooplanktivorous fish
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
!-----------------------------------------------------------------------
!  zooplanktivorous fish predated by piscivirious fish
!-----------------------------------------------------------------------
!  young_fish_consumption_by_Pisc_DW
!  PCLake_osis:sDPisc,aDFish,in g/m^2
   tDConsFiJvPisc = sDFiJv / aDFish * tDConsPisc
!  young_fish_consumption_by_Pisc
   tPConsFiJvPisc = rPDFiJv * tDConsFiJvPisc
!  young_fish_consumption_by_Pisc
   tNConsFiJvPisc = rNDFiJv * tDConsFiJvPisc
!-----------------------------------------------------------------------
!  benthivorous fish predated by piscivirious fish
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
!  Fish manipulation, if turned on, operated on g/m2/s base
!-----------------------------------------------------------------------
!  If benthivorous fish manipulation tured on
   If(self%Manipulate_FiAd) then
     tDManFiAd= sDFiAd*rManFiAd/secs_pr_day
   else
     tDManFiAd = 0.0_rk
   endif
!  If zooplanktivorous manipulation tured on
   If(self%Manipulate_FiJv) then
      tDManFiJv= sDFiJv*rManFiJv/secs_pr_day
   else
      tDManFiJv = 0.0_rk
   endif
!  If piscivorous fish manipulation tured on
   If(self%Manipulate_Pisc) then
      tDManPisc= sDPisc*rManPisc/secs_pr_day
   else
      tDManPisc = 0.0_rk
   endif
!  Change in N, P rate due to manipulation
   tNManFiAd = rNDFiAd*tDManFiAd
   tNManFiJv = rNDFiJv*tDManFiJv
   tPManFiAd = rPDFiAd*tDManFiAd
   tPManFiJv = rPDFiJv*tDManFiJv
!-----------------------------------------------------------------------
!  total flux of Fish change to state variables
!-----------------------------------------------------------------------
!  total_fish_flux_of_DW_in_Young_fish
   tDFiJv = tDMigrFiJv + tDReprFish + tDAssFiJv - tDRespFiJv - tDMortFiJv - tDConsFiJvPisc - tDAgeFish - tDManFiJv
!  temperal solution, vertial averaged
   wDFiJv=tDFiJv/sDepthW
!  total_fish_flux_of_P_in_Young_fish
   tPFiJv = tPMigrFiJv + tPReprFish  + tPAssFiJv - tPExcrFiJv - tPMortFiJv - tPConsFiJvPisc - tPAgeFish - tNManFiJv
!  temperal solution, vertial averaged
   wPFiJv = tPFiJv/sDepthW
!  total_fish_flux_of_N_in_Young_fish
   tNFiJv = tNMigrFiJv + tNReprFish + tNAssFiJv - tNExcrFiJv - tNMortFiJv - tNConsFiJvPisc - tNAgeFish - tPManFiJv
!  temperal solution, vertial averaged
   wNFiJv= tNFiJv/ sDepthW
!  total_fish_flux_of_DW_in_Adult_fish
   tDFiAd = tDMigrFiAd - tDRespFiAd - tDReprFish - tDConsFiAdPisc + tDAgeFish- tDManFiAd - tDMortFiAd
!  temperal solution, vertial averaged
   wDFiAd= tDFiAd/ sDepthW
!  total_fish_flux_of_P_in_Adult_fish
   tPFiAd = tPMigrFiAd  - tPExcrFiAd - tPMortFiAd - tPReprFish - tPConsFiAdPisc + tPAgeFish - tPManFiAd
!  temperal solution, vertial averaged
   wPFiAd= tPFiAd/ sDepthW
!  total_fish_flux_of_N_in_Adult_fish
   tNFiAd = tNMigrFiAd - tNExcrFiAd - tNMortFiAd - tNReprFish - tNConsFiAdPisc + tNAgeFish - tNManFiAd
!  temperal solution, vertial averaged
   wNFiAd= tNFiAd/ sDepthW
!  total_fish_flux_of_DW_in_predatory_fish
   tDPisc = tDMigrPisc + tDAssPisc - tDRespPisc - tDMortPisc - tDManPisc
   wDPisc=tDPisc/sDepthW
!=======================================================================
!  fish processes relating to other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in water
!-----------------------------------------------------------------------
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
!  total_fish_flux_of_N_in_ammonium_in_water_in_lake_water
   wNFishNH4W = (tNExcrFiJv + tNExcrFiAd + tNEgesFiJvNH4 + tNMortFishNH4 + &
   & tNExcrPisc + tNEgesPiscNH4 + tNMortPiscNH4)/sDepthW
!-----------------------------------------------------------------------
!  Update PO4 in water
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
!  total_fish_flux_of_P_in_SRP_in_water_in_lake_water
   tPFishPO4W= (tPExcrFiJv + tPExcrFiAd + tPEgesFiJvPO4 + tPMortFishPO4 + tPExcrPisc&
   & + tPEgesPiscPO4 + tPMortPiscPO4)
   wPFishPO4W = tPFishPO4W/sDepthW
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
!  total_fish_flux_of_DW_in_Detritus_in_lake_water
   wDFishDetW = (tDEgesFiJv + tDMortFishDet + tDEgesPisc + tDMortPiscDet)/sDepthW
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
!  total_fish_flux_of_N_in_Detritus_in_lake_water
   wNFishDetW = (tNEgesFiJvDet + tNMortFishDet + tNEgesPiscDet + tNMortPiscDet)/sDepthW
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
!  total_fish_flux_of_P_in_Detritus_in_lake_water
   wPFishDetW = (tPEgesFiJvDet + tPMortFishDet + tPEgesPiscDet + tPMortPiscDet)/sDepthW
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
   _SET_ODE_(self%id_sDFiJv,wDFiJv)
   _SET_ODE_(self%id_sPFiJv,wPFiJv)
   _SET_ODE_(self%id_sNFiJv,wNFiJv)
   _SET_ODE_(self%id_sDFiAd,wDFiAd)
   _SET_ODE_(self%id_sPFiAd,wPFiAd)
   _SET_ODE_(self%id_sNFiAd,wNFiAd)
   _SET_ODE_(self%id_sDPisc,wDPisc)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
!  update abiotic variables in water
   _SET_ODE_(self%id_NH4poolW,  wNFishNH4W)
   _SET_ODE_(self%id_PO4poolW,  wPFishPO4W)
   _SET_ODE_(self%id_DDetpoolW, wDFishDetW)
   _SET_ODE_(self%id_NDetpoolW, wNFishDetW)
   _SET_ODE_(self%id_PDetpoolW, wPFishDetW)
   _SET_ODE_(self%id_DFoodZoo,  -tDConsFiJv/sDepthW)
   _SET_ODE_(self%id_NFoodZoo,  -tNConsFiJv/sDepthW)
   _SET_ODE_(self%id_PFoodZoo,  -tPConsFiJv/sDepthW)
!-----------------------------------------------------------------------
!  output diagnostic variables for external links
!-----------------------------------------------------------------------
!  Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_aNPisc,aNPisc)
   _SET_DIAGNOSTIC_(self%id_aPPisc,aPPisc)

!  output diagnostic variables for modular fluxes
   _SET_DIAGNOSTIC_(self%id_wDFiJv,    wDFiJv*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPFiJv,    wPFiJv*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNFiJv,    wNFiJv*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDPisc,    wDPisc*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDFishZoo, -tDConsFiJv/sDepthW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNFishZoo, -tNConsFiJv/sDepthW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPFishZoo, -tPConsFiJv/sDepthW*86400.0_rk)
!  feh: This temperal solution
!  feh: all the variables are connected to tMortFiAd, which has dependency on
!  external variables, so can not be updated at the first time step, or 37*4
!  calc. step. therefore, set this little higher
   if (n .GE.1000) then
      _SET_DIAGNOSTIC_(self%id_wDFiAd,     wDFiAd*86400.0_rk)
      _SET_DIAGNOSTIC_(self%id_wPFiAd,     wPFiAd*86400.0_rk)
      _SET_DIAGNOSTIC_(self%id_wNFiAd,     wNFiAd*86400.0_rk)
      _SET_DIAGNOSTIC_(self%id_wNFishNH4W, wNFishNH4W*86400.0_rk)
      _SET_DIAGNOSTIC_(self%id_wPFishPO4W, tPFishPO4W/sDepthW*86400.0_rk)
      _SET_DIAGNOSTIC_(self%id_wDFishDetW, wDFishDetW*86400.0_rk)
      _SET_DIAGNOSTIC_(self%id_wNFishDetW, wNFishDetW*86400.0_rk)
      _SET_DIAGNOSTIC_(self%id_wPFishDetW, wPFishDetW*86400.0_rk)

   endif
   n=n+1
!-----------------------------------------------------------------------
!  Updated changed fish biomass for biomanipulation
!-----------------------------------------------------------------------
!  If benthivorous fish manipulation tured on
   If(self%Manipulate_FiAd) then
     _SET_ODE_(self%id_ChangedFiAd,tDManFiAd)
   endif
!  If zooplanktivorous manipulation tured on
   If(self%Manipulate_FiJv) then
     _SET_ODE_(self%id_ChangedFiJv,tDManFiJv)
   endif
!  If Piscivorous fish manipulation tured ontPEgesFiAdDet
   If(self%Manipulate_Pisc) then
     _SET_ODE_(self%id_ChangedPisc,tDManPisc)
   endif


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
   class (type_au_pclake_fish), intent(in) :: self
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


   end module au_pclake_fish

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

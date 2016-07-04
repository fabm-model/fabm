#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module au_pclake_zooplankton
!
! !DESCRIPTION:
!  The au_pclake_zooplankton module describes the state variables regarding zoo-
!  -planktons
!  Therefore, the state variables and their related its local processes are:
!  sDZoo,sPZoo,sNZoo,(zooplankton in dry-weight, nitrogen element and phosphorus
!  element respectively)
!  units: gDW/m**3, gP/m**3,gP/m**3
!  local processes:assimilation,respiration(only for sDZoo),excretion(only for sNZoo,sPZoo),
!  This module also describes the processes which influence the state variables registered in
!  other modules, including: (aPhytW stands for all groups of phytoplankton)
!  Zooplankton grazing on phytoplankton in water column: aDPhytW==>aDZoo,
!  aNPhytW==>aNZoo,aPPhytW==>aPZoo
!  Zooplankton grazing on detritus in the water column: sDDetW==>aDZoo,
!  sNDetW==>aNZoo,sPDetW==>aPZoo.
!  Detritus morted by zooplankton.
!  sDDetW<==sDZoo,sNDetW<==sNZoo,sPDetW<==sPZoo.
!  nutrients excreted by zooplankton.
!   sNH4W<==sNZoo,sPO4<==sPZoo,sNO3W<==0.0(special)
! !USES:
   use fabm_types
   use au_pclake_utility, ONLY:uFunTmBio

   implicit none

!  default: all is private.
      private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_au_pclake_zooplankton
!     local state variable identifiers
!     id_sDZoo,zooplankton concentration in dry-weight, gDW/m**3
!     id_sPZoo,zooplankton concentration in nitrogen element, gN/m**3
!     id_sNZoo,zooplankton concentration in phosphorus element, gP/m**3
      type (type_state_variable_id)            :: id_sDZoo,id_sPZoo,id_sNZoo
!     diagnostic variables for dependencies(without output)
!     diagnostic variables for modular fluxes
      type (type_diagnostic_variable_id)       :: id_wDZoo,id_wNZoo,id_wPZoo
      type (type_diagnostic_variable_id)       :: id_wNZooNO3W,id_wPZooPO4W,id_wDZooDetW
      type (type_diagnostic_variable_id)       :: id_wNZooDetW,id_wPZooDetW,id_wSiZooDetW
      type (type_diagnostic_variable_id)       :: id_wDZooDiatW,id_wNZooDiatW,id_wPZooDiatW
      type (type_diagnostic_variable_id)       :: id_wDZooGrenW,id_wNZooGrenW,id_wPZooGrenW
      type (type_diagnostic_variable_id)       :: id_wDZooBlueW,id_wNZooBlueW,id_wPZooBlueW
!     state dependencies identifiers
      type (type_state_variable_id)            :: id_DfoodDiat,id_DfoodGren,id_DfoodBlue,id_DDetpoolW
      type (type_state_variable_id)            :: id_NfoodDiat,id_NfoodGren,id_NfoodBlue,id_NDetpoolW
      type (type_state_variable_id)            :: id_PfoodDiat,id_PfoodGren,id_PfoodBlue,id_PDetpoolW,id_SiDetpoolW
      type (type_state_variable_id)            :: id_NH4poolW,id_NO3poolW,id_PO4poolW
!     environmental dependencies
      type (type_dependency_id)                :: id_uTm ,id_dz
      type (type_horizontal_dependency_id)     :: id_sDepthW
      type (type_global_dependency_id)         :: id_Day
!!    Model parameters
!     parameter for temperature
      real(rk)      :: cSigTmZoo,cTmOptZoo
!     parameters for zooplankton
      real(rk)      :: cDCarrZoo,kMortZoo,kDRespZoo,cPrefDiat
      real(rk)      :: cPrefGren,cPrefBlue,cPrefDet,hFilt
      real(rk)      :: fDAssZoo,cFiltMax
      real(rk)      :: cPDZooRef,cNDZooRef
      real(rk)      :: fDissEgesZoo,fDissMortZoo,cSiDDiat
!     nutrient ratios parameter
      real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
      real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax
!     Fish manipulation parameters, switch for turned on/off fish manipulation
      logical    :: Manipulate_FiAd, Manipulate_FiJv, Manipulate_Pisc
!     minimum state variable values
      real(rk)   :: cDZooMin
   contains

!     Model procedures
      procedure :: initialize
      procedure :: do
      !procedure :: do_bottom
      procedure :: get_light_extinction

   end type type_au_pclake_zooplankton

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk)           :: Pi=3.14159265358979_rk

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
   class (type_au_pclake_zooplankton), intent(inout), target :: self
   integer,                     intent(in)            :: configunit

!EOP
!-----------------------------------------------------------------------
!BOC

!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cSigTmZoo,     'cSigTmZoo',     '�C',        'temperature constant zooplankton(sigma in Gaussian curve)',                      default=13.0_rk)
   call self%get_parameter(self%cTmOptZoo,     'cTmOptZoo',     '�C',        'optimum temp. zooplankton',                                                      default=25.0_rk)
   call self%get_parameter(self%cDCarrZoo,     'cDCarrZoo',     'gm-3',      'carrying capacity of zooplankton',                                               default=25.0_rk)
   call self%get_parameter(self%kMortZoo,      'kMortZoo',      'd-1',       'mortality_constant_herb.zooplankton',                                            default=0.04_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespZoo,     'kDRespZoo',     'd-1',       'maintenance respiration constant herb.zooplankton',                              default=0.15_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cPrefDiat,     'cPrefDiat',     '[-]',       'selection factor for Diatoms',                                                   default=0.75_rk)
   call self%get_parameter(self%cPrefGren,     'cPrefGren',     '[-]',       'selection factor for Greens',                                                    default=0.75_rk)
   call self%get_parameter(self%cPrefBlue,     'cPrefBlue',     '[-]',       'selection factor for Bluegreens Cal.',                                           default=0.125_rk)
   call self%get_parameter(self%cPrefDet,      'cPrefDet',      '[-]',       'selection factor for detritus',                                                  default=0.25_rk)
   call self%get_parameter(self%hFilt,         'hFilt',         'gDW m-3',   'half-sat. food conc. for filtering',                                             default=1.0_rk)
   call self%get_parameter(self%fDAssZoo,      'fDAssZoo',      '[-]',       'DW-assimilation efficiency of herb. zooplankton',                                default=0.35_rk)
   call self%get_parameter(self%cFiltMax,      'cFiltMax',      'ltr/mgDW/d','maximum filtering rate',                                                         default=4.5_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cPDZooRef,     'cPDZooRef',     'mgP/mgDW',  'reference P/C-ratio herb. zooplankton',                                          default=0.01_rk)
   call self%get_parameter(self%cNDZooRef,     'cNDZooRef',     'mgN/mgDW',  'reference N/C-ratio herb. zooplankton',                                          default=0.07_rk)
   call self%get_parameter(self%fDissEgesZoo,  'fDissEgesZoo',  '[-]',       'soluble_nutrient_fraction_of_by_herb.zoopl._egested_foodcVSetDiat',              default=0.25_rk)
   call self%get_parameter(self%fDissMortZoo,  'fDissMortZoo',  '[-]',       'soluble_nutrient_fraction_of_died_zooplanktoncVSetGren',                         default=0.1_rk)
   call self%get_parameter(self%cSiDDiat,      'cSiDDiat',      'mgSi/mgDW', 'Si/DW_ratio_of_daitomscVSetBlue',                                                default=0.15_rk)
   call self%get_parameter(self%cNDDiatMin,    'cNDDiatMin',    'mgN/mgDW',  'minimum N/day ratio Diatoms',                                                    default=0.01_rk)
   call self%get_parameter(self%cPDDiatMin,    'cPDDiatMin',    'mgP/mgDW',  'minimum P/day ratio Diatoms',                                                    default=0.0005_rk)
   call self%get_parameter(self%cNDGrenMin,    'cNDGrenMin',    'mgN/mgDW',  'minimum N/day ratio greens',                                                     default=0.02_rk)
   call self%get_parameter(self%cPDGrenMin,    'cPDGrenMin',    'mgP/mgDW',  'minimum P/day ratio greens',                                                     default=0.0015_rk)
   call self%get_parameter(self%cNDBlueMin,    'cNDBlueMin',    'mgN/mgDW',  'minimum N/day ratio Bluegreens',                                                 default=0.03_rk)
   call self%get_parameter(self%cPDBlueMin,    'cPDBlueMin',    'mgP/mgDW',  'minimum P/day ratio Bluegreens',                                                 default=0.0025_rk)
   call self%get_parameter(self%cNDBlueMax,    'cNDBlueMax',    'mgN/mgDW',  'max. N/day ratio Bluegreens',                                                    default=0.15_rk)
   call self%get_parameter(self%cNDDiatMax,    'cNDDiatMax',    'mgN/mgDW',  'max. N/day ratio Diatoms',                                                       default=0.05_rk)
   call self%get_parameter(self%cNDGrenMax,    'cNDGrenMax',    'mgN/mgDW',  'max. N/day ratio greens',                                                        default=0.1_rk)
   call self%get_parameter(self%cPDBlueMax,    'cPDBlueMax',    'mgP/mgDW',  'max. P/day ratio blue-greens',                                                   default=0.025_rk)
   call self%get_parameter(self%cPDDiatMax,    'cPDDiatMax',    'mgP/mgDW',  'max. P/day ratio Diatoms',                                                       default=0.005_rk)
   call self%get_parameter(self%cPDGrenMax,    'cPDGrenMax',    'mgP/mgDW',  'max. P/day ratio greens',                                                        default=0.015_rk)
!  the user defined minumun value for state variables
   call self%get_parameter(self%cDZooMin,      'cDZooMin',      'gDW/m3',    'minimun zooplankton biomass in system',                                          default=0.00001_rk)
   
   
!  Register local state variable
!  zooplankton
   call self%register_state_variable(self%id_sDZoo,'sDZoo','g m-3','zooplankton biomass',     &
                                    initial_value=0.05_rk,minimum=self%cDZooMin,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPZoo,'sPZoo','g m-3','zooplankton phosphorus content',     &
                                    initial_value=0.0005_rk,minimum=self%cDZooMin * self%cPDZooRef,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNZoo,'sNZoo','g m-3','zooplankton nitrogen content',     &
                                    initial_value=0.0035_rk,minimum=self%cDZooMin * self%cNDZooRef,no_river_dilution=.TRUE.)
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNZoo)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPZoo)
!  register state variables dependencies
   call self%register_state_dependency(self%id_DfoodDiat,     'Diatom_as_food_DW',      'g m-3', 'Diatom_as_food_DW')
   call self%register_state_dependency(self%id_DfoodGren,     'Green_as_food_DW',       'g m-3', 'Green_as_food_DW')
   call self%register_state_dependency(self%id_DfoodBlue,     'Blue_as_food_DW',        'g m-3', 'Blue_as_food_DW')
   call self%register_state_dependency(self%id_NfoodDiat,     'Diatom_as_food_N',       'g m-3', 'Diatom_as_food_N')
   call self%register_state_dependency(self%id_NfoodGren,     'Green_as_food_N',        'g m-3', 'Green_as_food_N')
   call self%register_state_dependency(self%id_NfoodBlue,     'Blue_as_food_N',         'g m-3', 'Blue_as_food_N')
   call self%register_state_dependency(self%id_PfoodDiat,     'Diatom_as_food_P',       'g m-3', 'Diatom_as_food_P')
   call self%register_state_dependency(self%id_PfoodGren,     'Green_as_food_P',        'g m-3', 'Green_as_food_P')
   call self%register_state_dependency(self%id_PfoodBlue,     'Blue_as_food_P',         'g m-3', 'Blue_as_food_P')
   call self%register_state_dependency(self%id_DDetpoolW,     'Detritus_DW_pool_water', 'g m-3', 'Detritus_DW_pool_water')
   call self%register_state_dependency(self%id_NDetpoolW,     'Detritus_N_pool_water',  'g m-3', 'Detritus_N_pool_water')
   call self%register_state_dependency(self%id_PDetpoolW,     'Detritus_P_pool_water',  'g m-3', 'Detritus_P_pool_water')
   call self%register_state_dependency(self%id_SiDetpoolW,    'Detritus_Si_pool_water', 'g m-3', 'Detritus_Si_pool_water')
   call self%register_state_dependency(self%id_NH4poolW,      'NH4_pool_water',         'g m-3', 'NH4_pool_water')
   call self%register_state_dependency(self%id_NO3poolW,      'NO3_pool_water',         'g m-3', 'NO3_pool_water')
   call self%register_state_dependency(self%id_PO4poolW,      'PO4_pool_water',         'g m-3', 'PO4_pool_water')
!  register diagnostic variables for modular fluxes
   call self%register_diagnostic_variable(self%id_wDZoo,      'wDZoo',                   'g m-3 s-1', 'zooplankton_DZoo_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZoo,      'wNZoo',                   'g m-3 s-1', 'zooplankton_NZoo_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZoo,      'wPZoo',                   'g m-3 s-1', 'zooplankton_PZoo_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooNO3W,  'wNZooNO3W',               'g m-3 s-1', 'zooplankton_NO3W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooPO4W,  'wPZooPO4W',               'g m-3 s-1', 'zooplankton_PO4W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDZooDetW,  'wDZooDetW',               'g m-3 s-1', 'zooplankton_DDetW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooDetW,  'wNZooDetW',               'g m-3 s-1', 'zooplankton_NDetW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooDetW,  'wPZooDetW',               'g m-3 s-1', 'zooplankton_PDetW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wSiZooDetW, 'wSiZooDetW',              'g m-3 s-1', 'zooplankton_SiDetW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDZooDiatW, 'wDZooDiatW',              'g m-3 s-1', 'zooplankton_DDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooDiatW, 'wNZooDiatW',              'g m-3 s-1', 'zooplankton_NDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooDiatW, 'wPZooDiatW',              'g m-3 s-1', 'zooplankton_PDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDZooGrenW, 'wDZooGrenW',              'g m-3 s-1', 'zooplankton_DGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooGrenW, 'wNZooGrenW',              'g m-3 s-1', 'zooplankton_NGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooGrenW, 'wPZooGrenW',              'g m-3 s-1', 'zooplankton_PGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDZooBlueW, 'wDZooBlueW',              'g m-3 s-1', 'zooplankton_DBlue_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooBlueW, 'wNZooBlueW',              'g m-3 s-1', 'zooplankton_NBlue_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooBlueW, 'wPZooBlueW',              'g m-3 s-1', 'zooplankton_PBlue_change',  output=output_instantaneous)

!  register environmental dependencies
   call self%register_dependency(self%id_uTm,    standard_variables%temperature)
   call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_dz,     standard_variables%cell_thickness)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)

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
   class (type_au_pclake_zooplankton), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_
!  LOCAL VARIABLES:
!  Carriers for environment dependencies
   real(rk)     :: uTm,Day,dz,sDepthW
!  carriers for local state variables
   real(rk)      :: sDZoo,sNZoo,sPZoo
!  carriers for exteral link state variables
   real(rk)      :: sDDiatW,sDGrenW,sDBlueW,sDDetW
   real(rk)      :: sNDiatW,sNGrenW,sNBlueW,sNDetW
   real(rk)      :: sPDiatW,sPGrenW,sPBlueW,sPDetW
!  nutrient ratios variables
   real(rk)      :: rNDDiatW,rNDGrenW,rNDBlueW,rNDDetW
   real(rk)      :: rPDDiatW,rPDGrenW,rPDBlueW,rPDDetW
   real(rk)      :: rPDZoo,rNDZoo
!  variables for temperature function
   real(rk)      :: uFunTmZoo
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
   real(rk)      :: wDZoo,wDAssZoo,wDRespZoo,wDMortZoo
!  variables for Zooplankton N fluxes
   real(rk)      :: wNZoo,wNAssZoo,wNExcrZoo,wNMortZoo
!  variables for Zooplankton P fluxes
   real(rk)      ::wPZoo,wPAssZoo,wPExcrZoo,wPMortZoo
!  variables for exchange of NH4
!  PCLake_Osis, /m^2
   real(rk)     :: wNZooNH4W,wNEgesZooNH4,wNEgesZoo,wNMortZooNH4
!  variables for exchange of NH3
   real(rk)     :: wNZooNO3W
!  variables for exchange of PO4
!  PCLake_Osis, /m^2
   real(rk)     :: wPZooPO4W,wPEgesZooPO4,wPEgesZoo,wPMortZooPO4
!  variables for exchange of Detritus DW
   real(rk)     :: wDZooDetW,wDEgesZoo
!  variables for exchange of Detritus N
   real(rk)     :: wNZooDetW,wNEgesZooDet,wNMortZooDet
!  variables for exchange of detritus P
   real(rk)     :: wPZooDetW,wPEgesZooDet,wPMortZooDet
!  variables for exchange of detritus Si
   real(rk)     :: wSiZooDetW,wSiConsDiatZoo
!  variables for exchange of diatoms
   real(rk)     :: wDZooDiatW,wNZooDiatW,wPZooDiatW
!  variables for exchange of green algae
   real(rk)     :: wDZooGrenW,wNZooGrenW,wPZooGrenW
!  variables for exchange of green algae
   real(rk)     :: wDZooBlueW,wNZooBlueW,wPZooBlueW
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_(self%id_sDZoo,sDZoo)
   _GET_(self%id_sNZoo,sNZoo)
   _GET_(self%id_sPZoo,sPZoo)
!-----------------------------------------------------------------------
!  Retrieve dependencies  value
!-----------------------------------------------------------------------
!  Retrieve state dependencies value
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
   _GET_(self%id_dz,dz)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
!-------------------------------------------------------------------------
!  Current local nutrients ratios in zooplankton(check the current state)
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
!-----------------------------------------------------------------------
!  temperature function
!-----------------------------------------------------------------------
!  temp._function_of_zooplankton
   uFunTmZoo =uFunTmBio(uTm,self%cSigTmZoo,self%cTmOptZoo)
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
!  total flux of state variables
!-----------------------------------------------------------------------
!  total_flux_of_DW_in_Herbivorous_zooplankton
   wDZoo = wDAssZoo - wDRespZoo - wDMortZoo
!  total_flux_of_N_in_Herbivorous_zooplankton
   wNZoo = wNAssZoo - wNExcrZoo - wNMortZoo
!  total_flux_of_P_in_Herbivorous_zooplankton
   wPZoo = wPAssZoo - wPExcrZoo - wPMortZoo
!=======================================================================
!  zooplankton process relating to other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in water
!-----------------------------------------------------------------------
!  soluble_N_egestion
   wNEgesZooNH4 = self%fDissEgesZoo*wNEgesZoo
!  soluble_N_mortality
   wNMortZooNH4 = self%fDissMortZoo*wNMortZoo
!  total_Zoo_flux_of_N_in_ammonium_in_water_in_lake_water
  wNZooNH4W = wNExcrZoo + wNEgesZooNH4 + wNMortZooNH4
!-----------------------------------------------------------------------
!  Update NO3 in water   (no NO3????)
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_N_in_nitrate_in_water_in_lake_water
   wNZooNO3W = 0.0_rk
!-----------------------------------------------------------------------
!  Update PO4 in water
!-----------------------------------------------------------------------
!  soluble_P_egestion
   wPEgesZooPO4 = self%fDissEgesZoo*wPEgesZoo
!  soluble_P_mortality
   wPMortZooPO4 = self%fDissMortZoo * wPMortZoo
!  total_Zoo_flux_of_P_in_SRP_in_water_in_lake_water
  wPZooPO4W = wPExcrZoo + wPEgesZooPO4 + wPMortZooPO4
!-----------------------------------------------------------------------
!  Update detrital DW in water
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_DW_in_Detritus_in_lake_water
   wDZooDetW = - wDConsDetZoo + wDEgesZoo + wDMortZoo
!-----------------------------------------------------------------------
!  Update detrital N in water
!-----------------------------------------------------------------------
!  detrital_N_mortality
   wNMortZooDet = wNMortZoo - wNMortZooNH4
!  detrital_N_egestion
   wNEgesZooDet = wNEgesZoo - wNEgesZooNH4
!  total_Zoo_flux_of_N_in_Detritus_in_lake_water
   wNZooDetW = - wNConsDetZoo + wNEgesZooDet + wNMortZooDet
!-----------------------------------------------------------------------
!  Update detrital P in water
!-----------------------------------------------------------------------
!  detrital_P_mortality
   wPMortZooDet = wPMortZoo - wPMortZooPO4
!  detrital_P_egestion
   wPEgesZooDet = wPEgesZoo - wPEgesZooPO4
!  total_Zoo_flux_of_P_in_Detritus_in_lake_water
   wPZooDetW = - wPConsDetZoo + wPEgesZooDet + wPMortZooDet
!-----------------------------------------------------------------------
!  Update detrital Si in water
!-----------------------------------------------------------------------
!  consumption_of_diatoms
   wSiConsDiatZoo = self%cSiDDiat * wDConsDiatZoo
!  total_Zoo_flux_of_silica_in_lake_water_detritus
   wSiZooDetW = wSiConsDiatZoo
!-----------------------------------------------------------------------
!  Update diatom state variables
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_DW_in_Diatoms_in_lake_water
   wDZooDiatW = - wDConsDiatZoo
!  total_Zoo_flux_of_N_in_Diatoms_in_lake_water
   wNZooDiatW = - wNConsDiatZoo
!  total_Zoo_flux_of_P_in_Diatoms_in_lake_water
   wPZooDiatW = - wPConsDiatZoo
!-----------------------------------------------------------------------
!  Update green algae state variables
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_DW_in_Greens_in_lake_water
   wDZooGrenW = - wDConsGrenZoo
!  total_Zoo_flux_of_N_in_Greens_in_lake_water
   wNZooGrenW = - wNConsGrenZoo
!  total_Zoo_flux_of_P_in_Greens_in_lake_water
   wPZooGrenW = - wPConsGrenZoo
!-----------------------------------------------------------------------
!  Update blue algae state variables
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_DW_in_Blue-greens_in_lake_water
   wDZooBlueW = - wDConsBlueZoo
!  total_Zoo_flux_of_N_in_Blue-greens_in_lake_water
   wNZooBlueW = - wNConsBlueZoo
!  total_Zoo_flux_of_P_in_Blue-greens_in_lake_water
   wPZooBlueW = - wPConsBlueZoo
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
   _SET_ODE_(self%id_sDZoo,wDZoo)
   _SET_ODE_(self%id_sNZoo,wNZoo)
   _SET_ODE_(self%id_sPZoo,wPZoo)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
!  update abiotic variables in water
   _SET_ODE_(self%id_NH4poolW,   wNZooNH4W)
   _SET_ODE_(self%id_NO3poolW,   wNZooNO3W)
   _SET_ODE_(self%id_PO4poolW,   wPZooPO4W)
   _SET_ODE_(self%id_DDetpoolW,  wDZooDetW)
   _SET_ODE_(self%id_NDetpoolW,  wNZooDetW)
   _SET_ODE_(self%id_PDetpoolW,  wPZooDetW)
   _SET_ODE_(self%id_SiDetpoolW, wSiZooDetW)
!  update phytoplankton in water
   _SET_ODE_(self%id_DfoodDiat,  wDZooDiatW)
   _SET_ODE_(self%id_NfoodDiat,  wNZooDiatW)
   _SET_ODE_(self%id_PfoodDiat,  wPZooDiatW)
   _SET_ODE_(self%id_DfoodGren,  wDZooGrenW)
   _SET_ODE_(self%id_NfoodGren,  wNZooGrenW)
   _SET_ODE_(self%id_PfoodGren,  wPZooGrenW)
   _SET_ODE_(self%id_DfoodBlue,  wDZooBlueW)
   _SET_ODE_(self%id_NfoodBlue,  wNZooBlueW)
   _SET_ODE_(self%id_PfoodBlue,  wPZooBlueW)
!  output diagnostic variables for modular fluxes
   _SET_DIAGNOSTIC_(self%id_wDZoo,      wDZoo*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNZoo,      wNZoo*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPZoo,      wPZoo*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNZooNO3W,  wNZooNO3W*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPZooPO4W,  wPZooPO4W*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDZooDetW,  wDZooDetW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNZooDetW,  wNZooDetW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPZooDetW,  wPZooDetW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wSiZooDetW, wSiZooDetW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDZooDiatW, wDZooDiatW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNZooDiatW, wNZooDiatW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPZooDiatW, wPZooDiatW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDZooGrenW, wDZooGrenW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNZooGrenW, wNZooGrenW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPZooGrenW, wPZooGrenW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDZooBlueW, wDZooBlueW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNZooBlueW, wNZooBlueW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPZooBlueW, wPZooBlueW*86400.0_rk)



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
   class (type_au_pclake_zooplankton), intent(in) :: self
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


   end module au_pclake_zooplankton

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

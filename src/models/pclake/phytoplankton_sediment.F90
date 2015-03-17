#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_phytoplankton_sediment
!
! !DESCRIPTION:
! Module description
!  The pclake_abiotic_water module describes the processes related to settled phytoplankton in the sediment.
!  Three groups of phytoplankton are described here: Diatom, green algae and 
!  cyanobacteria(blue algae). Each group is described in three elements, dry-weight
!  nitrogen and phosphrus. Silica contration in diatom is not a state variables here
!  but a diagnostic instead, since model assumes diatom have fixed Si/D ration,i.e. 0.1
!  The state variables and their involved processes are:
!  State variables:           sDDiatS, sNDiatS, sPDiatS(diatom concentration in DW, N and P element)
!                             sDGrenS, sNGrenS, sPGrenS(grenn algae concentration in DW, N and P element)
!                             sDBlueS, sNBlueS, sPBlueS(cyanobacteria concentration in DW, N and P element)
!  units( for all the groups):gDW/m**2, gN/m**2, gP/m**2 in three elements respectively,
!  involved processes( for all the groups): 
!  respiration(only for sDDiatW,sDGrenW,sDBlueW), excretion(only for sNDiatW,sNGrenW,sNBlueW,
!  sPDiatW,sPGrenW,sPBlueW), mortality(for all state variables)
!  This module also discribes the processes which influence the state variables registered in
!  other modules, including:(aPhytS stands for all groups of settled phytoplankton)
!  nutrients excreted by settled phytoplankton: sNH4S<==aPhytS,sPO4S<==aPhytS,(no sNO3S excreted),sSiO2W<==oSiDiatS
!  detritus morted by settled phytoplankton: sDDetS<==aPhytS,sPDetS<==aPhytS,sNDetS<==aPhytS,sSiDetS<==oSiDiatS
!  This module also provide important diagnostic variable will be used in other modules, including:
!  Sediment detritus change, tDPrimDetS, used bymodule:auxilary
! !USES:
   use fabm_types
   use pclake_utility, ONLY:uFunTmBio

   
   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_phytoplankton_sediment
!     local state variable identifers
!     local state variable identifers
!     id_sDDiatS,id_sDGrenS,id_sDBlueS: phytoplankton concentration in dry-weight, gDW/m**2
!     id_sPDiatS,id_sPGrenS,id_sPBlueS: phytoplankton concentration in nitrogen element, gN/m**2
!     id_sNDiatS,id_sNGrenS,id_sNBlueS: phytoplankton concentration in phosphorus element, gP/m**2
      type (type_bottom_state_variable_id) :: id_sDDiatS,id_sPDiatS,id_sNDiatS
      type (type_bottom_state_variable_id) :: id_sDGrenS,id_sPGrenS,id_sNGrenS
      type (type_bottom_state_variable_id) :: id_sDBlueS,id_sPBlueS,id_sNBlueS
!     diagnostic variables for local output
!     id_oSiDiatS, diatom concentration in silica element, gSi/m**2
!     id_rPDPhytS, P/D ratio of settled phytoplankton
!     id_rNDPhytS, N/D rartio of settled phytoplankton
      type (type_horizontal_diagnostic_variable_id)       :: id_oSiDiatS
      type (type_horizontal_diagnostic_variable_id)       :: id_rPDPhytS,id_rNDPhytS
!     diagnostic variables for dependencies(without output)
      type (type_horizontal_diagnostic_variable_id)       :: id_tDPrimDetS
!     state dependencies identifers
     type (type_bottom_state_variable_id) :: id_PO4poolS,id_NO3poolS,id_NH4poolS
     type (type_bottom_state_variable_id) :: id_DDetpoolS,id_NDetpoolS,id_PDetpoolS,id_SiDetpoolS
     type (type_state_variable_id) :: id_SiO2poolW
!     environmental dependencies
      type (type_dependency_id)                :: id_uTm
!     Model parameters
!     temperature parameters
      real(rk)   :: cSigTmDiat,cTmOptDiat
      real(rk)   :: cSigTmBlue,cTmOptBlue,cSigTmGren,cTmOptGren
!     diatoms related parameters
      real(rk)   :: kMortDiatS,kDRespDiat,kMortGrenS,kDRespGren,kMortBlueS,kDRespBlue
!     nutrient ratios parameter
      real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
      real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax
!     exchange process related parameters
      real(rk)   :: fDissMortPhyt,cSiDDiat
      
      
   contains

     procedure initialize
     procedure do_bottom
     
     
     
      end type type_pclake_phytoplankton_sediment


!  private data memebers(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero=0.000000000000000000000000000000001_rk
!  Lowest phytoplankton value in sediment
   real(rk),parameter :: PhyZeroS=0.0000001_rk
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
! !INPUT PARAMETERS:
   class (type_pclake_phytoplankton_sediment), intent(inout), target :: self
   integer,                     intent(in)            :: configunit

!EOP
!-----------------------------------------------------------------------
!BOC

!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cSigTmDiat,   'cSigTmDiat',   '°C',       'temperature constant diatoms (sigma in Gaussian curve)',    default=20.0_rk)
   call self%get_parameter(self%cTmOptDiat,   'cTmOptDiat',   '°C',       'optimum temp. diatoms',                                     default=18.0_rk)
   call self%get_parameter(self%cSigTmBlue,   'cSigTmBlue',   '°C',       'temperature constant blue-greens (sigma in Gaussian curve)',default=12.0_rk)
   call self%get_parameter(self%cTmOptBlue,   'cTmOptBlue',   '°C',       'optimum temp. blue-greens',                                 default=25.0_rk)
   call self%get_parameter(self%cSigTmGren,   'cSigTmGren',   '°C',       'temperature constant greens (sigma in Gaussian curve)',     default=15.0_rk)
   call self%get_parameter(self%cTmOptGren,   'cTmOptGren',   '°C',       'optimum temp. of greens',                                   default=25.0_rk)
   call self%get_parameter(self%kMortDiatS,   'kMortDiatS',   'd-1',      'mortality constant of sed. Diatoms',                        default=0.05_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespDiat,   'kDRespDiat',   'd-1',      'maintenance respiration constant diatoms',                  default=0.1_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cNDDiatMin,   'cNDDiatMin',   'mgN/mgDW', 'minimum N/day ratio Diatoms',                               default=0.01_rk)
   call self%get_parameter(self%cPDDiatMin,   'cPDDiatMin',   'mgP/mgDW', 'minimum P/day ratio Diatoms',                               default=0.0005_rk)
   call self%get_parameter(self%kMortGrenS,   'kMortGrenS',   'd-1',      'mortality constant greens',                                 default=0.05_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespGren,   'kDRespGren',   'd-1',      'maintenance respiration constant greens ',                  default=0.075_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cNDGrenMin,   'cNDGrenMin',   'mgN/mgDW', 'minimum N/day ratio greens',                                default=0.02_rk)
   call self%get_parameter(self%cPDGrenMin,   'cPDGrenMin',   'mgP/mgDW', 'minimum P/day ratio greens',                                default=0.0015_rk)
   call self%get_parameter(self%kMortBlueS,   'kMortBlueS',   'd-1',      'mortality constant Bluegreens',                             default=0.2_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespBlue,   'kDRespBlue',   'd-1',      'maintenance respiration constant blue-greens',              default=0.03_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cNDBlueMin,   'cNDBlueMin',   'mgN/mgDW', 'minimum N/day ratio Bluegreens',                            default=0.03_rk)
   call self%get_parameter(self%cPDBlueMin,   'cPDBlueMin',   'mgP/mgDW', 'minimum P/day ratio Bluegreens',                            default=0.0025_rk)
   call self%get_parameter(self%fDissMortPhyt,'fDissMortPhyt','[-]',      'soluble nutrient fraction of died Algae',                   default=0.2_rk)
   call self%get_parameter(self%cSiDDiat,     'cSiDDiat',     'mgSi/mgDW','Si/D ratio of diatoms',                                     default=0.15_rk)
   call self%get_parameter(self%cNDBlueMax,   'cNDBlueMax',   'mgN/mgDW', 'max. N/day ratio Bluegreens',                               default=0.15_rk)
   call self%get_parameter(self%cNDDiatMax,   'cNDDiatMax',   'mgN/mgDW', 'max. N/day ratio Diatoms',                                  default=0.05_rk)
   call self%get_parameter(self%cNDGrenMax,   'cNDGrenMax',   'mgN/mgDW', 'max. N/day ratio greens',                                   default=0.1_rk)
   call self%get_parameter(self%cPDBlueMax,   'cPDBlueMax',   'mgP/mgDW', 'max. P/day ratio blue-greens',                              default=0.025_rk)
   call self%get_parameter(self%cPDDiatMax,   'cPDDiatMax',   'mgP/mgDW', 'max. P/day ratio Diatoms',                                  default=0.005_rk)
   call self%get_parameter(self%cPDGrenMax,   'cPDGrenMax',   'mgP/mgDW', 'max. P/day ratio greens',                                   default=0.015_rk)

!  Register local state variable
   call self%register_state_variable(self%id_sDDiatS,'sDDiatS','g m-2','diatom_D in sediment',     &
                                    initial_value=0.001_rk,minimum= PhyZeroS)
   call self%register_state_variable(self%id_sPDiatS,'sPDiatS','g m-2','diatom_P in sediment',     &
                                    initial_value=0.00001_rk,minimum= PhyZeroS)
   call self%register_state_variable(self%id_sNDiatS,'sNDiatS','g m-2','diatom_N in sediment',     &
                                    initial_value=0.0001_rk,minimum= PhyZeroS)
   call self%register_state_variable(self%id_sDGrenS,'sDGrenS','g m-2','green_D in sediment',     &
                                    initial_value=0.001_rk,minimum= PhyZeroS)
   call self%register_state_variable(self%id_sPGrenS,'sPGrenS','g m-2','green_P in sediment',     &
                                    initial_value=0.00001_rk,minimum= PhyZeroS)
   call self%register_state_variable(self%id_sNGrenS,'sNGrenS','g m-2','green_N in sediment',     &
                                    initial_value=0.0001_rk,minimum=PhyZeroS)
   call self%register_state_variable(self%id_sDBlueS,'sDBlueS','g m-2','blue_D in sediment',     &
                                    initial_value=0.001_rk,minimum=PhyZeroS)
   call self%register_state_variable(self%id_sPBlueS,'sPBlueS','g m-2','blue_P in sediment',     &
                                    initial_value=0.00001_rk,minimum=PhyZeroS)
   call self%register_state_variable(self%id_sNBlueS,'sNBlueS','g m-2','blue_N in sediment',     &
                                    initial_value=0.0001_rk,minimum=PhyZeroS)
!  register diagnostic variables
   call self%register_diagnostic_variable(self%id_oSiDiatS,  'oSiDiatS',  'g m-2 s-1','oSiDiatS',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDPhytS,  'rPDPhytS',  '[-]',      'rPDPhytS',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDPhytS,  'rNDPhytS',  '[-]',      'rNDPhytS',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDPrimDetS,'tDPrimDetS','g m-2 s-1','tDPrimDetS',output=output_none)
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNDiatS)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNGrenS)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNBlueS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPDiatS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPGrenS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPBlueS)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_oSiDiatS)
   
!  regirster state variables dependencies
   call self%register_state_dependency(self%id_PO4poolS,  'sPO4_pool_sediment',       'g m-2', 'sPO4_pool_sediment')
   call self%register_state_dependency(self%id_NO3poolS,  'sNO3_pool_sediment',       'g m-2', 'sNO3_pool_sediment')
   call self%register_state_dependency(self%id_NH4poolS,  'sNH4_pool_sediment',       'g m-2', 'sNH4_pool_sediment')
   call self%register_state_dependency(self%id_DDetpoolS, 'detritus_DW_pool_sediment','g m-2', 'detritus_DW_pool_sediment')
   call self%register_state_dependency(self%id_NDetpoolS, 'detritus_N_pool_sediment', 'g m-2', 'detritus_N_pool_sediment')
   call self%register_state_dependency(self%id_PDetpoolS, 'detritus_P_pool_sediment', 'g m-2', 'detritus_P_pool_sediment')
   call self%register_state_dependency(self%id_SiDetpoolS,'detritus_Si_pool_sediment','g m-2', 'detritus_Si_pool_sediment')
   call self%register_state_dependency(self%id_SiO2poolW, 'SiO2_pool_water',          'g m-3', 'SiO2_pool_water')


   call self%register_dependency(self%id_uTm,standard_variables%temperature)

      return

   end subroutine initialize
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_pclake_phytoplankton_sediment), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
! !LOCAL VARIABLES:
!  state variables value carriers
   real(rk)                   :: sDDiatS,sPDiatS,sNDiatS
   real(rk)                   :: sDGrenS,sPGrenS,sNGrenS
   real(rk)                   :: sDBlueS,sPBlueS,sNBlueS
   real(rk)                   :: oSiDiatS,rPDPhytS,rNDPhytS
!  environmental dependencies carrier
   real(rk)                   :: uTm
!  external links value carrier
   real(rk)                  :: sPO4S,sNO3S,sNH4S,sSiO2W
   real(rk)                  :: sDDetS,sNDetS,sPDetS,sSiDetS
   
!  Nutrients rations
   real(rk)   :: rPDDiatS,rNDDiatS
   real(rk)   :: rPDGrenS,rNDGrenS
   real(rk)   :: rPDBlueS,rNDBlueS
!  Temperature variables
   real(rk)   :: uFunTmDiat,uFunTmBlue,uFunTmGren
!  Diatom dry-weight variables
   real(rk)   :: tDPrimDiatS,tDMortDiatS,tDRespDiatS,ukDRespTmDiat
!  Diatom Nitrogen variables
   real(rk)   :: tNPrimDiatS,tNMortDiatS,tNExcrDiatS
!  Diatom Phosphrus variables
   real(rk)   :: tPPrimDiatS,tPMortDiatS,tPExcrDiatS
!  Green algae dry-weight variables
   real(rk)   :: tDPrimGrenS,tDMortGrenS,tDRespGrenS,ukDRespTmGren
!  Green algae Nitrogen variables
   real(rk)   :: tNPrimGrenS,tNMortGrenS,tNExcrGrenS
!  Green algae Phosphrus variables
   real(rk)   :: tPPrimGrenS,tPMortGrenS,tPExcrGrenS
!  Blue algae dry-weight variables
   real(rk)   :: tDPrimBlueS,tDMortBlueS,tDRespBlueS,ukDRespTmBlue
!  Blue algae Nitrogen variables
   real(rk)   :: tNPrimBlueS,tNMortBlueS,tNExcrBlueS
!  Blue algae Phosphrus variables
   real(rk)   :: tPPrimBlueS,tPMortBlueS,tPExcrBlueS
!  Exchange related variables
   real(rk)   :: tNPrimNH4S,tNExcrPhytS,tNMortPhytNH4S,tNMortPhytS
   real(rk)   :: tNPrimNO3S
   real(rk)   :: tPPrimPO4S,tPExcrPhytS,tPMortPhytPO4S,tPMortPhytS
   real(rk)   :: tDPrimDetS,tDMortPhytS
   real(rk)   :: tNPrimDetS,tNMortPhytDetS
   real(rk)   :: tPPrimDetS,tPMortPhytDetS
   real(rk)   :: tSiPrimDetS,tSiMortDiatS
   real(rk)   :: tSiExcrDiatS

   
!  Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_
!-----------------------------------------------------------------------
!  Retrieve current (local) state variable values.
!-----------------------------------------------------------------------
   _GET_HORIZONTAL_(self%id_sDDiatS,sDDiatS)
   _GET_HORIZONTAL_(self%id_sNDiatS,sNDiatS)
   _GET_HORIZONTAL_(self%id_sPDiatS,sPDiatS)
   _GET_HORIZONTAL_(self%id_sDGrenS,sDGrenS)
   _GET_HORIZONTAL_(self%id_sNGrenS,sNGrenS)
   _GET_HORIZONTAL_(self%id_sPGrenS,sPGrenS)
   _GET_HORIZONTAL_(self%id_sDBlueS,sDBlueS)
   _GET_HORIZONTAL_(self%id_sNBlueS,sNBlueS)
   _GET_HORIZONTAL_(self%id_sPBlueS,sPBlueS)
!-----------------------------------------------------------------------
!  Retrieve dependencis value
!-----------------------------------------------------------------------
!  Retrieve state dependencie value
   _GET_HORIZONTAL_(self%id_PO4poolS,sPO4S)
   _GET_HORIZONTAL_(self%id_NO3poolS,sNO3S)
   _GET_HORIZONTAL_(self%id_NH4poolS,sNH4S)
   _GET_HORIZONTAL_(self%id_DDetpoolS,sDDetS)
   _GET_HORIZONTAL_(self%id_NDetpoolS,sNDetS)
   _GET_HORIZONTAL_(self%id_PDetpoolS,sPDetS)
   _GET_HORIZONTAL_(self%id_SiDetpoolS,sSiDetS)
   _GET_(self%id_SiO2poolW,sSiO2W)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
!-----------------------------------------------------------------------
!  Current local nutrients ratios(check the curent state)
!-----------------------------------------------------------------------
!  Local status
!  P/D_ratio_of_Diatom
   rPDDiatS = sPDiatS /(sDDiatS+NearZero)
!  N/D_ratio_of_Diatom
   rNDDiatS = sNDiatS /(sDDiatS+NearZero)
!  P/D_ratio_of_Green_Algae
   rPDGrenS = sPGrenS /(sDGrenS+NearZero)
!  N/D_ratio_of_Green_Algae
   rNDGrenS = sNGrenS /(sDGrenS+NearZero)
!  P/D_ratio_of_Blue_Algae
   rPDBlueS = sPBlueS /(sDBlueS+NearZero)
!  N/D_ratio_of_Blue_Algae
   rNDBlueS = sNBlueS /(sDBlueS+NearZero)
!  Settled diatom s in silica element
   oSiDiatS= self%cSiDDiat*sDDiatS
!  P/D_ratio_of_settled_phytoplankton
   rPDPhytS= (sPDiatS+sPGrenS+sPBlueS)/(sDDiatS+sDGrenS+sDBlueS+NearZero)
   rNDPhytS= (sNDiatS+sNGrenS+sNBlueS)/(sDDiatS+sDGrenS+sDBlueS+NearZero)
!-----------------------------------------------------------------------
!  Temperature functions for pelagic phytoplankton
!-----------------------------------------------------------------------
!  temperature_function_of_Diatom
   uFunTmDiat = uFunTmBio(uTm,self%cSigTmDiat,self%cTmOptDiat)
!  temperature_function_of_Blue_algae
   uFunTmBlue = uFunTmBio(uTm,self%cSigTmBlue,self%cTmOptBlue)
!  temperature_function_of_Green_Algae
   uFunTmGren = uFunTmBio(uTm,self%cSigTmGren,self%cTmOptGren)
!-----------------------------------------------------------------------
!  Dry-weight change of Diatoms in sediment
!-----------------------------------------------------------------------
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmDiat = self%kDRespDiat * uFunTmDiat
!  respiration_of_sediment_Algae
   tDRespDiatS = ukDRespTmDiat * sDDiatS
!  mortality_in_sed.
   tDMortDiatS = self%kMortDiatS * sDDiatS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tDPrimDiatS = - tDMortDiatS - tDRespDiatS
!-----------------------------------------------------------------------
!  Nitrogen change of Diatoms in sediment
!-----------------------------------------------------------------------
!  N_excretion_of_algae_in_sediment
!   tNExcrDiatS = rNDDiatS /(self%cNDDiatMin + rNDDiatS) * rNDDiatS * tDRespDiatS ! V5.09
   tNExcrDiatS = (2.0_rk * rNDDiatS) /(self%cNDDiatMax + rNDDiatS) * rNDDiatS * tDRespDiatS  ! pl613
!  N_mortality_of_algae_in_sediment
   tNMortDiatS = self%kMortDiatS * sNDiatS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tNPrimDiatS = - tNMortDiatS - tNExcrDiatS
!-----------------------------------------------------------------------
!  Phosphrus change of Diatoms in sediment
!-----------------------------------------------------------------------
!  P_excretion_of_algae_in_sediment
!   tPExcrDiatS = rPDDiatS /(self%cPDDiatMin + rPDDiatS) * rPDDiatS * tDRespDiatS ! V5.09
   tPExcrDiatS = (2.0_rk * rPDDiatS) /(self%cPDDiatMax + rPDDiatS) * rPDDiatS * tDRespDiatS  !pl613
!  P_mortality_of_algae_in_sediment
   tPMortDiatS = self%kMortDiatS * sPDiatS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tPPrimDiatS = -tPMortDiatS - tPExcrDiatS
!-----------------------------------------------------------------------
!  Dry-weight change of green algae in sediment
!-----------------------------------------------------------------------
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmGren = self%kDRespGren * uFunTmGren
!  respiration_of_sediment_Algae
   tDRespGrenS = ukDRespTmGren * sDGrenS
!  mortality_in_sed.
   tDMortGrenS = self%kMortGrenS * sDGrenS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tDPrimGrenS = - tDMortGrenS - tDRespGrenS
!-----------------------------------------------------------------------
!  Nitrogen change of green algae in sediment
!-----------------------------------------------------------------------
!  N_excretion_of_algae_in_sediment
!   tNExcrGrenS = rNDGrenS /(self%cNDGrenMin + rNDGrenS) * rNDGrenS * tDRespGrenS  ! V5.09
   tNExcrGrenS = (2.0_rk * rNDGrenS) /(self%cNDGrenMax + rNDGrenS) * rNDGrenS * tDRespGrenS  !pl613
!  N_mortality_of_algae_in_sediment
   tNMortGrenS = self%kMortGrenS * sNGrenS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tNPrimGrenS = - tNMortGrenS - tNExcrGrenS
!-----------------------------------------------------------------------
!  Phosphrus change of green algae in sediment
!-----------------------------------------------------------------------
!  P_excretion_of_algae_in_sediment
!   tPExcrGrenS = rPDGrenS /(self%cPDGrenMin + rPDGrenS) * rPDGrenS * tDRespGrenS  ! v5.09
   tPExcrGrenS = (2.0_rk * rPDGrenS) /(self%cPDGrenMax + rPDGrenS) * rPDGrenS * tDRespGrenS  !pl613
!  P_mortality_of_algae_in_sediment
   tPMortGrenS = self%kMortGrenS * sPGrenS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tPPrimGrenS = - tPMortGrenS - tPExcrGrenS 
!-----------------------------------------------------------------------
!  Dry-weight change of blue algae in sediment
!-----------------------------------------------------------------------
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmBlue = self%kDRespBlue * uFunTmBlue
!  respiration_of_sediment_Algae
   tDRespBlueS = ukDRespTmBlue * sDBlueS
!  mortality_in_sed.
   tDMortBlueS = self%kMortBlueS * sDBlueS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tDPrimBlueS = - tDMortBlueS - tDRespBlueS
!-----------------------------------------------------------------------
!  Nitrogen change of blue algae in sediment
!-----------------------------------------------------------------------
!  N_excretion_of_algae_in_sediment
!   tNExcrBlueS = rNDBlueS /(self%cNDBlueMin + rNDBlueS) * rNDBlueS * tDRespBlueS !v5.09
   tNExcrBlueS = (2.0_rk * rNDBlueS) /(self%cNDBlueMax + rNDBlueS) * rNDBlueS * tDRespBlueS  !pl613
!  N_mortality_of_algae_in_sediment
   tNMortBlueS = self%kMortBlueS * sNBlueS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tNPrimBlueS = - tNMortBlueS - tNExcrBlueS
!-----------------------------------------------------------------------
!  Phosphrus change of blue algae in sediment
!-----------------------------------------------------------------------
!  P_excretion_of_algae_in_sediment
!   tPExcrBlueS = rPDBlueS /(self%cPDBlueMin + rPDBlueS) * rPDBlueS * tDRespBlueS ! v5.09
   tPExcrBlueS = (rPDBlueS *2.0_rk)/(self%cPDBlueMax + rPDBlueS) * rPDBlueS * tDRespBlueS !pl613
!  P_mortality_of_algae_in_sediment
   tPMortBlueS = self%kMortBlueS * sPBlueS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tPPrimBlueS = - tPMortBlueS - tPExcrBlueS
!-----------------------------------------------------------------------
!  exchange to other moudules
!-----------------------------------------------------------------------
!  sNH4S exhange
!-----------------------------------------------------------------------
!  total_phytoplankton_mortality
   tNMortPhytS = tNMortDiatS + tNMortGrenS + tNMortBlueS
!  ammonium_flux_from_died_Algae
   tNMortPhytNH4S = self%fDissMortPhyt * tNMortPhytS
!  total_N_excretion_sediment_phytoplankton
   tNExcrPhytS = tNExcrDiatS + tNExcrGrenS + tNExcrBlueS
!  Pore_water_ammonium
   tNPrimNH4S = tNExcrPhytS + tNMortPhytNH4S
!-----------------------------------------------------------------------
!  sNO3S exhange, the phytoplankton in the sediment are considered
!  without uptaking process
!-----------------------------------------------------------------------
!  Pore_water_nitrate
   tNPrimNO3S = 0.0_rk
!-----------------------------------------------------------------------
!  sPO4S exhange
!-----------------------------------------------------------------------
!  total_phytoplankton_mortality
   tPMortPhytS = tPMortDiatS + tPMortGrenS + tPMortBlueS
!  soluble_P_flux_from_died_Algae
   tPMortPhytPO4S = self%fDissMortPhyt * tPMortPhytS
!   total_P_excretion_sediment_phytoplankton
    tPExcrPhytS = tPExcrDiatS + tPExcrGrenS + tPExcrBlueS
!   Pore_water_P
    tPPrimPO4S = tPExcrPhytS + tPMortPhytPO4S
!-----------------------------------------------------------------------
!  sDDetS exhange
!-----------------------------------------------------------------------
!  mortality_of_algae_on_bottom
   tDMortPhytS = tDMortDiatS + tDMortGrenS + tDMortBlueS
!  Flux_to_sediment_detritus
   tDPrimDetS = tDMortPhytS
!-----------------------------------------------------------------------
!  sNDetS exhange
!-----------------------------------------------------------------------
!  detrital_N_flux_from_died_Algae
   tNMortPhytDetS = tNMortPhytS - tNMortPhytNH4S
!  Sediment_detritus
   tNPrimDetS = tNMortPhytDetS 
!-----------------------------------------------------------------------
!  sPDetS exhange
!-----------------------------------------------------------------------
!  detrital_P_flux_from_died_Algae
   tPMortPhytDetS = tPMortPhytS - tPMortPhytPO4S
!  Sediment_detritus
   tPPrimDetS = tPMortPhytDetS 
!-----------------------------------------------------------------------
!  sSiDetS exhange
!-----------------------------------------------------------------------
!  mortality_of_bottom_Algae
   tSiMortDiatS = self%cSiDDiat * tDMortDiatS
!  Sediment_detritus
   tSiPrimDetS = tSiMortDiatS
!-----------------------------------------------------------------------
!  sSiO2W exhange
!-----------------------------------------------------------------------
!  Si_excretion_of_bottom_Algae
   tSiExcrDiatS = self%cSiDDiat * tDRespDiatS
!---------------------------------------------------------------------------
!  Update local state variables
!---------------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_sDDiatS,tDPrimDiatS)
   _SET_ODE_BEN_(self%id_sNDiatS,tNPrimDiatS)
   _SET_ODE_BEN_(self%id_sPDiatS,tPPrimDiatS)
   _SET_ODE_BEN_(self%id_sDGrenS,tDPrimGrenS)
   _SET_ODE_BEN_(self%id_sNGrenS,tNPrimGrenS)
   _SET_ODE_BEN_(self%id_sPGrenS,tPPrimGrenS)
   _SET_ODE_BEN_(self%id_sDBlueS,tDPrimBlueS)
   _SET_ODE_BEN_(self%id_sNBlueS,tNPrimBlueS)
   _SET_ODE_BEN_(self%id_sPBlueS,tPPrimBlueS)
!-----------------------------------------------------------------------
!  Output local diagnostic variables
!-----------------------------------------------------------------------
!  currently none
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_PO4poolS,tPPrimPO4S)
   _SET_ODE_BEN_(self%id_NO3poolS,tNPrimNO3S)
   _SET_ODE_BEN_(self%id_NH4poolS,tNPrimNH4S)
   _SET_ODE_BEN_(self%id_DDetpoolS,tDPrimDetS)
   _SET_ODE_BEN_(self%id_NDetpoolS,tNPrimDetS)
   _SET_ODE_BEN_(self%id_PDetpoolS,tPPrimDetS)
   _SET_ODE_BEN_(self%id_SiDetpoolS,tSiPrimDetS)
   _SET_BOTTOM_EXCHANGE_(self%id_SiO2poolW,tSiExcrDiatS)
!-----------------------------------------------------------------------
!  Output denpendent diagnostic variables for other modules
!-----------------------------------------------------------------------
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDPrimDetS,tDPrimDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_oSiDiatS,oSiDiatS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rPDPhytS,rPDPhytS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rNDPhytS,rNDPhytS)

   _FABM_HORIZONTAL_LOOP_END_


   end subroutine do_bottom
!-----------------------------------------------------------------------
! Spatial loop end
!-----------------------------------------------------------------------
!
!EOP
!-----------------------------------------------------------------------

   end module pclake_phytoplankton_sediment

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

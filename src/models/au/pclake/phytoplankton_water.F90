#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
   module au_pclake_phytoplankton_water
!
! !DESCRIPTION:
!  The pclake_abiotic_water module describes the processes related to phytoplankton in water column.
!  Three groups of phytoplankton are described here: Diatom, green algae and
!  cyanobacteria(blue algae). Each group is described in three elements, dry-weight
!  nitrogen and phosphorus. Silica contration in diatom is not a state variables here
!  but a diagnostic instead, since model assumes diatom have fixed Si/D ration,i.e. 0.1
!  The state variables and their involved processes are:
!  State variables:           sDDiatW, sNDiatW, sPDiatW(diatom concentration in DW, N and P element)
!                             sDGrenW, sNGrenW, sPGrenW(green algae concentration in DW, N and P element)
!                             sDBlueW, sNBlueW, sPBlueW(cyanobacteria concentration in DW, N and P element)
!  units( for all the groups):gDW m-3, gN m-3, gP m-3 in three elements respectively,
!  involved processes( for all the groups):
!  assimilation(primary production,only for sDDiatW,sDGrenW,sDBlueW),nutrient uptake(only for
!  sNDiatW,sNGrenW,sNBlueW,sPDiatW,sPGrenW,sPBlueW), respiration(only for sDDiatW,sDGrenW,sDBlueW)
!  excretion(only for sNDiatW,sNGrenW,sNBlueW,sPDiatW,sPGrenW,sPBlueW), mortality(for all state
!  variables)
!  This module also describes the processes which influence the state variables registered in
!  other modules, including: (aPhytW stands for all groups of phytoplankton)
!  nutrients uptaken and excreted by phytoplankton in water column sNH4W<==>aNPhytW,sPO4<==>aPPhytW
!                           sNO3==>aPhytW(only take up),SiO2W<==sSiDiatW(only excretion)
!  detritus morted by phytoplankton: sDDetW<==aDPhytW,sNDetW<==aNPhytW,sPDetW<==aPPhytW,sSiDetW<==sSiDiatW
!  oxygen produced by primary production and consumed by respiration: sO2W<==>aDPhytW
!
! !USES:
   use fabm_types
   use fabm_expressions
   use fabm_standard_variables
   use au_pclake_utility, ONLY:uFunTmBio



   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_au_pclake_phytoplankton_water
!     local state variable identifiers
!     id_sDDiatW,id_sDGrenW,id_sDBlueW: phytoplankton concentration in dry-weight, gDW m-**3
!     id_sPDiatW,id_sPGrenW,id_sPBlueW: phytoplankton concentration in nitrogen element, gN m-**3
!     id_sNDiatW,id_sNGrenW,id_sNBlueW: phytoplankton concentration in phosphorus element, gP m-**3
      type (type_state_variable_id)            :: id_sDDiatW,id_sPDiatW,id_sNDiatW
      type (type_state_variable_id)            :: id_sDGrenW,id_sPGrenW,id_sNGrenW
      type (type_state_variable_id)            :: id_sDBlueW,id_sPBlueW,id_sNBlueW
!     diagnostic variables for local output
!     id_dPAR, photosynthetic active radiation
!     id_oSiDiatW, diatom concentration in silica element, gSi m-**3
!     id_wO2PrimW, primary production of O2
!     id_oChlaDiat,diatom concentration in chlorophyll a unit,mg m-3
!     id_oChlaGren,green algae concentration in chlorophyll a unit,mg m-3
!     id_oChlaBlue,cyanobacteria concentration in chlorophyll a unit,mg m-3
!     id_rPDDiatW,id_rPDGrenW,id_rPDBlueW,id_rPDPhytW, P/D ratio of phytoplankton
!     id_rNDDiatW,id_rNDGrenW,id_rNDBlueW,id_rNDPhytW, N/D rartio of phytoplankton
      type (type_diagnostic_variable_id)       :: id_oSiDiatW
      type (type_diagnostic_variable_id)       :: id_oChlaDiat,id_oChlaGren,id_oChlaBlue
      type (type_diagnostic_variable_id)       :: id_rPDDiatW,id_rPDGrenW,id_rPDBlueW,id_rPDPhytW
      type (type_diagnostic_variable_id)       :: id_rNDDiatW,id_rNDGrenW,id_rNDBlueW,id_rNDPhytW
!     diagnostic variables relating light limitation function
      type (type_diagnostic_variable_id)       :: id_extDiat,id_extGren,id_extBlue
      type (type_diagnostic_variable_id)       :: id_phypar,id_phytoextinction
      type (type_diagnostic_variable_id)       :: id_aLLimDiat,id_aLLimGren,id_aLLimBlue
      type (type_diagnostic_variable_id)       :: id_aNutLimDiat,id_aNutLimBlue,id_aNutLimGren
!     diagnostic variables for modular fluxes
      type (type_diagnostic_variable_id)       :: id_wSiPrimSiO2W,id_wPPrimPO4W,id_wNPrimNO3W
      type (type_diagnostic_variable_id)       :: id_wNPrimNH4W,id_wO2PrimW,id_wDPrimDetW
      type (type_diagnostic_variable_id)       :: id_wNPrimDetW,id_wPPrimDetW,id_wSiPrimDetW
      type (type_diagnostic_variable_id)       :: id_wDPrimDiatW,id_wPPrimDiatW,id_wNPrimDiatW
      type (type_diagnostic_variable_id)       :: id_wDPrimGrenW,id_wPPrimGrenW,id_wNPrimGrenW
      type (type_diagnostic_variable_id)       :: id_wDPrimBlueW,id_wPPrimBlueW,id_wNPrimBlueW
!     state dependencies identifiers
      type (type_state_variable_id)            :: id_SiO2poolW,id_PO4poolW,id_NO3poolW,id_NH4poolW  ! dissolved nutrient for uptaking
      type (type_state_variable_id)            :: id_O2poolW,id_DDetpoolW,id_NDetpoolW,id_PDetpoolW,id_SiDetpoolW
!     environmental dependencies
      type (type_global_dependency_id)         :: id_Day
      type (type_dependency_id)                :: id_uTm,id_dz
      type (type_dependency_id)                :: id_par,id_meanpar,id_extc
!     diagnostic dependency
      type (type_horizontal_dependency_id)     :: id_afCovSurfVeg
!     Model parameters
!     temperature parameters
      real(rk)   :: cSigTmDiat,cTmOptDiat
      real(rk)   :: cSigTmBlue,cTmOptBlue,cSigTmGren,cTmOptGren,cSigTmLoss
!     nutrient parameters
      real(rk)   :: cPDDiatMin,cPDDiatMax,cNDDiatMin,cNDDiatMax,hSiAssDiat  ! Diatoms
      real(rk)   :: cPDGrenMin,cPDGrenMax,cNDGrenMin,cNDGrenMax,hSiAssGren  ! Green algae
      real(rk)   :: cPDBlueMin,cPDBlueMax,cNDBlueMin,cNDBlueMax,hSiAssBlue  ! Blue algae
!     light function parameters
      integer    :: UseLightMethodGren,UseLightMethodBlue,UseLightMethodDiat
      real(rk)   :: hLRefGren,hLRefDiat,hLRefBlue
      real(rk)   :: cLOptRefGren,cLOptRefBlue,cLOptRefDiat
!     growth rate parameter
      real(rk)   :: cMuMaxDiat,cMuMaxGren,cMuMaxBlue
!     respiration parameters
      real(rk)   :: kDRespDiat,kDRespGren,kDRespBlue

!     mortality parameters
      real(rk)   :: kMortDiatW,kMortGrenW,kMortBlueW
!     uptaking parameters
      real(rk)   :: cVPUptMaxDiat,cVPUptMaxGren,cVPUptMaxBlue
      real(rk)   :: cAffPUptDiat,cAffPUptGren,cAffPUptBlue
      real(rk)   :: cVNUptMaxDiat,cVNUptMaxGren,cVNUptMaxBlue
      real(rk)   :: cAffNUptDiat,cAffNUptGren,cAffNUptBlue
!  nutrient exchange parameters
      real(rk)   :: fDissMortPhyt
      real(rk)   :: cCPerDW,hO2BOD
!     Si/D ratio in diatom
      real(rk)   :: cSiDDiat
!     chla related variables
      real(rk)   :: cChDDiatMin,cChDDiatMax,cChDGrenMin,cChDGrenMax,cChDBlueMin,cChDBlueMax
!     sinking parameters
      real(rk)   :: cVSetDiat,cVSetGren,cVSetBlue
!     parameter for specific light attenuation coefficient
      real(rk)   :: cExtSpDiat,cExtSpGren,cExtSpBlue
!     minimum state variable values
      real(rk)   :: cDBlueMinW,cDGrenMinW,cDDiatMinW

   contains

!     Model procedure
      procedure :: initialize
      procedure :: do
      procedure :: get_light_extinction
   end type type_au_pclake_phytoplankton_water

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: Pi=3.14159265358979_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
!  ratio of mol.weights, = 32/12 [gO2/gC],
   real(rk),parameter :: molO2molC = 2.6667_rk
!  mol_O2_formed_per_mol_NO3-_ammonified
   real(rk),parameter ::O2PerNO3 = 1.5_rk
!  ratio of mol.weights,32/14 [gO2/gN],
   real(rk),parameter :: molO2molN = 2.2857_rk
   real(rk),parameter :: mgPerg = 1000_rk
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
   class (type_au_pclake_phytoplankton_water), intent(inout), target :: self
   integer,                     intent(in)            :: configunit

!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cSigTmDiat,        'cSigTmDiat',        '�C',        'temperature constant diatoms (sigma in Gaussian curve)',                                                                   default=20.0_rk)
   call self%get_parameter(self%cTmOptDiat,        'cTmOptDiat',        '�C',        'optimum temp. diatoms',                                                                                                    default=18.0_rk)
   call self%get_parameter(self%cSigTmBlue,        'cSigTmBlue',        '�C',        'temperature constant blue-greens (sigma in Gaussian curve)',                                                               default=12.0_rk)
   call self%get_parameter(self%cTmOptBlue,        'cTmOptBlue',        '�C',        'optimum temp. blue-greens',                                                                                                default=25.0_rk)
   call self%get_parameter(self%cSigTmGren,        'cSigTmGren',        '�C',        'temperature constant greens (sigma in Gaussian curve)',                                                                    default=15.0_rk)
   call self%get_parameter(self%cTmOptGren,        'cTmOptGren',        '�C',        'optimum temp. of greens',                                                                                                  default=25.0_rk)
   call self%get_parameter(self%cPDDiatMin,        'cPDDiatMin',        'mgP/mgDW',  'minimum P/day ratio Diatoms',                                                                                              default=0.0005_rk)
   call self%get_parameter(self%cPDDiatMax,        'cPDDiatMax',        'mgP/mgDW',  'max. P/day ratio Diatoms',                                                                                                 default=0.005_rk)
   call self%get_parameter(self%cNDDiatMin,        'cNDDiatMin',        'mgN/mgDW',  'minimum N/day ratio Diatoms',                                                                                              default=0.01_rk)
   call self%get_parameter(self%cNDDiatMax,        'cNDDiatMax',        'mgN/mgDW',  'max. N/day ratio Diatoms',                                                                                                 default=0.05_rk)
   call self%get_parameter(self%hSiAssDiat,        'hSiAssDiat',        'mgSi/l',    'half-sat. Si for diatoms',                                                                                                 default=0.09_rk)
   call self%get_parameter(self%cPDGrenMin,        'cPDGrenMin',        'mgP/mgDW',  'minimum P/day ratio greens',                                                                                               default=0.0015_rk)
   call self%get_parameter(self%cPDGrenMax,        'cPDGrenMax',        'mgP/mgDW',  'max. P/day ratio greens',                                                                                                  default=0.015_rk)
   call self%get_parameter(self%cNDGrenMin,        'cNDGrenMin',        'mgN/mgDW',  'minimum N/day ratio greens',                                                                                               default=0.02_rk)
   call self%get_parameter(self%cNDGrenMax,        'cNDGrenMax',        'mgN/mgDW',  'max. N/day ratio greens',                                                                                                  default=0.1_rk)
   call self%get_parameter(self%hSiAssGren,        'hSiAssGren',        'mgSi/l',    'half-sat. Si conc. for growth of green algae',                                                                             default=0.0_rk)
   call self%get_parameter(self%cPDBlueMin,        'cPDBlueMin',        'mgP/mgDW',  'minimum P/day ratio Bluegreens',                                                                                           default=0.0025_rk)
   call self%get_parameter(self%cPDBlueMax,        'cPDBlueMax',        'mgP/mgDW',  'max. P/day ratio blue-greens',                                                                                             default=0.025_rk)
   call self%get_parameter(self%cNDBlueMin,        'cNDBlueMin',        'mgN/mgDW',  'minimum N/day ratio Bluegreens',                                                                                           default=0.03_rk)
   call self%get_parameter(self%cNDBlueMax,        'cNDBlueMax',        'mgN/mgDW',  'max. N/day ratio Bluegreens',                                                                                              default=0.15_rk)
   call self%get_parameter(self%hSiAssBlue,        'hSiAssBlue',        'mgSi/l',    'half-sat. Si conc. for growth of blue-greens',                                                                             default=0.0_rk)
   call self%get_parameter(self%hLRefGren,         'hLRefGren',         'W m-2',     'half-sat. PAR for green algae at 20 �C (Lehmann function)',                                                                default=17.0_rk)
   call self%get_parameter(self%cLOptRefBlue,      'cLOptRefBlue',      'W m-2',     'optimum PAR for blue-greens at 20 �C (Steele function)',                                                                   default=13.6_rk)
   call self%get_parameter(self%cLOptRefDiat,      'cLOptRefDiat',      'W m-2',     'optimum PAR for Diatoms at 20 �C (Steele function)',                                                                       default=54.0_rk)
   call self%get_parameter(self%cMuMaxBlue,        'cMuMaxBlue',        'd-1',       'maximum_growth_rate_Bluegreens',                                                                                           default=0.6_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cMuMaxGren,        'cMuMaxGren',        'd-1',       'maximum_growth_rate_greens',                                                                                               default=1.5_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cMuMaxDiat,        'cMuMaxDiat',        'd-1',       'maximum_growth_rate_Diatoms',                                                                                              default=2.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespDiat,        'kDRespDiat',        'd-1',       'maintenance respiration constant diatoms ',                                                                                default=0.1_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespGren,        'kDRespGren',        'd-1',       'maintenance respiration constant greens ',                                                                                 default=0.075_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespBlue,        'kDRespBlue',        'd-1',       'maintenance respiration constant blue-greens ',                                                                            default=0.03_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortDiatW,        'kMortDiatW',        'd-1',       'mortality constant of Diatoms in water',                                                                                   default=0.01_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortGrenW,        'kMortGrenW',        'd-1',       'mortality constant of Diatoms in water',                                                                                   default=0.01_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortBlueW,        'kMortBlueW',        'd-1',       'mortality constant of blue-greens in water',                                                                               default=0.01_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVPUptMaxDiat,     'cVPUptMaxDiat',     'mgP/mgDW/d',' maximum P uptake capacity of Diatoms',                                                                                    default=0.01_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVPUptMaxGren,     'cVPUptMaxGren',     'mgP/mgDW/d',' maximum P uptake capacity of greens',                                                                                     default=0.01_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVPUptMaxBlue,     'cVPUptMaxBlue',     'mgP/mgDW/d',' maximum P uptake capacity of Bluegreens',                                                                                 default=0.04_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cAffPUptDiat,      'cAffPUptDiat',      'l/mgDW/d',  'initial P uptake affinity Diatoms',                                                                                        default=0.2_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cAffPUptGren,      'cAffPUptGren',      'l/mgDW/d',  'initial P uptake affinity greens',                                                                                         default=0.2_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cAffPUptBlue,      'cAffPUptBlue',      'l/mgDW/d',  'initial P uptake affinity Bluegreens',                                                                                     default=0.8_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVNUptMaxDiat,     'cVNUptMaxDiat',     'mgN/mgDW/d',' maximum N uptake capacity of Diatoms',                                                                                    default=0.07_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVNUptMaxGren,     'cVNUptMaxGren',     'mgN/mgDW/d',' maximum N uptake capacity of greens',                                                                                     default=0.07_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVNUptMaxBlue,     'cVNUptMaxBlue',     'mgN/mgDW/d',' maximum N uptake capacity of Bluegreens',                                                                                 default=0.07_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cAffNUptDiat,      'cAffNUptDiat',      'l/mgDW/d',  'initial N uptake affinity Diatoms',                                                                                        default=0.2_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cAffNUptGren,      'cAffNUptGren',      'l/mgDW/d',  'initial N uptake affinity greensupling',                                                                                   default=0.2_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cAffNUptBlue,      'cAffNUptBlue',      'l/mgDW/d',  'initial N uptake affinity Bluegreens',                                                                                     default=0.2_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortPhyt,     'fDissMortPhyt',     '[-]',       'soluble_nutrient_fraction_of_died_Algae6',                                                                                 default=0.2_rk)
   call self%get_parameter(self%cSiDDiat,          'cSiDDiat',          'mgSi/mgDW', 'Si/D ratio of diatoms',                                                                                                    default=0.15_rk)
   call self%get_parameter(self%cCPerDW,           'cCPerDW',           'gC/gDW',    'C content of organic matter',                                                                                              default=0.4_rk)
   call self%get_parameter(self%hO2BOD,            'hO2BOD',            'mgO2/l',    'half-sat. oxygen conc. for BOD',                                                                                           default=1.0_rk)
   call self%get_parameter(self%cChDDiatMin,       'cChDDiatMin',       'mgChl/mgDW',' min._chlorophyll/C_ratio_Diatoms',                                                                                        default=0.004_rk)
   call self%get_parameter(self%cChDDiatMax,       'cChDDiatMax',       'mgChl/mgDW',' max._chlorophyll/C_ratio_Diatoms',                                                                                        default=0.012_rk)
   call self%get_parameter(self%cChDGrenMin,       'cChDGrenMin',       'mgChl/mgDW',' min._chlorophyll/C_ratio_greens',                                                                                         default=0.01_rk)
   call self%get_parameter(self%cChDGrenMax,       'cChDGrenMax',       'mgChl/mgDW',' max._chlorophyll/C_ratio_greens',                                                                                         default=0.02_rk)
   call self%get_parameter(self%cChDBlueMin,       'cChDBlueMin',       'mgChl/mgDW',' min._chlorophyll/C_ratio_Bluegreens',                                                                                     default=0.005_rk)
   call self%get_parameter(self%cChDBlueMax,       'cChDBlueMax',       'mgChl/mgDW',' max._chlorophyll/C_ratio_Bluegreens',                                                                                     default=0.015_rk)
   call self%get_parameter(self%cVSetDiat,         'cVSetDiat',         'm/d',       'sedimentation velocity Diatoms',                                                                                           default=-0.5_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVSetGren,         'cVSetGren',         'm/d',       'sedimentation velocity of greens',                                                                                         default=-0.2_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVSetBlue,         'cVSetBlue',         'm/d',       'sedimentation velocity Blue-greens',                                                                                       default=-0.06_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cExtSpDiat,        'cExtSpDiat',        'm2/gDW',    'specific_extinction_Diatoms"',                                                                                             default=0.25_rk)
   call self%get_parameter(self%cExtSpGren,        'cExtSpGren',        'm2/gDW',    'specific_extinction_greens"',                                                                                              default=0.25_rk)
   call self%get_parameter(self%cExtSpBlue,        'cExtSpBlue',        'm2/gDW',    'specific_extinction_Bluegreens"',                                                                                          default=0.35_rk)
   call self%get_parameter(self%hLRefDiat,         'hLRefDiat',         'W m-2',     'half-sat. PAR for diatoms at 20 �C (Lehmann function)',                                                                    default=6.5_rk)
   call self%get_parameter(self%hLRefBlue,         'hLRefBlue',         'W m-2',     'half-sat. PAR for blue algea at 20 �C (Lehmann function)',                                                                 default=34.0_rk)
   call self%get_parameter(self%cLOptRefGren,      'cLOptRefGren',      'W m-2',     'optimum PAR for diatoms at 20 �C (Steele function)',                                                                       default=30.0_rk)
   call self%get_parameter(self%UseLightMethodGren,'UseLightMethodGren','[-]',       'lightmethod,1without photoinhibition,Chalker(1980)model,2with photoinhibition,Klepperetal.(1988)/Ebenhohetal.(1997)model.',default=1)
   call self%get_parameter(self%UseLightMethodBlue,'UseLightMethodBlue','[-]',       'lightmethod,1without photoinhibition,Chalker(1980)model,2with photoinhibition,Klepperetal.(1988)/Ebenhohetal.(1997)model.',default=2)
   call self%get_parameter(self%UseLightMethodDiat,'UseLightMethodDiat','[-]',       'lightmethod,1without photoinhibition,Chalker(1980)model,2with photoinhibition,Klepperetal.(1988)/Ebenhohetal.(1997)model.',default=2)
!  the user defined minumun value for state variables
   call self%get_parameter(self%cDBlueMinW,           'cDBlueMinW',           'gDW/m3',             'minimun blue-green algae biomass in system',                                                             default=0.00001_rk)
   call self%get_parameter(self%cDGrenMinW,           'cDGrenMinW',           'gDW/m3',             'minimun green algae biomass in system',                                                                  default=0.00001_rk)
   call self%get_parameter(self%cDDiatMinW,           'cDDiatMinW',           'gDW/m3',             'minimun diatom biomass in system',                                                                       default=0.00001_rk)
!  Register local state variable
!  all phytoplankton has vertical movement activated,normally netgative, meaning settling.
   call self%register_state_variable(self%id_sDDiatW,'sDDiatW','g m-3','diatom_D in water',    &
                                    initial_value=0.5_rk,minimum=self%cDDiatMinW,vertical_movement= self%cVSetDiat,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPDiatW,'sPDiatW','g m-3','diatom_P in water',     &
                                    initial_value=0.005_rk,minimum=self%cDDiatMinW* self%cPDDiatMin,vertical_movement= self%cVSetDiat,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNDiatW,'sNDiatW','g m-3','diatom_N in water',     &
                                    initial_value=0.05_rk,minimum=self%cDDiatMinW * self%cNDDiatMin,vertical_movement= self%cVSetDiat,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sDGrenW,'sDGrenW','g m-3','green_D in water',    &
                                    initial_value=0.5_rk,minimum=NearZero,vertical_movement= self%cVSetGren,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPGrenW,'sPGrenW','g m-3','green_P in water',     &
                                    initial_value=0.005_rk,minimum=self%cDGrenMinW * self%cPDGrenMin,vertical_movement= self%cVSetGren,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNGrenW,'sNGrenW','g m-3','green_N in water',     &
                                    initial_value=0.05_rk,minimum=self%cDGrenMinW * self%cNDGrenMin,vertical_movement= self%cVSetGren,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sDBlueW,'sDBlueW','g m-3','blue_D in water',     &
                                    initial_value=3.0_rk,minimum=self%cDBlueMinW,vertical_movement= self%cVSetBlue,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPBlueW,'sPBlueW','g m-3','blue_P in water',     &
                                    initial_value=0.03_rk,minimum=self%cDBlueMinW * self%cPDBlueMin,vertical_movement= self%cVSetBlue,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNBlueW,'sNBlueW','g m-3','blue_N in water',     &
                                    initial_value=0.3_rk,minimum=self%cDBlueMinW * self%cNDBlueMin,vertical_movement= self%cVSetBlue,no_river_dilution=.TRUE.)
!     register diagnostic variables
   call self%register_diagnostic_variable(self%id_oSiDiatW,       'oSiDiatW',       'g m-3 ',   'diatom_Si in water',                           output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oChlaDiat,      'oChlaDiat',      'mg m-3',   'diatom_chla',                                  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oChlaGren,      'oChlaGren',      'mg m-3',   'green_algae_chla',                             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_oChlaBlue,      'oChlaBlue',      'mg m-3',   'blue_algae_chla',                              output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDDiatW,       'rPDDiatW',       '[-]',      'diatom_P/D_ration_wat',                        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDGrenW,       'rPDGrenW',       '[-]',      'green_P/D_ration_wat',                         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDBlueW,       'rPDBlueW',       '[-]',      'blue_P/D_ration_wat',                          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDDiatW,       'rNDDiatW',       '[-]',      'diatom_N/D_ration_wat',                        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDGrenW,       'rNDGrenW',       '[-]',      'green_N/D_ration_wat',                         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDBlueW,       'rNDBlueW',       '[-]',      'blue_N/D_ration_wat',                          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDPhytW,       'rPDPhytW',       '[-]',      'Phyt_P/D_ration_wat',                          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDPhytW,       'rNDPhytW',       '[-]',      'Phyt_N/D_ration_wat',                          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_extDiat,        'extDiat',        '[-]',      'extinction factor caused by diatom',           output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_extGren,        'extGren',        '[-]',      'extinction factor caused by green algae',      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_extBlue,        'extBlue',        '[-]',      'extinction factor caused by blue green algae', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_phytoextinction,'phytoextinction','[-]',      'phytoextinction',                              output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_phypar,         'phypar',         '[-]',      'local PAR',                                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aLLimDiat,      'aLLimDiat',      '[-]',      'diatoms light limitation factor',              output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aLLimGren,      'aLLimGren',      '[-]',      'green algae light limitation factor',          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aLLimBlue,      'aLLimBlue',      '[-]',      'green blue algae light limitation factor',     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aNutLimDiat,    'aNutLimDiat',    '[-]',      'diatoms nutrients limitation factor',          output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aNutLimGren,    'aNutLimGren',    '[-]',      'green algae  nutrients limitation factor',     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aNutLimBlue,    'aNutLimBlue',    '[-]',      'green blue algae nutrients limitation factor', output=output_instantaneous)
!  diagnostic variables for modular fluexes
!  modular fluxes regarding other modules, mainly abiotic water module
   call self%register_diagnostic_variable(self%id_wSiPrimSiO2W,   'wSiPrimSiO2W',   'g m-3 s-1', 'phytoplankton_water_SiO2_change',             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPPrimPO4W,     'wNPrimPO4W',     'g m-3 s-1', 'phytoplankton_water_PO4W_change',             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNPrimNH4W,     'wNPrimNH4W',     'g m-3 s-1', 'phytoplankton_water_NH4W_change',             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNPrimNO3W,     'wNPrimNO3W',     'g m-3 s-1', 'phytoplankton_water_NO3W_change',             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wO2PrimW,       'wO2PrimW',       'g m-3 s-1', 'phytoplankton_water_O2_change',               output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDPrimDetW,     'wDPrimDetW',     'g m-3 s-1', 'phytoplankton_water_DDetW_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNPrimDetW,     'wNPrimDetW',     'g m-3 s-1', 'phytoplankton_water_NDetW_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPPrimDetW,     'wPPrimDetW',     'g m-3 s-1', 'phytoplankton_water_PDetW_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wSiPrimDetW,    'wSiPrimDetW',    'g m-3 s-1', 'phytoplankton_water_SiDetW_change',           output=output_instantaneous)
!  modular fluxes regarding local state variables(algal concentration change)
   call self%register_diagnostic_variable(self%id_wDPrimDiatW,    'wDPrimDiatW',    'g m-3 s-1', 'phytoplankton_water_DDiat_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNPrimDiatW,    'wNPrimDiatW',    'g m-3 s-1', 'phytoplankton_water_NDiat_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPPrimDiatW,    'wPPrimDiatW',    'g m-3 s-1', 'phytoplankton_water_PDiat_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDPrimGrenW,    'wDPrimGrenW',    'g m-3 s-1', 'phytoplankton_water_DGren_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNPrimGrenW,    'wNPrimGrenW',    'g m-3 s-1', 'phytoplankton_water_NGren_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPPrimGrenW,    'wPPrimGrenW',    'g m-3 s-1', 'phytoplankton_water_PGren_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDPrimBlueW,    'wDPrimBlueW',    'g m-3 s-1', 'phytoplankton_water_DBlue_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNPrimBlueW,    'wNPrimBlueW',    'g m-3 s-1', 'phytoplankton_water_NBlue_change',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPPrimBlueW,    'wPPrimBlueW',    'g m-3 s-1', 'phytoplankton_water_PBlue_change',            output=output_instantaneous)


!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNDiatW)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNGrenW)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNBlueW)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPDiatW)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPGrenW)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPBlueW)
!   call self%add_to_aggregate_variable(standard_variables%total_silicate,              self%id_oSiDiatW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNDiatW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNGrenW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNBlueW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPDiatW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPGrenW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPBlueW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_chla'),self%id_oChlaDiat)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_chla'),self%id_oChlaGren)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_chla'),self%id_oChlaBlue)
!  register state variables dependencies
   call self%register_state_dependency(self%id_SiO2poolW, 'SiO2_pool_water',            'g m-3', 'SiO2_pool_water')
   call self%register_state_dependency(self%id_PO4poolW,  'PO4_pool_water',             'g m-3', 'PO4_pool_water')
   call self%register_state_dependency(self%id_NH4poolW,  'NH4_pool_water',             'g m-3', 'NH4_pool_water')
   call self%register_state_dependency(self%id_NO3poolW,  'NO3_pool_water',             'g m-3', 'NO3_pool_water')
   call self%register_state_dependency(self%id_O2poolW,   'Oxygen_pool_water',          'g m-3', 'Oxygen_pool_water)')
   call self%register_state_dependency(self%id_DDetpoolW, 'Detritus_DW_pool_water',     'g m-3', 'Detritus_DW_pool_water')
   call self%register_state_dependency(self%id_NDetpoolW, 'Detritus_N_pool_water',      'g m-3', 'Detritus_N_pool_water')
   call self%register_state_dependency(self%id_PDetpoolW, 'Detritus_P_pool_water',      'g m-3', 'Detritus_P_pool_water')
   call self%register_state_dependency(self%id_SiDetpoolW,'Detritus_Si_pool_water',     'g m-3', 'Detritus_Si_pool_water')
!  register diagnostic dependencies
   call self%register_dependency(self%id_afCovSurfVeg,    'Surface_vegetation_coverage','[-]',   'Surface_vegetation_coverage')
!  register environmental dependencies
   call self%register_dependency(self%id_uTm, standard_variables%temperature)
   call self%register_dependency(self%id_dz,  standard_variables%cell_thickness)
   call self%register_dependency(self%id_extc,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_Day, standard_variables%number_of_days_since_start_of_the_year)
!  environmental light dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_meanpar,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk))

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
! !INPUT PARAMETERS:
   class (type_au_pclake_phytoplankton_water), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_
! !LOCAL VARIABLES:
!  environmental dependency variables
   real(rk)   :: uTm,dz,extc,Day,meanpar
!  diagnostic dependency
   real(rk)   :: aCorO2BOD,afCovSurfVeg
!  state variable value carrier
   real(rk)   :: sDDiatW,sPDiatW,sNDiatW
   real(rk)   :: sDGrenW,sPGrenW,sNGrenW
   real(rk)   :: sDBlueW,sPBlueW,sNBlueW
!  external links variable carrier
   real(rk)   :: sSiO2W,sPO4W,sNH4W,sNO3W
   real(rk)   :: sO2W,sDDetW,sNDetW,sPDetW
!  Nutrients rations
   real(rk)   :: rPDDiatW,rNDDiatW
   real(rk)   :: rPDGrenW,rNDGrenW
   real(rk)   :: rPDBlueW,rNDBlueW
   real(rk)   :: rPDPhytW,rNDPhytW
!  Temperature variables
   real(rk)   :: uFunTmDiat,uFunTmBlue,uFunTmGren
!  Primary production variables
!  Primary production fluxes
   real(rk)   :: wDPrimDiatW,wDPrimGrenW,wDPrimBlueW
   real(rk)   :: wPPrimDiatW,wPPrimGrenW,wPPrimBlueW
   real(rk)   :: wNPrimDiatW,wNPrimGrenW,wNPrimBlueW
!  Nutrient limitation variables
   real(rk)   :: aPLimDiat,aNLimDiat,aSiLimDiat
   real(rk)   :: aPLimGren,aNLimGren,aSiLimGren
   real(rk)   :: aPLimBlue,aNLimBlue,aSiLimBlue
   real(rk)   :: aNutLimDiat,aNutLimGren,aNutLimBlue
!  Light function variables
   real(rk)   :: aLPARBot,ufDay,uLPARSurf
   real(rk)   :: uhLGren,uOptBlue,uOptDiat
   real(rk)   :: uOptGren,uhLBlue,uhLDiat
   real(rk)   :: aLLimDiat,aLLimBlue,aLLimGren
!  growth rate at current temp_and_light
   real(rk)   :: aMuTmLDiat,aMuTmLGren,aMuTmLBlue
!  total growth rate
   real(rk)   :: aMuDiat,aMuGren,aMuBlue
!  assimilation rate
   real(rk)   :: wDAssDiat,wDAssGren,wDAssBlue,wDAssPhyt
!  Respiration variables
   real(rk)   :: wDRespDiatW,wDRespGrenW,wDRespBlueW,wDRespPhytW
   real(rk)   :: ukDRespTmDiat,ukDRespTmGren,ukDRespTmBlue
!  mortality variables
   real(rk)   :: wDMortDiatW,wDMortGrenW,wDMortBlueW,wDMortPhytW
   real(rk)   :: wPMortDiatW,wPMortGrenW,wPMortBlueW,wPMortPhytW
   real(rk)   :: wNMortDiatW,wNMortGrenW,wNMortBlueW,wNMortPhytW
!  uptaking variables
   real(rk)   :: wPUptDiat,wPUptGren,wPUptBlue
   real(rk)   :: aVPUptDiat,aVPUptGren,aVPUptBlue
   real(rk)   :: aVPUptMaxCrDiat,aVPUptMaxCrGren,aVPUptMaxCrBlue
   real(rk)   :: aVNUptMaxCrDiat,aVNUptMaxCrGren,aVNUptMaxCrBlue
   real(rk)   :: ahNUptDiat,ahNUptGren,ahNUptBlue
   real(rk)   :: aVNUptDiat,aVNUptGren,aVNUptBlue
   real(rk)   :: wNUptDiat,wNUptGren,wNUptBlue
!  excretion variables
   real(rk)   :: wPExcrDiatW,wPExcrGrenW,wPExcrBlueW
   real(rk)   :: wNExcrDiatW,wNExcrGrenW,wNExcrBlueW
!  variables related to exchanging to other modules
   real(rk)   :: wNPrimNH4W,wNPrimNO3W,wPPrimPO4W,wSiPrimSiO2W,wO2PrimW
   real(rk)   :: wDPrimDetW,wNPrimDetW,wPPrimDetW,wSiPrimDetW
   real(rk)   :: wNUptNH4Diat,wNUptNH4Gren,wNUptNH4Blue
   real(rk)   :: afNH4UptDiat,afNH4UptGren,afNH4UptBlue
   real(rk)   :: wNUptNH4Phyt,wNUptNO3Phyt,wNExcrPhytW,wNMortPhytNH4W
   real(rk)   :: wNUptNO3Diat,wNUptNO3Gren,wNUptNO3Blue
   real(rk)   :: wPUptPhyt,wPExcrPhytW,wPMortPhytPO4W
   real(rk)   :: wSiExcrDiatW, wSiUptDiat
   real(rk)   :: wO2ProdPhyt,wO2RespPhytW,wO2UptNO3Phyt
   real(rk)   :: wNMortPhytDetW,wPMortPhytDetW
   real(rk)   :: wSiMortDiatW
!  Diatom_Si variables
   real(rk)    ::  oSiDiatW
!  chla calculation variables
   real(rk)    :: rChDDiat,rChDGren,rChDBlue
   real(rk)    :: oChlaDiat,oChlaGren,oChlaBlue

!EOP
!-----------------------------------------------------------------------
!BOC
!------------------------------------------------------------------------
!  Spatial loop
   _LOOP_BEGIN_
!-----------------------------------------------------------------------
!  Retrieve current (local) state variable values.
!-----------------------------------------------------------------------
   _GET_(self%id_sDDiatW,sDDiatW)
   _GET_(self%id_sPDiatW,sPDiatW)
   _GET_(self%id_sNDiatW,sNDiatW)
   _GET_(self%id_sDGrenW,sDGrenW)
   _GET_(self%id_sPGrenW,sPGrenW)
   _GET_(self%id_sNGrenW,sNGrenW)
   _GET_(self%id_sDBlueW,sDBlueW)
   _GET_(self%id_sPBlueW,sPBlueW)
   _GET_(self%id_sNBlueW,sNBlueW)
!-----------------------------------------------------------------------
!  Retrieve dependencies  value
!-----------------------------------------------------------------------
!  Retrieve state dependencie value
   _GET_(self%id_SiO2poolW,sSiO2W)
   _GET_(self%id_PO4poolW,sPO4W)
   _GET_(self%id_NH4poolW,sNH4W)
   _GET_(self%id_NO3poolW,sNO3W)
   _GET_(self%id_O2poolW,sO2W)
   _GET_(self%id_DDetpoolW,sDDetW)
   _GET_(self%id_NDetpoolW,sNDetW)
   _GET_(self%id_PDetpoolW,sPDetW)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_(self%id_extc,extc)
   _GET_(self%id_meanpar,meanpar)
   _GET_(self%id_dz,dz)
   _GET_GLOBAL_(self%id_Day,Day)
!  retrieve diagnostic dependency
   _GET_HORIZONTAL_(self%id_afCovSurfVeg,afCovSurfVeg)
!-----------------------------------------------------------------------
!   SiDiat
!-----------------------------------------------------------------------
   oSiDiatW=self%cSiDDiat*sDDiatW
!-----------------------------------------------------------------------
!   Current local nutrients ratios(check the current state)
!-----------------------------------------------------------------------
!  Local status
!   P/D_ratio_of_Diatom
   rPDDiatW = sPDiatW /(sDDiatW+NearZero)
!   N/D_ratio_of_Diatom
   rNDDiatW = sNDiatW /(sDDiatW+NearZero)
!  P/D_ratio_of_Green_Algae
   rPDGrenW = sPGrenW /(sDGrenW+NearZero)
!  N/D_ratio_of_Green_Algae
   rNDGrenW = sNGrenW /(sDGrenW+NearZero)
!  P/D_ratio_of_Blue_Algae
   rPDBlueW = sPBlueW /(sDBlueW+NearZero)
!  N/D_ratio_of_Blue_Algae
   rNDBlueW = sNBlueW /(sDBlueW+NearZero)
!  P/D_ratio_of_all_phytos
   rPDPhytW = (sPDiatW+sPGrenW+sPBlueW) /(sDDiatW+sDGrenW+sDBlueW+NearZero)
!  N/D_ratio_of_all_phytos
   rNDPhytW = (sNDiatW+sNGrenW+sNBlueW) /(sDDiatW+sDGrenW+sDBlueW+NearZero)
!-----------------------------------------------------------------------
!  Light and temperature influence for pelagic phytoplankton
!-----------------------------------------------------------------------
!--------------------------------------------------
!  Temperature functions for pelagic phytoplankton
!--------------------------------------------------
!  temperature_function_of_Diatom
   uFunTmDiat = uFunTmBio(uTm,self%cSigTmDiat,self%cTmOptDiat)
!  temperature_function_of_Blue_algae
   uFunTmBlue =uFunTmBio(uTm,self%cSigTmBlue,self%cTmOptBlue)
!  temperature_function_of_Green_Algae
   uFunTmGren = uFunTmBio(uTm,self%cSigTmGren,self%cTmOptGren)

!-----------------------------------------------------------------------
!  Light functions for pelagic phytoplankton
!-----------------------------------------------------------------------
!  Light function in PCLake include 2 motheds: 1)with depth depended
!  PAR(botthom par and top par), originally from PCLake
!  mannual with self-shading. For the monod type,we use the
!  functions from Jef Huisman & Franz J. 1994. ,equation 7.
!  The steele function is simplified. 2) Use layer-centered PAR
!  with Monod type for non-photpinhibiton and Klepper et al. (1988)
!   / Ebenhoh et al. (1997) model with photoinhibition
!     light_at_the_surface,in FABM it's at the top of the current layer
      uLPARSurf = meanpar / EXP( -extc * dz/2.0_rk )
!     light_at_the_bottom,in FABM it's at the bottom of the current layer
      aLPARBot = meanpar * EXP( -extc * dz/2.0_rk )
!  light limitation function for each group, each group has two different types
   select case(self%UseLightMethodGren)
!  case 1, without photoinhibition,Chalker (1980) model
   case(1)
!      half-saturation light intensity at current temperature
       uhLGren=self%hLRefGren*uFunTmGren
!      light limitation function for green algae, no photo-inhibition
       aLLimGren= 1.0_rk /(extc * dz) * log((1.0_rk + uLPARSurf / uhLGren) / (1.0_rk + aLPARBot /uhLGren))
!  case 2, Klepper et al. (1988) / Ebenhoh et al. (1997) model.
   case(2)
       uOptGren=uFunTmGren*self%cLOptRefGren
       aLLimGren = exp(1.0_rk) /(extc * dz) *(exp(- aLPARBot /uOptGren) - exp(- uLPARSurf /uOptGren))
   end select

   select case(self%UseLightMethodBlue)
!  case 1, without photoinhibition,Chalker (1980) model
      case(1)
!      half-saturation light intensity at current temperature
       uhLBlue=self%hLRefGren*uFunTmGren
!      light limitation function for green algae, no photo-inhibition
       aLLimBlue= 1.0_rk /(extc * dz) * log((1.0_rk + uLPARSurf / uhLBlue) / (1.0_rk + aLPARBot /uhLBlue))
!  case 2, Klepper et al. (1988) / Ebenhoh et al. (1997) model.
      case(2)
       uOptBlue=uFunTmGren*self%cLOptRefBlue
       aLLimBlue = exp(1.0_rk) /(extc * dz) *(exp(- aLPARBot /uOptBlue) - exp(- uLPARSurf /uOptBlue))
   end select

   select case(self%UseLightMethodDiat)
!  case 1, without photoinhibition,Chalker (1980) model
      case(1)
!      half-saturation light intensity at current temperature
       uhLDiat=self%hLRefDiat*uFunTmDiat
!      light limitation function for green algae, no photo-inhibition
       aLLimDiat= 1.0_rk /(extc * dz) * log((1.0_rk + uLPARSurf / uhLDiat) / (1.0_rk + aLPARBot /uhLDiat))
!  case 2, Klepper et al. (1988) / Ebenhoh et al. (1997) model.
   case(2)
       uOptDiat=uFunTmDiat*self%cLOptRefDiat
       aLLimDiat = exp(1.0_rk) /(extc * dz) *(exp(- aLPARBot /uOptDiat) - exp(- uLPARSurf /uOptDiat))
   end select
!  day_length
   ufDay = 0.5_rk - 0.2_rk * cos(2.0_rk*Pi*Day / 365.0_rk)
!  growth_rate_at_current_light_AND_temp.
   aMuTmLDiat = ufDay *(1.0_rk - afCovSurfVeg)*aLLimDiat * uFunTmDiat*self%cMuMaxDiat
!  growth_rate_at_current_light_AND_temp.
   aMuTmLGren = ufDay *(1.0_rk - afCovSurfVeg)*aLLimGren * uFunTmGren*self%cMuMaxGren
!  growth_rate_at_current_light_AND_temp.
   aMuTmLBlue = ufDay *(1.0_rk - afCovSurfVeg)*aLLimBlue * uFunTmBlue*self%cMuMaxBlue
!-----------------------------------------------------------------------
!  Nutrient limitation functions
!-----------------------------------------------------------------------
!  Nutrient functions for diatom
!  Droop_function(P)_for_Diatom
   aPLimDiat = max(0.0_rk,(1.0_rk - self%cPDDiatMin / rPDDiatW) * self%cPDDiatMax /(self%cPDDiatMax - self%cPDDiatMin))
!  Droop_function(N)_for_Diatom
   aNLimDiat = max(0.0_rk,(1.0_rk - self%cNDDiatMin / rNDDiatW) * self%cNDDiatMax /(self%cNDDiatMax - self%cNDDiatMin))
!  silica_dependence_of_growth_rate
   aSiLimDiat = sSiO2W /(self%hSiAssDiat + sSiO2W)
!  nutrient_limitation_function_of_Diatom
!   aNutLimDiat = min(aPLimDiat,aNLimDiat,aSiLimDiat)  ! v5.09
   aNutLimDiat = min(aPLimDiat,(min(aNLimDiat,aSiLimDiat)))   ! pl613
!  Nutrient functions for Green algae
!  Droop_function(P)_for_Gren_Algae
   aPLimGren = max(0.0_rk,(1.0_rk - self%cPDGrenMin / rPDGrenW) * self%cPDGrenMax /(self%cPDGrenMax - self%cPDGrenMin))
!  Droop_function(N)_for_Algae
   aNLimGren = max(0.0_rk,(1.0_rk - self%cNDGrenMin / rNDGrenW) * self%cNDGrenMax /(self%cNDGrenMax - self%cNDGrenMin))
!  silica_dependence_of_growth_rate
   aSiLimGren = sSiO2W /(self%hSiAssGren + sSiO2W)
!  nutrient_limitation_function_of_Gren_Algae
!   aNutLimGren = min(aPLimGren,aNLimGren,aSiLimGren) ! v5.09
   aNutLimGren = min(aPLimGren,(min(aNLimGren,aSiLimGren))) ! pl613
!  Nutrient functions for Blue algae
!  Droop_function(P)_for_blue_Algae
   aPLimBlue = max(0.0_rk,(1.0_rk - self%cPDBlueMin / rPDBlueW) * self%cPDBlueMax /(self%cPDBlueMax - self%cPDBlueMin))
!  Droop_function(N)_for_Algae
   aNLimBlue = max(0.0_rk,(1.0_rk - self%cNDBlueMin / rNDBlueW) * self%cNDBlueMax /(self%cNDBlueMax - self%cNDBlueMin))
!  silica_dependence_of_growth_rate
   aSiLimBlue = sSiO2W /(self%hSiAssBlue + sSiO2W)
!  nutrient_limitation_function_of_Algae
   aNutLimBlue = min(aPLimBlue,(min(aNLimBlue,aSiLimBlue)))
!-----------------------------------------------------------------------
!  Algae growth_DW
!-----------------------------------------------------------------------
!  growth_rate_diatoms
   aMuDiat = aNutLimDiat * aMuTmLDiat
!  assimilation_Algae
   wDAssDiat = aMuDiat*sDDiatW
!  growth_rate_green_algae
   aMuGren = aNutLimGren * aMuTmLGren
!  assimilation_green-Algae
   wDAssGren = aMuGren*sDGrenW
!  growth_rate_blue_algae
   aMuBlue = aNutLimBlue * aMuTmLBlue
!  assimilation_blue_Algae
   wDAssBlue = aMuBlue*sDBlueW
!  total_algal_growth
   wDAssPhyt = wDAssDiat + wDAssGren + wDAssBlue    ! used in O2 exchange part
!-----------------------------------------------------------------------
!  Algae respiration_DW
!-----------------------------------------------------------------------
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmDiat = self%kDRespDiat * uFunTmDiat
!  respiration_of_Algae_in_water
   wDRespDiatW = ukDRespTmDiat * sDDiatW
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmGren = self%kDRespGren * uFunTmGren
!  respiration_of_Algae_in_water
   wDRespGrenW = ukDRespTmGren * sDGrenW
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmBlue = self%kDRespBlue * uFunTmBlue
!  respiration_of_Algae_in_water
   wDRespBlueW = ukDRespTmBlue * sDBlueW
!  total_algal_respiration_in_water
   wDRespPhytW = wDRespDiatW + wDRespGrenW + wDRespBlueW   ! Used in O2 exchange
!-----------------------------------------------------------------------
!  Algae mortality loss_DW
!-----------------------------------------------------------------------
!  mortality_in_water
   wDMortDiatW = self%kMortDiatW * sDDiatW
!  mortality_in_water
   wDMortGrenW = self%kMortGrenW * sDGrenW
!  mortality_in_water
   wDMortBlueW = self%kMortBlueW * sDBlueW
!-----------------------------------------------------------------------
!  Algae mortality loss_P
!-----------------------------------------------------------------------
!  mortality_Algae_in_water
   wPMortDiatW = self%kMortDiatW * sPDiatW
!  mortality_Algae_in_water
   wPMortGrenW = self%kMortGrenW * sPGrenW
!  mortality_Algae_in_water
   wPMortBlueW = self%kMortBlueW * sPBlueW
!-------------------------------------------------------------------------------------------------------------
!  Algae mortality loss_N
!-------------------------------------------------------------------------------------------------------------
!  mortality_Algae_in_water
   wNMortDiatW = self%kMortDiatW * sNDiatW
!  mortality_Algae_in_water
   wNMortGrenW = self%kMortGrenW * sNGrenW
!  mortality_Algae_in_water
   wNMortBlueW = self%kMortBlueW * sNBlueW
!  total_algal_mortality_in_water
   wDMortPhytW = wDMortDiatW + wDMortGrenW + wDMortBlueW   ! used in detritues exchange
!-------------------------------------------------------------------------------------------------------------
!  Algae uptake_P
!-------------------------------------------------------------------------------------------------------------
! maximum_P_uptake_rate_of_Algae,corrected_for_P/D_ratio
   aVPUptMaxCrDiat = max(0.0_rk,self%cVPUptMaxDiat * uFunTmDiat *(self%cPDDiatMax - rPDDiatW)&
                     &/(self%cPDDiatMax - self%cPDDiatMin))
! P_uptake_rate_of_Algae
   aVPUptDiat = sPO4W * aVPUptMaxCrDiat /(aVPUptMaxCrDiat / self%cAffPUptDiat + sPO4W)
! P_uptake_Algae
   wPUptDiat = aVPUptDiat * sDDiatW
! maximum_P_uptake_rate_of_Algae,corrected_for_P/D_ratio
   aVPUptMaxCrGren = max(0.0_rk,self%cVPUptMaxGren * uFunTmGren *(self%cPDGrenMax - rPDGrenW)&
                      & /(self%cPDGrenMax - self%cPDGrenMin))
! P_uptake_rate_of_Algae
   aVPUptGren = sPO4W * aVPUptMaxCrGren /(aVPUptMaxCrGren / self%cAffPUptGren + sPO4W)
! P_uptake_Algae
   wPUptGren = aVPUptGren * sDGrenW
! maximum_P_uptake_rate_of_Algae,corrected_for_P/D_ratio
   aVPUptMaxCrBlue = max(0.0_rk,self%cVPUptMaxBlue * uFunTmBlue *(self%cPDBlueMax - rPDBlueW)&
                     &/(self%cPDBlueMax - self%cPDBlueMin))
! P_uptake_rate_of_Algae
   aVPUptBlue = sPO4W * aVPUptMaxCrBlue /(aVPUptMaxCrBlue / self%cAffPUptBlue + sPO4W)
! P_uptake_Algae
   wPUptBlue = aVPUptBlue * sDBlueW
!-----------------------------------------------------------------------
!  Algae uptake_N
!-----------------------------------------------------------------------
!  maximum_N_uptake_rate_of_Algae,corrected_for_N/D_ratio
   aVNUptMaxCrDiat = max(0.0_rk,self%cVNUptMaxDiat * uFunTmDiat * (self%cNDDiatMax - rNDDiatW) &
                    & /(self%cNDDiatMax - self%cNDDiatMin))
!  half-sat._NDissW_for_uptake_by_Algae
   ahNUptDiat = aVNUptMaxCrDiat / self%cAffNUptDiat
!  N_uptake_rate_of_Algae
   aVNUptDiat = (sNO3W + sNH4W) * aVNUptMaxCrDiat /(ahNUptDiat + sNO3W + sNH4W)
!  N_uptake_Algae
   wNUptDiat = aVNUptDiat * sDDiatW
!  maximum_N_uptake_rate_of_Algae,corrected_for_N/D_ratio
   aVNUptMaxCrGren = max(0.0_rk,self%cVNUptMaxGren * uFunTmGren * (self%cNDGrenMax - rNDGrenW)&
                    & /(self%cNDGrenMax - self%cNDGrenMin))
!  half-sat._NDissW_for_uptake_by_Algae
   ahNUptGren = aVNUptMaxCrGren / self%cAffNUptGren
!  N_uptake_rate_of_Algae
   aVNUptGren = (sNO3W+sNH4W) * aVNUptMaxCrGren /(ahNUptGren + sNO3W + sNH4W)
!  N_uptake_Algae
   wNUptGren = aVNUptGren * sDGrenW
!  maximum_N_uptake_rate_of_Algae,corrected_for_N/D_ratio
   aVNUptMaxCrBlue = max(0.0_rk,self%cVNUptMaxBlue * uFunTmBlue * (self%cNDBlueMax - rNDBlueW)&
                     & /(self%cNDBlueMax - self%cNDBlueMin))
!  half-sat._NDissW_for_uptake_by_Algae
   ahNUptBlue = aVNUptMaxCrBlue / self%cAffNUptBlue
!  N_uptake_rate_of_Algae
   aVNUptBlue = (sNO3W+sNH4W) * aVNUptMaxCrBlue /(ahNUptBlue + sNO3W+sNH4W)
!  N_uptake_Algae
   wNUptBlue = aVNUptBlue * sDBlueW
!-----------------------------------------------------------------------
!  Algae excretion_P
!-----------------------------------------------------------------------
!  P_excretion_Algae_in_water
!   wPExcrDiatW = rPDDiatW /(self%cPDDiatMin + rPDDiatW) * rPDDiatW * wDRespDiatW ! v5.09
   wPExcrDiatW = (2.0_rk * rPDDiatW) /(self%cPDDiatMax + rPDDiatW) * rPDDiatW * wDRespDiatW
!  P_excretion_Algae_in_water
!   wPExcrGrenW = rPDGrenW /(self%cPDGrenMin + rPDGrenW) * rPDGrenW * wDRespGrenW ! V5.09
   wPExcrGrenW = (2.0_rk *rPDGrenW) /(self%cPDGrenMax + rPDGrenW) * rPDGrenW * wDRespGrenW
!  P_excretion_Algae_in_water
!   wPExcrBlueW = rPDBlueW /(self%cPDBlueMin + rPDBlueW) * rPDBlueW * wDRespBlueW  ! V5.09
   wPExcrBlueW = (rPDBlueW * 2.0_rk )/(self%cPDBlueMax + rPDBlueW) * rPDBlueW * wDRespBlueW ! pl613
!-----------------------------------------------------------------------
!  Algae excretion_N
!-----------------------------------------------------------------------
!  N_excretion_Algae_in_water
!   wNExcrDiatW = rNDDiatW /(self%cNDDiatMin + rNDDiatW) * rNDDiatW * wDRespDiatW ! V5.09
   wNExcrDiatW = (2.0_rk * rNDDiatW) /(self%cNDDiatMax + rNDDiatW) * rNDDiatW * wDRespDiatW ! pl613
!  N_excretion_Algae_in_water
!   wNExcrGrenW = rNDGrenW /(self%cNDGrenMin + rNDGrenW) * rNDGrenW * wDRespGrenW ! v5.09
   wNExcrGrenW = (2.0_rk * rNDGrenW) /(self%cNDGrenMax + rNDGrenW) * rNDGrenW * wDRespGrenW  ! pl613
!  N_excretion_Algae_in_water
!   wNExcrBlueW = rNDBlueW /(self%cNDBlueMin + rNDBlueW) * rNDBlueW * wDRespBlueW  ! v5.09
   wNExcrBlueW = (rNDBlueW *2.0_rk)/(self%cNDBlueMax + rNDBlueW) * rNDBlueW * wDRespBlueW   ! pl613
!-----------------------------------------------------------------------
!  total_PRIM_flux_to_algae_in_water
!-----------------------------------------------------------------------
!  total_PRIM_flux_to_algae_in_water_DW
   wDPrimDiatW = wDAssDiat - wDRespDiatW - wDMortDiatW
   wDPrimGrenW = wDAssGren - wDRespGrenW - wDMortGrenW
   wDPrimBlueW = wDAssBlue - wDRespBlueW - wDMortBlueW
!  total_PRIM_flux_to_algae_in_water_P
   wPPrimDiatW = wPUptDiat - wPExcrDiatW  - wPMortDiatW
   wPPrimGrenW = wPUptGren - wPExcrGrenW  - wPMortGrenW
   wPPrimBlueW = wPUptBlue - wPExcrBlueW  - wPMortBlueW
!  total_PRIM_flux_to_algae_in_water_N
   wNPrimDiatW = wNUptDiat- wNExcrDiatW - wNMortDiatW
   wNPrimGrenW = wNUptGren- wNExcrGrenW - wNMortGrenW
   wNPrimBlueW = wNUptBlue- wNExcrBlueW - wNMortBlueW
!-----------------------------------------------------------------------
!  Nutrient exchange with abiotic module
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  NH4 exchange with abiotic module
!-----------------------------------------------------------------------
!  fraction_ammonium_uptake_by_Diatoms
   afNH4UptDiat = sNH4W * sNO3W /((ahNUptDiat + sNH4W) *(ahNUptDiat + sNO3W)) + sNH4W &
   & * ahNUptDiat /((sNH4W + sNO3W) *(ahNUptDiat + sNO3W))
!  ammonium_uptake_by_Algae
   wNUptNH4Diat = afNH4UptDiat * wNUptDiat
!  fraction_ammonium_uptake_by_Algae
   afNH4UptGren = sNH4W * sNO3W /((ahNUptGren + sNH4W) *(ahNUptGren + sNO3W)) + sNH&
    &4W * ahNUptGren /((sNH4W + sNO3W) *(ahNUptGren + sNO3W))
!  ammonium_uptake_by_Algae
   wNUptNH4Gren = afNH4UptGren * wNUptGren
!  fraction_ammonium_uptake_by_Algae
   afNH4UptBlue = sNH4W * sNO3W /((ahNUptBlue + sNH4W) *(ahNUptBlue + sNO3W)) + sNH&
   &4W * ahNUptBlue /((sNH4W + sNO3W) *(ahNUptBlue + sNO3W))
!  ammonium_uptake_by_Algae
   wNUptNH4Blue = afNH4UptBlue * wNUptBlue
!  total_ammonium-N_uptake_phytoplankton
!   wNUptNH4Phyt = wNUptNH4Diat + wNUptNH4Gren + wNUptNH4Blue
!  total_ammonium-N_uptake_phytoplankton
   wNUptNH4Phyt = wNUptNH4Diat + wNUptNH4Gren + wNUptNH4Blue
!  total_N_excretion_phytoplankton_in_water
   wNExcrPhytW = wNExcrDiatW + wNExcrGrenW + wNExcrBlueW
!   total_N_grazing_loss
!   wNLossPhyt = wNLossDiat + wNLossGren + wNLossBlue
!   NH4-N_grazing_loss
!   wNLossPhytNH4 = self%fDissLoss * wNLossPhyt
!  total_N_mortality_phytoplankton_in_water
   wNMortPhytW = wNMortDiatW + wNMortGrenW + wNMortBlueW
!  ammonium_flux_from_died_Algae
   wNMortPhytNH4W = self%fDissMortPhyt * (wNMortDiatW + wNMortGrenW + wNMortBlueW)
!  ammonium_in_water
   wNPrimNH4W = - wNUptNH4Phyt + wNExcrPhytW + wNMortPhytNH4W   ! + wNLossPhytNH4
!-----------------------------------------------------------------------
!  NO3 exchange with abiotic module
!-----------------------------------------------------------------------
!   nitrate_uptake_by_Algae
    wNUptNO3Blue = wNUptBlue - wNUptNH4Blue
!   nitrate_uptake_by_Algae
    wNUptNO3Gren = wNUptGren - wNUptNH4Gren
!   nitrate_uptake_by_Algae
    wNUptNO3Diat = wNUptDiat - wNUptNH4Diat
!   total_nitrate-N_uptake_phytoplankton
    wNUptNO3Phyt = wNUptNO3Diat + wNUptNO3Gren + wNUptNO3Blue
!   nitrate_in_water
    wNPrimNO3W = - wNUptNO3Phyt
!-----------------------------------------------------------------------
!  SRP(PO4) exchange with abiotic module
!-----------------------------------------------------------------------
!  total_N_mortality_phytoplankton_in_water
   wPMortPhytW = wPMortDiatW + wPMortGrenW + wPMortBlueW
!  soluble_P_flux_from_died_Algae
   wPMortPhytPO4W = self%fDissMortPhyt * wPMortPhytW
!  total_grazing_loss
!   wPLossPhyt = wPLossDiat + wPLossGren + wPLossBlue
!   soluble_P_grazing_loss
!   wPLossPhytPO4 = self%fDissLoss * wPLossPhyt
!  total_P_excretion_phytoplankton_in_water
   wPExcrPhytW = wPExcrDiatW + wPExcrGrenW + wPExcrBlueW
!  total_P_uptake_phytoplankton
   wPUptPhyt = wPUptDiat + wPUptGren + wPUptBlue
!  SRP_in_water
   wPPrimPO4W = - wPUptPhyt + wPExcrPhytW + wPMortPhytPO4W !+ wPLossPhytPO4
!-----------------------------------------------------------------------
!   Si exchange with abiotic module
!-----------------------------------------------------------------------
!  Diatoms_silica_uptake
   wSiUptDiat = self%cSiDDiat * wDAssDiat
!  Si_excretion
   wSiExcrDiatW = self%cSiDDiat * wDRespDiatW
!  total_Si_flux_to_SiO2_in_PRIM_module
   wSiPrimSiO2W = wSiExcrDiatW - wSiUptDiat
!-----------------------------------------------------------------------
!  O2 exchange with abiotic module
!-----------------------------------------------------------------------
!  O2_production_by_phytoplankton
   wO2ProdPhyt = molO2molC * self%cCPerDW * wDAssPhyt
!  correction_of_O2_demand_in_water_at_low_oxygen_conc.
   aCorO2BOD = sO2W / (self%hO2BOD + sO2W)
!  O2_production_by_phytoplankton
   wO2RespPhytW = molO2molC * self%cCPerDW * wDRespPhytW * aCorO2BOD
!  O2_production_due_to_NO3_uptake_by_phytopl.
   wO2UptNO3Phyt = O2PerNO3 * molO2molN * wNUptNO3Phyt
!  O2_flux_by_water_algae
   wO2PrimW = wO2ProdPhyt - wO2RespPhytW + wO2UptNO3Phyt
!-----------------------------------------------------------------------
!  Detritus exchange with abiotic module
!-----------------------------------------------------------------------
!  Flux_to_water_detritus
   wDPrimDetW = wDMortPhytW  !+ wDLossPhyt
!   detrital_N_grazing_loss
!   wNLossPhytDet = wNLossPhyt - wNLossPhytNH4
!  detrital_N_flux_from_died_Algae
   wNMortPhytDetW = wNMortPhytW - wNMortPhytNH4W
!  Detritus_in_water
   wNPrimDetW =wNMortPhytDetW  !+ wNLossPhytDet
!   detrital_P_grazing_loss
!   wPLossPhytDet = wPLossPhyt - wPLossPhytPO4
!  detrital_P_flux_from_died_Algae
   wPMortPhytDetW = wPMortPhytW - wPMortPhytPO4W
!  Detritus_in_water
   wPPrimDetW = wPMortPhytDetW  !+wPLossPhytDet
!   diatom_grazing_loss
!   wSiLossDiat = self%cSiDDiat * wDLossDiat
!  Diatoms_mortality_in_water
   wSiMortDiatW = self%cSiDDiat * wDMortDiatW
!  total_Si_flux_to_sed._detritus_in_PRIM_module
   wSiPrimDetW = wSiMortDiatW  !+ wSiLossDiat
!-----------------------------------------------------------------------
!  chla calculation, for output purpose
!-----------------------------------------------------------------------
! chlorophyll-a/DW_ratio_Algae
   rChDDiat = self%cChDDiatMax -(self%cChDDiatMax - self%cChDDiatMin) * aLLimDiat
! chlorophyll-a_conc.
   oChlaDiat = mgPerg * rChDDiat * sDDiatW
! chlorophyll-a/DW_ratio_Algae
   rChDGren = self%cChDGrenMax -(self%cChDGrenMax - self%cChDGrenMin) * aLLimGren
! chlorophyll-a_conc.
   oChlaGren = mgPerg * rChDGren * sDGrenW
! chlorophyll-a/DW_ratio_Algae
   rChDBlue = self%cChDBlueMax -(self%cChDBlueMax - self%cChDBlueMin) * aLLimBlue
! chlorophyll-a_conc.
   oChlaBlue = mgPerg * rChDBlue * sDBlueW
!-----------------------------------------------------------------------
!  chla calculation, for output purpose
!-----------------------------------------------------------------------
! chlorophyll-a/DW_ratio_Algae
   rChDDiat = self%cChDDiatMax -(self%cChDDiatMax - self%cChDDiatMin) * aLLimDiat
! chlorophyll-a_conc.
   oChlaDiat = mgPerg * rChDDiat * sDDiatW
! chlorophyll-a/DW_ratio_Algae
   rChDGren = self%cChDGrenMax -(self%cChDGrenMax - self%cChDGrenMin) * aLLimGren
! chlorophyll-a_conc.
   oChlaGren = mgPerg * rChDGren * sDGrenW
! chlorophyll-a/DW_ratio_Algae
   rChDBlue = self%cChDBlueMax -(self%cChDBlueMax - self%cChDBlueMin) * aLLimBlue
! chlorophyll-a_conc.
   oChlaBlue = mgPerg * rChDBlue * sDBlueW
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
!  set out diatom variables
   _SET_ODE_(self%id_sDDiatW, wDPrimDiatW)
   _SET_ODE_(self%id_sPDiatW, wPPrimDiatW)
   _SET_ODE_(self%id_sNDiatW, wNPrimDiatW)
!   set out green algae variables
   _SET_ODE_(self%id_sDGrenW, wDPrimGrenW)
   _SET_ODE_(self%id_sPGrenW, wPPrimGrenW)
   _SET_ODE_(self%id_sNGrenW, wNPrimGrenW)
!set out blue algae variables
   _SET_ODE_(self%id_sDBlueW, wDPrimBlueW)
   _SET_ODE_(self%id_sPBlueW, wPPrimBlueW)
   _SET_ODE_(self%id_sNBlueW, wNPrimBlueW)
!-----------------------------------------------------------------------
!  Output local diagnostic variables
!-----------------------------------------------------------------------
   _SET_DIAGNOSTIC_(self%id_oSiDiatW, oSiDiatW)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
   _SET_ODE_(self%id_SiO2poolW,  wSiPrimSiO2W)
   _SET_ODE_(self%id_PO4poolW,   wPPrimPO4W)
   _SET_ODE_(self%id_NO3poolW,   wNPrimNO3W)
   _SET_ODE_(self%id_NH4poolW,   wNPrimNH4W)
   _SET_ODE_(self%id_O2poolW,    wO2PrimW)
   _SET_ODE_(self%id_DDetpoolW,  wDPrimDetW)
   _SET_ODE_(self%id_NDetpoolW,  wNPrimDetW)
   _SET_ODE_(self%id_PDetpoolW,  wPPrimDetW)
   _SET_ODE_(self%id_SiDetpoolW, wSiPrimDetW)
!-----------------------------------------------------------------------
!  Output dependent diagnostic variables for other modules
!-----------------------------------------------------------------------
   _SET_DIAGNOSTIC_(self%id_oChlaDiat,   oChlaDiat)
   _SET_DIAGNOSTIC_(self%id_oChlaGren,   oChlaGren)
   _SET_DIAGNOSTIC_(self%id_oChlaBlue,   oChlaBlue)
   _SET_DIAGNOSTIC_(self%id_rPDDiatW,    rPDDiatW)
   _SET_DIAGNOSTIC_(self%id_rPDGrenW,    rPDGrenW)
   _SET_DIAGNOSTIC_(self%id_rPDBlueW,    rPDBlueW)
   _SET_DIAGNOSTIC_(self%id_rNDDiatW,    rNDDiatW)
   _SET_DIAGNOSTIC_(self%id_rNDGrenW,    rNDGrenW)
   _SET_DIAGNOSTIC_(self%id_rNDBlueW,    rNDBlueW)
   _SET_DIAGNOSTIC_(self%id_rPDPhytW,    rPDPhytW)

   _SET_DIAGNOSTIC_(self%id_aLLimDiat,   aLLimDiat)
   _SET_DIAGNOSTIC_(self%id_aLLimGren,   aLLimGren)
   _SET_DIAGNOSTIC_(self%id_aLLimBlue,   aLLimBlue)
   _SET_DIAGNOSTIC_(self%id_aNutLimDiat, aNutLimDiat)
   _SET_DIAGNOSTIC_(self%id_aNutLimGren, aNutLimGren)
   _SET_DIAGNOSTIC_(self%id_aNutLimBlue, aNutLimBlue)
   _SET_DIAGNOSTIC_(self%id_phypar,      meanpar)
   _SET_DIAGNOSTIC_(self%id_phytoextinction, extc)

!  output diagnostic variables for modular fluxes
!  fluxes for other modules, mainly abiotic water module
   _SET_DIAGNOSTIC_(self%id_wSiPrimSiO2W,wSiPrimSiO2W*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPPrimPO4W,  wPPrimPO4W*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNPrimNH4W,  wNPrimNH4W*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNPrimNO3W,  wNPrimNO3W*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wO2PrimW,    wO2PrimW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDPrimDetW,  wDPrimDetW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNPrimDetW,  wNPrimDetW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPPrimDetW,  wPPrimDetW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wSiPrimDetW, wSiPrimDetW*86400.0_rk)
!  fluxes for local state variables
   _SET_DIAGNOSTIC_(self%id_wDPrimDiatW, wDPrimDiatW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNPrimDiatW, wNPrimDiatW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPPrimDiatW, wPPrimDiatW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDPrimGrenW, wDPrimGrenW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNPrimGrenW, wNPrimGrenW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPPrimGrenW, wPPrimGrenW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wDPrimBlueW, wDPrimBlueW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wNPrimBlueW, wNPrimBlueW*86400.0_rk)
   _SET_DIAGNOSTIC_(self%id_wPPrimBlueW, wPPrimBlueW*86400.0_rk)

   _LOOP_END_

!  Spatial loop end

   end subroutine do
!EOC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_au_pclake_phytoplankton_water), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk) :: sDDiatW,sDGrenW,sDBlueW
   real(rk) :: extDiat,extGren,extBlue
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_sDDiatW,sDDiatW)
   _GET_(self%id_sDGrenW,sDGrenW)
   _GET_(self%id_sDBlueW,sDBlueW)

   extDiat=self%cExtSpDiat*sDDiatW
   extGren=self%cExtSpGren*sDGrenW
   extBlue=self%cExtSpBlue*sDBlueW


   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(extDiat+extGren+extBlue)



   _SET_DIAGNOSTIC_(self%id_extDiat,extDiat)
   _SET_DIAGNOSTIC_(self%id_extGren,extGren)
   _SET_DIAGNOSTIC_(self%id_extBlue,extBlue)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC
   end module au_pclake_phytoplankton_water

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

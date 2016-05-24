#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module au_pclake_abiotic_sediment
!
! !DESCRIPTION:
!
!  The au_pclake_abiotic_sediment module describes all the state variables
!  which are related to abiotic processes in the sediment, including:
!  inorganic matter(IM), organic matters(detritus),dissolved nutirents
!  (ammonia, nitrate,phosphate) immobilized phosphrus(absorbed phosphrus).
!  Each state variable and its related local processes are:
!  Inorganic matter: sDIMS, processes: none local processes
!  Organic matter: sDDetS,sNDetS,sPDetS, processes: mineralization
!  Dissolved nutrients: sNH4S,processes:mineralisation,nitrification,
!                       diffusion to water column
!  Dissolved nutrients:sNO3W,processes: nitrification,denitrification,
!                      diffusion to water column
!  Dissolved nutrients: sPO4S,processes: mineralisation,phosphrus
!                       absorption,difusion to water column
!  Absorbed_P: sPAIMW, processes:phosphrus absorption
!  This module also discribes the processes which influence the state
!  variables registered in other modules, including:
!  Diffusion, influences ammonia,nitrate,phosphate in water column,
!  sNH4S<==>sNH4W, sNO3S<==>sNO3W,sPO4S<==>sPO4W
!  Sediment oxygen consumption, influences oxygen in water column,
!   ==>sO2W
!  Organic silica mineralization, influences silica dioxide in water
!  column, ==>sSiO2W
!  This module also provides important diagnostic variable which will be
!   used in other modules, including:
!  Sediment oxic layer fraction, afOxySed, used by module:
!                                            macrophytes module
!  Sediment detritus change, tDAbioDetS, used bymodule:auxilary
!
!  feh: April 13th, major change: added all the modular fluxes as diagostic variables.
! ! END OF DESCRIPTION

! !USES:
   use fabm_types
   use au_pclake_utility, ONLY: uFunTmAbio
   implicit none
!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_au_pclake_abiotic_sediment
!     local state variable identifers
!     sDIMS: inorganic matter concentration, in dry-weight,gDW/m**2
!     sDDetS, detritus concentration, in dry-weight,gDW/m**2
!     sNDetS, detritus concentration, in nitrogen element,gN/m**2
!     sPDetS, detritus concentration, in phosphorus element,gP/m**2
!     sNH4S,  ammonia concentration, in nitrogen element,gN/m**2
!     sNO3S, nitrate concentration, in nitrogen element,gN/m**2
!     sPO4S, phosphate concentration, in phosphorus element,gP/m**2
!     sPAIMS, absorbed phosphorus concentration, in phosphorus element,gP/m**2
      type (type_bottom_state_variable_id) :: id_sDIMS,id_sDDetS,id_sNDetS
      type (type_bottom_state_variable_id) :: id_sPDetS,id_sSiDetS
      type (type_bottom_state_variable_id) :: id_sPO4S,id_sPAIMS,id_sNH4S
      type (type_bottom_state_variable_id) :: id_sNO3S
      type (type_bottom_state_variable_id) :: id_sDHumS,id_sNHumS,id_sPHumS
!     diagnostic variables for local output
!     rPDDetS: P/D ratio of detritus
!     rNDDetS: N/D ratio of detritus
!     aPEqIMS: equilibrium absorbed phosphrus concentration
      type (type_horizontal_diagnostic_variable_id) :: id_rPDDetS,id_rNDDetS
      type (type_horizontal_diagnostic_variable_id) :: id_aPEqIMS,id_afOxySed
!     diagnostic variable for external dependencies: 
      type (type_horizontal_diagnostic_variable_id) :: id_tDAbioHumS,id_tDAbioDetS
!   diagnostic variables for modular fluxes for each module
      type (type_horizontal_diagnostic_variable_id) :: id_tDAbioIMS,id_tPAbioDetS
      type (type_horizontal_diagnostic_variable_id) :: id_tNAbioDetS,id_tSiAbioDetS,id_tNAbioNH4S
      type (type_horizontal_diagnostic_variable_id) :: id_tNAbioNO3S,id_tPAbioPO4S,id_tPAbioAIMS
      type (type_horizontal_diagnostic_variable_id) :: id_tPAbioHumS,id_tNAbioHumS
      type (type_horizontal_diagnostic_variable_id) :: id_tNdifNH4,id_tNdifNO3,id_tPdifPO4
      type (type_horizontal_diagnostic_variable_id) :: id_tSiAbioSiO2S,id_tDAbioO2S
!     detritus and humus fluxes will be named differently due to they are used by 
!     external dependencies
      type (type_horizontal_diagnostic_variable_id) :: id_tDAbioHumSflux,id_tDAbioDetSflux

!     state dependencies identifiers
!     MinSiO2Sed: Mineralization generated SiO2 from sediment
!     O2ConsumpSed: O2 consumption in sediment
!     diff+nut: diffusion fluxes of nutrients between water and sediment
      type (type_state_variable_id) :: id_MinSiO2Sed ,id_O2ConsumpSed
      type (type_state_variable_id) :: id_diffNH4,id_diffNO3,id_diffPO4
!     environmental dependencies
!     April 15th,2016, added dz, cell_thickness(feh) 
      type (type_dependency_id)                :: id_uTm,id_dz
      type (type_horizontal_dependency_id)     :: id_depth
!     Model parameters
!     Model scale parameters
      real(rk)                   :: cDepthS,cCPerDW,O2PerNH4

!     sediment properties parameters
      real(rk)                   :: bPorS,bPorCorS

!     P-sorption parameters
      real(rk)                   :: kPSorp,cRelPAdsD
      real(rk)                   :: cRelPAdsFe,fFeDIM,cRelPAdsAl,fAlDIM
      real(rk)                   :: fRedMax,cKPAdsOx,kPChemPO4,coPO4Max
!     denitrification parameters
      real(rk)                   :: NO3PerC,hNO3Denit
!     detritus related parameters
      real(rk)                   :: fRefrDetS,cThetaMinS,kDMinDetS
      real(rk)                   :: kNitrS,cThetaNitr
!     Humus related parameters
      real(rk)                   :: kDMinHum
!     diffusion parameters
      real(rk)                   :: fDepthDifS,cThetaDif,cTurbDifNut
      real(rk)                   :: kNDifNH4,kNDifNO3,kPDifPO4
      real(rk)                   :: kO2Dif,cTurbDifO2


!
      contains
!
      procedure initialize
      procedure do_bottom
   end type type_au_pclake_abiotic_sediment


!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero=0.000000000000000000000000000000001_rk
!  ratio of mol.weights,=32/12 [gO2/gC],
   real(rk),parameter :: molO2molC=2.6667_rk
!  ratio of mol.weights,32/14 [gO2/gN],
   real(rk),parameter :: molO2molN=2.2857_rk
!  ratio of mol.weights,14/12 [gN/gC],
   real(rk),parameter :: molNmolC=1.1667_rk
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the pclake_abioti_sediment module
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the au_pclake_abiotic_sediment namelist is read and the variables
!  are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_au_pclake_abiotic_sediment), intent(inout), target :: self
   integer,                              intent(in)            :: configunit


!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cDepthS,    'cDepthS',     'm',                    'sediment depth',                                           default=0.1_rk)
   call self%get_parameter(self%fRefrDetS,  'fRefrDetS',   '[-]',                  'nutrient diffusion distance as fraction of sediment depth',default=0.15_rk)
   call self%get_parameter(self%kDMinDetS,  'kDMinDetS',   'd-1',                  'decomposition constant of sediment detritus',              default=0.002_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cThetaMinS, 'cThetaMinS',  '[-]',                  'expon. temp. constant of sediment mineralization',         default=1.07_rk)
   call self%get_parameter(self%cCPerDW,    'cCPerDW',     'gC/gDW',               'C content of organic matter',                              default=0.4_rk)
   call self%get_parameter(self%O2PerNH4,   'O2PerNH4',    'mol',                  'O2 used per mol NH4+ nitrified',                           default=2.0_rk)
   call self%get_parameter(self%kNitrS,     'kNitrS',      'd-1',                  'nitrification rate constant in sediment',                  default=1.0_rk,    scale_factor =1.0_rk/secs_pr_day)
   call self%get_parameter(self%cThetaNitr, 'cThetaNitr',  '[-]',                  'temperature coefficient for nitrification',                default=1.08_rk)
   call self%get_parameter(self%NO3PerC,    'NO3PerC',     '[-]',                  'NO3 denitrified per mol C mineralised',                    default=0.8_rk)
   call self%get_parameter(self%hNO3Denit,  'hNO3Denit',   'mgN/l',                'quadratic half-sat. NO3 conc. for denitrification',        default=2.0_rk)
   call self%get_parameter(self%kPSorp,     'kPSorp',      'd-1',                  'P sorption rate constant not too high -> model speed',     default=0.05_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cRelPAdsD,  'cRelPAdsD',   'gP/gD',                'max. P adsorption per g DW',                               default=0.00003_rk)
   call self%get_parameter(self%cRelPAdsFe, 'cRelPAdsFe',  'gP/gFe',               'max. P adsorption per g Fe',                               default=0.065_rk)
   call self%get_parameter(self%fFeDIM,     'fFeDIM',      'gFe/gD',               'Fe content of inorg. matter',                              default=0.01_rk)
   call self%get_parameter(self%cRelPAdsAl, 'cRelPAdsAl',  'gP/gAl',               'max. P adsorption per g Al',                               default=0.134_rk)
   call self%get_parameter(self%fAlDIM,     'fAlDIM',      'gAl/gD',               'Al content of inorg. matter',                              default=0.01_rk)
   call self%get_parameter(self%fRedMax,    'fRedMax',     '[-]',                  'max. reduction factor of P adsorption affinity',           default=0.9_rk)
   call self%get_parameter(self%cKPAdsOx,   'cKPAdsOx',    'm3/gP',                'P adsorption affinity at oxidized conditions',             default=0.6_rk)
   call self%get_parameter(self%kPChemPO4,  'kPChemPO4',   'd-1',                  'chem. PO4 loss rate',                                      default=0.03_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%coPO4Max,   'coPO4Max',    'mgP/l',                'max. SRP conc. in pore water',                             default=1.0_rk)
   call self%get_parameter(self%bPorS,      'bPorS',       'm3 water m-3 sediment','sediment porosity',                                        default=0.847947_rk)
   call self%get_parameter(self%cThetaDif,  'cThetaDif',   '[-]',                  'temperature coefficient for diffusion',                    default=1.02_rk)
   call self%get_parameter(self%fDepthDifS, 'fDepthDifS',  '[-]',                  'utrient diffusion distance as fraction of sediment depth', default=0.5_rk)
   call self%get_parameter(self%kNDifNH4,   'kNDifNH4',    'm2/day',               'mol. NH4 diffusion constant',                              default=0.000112_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cTurbDifNut,'cTurbDifNut', '[-]',                  'bioturbation factor for diffusion',                        default=5.0_rk)
   call self%get_parameter(self%bPorCorS,   'bPorCorS',    'm3 water m-3 sediment','sediment porosity, corrected for tortuosity',              default=0.737275_rk)
   call self%get_parameter(self%kNDifNO3,   'kNDifNO3',    'm2/day',               'mol. NO3 diffusion constant',                              default=0.000086_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPDifPO4,   'kPDifPO4',    'm2/day',               'mol. PO4 diffusion constant',                              default=0.000072_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kO2Dif,     'kO2Dif',      'm2/day',               'mol. O2 diffusion constant',                               default=0.000026_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cTurbDifO2, 'cTurbDifO2',  '[-]',                  'bioturbation factor for diffusion',                        default=5.0_rk)
   call self%get_parameter(self%kDMinHum,   'kDMinHum',    'd-1',                  'maximum_decomposition_constant_of_humic_material_(1D-5)',  default=0.00001_rk , scale_factor=1.0_rk/secs_pr_day)

!  Register local state variable
!  Inorganic matter
   call self%register_state_variable(self%id_sDIMS,  'sDIMS',  'g m-2','sediment inorg.Matter',  initial_value=39611.3_rk, minimum=_ZERO_)
!  Detritus
   call self%register_state_variable(self%id_sDDetS, 'sDDetS', 'g m-2','sediment detritus DW',   initial_value=181.7_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sNDetS, 'sNDetS', 'g m-2','sediment detritus N',    initial_value=4.54_rk,    minimum=_ZERO_)
   call self%register_state_variable(self%id_sPDetS, 'sPDetS', 'g m-2','sediment detritus P',    initial_value=0.454_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sSiDetS,'sSiDetS','g m-2','sediment detritus Si',   initial_value=1.82_rk,    minimum=_ZERO_)
!  Dissolved nutrients
   call self%register_state_variable(self%id_sPO4S,  'sPO4S',   'g m-2','Sediment Phosphate',    initial_value=0.182_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sPAIMS, 'sPAIMS',  'g m-2','SED_Absorbed Phosphate',initial_value=17.99_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sNH4S,  'sNH4S',   'g m-2','Sediment Amonia',       initial_value=0.02_rk,    minimum=_ZERO_)
   call self%register_state_variable(self%id_sNO3S,  'sNO3S',   'g m-2','Sediment Nitrates',     initial_value=0.002_rk,   minimum=_ZERO_)
!  Humus
   call self%register_state_variable(self%id_sDHumS,'sDHumS','g m-2','sediment Humus DW',        initial_value=3452.34_rk, minimum=_ZERO_)
   call self%register_state_variable(self%id_sNHumS,'sNHumS','g m-2','sediment Humus N',         initial_value=172.62_rk,  minimum=_ZERO_)
   call self%register_state_variable(self%id_sPHumS,'sPHumS','g m-2','sediment Humus P',         initial_value=17.26_rk,   minimum=_ZERO_)
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNH4S)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNO3S)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNDetS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPO4S)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPAIMS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPDetS)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_sSiDetS)
!---------------------------------------------------------------------------------------------------------------
!  register state variables dependencies (2 steps)
!---------------------------------------------------------------------------------------------------------------
!   Register dependencies on external state variables
   call self%register_state_dependency(self%id_O2ConsumpSed, 'oxygen_pool_water',               'g m-3', 'oxygen_pool_water')
   call self%register_state_dependency(self%id_MinSiO2Sed,   'SiO2_generated_by_mineralization','g m-3', 'SiO2_generated_by_mineralization')
   call self%register_state_dependency(self%id_diffNH4,      'NH4_diffusion_flux',              'g m-3', 'NH4_diffusion_flux')
   call self%register_state_dependency(self%id_diffNO3,      'NO3_diffusion_flux',              'g m-3', 'NO3_diffusion_flux')
   call self%register_state_dependency(self%id_diffPO4,      'PO4_diffusion_flux',              'g m-3', 'PO4_diffusion_flux')
!  Register diagnostic variables
   call self%register_diagnostic_variable(self%id_afOxySed,       'afOxySed',       '[-]',       'fraction of aerobic sediment',        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDDetS,        'rPDDetS',        '[-]',       'detritus_P/D_ration_sed',             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDDetS,        'rNDDetS',        '[-]',       'detritus_N/D_ration_sed',             output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aPEqIMS,        'aPEqIMS',        '[-]',       'equilibrium_absorbed_PO4',            output=output_instantaneous)
!  Register diagnostic variables for modular fluxes
!  Total fluxes to local state variables
   call self%register_diagnostic_variable(self%id_tDAbioIMS,      'tDAbioIMS',      'g m-2 s-1', 'abiotic_sediment_DIMS_change',        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAbioDetSflux, 'tDAbioDetSflux', 'g m-2 s-1', 'abiotic_sediment_DDetS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAbioDetS,     'tNAbioDetS',     'g m-2 s-1', 'abiotic_sediment_NDetS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAbioDetS,     'tPAbioDetS',     'g m-2 s-1', 'abiotic_sediment_PDetS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiAbioDetS,    'tSiAbioDetS',    'g m-2 s-1', 'abiotic_sediment_SiDetS_change',      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAbioPO4S,     'tPAbioPO4S',     'g m-2 s-1', 'abiotic_sediment_PAIMS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAbioNH4S,     'tNAbioNH4S',     'g m-2 s-1', 'abiotic_sediment_NH4S_change',        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAbioNO3S,     'tNAbioNO3S',     'g m-2 s-1', 'abiotic_sediment_NO3S_change',        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAbioAIMS,     'tPAbioAIMS',     'g m-2 s-1', 'abiotic_sediment_PAIMS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAbioHumSflux, 'tDAbioHumSflux', 'g m-2 s-1', 'abiotic_sediment_DHumS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAbioHumS,     'tNAbioHumS',     'g m-2 s-1', 'abiotic_sediment_NHumS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAbioHumS,     'tPAbioHumS',     'g m-2 s-1', 'abiotic_sediment_NHumS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPdifPO4,       'tPdifPO4',       'g m-2 s-1', 'abiotic_sediment_PO4W_diffusion',     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNdifNH4,       'tNdifNH4',       'g m-2 s-1', 'abiotic_sediment_NH4W_diffusion',     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNdifNO3,       'tNdifNO3',       'g m-2 s-1', 'abiotic_sediment_NO3W_diffusion',     output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiAbioSiO2S,   'tSiAbioSiO2S',   'g m-2 s-1', 'abiotic_sediment_SiO2_change',        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAbioO2S,      'tDAbioO2S',      'g m-2 s-1', 'abiotic_sediment_O2_consumption',     output=output_instantaneous)
!  Register diagnostic variable for external dependencies
   call self%register_diagnostic_variable(self%id_tDAbioDetS,     'tDAbioDetS',     'g m-2 s-1', 'Detritus change in abiotic sediment', output=output_none)
   call self%register_diagnostic_variable(self%id_tDAbioHumS,     'tDAbioHumS',     'g m-2 s-1', 'Humus change in abiotic sediment',    output=output_none)
!  Register environmental dependencies
   call self%register_dependency(self%id_uTm,standard_variables%temperature)
   call self%register_dependency(self%id_depth,standard_variables%bottom_depth)
   call self%register_dependency(self%id_dz,standard_variables%cell_thickness)


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
! ! INPUT PARAMETERS:
   class (type_au_pclake_abiotic_sediment), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
!  carriers for local state variables
   real(rk)                   :: sNH4S,sNO3S,sPO4S,sPAIMS,sDIMS
   real(rk)                   :: sDDetS,sNDetS,sPDetS,sSiDetS
   real(rk)                   :: sDHumS,sNHumS,sPHumS
!  nutrients ratios
   real(rk)                   :: rPDDetS,rNDDetS
!  carriers for environmental dependencies
   real(rk)                   :: uTm,depth,dz
!  carriers for diagnostic dependencies
!  in abiotic water column module
   real(rk)                   :: sO2W,sNH4W,sNO3W,sPO4W
!  variables for local processes
   real(rk)                   :: kPMinDetS,kNMinDetS,kSiMinDetS
   real(rk)                   :: afOxySed, aDepthOxySed
   real(rk)                   :: tSOD,tDMinDetS,uFunTmNitr,tNNitrS
   real(rk)                   :: oNH4S,oNO3S,oPO4S
   real(rk)                   :: tO2MinDetS, tDMinOxyDetS,tO2NitrS
   real(rk)                   :: tNAbioNO3S,tNDenitS,tDDenitS
   real(rk)                   :: tNAbioNH4S,tNMinDetS,uFunTmMinS
   real(rk)                   :: tPAbioPO4S,tPMinDetS,tPSorpIMS,aPEqIMS
   real(rk)                   :: aPIsoAdsS,aPAdsMaxS,aKPAdsS,tPChemPO4,tPAbioAIMS
   real(rk)                   :: tSiMinDetS,tDAbioDetS
   real(rk)                   :: tNAbioDetS
   real(rk)                   :: tPAbioDetS
   real(rk)                   :: tSiAbioDetS
   real(rk)                   :: tDAbioIMS
!  variables for diffusion(in the order of apperance)
   real(rk)                   :: aDepthDif,tNDifNH4,tNDifNO3,tPDifPO4
   real(rk)                   :: akO2DifCor,tO2Dif,uFunTmDif
!  Variables for humus
   real(rk)    :: tDMinHumS,tNMinHumS,tPMinHumS
   real(rk)    :: tDAbioHumS,tNAbioHumS,tPAbioHumS
!EOP
!-----------------------------------------------------------------------
!BOC

!  Enter spatial loops
   _FABM_HORIZONTAL_LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_HORIZONTAL_(self%id_sDIMS,sDIMS)
   _GET_HORIZONTAL_(self%id_sDDetS,sDDetS)
   _GET_HORIZONTAL_(self%id_sNDetS,sNDetS)
   _GET_HORIZONTAL_(self%id_sPDetS,sPDetS)
   _GET_HORIZONTAL_(self%id_sSiDetS,sSiDetS)
   _GET_HORIZONTAL_(self%id_sPO4S,sPO4S)
   _GET_HORIZONTAL_(self%id_sPAIMS,sPAIMS)
   _GET_HORIZONTAL_(self%id_sNH4S,sNH4S)
   _GET_HORIZONTAL_(self%id_sNO3S,sNO3S)

!  Humus
   _GET_HORIZONTAL_(self%id_sDHumS,sDHumS)
   _GET_HORIZONTAL_(self%id_sNHumS,sNHumS)
   _GET_HORIZONTAL_(self%id_sPHumS,sPHumS)
!     Retrieve dependencis value
!  from abiotic water module
   _GET_(self%id_O2ConsumpSed,sO2W)
   _GET_(self%id_diffNH4,sNH4W)
   _GET_(self%id_diffNO3,sNO3W)
   _GET_(self%id_diffPO4,sPO4W)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_HORIZONTAL_(self%id_depth,depth)
   _GET_(self%id_dz,dz)
   
!-----------------------------------------------------------------------
!  Current local nutrients ratios(check the curent state)
!-----------------------------------------------------------------------
   rPDDetS=sPDetS/(sDDetS+NearZero)
   rNDDetS=sNDetS/(sDDetS+NearZero)
!-----------------------------------------------------------------------
!  Temperature functions for sediment abiotic process
!-----------------------------------------------------------------------
!  temp._function_of_mineralization
   uFunTmMinS=uFunTmAbio(uTm,self%cThetaMinS)
!  Temperature_dependence_for_nitrification
   uFunTmNitr=uFunTmAbio(uTm,self%cThetaNitr)
!  temperature_function_of_diffusion
   uFunTmDif= uFunTmAbio(uTm,self%cThetaDif)
!-----------------------------------------------------------------------
!  dissolved nutrients concentration in sediment(converting)
!-----------------------------------------------------------------------
!  conc._dissolved_N-NO3_in_interstitial_water
   oNO3S = sNO3S / self%cDepthS / self%bPorS
!  conc._dissolved_N-NH4_in_interstitial_water
   oNH4S=sNH4S/self%cDepthS/self%bPorS
!  conc._dissolved_P_in_interstitial_water
   oPO4S = sPO4S / self%cDepthS / self%bPorS
!-----------------------------------------------------------------------
!  Mineralization functions
!-----------------------------------------------------------------------
!  P_mineralisation_constant_in_sed.
   kPMinDetS=self%kDMinDetS
!  N_mineralisation_constant_in_sed.
   kNMinDetS=self%kDMinDetS
!  Si_mineralisation_constant_in_sed.
   kSiMinDetS=self%kDMinDetS
!  decomposition_of_upper_sediment
   tDMinDetS=self%kDMinDetS*uFunTmMinS*sDDetS
!  mineralization_of_P_in_upper_sediment
   tPMinDetS=kPMinDetS*uFunTmMinS*sPDetS
!  mineralization_of_N_in_upper_sediment
   tNMinDetS=kNMinDetS*uFunTmMinS*sNDetS
!  mineralization_of_Si_in_upper_sediment
   tSiMinDetS=kSiMinDetS*uFunTmMinS*sSiDetS
!-----------------------------------------------------------------------
!  diffusion process
!-----------------------------------------------------------------------
!  average_diffusion_distance
   aDepthDif=self%fDepthDifS*self%cDepthS
!  diffusion_flux_of_NH4_from_sediment_to_water
   tNDifNH4=self%kNDifNH4*uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oNH4S-sNH4W)/aDepthDif
!  diffusion_flux_of_NO3_from_sediment_to_water
   tNDifNO3=self%kNDifNO3*uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oNO3S-sNO3W)/aDepthDif
!  diffusion_flux_of_dissolved_P_from_sediment_to_water
   tPDifPO4=self%kPDifPO4*uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oPO4S-sPO4W)/aDepthDif
!  corrected_O2_diffusion_coefficient
   akO2DifCor=self%kO2Dif*uFunTmDif*self%cTurbDifO2*self%bPorCorS
!  O2_diffusion_(water_->_sediment)
   tO2Dif= akO2DifCor*sO2W/aDepthDif
!-----------------------------------------------------------------------
!  Oxygen conditions in sediment
!-----------------------------------------------------------------------
!  sediment_oxygen_demand
   tSOD=(molO2molC*self%cCPerDW*(1.0_rk-self%fRefrDetS)*tDMinDetS+self%O2PerNH4*molO2molN*self%kNitrS*uFunTmNitr*sNH4S)/self%cDepthS
!  oxygen_penetration_depth
   aDepthOxySed=(((2.0_rk * sO2W * akO2DifCor / tSOD) )** (0.5_rk))
!  fraction_aerobic_sediment
   afOxySed=aDepthOxySed/self%cDepthS
!  aerobic_mineralisation
   tDMinOxyDetS=afOxySed*(1.0_rk-self%fRefrDetS)*tDMinDetS
!  sediment_oxygen_demand
   tO2MinDetS=molO2molC*self%cCPerDW*tDMinOxyDetS
!-----------------------------------------------------------------------
!  denitrification flux
!-----------------------------------------------------------------------
!  mineralisation_flux_by_denitrification
   tDDenitS=oNO3S*oNO3S/(self%hNO3Denit*self%hNO3Denit+oNO3S*oNO3S)*(1.0-afOxySed)*(1.0_rk-self%fRefrDetS)*tDMinDetS
!  Denitrification_flux
   tNDenitS=self%NO3PerC*molNmolC*self%cCPerDW*tDDenitS
!-----------------------------------------------------------------------
!  nitrification flux
!-----------------------------------------------------------------------
!  nitrification_flux
   tNNitrS=afOxySed*self%kNitrS*uFunTmNitr*sNH4S
!  O2_flux_due_to_nitrification
   tO2NitrS=self%O2PerNH4*molO2molN*tNNitrS
!-----------------------------------------------------------------------
!  absorbed P in sediment,oxygen dependent
!-----------------------------------------------------------------------
!  max._P_adsorption_per_g_inorg._matter_in_sediment
   aPAdsMaxS =self%cRelPAdsD+afOxySed*self%cRelPAdsFe*self%fFeDIM+self%cRelPAdsAl*self%fAlDIM
!  P_adsorption_affinity,_corrected_for_redox_conditions
   aKPAdsS=(1.0_rk-self%fRedMax*(1.0_rk-afOxySed))*self%cKPAdsOx
!  P_adsorption_isotherm_onto_inorg._matter_in_sediment
   aPIsoAdsS=aPAdsMaxS*aKPAdsS*oPO4S/(1.0_rk+aKPAdsS*oPO4S)
!  equilibrium_amount
   aPEqIMS = aPIsoAdsS * sDIMS
!  sorption
   tPSorpIMS=self%kPSorp*(aPEqIMS-sPAIMS)
!  chem._loss_of_dissolved_P_from_pore_water
   tPChemPO4=max( 0.0_rk,self%kPChemPO4*(oPO4S-self%coPO4Max))
!  decomposition_of_upper_sediment_humus
   tDMinHumS = self%kDMinHum * uFunTmMinS * afOxySed * sDHumS
!  mineralization_of_P_in_upper_sediment_humus
   tPMinHumS = self%kDMinHum * uFunTmMinS * afOxySed * sPHumS
!  mineralization_of_N_in_upper_sediment_humus
   tNMinHumS = self%kDMinHum * uFunTmMinS * afOxySed * sNHumS
!-----------------------------------------------------------------------
!  total abiotic flux for each state variable in sediment
!-----------------------------------------------------------------------
!  total_abiotic/microbial_DW_inorganic_matter_flux_in_sediment
   tDAbioIMS=0.0_rk
!  total_abiotic/microbial_DW_detritus_flux_in_sediment
   tDAbioDetS=-tDMinDetS
!  total_abiotic/microbial_P_detritus_flux_in_sediment
   tPAbioDetS =-tPMinDetS
!  total_abiotic/microbial_dissolved_P_flux_in_sediment
   tPAbioPO4S= (1.0_rk-self%fRefrDetS)*tPMinDetS + tPMinHumS -tPDifPO4-tPSorpIMS  -tPChemPO4
!  total_abiotic/microbial_P_absorbed_onto_inorganic_matter_flux_in_sediment
   tPAbioAIMS=tPSorpIMS
!  total_abiotic/microbial_N_NH4_flux_in_sediment
   tNAbioNH4S=(1.0_rk-self%fRefrDetS)*tNMinDetS +tNMinHumS -tNDifNH4 -tNNitrS
!  total_abiotic/microbial_N_NO3_flux_in_sediment
   tNAbioNO3S= tNNitrS-tNDenitS-tNDifNO3
!  total_abiotic/microbial_N_detritus_flux_in_sediment
   tNAbioDetS =-tNMinDetS
!  total_abiotic/microbial_Si_detritus_flux_in_sediment
   tSiAbioDetS =-tSiMinDetS
!  Humus process
!  total_abiotic/microbial_DW_humus_flux_in_sediment
   tDAbioHumS = self%fRefrDetS * tDMinDetS - tDMinHumS
!  total_abiotic/microbial_N_humus_flux_in_sediment
   tNAbioHumS = self%fRefrDetS * tNMinDetS - tNMinHumS
!  total_abiotic/microbial_P_humus_flux_in_sediment
   tPAbioHumS = self%fRefrDetS * tPMinDetS - tPMinHumS
!-----------------------------------------------------------------------
!  update the state variables
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_sDIMS,tDAbioIMS)
   _SET_ODE_BEN_(self%id_sDDetS,tDAbioDetS)
   _SET_ODE_BEN_(self%id_sPDetS,tPAbioDetS)
   _SET_ODE_BEN_(self%id_sNDetS,tNAbioDetS)
   _SET_ODE_BEN_(self%id_sSiDetS,tSiAbioDetS)
   _SET_ODE_BEN_(self%id_sNH4S,tNAbioNH4S)
   _SET_ODE_BEN_(self%id_sNO3S,tNAbioNO3S)
   _SET_ODE_BEN_(self%id_sPO4S,tPAbioPO4S)
   _SET_ODE_BEN_(self%id_sPAIMS,tPAbioAIMS)
   _SET_ODE_BEN_(self%id_sDHumS,tDAbioHumS)
   _SET_ODE_BEN_(self%id_sPHumS,tPAbioHumS)
   _SET_ODE_BEN_(self%id_sNHumS,tNAbioHumS)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
   _SET_BOTTOM_EXCHANGE_(self%id_MinSiO2Sed,(1.0_rk-self%fRefrDetS)*tSiMinDetS)
   _SET_BOTTOM_EXCHANGE_(self%id_diffNH4, tNdifNH4)
   _SET_BOTTOM_EXCHANGE_(self%id_diffNO3,tNdifNO3)
   _SET_BOTTOM_EXCHANGE_(self%id_diffPO4,tPdifPO4)
   _SET_BOTTOM_EXCHANGE_(self%id_O2ConsumpSed,-tO2MinDetS - tO2NitrS)
!-----------------------------------------------------------------------
!  Output dependent diagnostic variables for other modules
!-----------------------------------------------------------------------
!  output local diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_afOxySed,afOxySed)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rPDDetS,rPDDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rNDDetS,rNDDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aPEqIMS,aPEqIMS)
!  output diagnostic values for external usage
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioDetS,tDAbioDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioHumS,tDAbioHumS)

!  output diagnostic variables for modular fluxes
!  total fluxes for local sediment state variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioIMS,tDAbioIMS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioDetSflux,tDAbioDetS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAbioDetS,tNAbioDetS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAbioDetS,tPAbioDetS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiAbioDetS,tSiAbioDetS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAbioAIMS,tPAbioAIMS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioHumSflux,tDAbioHumS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAbioHumS,tNAbioHumS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAbioHumS,tPAbioHumS*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAbioNH4S,tNAbioNH4S*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAbioNO3S,tNAbioNO3S*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAbioPO4S,tPAbioPO4S*86400.0_rk)
!  total fluxes for abiotic water state variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPdifPO4,tPdifPO4/dz*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNdifNH4,tNdifNH4/dz*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNdifNO3,tNdifNO3/dz*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioO2S,(-tO2MinDetS - tO2NitrS)*86400.0_rk)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiAbioSiO2S,(1.0_rk-self%fRefrDetS)*tSiMinDetS*86400.0_rk)


   _FABM_HORIZONTAL_LOOP_END_
!-----------------------------------------------------------------------
! Spatial loop end
!-----------------------------------------------------------------------

   end subroutine do_bottom

!EOC
!-----------------------------------------------------------------------

   end module au_pclake_abiotic_sediment

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

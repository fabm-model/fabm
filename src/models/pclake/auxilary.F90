#include "fabm_driver.h"
!
!BOP
!
! !INTERFACE:
   module pclake_auxilary
!
! !DESCRIPTION:
!  pclake_auxilary is created for the purpose of computing resuspension and sedimentation,
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
   use pclake_utility, ONLY: uFunTmAbio,uFunTmBio

   implicit none

!  default: all is private.
      private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_auxilary
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

   contains

!     Model procedures
      procedure :: initialize
      procedure :: do_bottom
      procedure :: do
   end type type_pclake_auxilary
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
   class (type_pclake_auxilary), intent(inout), target :: self
   integer,                     intent(in)            :: configunit

!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cDepthS,      'cDepthS',      'm',                 'sediment depth',                                         default=0.1_rk)
   call self%get_parameter(self%cThetaSet,    'cThetaSet',    '1/e^oC',            'temp. parameter of sedimentation',                       default=1.01_rk)
   call self%get_parameter(self%kVegResus,    'kVegResus',    'm2/gDW',            'rel. resuspension reduction per g vegetation',           default=0.01_rk)
   call self%get_parameter(self%kTurbFish,    'kTurbFish',    'g/gfish/d',         'relative resuspension by adult fish browsing',           default=1.0_rk)
!  don't convert here,due to PCLake imperical function. kTurbFish
   call self%get_parameter(self%cSuspRef,     'cSuspRef',     '[-]',               'reference suspended matter function [-]',                default=0.0_rk)
   call self%get_parameter(self%cSuspMin,     'cSuspMin',     '[-]',               'minimum value of logistic function',                     default=6.1_rk)
   call self%get_parameter(self%cSuspMax,     'cSuspMax',     '[-]',               'maximum value of logistic function',                     default=25.2_rk)
   call self%get_parameter(self%cSuspSlope,   'cSuspSlope',   '[-]',               'slope of logistic function',                             default=2.1_rk)
   call self%get_parameter(self%hDepthSusp,   'hDepthSusp',   '[-]',               'half-sat. value of depth in logistic function',          default=2.0_rk)
   call self%get_parameter(self%cFetchRef,    'cFetchRef',    'm',                 'reference fetch',                                        default=1000.0_rk)
   call self%get_parameter(self%cFetch,       'cFetch',       'm',                 'the length of the lake in the prevailing wind direction',default=1000.0_rk)
   call self%get_parameter(self%fLutum ,      'fLutum',       '[-]',               'lutum content of inorg. matter',                         default=0.1_rk)
   call self%get_parameter(self%fLutumRef,    'fLutumRef',    '[-]',               'reference lutum fraction',                               default=0.2_rk)
   call self%get_parameter(self%bPorS,        'bPorS',        'm3waterm-3sediment','porosity',                                               default=0.847947_rk)
   call self%get_parameter(self%cResusPhytExp,'cResusPhytExp','(gDW m-2 d-1)-1',   'exp._par. for phytopl. resuspension',                    default=-0.379_rk)
   call self%get_parameter(self%kResusPhytMax,'kResusPhytMax','d-1',               'max. phytopl. resuspension',                             default=0.25_rk)
   call self%get_parameter(self%cVSetIM,      'cVSetIM',      'm d-1',             'max. sedimentation velocity of inert org. matter',       default=-1.0_rk)
   call self%get_parameter(self%cVSetDet,     'cVSetDet',     'm d-1',             'max. sedimentation velocity of detritus',                default=-0.25_rk)
   call self%get_parameter(self%cVSetDiat,    'cVSetDiat',    'm d-1',             'sedimentation velocity Diatoms',                         default=0.5_rk)
   call self%get_parameter(self%cVSetGren,    'cVSetGren',    'm d-1',             'sedimentation velocity of greens',                       default=0.2_rk)
   call self%get_parameter(self%cVSetBlue ,   'cVSetBlue',    'm d-1',             'sedimentation velocity Blue-greens',                     default=0.06_rk)
   call self%get_parameter(self%cRhoIM,       'cRhoIM',       'g m-3',             'density of sediment IM',                                 default=2500000.0_rk)
   call self%get_parameter(self%cRhoOM,       'cRhoOM',       'g m-3',             'density of sediment detritus',                           default=1400000.0_rk)
   call self%get_parameter(self%fDOrgSoil,    'fDOrgSoil',    '[-]',               'fraction soil organic matter',                           default=0.1_rk)
   call self%get_parameter(self%cPO4Ground,   'cPO4Ground',   'mgP l-1',           'PO4 cone in groundwater',                                default=0.1_rk)
   call self%get_parameter(self%cNH4Ground,   'cNH4Ground',   'mgN l-1',           'NH4 cone in groundwater',                                default=1.0_rk)
   call self%get_parameter(self%cNO3Ground,   'cNO3Ground',   'mgN l-1',           'NO3 cone in groundwater',                                default=0.1_rk)
   call self%get_parameter(self%cTmOptFish,   'cTmOptFish',   '°C',                'optimum temp. of fish',                                  default=25.0_rk)
   call self%get_parameter(self%cSigTmFish,   'cSigTmFish',   '°C',                'temperature constant of fish(sigma_in_Gaussian_curve)',  default=10.0_rk)
   call self%get_parameter(self%cNDDiatMin,   'cNDDiatMin',   'mgN/mgDW',          'minimum N/day ratio Diatoms',                            default=0.01_rk)
   call self%get_parameter(self%cPDDiatMin,   'cPDDiatMin',   'mgP/mgDW',          'minimum P/day ratio Diatoms',                            default=0.0005_rk)
   call self%get_parameter(self%cNDGrenMin,   'cNDGrenMin',   'mgN/mgDW',          'minimum N/day ratio greens',                             default=0.02_rk)
   call self%get_parameter(self%cPDGrenMin,   'cPDGrenMin',   'mgP/mgDW',          'minimum P/day ratio greens',                             default=0.0015_rk)
   call self%get_parameter(self%cNDBlueMin,   'cNDBlueMin',   'mgN/mgDW',          'minimum N/day ratio Bluegreens',                         default=0.03_rk)
   call self%get_parameter(self%cPDBlueMin,   'cPDBlueMin',   'mgP/mgDW',          'minimum P/day ratio Bluegreens',                         default=0.0025_rk)
   call self%get_parameter(self%cNDBlueMax,   'cNDBlueMax',   'mgN/mgDW',          'max. N/day ratio Bluegreens',                            default=0.15_rk)
   call self%get_parameter(self%cNDDiatMax,   'cNDDiatMax',   'mgN/mgDW',          'max. N/day ratio Diatoms',                               default=0.05_rk)
   call self%get_parameter(self%cNDGrenMax,   'cNDGrenMax',   'mgN/mgDW',          'max. N/day ratio greens',                                default=0.1_rk)
   call self%get_parameter(self%cPDBlueMax,   'cPDBlueMax',   'mgP/mgDW',          'max. P/day ratio blue-greens',                           default=0.025_rk)
   call self%get_parameter(self%cPDDiatMax,   'cPDDiatMax',   'mgP/mgDW',          'max. P/day ratio Diatoms',                               default=0.005_rk)
   call self%get_parameter(self%cPDGrenMax,   'cPDGrenMax',   'mgP/mgDW',          'max. P/day ratio greens',                                default=0.015_rk)
!  loadings
   call self%get_parameter(self%cLoadPO4,     'cLoadPO4',     'gP m-2 d-1',         'constant, different from cLoadIM',                       default=0.005_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cLoadNO3,     'cLoadNO3',     'gN m-2 d-1',         'constant',                                               default=0.05_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%uQIn,         'uQIn',         'mm d-1',             'net water load rate, constant',                          default=20.0_rk,   scale_factor=1.0_rk/secs_pr_day)
!  erosion parameter
   call self%get_parameter(self%fSedErosIM,   'fSedErosIM',   '[-]',                'instantly_sedimentating_fraction_of_IM',                 default=0.95_rk)
   call self%get_parameter(self%cDErosTot,    'cDErosTot',    'g m-2 d-1',          'erosion input(tentative)',                               default=0.1_rk,    scale_factor=1.0_rk/secs_pr_day)
!  Register dependencies to abiotic sediment module
   call self%register_state_dependency(self%id_WSNH4,   'Amonia_pool_in_sediment',       'g m-2', 'Amonia_pool_in_sediment')
   call self%register_state_dependency(self%id_WSNO3,   'Nitrates_pool_in_sediment',     'g m-2', 'Nitrates_pool_in_sediment')
   call self%register_state_dependency(self%id_WSPO4,   'Phosphate_pool_in_sediment',    'g m-2', 'Phosphate_pool_in_sediment')
   call self%register_state_dependency(self%id_WSPAIM,  'Absorbed_phosphrus_in_sediment','g m-2', 'Absorbed_phosphrus_in_sediment')
   call self%register_state_dependency(self%id_WSDIM,   'Inorg_pool_in_sediment',        'g m-2', 'Inorg_pool_in_sediment')
   call self%register_state_dependency(self%id_WSDDet,  'Detritus_DW_in_sediment',       'g m-2', 'Detritus_DW_in_sediment')
   call self%register_state_dependency(self%id_WSNDet,  'Detritus_N_in_sediment',        'g m-2', 'Detritus_N_in_sediment')
   call self%register_state_dependency(self%id_WSPDet,  'Detritus_P_in_sediment',        'g m-2', 'Detritus_P_in_sediment')
   call self%register_state_dependency(self%id_WSSiDet, 'Detritus_Si_in_sediment',       'g m-2', 'Detritus_Si_in_sediment')
   call self%register_state_dependency(self%id_WSDHum,  'Humus_DW_in_sediment',          'g m-2', 'Humus_DW_in_sediment')
   call self%register_state_dependency(self%id_WSNHum,  'Humus_N_in_sediment',           'g m-2', 'Humus_N_in_sediment')
   call self%register_state_dependency(self%id_WSPHum,  'Humus_P_in_sediment',           'g m-2', 'Humus_P_in_sediment')
!  Register dependencies to phytoplankton in water column
   call self%register_state_dependency(self%id_SWDDiat, 'Diatom_DW_in_water',            'g m-3', 'Diatom_DW_in_water')
   call self%register_state_dependency(self%id_SWDGren, 'Green_DW_in_water',             'g m-3', 'Green_DW_in_water')
   call self%register_state_dependency(self%id_SWDBlue, 'Blue_DW_in_water',              'g m-3', 'Blue_DW_in_water')
   call self%register_state_dependency(self%id_SWNDiat, 'Diatom_N_in_water',             'g m-3', 'Diatom_N_in_water')
   call self%register_state_dependency(self%id_SWNGren, 'Green_N_in_water',              'g m-3', 'Green_N_in_water')
   call self%register_state_dependency(self%id_SWNBlue, 'Blue_N_in_water',               'g m-3', 'Blue_N_in_water')
   call self%register_state_dependency(self%id_SWPDiat, 'Diatom_P_in_water',             'g m-3', 'Diatom_P_in_water')
   call self%register_state_dependency(self%id_SWPGren, 'Green_P_in_water',              'g m-3', 'Green_P_in_water')
   call self%register_state_dependency(self%id_SWPBlue, 'Blue_P_in_water',               'g m-3', 'Blue_P_in_water')
!  Register dependencies to phytoplankton in sediment
   call self%register_state_dependency(self%id_WSDDiat, 'Diatom_DW_in_sediment',         'g m-2', 'Diatom_DW_in_sediment')
   call self%register_state_dependency(self%id_WSDGren, 'Green_DW_in_sediment',          'g m-2', 'Green_DW_in_sediment')
   call self%register_state_dependency(self%id_WSDBlue, 'Blue_DW_in_sediment',           'g m-2', 'Blue_DW_in_sediment')
   call self%register_state_dependency(self%id_WSNDiat, 'Diatom_N_in_sediment',          'g m-2', 'Diatom_N_in_sediment')
   call self%register_state_dependency(self%id_WSNGren, 'Green_N_in_sediment',           'g m-2', 'Green_N_in_sediment')
   call self%register_state_dependency(self%id_WSNBlue, 'Blue_N_in_sediment',            'g m-2', 'Blue_N_in_sediment')
   call self%register_state_dependency(self%id_WSPDiat, 'Diatom_P_in_sediment',          'g m-2', 'Diatom_P_in_sediment')
   call self%register_state_dependency(self%id_WSPGren, 'Green_P_in_sediment',           'g m-2', 'Green_P_in_sediment')
   call self%register_state_dependency(self%id_WSPBlue, 'Blue_P_in_sediment',            'g m-2', 'Blue_P_in_sediment')
!  register vegetation and fish for resuspension dependency
   call self%register_state_dependency(self%id_DragVeg, 'vegetation_DW',                 'g m-2', 'vegetation_DW')
   call self%register_state_dependency(self%id_TurbFish,'adult_fish_DW',                 'g m-3', 'adult_fish_DW')
!  register zooplankton for transport purpose
   call self%register_state_dependency(self%id_DTranZoo,'Zooplankton_DW',                'g m-3', 'Zooplankton_DW')
   call self%register_state_dependency(self%id_PTranZoo,'Zooplankton_P',                 'g m-3', 'Zooplankton_P')
   call self%register_state_dependency(self%id_NTranZoo,'Zooplankton_N',                 'g m-3', 'Zooplankton_N')
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,    standard_variables%temperature)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
   call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_shear,  standard_variables%bottom_stress)
!  step1, Register dependencies on external diagnostic variables
   call self%register_dependency(self%id_tDAbioDetS, 'Detritus_abiotic_update', 'g m-2 s-1', 'Detritus_abiotic_update')
   call self%register_dependency(self%id_tDAbioHumS, 'Humus_abiotic_update',    'g m-2 s-1', 'Humus_abiotic_update')
   call self%register_dependency(self%id_tDPrimDetS, 'Detritus_from_algae',     '[-]',       'Detritus_from_algae')
   call self%register_dependency(self%id_tDWebDetS,  'Detritus_from_foodweb',   '[-]',       'Detritus_from_foodweb')
   call self%register_dependency(self%id_tDBedDetS,  'Detritus_from_vegetation','[-]',       'Detritus_from_vegetation')
!  Regirster diagnostic variables
   call self%register_diagnostic_variable(self%id_tDBurIM,     'tDBurIM',     'g m-2 s-1','tDBurIM',     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_shearstress, 'shearstress', 'N/m**2',   'shearstress', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aFunDimSusp, 'aFunDimSusp', '[-]',      'aFunDimSusp', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDResusDead, 'tDResusDead', '[-]',      'tDResusDead', output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_aFunTauSet,  'aFunTauSet',  '[-]',      'basic settling rate',  output=output_time_step_integrated)

   return

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
   class (type_pclake_auxilary), intent(in)    :: self
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
   if (uTm >= 0.1_rk) then
       aFunDimSusp=self%cSuspRef*((self%cSuspMin+self%cSuspMax/(1.0_rk+exp(self%cSuspSlope*&
       & (sDepthW-self%hDepthSusp))))*((((self%cFetch +NearZero)/ self%cFetchRef) )** (0.5_rk)))
   else
       aFunDimSusp=0.0_rk
   endif
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
!  update diagnostic variables
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
   class (type_pclake_auxilary), intent(in)    :: self
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


   end module pclake_auxilary

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

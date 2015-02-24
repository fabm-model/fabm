#include "fabm_driver.h"
   module au_pclake_macrophytes
!
! !DESCRIPTION:
!  This module describes the submerged macrophytes group, and implemented as 0d benthic
!  module keeping most of the original PCLake methods.
!  The state variables include: sDVeg,sNVeg,sPVeg(vegetation mass in dry-weight, nitrogen
!  and phosphorus. The related processes are:
!  Assimilation(only for sDVeg), nutrient updake(only for sNVeg and sPVeg),respiration
!  (only for sDVeg), excretion(only for sNVeg and sPVeg), mortality, migration.
!  This module also discribes the processes which influence the state variables registered in
!  other modules, including:
!  nutrients taken up and excreted by shoots and roots: sNH4W&sNH4S<==>sNVeg,sPO4W&sPO4S<==>sPVeg
!                                                       sNO3W&sNO3W==>sNVeg(only taken up)
! detritus morted by shoots and roots: sDVeg==>sDDetW&sDDetS,sNVeg==>sNDetW&sNDetS,sPVeg==>sPDetW&sPDetS
!  oxygen produced by primary production and consumed by respiration: sO2W<==>sDVeg
!  This module also provide important diagnostic variable will be used in other modules, including:
!  Submerged macrophytes biomass, aDSubVeg, used by module: foodweb water module
!  macrophytes coverage percentage,aCovVeg, used by module: foodweb water module
!  Sediment detritus change, tDBedDetS, used bymodule:auxilary
! !USES:
   use fabm_types
   use fabm_expressions
   use fabm_standard_variables
   use au_pclake_utility, ONLY:uFunTmVeg

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_au_pclake_macrophytes
!     local state variable identifers
!     id_sDVeg,macrophytes in dry-weight, gDW/m**2
!     id_sPVeg,macrophytes in nitrogen element, gN/m**2
!     id_sPVeg,macrophytes in phosphorus element, gP/m**2
      type (type_bottom_state_variable_id)            :: id_sDVeg,id_sNVeg,id_sPVeg
!     diagnostic variables for local output
!     id_tO2BedW,oxygen production from macrophytes
      type (type_horizontal_diagnostic_variable_id)       :: id_tO2BedW
!     diagnostic variables for dependencies(without output)
      type (type_horizontal_diagnostic_variable_id)       :: id_aDSubVeg,id_aCovVeg,id_tDBedDetS,id_afCovSurfVeg

!     diagonostic variables for light attenuation coefficient for plant.
      type (type_horizontal_diagnostic_variable_id)       :: id_aDayInitVeg,id_tDBedVeg
!     state dependencies identifers
      type (type_state_variable_id)                :: id_NH4poolW,id_NO3poolW,id_PO4poolW,id_O2poolW
      type (type_state_variable_id)                :: id_DDetpoolW,id_DNetpoolW,id_DPetpoolW
      type (type_bottom_state_variable_id)         :: id_NH4poolS,id_NO3poolS,id_PO4poolS
      type (type_bottom_state_variable_id)         :: id_DDetpoolS,id_DNetpoolS,id_DPetpoolS
!     environmental dependencies
      type (type_global_dependency_id)         :: id_Day
      type (type_dependency_id)                :: id_uTm,id_extc,id_dz
      type (type_dependency_id)                :: id_par,id_meanpar
      type (type_horizontal_dependency_id)     :: id_sDepthW
!     output light variables
      type (type_horizontal_diagnostic_variable_id)       :: id_aLPAR1Veg,id_aLPAR2Veg
      type (type_horizontal_diagnostic_variable_id)       :: id_allimveg,id_aNutLimVeg
      type (type_horizontal_diagnostic_variable_id)       :: id_macroextinction
!     diagnostic dependencies
      type (type_horizontal_dependency_id)     :: id_afOxySed
!    Model parameters
!    Primary production parameters
     real(rk)    :: cDVegIn,kMigrVeg,cMuMaxVeg,cDCarrVeg,kDRespVeg
     real(rk)    :: hLRefVeg,fDepth1Veg,fDepth2Veg
!    nutrient ratio parameters
     real(rk)    :: cNDVegMin,cNDVegMax,cPDVegMin,cPDVegMax
!    vegetation shoots and roots allocation parameters
     real(rk)    :: cDayWinVeg,cLengAllo
     real(rk)    :: fRootVegWin,fRootVegSum,cTmInitVeg
     real(rk)    :: fEmergVeg,fFloatVeg,cDLayerVeg,cCovSpVeg
     real(rk)    :: kMortVegSum,cLengMort,fWinVeg
!    temperature function parameters
     real(rk)    :: cQ10ProdVeg,cQ10RespVeg
!    parameters for Nitrogen and Phosphrus equations
     real(rk)    :: cPDVeg0,cNDVeg0
     real(rk)    :: fSedUptVegMax,fSedUptVegCoef,fSedUptVegExp,cAffNUptVeg,cVNUptMaxVeg
     integer     :: UseEmpUpt
     real(rk)    :: cVPUptMaxVeg,cAffPUptVeg,fDissMortVeg,fDetWMortVeg
!    paremters for sediment properties(pore water concentration)
     real(rk)   :: cDepthS,bPorS,cCPerDW,hO2BOD
     real(rk)   :: cExtSpVeg  !,host
!    plant height
     real(rk)   :: cHeightVeg
   contains

!     Model procedure
      procedure :: initialize
      procedure :: do_bottom
      procedure :: get_light_extinction
   end type type_au_pclake_macrophytes
   
!  private data memebers(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk),parameter :: Pi=3.14159265358979_rk
!  ratio of mol.weights, = 32/12 [gO2/gC], 
   real(rk),parameter :: molO2molC = 2.6667_rk
!  mol_O2_formed_per_mol_NO3-_ammonified
   real(rk),parameter ::O2PerNO3 = 1.5_rk
!  ratio of mol.weights,32/14 [gO2/gN], 
   real(rk),parameter :: molO2molN = 2.2857_rk
!   lowest state variable value
    real(rk),parameter :: VegZero=0.0001_rk
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
   class (type_au_pclake_macrophytes), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!----------------------------------------------------------------------------
!  LOCAL VARIABLES:
!----------------------------------------------------------------------------
   real(rk)                  :: sDVeg_initial
   real(rk)                  :: sNVeg_initial
   real(rk)                  :: sPVeg_initial
   integer                   :: UseEmpUpt
   real(rk)                  :: cDVegIn
   real(rk)                  :: kMigrVeg
   real(rk)                  :: cNDVegMin
   real(rk)                  :: cNDVegMax
   real(rk)                  :: cPDVegMin
   real(rk)                  :: cPDVegMax
   real(rk)                  :: cMuMaxVeg
   real(rk)                  :: cDCarrVeg
   real(rk)                  :: kDRespVeg
   real(rk)                  :: cDayWinVeg
   real(rk)                  :: cLengAllo
   real(rk)                  :: fRootVegWin
   real(rk)                  :: fRootVegSum
   real(rk)                  :: cTmInitVeg
   real(rk)                  :: fEmergVeg
   real(rk)                  :: fFloatVeg
   real(rk)                  :: cDLayerVeg
   real(rk)                  :: cCovSpVeg
   real(rk)                  :: hLRefVeg
   real(rk)                  :: cQ10ProdVeg
   real(rk)                  :: cQ10RespVeg
   real(rk)                  :: kMortVegSum
   real(rk)                  :: cLengMort
   real(rk)                  :: fWinVeg
   real(rk)                  :: cPDVeg0
   real(rk)                  :: cNDVeg0
   real(rk)                  :: fSedUptVegMax
   real(rk)                  :: fSedUptVegCoef
   real(rk)                  :: fSedUptVegExp
   real(rk)                  :: cAffNUptVeg
   real(rk)                  :: cVNUptMaxVeg
   real(rk)                  :: cDepthS
   real(rk)                  :: bPorS
   real(rk)                  :: cVPUptMaxVeg
   real(rk)                  :: cAffPUptVeg
   real(rk)                  :: fDissMortVeg
   real(rk)                  :: cCPerDW
   real(rk)                  :: fDetWMortVeg
   real(rk)                  :: fDepth1Veg
   real(rk)                  :: fDepth2Veg
   real(rk)                  :: hO2BOD
   real(rk)                  :: cExtSpVeg
   real(rk)                  :: cHeightVeg
   character(len=64)         :: ammonia_pool_water
   character(len=64)         :: nitrate_pool_water
   character(len=64)         :: phosphate_pool_water
   character(len=64)         :: ammonia_pool_sediment
   character(len=64)         :: nitrate_pool_sediment
   character(len=64)         :: phosphate_pool_sediment
   character(len=64)         :: oxygen_pool_water
   character(len=64)         :: detritus_DW_pool_water
   character(len=64)         :: detritus_N_pool_water
   character(len=64)         :: detritus_P_pool_water
   character(len=64)         :: detritus_DW_pool_sediment
   character(len=64)         :: detritus_N_pool_sediment
   character(len=64)         :: detritus_P_pool_sediment
   character(len=64)         :: oxic_layer_value
!----------------------------------------------------------------------------
!  create namelist
!----------------------------------------------------------------------------
   namelist /au_pclake_macrophytes/ sDVeg_initial,sNVeg_initial,sPVeg_initial,UseEmpUpt,cDVegIn,kMigrVeg,cNDVegMin,cNDVegMax,cPDVegMin,cPDVegMax, &
                              & cMuMaxVeg,cDCarrVeg,kDRespVeg,cDayWinVeg,cLengAllo,fRootVegWin,fRootVegSum,cTmInitVeg, &
                              & fEmergVeg,fFloatVeg,cDLayerVeg,cCovSpVeg,hLRefVeg,cQ10ProdVeg,cQ10RespVeg,kMortVegSum,cLengMort,&
                              & fWinVeg,cPDVeg0,cNDVeg0,fSedUptVegMax,fSedUptVegCoef,fSedUptVegExp,cAffNUptVeg,cVNUptMaxVeg,&
                              & cDepthS,bPorS,cVPUptMaxVeg,cAffPUptVeg,fDissMortVeg,cCPerDW,fDetWMortVeg, &
                              & fDepth1Veg,fDepth2Veg ,hO2BOD, cExtSpVeg,cHeightVeg,&
                              & ammonia_pool_water,nitrate_pool_water,phosphate_pool_water,&
                              & ammonia_pool_sediment,nitrate_pool_sediment,phosphate_pool_sediment,&
                              & oxygen_pool_water,detritus_DW_pool_water,detritus_N_pool_water,detritus_P_pool_water,detritus_DW_pool_sediment,&
                              & detritus_N_pool_sediment,detritus_P_pool_sediment,oxic_layer_value
!EOP
!-----------------------------------------------------------------------
!BOC
!  initialize the parameters
   sDVeg_initial=1.0_rk
   sNVeg_initial=0.02_rk
   sPVeg_initial=0.002_rk
   UseEmpUpt= 0
   cDVegIn=1.0_rk
   kMigrVeg=0.001_rk
   cNDVegMin=0.01_rk
   cNDVegMax=0.035_rk
   cPDVegMin=0.0008_rk
   cPDVegMax=0.0035_rk
   cMuMaxVeg= 0.2_rk
   cDCarrVeg=400.0_rk
   kDRespVeg=0.02_rk
   cDayWinVeg=259.0_rk
   cLengAllo=15.0_rk  
   fRootVegWin=0.6_rk 
   fRootVegSum=0.1_rk 
   cTmInitVeg=9.0_rk  
   fEmergVeg=0.0_rk 
   fFloatVeg=0.0_rk 
   cDLayerVeg=0.0_rk
   cCovSpVeg=0.5_rk
   hLRefVeg=17.0_rk
   cQ10ProdVeg=1.2_rk
   cQ10RespVeg=2.0_rk
   kMortVegSum=0.005_rk
   cLengMort=42.0_rk
   fWinVeg=0.3_rk
   cPDVeg0=0.002_rk
   cNDVeg0= 0.02_rk
   fSedUptVegMax=0.998_rk 
   fSedUptVegCoef= 2.66_rk
   fSedUptVegExp =-0.83_rk
   cAffNUptVeg=0.2_rk
   cVNUptMaxVeg=0.1_rk
   cVPUptMaxVeg=0.01_rk 
   cAffPUptVeg=0.2_rk
   fDissMortVeg=0.25_rk
   cCPerDW=0.4_rk
   fDetWMortVeg=0.1_rk
   fDepth1Veg=0.0_rk
   fDepth2Veg=1.0_rk
   hO2BOD = 1.0_rk
   cExtSpVeg=0.01_rk
   cHeightVeg=0.1_rk
   ammonia_pool_water=''
   nitrate_pool_water=''
   phosphate_pool_water=''
   ammonia_pool_sediment=''
   nitrate_pool_sediment=''
   phosphate_pool_sediment=''
   oxygen_pool_water=''
   detritus_DW_pool_water=''
   detritus_N_pool_water=''
   detritus_P_pool_water=''
   detritus_DW_pool_sediment=''
   detritus_N_pool_sediment=''
   detritus_P_pool_sediment=''
   oxic_layer_value=''
!  Read parameters namelist
   if (configunit>0) read(configunit,nml=au_pclake_macrophytes,err=99,end=100)
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cDVegIn,'cDVegIn',default=cDVegIn)
   call self%get_parameter(self%kMigrVeg,'kMigrVeg',default=kMigrVeg,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cNDVegMin,'cNDVegMin',default=cNDVegMin)
   call self%get_parameter(self%cNDVegMax,'cNDVegMax',default=cNDVegMax)
   call self%get_parameter(self%cPDVegMin,'cPDVegMin',default=cPDVegMin)
   call self%get_parameter(self%cPDVegMax,'cPDVegMax',default=cPDVegMax)
   call self%get_parameter(self%cMuMaxVeg,'cMuMaxVeg',default=cMuMaxVeg,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrVeg,'cDCarrVeg',default=cDCarrVeg)
   call self%get_parameter(self%kDRespVeg,'kDRespVeg',default=kDRespVeg,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDayWinVeg,'cDayWinVeg',default=cDayWinVeg)
   call self%get_parameter(self%cLengAllo,'cLengAllo',default=cLengAllo)
   call self%get_parameter(self%fRootVegWin,'fRootVegWin',default=fRootVegWin)
   call self%get_parameter(self%fRootVegSum,'fRootVegSum',default=fRootVegSum)
   call self%get_parameter(self%cTmInitVeg,'cTmInitVeg',default=cTmInitVeg)
   call self%get_parameter(self%fEmergVeg,'fEmergVeg',default=fEmergVeg)
   call self%get_parameter(self%fFloatVeg,'fFloatVeg',default=fFloatVeg)
   call self%get_parameter(self%cDLayerVeg,'cDLayerVeg',default=cDLayerVeg)
   call self%get_parameter(self%cCovSpVeg,'cCovSpVeg',default=cCovSpVeg)
   call self%get_parameter(self%hLRefVeg,'hLRefVeg',default=hLRefVeg)
   call self%get_parameter(self%cQ10ProdVeg,'cQ10ProdVeg',default=cQ10ProdVeg)
   call self%get_parameter(self%cQ10RespVeg,'cQ10RespVeg',default=cQ10RespVeg)
   call self%get_parameter(self%kMortVegSum,'kMortVegSum',default=kMortVegSum,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cLengMort,'cLengMort',default=cLengMort)
   call self%get_parameter(self%fWinVeg,'fWinVeg',default=fWinVeg)
   call self%get_parameter(self%cPDVeg0,'cPDVeg0',default=cPDVeg0)
   call self%get_parameter(self%cNDVeg0,'cNDVeg0',default=cNDVeg0)
   call self%get_parameter(self%fSedUptVegMax,'fSedUptVegMax',default=fSedUptVegMax)
   call self%get_parameter(self%fSedUptVegCoef,'fSedUptVegCoef',default=fSedUptVegCoef)
   call self%get_parameter(self%fSedUptVegExp,'fSedUptVegExp',default=fSedUptVegExp)
   call self%get_parameter(self%cAffNUptVeg,'cAffNUptVeg',default=cAffNUptVeg,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVNUptMaxVeg,'cVNUptMaxVeg',default=cVNUptMaxVeg,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDepthS,'cDepthS',default=cDepthS)
   call self%get_parameter(self%bPorS,'bPorS',default=bPorS)
   call self%get_parameter(self%UseEmpUpt,'UseEmpUpt',default=UseEmpUpt)
   call self%get_parameter(self%cAffPUptVeg,'cAffPUptVeg',default=cAffPUptVeg,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVPUptMaxVeg,'cVPUptMaxVeg',default=cVPUptMaxVeg,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortVeg,'fDissMortVeg',default=fDissMortVeg)
   call self%get_parameter(self%cCPerDW,'cCPerDW',default=cCPerDW)
   call self%get_parameter(self%fDetWMortVeg,'fDetWMortVeg',default=fDetWMortVeg)
   call self%get_parameter(self%fDepth1Veg,'fDepth1Veg',default=fDepth1Veg)
   call self%get_parameter(self%fDepth2Veg,'fDepth2Veg',default=fDepth2Veg)
   call self%get_parameter(self%hO2BOD,'hO2BOD',default=hO2BOD)
   call self%get_parameter(self%cHeightVeg,'cHeightVeg',default=cHeightVeg)
!  for vegetation light attenuation
   call self%get_parameter(self%cExtSpVeg,'cExtSpVeg',default=cExtSpVeg)
!  Register local state variable
   call self%register_state_variable(self%id_sDVeg,'sDVeg','g/m**2','vegetation_dry_weight',     &
                                    sDVeg_initial,minimum=VegZero)
   call self%register_state_variable(self%id_sNVeg,'sNVeg','g/m**2','vegetation_Nitrogen',     &
                                    sNVeg_initial,minimum=VegZero)
   call self%register_state_variable(self%id_sPVeg,'sPVeg','g/m**2','vegetation_Phosphrus',     &
                                    sPVeg_initial,minimum=VegZero)
!  register diagnostic variables for local output(with different ways of output)
   call self%register_diagnostic_variable(self%id_tO2BedW,'tO2BedW','g/m**2/s','vegetation_O2_change',           &
                 output=output_none)
!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_aDSubVeg,'aDSubVeg','g/m**2','aDSubVeg',           &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_aCovVeg,'aCovVeg','%','aCovVeg',           &
                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_tDBedDetS,'tDBedDetS','g/m**2/s','tDBedDetS',           &
                     output=output_none)
   call self%register_diagnostic_variable(self%id_afCovSurfVeg,'afCovSurfVeg','--','afCovSurfVeg',           &
                     output=output_none)
   call self%register_diagnostic_variable(self%id_aDayInitVeg,'aDayInitVeg','day','aDayInitVeg',           &
                     output=output_none)
   call self%register_diagnostic_variable(self%id_tDBedVeg,'tDBedVeg','g/m**2/s','tDBedVeg',           &
                     output=output_none)
   call self%register_diagnostic_variable(self%id_aNutLimVeg,'aNutLimVeg','--','aNutLimVeg',  &
                  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_macroextinction,'macroextinction','--','macroextinction',  &
                  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_allimveg,'allimveg','--','allimveg',  &
                  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aLPAR1Veg,'aLPAR1Veg','--','aLPAR1Veg',  &
                  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aLPAR2Veg,'aLPAR2Veg','--','aLPAR2Veg',  &
                  output=output_instantaneous)
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNVeg)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPVeg)
!  regirster state variables dependencies (2 steps) (the external link to sediment is unsolved)
!  step1, Register dependencies on external state variables
   call self%register_state_dependency(self%id_NH4poolW, 'ammonia_pool_water','g/m**3','Amonia pool for nutrient uptake')
   call self%register_state_dependency(self%id_NO3poolW, 'nitrate_pool_water','g/m**3','Nitrate pool for nutrient uptake')
   call self%register_state_dependency(self%id_PO4poolW, 'phosphate_pool_water','g/m**3','Phosphate pool for nutrient uptake')
   call self%register_state_dependency(self%id_O2poolW, 'oxygen_pool_water','g/m**3','oxygen pool in water')
   call self%register_state_dependency(self%id_DDetpoolW, 'detritus_DW_pool_water','g/m**3','detritus_DW_pool_water')
   call self%register_state_dependency(self%id_DNetpoolW, 'detritus_N_pool_water','g/m**3','detritus_N_pool_water')
   call self%register_state_dependency(self%id_DPetpoolW, 'detritus_P_pool_water','g/m**3','detritus_P_pool_water')
   call self%register_state_dependency(self%id_NH4poolS, 'ammonia_pool_sediment','g/m**2','Amonia pool for nutrient uptake')
   call self%register_state_dependency(self%id_NO3poolS, 'nitrate_pool_sediment','g/m**2','Nitrate pool for nutrient uptake')
   call self%register_state_dependency(self%id_PO4poolS, 'phosphate_pool_sediment','g/m**2','Phosphate pool for nutrient uptake')
   call self%register_state_dependency(self%id_DDetpoolS, 'detritus_DW_pool_sediment','g/m**2','detritus_DW_pool_sediment')
   call self%register_state_dependency(self%id_DNetpoolS, 'detritus_N_pool_sediment','g/m**2','detritus_N_pool_sediment')
   call self%register_state_dependency(self%id_DPetpoolS, 'detritus_P_pool_sediment','g/m**2','detritus_P_pool_sediment')
!  step 2, Automatically couple dependencies if target variables have been specified.
   if (ammonia_pool_water/='') call self%request_coupling(self%id_NH4poolW,ammonia_pool_water)
   if (nitrate_pool_water/='') call self%request_coupling(self%id_NO3poolW,nitrate_pool_water)
   if (phosphate_pool_water/='') call self%request_coupling(self%id_PO4poolW,phosphate_pool_water)
   if (oxygen_pool_water/='') call self%request_coupling(self%id_O2poolW,oxygen_pool_water)
   if (detritus_DW_pool_water/='') call self%request_coupling(self%id_DDetpoolW,detritus_DW_pool_water)
   if (detritus_N_pool_water/='') call self%request_coupling(self%id_DNetpoolW,detritus_N_pool_water)
   if (detritus_P_pool_water/='') call self%request_coupling(self%id_DPetpoolW,detritus_P_pool_water)
   if (ammonia_pool_sediment/='') call self%request_coupling(self%id_NH4poolS,ammonia_pool_sediment)
   if (nitrate_pool_sediment/='') call self%request_coupling(self%id_NO3poolS,nitrate_pool_sediment)
   if (phosphate_pool_sediment/='') call self%request_coupling(self%id_PO4poolS,phosphate_pool_sediment)
   if (detritus_DW_pool_sediment/='') call self%request_coupling(self%id_DDetpoolS,detritus_DW_pool_sediment)
   if (detritus_N_pool_sediment/='') call self%request_coupling(self%id_DNetpoolS,detritus_N_pool_sediment)
   if (detritus_P_pool_sediment/='') call self%request_coupling(self%id_DPetpoolS,detritus_P_pool_sediment)
!------------------------------------------------------------------------------------------------------------
!  register environmental dependencies
!------------------------------------------------------------------------------------------------------------
   call self%register_dependency(self%id_uTm,standard_variables%temperature)
   call self%register_dependency(self%id_extc,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_Day,standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
   call self%register_dependency(self%id_dz,standard_variables%cell_thickness)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_meanpar,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk))
!------------------------------------------------------------------------------------------------------------
!  register diagnostic dependencies, 2 steps
!------------------------------------------------------------------------------------------------------------
!  step1, Register dependencies on external diagnostic variables
   call self%register_dependency(self%id_afOxySed, 'oxic_layer_value','--','oxic_layer_value')
  
!  step 2, Automatically couple dependencies if target variables have been specified.
   if (oxic_layer_value/='') call self%request_coupling(self%id_afOxySed,oxic_layer_value)
   return

99  call self%fatal_error('au_pclake_macrophytes_init','Error reading namelist au_pclake_macrophytes')

100 call self%fatal_error('au_pclake_macrophytes_init','Namelist au_pclake_macrophytes was not found.')


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
   class (type_au_pclake_macrophytes), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
! !LOCAL VARIABLES:
!  Carriers for environment dependencies
   real(rk)     :: uTm,extc,meanpar,dz
   real(rk)     :: day,sDepthW
!  state variable value carriers
   real(rk)     :: sDVeg,sNVeg,sPVeg
!  external state variable carriers
   real(rk)     :: sNH4W,sNO3W,sPO4W,sNH4S,sNO3S,sPO4S
   real(rk)     :: sO2W
   real(rk)     :: aCorO2BOD,afOxySed
!  variables with o- prefix
   real(rk)     :: oNH4S,oNO3S,oNDissS,oNDissW,oPO4S
!  varaibles for re-allocation of the roots and shoots
   real(rk)     :: bfRootVeg,bfShootVeg,aDRootVeg,aDShootVeg
   real(rk)     :: aDEmergVeg,aDFloatVeg,bfSubVeg,aDSubVeg
   real(rk)     :: bkMortVeg
   real(rk),save :: aDayInitVeg
!  variables for coverage of vegetation(light function and fish feeding)
   real(rk)     :: afCovEmergVeg,aCovVeg
!  temperature function variables
   real(rk)    ::  uFunTmProdVeg,uFunTmRespVeg
!  light funtion variables
   real(rk)    :: afCovSurfVeg,aLLimShootVeg
   real(rk)    :: uhLVeg,aMuTmLVeg,ufDay,par_bott
   real(rk)    :: aDepth1Veg,aDepth2Veg,aLPAR1Veg,aLPAR2Veg
!  nutrient ratio variables
   real(rk)    :: rPDVeg,rNDVeg
!  variables for dry-weight change of Vegetation
   real(rk)    :: tDBedVeg,tDMigrVeg,tDProdVeg,tDRespVeg
   real(rk)    :: aMuVeg,aNutLimVeg,aPLimVeg,aNLimVeg
   real(rk)    :: tDEnvProdVeg,akDIncrVeg,tDEnvVeg
   real(rk)    :: tDMortVeg,tDEnvMortVeg
!  for vegetation light attenuation
   real(rk)    :: wDBedVeg
!  variables for nitrogen change of Vegetation
   real(rk)    :: tNBedVeg,tNMigrVeg,tNUptVeg
   real(rk)    :: tNUptVegW,tNUptVegS,afNUptVegS,ahNUptVeg
   real(rk)    :: aVNUptMaxCrVeg
   real(rk)    :: aVNUptVegW,aVNUptVegS
   real(rk)    :: tNExcrVeg,tNMortVeg
!  variables for phosphrus change of Vegetation
   real(rk)   :: tPBedVeg,tPMigrVeg,tPUptVeg
   real(rk)   :: tPUptVegW,tPUptVegS,aVPUptMaxCrVeg,ahPUptVeg
   real(rk)   :: aVPUptVegW,aVPUptVegS,afPUptVegS
   real(rk)   :: tPExcrVeg,tPMortVeg
!  variables for external links
!  NH4W & NO3W
   real(rk)   :: wNBedNH4W,tNUptNH4VegW,afNH4UptVegW,tNMortVegNH4W
   real(rk)   :: tNMortVegNH4,wNBedNO3W,tNUptNO3VegW,tNExcrVegS,tNExcrVegW
!  NH4S & NO3S
   real(rk)   :: tNBedNH4S,tNUptNH4VegS,afNH4UptVegS,tNUptNO3VegS
   real(rk)   :: tNMortVegNH4S,tNBedNO3S
!  PO4W&PO4S
   real(rk)   :: wPBedPO4W,tPExcrVegS,tPExcrVegW,tPBedPO4S
   real(rk)   :: tPMortVegPO4W,tPMortVegPO4,tPMortVegPO4S
!  O2W
   real(rk)   :: tO2BedW,tO2ProdVegW,tO2ProdVeg,tO2ProdVegS
   real(rk)   :: tO2RespVegS,tO2RespVegW,tO2UptNO3VegW
!  Detritus in water,DW,Nitrogen,Phosphrus
   real(rk)   :: wDBedDetW,tDMortVegW
   real(rk)   :: wNBedDetW,tNMortVegDetW,tNMortVegDet
   real(rk)   :: wPBedDetW,tPMortVegDetW,tPMortVegDet
!  Detritus in sediment,DW,Nitrogen,Phosphrus
   real(rk)   :: tDBedDetS,tDMortVegS
   real(rk)   :: tNBedDetS,tNMortVegS,tNMortVegDetS
   real(rk)   :: tPBedDetS,tPMortVegS,tPMortVegDetS

! 
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _FABM_HORIZONTAL_LOOP_BEGIN_
!-----------------------------------------------------------------------
!     Retrieve current (local) state variable values.
!-----------------------------------------------------------------------
   _GET_HORIZONTAL_(self%id_sDVeg,sDVeg)
   _GET_HORIZONTAL_(self%id_sNVeg,sNVeg)
   _GET_HORIZONTAL_(self%id_sPVeg,sPVeg)
!-----------------------------------------------------------------------
!  Retrieve dependencis value
!-----------------------------------------------------------------------
!  Retrieve state dependencie value
   _GET_(self%id_NH4poolW,sNH4W)
   _GET_(self%id_NO3poolW,sNO3W)
   _GET_(self%id_PO4poolW,sPO4W)
   _GET_(self%id_O2poolW,sO2W)
   _GET_HORIZONTAL_(self%id_NH4poolS,sNH4S)
   _GET_HORIZONTAL_(self%id_NO3poolS,sNO3S)
   _GET_HORIZONTAL_(self%id_PO4poolS,sPO4S)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_(self%id_dz,dz)
!  current fabm0d, benthic retrieve meanpar in the center of water column
!  benthic par retrieve par at the bottom of the water column
!  this is up to December 3rd, 2014
   _GET_(self%id_meanpar,meanpar)
   _GET_GLOBAL_(self%id_Day,Day)
   _GET_(self%id_extc,extc)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
!  retrieve diagnostic denpendency
   _GET_HORIZONTAL_(self%id_afOxySed,afOxySed)
!  auxiliaries for nutrient variables
   oNDissW = sNO3W + sNH4W
   oNDissS = (sNO3S + sNH4S)/self%cDepthS / self%bPorS
   oNH4S=sNH4S/self%cDepthS / self%bPorS
   oNO3S=sNO3S/self%cDepthS / self%bPorS
   oPO4S=sPO4S/self%cDepthS / self%bPorS
!-----------------------------------------------------------------------
!  Current local nutrients ratios in phytoplankton(check the curent state)
!-----------------------------------------------------------------------
!  P/D_ratio_of_Vegetation
   rPDVeg=sPVeg/(sDVeg+NearZero)
!  N/D_ratio_of_Diatom
   rNDVeg = sNVeg /(sDVeg+NearZero)
!-----------------------------------------------------------------------
!  temperature functions of vegetation
!-----------------------------------------------------------------------
!  temperature_function_of_vegetation_production
   uFunTmProdVeg =  uFunTmVeg(uTm,self%cQ10ProdVeg )
!  temperature_function_of_vegetation_respiration
   uFunTmRespVeg = uFunTmVeg(uTm,self%cQ10RespVeg)
!-----------------------------------------------------------------------
!  the germination, allocation and reallocation process 
!-----------------------------------------------------------------------
!    Initial_growth_only_once_a_year
      if (Day < 1) then
         aDayInitVeg=367
      else if (uTm >= self%cTmInitVeg .and. aDayInitVeg > 366) then 
         aDayInitVeg = Day
      else
         aDayInitVeg=aDayInitVeg
      endif
      

!    setting_root_fration
      if (Day < aDayInitVeg) then
      bfRootVeg = self%fRootVegWin 
      else if (Day < aDayInitVeg + self%cLengAllo) then
      bfRootVeg = 0.5*(self%fRootVegWin + self%fRootVegSum) + 0.5*(self%fRootVegWin - self%fRootVegSum) * &
      &cos(Pi/self%cLengAllo * (Day - aDayInitVeg))
      else if (Day < self%cDayWinVeg) then
      bfRootVeg = self%fRootVegSum 
      else if (Day < self%cDayWinVeg + self%cLengAllo) then
      bfRootVeg = 0.5*(self%fRootVegWin + self%fRootVegSum) - 0.5*(self%fRootVegWin - self%fRootVegSum) * &
      &cos(Pi/self%cLengAllo * (Day - self%cDayWinVeg))
      else
      bfRootVeg = self%fRootVegWin 
      endif
!  mortality_constant
   if (Day < self%cDayWinVeg) then
   bkMortVeg = self%kMortVegSum 
   else if (Day < self%cDayWinVeg + self%cLengMort) then
   bkMortVeg = - log(self%fWinVeg) / self%cLengMort/secs_pr_day
   else
   bkMortVeg = self%kMortVegSum
   endif
!-----------------------------------------------------------------------
!  fractions of roots and shoots
!-----------------------------------------------------------------------
!  shoot_fraction
   bfShootVeg = 1.0_rk - bfRootVeg
!  root_biomass
   aDRootVeg = bfRootVeg * sDVeg
!  shoot_biomass
   aDShootVeg = bfShootVeg * sDVeg
!  emergent_biomass
   aDEmergVeg = self%fEmergVeg * aDShootVeg
!  floating_biomass
   aDFloatVeg = self%fFloatVeg * aDShootVeg
!  submerged_fraction_of_shoot
   bfSubVeg = 1.0_rk - self%fFloatVeg - self%fEmergVeg
!  submerged_biomass
   aDSubVeg = bfSubVeg * aDShootVeg
!-----------------------------------------------------------------------
!  vegetation migration
!-----------------------------------------------------------------------
!  migration_flux
   tDMigrVeg = self%kMigrVeg * (self%cDVegIn - sDVeg)
!  net_migration_flux
   tPMigrVeg = self%kMigrVeg * (self%cPDVeg0*self%cDVegIn - sPVeg)
!  net_migration_flux
   tNMigrVeg = self%kMigrVeg * (self%cNDVeg0*self%cDVegIn - sNVeg)
!-----------------------------------------------------------------------
!  water coverage by vegetation
!-----------------------------------------------------------------------
!  fraction_of_water_SURFACE_covered_by_plant_species
   afCovSurfVeg = min(1.0_rk, max(aDFloatVeg / (self%cDLayerVeg + NearZero), aDEmergVeg / (&
   &self%fEmergVeg * self%cDCarrVeg + NearZero) ))
!  fraction_emergent_coverage
   afCovEmergVeg = min(1.0_rk, 0.01_rk * self%cCovSpVeg * aDEmergVeg)
!  percent_cover
   aCovVeg = min(100.0_rk, self%cCovSpVeg * aDShootVeg)
!=======================================================================
!  the Assimilation part!
!=======================================================================
!-----------------------------------------------------------------------
!  Light function of vegetation  , standard PCLake method
!-----------------------------------------------------------------------
!   feh
!   introduced plant height to subtitute the aDpeht1Veg,if fDepth1Veg=0.5, then plant height
!   is sDepthW-sDepthW*fDepth1Veg
!   Since the bottom light can be retrieved, then light on top of the plant can be calculated.
!  half-sat._light_for_vegetation_production_at_current_temp.
   uhLVeg = self%hLRefVeg * uFunTmProdVeg
   par_bott= meanpar
   aLPAR2Veg = par_bott*exp(- extc *dz/2.0_rk)
   aLPAR1Veg = aLPAR2Veg / exp(- extc * self%cHeightVeg)
   aLLimShootVeg = self%fEmergVeg + self%fFloatVeg * (1.0 - afCovEmergVeg) + bfSubVeg * (1.0 &
   &- afCovSurfVeg) * 1.0 / (extc * sDepthW) * log( (1.0 + aLPAR1Veg / uhLVeg) /& 
   & (1.0 + aLPAR2Veg / uhLVeg))
!=======================================================================
   ufDay = 0.5_rk - 0.2_rk * cos(2.0_rk*Pi*Day / 365.0_rk)
!  max._growth_rate_at_current_temp._AND_light
   aMuTmLVeg =ufDay * bfShootVeg * aLLimShootVeg * uFunTmProdVeg * self%cMuMaxVeg
!-----------------------------------------------------------------------
!  Nutrient limitation functions
!-----------------------------------------------------------------------
!  Droop_function_(P)_for_vegetation
   aPLimVeg = max(0.0_rk, (1.0_rk - self%cPDVegMin / rPDVeg) * self%cPDVegMax / (self%cPDVegMax - self%cPDVegMin) )
!  Droop_function_(N)_for_vegetation
   aNLimVeg = max(0.0_rk, (1.0_rk - self%cNDVegMin / rNDVeg) * self%cNDVegMax / (self%cNDVegMax - self%cNDVegMin) )
!  nutrient_limitation_function_of_vegetation
   aNutLimVeg = min( aPLimVeg,aNLimVeg)
!  actual_growth_rate_of_vegetation
   aMuVeg = aMuTmLVeg * aNutLimVeg
!-----------------------------------------------------------------------
!  vegetation growth rate adjust and correction 
!-----------------------------------------------------------------------
!  intrinsic_net_increase_rate_of_vegetation
   akDIncrVeg = aMuTmLVeg - self%kDRespVeg * uFunTmRespVeg -bkMortVeg
!  logistic_correction_of_vegetation
   tDEnvVeg = max(0.0_rk, sDVeg**2.0_rk*akDIncrVeg / (self%cDCarrVeg+NearZero) )
!  logistic_correction_of_production
   tDEnvProdVeg = aMuVeg / self%cMuMaxVeg * tDEnvVeg
!  vegetation_production
   tDProdVeg = max(0.0_rk, aMuVeg * sDVeg - tDEnvProdVeg)
!-----------------------------------------------------------------------
!  vegetation nutrient uptake ---Nitrogen
!-----------------------------------------------------------------------
!  fraction_of_N_uptake_from_sediment
   if (self%UseEmpUpt==0) then
       afNUptVegS = 0.0_rk
   elseif (bfRootVeg <= NearZero) then
      afNUptVegS = 0.0_rk
   else if (self%fFloatVeg + bfSubVeg <= NearZero) then
      afNUptVegS = 1.0_rk
   else
      afNUptVegS = self%fSedUptVegMax / (1.0_rk + self%fSedUptVegCoef * ((((oNDissS+NearZero) / (oN&
      &DissW+NearZero)) )** self%fSedUptVegExp))
   endif
!  fraction_of_P_uptake_from_sediment
   if (self%UseEmpUpt==0) then
      afPUptVegS = 0.0_rk
   elseif (bfRootVeg <= NearZero) then
      afPUptVegS = 0.0_rk
   else if (self%fFloatVeg + bfSubVeg <= NearZero) then
      afPUptVegS = 1.0_rk
   else
      afPUptVegS = self%fSedUptVegMax / (1.0_rk + self%fSedUptVegCoef * ((((oPO4S+NearZero) / (sPO4&
      &W+NearZero)) )** self%fSedUptVegExp)) 
   endif
!  maximum_P_uptake_rate_of_vegetation,_corrected_for_P/D_ratio
   aVPUptMaxCrVeg = max( 0.0_rk, self%cVPUptMaxVeg * uFunTmProdVeg * (self%cPDVegMax-rPDVeg) / (&
   &self%cPDVegMax-self%cPDVegMin))
!    P_uptake_RATE_by_subm_AND_floating_parts
   if (self%UseEmpUpt==0) then
      aVPUptVegW = sPO4W * aVPUptMaxCrVeg / (aVPUptMaxCrVeg / self%cAffPUptVeg + sPO4W) 
   else
      aVPUptVegW = 0.0_rk 
   endif
!    P_uptake_rate_by_roots
   if  (self%UseEmpUpt==0) then
      aVPUptVegS = oPO4S * aVPUptMaxCrVeg / (aVPUptMaxCrVeg / self%cAffPUptVeg + oPO4S) 
   else
      aVPUptVegS = 0.0_rk
   endif
!  P_uptake_from_water
      if (self%UseEmpUpt==0) then
      tPUptVegW = aVPUptVegW * (aDSubVeg + aDFloatVeg) 
      else
      tPUptVegW = (1.0 - afPUptVegS) * aVPUptMaxCrVeg * sPO4W / (aVPUptMaxCrVeg / self%cAff&
      &PUptVeg + sPO4W) * sDVeg 
      endif
!    P_uptake_from_pore_water_(by_root_fraction)
      if (self%UseEmpUpt==0) then
      tPUptVegS = aVPUptVegS * aDRootVeg 
      else
      tPUptVegS = afPUptVegS * aVPUptMaxCrVeg * oPO4S / (aVPUptMaxCrVeg / self%cAffPUptVeg &
      &+ oPO4S) * sDVeg 
      endif

!    total_P_uptake_vegetation
      tPUptVeg = tPUptVegW + tPUptVegS 

!  maximum_N_uptake_rate_of_vegetation,_corrected_for_N/D_ratio
   aVNUptMaxCrVeg = max( 0.0_rk, self%cVNUptMaxVeg * uFunTmProdVeg * (self%cNDVegMax - rNDVeg) /&
   & (self%cNDVegMax - self%cNDVegMin))
!  half-sat._'constant'_for_N_uptake
   ahNUptVeg = aVNUptMaxCrVeg / self%cAffNUptVeg
!  N_uptake_RATE_by_subm_AND_floating_parts
   if (self%UseEmpUpt==0) then
      aVNUptVegW = oNDissW * aVNUptMaxCrVeg / (ahNUptVeg + oNDissW) 
   else
      aVNUptVegW = 0.0_rk
   endif
!  N_uptake_from_water_(by_shoots)
   if (self%UseEmpUpt==0) then
      tNUptVegW = aVNUptVegW * (aDSubVeg + aDFloatVeg) 
   else
      tNUptVegW = (1.0_rk - afNUptVegS) * aVNUptMaxCrVeg * oNDissW / (aVNUptMaxCrVeg / self%cA&
      &ffNUptVeg + oNDissW) * sDVeg 
   endif
!  N_uptake_RATE_of_roots
   if (self%UseEmpUpt==0) then
      aVNUptVegS = oNDissS * aVNUptMaxCrVeg / (ahNUptVeg + oNDissS) 
   else
      aVNUptVegS = 0.0 
   endif
!  N_uptake_from_pore_water_(by_roots)
   if (self%UseEmpUpt==0)  then
      tNUptVegS = aVNUptVegS * aDRootVeg 
   else
      tNUptVegS = afNUptVegS * aVNUptMaxCrVeg * oNDissS / (aVNUptMaxCrVeg / self%cAffNUptVeg&
      & + oNDissS) * sDVeg 
   endif
!  total Nitrogen uptake of water column and sediment
   tNUptVeg = tNUptVegW+ tNUptVegS

!-----------------------------------------------------------------------
!  vegetation nutrient uptake ---Phosphrus
!-----------------------------------------------------------------------
!  half-sat._'constant'_for_P_uptake
   ahPUptVeg = aVPUptMaxCrVeg / self%cAffPUptVeg
!  P_uptake_RATE_by_subm._AND_floating_parts
   aVPUptVegW = sPO4W * aVPUptMaxCrVeg / (ahPUptVeg + sPO4W)
!  P_uptake_rate_by_roots
   aVPUptVegS = oPO4S * aVPUptMaxCrVeg / (ahPUptVeg + oPO4S) 
!=======================================================================
!  the Dissimilation part!
!=======================================================================
!-----------------------------------------------------------------------
!  vegetation respiration
!-----------------------------------------------------------------------
!  dark_respiration_of_vegetation
   tDRespVeg = self%kDRespVeg * uFunTmRespVeg * sDVeg
!-----------------------------------------------------------------------
!  vegetation excretion,Nitroge
!-----------------------------------------------------------------------
!  N_excretion_by_vegetation
!   tNExcrVeg = rNDVeg / (self%cNDVegMin + rNDVeg) * rNDVeg * tDRespVeg, v5.09
   tNExcrVeg = (2.0_rk * rNDVeg) / (self%cNDVegMax + rNDVeg) * rNDVeg * tDRespVeg ! pl613
!-----------------------------------------------------------------------
!  vegetation excretion,Phosphrus
!-----------------------------------------------------------------------
!  P_excretion_by_vegetation
   tPExcrVeg = rPDVeg / (self%cPDVegMin + rPDVeg) * rPDVeg * tDRespVeg
!-----------------------------------------------------------------------
!  vegetation mortality,Dry-weight
!-----------------------------------------------------------------------
!  logistic_correction_of_mortality
   tDEnvMortVeg = tDEnvVeg - tDEnvProdVeg
!  total_mortality_flux_DW_vegetation
   tDMortVeg = bkMortVeg * sDVeg + tDEnvMortVeg
!-----------------------------------------------------------------------
!  vegetation Mortality,Nitrogen
!-----------------------------------------------------------------------
!  N_mortality_flux_of_vegetation
   tNMortVeg = rNDVeg * tDMortVeg
!-----------------------------------------------------------------------
!  vegetation Mortality,Phosphrus
!-----------------------------------------------------------------------
!  P_mortality_flux_of_vegetation
   tPMortVeg = rPDVeg * tDMortVeg
!-----------------------------------------------------------------------
!  derivative_of_vegetation,the bird grazing part and management part is left to discussion
!-----------------------------------------------------------------------
!  derivative_of_vegetation_biomass
   tDBedVeg = tDMigrVeg + tDProdVeg- tDMortVeg - tDRespVeg
!  total_vegetation_N_flux_in_bed_module
   tNBedVeg = tNMigrVeg + tNUptVeg - tNExcrVeg - tNMortVeg 
!  total_vegetation_P_flux_in_bed_module
   tPBedVeg = tPMigrVeg + tPUptVeg - tPExcrVeg - tPMortVeg
!=======================================================================
!  vegetation part relating to other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in water
!-----------------------------------------------------------------------
!  fraction_ammonium_uptake_from_water_column_(from_WASP_model,_EPA)
   afNH4UptVegW = sNH4W * sNO3W / ((ahNUptVeg + sNH4W) * (ahNUptVeg + sNO3W + NearZero))&
   & + sNH4W * ahNUptVeg / ((sNH4W + sNO3W + NearZero) * (ahNUptVeg + sNO3W + NearZero))
!  NH4_uptake_of_vegetation_from_water
   tNUptNH4VegW = afNH4UptVegW * tNUptVegW
!  N_excretion_by_vegetation_to_sediment
   tNExcrVegS = bfRootVeg * tNExcrVeg
!  N_excretion_by_vegetation_to_water
   tNExcrVegW = tNExcrVeg - tNExcrVegS
!  mortality_flux_of_vegetation_becoming_dissolved_N
   tNMortVegNH4 = self%fDissMortVeg * tNMortVeg
!  mortality_flux_of_vegetation_becoming_dissolved_N_in_sediment
   tNMortVegNH4S = bfRootVeg * tNMortVegNH4
!  mortality_flux_of_vegetation_becoming_dissolved_N_in_water
   tNMortVegNH4W = tNMortVegNH4 - tNMortVegNH4S
!  total_N_flux_from_Vegetation_module_to_NH4_in_water
   wNBedNH4W = - tNUptNH4VegW + tNExcrVegW + tNMortVegNH4W
!-----------------------------------------------------------------------
!  Update NH4 in sediment,still under discussion, Jan 17th, 2014
!-----------------------------------------------------------------------
!  fraction_ammonium_uptake_from_pore_water_(from_WASP_model,_EPA)
   afNH4UptVegS = oNH4S * oNO3S / ((ahNUptVeg + oNH4S +NearZero) * (ahNUptVeg + oNO&
   &3S +NearZero)) + oNH4S * ahNUptVeg / ((oNH4S + oNO3S+NearZero) * (ahNUptVeg + oN&
   &O3S+NearZero)) 
!  NH4_uptake_of_vegetation_from_sediment
   tNUptNH4VegS = afNH4UptVegS * tNUptVegS
!  total_N_flux_from_Vegetation_module_to_NH4_in_pore_water
   tNBedNH4S = - tNUptNH4VegS + tNExcrVegS + tNMortVegNH4S
!-----------------------------------------------------------------------
!  Update NO3 in water
!-----------------------------------------------------------------------
!  NO3_uptake_of_vegetation_from_water
   tNUptNO3VegW = tNUptVegW - tNUptNH4VegW
!  total_N_flux_from_Vegetation_module_to_NO3_in_water
   wNBedNO3W = - tNUptNO3VegW
!-----------------------------------------------------------------------
!  Update NO3 in sediment
!-----------------------------------------------------------------------
! NO3_uptake_of_vegetation_from_sediment
   tNUptNO3VegS = tNUptVegS - tNUptNH4VegS 
! total_N_flux_from_Vegetation_module_to_NO3_in_pore_water
   tNBedNO3S = - tNUptNO3VegS
!-----------------------------------------------------------------------
!  Update PO4 in water
!-----------------------------------------------------------------------
!  P_excretion_by_vegetation_in_sediment
   tPExcrVegS = bfRootVeg * tPExcrVeg
!  P_excretion_by_vegetation_in_water
   tPExcrVegW = tPExcrVeg - tPExcrVegS
!  mortality_flux_of_vegetation_becoming_dissolved_P
   tPMortVegPO4 = self%fDissMortVeg * tPMortVeg
!  mortality_flux_of_vegetation_becoming_dissolved_P_in_sediment
   tPMortVegPO4S = bfRootVeg * tPMortVegPO4
!  mortality_flux_of_vegetation_becoming_dissolved_P_in_water
   tPMortVegPO4W = tPMortVegPO4 - tPMortVegPO4S
!  total_P_flux_from_Vegetation_module_to_PO4_in_water
   wPBedPO4W = - tPUptVegW + tPExcrVegW + tPMortVegPO4W
!-----------------------------------------------------------------------
!  Update PO4 in sediment
!-----------------------------------------------------------------------
! total_P_flux_from_Vegetation_module_to_pore_water_PO4
   tPBedPO4S = - tPUptVegS + tPExcrVegS + tPMortVegPO4S
!-----------------------------------------------------------------------
!  Update O2 in water
!-----------------------------------------------------------------------
!  O2_production_to_water_due_to_NO3_uptake_by_macrophytes
   tO2UptNO3VegW = O2PerNO3 * molO2molN * bfSubVeg * tNUptNO3VegW
!  correction_of_O2_demand_in_water_at_low_oxygen_conc.
   aCorO2BOD = sO2W / (self%hO2BOD + sO2W)
!  submerged_O2_respiration
   tO2RespVegW = molO2molC * self%cCPerDW * bfSubVeg * tDRespVeg * aCorO2BOD
!  root_O2_respiration
   tO2RespVegS = molO2molC * self%cCPerDW * bfRootVeg * tDRespVeg * afOxySed
!   for first time step output
!   tO2RespVegS = molO2molC * self%cCPerDW * bfRootVeg * tDRespVeg * 0.323292
!  vegetation_O2_production
   tO2ProdVeg = molO2molC * self%cCPerDW * tDProdVeg
!  O2_transport_to_roots
   tO2ProdVegS = min (tO2RespVegS,tO2ProdVeg)
!  O2_used_for_vegetation_production
   tO2ProdVegW = min( tO2ProdVeg - tO2ProdVegS, bfSubVeg * tO2ProdVeg)
!  total_water_O2_flux_in_vegetation_module
   tO2BedW = tO2ProdVegW - tO2RespVegW + tO2UptNO3VegW

!-----------------------------------------------------------------------
!  Update detritus  in water, DW, N and P
!-----------------------------------------------------------------------
!  mortality_flux_becoming_water_detritus
   tDMortVegW = self%fDetWMortVeg * (1.0_rk - bfRootVeg) * tDMortVeg
!  total_DW_flux_from_Vegetation_module_to_water_detritus
   wDBedDetW = tDMortVegW
!  mortality_flux_of_vegetation_becoming_detritus_N
   tNMortVegDet = tNMortVeg - tNMortVegNH4
!  mortality_flux_of_vegetation_becoming_detritus_N_in_water
   tNMortVegDetW = self%fDetWMortVeg * (1.0_rk - bfRootVeg) * tNMortVegDet
!  total_N_flux_from_Vegetation_module_to_water_detritus
   wNBedDetW = tNMortVegDetW
!  mortality_flux_of_vegetation_becoming_detritus_P
   tPMortVegDet = tPMortVeg - tPMortVegPO4
!  mortality_flux_of_vegetation_becoming_detritus_P_in_water
   tPMortVegDetW = self%fDetWMortVeg * (1.0_rk - bfRootVeg) * tPMortVegDet
!  total_P_flux_from_Vegetation_module_to_water_detritus
   wPBedDetW = tPMortVegDetW
!---------------------------------------------------------------------------
!  Update detritus  in sediment, DW, N and P
!---------------------------------------------------------------------------
!  mortality_flux_becoming_sediment_detritus
   tDMortVegS = tDMortVeg - tDMortVegW
!  total_DW_flux_from_Vegetation_module_to_sediment_detritus
   tDBedDetS = tDMortVegS
!  mortality_flux_of_vegetation_becoming_detritus_N_in_sediment
   tNMortVegDetS = tNMortVegDet - tNMortVegDetW
!  total_N_flux_from_Vegetation_module_to_sediment_detritus
   tNBedDetS = tNMortVegDetS
!  mortality_flux_of_vegetation_becoming_detritus_P_in_sediment
   tPMortVegDetS = tPMortVegDet - tPMortVegDetW
!  total_P_flux_from_Vegetation_module_to_sediment_detritus
   tPBedDetS = tPMortVegDetS
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_sDVeg,tDBedVeg)
   _SET_ODE_BEN_(self%id_sNVeg,tNBedVeg)
   _SET_ODE_BEN_(self%id_sPVeg,tPBedVeg)
!-----------------------------------------------------------------------
!  Output local diagnostic variables
!-----------------------------------------------------------------------
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tO2BedW,tO2BedW)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
   _SET_BOTTOM_EXCHANGE_(self%id_NH4poolW,wNBedNH4W)
   _SET_BOTTOM_EXCHANGE_(self%id_NO3poolW,wNBedNO3W)
   _SET_BOTTOM_EXCHANGE_(self%id_PO4poolW,wPBedPO4W)
   _SET_BOTTOM_EXCHANGE_(self%id_O2poolW,tO2BedW)
   _SET_BOTTOM_EXCHANGE_(self%id_DDetpoolW,wDBedDetW)
   _SET_BOTTOM_EXCHANGE_(self%id_DNetpoolW,wNBedDetW)
   _SET_BOTTOM_EXCHANGE_(self%id_DPetpoolW,wPBedDetW)
!  in the sediment
   _SET_ODE_BEN_(self%id_NH4poolS,tNBedNH4S)
   _SET_ODE_BEN_(self%id_NO3poolS,tNBedNO3S)
   _SET_ODE_BEN_(self%id_PO4poolS,tPBedPO4S)
   _SET_ODE_BEN_(self%id_DDetpoolS,tDBedDetS)
   _SET_ODE_BEN_(self%id_DNetpoolS,tNBedDetS)
   _SET_ODE_BEN_(self%id_DPetpoolS,tPBedDetS)
!-----------------------------------------------------------------------
!  output diagnostic variables for external links
!-----------------------------------------------------------------------
!  Export diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aDSubVeg,aDSubVeg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aCovVeg,aCovVeg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBedDetS,tDBedDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_afCovSurfVeg,afCovSurfVeg)
   
!  for vegetation light attenutaion output
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aDayInitVeg,aDayInitVeg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBedVeg,tDBedVeg)
   
!  light output
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_macroextinction,extc)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aLPAR1Veg,aLPAR1Veg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aLPAR2Veg,aLPAR2Veg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_allimveg,aLLimShootVeg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aNutLimVeg,aNutLimVeg)

   _FABM_HORIZONTAL_LOOP_END_
!
! Spatial loop end
!
!EOP
!-----------------------------------------------------------------------

    end subroutine do_bottom

!EOC
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
   class (type_au_pclake_macrophytes), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk) :: sDVeg,uTm,sDepthW
   real(rk) :: Day
   real(rk),save :: aDayInitVeg
   real(rk) :: bfRootVeg,bfShootVeg,aDShootVeg
   real(rk) :: bfSubVeg,aDSubVeg
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_HORIZONTAL_(self%id_sDVeg,sDVeg)
   _GET_GLOBAL_(self%id_Day,Day)
   _GET_(self%id_uTm,uTm)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)

!  Self-shading with explicit contribution from background phytoplankton concentration.
!    Initial_growth_only_once_a_year
   if (Day < 1) then
      aDayInitVeg=367
   else if (uTm >= self%cTmInitVeg .and. aDayInitVeg > 366) then 
      aDayInitVeg = Day
   else
      aDayInitVeg=aDayInitVeg
   endif
 
! setting_root_fration
   if (Day < aDayInitVeg) then
   bfRootVeg = self%fRootVegWin 
   else if (Day < aDayInitVeg + self%cLengAllo) then
   bfRootVeg = 0.5*(self%fRootVegWin + self%fRootVegSum) + 0.5*(self%fRootVegWin - self%fRootVegSum) * &
   &cos(Pi/self%cLengAllo * (Day - aDayInitVeg))
   else if (Day < self%cDayWinVeg) then
   bfRootVeg = self%fRootVegSum 
   else if (Day < self%cDayWinVeg + self%cLengAllo) then
   bfRootVeg = 0.5*(self%fRootVegWin + self%fRootVegSum) - 0.5*(self%fRootVegWin - self%fRootVegSum) * &
   &cos(Pi/self%cLengAllo * (Day - self%cDayWinVeg)) 
   else
   bfRootVeg = self%fRootVegWin 
   endif
!-----------------------------------------------------------------------
!  fractions of roots and shoots
!-----------------------------------------------------------------------
!  shoot_fraction
   bfShootVeg = 1.0_rk - bfRootVeg
!  shoot_biomass
   aDShootVeg = bfShootVeg * sDVeg
!  submerged_fraction_of_shoot
   bfSubVeg = 1.0_rk - self%fFloatVeg - self%fEmergVeg
!  submerged_biomass
   aDSubVeg = bfSubVeg * aDShootVeg
!  convert unit from area loading to concentration
   aDSubVeg=aDSubVeg/sDepthW

   _SET_EXTINCTION_(self%cExtSpVeg*aDSubVeg)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------
!
!EOP
!-----------------------------------------------------------------------

   end module au_pclake_macrophytes

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

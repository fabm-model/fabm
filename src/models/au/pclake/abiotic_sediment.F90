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
!  Each state variable and its related its local processes are:
!  Inorganic matter: sDIMS, processes: none local processes
!  Organic mattter: sDDetS,sNDetS,sPDetS, processes: mineralization
!  Dissolved nutrients: sNH4S,process:mineralisation,nitrification,
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
!  This module also provide important diagnostic variable will be used 
!  in other modules, including:
!  Sediment aerobic layer fraction, afOxySed, used by module: 
!                                              macrophytes module
!  Sediment detritus change, tDAbioDetS, used bymodule:auxilary
!
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
!     tDAbioO2S: abiotic sediment oxygen consumption
!     aPEqIMS: equilibrium absorped phosphrus concentration
      type (type_horizontal_diagnostic_variable_id) :: id_rPDDetS,id_rNDDetS
      type (type_horizontal_diagnostic_variable_id) :: id_tDAbioO2S
      type (type_horizontal_diagnostic_variable_id) :: id_aPEqIMS
      
!     state dependencies identifers
!     MinSiO2Sed: Mineralization generated SiO2 from sediment
!     O2ConsumpSed: O2 consumption in sediment
!     diff+nut: diffusion fluxes of nutrients between water and sediment
      type (type_state_variable_id) :: id_MinSiO2Sed ,id_O2ConsumpSed
      type (type_state_variable_id) :: id_diffNH4,id_diffNO3,id_diffPO4
!     diagnostic variables for dependencies(without output)
      type (type_horizontal_diagnostic_variable_id) :: id_tDAbioDetS,id_afOxySed,id_tDAbioHumS
!     environmental dependencies
      type (type_dependency_id)                :: id_uTm
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
!     detritus related paramters
      real(rk)                   :: fRefrDetS,cThetaMinS,kDMinDetS
      real(rk)                   :: kNitrS,cThetaNitr
!     Humus related paramters
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


!  private data memebers(API0.92)
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
!
! !local variables
   real(rk)                  :: sNH4S_initial
   real(rk)                  :: sNO3S_initial
   real(rk)                  :: sPO4S_initial
   real(rk)                  :: sPAIMS_initial
   real(rk)                  :: sDIMS_initial
   real(rk)                  :: sDDetS_initial
   real(rk)                  :: sPDetS_initial
   real(rk)                  :: sNDetS_initial
   real(rk)                  :: sSiDetS_initial
   real(rk)                  :: sDHumS_initial
   real(rk)                  :: sPHumS_initial
   real(rk)                  :: sNHumS_initial
   real(rk)                  :: cDepthS
   real(rk)                  :: fRefrDetS
   real(rk)                  :: kDMinDetS
   real(rk)                  :: cThetaMinS
   real(rk)                  :: cCPerDW
   real(rk)                  :: O2PerNH4
   real(rk)                  :: kNitrS
   real(rk)                  :: cThetaNitr
   real(rk)                  :: NO3PerC
   real(rk)                  :: hNO3Denit
   real(rk)                  :: kPSorp
   real(rk)                  :: cRelPAdsD
   real(rk)                  :: cRelPAdsFe
   real(rk)                  :: fFeDIM
   real(rk)                  :: cRelPAdsAl
   real(rk)                  :: fAlDIM
   real(rk)                  :: fRedMax
   real(rk)                  :: cKPAdsOx
   real(rk)                  :: kPChemPO4
   real(rk)                  :: coPO4Max
   real(rk)                  :: bPorS
   real(rk)                  :: cThetaDif
   real(rk)                  :: fDepthDifS
   real(rk)                  :: kNDifNH4
   real(rk)                  :: cTurbDifNut
   real(rk)                  :: bPorCorS
   real(rk)                  :: kNDifNO3
   real(rk)                  :: kPDifPO4
   real(rk)                  :: kO2Dif
   real(rk)                  :: cTurbDifO2
   real(rk)                  :: kDMinHum
   character(len=64)         :: oxygen_pool_water
   character(len=64)         :: SiO2_generated_by_mineralization
   character(len=64)         :: NH4_diffusion_flux
   character(len=64)         :: NO3_diffusion_flux
   character(len=64)         :: PO4_diffusion_flux
!  create namelist
   namelist /au_pclake_abiotic_sediment/ sDHumS_initial,sNHumS_initial,sPHumS_initial, &
                      & sNH4S_initial,sNO3S_initial,sPO4S_initial,sPAIMS_initial,sDIMS_initial, &
                      & sDDetS_initial,sNDetS_initial,sPDetS_initial,sSiDetS_initial,&
                      & cDepthS,fRefrDetS,kDMinDetS, &
                      & cThetaMinS,cCPerDW,O2PerNH4,kNitrS,cThetaNitr,NO3PerC,hNO3Denit,kPSorp,&
                      & cRelPAdsD,cRelPAdsFe,fFeDIM,cRelPAdsAl,fAlDIM,fRedMax,cKPAdsOx,kPChemPO4,coPO4Max,bPorS,&
                      & cThetaDif,fDepthDifS,kNDifNH4,cTurbDifNut,bPorCorS, kNDifNO3,kPDifPO4,kO2Dif,cTurbDifO2,kDMinHum,&
                      & oxygen_pool_water,SiO2_generated_by_mineralization,NH4_diffusion_flux,NO3_diffusion_flux, &
                      & PO4_diffusion_flux
!EOP
!-----------------------------------------------------------------------
!BOC
!
!  initialize the parameters
   sNH4S_initial=0.02_rk
   sNO3S_initial=0.002_rk
   sPO4S_initial=0.182_rk
   sPAIMS_initial=17.99_rk
   sDIMS_initial=39611.3_rk
   sDDetS_initial=181.7_rk
   sNDetS_initial=4.54_rk
   sPDetS_initial=0.454_rk
   sSiDetS_initial=1.82_rk
   sDHumS_initial=3452.34_rk
   sNHumS_initial=172.62_rk
   sPHumS_initial=17.26_rk
   cDepthS=0.1_rk
   fRefrDetS=0.15_rk
   kDMinDetS=0.002_rk
   cThetaMinS=1.07_rk
   cCPerDW=0.4_rk
   O2PerNH4=2.0_rk
   kNitrS=1.0_rk
   cThetaNitr=1.08_rk
   NO3PerC=0.8_rk
   hNO3Denit=2.0_rk
   kPSorp=0.05_rk
   cRelPAdsD=0.00003_rk
   cRelPAdsFe=0.065_rk
   fFeDIM=0.01_rk
   cRelPAdsAl=0.134_rk
   fAlDIM=0.01_rk
   fRedMax=0.9_rk
   cKPAdsOx=0.6_rk
   kPChemPO4=0.03_rk
   coPO4Max=2.0_rk
   bPorS= 0.85_rk
   cThetaDif = 1.02_rk
   fDepthDifS = 0.5_rk
   kNDifNH4 = 0.000112_rk
   cTurbDifNut =5.0_rk
   bPorCorS = 0.975_rk
   kNDifNO3 = 0.000086_rk
   kPDifPO4=0.000072_rk
   kO2Dif = 0.000026_rk
   cTurbDifO2 = 5.0_rk
   kDMinHum=0.00001_rk
   oxygen_pool_water=''
   SiO2_generated_by_mineralization=''
   NH4_diffusion_flux=''
   NO3_diffusion_flux=''
   PO4_diffusion_flux=''
!
!  Read parameters namelist
   if (configunit>0) read(configunit,nml=au_pclake_abiotic_sediment,err=99,end=100)
!
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cDepthS,'cDepthS',default=cDepthS)
   call self%get_parameter(self%fRefrDetS,'fRefrDetS',default=fRefrDetS)
   call self%get_parameter(self%kDMinDetS,'kDMinDetS',default=kDMinDetS,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cThetaMinS,'cThetaMinS',default=cThetaMinS)
   call self%get_parameter(self%cCPerDW,'cCPerDW',default=cCPerDW)
   call self%get_parameter(self%O2PerNH4,'O2PerNH4',default=O2PerNH4)
   call self%get_parameter(self%kNitrS,'kNitrS',default=kNitrS,scale_factor =1.0_rk/secs_pr_day)
   call self%get_parameter(self%cThetaNitr,'cThetaNitr',default=cThetaNitr)
   call self%get_parameter(self%NO3PerC,'NO3PerC',default=NO3PerC)
   call self%get_parameter(self%hNO3Denit,'hNO3Denit',default=hNO3Denit)
   call self%get_parameter(self%kPSorp,'kPSorp',default=kPSorp,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cRelPAdsD,'cRelPAdsD',default=cRelPAdsD)
   call self%get_parameter(self%cRelPAdsFe,'cRelPAdsFe',default=cRelPAdsFe)
   call self%get_parameter(self%fFeDIM,'fFeDIM',default=fFeDIM)
   call self%get_parameter(self%cRelPAdsAl,'cRelPAdsAl',default=cRelPAdsAl)
   call self%get_parameter(self%fAlDIM,'fAlDIM',default=fAlDIM)
   call self%get_parameter(self%fRedMax,'fRedMax',default=fRedMax)
   call self%get_parameter(self%cKPAdsOx,'cKPAdsOx',default=cKPAdsOx)
   call self%get_parameter(self%kPChemPO4,'kPChemPO4',default=kPChemPO4,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%coPO4Max,'coPO4Max',default=coPO4Max)
   call self%get_parameter(self%bPorS,'bPorS',default=bPorS)
   call self%get_parameter(self%cThetaDif,'cThetaDif',default=cThetaDif)
   call self%get_parameter(self%fDepthDifS,'fDepthDifS',default=fDepthDifS)
   call self%get_parameter(self%kNDifNH4,'kNDifNH4',default=kNDifNH4, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cTurbDifNut,'cTurbDifNut',default=cTurbDifNut)
   call self%get_parameter(self%bPorCorS,'bPorCorS',default=bPorCorS)
   call self%get_parameter(self%kNDifNO3,'kNDifNO3',default=kNDifNO3,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPDifPO4,'kPDifPO4',default=kPDifPO4, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kO2Dif,'kO2Dif',default=kO2Dif,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cTurbDifO2,'cTurbDifO2',default=cTurbDifO2)
   call self%get_parameter(self%kDMinHum,'kDMinHum',default=kDMinHum ,scale_factor=1.0_rk/secs_pr_day)
   
!
!  Register local state variable
!
   call self%register_state_variable(self%id_sDIMS,'sDIMS','g/m**2','sediment inorg.Matter',     &
                                     sDIMS_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sDDetS,'sDDetS','g/m**2','sediment detritus DW',     &
                                     sDDetS_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sNDetS,'sNDetS','g/m**2','sediment detritus N',     &
                                     sNDetS_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sPDetS,'sPDetS','g/m**2','sediment detritus P',     &
                                     sPDetS_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sSiDetS,'sSiDetS','g/m**2','sediment detritus Si',     &
                                     sSiDetS_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sPO4S,'sPO4S','g/m**2','Sediment Phosphate',     &
                                     sPO4S_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sPAIMS,'sPAIMS','g/m**2','SED_Absorbed Phosphate',     &
                                     sPAIMS_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sNH4S,'sNH4S','g/m**2','Sediment Amonia',     &
                                     sNH4S_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sNO3S,'sNO3S','g/m**2','Sediment Nitrates',     &
                                     sNO3S_initial,minimum=_ZERO_)
!  Humus
   call self%register_state_variable(self%id_sDHumS,'sDHumS','g/m**2','sediment Humus DW',     &
                                     sDHumS_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sNHumS,'sNHumS','g/m**2','sediment Humus N',     &
                                     sNHumS_initial,minimum=_ZERO_)
   call self%register_state_variable(self%id_sPHumS,'sPHumS','g/m**2','sediment Humus P',     &
                                     sPHumS_initial,minimum=_ZERO_)
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNH4S)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNO3S)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNDetS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPO4S)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPAIMS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPDetS)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,self%id_sSiDetS)

!---------------------------------------------------------------------------------------------------------------
!  regirster state variables dependencies (2 steps) 
!---------------------------------------------------------------------------------------------------------------
!  step1, Register dependencies on external state variables
   call self%register_state_dependency(self%id_O2ConsumpSed, 'oxygen_pool_water','g/m**3','oxygen_pool_water')
   call self%register_state_dependency(self%id_MinSiO2Sed, 'SiO2_generated_by_mineralization','g/m**3','SiO2_generated_by_mineralization')
   call self%register_state_dependency(self%id_diffNH4, 'NH4_diffusion_flux','g/m**3','NH4_diffusion_flux')
   call self%register_state_dependency(self%id_diffNO3, 'NO3_diffusion_flux','g/m**3','NO3_diffusion_flux')
   call self%register_state_dependency(self%id_diffPO4, 'PO4_diffusion_flux','g/m**3','PO4_diffusion_flux')
!  step 2, Automatically couple dependencies if target variables have been specified.
   if (oxygen_pool_water/='') call self%request_coupling(self%id_O2ConsumpSed,oxygen_pool_water)
   if (SiO2_generated_by_mineralization/='') call self%request_coupling(self%id_MinSiO2Sed,SiO2_generated_by_mineralization)
   if (NH4_diffusion_flux/='') call self%request_coupling(self%id_diffNH4,NH4_diffusion_flux)
   if (NO3_diffusion_flux/='') call self%request_coupling(self%id_diffNO3,NO3_diffusion_flux)
   if (PO4_diffusion_flux/='') call self%request_coupling(self%id_diffPO4,PO4_diffusion_flux)!
!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_afOxySed,'afOxySed','-- ','afOxySed',output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_tDAbioDetS,'tDAbioDetS','g/m**2/s ','abiotic_sediment_DDet_change', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_tDAbioHumS,'tDAbioHumS','g/m**2/s ','abiotic_sediment_DHum_change', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_tDAbioO2S,'tDAbioO2S','g/m**2/s ','abiotic_sediment_O2_change', &
                                          output=output_none)
   call self%register_diagnostic_variable(self%id_rPDDetS,'rPDDetS','-- ','detritus_P/D_ration_sed', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rNDDetS,'rNDDetS','-- ','detritus_N/D_ration_sed', &
                                          output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_aPEqIMS,'aPEqIMS','-- ','equilibrium_absorbed_PO4', &
                                          output=output_time_step_averaged)

!  register environmental dependencies
   call self%register_dependency(self%id_uTm,standard_variables%temperature)
   return

99 call self%fatal_error('au_pclake_abiotic_sediment_init','Error reading namelist au_pclake_abiotic_sediment')

100 call self%fatal_error('au_pclake_abiotic_sediment_init','Namelist au_pclake_abiotic_sediment was not found.')
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
   real(rk)                   :: uTm
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
!
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
   _SET_BOTTOM_EXCHANGE_(self%id_diffNH4,tNdifNH4)
   _SET_BOTTOM_EXCHANGE_(self%id_diffNO3,tNdifNO3)
   _SET_BOTTOM_EXCHANGE_(self%id_diffPO4,tPdifPO4)
!  update O2 in water column
   _SET_BOTTOM_EXCHANGE_(self%id_O2ConsumpSed,-tO2MinDetS - tO2NitrS)
!-----------------------------------------------------------------------
!  Output denpendent diagnostic variables for other modules
!-----------------------------------------------------------------------
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_afOxySed,afOxySed)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioDetS,tDAbioDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioO2S,-tO2MinDetS - tO2NitrS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rPDDetS,rPDDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rNDDetS,rNDDetS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioHumS,tDAbioHumS)
   
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aPEqIMS,aPEqIMS)

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

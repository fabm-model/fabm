#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_abiotic_water
!
! !DESCRIPTION:
!
!  The pclake_abiotic_water module describes the state variables which are related to abiotic
!  processes in water column, including: inorganic matter(IM), organic matters(detritus),
!  dissolved nutirents(ammonia, nitrate,phosphate and dissolved sicica dioxide), immobilized phosphrus
!  (absorbed phosphrus), dissolved oxygen. Each state variable and its related processes are listed as
!  followed:
!  Inorganic matter: sDIMW, no local processes
!  Organic matter: sDDetW,sPDetW,sNDetW,sSiDetW, processes:mineralisation
!  Dissolved nutrients: sNH4W,process:mineralisation,nitrification
!  Dissolved nutrients:sNO3W,processes: nitrification,denitrification
!  Dissolved nutrients: sPO4W,processes: mineralisation,phosphrus absorption
!  Dissolved nutrients: sSiO2W, processes: mineralisation
!  Dissolved oxygen: sO2W,processes reaeration,mineralisation oxygen consumption, nitrification oxygen consumption
!  Absorbed_P: sPAIMW, processes:phosphrus absorption
! !USES:
   use fabm_types
   use fabm_standard_variables
   use pclake_utility, ONLY: uFunTmAbio

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_abiotic_water
!  Module description, variables
!     local state variable identifers
!     sDIMW: inorganic matter concentration, in dry-weight, gDW/m**3
!     sDDetW, detritus concentration, in dry-weight, gDW/m**3
!     sNDetW, detritus concentration, in nitrogen element, gN/m**3
!     sPDetW, detritus concentration, in phosphorus element, gP/m**3
!     sNH4W,  ammonia concentration, in nitrogen element,gN/m**3
!     sNO3W,  nitrate concentration, in nitrogen element,gN/m**3
!     sPO4W,  phosphate concentration, in phosphorus element,gP/m**3
!     sPAIMW, absorbed phosphorus concentration, in phosphorus element,gP/m**3
!     sSiO2W, silica dioxide concentration, in silica element,gSi/m**3
!     sO2W,   dissolved oxygen, in O2 molecule,gO2/m**3
      type (type_state_variable_id)   :: id_sDIMW,id_sDDetW,id_sNDetW,id_sPDetW,id_sSiDetW
      type (type_state_variable_id)   :: id_sPO4W,id_sPAIMW,id_sNH4W,id_sNO3W
      type (type_state_variable_id)   :: id_sSiO2W,id_sO2W

!     diagnostic variables for local output
!     rPDDetW: P/D ratio of detritus
!     rNDDetW: N/D ratio of detritus
!     wO2AbioW: abiotic water column oxygen consumption
!     tO2Aer: O2 reaeration rate at the water surface
      type (type_diagnostic_variable_id)           :: id_rPDDetW,id_rNDDetW,id_wO2AbioW
      type (type_diagnostic_variable_id)           :: id_extIM,id_extDet
      type (type_horizontal_diagnostic_variable_id):: id_tO2Aer,id_wind
!     environmental dependencies
      type (type_dependency_id)                :: id_uTm
      type (type_horizontal_dependency_id)     :: id_uVWind
!     Model parameters
      real(rk)                   :: kDMinDetW,cThetaMinW
      real(rk)                   :: kNitrW,cThetaNitr,hO2Nitr
      real(rk)                   :: cThetaAer,cCPerDW,hO2BOD,O2PerNH4
      real(rk)                   :: NO3PerC,hNO3Denit,kPSorp,cRelPAdsD
      real(rk)                   :: cRelPAdsFe,fFeDIM,cRelPAdsAl,fAlDIM
      real(rk)                   :: fRedMax,cKPAdsOx
!     sinking parameter
      real(rk)                   :: cVSetIM,cVSetDet
!     parameter for specific light attenuation coefficient
      real(rk)                   :: cExtSpIM,cExtSpDet

   contains
!     Model procedures
      procedure :: initialize
      procedure :: do
      procedure :: get_light_extinction
      procedure :: do_surface


   end type type_pclake_abiotic_water

!  private data memebers
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
!  ratio of mol.weights, = 32/12 [gO2/gC],
   real(rk),parameter :: molO2molC = 2.6667_rk
!  ratio of mol.weights,32/14 [gO2/gN],
   real(rk),parameter :: molO2molN = 2.2857_rk
!  ratio of mol.weights,14/12 [gN/gC],
   real(rk),parameter :: molNmolC = 1.1667_rk
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the inorganic and organic matter model in water
!
! !INTERFACE:

   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the pclake_abiotic_water namelist is read and the variables
!   are registered with FABM.They should be initialised as
!  concentrations.The detritus P and N variables are initialized according
!  to the initial P/D and N/D ratios. Adsorbed P, usually a minor component
!  in the water column,is initialised at 0 gP/m**3.
!  Adsorption is calculated, however, during the run.
!
! !INPUT PARAMETERS:
   class (type_pclake_abiotic_water), intent(inout), target  :: self
   integer,                          intent(in)            :: configunit

!EOP
!-----------------------------------------------------------------------------
!BOC

!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cCPerDW,   'cCPerDW',   'gC/gDW', 'C content of organic matter',                              default=0.4_rk)
   call self%get_parameter(self%cExtSpDet, 'cExtSpDet', 'm2/gDW', 'specific extinction detritus',                             default=0.15_rk)
   call self%get_parameter(self%cExtSpIM,  'cExtSpIM',  'm2/gDW', 'specific extinction inert matter',                         default=0.05_rk)
   call self%get_parameter(self%cKPAdsOx,  'cKPAdsOx',  'm3/gP',  'P adsorption affinity at oxidized conditions',             default=0.6_rk)
   call self%get_parameter(self%cRelPAdsAl,'cRelPAdsAl','gP/gAl', 'max. P adsorption per g Al',                               default=0.134_rk)
   call self%get_parameter(self%cRelPAdsD, 'cRelPAdsD', 'gP/gD',  'max. P adsorption per g DW',                               default=0.00003_rk)
   call self%get_parameter(self%cRelPAdsFe,'cRelPAdsFe','gP/gFe', 'max. P adsorption per g Fe',                               default=0.065_rk)
   call self%get_parameter(self%cThetaAer, 'cThetaAer', '1/e^oC', 'temperature coeff. for reaeration',                        default=1.024_rk)
   call self%get_parameter(self%cThetaMinW,'cThetaMinW','[-]',    'expon. temp. constant of mineralization in water',         default=1.07_rk)
   call self%get_parameter(self%cThetaNitr,'cThetaNitr','[-]',    'temperature coefficient for nitrification',                default=1.08_rk)
   call self%get_parameter(self%cVSetDet,  'cVSetDet',  'm/day',  'max. sedimentation velocity of detritus',                  default=-0.25_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVSetIM,   'cVSetIM',   'm/day',  'max. sedimentation velocity of inert org. matter',         default=-1.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fAlDIM,    'fAlDIM',    'gAl/gD', 'Al content of inorg. matter',                              default=0.01_rk)
   call self%get_parameter(self%fFeDIM,    'fFeDIM',    'gFe/gD', 'Fe content of inorg. matter',                              default=0.01_rk)
   call self%get_parameter(self%fRedMax,   'fRedMax',   '[-]',    'max. reduction factor of P adsorption affinity',           default=0.9_rk)
   call self%get_parameter(self%hNO3Denit, 'hNO3Denit', 'mgN/l',  'quadratic half-sat. NO3 conc. for denitrification',        default=2.0_rk)
   call self%get_parameter(self%hO2BOD,    'hO2BOD',    'mgO2/l', 'half-sat. oxygen conc. for BOD',                           default=1.0_rk)
   call self%get_parameter(self%hO2Nitr,   'hO2Nitr',   'mgO2/l', 'half sat.oxygen conc. for for nitrogen',                   default=2.0_rk)
   call self%get_parameter(self%kDMinDetW, 'kDMinDetW', 'day-1',  'decomposition constant of detritus',                       default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kNitrW,    'kNitrW',    'day-1',  'nitrification rate constant in water',                     default=0.1_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPSorp,    'kPSorp',    'day-1',  'P sorption rate constant not too high -> model speed day', default=0.05_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%NO3PerC,   'NO3PerC',   'molNO3', 'denitrified per mol C mineralised',                        default=0.8_rk)
   call self%get_parameter(self%O2PerNH4,  'O2PerNH4',  'molO2',  'used per mol NH4+ nitrified',                              default=2.0_rk)



!  Register local state variable
!  particles, including inorganic matter(sDIM) and organic matter(sDDetW,sNDetW,sPDetW,sSiDetW) have
!  vertical movement, usually settling(negative values)
   call self%register_state_variable(self%id_sDIMW,'sDIMW','g m-3','Inorg. matter in water',           &
                                    initial_value=5.0_rk,  minimum=_ZERO_, vertical_movement= self%cVSetIM,no_river_dilution=.FALSE.)
!  detritus
   call self%register_state_variable(self%id_sDDetW,'sDDetW','g m-3','detritus dry-weight in water',    &
                                    initial_value=2.0_rk,  minimum=_ZERO_, vertical_movement= self%cVSetDet,no_river_dilution=.FALSE.)
   call self%register_state_variable(self%id_sPDetW,'sPDetW','g m-3','detritus phosphrus in water',     &
                                    initial_value=0.005_rk,minimum=_ZERO_,vertical_movement= self%cVSetDet,no_river_dilution=.FALSE.)
   call self%register_state_variable(self%id_sNDetW,'sNDetW','g m-3','detritus nitrogen in water',      &
                                    initial_value=0.05_rk ,minimum=_ZERO_,vertical_movement=self%cVSetDet,no_river_dilution=.FALSE.)
   call self%register_state_variable(self%id_sSiDetW,'sSiDetW','g m-3','detritus silica in water',      &
                                    initial_value=0.02_rk, minimum=_ZERO_,vertical_movement= self%cVSetDet,no_river_dilution=.FALSE.)
!  dissolved nutrients
   call self%register_state_variable(self%id_sPO4W,  'sPO4W',   'g m-3','Phosphate in water',     &
                                    initial_value=0.01_rk, minimum=_ZERO_,no_river_dilution=.FALSE.)
   call self%register_state_variable(self%id_sPAIMW,'sPAIMW','g m-3','Absorbed_P in water',     &
                                    initial_value=0.0_rk,minimum=_ZERO_,vertical_movement= self%cVSetIM,no_river_dilution=.FALSE.)
   call self%register_state_variable(self%id_sNH4W,'sNH4W','g m-3','Amonia in water',     &
                                    initial_value=0.1_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNO3W,'sNO3W','g m-3','Nitrates in water',     &
                                    initial_value=0.1_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sSiO2W,'sSiO2W','g m-3','Silica dioxide in water',     &
                                    initial_value=3.0_rk,minimum=_ZERO_,no_river_dilution=.FALSE.)
!  oxygen
   call self%register_state_variable(self%id_sO2W,'sO2W','g m-3','oxygen in water',     &
                                    initial_value=10.0_rk,minimum=_ZERO_,no_river_dilution=.FALSE.)
!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_wO2AbioW,'wO2AbioW',  'g m-3 s-1','abiotic_water_O2_change',   output=output_none)
   call self%register_diagnostic_variable(self%id_tO2Aer,  'tO2Aer',    'g m-2 s-1','O2_reareation',             output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_rPDDetW, 'rPDDetW',   '[-]',       'detritus_P/D_ratio_wat',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDDetW, 'rNDDetW',   '[-]',       'detritus_N/D_ratio_wat',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wind,    'wind',      'm/s',       'windspeed',                output=output_none)
   call self%register_diagnostic_variable(self%id_extIM,   'extIM',     '[-]',       'extIM',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_extDet,  'extDet',    '[-]',       'extDet',                   output=output_instantaneous)


!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNDetW)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNH4W)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNO3W)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPO4W)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPAIMW)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPDetW)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,              self%id_sSiO2W)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,              self%id_sSiDetW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNDetW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNH4W)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNO3W)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPO4W)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPAIMW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPDetW)
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,   standard_variables%temperature)
   call self%register_dependency(self%id_uVWind,standard_variables%wind_speed)

!  register diagnostic dependencies

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
   class (type_pclake_abiotic_water), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_
! !LOCAL VARIABLES:
!  carriers for local state variables
   real(rk)     :: sNH4W,sNO3W,sPO4W,sPAIMW,sSiO2W,sO2W
   real(rk)     :: sDDetW,sNDetW,sPDetW,sDIMW,sSiDetW
!  Nutrients ratio value
   real(rk)     :: rPDDetW,rNDDetW
!  carriers for environmental dependencies
   real(rk)     :: uTm
!  local variables for processes
!  Temperature function variables
   real(rk)     :: uFunTmMinW,uFunTmNitr
!  O2 functions variables
   real(rk)     :: aCorO2BOD,aCorO2NitrW,wO2MinDetW,wO2NitrW
!  mineralisation function variables
   real(rk)     :: kNMinDetW,kPMinDetW,kSiMinDetW
   real(rk)     :: wDMinDetW,wNMinDetW,wPMinDetW,wSiMinDetW
!  sorption funtion variables
   real(rk)     :: wPSorpIMW,aPEqIMW,aPIsoAdsW,aPAdsMaxW,aKPAdsW
!  nitrification funtions variables
   real(rk)     :: wNNitrW
!  denitrification functions variables
   real(rk)     :: wDDenitW,wNDenitW
!  total fluxes
   real(rk)     :: wDAbioIMW,wDAbioDetW,wNAbioDetW,wPAbioDetW,wSiAbioDetW
   real(rk)     :: wNAbioNH4W,wNAbioNO3W,wPAbioPO4W,wPAbioAIMW
   real(rk)     :: wSiAbioSiO2W,wO2AbioW

!
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_(self%id_sDIMW,sDIMW)
   _GET_(self%id_sDDetW,sDDetW)
   _GET_(self%id_sPDetW,sPDetW)
   _GET_(self%id_sNDetW,sNDetW)
   _GET_(self%id_sSiDetW,sSiDetW)
   _GET_(self%id_sPO4W,sPO4W)
   _GET_(self%id_sPAIMW,sPAIMW)
   _GET_(self%id_sNH4W,sNH4W)
   _GET_(self%id_sNO3W,sNO3W)
   _GET_(self%id_sO2W,sO2W)
   _GET_(self%id_sSiO2W,sSiO2W)

!  retrieve current environmental dependencies
   _GET_(self%id_uTm,uTm)

!  Nutrients ratio of detritus
  rPDDetW=sPDetW/(sDDetW+NearZero)
  rNDDetW= sNDetW/(sDDetW+NearZero)
!-----------------------------------------------------------------------
!  Temperature functions
!-----------------------------------------------------------------------
!  temp._function_of_mineralization_in_water
   uFunTmMinW = uFunTmAbio(uTm,self%cThetaMinW)
!  Temperature_dependence_for_nitrification
   uFunTmNitr = uFunTmAbio(uTm,self%cThetaNitr)
!-----------------------------------------------------------------------
!  mineralization functions
!-----------------------------------------------------------------------
!  P_mineralisation_constant_in_water
   kPMinDetW = self%kDMinDetW
!  N_mineralisation_constant_in_water
   kNMinDetW = self%kDMinDetW
!  Si_mineralisation_constant_in_water
   kSiMinDetW = self%kDMinDetW
!  decomposition
   wDMinDetW = self%kDMinDetW * uFunTmMinW * sDDetW
!  detrital_P mineralization
   wPMinDetW = kPMinDetW * uFunTmMinW * sPDetW
!  detrital_N mineralization
   wNMinDetW = kNMinDetW * uFunTmMinW * sNDetW
!  detrital_Si mineralization
   wSiMinDetW = kSiMinDetW * uFunTmMinW * sSiDetW
!-----------------------------------------------------------------------
!   Phosphrus sorption
!-----------------------------------------------------------------------
!  correction_of_O2_demand_in_water_at_low_oxygen_conc.
   aCorO2BOD = sO2W / (self%hO2BOD + sO2W)
!  max._P_adsorption_per_g_inorg._matter_in_water
   aPAdsMaxW = self%cRelPAdsD+aCorO2BOD*self%cRelPAdsFe*self%fFeDIM+ &
   & self%cRelPAdsAl*self%fAlDIM
!  P_adsorption_affinity,_corrected_for_redox_conditions
   aKPAdsW = (1.0_rk - self%fRedMax * (1.0_rk-aCorO2BOD)) * self%cKPAdsOx
!  P_adsorption_isotherm_onto_inorg._matter_in_sediment
   aPIsoAdsW = aPAdsMaxW * aKPAdsW * sPO4W / (1.0_rk + aKPAdsW * sPO4W)
!  equilibrium_conc.
   aPEqIMW = aPIsoAdsW * sDIMW
!  sorption_flux_in_water
   wPSorpIMW = self%kPSorp * (aPEqIMW - sPAIMW)
!-----------------------------------------------------------------------
!  nitrification functions
!-----------------------------------------------------------------------
   aCorO2NitrW = sO2W*sO2W / (self%hO2Nitr*self%hO2Nitr + sO2W*sO2W)
!  nitrification_flux
   wNNitrW = self%kNitrW * uFunTmNitr * aCorO2NitrW * sNH4W
!-----------------------------------------------------------------------
!  Denitrification functions
!-----------------------------------------------------------------------
!  mineralisation_flux_by_denitrification
   wDDenitW = sNO3W*sNO3W/(self%hNO3Denit*self%hNO3Denit+sNO3W*sNO3W)* &
   & (1.0_rk-aCorO2BOD)*wDMinDetW
!  Denitrification_flux
   wNDenitW=self%NO3PerC*molNmolC*self%cCPerDW*wDDenitW
!-----------------------------------------------------------------------
!  O2 dynamics
!-----------------------------------------------------------------------
!  O2_flux_due_to_mineralization_of_detritus
   wO2MinDetW = molO2molC * self%cCPerDW * aCorO2BOD * wDMinDetW
!  O2_flux_due_to_nitrification
   wO2NitrW = self%O2PerNH4 * molO2molN * wNNitrW
!-----------------------------------------------------------------------
!  total abiotic flux for each state variable in water
!   the abiot auxilaries have different meanings of the original one
!   just changed to its right hand side part, for benchmark test
!   but still keep the formula,fen, Sep17
!-----------------------------------------------------------------------
!  total_abiotic/microbial_DW_inorganic_matter_flux_in_water
   wDAbioIMW=0.0_rk
!  total_abiotic/microbial_DW_detritus_flux_in_water
   wDAbioDetW=-wDMinDetW
!  total_abiotic/microbial_N_detritus_flux_in_water
   wNAbioDetW =-wNMinDetW
!  total_abiotic/microbial_P_detritus_flux_in_water
   wPAbioDetW=-wPMinDetW
!  total_abiotic/microbial_Si_detritus_flux_in_water
   wSiAbioDetW=-wSiMinDetW
!  total_abiotic/microbial_dissolved_P_flux_in_water
   wPAbioPO4W=wPMinDetW-wPSorpIMW
!  total_abiotic/microbial_P_absorbed_onto_inorganic_matter_flux_in_water
   wPAbioAIMW=wPSorpIMW
!  total_abiotic/microbial_N_NH4_flux_in_water
   wNAbioNH4W = wNMinDetW-wNNitrW
!  total_abiotic/microbial_N_NO3_flux_in_water
   wNAbioNO3W = wNNitrW - wNDenitW
!  total_abiotic/microbial_Si_SiO2_flux_in_water
   wSiAbioSiO2W=wSiMinDetW
!  total_abiotic/microbial_O2_flux_in_water
   wO2AbioW = - wO2MinDetW - wO2NitrW
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
!  Set temporal derivatives
   _SET_ODE_(self%id_sDIMW,wDAbioIMW)
   _SET_ODE_(self%id_sDDetW,wDAbioDetW)
   _SET_ODE_(self%id_sPDetW,wPAbioDetW)
   _SET_ODE_(self%id_sNDetW,wNAbioDetW)
   _SET_ODE_(self%id_sSiDetW,wSiAbioDetW)
   _SET_ODE_(self%id_sPO4W,wPAbioPO4W)
   _SET_ODE_(self%id_sPAIMW,wPAbioAIMW)
   _SET_ODE_(self%id_sNH4W,wNAbioNH4W)
   _SET_ODE_(self%id_sNO3W,wNAbioNO3W)
   _SET_ODE_(self%id_sSiO2W,wSiAbioSiO2W)
   _SET_ODE_(self%id_sO2W,wO2AbioW)
!-----------------------------------------------------------------------
!  Output local diagnostic variables
!-----------------------------------------------------------------------
   _SET_DIAGNOSTIC_(self%id_wO2AbioW,wO2AbioW)
   _SET_DIAGNOSTIC_(self%id_rPDDetW,rPDDetW)
   _SET_DIAGNOSTIC_(self%id_rNDDetW,rNDDetW)


   _LOOP_END_
!-----------------------------------------------------------------------
!  Spatial loop end
!-----------------------------------------------------------------------
   end subroutine do

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
   class (type_pclake_abiotic_water), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s):
! !LOCAL VARIABLES:
   real(rk) :: sDIMW,sDDetW
   real(rk) :: extIM,extDet
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_sDIMW,sDIMW)
   _GET_(self%id_sDDetW,sDDetW)

   extIM=self%cExtSpIM*sDIMW
   extDet=self%cExtSpDet*sDDetW

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(extIM+extDet)

   _SET_DIAGNOSTIC_(self%id_extIM,extIM)
   _SET_DIAGNOSTIC_(self%id_extDet,extDet)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for Oxygen reaeration from the air
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_pclake_abiotic_water),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!  local variables
   real(rk)                   :: uVWind,sO2W,uTm
   real(rk)                   :: tO2Aer,uO2Sat,aFunLemnAer,kAer,uFunTmAer
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_
!  Retrive environmental dependencies
   _GET_HORIZONTAL_(self%id_uVWind,uVWind)
   _GET_(self%id_uTm,uTm)
!  Retrive state variable values
   _GET_(self%id_sO2W,sO2W)
!-----------------------------------------------------------------------
!     oxygen_reaeration functions
!-----------------------------------------------------------------------
!  temperature_function_of_reaeration
   uFunTmAer =  uFunTmAbio(uTm,self%cThetaAer)
!  oxygen_saturation_concentration
   uO2Sat = 14.652_rk - 0.41022 * uTm + 0.007991_rk * uTm*uTm - 0.000077774_rk * uTm*uTm*uTm
!  reaeration_coefficient

   kAer=(0.727_rk*((uVWind)**(0.5_rk))-0.371_rk*uVWind+0.0376_rk*uVWind*uVWind)
!  duckweed_function_of_reaeration
   aFunLemnAer = 1.0_rk
!  reaeration_flux_of_O2_into_the_water
   tO2Aer = kAer * uFunTmAer * (uO2Sat - sO2W) * aFunLemnAer
!-----------------------------------------------------------------------
!  update oxygen in rearation process
!-----------------------------------------------------------------------
!  for PCLake Benchmark
!  convert daily rates to sedonds


!  keep this, t is in day
   _SET_SURFACE_EXCHANGE_(self%id_sO2W,tO2Aer/secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tO2Aer,tO2Aer/secs_pr_day)
!   gotm output for 0d input
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wind,uVWind)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------
   end module pclake_abiotic_water
!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

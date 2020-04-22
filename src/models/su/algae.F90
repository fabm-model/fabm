#include "fabm_driver.h"

module su_algae

  use fabm_types

  implicit none

  private

  type (type_bulk_standard_variable), parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)

  type, extends(type_base_model), public :: type_su_algae
    ! Variables identifiers
    type(type_state_variable_id)      :: id_aC, id_aN, id_aP, id_aChl
    type(type_state_variable_id)      :: id_NH4, id_NO3, id_DIP, id_DIC, id_sDOMC, id_sDOMN, id_sDOMP, id_DOC 
    type(type_diagnostic_variable_id) :: id_aNC, id_aPC, id_aChlC, id_aPqmM
    type(type_dependency_id)          :: id_PAR, id_TChl

    ! Parameters
    real(rk) :: aM
    real(rk) :: LD
    real(rk) :: mix_depth
    real(rk) :: water_atten
    real(rk) :: ChlAtt
    real(rk) :: aalpha
    real(rk) :: aChlm
    real(rk) :: aBR
    real(rk) :: aNpref
    real(rk) :: aApref
    real(rk) :: aNKt
    real(rk) :: aAKt
    real(rk) :: aNCm
    real(rk) :: aNCmin
    real(rk) :: aNCabs
    real(rk) :: aPKu
    real(rk) :: aPCmax
    real(rk) :: aPCmin
    real(rk) :: aPCabs
    real(rk) :: aKQN
    real(rk) :: aKQP
    real(rk) :: aQh
    real(rk) :: aKxi
    real(rk) :: abeta
    real(rk) :: redco_1
    real(rk) :: AAsyn_1
    real(rk) :: aUmax
    real(rk) :: aDOCpc
    real(rk) :: aMetMult
    real(rk) :: aPpref
  contains
    ! Add model procedures here
    procedure :: initialize
    procedure :: do
  end type type_su_algae

contains

  subroutine initialize(self,configunit)                   
    class(type_su_algae), intent(inout), target :: self
    integer,              intent(in)            :: configunit

    ! Register variables in FABM
    call self%register_state_variable(self%id_aC,      'aC',    'ugC L-1',   'core algae C-biomass')
    call self%register_state_variable(self%id_aN,      'aN',    'ugN L-1',   'core algae N-biomass')
    call self%register_state_variable(self%id_aP,      'aP',    'ugP L-1',   'core algae P-biomass')
    call self%register_state_variable(self%id_aChl,    'aChl',  'ugChl L-1', 'core algae Chl-biomass')

    ! Register links to external state variables
    call self%register_dependency(self%id_TChl, total_chlorophyll)

    ! Register links to external nutrient pools
    call self%register_state_dependency(self%id_NH4,   'NH4',   'ugN L-1',   'ammonium')
    call self%register_state_dependency(self%id_NO3,   'NO3',   'ugN L-1',   'nitrate')
    call self%register_state_dependency(self%id_DIP,   'DIP',   'ugP L-1',   'dissolved inorganic phosphorus')
    call self%register_state_dependency(self%id_DIC,   'DIC',   'ugC L-1',   'dissolved inorganic carbon')
    call self%register_state_dependency(self%id_sDOMC, 'sDOMC', 'ugC L-1',   'semi-labile dissolved organic material-carbon')
    call self%register_state_dependency(self%id_sDOMN, 'sDOMN', 'ugN L-1',   'semi-labile dissolved organic material-nitrogen')
    call self%register_state_dependency(self%id_sDOMP, 'sDOMP', 'ugP L-1',   'semi-labile dissolved organic material-phosphorus')
    call self%register_state_dependency(self%id_DOC,   'DOC',   'ugC L-1',   'dissolved organic carbon')

    ! Register diagnostic variables
    call self%register_diagnostic_variable(self%id_aNC,    'aNC',   'gN (gC)-1',         'core algae N:C')
    call self%register_diagnostic_variable(self%id_aPC,    'aPC',   'gP (gC)-1',         'core algae P:C')
    call self%register_diagnostic_variable(self%id_aChlC,  'aChlC', 'gChl (gC)-1',       'core algae Chl:C')
    call self%register_diagnostic_variable(self%id_aPqmM,  'aPqmM', 'gC (gC)-1 (day)-1', 'mixotroph kleptochloroplastic maximum possible Pmax for prey')

    ! Register links to external dependencies 
    call self%register_dependency(self%id_PAR,standard_variables%downwelling_photosynthetic_radiative_flux)

    ! Register the contribution of variables to total N,P, and C
    call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_aN, scale_factor=1._rk/14._rk)
    call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_aP, scale_factor=1._rk/31._rk)
    call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_aC, scale_factor=1._rk/12._rk)
    call self%add_to_aggregate_variable(total_chlorophyll,                   self%id_aChl) ! no need for scale_factor (ug/L = mg/m3)

    ! Register parameters in FABM
    call self%get_parameter(self%aM,          'aM',          'dl',                'algal scalar for controlling photoacclimation',                                  default=3._rk)
    call self%get_parameter(self%LD,          'LD',          'dl',                         'fraction of day illuminated',                                       default=1._rk)
    call self%get_parameter(self%mix_depth,   'mix_depth',   'm',                          'mixed layer depth',                                             default=10._rk)
    call self%get_parameter(self%water_atten, 'water_atten', 'm-1',                        'attenuation of light by water',                                             default=0.032323_rk)
    call self%get_parameter(self%ChlAtt,      'ChlAtt',      'm2 (mgChl)-1',               'attenuation of light by Chl',                                               default=0.02_rk)
    call self%get_parameter(self%aalpha,      'aalpha',      '(m2/gChla)*(gC/umolphoton)', 'algae Chl-specific initial slope to PI curve',                                             default=0.000005_rk)
    call self%get_parameter(self%aChlm,       'aChlm',       'gChla (gC)-1',                         'algal maximum possible ChlC',                                              default=0.06_rk)
    call self%get_parameter(self%aBR,         'aBR',         'gC (gC)-1 (day)-1',          'algae basal respiration rate',                                              default=0.05_rk)
    call self%get_parameter(self%aNpref,      'aNpref',      'dl',                         'algae relative preference for NO3 usage',                                             default=1._rk)
    call self%get_parameter(self%aApref,      'aApref',      'dl',                         'algae relative preference for NH4 usage',                                             default=2._rk)
    call self%get_parameter(self%aNKt,        'aNKt',        'ugN L-1',                    'algae half saturation constant for NO3 transport',                                         default=14._rk)
    call self%get_parameter(self%aAKt,        'aAKt',        'ugN L-1',                    'algae half saturation constant for NH4 transport',                                         default=14._rk)
    call self%get_parameter(self%aNCm,        'aNCm',        'gN (gC)-1',                  'algae maximum N:C quota',                                             default=0.25_rk)
    call self%get_parameter(self%aNCmin,      'aNCmin',      'gN (gC)-1',                  'algae minimum N:C quota',                                             default=0.06_rk)
    call self%get_parameter(self%aNCabs,      'aNCabs',      'gN (gC)-1',                  'algae absolute maximum N:C quota',                                             default=0.25_rk)
    call self%get_parameter(self%aPKu,        'aPKu',        'ugP L-1',                    'algae half-saturation constant for DIP uptake',                                            default=31._rk)
    call self%get_parameter(self%aPCmax,      'aPCmax',      'gP (gC)-1',                  'algae maximum P:C quota that affects growth',                                            default=0.02_rk)
    call self%get_parameter(self%aPCmin,      'aPCmin',      'gP (gC)-1',                  'algae minimum P:C quota',                                             default=0.005247_rk)
    call self%get_parameter(self%aPCabs,      'aPCabs',      'gP (gC)-1',                  'algae absolute maximum P:C quota',                                             default=0.04_rk)
    call self%get_parameter(self%aKQN,        'aKQN',        'dl',                         'algae half-saturation for N quota curve',                                             default=10._rk)
    call self%get_parameter(self%aKQP,        'aKQP',        'dl',                         'algae half-saturation for P quota curve',                                             default=0.1_rk)
    call self%get_parameter(self%aQh,         'aQh',         'dl',                         'algae Hill number for control of nutrient uptake',                                            default=4._rk)
    call self%get_parameter(self%aKxi,        'aKxi',        'dl',                         'algae half-saturation constant for control of nutrient uptake',                                            default=0.1_rk)
    call self%get_parameter(self%abeta,       'abeta',       'dl',                         'algae power for controlling uptake of non-limiting nutrient',                                          default=0.05_rk)
    call self%get_parameter(self%redco_1,     'redco_1',     'gC (gN)-1',                  'algae C respired to support nitrate reduction to intracellular NH4',                                 default=1.71_rk)
    call self%get_parameter(self%AAsyn_1,     'AAsyn_1',     'gC (gN)-1',                  'algae amino acid synthesis cost',                                              default=1.5_rk)
    call self%get_parameter(self%aUmax,       'aUmax',       'gC (gC)-1 (day)-1',          'algae maximum rate of photosynthesis-driven growth',                                            default=1._rk)
    call self%get_parameter(self%aDOCpc,      'aDOCpc',      '%',                          'algae % of Cfix released as DOC',                                               default=0.1_rk)
    call self%get_parameter(self%aMetMult,    'aMetMult',    'dl',                         'algae metabolic cost multiplier to achieve required GGE for dinos',                                             default=3._rk)
    call self%get_parameter(self%aPpref,      'aPpref',      'dl',                         'algae preference for DIP',                                               default=5._rk)
    self%dt = 86400._rk          ! Change time unit of all parameters to meet FABM default (seconds)
  end subroutine initialize

  subroutine do(self, _ARGUMENTS_DO_)
      class (type_su_algae), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: aC, aN, aP, aChl                             ! internal state variables
      real(rk) :: NH4, NO3, DIP, DIC, sDOMC, sDOMN, sDOMP, DOC ! external state variables
      real(rk) :: aUmpN, aNVP, aAVP, aAV, aNV, afrat           ! Nitrogen transport and uptake
      real(rk) :: aNC, aNCu, aPC, aPCu, aNPCu                  ! Internal N and P quotas
      real(rk) :: frNC, aupN, aupN_aC, aNH4up, aNO3up, frPC, aupP, aupP_aC ! Actual uptake of N and P
      real(rk) :: aBRop, aPqm, aRphot, PFD, TChl               ! Respiration and photosynthesis
      real(rk) :: attenuation, exat, aChlC, aPyt, aPS          ! Respiration and photosynthesis
      real(rk) :: aCu, daChlC, daChlC_aC, agro                 ! Algae growth and chlorophyll synthesis 
      real(rk) :: aICout, aICup, aDOCout                       ! Algae contributions to DIC/DOC pool       
      real(rk) :: aPqmM                                        ! To use in the mixotroph model

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

      ! Retrieve internal state variables values
      _GET_(self%id_aC,   aC)
      _GET_(self%id_aN,   aN)
      _GET_(self%id_aP,   aP)
      _GET_(self%id_aChl, aChl)

      ! Retrieve external state variables values
      _GET_(self%id_TChl, TChl)
      _GET_(self%id_NH4,   NH4)
      _GET_(self%id_NO3,   NO3)
      _GET_(self%id_DIP,   DIP)
      _GET_(self%id_DIC,   DIC)
      _GET_(self%id_sDOMC, sDOMC)
      _GET_(self%id_sDOMN, sDOMN)
      _GET_(self%id_sDOMP, sDOMP)
      _GET_(self%id_DOC,   DOC)

      ! Retrieve external dependencies
      _GET_(self%id_PAR,PFD) 
      
      ! Calculations:

!-------------------------------------
!-----------Nutrient uptake-----------
!-------------------------------------

      ! Nitrogen (NO3 and NH4) transport and uptake:
      ! Total algal N uptake (gN/gC/day):
      aUmpN = self%aUmax * self%aNCm

      ! Algal potential relative NO3 trasnport (gN/gC/day):
      if (NO3 > 0._rk) then
      aNVP = aUmpN * self%aNpref * NO3/(NO3 + self%aNKt)
      else
      aNVP = 0._rk
      end if

      ! Algal potential relative NH4 trasnport (gN/gC/day):
      if (NH4 > 0._rk) then
      aAVP = aUmpN * self%aApref * NH4/(NH4 + self%aAKt)
      else
      aAVP = 0._rk
      end if

      ! Algal NH4 uptake (gN/gC/day):
      if (aUmpN < aAVP) then
      aAV = aUmpN
      else
      aAV = aAVP
      end if 

      ! Algal NO3 uptake (gN/gC/day):
      if (aAVP < aUmpN) then
      	 if (aAVP + aNVP < aUmpN) then
         aNV = aNVP
         else
         aNV = aUmpN - aAVP
         end if
      else
      aNV = 0._rk
      end if

      ! Algal f-ratio (dl):
      afrat = aNV/(aNV + aAV + 1e-12_rk)

      ! NC internal quota (gN/gC):
      aNC = aN/aC

      ! PC internal quota (gP/gC):
      aPC = aP/aC

      ! Algae normalised NC quota description (dl): 
      if (aNC <= self%aNCm) then
         if (aNC >= self%aNCmin) then
         aNCu = (1._rk + self%aKQN)*(aNC - self%aNCmin)/((aNC - self%aNCmin) + self%aKQN*(self%aNCm - self%aNCmin))
         else
         aNCu = 0._rk
         end if
      else
      aNCu = 1._rk
      end if

      ! Algae normalised PC quota description (dl): 
      if (aPC < self%aPCmax) then
         if (aPC > self%aPCmin) then
         aPCu = (1._rk + self%aKQP)*(aPC - self%aPCmin)/((aPC - self%aPCmin) + self%aKQP*(self%aPCmax - self%aPCmin))
         else
         aPCu = 0._rk
         end if
      else
      aPCu = 1._rk
      end if

      ! Algal threshold selection of phototrophic growth control (by N or P status) (dl):
      aNPCu = MIN(1._rk,MIN(aNCu, aPCu))

      ! Normalised feedback from NC quota (dl)
      frNC = (1._rk + self%aKxi**self%aQh)*(1._rk - aNC/self%aNCabs)**self%aQh/((1._rk - aNC/self%aNCabs)**self%aQh + self%aKxi**self%aQh)

      ! Algal total DIN uptake (gN/gC/day):
      if (aNC < self%aNCabs - 0.002_rk) then 
         if (aNCu > aNPCu) then
         aupN = (aAV + aNV)*aNPCu**self%abeta*frNC
         else 
         aupN = (aAV + aNV)*frNC
         end if
      else
      aupN = 0._rk
      end if

      ! Algae population total DIN uptake (ugN/L/day):
      aupN_aC = aupN * aC 

      ! Algae population total NH4 uptake (ugN/L/day):
      aNH4up = aupN_aC*(1._rk - afrat) 

      ! Algae population total NO3 uptake (ugN/L/day):
      aNO3up = aupN_aC*afrat

      ! Normalised feedback from PC quota (dl)
      frPC = DIP/(DIP + self%aPKu)*(1._rk + self%aKxi**self%aQh)*(1._rk - aPC/self%aPCabs)**self%aQh/((1._rk - aPC/self%aPCabs)**self%aQh + self%aKxi**self%aQh)

      ! Algal total DIP uptake (gP/gC/day):
      if (DIP > 0._rk .and. aPC < self%aPCabs) then 
         if (aPCu > aNPCu) then
         aupP = self%aUmax*self%aPpref*self%aPCmax*aNPCu**self%abeta*frPC
         else 
         aupP = self%aUmax*self%aPpref*self%aPCmax*frPC
         end if
      else
      aupP = 0._rk
      end if

      ! Algae population total DIP uptake (ugP/L/day):
      aupP_aC = aupP * aC 

!--------------------------------------
!----Photosynthesis and respiration----
!--------------------------------------

      ! Algal basal respiration rate (halting respiration at high NC) (gC/gC/day):
      if (aNC < self%aNCabs) then
      aBRop = self%aUmax*self%aBR*1.01_rk*((self%aNCabs - aNC)/(self%aNCabs - self%aNCmin))/((self%aNCabs - aNC)/(self%aNCabs - self%aNCmin)+0.01_rk)
      else
      aBRop = 0._rk
      end if

      ! Algae potential maximum photosynthesis rate (gC/gC/day):
      aPqm = (self%aUmax+aBRop+self%aNCm*self%aUmax*(self%redco_1+self%AAsyn_1*self%aMetMult))*aNPCu+1e-6_rk

      ! Algae phototrophic driven respiration (gC/gC/day):
      aRphot = self%redco_1*aupN*afrat + aupN*self%AAsyn_1*self%aMetMult

      ! Attenuation of light by water and chlorophyll (ugChl/L):
      attenuation = self%mix_depth * (self%water_atten + TChl*self%ChlAtt)
      exat = exp(-attenuation)

      ! Algae Chl:C quota (gChl/gC):
      aChlC = aChl/aC

      ! Algal photosynthesis rate according to the Smith equation (gC/gC/day):
      aPyt = (self%aalpha*aChlC*PFD*24._rk*60._rk*60._rk)/aPqm

      ! Algal depth integrated photosynthesis rate (gC/gC/day):
      aPS = aPqm*(log(aPyt + SQRT(1._rk + aPyt**2._rk)) - log(aPyt*exat + SQRT(1._rk + (aPyt*exat)**2._rk)))/attenuation

      ! Net growth rate (gC/gC/day):
      aCu = aPS - aBRop - aRphot

!---------------------------------------
!---------Chlorophyll synthesis---------
!---------------------------------------

      ! Algal change in ChlC with synthesis and degradation with C growth (gChl/gC/day):
      if (aChlC < self%aChlm) then
      daChlC =  self%aChlm*aNPCu*self%aM*self%aUmax*(1._rk - MIN(1._rk,aPS/aPqm))**0.5_rk*(1._rk + 0.05_rk)*(1._rk - aChlC/self%aChlm)/(1._rk - aChlC/self%aChlm + 0.05_rk) - aChlC*((1._rk - aNCu)*aCu)
      else
      daChlC = 0._rk
      end if

      ! daChlC in units of ugChl/L/day:
      daChlC_aC = daChlC*aC

!----------------------------------
!---------Population rates---------
!----------------------------------

      ! Algae population growth (ugC/L/day):
      agro = aCu*aC

      ! Algae contribution to DIC pool (ugC/L/day):
      aICout = (aBRop + aRphot)*aC

      ! Total algae DIC usage (ugC/L/day):
      aICup = aPS*aC
      aDOCout = aICup*self%aDOCpc

      !Mixotroph kleptochloroplastic maximum possible Pmax for prey; assumes that on ingestion the nutrient status within the kleptochloroplasts is enhanced to maximum (gC/gC/d):
      aPqmM = (self%aUmax + self%aBR + self%aNCm * self%aUmax * (self%redco_1 + self%AAsyn_1 * self%aMetMult)) + 1e-6_rk

      ! Set ODEs:
      _ADD_SOURCE_(self%id_aC,   + agro)
      _ADD_SOURCE_(self%id_aN,   + aupN_aC)
      _ADD_SOURCE_(self%id_aP,   + aupP_aC)
      _ADD_SOURCE_(self%id_aChl, + daChlC_aC)

      _ADD_SOURCE_(self%id_NH4,  - aNH4up)
      _ADD_SOURCE_(self%id_NO3,  - aNO3up)
      _ADD_SOURCE_(self%id_DIP,  - aupP_aC)
      _ADD_SOURCE_(self%id_DIC,  + aICout - aICup - aDOCout)
      _ADD_SOURCE_(self%id_DOC,  + aDOCout)

      ! Set diagnostic variables:
      _SET_DIAGNOSTIC_(self%id_aPqmM, aPqmM)
      _SET_DIAGNOSTIC_(self%id_aNC, aNC)
      _SET_DIAGNOSTIC_(self%id_aPC, aPC)
      _SET_DIAGNOSTIC_(self%id_aChlC, aChlC)

      ! Leave spatial loops (if any)
      _LOOP_END_

  end subroutine do

end module su_algae

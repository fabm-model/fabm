#include "fabm_driver.h"

module su_bacteria

  use fabm_types

  implicit none

  private

  type, extends(type_base_model), public :: type_su_bacteria
    ! Variables identifiers
    type(type_state_variable_id)      :: id_bC
    type(type_state_variable_id)      :: id_NH4, id_NO3, id_DIP, id_DIC, id_DOC, id_sDOMC, id_sDOMN, id_sDOMP, id_cDOMC, id_cDOMN, id_cDOMP
    type(type_diagnostic_variable_id) :: id_bOCV, id_bsOMCV, id_bbrelN, id_bbrelP, id_bINV, id_bsOMNV, id_bVm, id_btotCV, id_btotNV, id_btotPV, id_bgroCu, id_bIPV, id_bsOMPV, id_bCup, id_bICout !,id_bN, id_bP, 

    ! Parameters
    real(rk) :: brbasN
    real(rk) :: brbasP
    real(rk) :: brbasC
    real(rk) :: brelpref
    real(rk) :: brelcXpref
    real(rk) :: brelsCpref
    real(rk) :: relcCpref
    real(rk) :: bKN
    real(rk) :: bKP
    real(rk) :: bKcN
    real(rk) :: bKcP
    real(rk) :: bKOC
    real(rk) :: bKsOMC
    real(rk) :: bKcOMC
    real(rk) :: bNC
    real(rk) :: bmNC
    real(rk) :: bPC
    real(rk) :: bmPC
    real(rk) :: bUm
    real(rk) :: bINrs
    real(rk) :: bONrs
    real(rk) :: bcM
  contains
    ! Add model procedures here
    procedure :: initialize
    procedure :: do
  end type type_su_bacteria

contains

  subroutine initialize(self,configunit)                   
    class(type_su_bacteria), intent(inout), target :: self
    integer,                 intent(in)            :: configunit

    ! Register variables in FABM
    call self%register_state_variable(self%id_bC,      'bC',    'ugC L-1', 'core bacterial C-biomass')

    ! Register links to external state variables (mixotroph and nutrient pools)  
    call self%register_state_dependency(self%id_NH4,   'NH4',   'ugN L-1', 'ammonium')
    call self%register_state_dependency(self%id_NO3,   'NO3',   'ugN L-1', 'nitrate')
    call self%register_state_dependency(self%id_DIP,   'DIP',   'ugP L-1', 'dissolved inorganic phosphorus')
    call self%register_state_dependency(self%id_DIC,   'DIC',   'ugC L-1', 'dissolved inorganic carbon')
    call self%register_state_dependency(self%id_DOC,   'DOC',   'ugC L-1', 'dissolved organic carbon')
    call self%register_state_dependency(self%id_sDOMC, 'sDOMC', 'ugC L-1', 'semi-labile dissolved organic material-carbon')
    call self%register_state_dependency(self%id_sDOMN, 'sDOMN', 'ugN L-1', 'semi-labile dissolved organic material-nitrogen')
    call self%register_state_dependency(self%id_sDOMP, 'sDOMP', 'ugP L-1', 'semi-labile dissolved organic material-phosphorus')
    call self%register_state_dependency(self%id_cDOMC, 'cDOMC', 'ugC L-1', 'combined dissolved organic material-carbon')
    call self%register_state_dependency(self%id_cDOMN, 'cDOMN', 'ugN L-1', 'combined dissolved organic material-nitrogen')
    call self%register_state_dependency(self%id_cDOMP, 'cDOMP', 'ugP L-1', 'combined dissolved organic material-phosphorus')

    ! Register diagnostic variables
    !call self%register_diagnostic_variable(self%id_bN,     'bN',    'ugN L-1', 'bacterial biomass in terms of nitrogen per water volume')
    !call self%register_diagnostic_variable(self%id_bP,     'bP',    'ugP L-1', 'bacterial biomass in terms of phosphorus per water volume')
    call self%register_diagnostic_variable(self%id_bOCV,   'bOCV',    '', '')
    call self%register_diagnostic_variable(self%id_bsOMCV, 'bsOMCV',    '', '')
    call self%register_diagnostic_variable(self%id_bbrelN, 'bbrelN',    '', '')
    call self%register_diagnostic_variable(self%id_bbrelP, 'bbrelP',    '', '')
    call self%register_diagnostic_variable(self%id_bINV,   'bINV',    '', '')
    call self%register_diagnostic_variable(self%id_bsOMNV, 'bsOMNV',    '', '')
    call self%register_diagnostic_variable(self%id_bVm,    'bVm',    '', '')
    call self%register_diagnostic_variable(self%id_btotCV, 'btotCV',    '', '')
    call self%register_diagnostic_variable(self%id_btotNV, 'btotNV',    '', '')
    call self%register_diagnostic_variable(self%id_btotPV, 'btotPV',    '', '')
    call self%register_diagnostic_variable(self%id_bgroCu, 'bgroCu',    '', '')
    call self%register_diagnostic_variable(self%id_bIPV,   'bIPV',    '', '')
    call self%register_diagnostic_variable(self%id_bsOMPV, 'bsOMPV',    '', '')
    call self%register_diagnostic_variable(self%id_bCup,   'bCup',    '', '')
    call self%register_diagnostic_variable(self%id_bICout, 'bICout',    '', '')

    ! Register parameters in FABM
    call self%get_parameter(self%brbasN,     'brbasN',     'gC (gC)-1 (day)-1','bacterial basal respiration rate when N is limiting',   default=0.1_rk)
    call self%get_parameter(self%brbasP,     'brbasP',     'gC (gC)-1 (day)-1','bacterial basal respiration rate when P is limiting',   default=0.2_rk)
    call self%get_parameter(self%brbasC,     'brbasC',     'gC (gC)-1 (day)-1','bacterial basal respiration rate when just C is limiting',                                                                                                                              default=0.05_rk)
    call self%get_parameter(self%brelpref,   'brelpref',   'dl','bacterial relative preference of DOMX over DIX; must be >=1; if preference is for DIX then this index must be pointed to DIX',                                                                                                default=1._rk)
    call self%get_parameter(self%brelcXpref, 'brelcXpref', 'dl','bacterial relative C preference',                                      default=1._rk)
    call self%get_parameter(self%brelsCpref, 'brelsCpref', 'dl','bacterial relative preference of DOM-C over DOC; must be >=1; if preference is for DOC then this index must be pointed to DOC',                                                                                                default=1._rk)
    call self%get_parameter(self%relcCpref,  'relcCpref',  'dl',                'bacterial relative preference for combined C',         default=0.25_rk)
    call self%get_parameter(self%bKN,        'bKN',        'ugN L-1',           'bacterial half saturation constant for DIN uptake',    default=0.14_rk)
    call self%get_parameter(self%bKP,        'bKP',        'ugP L-1',           'bacterial half saturation constant for DIP uptake',    default=0.31_rk)
    call self%get_parameter(self%bKcN,       'bKcN',       'ugN L-1',           'bacterial half saturation constant for combined DOM-N',default=14._rk)
    call self%get_parameter(self%bKcP,       'bKcP',       'ugP L-1',           'bacterial half saturation constant for combined DOM-P',default=0.31_rk)
    call self%get_parameter(self%bKOC,       'bKOC',       'ugC L-1',           'bacterial half saturation constant for DOC-C',         default=12._rk)
    call self%get_parameter(self%bKsOMC,     'bKsOMC',     'ugC L-1',           'bacterial half saturation constant for semi-labile DOM-C',default=12._rk)
    call self%get_parameter(self%bKcOMC,     'bKcOMC',     'ugC L-1',           'bacterial half saturation constant for combined DOM-C',default=24._rk)
    call self%get_parameter(self%bNC,        'bNC',        'gN (gC)-1',         'bacterial N:C quota',                                  default=0.2_rk)
    call self%get_parameter(self%bmNC,       'bmNC',       'gN (gC)-1',         'bacterial maximum N:C for usage',                      default=0.5_rk)
    call self%get_parameter(self%bPC,        'bPC',        'gP (gC)-1',         'bacterial P:C quota',                                  default=0.02_rk)
    call self%get_parameter(self%bmPC,       'bmPC',       'gP (gC)-1',         'bacterial maximum P:C for usage',                      default=0.05_rk)
    call self%get_parameter(self%bUm,        'bUm',        'gC (gC)-1 (day)-1', 'bacterial maximum net growth rate',                    default=3._rk)
    call self%get_parameter(self%bINrs,      'bINrs',      'gC (gN)-1',         'bacterial respiratory cost for the assimilation of ammonium; gC respired per gDIN incorporated (E.coli data assuming N is as NH4 and C is sugars)',                                                     default=1.46_rk)
    call self%get_parameter(self%bONrs,      'bONrs',      'gC (gN)-1',         'bacterial respiratory cost for the assimilation of amino acids; gC respired per gDON assimilated, assumes DON is amino acids BUT that C is as sugars',                                                    default=0.857_rk)
    call self%get_parameter(self%bcM,        'bcM',        'dl',                'bacterial cost multiplier to alter the respiration rate of the cell',                                                                                                                                 default=2._rk)

    ! Register the contribution of variables to total N,P, and C (system balance)
    call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_bC, scale_factor=self%bNC/14._rk)
    call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_bC, scale_factor=self%bPC/31._rk)
    call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_bC, scale_factor=1._rk/12._rk)

    self%dt = 86400._rk      ! Change time unit of all parameters to meet FABM default (seconds)
  end subroutine initialize

  subroutine do(self,_ARGUMENTS_DO_)
      class (type_su_bacteria), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: bC                                                ! internal state variables
      real(rk) :: NH4, NO3, DIP, DIC, DOC                           ! external state variables
      real(rk) :: sDOMC, sDOMN, sDOMP, cDOMC, cDOMN, cDOMP          ! external state variables
      real(rk) :: bVm, bOCV, bsOMCV, bcOMCV, btotCV                 ! auxiliaries for C usage
      real(rk) :: bINV, bsOMNV, bcOMNV, btotNV, bbrelN              ! auxiliaries for N usage
      real(rk) :: bIPV, bsOMPV, bcOMPV, btotPV, bbrelP              ! auxiliaries for P usage
      real(rk) :: bgroCu, bOCuse, bsOMu, bcOMu                      ! rates of C-sources usage
      real(rk) :: bNrsyn, bgroNu, bINu, bsONu, bcONu                ! Rates of N-sources use
      real(rk) :: bgroPu, bIPu, bsOPu, bcOPu                        ! Rates of P-sources use
      real(rk) :: bDOCup, bcOMCup, bsOMCup, bNH4up, bsOMNup, bcOMNup, bDIPup, bsOMPup, bcOMPup ! Bact pop
      real(rk) :: brN, brP, brbas, brelU, bBR, bmCres, bCRes, bCu   ! Bacterial respiration and growth
      real(rk) :: bsONC, bcONC, bsOPC, bcOPC                        ! External NC and PC quotas
      real(rk) :: bConsON, bConcON, bConsOP, bConcOP, bregN, bregP  ! Control of N and P usage and reg
      real(rk) :: bN, bP, bCup, bICout                              ! N and P biomass; bC flows

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

      ! Retrieve internal state variables values
      _GET_(self%id_bC,    bC)

      ! Retrieve external state variables values
      _GET_(self%id_NH4,   NH4)
      _GET_(self%id_NO3,   NO3)
      _GET_(self%id_DIP,   DIP)
      _GET_(self%id_DIC,   DIC)
      _GET_(self%id_DOC,   DOC)
      _GET_(self%id_sDOMC, sDOMC)
      _GET_(self%id_sDOMN, sDOMN)
      _GET_(self%id_sDOMP, sDOMP)
      _GET_(self%id_cDOMC, cDOMC)
      _GET_(self%id_cDOMN, cDOMN)
      _GET_(self%id_cDOMP, cDOMP)

      ! Calculations:

      ! Gross growth rate (gC/gC/day):
      bVm = self%bUm + self%bcM * self%bUm * self%bONrs * self%bNC

      ! Quotient for uptake of DOC (dl):
      if (DOC > 0._rk) then
         bOCV = DOC/(DOC + self%bKOC)
      else
         bOCV = 0._rk
      end if

      ! Bacterial sDOM-C availability; quotient for uptake of C from sDOM, with preference scalar (dl):
      if (sDOMC> 0._rk) then
         bsOMCV = self%brelsCpref * sDOMC/(sDOMC + self%bKsOMC)
      else
         bsOMCV = 0._rk
      end if

      ! Bacterial cDOM-C availability; quotient for uptake of C from cDOM, with preference scalar (dl):
      if (cDOMC> 0._rk) then
         bcOMCV = self%relcCpref * cDOMC/(cDOMC + self%bKcOMC)
      else
         bcOMCV = 0._rk
      end if

      ! Actual rate of C-uptake (gC/gC/day):
      btotCV = MIN(bVm, bVm*(bOCV + bsOMCV + bcOMCV))

      ! Bacterial NH4 availability (dl)
      if (NH4> 0._rk) then
         bINV = NH4/(NH4 + self%bKN)
      else
         bINV = 0._rk
      end if
      
      ! Bacterial sDOM-N availability (dl):
      if (sDOMN> 0._rk) then
         bsOMNV = self%brelpref * sDOMN/(sDOMN + self%bKN)
      else
         bsOMNV = 0._rk
      end if

      ! Bacterial cDOM-N availability (dl):
      if (cDOMN> 0._rk) then
         bcOMNV = self%brelcXpref * cDOMN/(cDOMN + self%bKcN)
      else
         bcOMNV = 0._rk
      end if

      ! Actual rate of N-uptake (gN/gC/day):
      btotNV = MIN(bVm*self%bNC, bVm*self%bNC*(bINV + bsOMNV + bcOMNV))

      ! N:C stoichiometric balance (dl):
      if(btotCV == 0._rk) then
      bbrelN = 0._rk
      else
      bbrelN = (btotNV/btotCV)/self%bNC  ! if > 1, no N limitation; if < 1 N-limited.
      end if
      
      ! Bacterial DIP availability (dl)
      if (DIP> 0._rk) then
         bIPV = DIP/(DIP + self%bKP)
      else
         bIPV = 0._rk
      end if
      
      ! Bacterial sDOM-P availability (dl):
      if (sDOMP> 0._rk) then
         bsOMPV = self%brelpref * sDOMP/(sDOMP + self%bKP)
      else
         bsOMPV = 0._rk
      end if

      ! Bacterial cDOM-P availability (dl):
      if (cDOMP> 0._rk) then
         bcOMPV = self%brelcXpref * cDOMP/(cDOMP + self%bKcP)
      else
         bcOMPV = 0._rk
      end if

      ! Actual rate of P-uptake (gP/gC/day):
      btotPV = MIN(bVm*self%bPC, bVm*self%bPC*(bIPV + bsOMPV + bcOMPV))

      ! P:C stoichiometric balance (dl):
      if(btotCV == 0._rk) then
      bbrelP = 0._rk
      else
      bbrelP = (btotPV/btotCV)/self%bPC  ! if > 1, no P limitation; if < 1 P-limited.
      end if

      ! Total need for C by bacteria (gC/gC/day):
      bgroCu = MIN(bbrelN,bbrelP,1._rk)*btotCV

      ! Bacterial use of DOC (gC/gC/day):
      bOCuse = bgroCu*bOCV/(bOCV + bsOMCV + bcOMCV)

      ! Bacterial use of C from semi-labile DOM (gC/gC/day):
      bsOMu = bgroCu*bsOMCV/(bOCV + bsOMCV + bcOMCV)

      ! Bacterial use of C from combined DOM (gC/gC/day):
      bcOMu = bgroCu - bOCuse - bsOMu

      ! Bacterial cost of N-assimilation (gC/gN):
      bNrsyn = self%bcM*(self%bINrs*bINV + self%bONrs*(bsOMNV + bcOMNV))/(bINV + bsOMNV + bcOMNV + 1e-23_rk) 

      ! Total need for N by bacteria (gN/gC/day):
      bgroNu = bgroCu/(1._rk/self%bNC + bNrsyn)

      ! Bacterial DIN use (gN/gC/day):
      bINu = bgroNu*bINV/(bINV + bsOMNV + bcOMNV + 1e-23_rk)

      ! Bacterial sDOM-N use (gN/gC/day):
      bsONu = bgroNu*bsOMNV/(bINV + bsOMNV + bcOMNV + 1e-23_rk)

      ! Bacterial cDOM_N use (gN/gC/day):
      bcONu = bgroNu - bINu - bsONu

      ! Total need for P by bacteria (gP/gC/day):
      bgroPu = bgroNu*self%bPC/self%bNC

      ! Bacterial DIP use (gP/gC/day):
      bIPu = bgroPu*bIPV/(bIPV + bsOMPV + bcOMPV + 1e-23_rk)

      ! Bacterial sDOM-P use (gP/gC/day):
      bsOPu = bgroPu*bsOMPV/(bIPV + bsOMPV + bcOMPV + 1e-23_rk)

      ! Bacterial cDOM-P use (gP/gC/day):
      bcOPu = bgroPu - bIPu - bsOPu

      ! Rates of uptake of C,N, and P sources by bacterial population (ugX/L/day)
      bDOCup = bOCuse*bC

      bcOMCup = bcOMu*bC

      bsOMCup = bsOMu*bC

      if (NH4 > 0._rk) then
      bNH4up = bINu*bC
      else
      bNH4up = 0._rk
      end if

      bsOMNup = bsONu*bC

      bcOMNup = bcONu*bC

      if (DIP > 0._rk) then
      bDIPup = bIPu*bC
      else
      bDIPup = 0._rk
      end if

      bsOMPup = bsOPu*bC

      bcOMPup = bcOPu*bC

      ! Basal respiration rate control; contains a switch to set basal respiration depending on whether N or P is limiting (gC/gC/day):
      if(bbrelN > 1._rk) then
      brN = 0._rk
      else
      brN = (1._rk-bbrelN)*(self%brbasN - self%brbasC)
      end if

      if(bbrelP > 1._rk) then
      brP = 0._rk
      else
      brP = (1._rk - bbrelP)*(self%brbasP - self%brbasC)
      end if

      brbas = self%brbasC + MAX(brN, brP)

      ! Bacterial relative growth rate (dl):
      brelU = bgroCu/bVm

      ! Bacterial basal respiration (gC/gC/day):
      bBR = brbas*self%bUm*(1._rk - brelU)

      ! Bacterial metabolic respiration (gC/gC/day):
      bmCres = ((bsONu*self%bONrs+bINu*self%bINrs)+bcONu*self%bONrs)*self%bcM

      ! Respiration of C from bacteria associated with basal and metabolic respiration (gC/gC/day):
      bCRes = bmCres + bBR

      ! Bacterial net growth rate (gC/gC/day):
      bCu = bOCuse + bsOMu + bcOMu - bCRes

      ! External N:C (gN/gC) and P:C (gP/gC) quotas:
      if (sDOMN > 0._rk) then
      bsONC = sDOMN/sDOMC
      else
      bsONC = 0._rk
      end if

      if (cDOMN > 0._rk) then
      bcONC = cDOMN/cDOMC
      else
      bcONC = 0._rk
      end if

      if (sDOMP > 0._rk) then
      bsOPC = sDOMP/sDOMC
      else
      bsOPC = 0._rk
      end if

      if (cDOMP > 0._rk) then
      bcOPC = cDOMP/cDOMC
      else
      bcOPC = 0._rk
      end if

      ! Control of N and P usage based on internal and external XC quotas (dl):
      if (bsONC > self%bmNC .and. sDOMN > 0._rk) then
      bConsON = bsONC/self%bmNC*2._rk
      else
      bConsON = 0._rk
      end if

      if (bcONC > self%bmNC .and. cDOMN > 0._rk) then
      bConcON = bcONC/self%bmNC*2._rk
      else
      bConcON = 0._rk
      end if

      if(bsOPC > self%bmPC .and. sDOMP > 0._rk) then
      bConsOP = bsOPC/self%bmPC*2._rk
      else
      bConsOP = 0._rk
      end if

      if(bcOPC > self%bmPC .and. cDOMP > 0._rk) then
      bConcOP = bcOPC/self%bmPC*2._rk
      else
      bConcOP = 0._rk
      end if

      ! Regeneration rate of NH4 and DIP by bacterial population (ugX/L/day): 
      bregN = bBR*bC*self%bNC + bConsON + bConcON
      bregP = bBR*bC*self%bPC + bConsOP + bConcOP

      ! Nitrogen and phosphorous biomass of bacterial population (ugX/L):
      bN = bC*self%bNC
      bP = bC*self%bPC

      ! Total gross intake of C to bacteria (ugC/L/day):
      bCup = bC*(bOCuse + bsOMu + bcOMu)

      ! Bacterial population DIC output (ugC/L/day):
      bICout = bC*bCRes

      ! Set ODEs:
      _ADD_SOURCE_(self%id_bC,    + bCup - bICout)

      _ADD_SOURCE_(self%id_NH4,   - bNH4up + bregN)
      _ADD_SOURCE_(self%id_DIP,   - bDIPup + bregP)
      _ADD_SOURCE_(self%id_DIC,   + bICout)
      _ADD_SOURCE_(self%id_DOC,   - bDOCup)
      _ADD_SOURCE_(self%id_sDOMC, - bsOMCup)
      _ADD_SOURCE_(self%id_sDOMN, - bsOMNup - bConsON)
      _ADD_SOURCE_(self%id_sDOMP, - bsOMPup - bConsOP)

      ! Set diagnostic variables:
      !_SET_DIAGNOSTIC_(self%id_bN, bN)
      !_SET_DIAGNOSTIC_(self%id_bP, bP)
      _SET_DIAGNOSTIC_(self%id_bOCV, bOCV)
      _SET_DIAGNOSTIC_(self%id_bsOMCV, bsOMCV)
      _SET_DIAGNOSTIC_(self%id_bbrelN, bbrelN)
      _SET_DIAGNOSTIC_(self%id_bbrelP, bbrelP)
      _SET_DIAGNOSTIC_(self%id_bINV, bINV)
      _SET_DIAGNOSTIC_(self%id_bsOMNV, bsOMNV)
      _SET_DIAGNOSTIC_(self%id_bVm, bVm)
      _SET_DIAGNOSTIC_(self%id_btotCV, btotCV)
      _SET_DIAGNOSTIC_(self%id_btotNV, btotNV)
      _SET_DIAGNOSTIC_(self%id_btotPV, btotPV)
      _SET_DIAGNOSTIC_(self%id_bgroCu, bgroCu)
      _SET_DIAGNOSTIC_(self%id_bIPV, bIPV)
      _SET_DIAGNOSTIC_(self%id_bsOMPV, bsOMPV)
      _SET_DIAGNOSTIC_(self%id_bCup, bCup)
      _SET_DIAGNOSTIC_(self%id_bICout, bICout)

      ! Leave spatial loops (if any)
      _LOOP_END_

  end subroutine do

end module su_bacteria

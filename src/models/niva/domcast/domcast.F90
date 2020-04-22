#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: niva_domcast
!
! !INTERFACE:
   module niva_domcast
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   implicit none
!
!  default: all is private.
   private
   ! !REVISION HISTORY:
!  Original author(s): Members of the PROGNOS project
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_domcast
      ! Variable identifiers
      type (type_state_variable_id)        :: id_DOMa, id_DOMb
      type (type_state_variable_id)        :: id_dic
      type (type_dependency_id)            :: id_swf,id_temp !shortwave flux and temperature at depth
      type (type_horizontal_dependency_id) :: id_I_0
      type (type_state_variable_id)        :: id_eo2, id_enn, id_edd, id_epp, id_ebb, id_eff, id_eaa  !ERGOM dependencies
      type (type_diagnostic_variable_id)   :: id_rate11, id_rate12, id_rate21, id_rate22, id_rate3, id_rate5a, id_rate5b, id_light !reaction rates
      type (type_state_variable_id)        :: id_Kd_DOC  !light extinction coefficient

      ! Model parameters
      ! Half saturation coefficients
      real(rk) :: Km_O2, Km_NO3
      ! Inhibition coefficients
      real(rk) :: Kin_O2
      ! Rate constants
      real(rk) :: kOM1, kOM2
      ! Physical constants
      real(rk) :: oc_DOM, qy_DOM, f_par, e_par, theta, k_floc, frac
      ! Carbon to nitrogen ratio
      real(rk) :: rfc 
   contains
      procedure :: initialize
      procedure :: do
      !procedure :: get_light_extinction
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the CDOM model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
   
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_niva_domcast), intent(inout), target :: self
   integer,                  intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Members of the PROGNOS project
!
! !LOCAL VARIABLES:
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk), parameter :: y_per_s = d_per_s/365.25_rk
   character(len=64)   :: dic_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Model parameters
   call self%get_parameter(self%Km_O2,  'Km_O2',  'mmol m-3',     'respiration',                           default=1.23e-2_rk)
   call self%get_parameter(self%Km_NO3, 'Km_NO3', 'mmol m-3',     'denitrification',                       default=1.e-2_rk)
   call self%get_parameter(self%Kin_O2, 'Kin_O2', 'mmol m-3',     'inhibitation of denitrification by O2', default=0.33_rk)
   call self%get_parameter(self%kOM1,   'kOM1',   'year-1',       'OM1 degradation',                       default=1.0_rk,      scale_factor=y_per_s)
   call self%get_parameter(self%kOM2,   'kOM2',   'year-1',       'OM2 degradation',                       default=0.1_rk,      scale_factor=y_per_s)
   call self%get_parameter(self%oc_DOM, 'oc_DOM', 'm2.mg/DOM',    'Optical cross section of DOM',          default=0.01_rk)
   call self%get_parameter(self%qy_DOM, 'qy_DOM', 'unitless',     'Quantum yield',                         default=0.1_rk)
   call self%get_parameter(self%f_par,  'f_par',  'unitless',     'PAR fraction',                          default=0.45_rk)
   call self%get_parameter(self%e_par,  'e_par',  'J mol-1',      'Energy of PAR photons',                 default=240800._rk)
   call self%get_parameter(self%theta,  'theta',  'J mol-1',      'Temperature adjustment coefficient',    default=1.047_rk)
   call self%get_parameter(self%k_floc, 'k_floc', 'mmol m-3 d-1', 'floculation rate',                      default=6.8e-9_rk)
   call self%get_parameter(self%rfc,    'rfc',    '-',            'carbon | nitrogen ratio',               default=6.625_rk)
   call self%get_parameter(self%frac,  'frac',    '-',            'Fraction of DOMb that is refractory',   default=0.0_rk)
   ! Register state variables
   call self%register_state_variable(self%id_DOMa, 'DOMa', 'mg m**-3', 'labile',     0.0_rk, minimum=0.0_rk)
   call self%register_state_variable(self%id_DOMb, 'DOMb', 'mg m**-3', 'semi-labile', 0.0_rk, minimum=0.0_rk)
   
   ! Register dependencies on external state variables
   call self%register_state_dependency(self%id_enn,'enn','mmol n/m**3', 'nitrate')
   call self%register_state_dependency(self%id_eo2,'eo2','mmol o2/m**3','oxygen')
   call self%register_state_dependency(self%id_edd,'edd','mmol n/m**3', 'detritus')
   call self%register_state_dependency(self%id_epp,'epp','mmol n/m**3', 'diatoms')
   call self%register_state_dependency(self%id_ebb,'ebb','mmol n/m**3', 'cyanobacteria')
   call self%register_state_dependency(self%id_eff,'eff','mmol n/m**3', 'flagellates')
   call self%register_state_dependency(self%id_eaa,'eaa','mmol n/m**3', 'ammonium')
   
   ! register diagnostic variable -> rates included only for debugging purposes
   call self%register_diagnostic_variable(self%id_rate11, 'rate11', '----', 'DOMa_O2')
   call self%register_diagnostic_variable(self%id_rate12, 'rate12', '----', 'DOMb_O2')
   call self%register_diagnostic_variable(self%id_rate21, 'rate21', '----', 'DOMa_N')
   call self%register_diagnostic_variable(self%id_rate22, 'rate22', '----', 'DOMb_N')
   call self%register_diagnostic_variable(self%id_rate3,  'rate3',  '----', 'PhotoOx')
   call self%register_diagnostic_variable(self%id_rate5a, 'rate5a', '----', 'DOMa_floc')
   call self%register_diagnostic_variable(self%id_rate5b, 'rate5b', '----', 'DOMb_floc')
   call self%register_diagnostic_variable(self%id_light,  'light', '----', 'light')
#if 0
   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   call self%register_state_dependency(self%id_dic,'dic','mmol m-3','total dissolved inorganic carbon',required=.false.)
   if (dic_variable/='') call self%request_coupling(self%id_dic,dic_variable)
#endif
   ! Register diagnostic variables

   ! Register environmental dependencies
   call self%register_dependency(self%id_swf,  standard_variables%downwelling_shortwave_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)

   return
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_niva_domcast),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: DOMa, DOMb 
   real(rk)                   :: swfz, temperature, dummy !dummy is a temporary variable for the rate computations
   real(rk)                   :: R11, R12, R21, R22, R3, R5a, R5b
   real(rk)                   :: POMa !Sum of cyanobacteria, diatoms and flagellates from ERGOM
   
   real(rk) :: NO3,O2,POMb,cyanobacteria,diatoms,flagellates,ammonium
   integer :: cnt = 0
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
   cnt = cnt + 1

   ! Retrieve current (local) state variable values.
   _GET_(self%id_DOMa, DOMa) ! Labile DOM
   _GET_(self%id_DOMb, DOMb) ! Semi-labile DOM

#if 1
   ! Retrieve current environmental conditions.
   _GET_(self%id_swf,swfz)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temperature)    ! temperature
#endif
   
   !retrieve domcast variables
   _GET_(self%id_enn,NO3)
   _GET_(self%id_eo2,O2)
   _GET_(self%id_edd,POMb)
   _GET_(self%id_ebb,cyanobacteria)
   _GET_(self%id_epp,diatoms)
   _GET_(self%id_eff,flagellates)
   _GET_(self%id_eaa,ammonium)
   
   ! for debugging
   !   if (cnt < 5) then
   !       Print *, cnt, DOMa, DOMb, swfz
   !   end if
   DOMa = self%frac * DOMb
   DOMb = (1 - self%frac) * DOMb
   POMa = (cyanobacteria + diatoms + flagellates) * self%rfc !Not currently used.    
   
   temperature = self%theta ** (temperature - 20.0_rk) !temperature adjustment factor  

   ! DOC mineralization by bacteria in the water colum
   dummy = O2 / (self%Km_O2 + O2) *  temperature
   R11 = self%kOM1 * dummy * DOMa / 24.0_rk
   R12 = self%kOM2 * dummy * DOMb / 24.0_rk  

   ! DOC mineralization by bacteria in the water colum
   dummy = NO3 / (self%Km_NO3 + NO3) * self%Kin_O2 / (self%Kin_O2 + O2) * temperature
   R21 = self%kOM1 * dummy * DOMa / 24.0_rk
   R22 = self%kOM2 * dummy * DOMb / 24.0_rk

   !Photo-oxidation and photo-mineralization
   R3 = +1.0 * self%oc_DOM * self%qy_DOM * self%f_par * 1._rk/self%e_par * swfz !R3,4

   !Flocculation
   R5a = self%k_floc * DOMa
   R5b = self%k_floc * DOMb * swfz
   
   ! All processes degrade DOM pools
   _ADD_SOURCE_(self%id_DOMa, -((R11 + R21) * 24.0_rk + R3 + R5a) )
   _ADD_SOURCE_(self%id_DOMb, -((R12 + R22) * 24.0_rk + R3 + R5b) )

   !O2 is consumed for DOM bacteria mineralization of both DOM pools
   _ADD_SOURCE_(self%id_eo2,  -(R11+R12) )

   !Nitrate is consumed for DOM bacteria mineralization of both DOM pools
   _ADD_SOURCE_(self%id_enn,  -(R21+R22)*0.8_rk)

   !Ammonium is created during bacteria mineralization
   _ADD_SOURCE_(self%id_eaa,  (R21+R22)*0.8_rk)
   
   !POMb is produced during flocculation (add to detritus in selma)
   _ADD_SOURCE_(self%id_edd,  (R5a / 5.6_rk + R5b / 10.0_rk ))
   
   ! Export diagnostic variables -> rates included only for debugging purposes
   _SET_DIAGNOSTIC_(self%id_rate11, R11)
   _SET_DIAGNOSTIC_(self%id_rate12, R12)
   _SET_DIAGNOSTIC_(self%id_rate21, R21)
   _SET_DIAGNOSTIC_(self%id_rate22, R22)
    !WRITE(UNIT=1, FMT=*) cnt, R3
   _SET_DIAGNOSTIC_(self%id_rate3,  R3)
   _SET_DIAGNOSTIC_(self%id_rate5a, R5a)
   _SET_DIAGNOSTIC_(self%id_rate5b, R5b)
   _SET_DIAGNOSTIC_(self%id_light, swfz)
 
   _LOOP_END_
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
   class (type_niva_domcast), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Members of the PROGNOS project
!
! !LOCAL VARIABLES:
   real(rk)                     :: DOMb
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_DOMb, DOMb) !
   
   ! Self-shading with explicit contribution colored dissolved organic carbon
   _SET_EXTINCTION_(0.1_rk*DOMb**1.22_rk)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
!   pure real(rk) function fpz(self,p,z)
!!
!! !DESCRIPTION:
!! Here, the classical Ivlev formulation for zooplankton grazing on
!! phytoplankton is formulated.
!!
!! !INPUT PARAMETERS:
!   type (type_gotm_npzd), intent(in) :: self
!   real(rk), intent(in)         :: p,z
!!
!! !REVISION HISTORY:
!!  Original author(s): Members of the PROGNOS project
!!
!!EOP
!!-----------------------------------------------------------------------
!!BOC
!   fpz = self%gmax*(1.0_rk-exp(-self%iv*self%iv*p*p))*(z+self%z0)
!
!   end function fpz
!EOC

!-----------------------------------------------------------------------

   end module niva_domcast

!-----------------------------------------------------------------------
! Copyright by the PROGNOS team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

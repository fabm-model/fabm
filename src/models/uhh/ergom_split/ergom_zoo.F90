#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_ergom_split_diatoms - ERGOM zooplankton compartment
!
! !INTERFACE:
   module fabm_uhh_ergom_split_zoo
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   use fabm_uhh_ergom_split_utilities

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_ergom_split_zoo
!     Variable identifiers
      type (type_state_variable_id)        :: id_zoo
      type (type_state_variable_id)        :: id_dia,id_din
      type (type_state_variable_id)        :: id_detritus
      type (type_state_variable_id)        :: id_ammonium, id_phosphate
      type (type_state_variable_id)        :: id_oxygen
      type (type_dependency_id)            :: id_temp
      type (type_diagnostic_variable_id)   :: id_secprod_total  
      type (type_diagnostic_variable_id)   :: id_secprod_dia    
      type (type_diagnostic_variable_id)   :: id_secprod_din 

!     Model parameters
      real(rk) :: topt
      real(rk) :: gmax_dia  ! maximum grazing rate diatoms
      real(rk) :: gmax_din  ! maximum grazing rate dinoflagellates
      real(rk) :: p_to_n    ! P:N ratio of zooplankton
      real(rk) :: zoo0      ! background concentration
      real(rk) :: alpha
      real(rk) :: iv
      real(rk) :: s2        ! stoichiometric factor
      real(rk) :: s3        ! stoichiometric factor
      real(rk) :: lza
      real(rk) :: lzd
      real(rk) :: slopf     ! sloppy feeding factor [0..1]
      logical  :: graz_dia, graz_din 
      logical  :: use_oxy, use_pho

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the ERGOM diatoms model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the uhh_ergom_split_diatoms namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_uhh_ergom_split_zoo), intent(inout), target :: self
   integer,                        intent(in)            :: configunit

!  LOCAL VARIABLES:
   real(rk)                  :: zoo0 !background_concentration
   real(rk)                  :: topt
   real(rk)                  :: gmax_dia  ! maximum grazing rate diatoms
   real(rk)                  :: gmax_din  ! maximum grazing rate dinoflagellates
   real(rk)                  :: p_to_n    ! P:N ratio of zooplankton
   real(rk)                  :: excretion_rate ! lza
   real(rk)                  :: mortality_rate ! lzd
   real(rk)                  :: s2, s3, alpha, iv
   real(rk)                  :: slopf
   character(len=64)         :: oxygen_variable
   character(len=64)         :: diatoms_variable
   character(len=64)         :: dinoflagellates_variable
   character(len=64)         :: ammonium_variable
   character(len=64)         :: phosphate_variable
   character(len=64)         :: detritus_variable

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   namelist /uhh_ergom_split_zoo/ zoo0, iv, &
             gmax_dia, gmax_din, p_to_n, s2, &
             s3, topt,  &
             mortality_rate, excretion_rate, slopf, &
             ammonium_variable,& 
             phosphate_variable, detritus_variable, &
             oxygen_variable, diatoms_variable, &
             dinoflagellates_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   zoo0 = 0.0225_rk
   gmax_dia     = -1.0_rk
   gmax_din     = 0.03_rk
   alpha        = 0.3_rk
   p_to_n       = 0.0625_rk
   s2           = 6.625_rk
   s3           = 8.625_rk
   iv           = 0.24444444_rk
   topt         = 20.0_rk
   mortality_rate = 0.02_rk ! lzd [1/d]
   excretion_rate = 0.01_rk ! lza [1/d]
   slopf        = 1.0_rk ! default: no sloppy feeding
   ammonium_variable = 'uhh_ergom_split_base_amm'
   phosphate_variable = 'uhh_ergom_split_base_pho'
   detritus_variable = 'uhh_ergom_split_base_det'
   oxygen_variable = 'uhh_ergom_split_base_oxy'
   dinoflagellates_variable= 'uhh_dinoflag_veg'
   diatoms_variable= 'uhh_diatoms_veg'
   
   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_ergom_split_zoo,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%zoo0,   'zoo0', default=zoo0)
   call self%get_parameter(self%p_to_n,   'p_to_n',   default=p_to_n)
   call self%get_parameter(self%s2,       's2',   default=s2)
   call self%get_parameter(self%s3,       's3',   default=s3)
   call self%get_parameter(self%iv,       'iv',   default=iv)
   call self%get_parameter(self%slopf,    'slopf',default=slopf)
   call self%get_parameter(self%topt,     'topt', default=topt)
   call self%get_parameter(self%gmax_dia, 'gmax_dia', default=gmax_dia, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%gmax_din, 'gmax_din', default=gmax_din, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%lza,    'excretion_rate', default=excretion_rate,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%lzd,    'mortality_rate', default=mortality_rate,  scale_factor=1.0_rk/secs_pr_day)
   
   ! Register state variables
   call self%register_state_variable(self%id_zoo,'zoo','mmol n/m**3','zooplankton')
   
   ! Register dependencies on external standard variables
   call self%register_state_dependency(self%id_ammonium, 'ammonium_target', 'mmol/m**3','sink for excreted matter and ammonium source')
   self%use_pho = (phosphate_variable /= '')
   if (self%use_pho) &
     call self%register_state_dependency(self%id_phosphate, 'phosphate_target',  'mmol/m**3','phosphate source')

   
   ! Register external state dependencies
   call self%register_state_dependency(self%id_detritus, 'mortality_target','mmol/m**3','sink for dead matter')
   self%use_oxy = (oxygen_variable /= '')
   if (self%use_oxy) &
     call self%register_state_dependency(self%id_oxygen,   'oxygen_target'   ,'mmol-O2/m**3','dissolved oxygen pool')
   self%graz_dia = (diatoms_variable /= '')
   self%graz_din = (dinoflagellates_variable /= '')


   if (self%graz_dia) &
      call self%register_state_dependency(self%id_dia, 'diatoms_target','mmol/m**3','diatom prey')
   if (self%graz_din) &
      call self%register_state_dependency(self%id_din, 'dinoflagellates_target','mmol/m**3','dinoflagellate prey')

   call self%request_coupling(self%id_ammonium, ammonium_variable)
   if (self%use_pho) &
      call self%request_coupling(self%id_phosphate, phosphate_variable)
   call self%request_coupling(self%id_detritus,detritus_variable)
   if (self%use_oxy) &
      call self%request_coupling(self%id_oxygen, oxygen_variable)
   if (self%graz_dia) &
      call self%request_coupling(self%id_dia,diatoms_variable)
   if (self%graz_din) &
      call self%request_coupling(self%id_din,dinoflagellates_variable)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

 ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_secprod_total,'secprod_total','mmol-N/m**3/d', &
         'total secondary production rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_secprod_dia,'secprod_dia','mmol-N/m**3/d', &
         'dia secondary production rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_secprod_din,'secprod_din','mmol-N/m**3/d', &
         'din secondary production rate', output=output_instantaneous)
   return

99 call self%fatal_error('fabm_uhh_ergom_split_zoo','Error reading namelist uhh_ergom_split_zoo')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of ERGOM-split zooplankton model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_uhh_ergom_split_zoo), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: dia,din,zoo,temp
   real(rk)                   :: psum, gross_zoo
   real(rk)                   :: mortality, excretion
   real(rk)                   :: graz_dia,graz_din
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   dia=0.0_rk; graz_dia=0.0_rk
   din=0.0_rk; graz_din=0.0_rk
   _GET_(self%id_zoo,zoo)                    ! zooplankton
   if (self%graz_dia) _GET_(self%id_dia,dia) ! diatoms
   if (self%graz_din) _GET_(self%id_din,din) ! dinoflagellates

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)      ! local temperature
   
   psum = dia  + din  + self%zoo0
   gross_zoo = zoo + self%zoo0
   if (self%graz_dia) then
      graz_dia = fpz(self,self%gmax_dia,temp,self%topt,psum) * &
                 dia/psum * gross_zoo
      _ADD_SOURCE_(self%id_dia,-graz_dia)
   end if
   
   if (self%graz_din) then
      graz_din = fpz(self,self%gmax_din,temp,self%topt,psum) * &
                 din/psum * gross_zoo
      _ADD_SOURCE_(self%id_din,-graz_din)   
   end if
   excretion = self%lza * zoo * gross_zoo
   mortality = self%lzd * zoo * gross_zoo

   ! Set temporal derivatives
   _ADD_SOURCE_(self%id_zoo,self%slopf*(graz_dia + graz_din ) - mortality - excretion)
   _ADD_SOURCE_(self%id_ammonium,excretion)
   _ADD_SOURCE_(self%id_detritus,mortality + (1.0_rk-self%slopf)*(graz_dia  + graz_din ))
   if (self%use_pho) then
   _ADD_SOURCE_(self%id_phosphate, self%p_to_n * excretion)
   end if
   if (self%use_oxy) then
     _ADD_SOURCE_(self%id_oxygen, -self%s2 * excretion)
   end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_secprod_total,(self%slopf*(graz_dia  + graz_din ))*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_secprod_dia,self%slopf*graz_dia*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_secprod_din,self%slopf*graz_din*secs_pr_day)
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: PP/DD of ERGOM-split zooplankton model
!
! !INTERFACE:
   subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!
! !INPUT PARAMETERS:
   class (type_uhh_ergom_split_zoo), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
! !LOCAL VARIABLES:
   real(rk)                   :: dia,din,zoo,temp
   real(rk)                   :: psum, gross_zoo
   real(rk)                   :: mortality, excretion
   real(rk)                   :: graz_dia,graz_din,graz_cya
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   !! Retrieve current (local) state variable values.
   dia=0.0_rk; graz_dia=0.0_rk
   din=0.0_rk; graz_din=0.0_rk
   _GET_(self%id_zoo,zoo)                    ! zooplankton
   if (self%graz_dia) _GET_(self%id_dia,dia) ! diatoms
   if (self%graz_din) _GET_(self%id_din,din) ! dinoflagellates

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)      ! local temperature
   
   psum = dia  + din  + self%zoo0
   gross_zoo = zoo + self%zoo0
   if (self%graz_dia) then
      graz_dia = fpz(self,self%gmax_dia,temp,self%topt,psum) * &
                 dia/psum * gross_zoo
      _SET_DD_SYM_(self%id_dia,self%id_zoo,graz_dia)
   end if
   if (self%graz_din) then
      graz_din = fpz(self,self%gmax_din,temp,self%topt,psum) * &
                 din/psum * gross_zoo
      _SET_DD_SYM_(self%id_din,self%id_zoo,graz_din)
   end if
   excretion = self%lza * zoo * gross_zoo
   mortality = self%lzd * zoo * gross_zoo

   ! Set temporal derivatives
   if (self%use_pho) then
   _SET_PP_(self%id_phosphate,self%id_phosphate, self%p_to_n * excretion)
   end if
   if (self%use_oxy) then
     _SET_DD_(self%id_oxygen,self%id_oxygen, self%s2 * excretion)
   end if
   _SET_DD_SYM_(self%id_zoo,self%id_detritus,mortality + (1.0_rk-self%slopf)*( graz_din  + graz_dia))
   _SET_DD_SYM_(self%id_zoo,self%id_ammonium,excretion)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_secprod_total,(self%slopf*(graz_dia  + graz_din ))*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_secprod_dia,self%slopf*graz_dia*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_secprod_din,self%slopf*graz_din*secs_pr_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure real(rk) function fpz(self,g,t,topt,psum)
!
! !DESCRIPTION:
! The Ivlev formulation for zooplankton grazing on the three phytoplankton
! classes $c_1$, $c_2$, and $c_3$ is given here as a function:
! \begin{equation}\label{neu_di4}
! d_{i,4}=g_i^{\max}\left(1+\frac{T^2}{T_{opt}^2}\exp
! \left(1-\frac{2T}{T_{opt}} \right)\right)
! \left( 1-\exp\left(-I_v^2 \left( \sum_{
! j=1}^3c_j \right)^2\right)  \right)
! \frac{c_i}{\sum_{j=1}^3c_j}\left( c_4+c_4^{\min} \right)
! \end{equation}
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   class (type_uhh_ergom_split_zoo), intent(in) :: self
   real(rk), intent(in)                                 :: g,t,topt,psum
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz=g*(1.0_rk+t**2/topt**2*exp(1.0_rk-2.0_rk*t/topt))*               &
        (1.0_rk-exp(-self%iv**2*psum**2))
   return
   end function fpz
!EOC

   end module fabm_uhh_ergom_split_zoo

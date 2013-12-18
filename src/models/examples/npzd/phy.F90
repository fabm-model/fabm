#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_npzd_phy - Fennel & Neumann 1996 NPZD model - phytoplankton component
!
! !INTERFACE:
   module fabm_examples_npzd_phy
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_examples_npzd_phy
!     Variable identifiers
      type (type_state_variable_id)        :: id_p
      type (type_state_variable_id)        :: id_exctarget,id_morttarget,id_upttarget
      type (type_dependency_id)            :: id_par
      type (type_horizontal_dependency_id) :: id_I_0
      type (type_diagnostic_variable_id)   :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR

!     Model parameters
      real(rk) :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rpdu,rpdl
      real(rk) :: dic_per_n

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_ppdd
      procedure :: get_light_extinction
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples_npzd_phy namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_phy), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk)                  :: p_initial
   real(rk)                  :: p0
   real(rk)                  :: w_p
   real(rk)                  :: kc
   real(rk)                  :: i_min
   real(rk)                  :: rmax
   real(rk)                  :: alpha
   real(rk)                  :: rpn
   real(rk)                  :: rpdu
   real(rk)                  :: rpdl
   character(len=64)         :: excretion_target_variable
   character(len=64)         :: mortality_target_variable
   character(len=64)         :: uptake_target_variable

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   namelist /examples_npzd_phy/ &
             p_initial,p0,w_p,kc,i_min,rmax,alpha,rpn,rpdu,rpdl,     &
             excretion_target_variable,mortality_target_variable, &
             uptake_target_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   p_initial = 0.0_rk
   p0        = 0.0225_rk
   w_p       = -1.0_rk
   kc        = 0.03_rk
   i_min     = 25.0_rk
   rmax      = 1.0_rk
   alpha     = 0.3_rk
   rpn       = 0.01_rk
   rpdu      = 0.02_rk
   rpdl      = 0.1_rk
   excretion_target_variable = ''
   mortality_target_variable = ''
   uptake_target_variable    = ''

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=examples_npzd_phy,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%p0,   'p0',   default=p0)
   call self%get_parameter(self%kc,   'kc',   default=kc)
   call self%get_parameter(self%i_min,'i_min',default=i_min)
   call self%get_parameter(self%rmax, 'rmax', default=rmax, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%alpha,'alpha',default=alpha)
   call self%get_parameter(self%rpn,  'rpn',  default=rpn,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%rpdu, 'rpdu', default=rpdu, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%rpdl, 'rpdl', default=rpdl, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(w_p,       'w_p',  default=w_p,  scale_factor=1.0_rk/secs_pr_day)

   ! Register state variables
   call self%register_state_variable(self%id_p,'phy','mmol/m**3','phytoplankton', &
                                p_initial,minimum=0.0_rk,vertical_movement=w_p)

   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_p)

   ! Register dependencies on external state variables
   call self%register_state_dependency(self%id_exctarget, 'excretion_target','mmol/m**3','sink for excreted matter')
   call self%register_state_dependency(self%id_morttarget,'mortality_target','mmol/m**3','sink for dead matter')
   call self%register_state_dependency(self%id_upttarget, 'uptake_target',   'mmol/m**3','nutrient source')

   ! Automatically couple dependencies if target variables have been specified.
   if (excretion_target_variable/='') call self%request_coupling(self%id_exctarget, excretion_target_variable)
   if (mortality_target_variable/='') call self%request_coupling(self%id_morttarget,mortality_target_variable)
   if (uptake_target_variable   /='') call self%request_coupling(self%id_upttarget, uptake_target_variable)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_GPP, 'GPP','mmol/m**3',  'gross primary production',           &
                                     time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_NCP, 'NCP','mmol/m**3',  'net community production',           &
                                     time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_PPR, 'PPR','mmol/m**3/d','gross primary production rate',      &
                                     time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_NPR, 'NPR','mmol/m**3/d','net community production rate',      &
                                     time_treatment=time_treatment_averaged)
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2',     'photosynthetically active radiation',&
                                     time_treatment=time_treatment_averaged)

   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

   return

99 call self%fatal_error('fabm_examples_npzd_phy','Error reading namelist examples_npzd_phy')

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
   class (type_examples_npzd_phy), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: n,p,par,I_0
   real(rk)                   :: iopt,rpd,primprod
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p)         ! phytoplankton
   _GET_(self%id_upttarget,n) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,self%I_min)

   ! Loss rate of phytoplankton to detritus depends on local light intensity.
   if (par .ge. self%I_min) then
      rpd = self%rpdu
   else
      rpd = self%rpdl
   end if

   ! Define some intermediate quantities that will be reused multiple times.
   primprod = fnp(self,n,p,par,iopt)

   ! Set temporal derivatives
   _SET_ODE_(self%id_p,primprod - self%rpn*p - rpd*p)

   ! If an externally maintained ...
   _SET_ODE_(self%id_upttarget,-primprod)
   _SET_ODE_(self%id_morttarget,rpd*p)
   _SET_ODE_(self%id_exctarget,self%rpn*p)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP ,primprod)
   _SET_DIAGNOSTIC_(self%id_NCP ,primprod - self%rpn*p)
   _SET_DIAGNOSTIC_(self%id_PPR ,primprod*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR ,(primprod - self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
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
   class (type_examples_npzd_phy), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: p
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*(self%p0+p))

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_phy), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                   :: n,p,par,I_0
   real(rk)                   :: iopt,rpd,primprod
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p)         ! phytoplankton
   _GET_(self%id_upttarget,n) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)     ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,self%I_min)

   ! Loss rate of phytoplankton to detritus depends on local light intensity.
   if (par .ge. self%I_min) then
      rpd = self%rpdu
   else
      rpd = self%rpdl
   end if

   ! Rate of primary production will be reused multiple times - calculate it once.
   primprod = fnp(self,n,p,par,iopt)

   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_ [as opposed to _SET_DD_], assignments to dd(i,j)
   ! are automatically assigned to pp(j,i) as well.
   _SET_DD_SYM_(self%id_upttarget,self%id_p,primprod)
   _SET_DD_SYM_(self%id_p,self%id_exctarget,self%rpn*p)
   _SET_DD_SYM_(self%id_p,self%id_morttarget,rpd*p)

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP,primprod)
   _SET_DIAGNOSTIC_(self%id_NCP,primprod-self%rpn*p)
   _SET_DIAGNOSTIC_(self%id_PPR,primprod*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR,(primprod-self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   pure real(rk) function fnp(self,n,p,par,iopt)
!
! !DESCRIPTION:
! Here, the classical Michaelis-Menten formulation for nutrient uptake
! is formulated.
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_phy), intent(in) :: self
   real(rk), intent(in)                       :: n,p,par,iopt
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fnp = self%rmax*par/iopt*exp(1.0_rk-par/iopt)*n/(self%alpha+n)*(p+self%p0)

   end function fnp
!EOC

!-----------------------------------------------------------------------

   end module fabm_examples_npzd_phy

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

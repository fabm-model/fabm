#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: examples_npzd_phy - Fennel & Neumann 1996 NPZD model - phytoplankton component
!
! !INTERFACE:
   module examples_npzd_phy
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
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk)            :: w_p
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day and are converted here to values per second.
   call self%get_parameter(self%p0,   'p0',   'mmol m-3', 'background concentration ',               default=0.0225_rk)
   call self%get_parameter(self%kc,   'kc',   'm2 mmol-1','specific light extinction',               default=0.03_rk)
   call self%get_parameter(self%i_min,'i_min','W m-2',    'minimum light intensity in euphotic zone',default=25.0_rk)
   call self%get_parameter(self%rmax, 'rmax', 'd-1',      'maximum specific growth rate',            default=1.0_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%alpha,'alpha','mmol m-3', 'half-saturation nutrient concentration',  default=0.3_rk)
   call self%get_parameter(self%rpn,  'rpn',  'd-1',      'excretion rate',                          default=0.01_rk, scale_factor=d_per_s)
   call self%get_parameter(self%rpdu, 'rpdu', 'd-1',      'mortality in euphotic zone',              default=0.02_rk, scale_factor=d_per_s)
   call self%get_parameter(self%rpdl, 'rpdl', 'd-1',      'mortality below euphotic zone',           default=0.1_rk,  scale_factor=d_per_s)
   call self%get_parameter(w_p,       'w_p',  'm d-1',    'vertical velocity (<0 for sinking)',      default=-1.0_rk, scale_factor=d_per_s)

   ! Register state variables
   call self%register_state_variable(self%id_p,'c','mmol m-3','concentration',0.0_rk,minimum=0.0_rk,vertical_movement=w_p)

   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_p)

   ! Register dependencies on external state variables
   call self%register_state_dependency(self%id_upttarget, 'uptake_target',   'mmol m-3','nutrient source')
   call self%register_state_dependency(self%id_exctarget, 'excretion_target','mmol m-3','sink for excreted matter')
   call self%register_state_dependency(self%id_morttarget,'mortality_target','mmol m-3','sink for dead matter')

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_GPP, 'GPP','mmol m-3',    'gross primary production',           &
                                     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_NCP, 'NCP','mmol m-3',    'net community production',           &
                                     output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_PPR, 'PPR','mmol m-3 d-1','gross primary production rate',      &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_NPR, 'NPR','mmol m-3 d-1','net community production rate',      &
                                     output=output_time_step_averaged)
   call self%register_diagnostic_variable(self%id_dPAR,'PAR','W m-2',       'photosynthetically active radiation',&
                                     output=output_time_step_averaged)

   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

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
   if (par>=self%I_min) then
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
   if (par>=self%I_min) then
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

   end module examples_npzd_phy

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

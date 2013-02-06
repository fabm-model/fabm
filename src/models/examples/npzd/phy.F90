#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_npzd_phy --- phytoplankton biogeochemical model
!
! !INTERFACE:
   module fabm_examples_npzd_phy
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   use fabm_driver

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_examples_npzd_phy, &
          examples_npzd_phy_init, &
          examples_npzd_phy_do,   &
          aed_p_do_ppdd, &
          aed_p_get_light_extinction, &
          examples_npzd_phy_get_conserved_quantities
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type type_examples_npzd_phy
!     Variable identifiers
      type (type_state_variable_id)        :: id_p
      type (type_state_variable_id)        :: id_exctarget,id_morttarget,id_upttarget
      type (type_dependency_id)            :: id_par
      type (type_horizontal_dependency_id) :: id_I_0
      type (type_diagnostic_variable_id)   :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR
      type (type_conserved_quantity_id)    :: id_totN

!     Model parameters
      real(rk) :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      real(rk) :: dic_per_n
      logical  :: do_exc,do_mort,do_upt
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
   subroutine examples_npzd_phy_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_phy), intent(out)   :: self
   _CLASS_ (type_model_info),     intent(inout) :: modelinfo
   integer,                       intent(in)    :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                  :: p_initial=0.
   real(rk)                  :: p0=0.0225
   real(rk)                  :: w_p=-1.157407e-05
   real(rk)                  :: i_min=25.
   real(rk)                  :: rmax=1.157407e-05
   real(rk)                  :: alpha=0.3
   real(rk)                  :: rpn=1.157407e-07
   real(rk)                  :: rpdu=2.314814e-07
   real(rk)                  :: rpdl=1.157407e-06
   character(len=64)         :: excretion_target_variable=''
   character(len=64)         :: mortality_target_variable=''
   character(len=64)         :: uptake_target_variable=''

   real(rk), parameter :: secs_pr_day = 86400.
   namelist /examples_npzd_phy/ &
             p_initial,p0,w_p,i_min,rmax,alpha,rpn,rpdu,rpdl,     &
             excretion_target_variable,mortality_target_variable, &
             uptake_target_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=examples_npzd_phy,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%p0    = p0
   self%i_min = i_min
   self%rmax  = rmax/secs_pr_day
   self%alpha = alpha
   self%rpn  = rpn /secs_pr_day
   self%rpdu = rpdu/secs_pr_day
   self%rpdl = rpdl/secs_pr_day

   ! Register state variables
   call register_state_variable(modelinfo,self%id_p,'phy','mmol/m**3','phytoplankton', &
                                    p_initial,minimum=_ZERO_,vertical_movement=w_p/secs_pr_day)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   self%do_exc = excretion_target_variable.ne.''
   if (self%do_exc) call register_state_dependency(modelinfo,self%id_exctarget,excretion_target_variable)
   self%do_mort = mortality_target_variable.ne.''
   if (self%do_mort) call register_state_dependency(modelinfo,self%id_morttarget,mortality_target_variable)
   self%do_upt = uptake_target_variable.ne.''
   if (self%do_upt) call register_state_dependency(modelinfo,self%id_upttarget,uptake_target_variable)


   ! Register diagnostic variables
   call register_diagnostic_variable(modelinfo,self%id_GPP,'GPP','mmol/m**3',  'gross primary production',           &
                     time_treatment=time_treatment_step_integrated)
   call register_diagnostic_variable(modelinfo,self%id_NCP,'NCP','mmol/m**3',  'net community production',           &
                     time_treatment=time_treatment_step_integrated)
   call register_diagnostic_variable(modelinfo,self%id_PPR,'PPR','mmol/m**3/d','gross primary production rate',      &
                     time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_NPR,'NPR','mmol/m**3/d','net community production rate',      &
                     time_treatment=time_treatment_averaged)
   call register_diagnostic_variable(modelinfo,self%id_dPAR,'PAR','W/m**2',     'photosynthetically active radiation',&
                     time_treatment=time_treatment_averaged)

   ! Register conserved quantities
!KB   self%id_totN = register_conserved_quantity(modelinfo,'N','mmol/m**3','nitrogen')

   ! Register environmental dependencies
   call register_dependency(modelinfo, self%id_par, varname_par)
   call register_dependency(modelinfo, self%id_I_0, varname_par_sf)

   return

99 call fatal_error('fabm_examples_npzd_phy','Error reading namelist npzd')

   end subroutine examples_npzd_phy_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine examples_npzd_phy_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_phy), intent(in)     :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                   :: n,p,par,I_0
   real(rk)                   :: iopt,rpd,primprod,dn
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_upttarget,n) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_   (self%id_par,par)  ! local photosynthetically active radiation
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
   if (self%do_upt) then
      _SET_ODE_(self%id_upttarget,-primprod)
   end if
   if (self%do_mort) then
      _SET_ODE_(self%id_morttarget,rpd*p)
   end if
   if (self%do_exc) then
      _SET_ODE_(self%id_exctarget,self%rpn*p)
   end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP ,primprod)
   _SET_DIAGNOSTIC_(self%id_NCP ,primprod - self%rpn*p)
   _SET_DIAGNOSTIC_(self%id_PPR ,primprod*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR ,(primprod - self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine examples_npzd_phy_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   pure subroutine aed_p_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_phy), intent(in)     :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                     :: p
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(0.0*p)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_p_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   pure subroutine examples_npzd_phy_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_phy), intent(in)     :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk)                     :: p
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton

   ! Total nutrient is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totN,p)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine examples_npzd_phy_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine aed_p_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_phy), intent(in)     :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   real(rk)                   :: n,p,par,I_0
   real(rk)                   :: iopt,rpd,dn,primprod
   real(rk), parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_upttarget,n) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_   (self%id_par,par)  ! local photosynthetically active radiation

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
!  _SET_DD_SYM_(self%id_p,self%id_z,fpz(self,p,z))          ! spz
!  _SET_DD_SYM_(self%id_p,self%id_n,self%rpn*p)             ! spn
!  _SET_DD_SYM_(self%id_p,self%id_d,rpd*p)                  ! spd

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
!  dn = - fnp(self,n,p,par,iopt) + self%rpn*p + self%rzn*z + self%rdn*d
!  if (self%use_dic) _SET_PP_(self%id_dic,self%id_dic,self%dic_per_n*dn)

   _SET_DD_(self%id_p,self%id_p,primprod)
   _SET_DD_(self%id_p,self%id_p,- self%rpn*p)
   _SET_DD_(self%id_p,self%id_p,- rpd*p)

   ! If an externally maintained ...
   if (self%do_upt) then
      _SET_PP_(self%id_upttarget,self%id_upttarget,-primprod)
   end if
   if (self%do_mort) then
      _SET_PP_(self%id_morttarget,self%id_morttarget,rpd*p)
   end if
   if (self%do_exc) then
      _SET_PP_(self%id_exctarget,self%id_exctarget,self%rpn*p)
   end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP,primprod)
   _SET_DIAGNOSTIC_(self%id_NCP,primprod-self%rpn*p)
   _SET_DIAGNOSTIC_(self%id_PPR,primprod*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR,(primprod-self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_p_do_ppdd
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
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_phy), intent(in)     :: self
   real(rk), intent(in)         :: n,p,par,iopt
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fnp = self%rmax*par/iopt*exp(_ONE_-par/iopt)*n/(self%alpha+n)*(p+self%p0)

   end function fnp
!EOC

!-----------------------------------------------------------------------

   end module fabm_examples_npzd_phy

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_npzd_zoo --- Zooplankton biogeochemical model
!
! !INTERFACE:
   module fabm_examples_npzd_zoo
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
   public type_examples_npzd_zoo, &
          examples_npzd_zoo_init, &
          examples_npzd_zoo_do,   &
          aed_z_do_ppdd, &
          aed_z_get_light_extinction, aed_z_get_conserved_quantities
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
! !PUBLIC DERIVED TYPES:
   type type_examples_npzd_zoo
!     Variable identifiers
      type (type_state_variable_id)      :: id_z
      type (type_state_variable_id)      :: id_exctarget,id_morttarget,id_grztarget
      type (type_dependency_id)          :: id_temp

!     Model parameters
      REALTYPE :: z0,gmax,iv,rzn,rzd
      logical  :: do_exc,do_mort,do_grz
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
   subroutine examples_npzd_zoo_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_zoo), intent(out)   :: self
   _CLASS_ (type_model_info),     intent(inout) :: modelinfo
   integer,                       intent(in)    :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                  :: z_initial=0.
   REALTYPE                  :: z0=0.0225
   REALTYPE                  :: gmax=5.787037e-06
   REALTYPE                  :: iv=1.1
   REALTYPE                  :: rzn=1.157407e-07
   REALTYPE                  :: rzd=2.314814e-07
   character(len=64)         :: excretion_target_variable=''
   character(len=64)         :: mortality_target_variable=''
   character(len=64)         :: grazing_target_variable=''

   REALTYPE, parameter :: secs_pr_day = 86400.
   namelist /examples_npzd_zoo/ &
            z_initial,z0,gmax,iv,rzn,rzd,excretion_target_variable, &
            mortality_target_variable,grazing_target_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=examples_npzd_zoo,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%z0    = z0
   self%gmax  = gmax/secs_pr_day
   self%iv    = iv
   self%rzn  = rzn /secs_pr_day
   self%rzd  = rzd /secs_pr_day

   ! Register state variables
   call register_state_variable(modelinfo,self%id_z,'zoo','mmol/m**3','zooplankton', &
                                    z_initial,minimum=_ZERO_)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   self%do_exc = excretion_target_variable.ne.''
   if (self%do_exc) call register_state_dependency(modelinfo,self%id_exctarget,excretion_target_variable)
   self%do_mort = mortality_target_variable.ne.''
   if (self%do_mort) call register_state_dependency(modelinfo,self%id_morttarget,mortality_target_variable)
   self%do_grz = grazing_target_variable.ne.''
   if (self%do_grz) call register_state_dependency(modelinfo,self%id_grztarget,grazing_target_variable)
   ! Register diagnostic variables

   ! Register conserved quantities

   ! Register environmental dependencies
   call register_dependency(modelinfo,self%id_temp,varname_temp)

   return

99 call fatal_error('examples_npzd_zoo_init','Error reading namelist aed_z')

   end subroutine examples_npzd_zoo_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine examples_npzd_zoo_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_zoo), intent(in)     :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: p,z,temp
   REALTYPE, parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_z,z) ! zooplankton
   _GET_(self%id_grztarget,p) ! phytoplankton

   ! Retrieve current environmental conditions.
   _GET_   (self%id_temp,temp)  ! local photosynthetically active radiation

   ! Loss rate of zooplankton to detritus depends on mortality


   ! Set temporal derivatives
    _SET_ODE_(self%id_z,fpz(self,p,z) - self%rzn*z - self%rzd*z)

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   if (self%do_grz) then
      _SET_ODE_(self%id_grztarget,-fpz(self,p,z))
   end if
   if (self%do_mort) then
      _SET_ODE_(self%id_morttarget,self%rzd*z)
   end if
   if (self%do_exc) then
      _SET_ODE_(self%id_exctarget,self%rzn*z)
   end if

   ! Export diagnostic variables

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine examples_npzd_zoo_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   pure subroutine aed_z_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_zoo), intent(in)     :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: p,d
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_z_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   pure subroutine aed_z_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_zoo), intent(in)     :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: n,p,z,d
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_z,z) ! zooplankton

   ! Total nutrient is simply the sum of all variables.
!#   _SET_CONSERVED_QUANTITY_(self%id_totN,n+p+z+d)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_z_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine aed_z_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_zoo), intent(in)     :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: p,z,temp
   REALTYPE                   :: iopt,rpd,dn,primprod
   REALTYPE, parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_z,z) ! zooplankton

   ! Retrieve current environmental conditions.
   _GET_   (self%id_temp,temp)  ! local photosynthetically active radiation


   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_ [as opposed to _SET_DD_], assignments to dd(i,j)
   ! are automatically assigned to pp(j,i) as well.
!  _SET_DD_SYM_(self%id_z,self%id_n,self%rzn*z)             ! szn
!  _SET_DD_SYM_(self%id_z,self%id_d,self%rzd*z)             ! szd

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
!   dn = 1.0
!   if (self%use_phy) _SET_PP_(self%id_phy,self%id_phy,self%phy_per_n*dn)

   ! Set temporal derivatives
    _SET_DD_(self%id_z,self%id_z, fpz(self,p,z) )
    _SET_DD_(self%id_z,self%id_z,-self%rzn*z    )
    _SET_DD_(self%id_z,self%id_z, self%rzd*z    )

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   if (self%do_grz) then
      _SET_PP_(self%id_grztarget,self%id_grztarget,-fpz(self,p,z))
   end if
   if (self%do_mort) then
      _SET_PP_(self%id_morttarget,self%id_grztarget,self%rzd*z)
   end if
   if (self%do_exc) then
      _SET_PP_(self%id_exctarget,self%id_grztarget,self%rzn*z)
   end if


   ! Export diagnostic variables

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_z_do_ppdd
!EOC



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure REALTYPE function fpz(self,p,z)
!
! !DESCRIPTION:
! Here, the classical Ivlev formulation for zooplankton grazing on
! phytoplankton is formulated.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_zoo), intent(in)     :: self
   REALTYPE, intent(in)          :: p,z
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz = self%gmax*(_ONE_-exp(-self%iv*self%iv*p*p))*(z+self%z0)

   end function fpz
!EOC
!-----------------------------------------------------------------------

   end module fabm_examples_npzd_zoo

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

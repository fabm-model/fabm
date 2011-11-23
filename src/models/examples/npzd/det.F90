#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_npzd_det --- AED detritus biogeochemical model
!
! !INTERFACE:
   module fabm_examples_npzd_det
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
   public type_examples_npzd_det, &
          examples_npzd_det_init, &
          examples_npzd_det_do,   &
          aed_d_do_ppdd, &
          aed_d_get_light_extinction, aed_d_get_conserved_quantities
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
! !PUBLIC DERIVED TYPES:
   type type_examples_npzd_det
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_d
      _TYPE_STATE_VARIABLE_ID_      :: id_zoo, id_phy, id_mintarget
      _TYPE_DEPENDENCY_ID_          :: id_temp

!     Model parameters
      REALTYPE :: kc,rdn
      REALTYPE :: zoo_per_n
      logical  :: use_zoo
      logical  :: use_phy
      logical  :: do_min
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the Detritus model
!
! !INTERFACE:
   subroutine examples_npzd_det_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_det), intent(out)    :: self
   type (type_model_info),intent(inout)          :: modelinfo
   integer,               intent(in)             :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                  :: d_initial=4.5
   REALTYPE                  :: w_d=-5.787037e-05
   REALTYPE                  :: rdn=3.472222e-08
   REALTYPE                  :: kc=0.03
   character(len=64)         :: mineralisation_target_variable=''

   REALTYPE, parameter :: secs_pr_day = 86400.
   namelist /examples_npzd_det/ &
            d_initial, w_d, rdn, kc, mineralisation_target_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=examples_npzd_det,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%rdn  = rdn /secs_pr_day

   ! Register state variables
   self%id_d = register_state_variable(modelinfo,'det','mmol/m**3','detritus', &
                                    d_initial,minimum=_ZERO_,vertical_movement=w_d/secs_pr_day)

   ! Register external state variable dependencies
   self%do_min = mineralisation_target_variable .ne. ''
   if (self%do_min) self%id_mintarget = register_state_dependency(modelinfo,mineralisation_target_variable)

   ! Register diagnostic variables

   ! Register conserved quantities

   ! Register environmental dependencies
   self%id_temp = register_dependency(modelinfo, varname_temp)

   return

99 call fatal_error('examples_npzd_det_init','Error reading namelist npzd')

   end subroutine examples_npzd_det_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Detritus model
!
! !INTERFACE:
   subroutine examples_npzd_det_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_det), intent(in)     :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: d,temp
   REALTYPE                   :: rpd
   REALTYPE, parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_d,d) ! detritus

   ! Retrieve current environmental conditions.
   _GET_DEPENDENCY_(self%id_temp,temp)  ! temperature

   ! Set temporal derivatives
   _SET_ODE_(self%id_d,-self%rdn*d)


   ! If an externally maintained NUT pool is present, add mineralisation to it
   if (self%do_min) then
     _SET_ODE_(self%id_mintarget, self%rdn*d)
   end if

   ! Export diagnostic variables

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine examples_npzd_det_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   pure subroutine aed_d_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_det), intent(in)     :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: d
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_d,d) ! detritus

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(self%kc*d)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_d_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   pure subroutine aed_d_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_det), intent(in)     :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: d
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_d,d) ! detritus

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_d_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Detritus model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine aed_d_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_det), intent(in)     :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: d,temp
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_d,d) ! detritus

   ! Retrieve current environmental conditions.
   _GET_DEPENDENCY_(self%id_temp,temp)

   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_ [as opposed to _SET_DD_], assignments to dd(i,j)
   ! are automatically assigned to pp(j,i) as well.
  ! _SET_DD_SYM_(self%id_d,self%id_d,self%rdn*d)             ! sdn  ?????????????????????
   _SET_DD_(self%id_d,self%id_d,self%rdn*d)             ! sdn  ?????????????????????

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   if (self%do_min) _SET_PP_(self%id_mintarget,self%id_mintarget,self%rdn*d)

   ! Export diagnostic variables

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_d_do_ppdd
!EOC

!-----------------------------------------------------------------------

   end module fabm_examples_npzd_det

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

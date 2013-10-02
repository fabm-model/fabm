#ifdef _FABM_F2003_

#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_npzd_det - Fennel & Neumann 1996 NPZD model - detritus component
!
! !INTERFACE:
   module fabm_examples_npzd_det
!
! !DESCRIPTION:
! This model features a single detritus variable, characterized by a rate of decay (rdn)
! and a sinking rate. Mineralized detritus feeds into a dissolved mineral pool that must
! be provided by an external model (e.g., fabm_examples_npzd_nut).
!
! !USES:
   use fabm_types
   use fabm_driver
   
   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_examples_npzd_det
!     Variable identifiers
      type (type_state_variable_id)     :: id_d
      type (type_state_variable_id)     :: id_mintarget
      type (type_conserved_quantity_id) :: id_totN

!     Model parameters
      real(rk) :: rdn
      logical  :: do_min

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
! !IROUTINE: Initialise the Detritus model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples_npzd_det namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_det), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk)          :: d_initial
   real(rk)          :: w_d
   real(rk)          :: rdn
   real(rk)          :: kc
   character(len=64) :: mineralisation_target_variable

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   namelist /examples_npzd_det/ d_initial, w_d, rdn, kc, mineralisation_target_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   d_initial = 4.5_rk
   w_d       = -5.0_rk
   kc        = 0.03_rk
   rdn       = 0.003_rk
   mineralisation_target_variable = ''

   ! Read the namelist
   read(configunit,nml=examples_npzd_det,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%rdn = rdn/secs_pr_day

   ! Register state variables
   call self%register_state_variable(self%id_d,'det','mmol/m**3','detritus',           &
                                d_initial,minimum=0.0_rk,vertical_movement=w_d/secs_pr_day, &
                                specific_light_extinction=kc)

   ! Register external state variable dependencies
   self%do_min = mineralisation_target_variable/=''
   if (self%do_min) call self%register_state_dependency(self%id_mintarget,mineralisation_target_variable)

   call self%register_conserved_quantity(self%id_totN,standard_variables%total_nitrogen)
   call self%add_conserved_quantity_component(self%id_totN,self%id_d)

   return

99 call fatal_error('examples_npzd_det_init','Error reading namelist examples_npzd_det')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Detritus model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_det), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: d
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_d,d) ! detritus

   ! Set temporal derivatives
   _SET_ODE_(self%id_d,-self%rdn*d)

   ! If an externally maintained NUT pool is present, add mineralisation to it
   if (self%do_min) then
     _SET_ODE_(self%id_mintarget, self%rdn*d)
   end if

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of Detritus model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_det), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
! !LOCAL VARIABLES:
   real(rk)                   :: d
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_d,d) ! detritus

   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_ [as opposed to _SET_DD_], assignments to dd(i,j)
   ! are automatically assigned to pp(j,i) as well.
  ! _SET_DD_SYM_(self%id_d,self%id_d,self%rdn*d)             ! sdn  ?????????????????????
   _SET_DD_(self%id_d,self%id_d,self%rdn*d)             ! sdn  ?????????????????????

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   if (self%do_min) _SET_PP_(self%id_mintarget,self%id_mintarget,self%rdn*d)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd
!EOC

!-----------------------------------------------------------------------

   end module fabm_examples_npzd_det

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#endif
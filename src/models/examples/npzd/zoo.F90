#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_npzd_zoo - Fennel & Neumann 1996 NPZD model - zooplankton component
!
! !INTERFACE:
   module fabm_examples_npzd_zoo
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none

   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_examples_npzd_zoo
!     Variable identifiers
      type (type_state_variable_id)     :: id_z
      type (type_state_variable_id)     :: id_exctarget,id_morttarget,id_grztarget

!     Model parameters
      real(rk) :: z0,gmax,iv,rzn,rzd

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
! !IROUTINE: Initialise the NPZD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples_npzd_zoo namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_zoo), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!
! !LOCAL VARIABLES:
   real(rk)                  :: z_initial
   real(rk)                  :: z0
   real(rk)                  :: gmax
   real(rk)                  :: iv
   real(rk)                  :: rzn
   real(rk)                  :: rzd
   character(len=64)         :: excretion_target_variable
   character(len=64)         :: mortality_target_variable
   character(len=64)         :: grazing_target_variable

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   namelist /examples_npzd_zoo/ &
            z_initial,z0,gmax,iv,rzn,rzd,excretion_target_variable, &
            mortality_target_variable,grazing_target_variable
!EOP
!-----------------------------------------------------------------------
!BOC
   z_initial = 0.0_rk
   z0        = 0.0225_rk
   gmax      = 0.5_rk
   iv        = 1.1_rk
   rzn       = 0.01_rk
   rzd       = 0.02_rk
   excretion_target_variable = ''
   mortality_target_variable = ''
   grazing_target_variable   = ''

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=examples_npzd_zoo,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%z0,  'z0',  default=z0)
   call self%get_parameter(self%gmax,'gmax',default=gmax, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%iv,  'iv',  default=iv)
   call self%get_parameter(self%rzn, 'rzn', default=rzn, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%rzd, 'rzd', default=rzd, scale_factor=1.0_rk/secs_pr_day)

   ! Register state variables
   call self%register_state_variable(self%id_z,'zoo','mmol/m**3','zooplankton',z_initial,minimum=0.0_rk)

   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_z)

   ! Register dependencies on external state variables.
   call self%register_state_dependency(self%id_exctarget, 'excretion_target','mmol/m**3','sink for excreted matter')
   call self%register_state_dependency(self%id_morttarget,'mortality_target','mmol/m**3','sink for dead matter')
   call self%register_state_dependency(self%id_grztarget, 'grazing_target',  'mmol/m**3','prey source')

   ! Automatically couple dependencies if target variables have been specified.
   if (excretion_target_variable/='') call self%request_coupling(self%id_exctarget, excretion_target_variable)
   if (mortality_target_variable/='') call self%request_coupling(self%id_morttarget,mortality_target_variable)
   if (grazing_target_variable  /='') call self%request_coupling(self%id_grztarget, grazing_target_variable)

   return

99 call self%fatal_error('examples_npzd_zoo::initialize','Error reading namelist examples_npzd_zoo')

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
   class (type_examples_npzd_zoo), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: p,z
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_z,z)         ! zooplankton
   _GET_(self%id_grztarget,p) ! phytoplankton

   ! Set temporal derivatives
    _SET_ODE_(self%id_z,fpz(self,p,z) - self%rzn*z - self%rzd*z)

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)
   _SET_ODE_(self%id_grztarget,-fpz(self,p,z))
   _SET_ODE_(self%id_morttarget,self%rzd*z)
   _SET_ODE_(self%id_exctarget,self%rzn*z)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
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
   class (type_examples_npzd_zoo), intent(in)     :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
! !LOCAL VARIABLES:
   real(rk)                   :: p,z
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_z,z)         ! zooplankton
   _GET_(self%id_grztarget,p) ! phytoplankton

   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_ [as opposed to _SET_DD_], assignments to dd(i,j)
   ! are automatically assigned to pp(j,i) as well.
    _SET_DD_SYM_(self%id_grztarget,self%id_z,fpz(self,p,z))
    _SET_DD_SYM_(self%id_z,self%id_exctarget,self%rzn*z)
    _SET_DD_SYM_(self%id_z,self%id_morttarget,self%rzd*z)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure real(rk) function fpz(self,p,z)
!
! !DESCRIPTION:
! Here, the classical Ivlev formulation for zooplankton grazing on
! phytoplankton is formulated.
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_zoo), intent(in) :: self
   real(rk),                       intent(in) :: p,z
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   fpz = self%gmax*(1.0_rk-exp(-self%iv*self%iv*p*p))*(z+self%z0)

   end function fpz
!EOC
!-----------------------------------------------------------------------

   end module fabm_examples_npzd_zoo

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: examples_npzd_zoo - Fennel & Neumann 1996 NPZD model - zooplankton component
!
! !INTERFACE:
   module examples_npzd_zoo
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
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day and are converted here to values per second.
   call self%get_parameter(self%z0,  'z0',  'mmol m-3', 'background concentration',     default=0.0225_rk)
   call self%get_parameter(self%gmax,'gmax','d-1',      'maximum specific grazing rate',default=0.5_rk,  scale_factor=d_per_s)
   call self%get_parameter(self%iv,  'iv',  'm3 mmol-1','Ivlev grazing constant',       default=1.1_rk)
   call self%get_parameter(self%rzn, 'rzn', 'd-1',      'excretion rate',               default=0.01_rk, scale_factor=d_per_s)
   call self%get_parameter(self%rzd, 'rzd', 'd-1',      'mortality',                    default=0.02_rk, scale_factor=d_per_s)

   ! Register state variables
   call self%register_state_variable(self%id_z,'c','mmol m-3','concentration',0.0_rk,minimum=0.0_rk)

   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_z)

   ! Register dependencies on external state variables.
   call self%register_state_dependency(self%id_grztarget, 'grazing_target',  'mmol m-3','prey source')
   call self%register_state_dependency(self%id_exctarget, 'excretion_target','mmol m-3','sink for excreted matter')
   call self%register_state_dependency(self%id_morttarget,'mortality_target','mmol m-3','sink for dead matter')

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

   end module examples_npzd_zoo

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

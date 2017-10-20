#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: examples_npzd_nut - Fennel & Neumann 1996 NPZD model - nutrient component
!
! !INTERFACE:
   module examples_npzd_nut
!
! !DESCRIPTION:
! This is a general nutrient (passive non-sinking, non-floating tracer), characterized by
! an initial concentration only.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_examples_npzd_nut
!     Variable identifiers
      type (type_state_variable_id) :: id_n
   contains
      procedure :: initialize
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the nutrient component
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the examples_npzd_nut namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_examples_npzd_nut), intent(inout), target :: self
   integer,                        intent(in)            :: configunit
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Register state variables
   call self%register_state_variable(self%id_n,'c','mmol m-3','concentration',     &
                                1.0_rk,minimum=0.0_rk,no_river_dilution=.true.)

   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------

   end module examples_npzd_nut

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

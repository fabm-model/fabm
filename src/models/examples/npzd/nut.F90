#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_npzd_nut - Fennel & Neumann 1996 NPZD model - nutrient component
!
! !INTERFACE:
   module fabm_examples_npzd_nut
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
      type (type_state_variable_id)     :: id_n

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
!
! !LOCAL VARIABLES:
   real(rk)                  :: n_initial
   namelist /examples_npzd_nut/ n_initial
!EOP
!-----------------------------------------------------------------------
!BOC
   n_initial = 4.5_rk

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=examples_npzd_nut,err=99)

   ! Register state variables
   call self%register_state_variable(self%id_n,'nut','mmol/m**3','nutrients',     &
                                n_initial,minimum=0.0_rk,no_river_dilution=.true.)

   ! Register contribution of state to global aggregate variables.
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_n)

   return

99 call self%fatal_error('examples_npzd_nut::initialize','Error reading namelist examples_npzd_nut')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------

   end module fabm_examples_npzd_nut

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

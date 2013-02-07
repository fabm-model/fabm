#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_examples_npzd_nut --- nutrients biogeochemical model
!
! !INTERFACE:
   module fabm_examples_npzd_nut
!
! !DESCRIPTION:
! The NPZD (nutrient-phytoplankton-zooplankton-detritus) model described here
! consists of $I=4$ state variables.
! Nutrient uptake (phytoplankton growth) is limited by light and nutrient
! availability, the latter of which is modelled by means
! of Michaelis-Menten kinetics, see eq.\ (\ref{dnp}).
! The half-saturation nutrient concentration $\alpha$ used in this
! formulation has typically a value between 0.2 and 1.5 mmol N\, m$^{-3}$.
! Zooplankton grazing which is limited by the phytoplankton standing stock
! is modelled by means of an Ivlev formulation, see eq.\ (\ref{dpz}).
! All other processes are based on linear first-order kinematics,
! see eqs.\ (\ref{dpn}) - (\ref{dzd}).
! For all details of the NPZD model implemented here,
! see \cite{Burchardetal2005b}.
!
! !USES:
   use fabm_types
   use fabm_driver
   
   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_examples_npzd_nut, examples_npzd_nut_init
!
! !PUBLIC DERIVED TYPES:
   type type_examples_npzd_nut
!     Variable identifiers
      type (type_state_variable_id)      :: id_n
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the nutrient componet
!
! !INTERFACE:
   subroutine examples_npzd_nut_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the examples_npzd_nut namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_nut), intent(out)   :: self
   _CLASS_ (type_model_info),     intent(inout) :: modelinfo
   integer,                       intent(in)    :: namlst
!
! !LOCAL VARIABLES:
   real(rk)                  :: n_initial
   namelist /examples_npzd_nut/ n_initial
!EOP
!-----------------------------------------------------------------------
!BOC
   n_initial = 4.5_rk

   ! Read the namelist
   read(namlst,nml=examples_npzd_nut,err=99)

   ! Register state variables
   call register_state_variable(modelinfo,self%id_n,'nut','mmol/m**3','nutrients',     &
                                n_initial,minimum=0.0_rk,no_river_dilution=.true.)

   return

99 call fatal_error('examples_npzd_nut_init','Error reading namelist examples_npzd_nut')

   end subroutine examples_npzd_nut_init
!EOC

!-----------------------------------------------------------------------

   end module fabm_examples_npzd_nut

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

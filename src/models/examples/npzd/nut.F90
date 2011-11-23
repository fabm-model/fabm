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

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_examples_npzd_nut, &
          examples_npzd_nut_init, &
          examples_npzd_nut_do,   &
          aed_n_do_ppdd, &
          aed_n_get_light_extinction, aed_n_get_conserved_quantities
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!
! !PUBLIC DERIVED TYPES:
   type type_examples_npzd_nut
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_n

!     Model parameters
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
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_nut), intent(out)    :: self
   type (type_model_info),intent(inout)          :: modelinfo
   integer,               intent(in)             :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                  :: n_initial=4.5
   REALTYPE, parameter       :: secs_pr_day=86400.
   namelist /examples_npzd_nut/ n_initial
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=examples_npzd_nut,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   ! Register state variables
   self%id_n = register_state_variable(modelinfo,'nut','mmol/m**3','nutrients',     &
                                    n_initial,minimum=_ZERO_,no_river_dilution=.true.)

   ! Register link to external pools

   ! Register diagnostic variables

   ! Register conserved quantities

   ! Register environmental dependencies

   return

99 call fatal_error('aed_n_init','Error reading namelist npzd')

   end subroutine examples_npzd_nut_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine examples_npzd_nut_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_nut), intent(in)     :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: n
   REALTYPE                   :: dn
   REALTYPE, parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_n,n) ! nutrient

   ! Set temporal derivatives
   dn = _ZERO_
   _SET_ODE_(self%id_n,dn)

   ! If an externally maintained pool is present, change the  pool according

   ! Export diagnostic variables

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine examples_npzd_nut_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   pure subroutine aed_n_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_nut), intent(in)     :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   _FABM_LOOP_END_
   end subroutine aed_n_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of conserved quantities (currently only nitrogen)
!
! !INTERFACE:
   pure subroutine aed_n_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_nut), intent(in)     :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE                     :: n
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_n,n) ! nutrient

   ! Total nutrient is simply the sum of all variables.
!   _SET_CONSERVED_QUANTITY_(self%id_totN,n+p+z+d)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_n_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model exporting production/destruction matrices
!
! !INTERFACE:
   subroutine aed_n_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_examples_npzd_nut), intent(in)     :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE                   :: n,p,z,d,par,I_0
   REALTYPE                   :: iopt,rpd,dn,primprod
   REALTYPE, parameter        :: secs_pr_day = 86400.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_n,n) ! nutrient

   ! Retrieve current environmental conditions.


   ! Export diagnostic variables

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine aed_n_do_ppdd
!EOC

!-----------------------------------------------------------------------

   end module fabm_examples_npzd_nut

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

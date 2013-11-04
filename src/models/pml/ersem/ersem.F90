#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_pml_ersem --- PML's ERSEM biogeochemical model ,
! adapted for FABM by momm@pml.ac.uk
!
! !INTERFACE:
   module fabm_pml_ersem
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types
   use fabm_driver
#ifdef FABM_PMLERSEM
   use global_declarations
   use ersem, ONLY:allocate_ersem,init_ersem,ersem_loop
   use allocationHelpers, ONLY: allocerr
   use ncdfRestartErsem, ONLY: readErsemRestart
#endif

   implicit none

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_pml_ersem, pml_ersem_init, pml_ersem_set_domain, &
          pml_ersem_do, pml_ersem_get_light_extinction
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!

!  Original author(s): Momme Butenschön
!
!
! !PUBLIC DERIVED TYPES:
   type type_pml_ersem
!     Variable identifiers
      type (type_state_variable_id),allocatable      :: id_ccc(:),id_ccb(:)
      type (type_dependency_id)          :: id_EIR,id_ETW,id_x1X,id_EPW
!      type (type_dependency_id)          :: id_Water,id_SeaSurface,id_BoxDepth
!      type (type_dependency_id)          :: id_SeaFloor,id_BoxFaces,id_Bathymetry

!     Model parameters
!      real(rk) :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
!      real(rk) :: dic_per_n
   end type
   logical  :: bioshade_feedback
   integer :: nbudget
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the PML-ERSEM model
!
! !INTERFACE:
   subroutine pml_ersem_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_pml_ersem),    intent(out)   :: self
   _CLASS_ (type_model_info),intent(inout) :: modelinfo
   integer,                  intent(in)    :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschön
!
! !LOCAL VARIABLES:
   integer :: n,ialloc
   character(len=255) :: ncdfErsemFile,ncdfErsemTitle

!   real(rk), parameter :: secs_pr_day = 86400.
#ifdef FABM_PMLERSEM
   namelist /pml_ersem/ ncdfErsemFile,ncdfErsemTitle, &
                        bioshade_feedback, &
                        nbudget,readErsemRestart, &
                        ncdfInstOut,ncdfDailyOut,ncdfWeeklyOut,ncdfMonthlyOut

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   !read(namlst,nml=pml_ersem,err=99,end=100)

   ! CURRENTLY BROKEN!
   ! ERSEM allocates arrays that describe the model (cccstr, ccbstr) and arrays that
   ! have a spatial dimension simultaneously (call to allocate_ersem).
   ! This is no longer allowed in FABM - non-spatial model information must be set during
   ! initialization, and spatial information is provided in a call to set_domain.
   ! Currently I_STATE, I_STATEBEN, cccstr and ccbstr will not be available,
   ! causing the code below to break.
   !call allocate_ersem()

   allocate(self%id_ccc(I_STATE),stat=ialloc)
   call allocerr ('id_ccc',1,i_state,ialloc)
   allocate(self%id_ccb(I_STATEBEN),stat=ialloc)
   call allocerr ('id_ccb',1,i_stateben,ialloc)
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   ! wdepth = WILL NEED TO BE ASSIGNED FOR ANALYTIC BEN INIT TO WORK

   call init_ersem()
! Register state variables
   do n=1,I_STATE
        self%id_ccc(n) = register_state_variable(modelinfo,cccstr(n),'undefined',cccstr(n),     &
                                    ccc(1,n)) ! EVENTUALLY SET A MINIMUM AT LATER STAGE, replace 2nd cccstr(n) with real long name
   enddo
   do n=1,I_STATEBEN
        self%id_ccb(n) = register_state_variable(modelinfo,ccbstr(n),'undefined',ccbstr(n),     &
                                    ccb(0,n),benthic=.true.) ! EVENTUALLY SET A MINIMUM AT LATER STAGE, replace 2nd ccbstr(n) with real long name
   enddo
#endif

   ! Register diagnostic variables
   ! NONE FOR NOW
   ! Register conserved quantities
!   self%id_totN = register_conserved_quantity(modelinfo,'N','mmol/m**3','nitrogen') LATER

   ! Register environmental dependencies
   call register_dependency(modelinfo, self%id_ETW, standard_variables%temperature)
   call register_dependency(modelinfo, self%id_EIR, standard_variables%downwelling_shortwave_flux)
   call register_dependency(modelinfo, self%id_EPW, standard_variables%pressure)
   call register_dependency(modelinfo, self%id_x1X, standard_variables%practical_salinity)

   return

99 call fatal_error('pml_ersem_init','Error reading namelist pml_ersem')
100 call fatal_error('pml_ersem_init','Namelist pml_ersem was not found')

   end subroutine pml_ersem_init
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Specify the spatial domain
!
! !INTERFACE:
   subroutine pml_ersem_set_domain(self,domainsize)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_pml_ersem), intent(in)    :: self
   integer,               intent(in)    :: domainsize
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschön
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef FABM_PMLERSEM
   N_COMP = domainsize
   call allocate_ersem()
#endif

   end subroutine pml_ersem_set_domain
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate sink and source terms (rates) of the pelagic component of ERSEM.
!
! !INTERFACE:
   subroutine pml_ersem_do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   type (type_pml_ersem),       intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschön
!
! !LOCAL VARIABLES:
   real(rk), parameter        :: secs_pr_day = 86400.
   integer :: n
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_1D_

#ifdef FABM_PMLERSEM
   ! Retrieve current (local) state variable values.
   DO n=1,I_STATE
        _GET_STATE_1D_(self%id_ccc(n),ccc(_DOMAIN_1D_,n)) ! pelagic
   ENDDO

   ! Retrieve current environmental conditions.
   _GET_DEPENDENCY_1D_(self%id_EIR,EIR(_DOMAIN_1D_))  ! local short wave radiation
   _GET_DEPENDENCY_1D_(self%id_ETW,ETW(_DOMAIN_1D_))  ! local temperature
   _GET_DEPENDENCY_1D_(self%id_x1X,x1X(_DOMAIN_1D_))  ! local salinity
!   _GET_DEPENDENCY_1D_(self%id_EPW,EPW(_DOMAIN_1D_))  ! local pressure

   ! TODO:
   ! (1) Benthos should be disabled in ersem_loop; it must be called separately instead.
   ! (2) Inclusion of subsidence in sink and source terms may be disabled (no call to calc_subsidence)
   !     See discussion in pml_ersem_get_vertical_movement below.
   call ersem_loop()

   ! Set temporal derivatives
   do n=1,I_STATE
     _SET_ODE_1D_(self%id_ccc(n),sccc(_DOMAIN_1D_,n)/secs_pr_day)
   enddo

   ! Export diagnostic variables
!   _SET_DIAGNOSTIC_1D_(self%id_dPAR,par)
   call reset_ersem()

#endif

   ! Leave spatial loops (if any)
   _LOOP_END_1D_

   end subroutine pml_ersem_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate pelagic bottom fluxes and benthic sink and source terms of ERSEM.
! Everything in units per surface area (not volume!) per time.
!
! !INTERFACE:
   subroutine pml_ersem_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   type (type_pml_ersem),       intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   real(rk), parameter        :: secs_pr_day = 86400.
   integer                    :: n
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_1D_

#ifdef FABM_PMLERSEM
   ! Retrieve current (local) state variable values for the bottom pelagic layer.
   DO n=1,I_STATE
        _GET_STATE_HZ_1D_(self%id_ccc(n),ccc(_INDEX_HZ_1D_,n))
   ENDDO

   ! Retrieve current (local) state variable values for the benthos.
   DO n=1,I_STATEBEN
        _GET_STATE_BEN_1D_(self%id_ccb(n),ccb(_DOMAIN_HZ_1D_,n))
   ENDDO

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_DEPENDENCY_HZ_1D_(self%id_EIR,EIR(_INDEX_HZ_1D_))  ! local short wave radiation
   _GET_DEPENDENCY_HZ_1D_(self%id_ETW,ETW(_INDEX_HZ_1D_))  ! local temperature
   _GET_DEPENDENCY_HZ_1D_(self%id_x1X,x1X(_INDEX_HZ_1D_))  ! local salinity
!   _GET_DEPENDENCY_HZ_1D_(self%id_EPW,EPW(_INDEX_HZ_1D_))  ! local pressure

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   do n=1,I_STATE
     _SET_BOTTOM_FLUX_1D_(self%id_ccc(n),0.0)
   enddo

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   do n=1,I_STATEBEN
     _SET_ODE_BEN_1D_(self%id_ccb(n),sccb(_DOMAIN_HZ_1D_,n))
   enddo
#endif

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_1D_

   end subroutine pml_ersem_do_benthos
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine pml_ersem_get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   type (type_pml_ersem), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                     :: n
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_1D_

#ifdef FABM_PMLERSEM
   ! Retrieve current (local) state variable values.
   DO n=1,I_STATE
        _GET_STATE_1D_(self%id_ccc(n),ccc(_DOMAIN_1D_,n)) ! pelagic
   ENDDO
   call calculate_extinction()
   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_1D_(xEPS(_DOMAIN_1D_))
#endif

   ! Leave spatial loops (if any)
   _LOOP_END_1D_

   end subroutine pml_ersem_get_light_extinction
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical velocity of pelagic biogeochemical variables
! [not needed as long as ERSEM handles subsidence as part of its sink and source terms]
!
! !INTERFACE:
   subroutine pml_ersem_get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
!
! !INPUT PARAMETERS:
   type (type_pml_ersem), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                    :: n
   real(rk), parameter        :: secs_pr_day = 86400.
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_1D_

   ! Here a call would be needed to ensure sdCCC contains the correct velocities
   ! for the current point in horizontal space. This may be non-trivial, though,
   ! as ERSEM seems to calculate these velocities in the same code that handles sink and source
   ! terms.
   !
   ! Solutions:
   ! (1) let ERSEM handle subsidence as part of its sinks and source terms
   !     (all vertical velocities would be set to zero for FABM). This is the fastest solution.
   ! (2) Cache the full 3D field for sdCCC while calcualting sink and source terms
   !     per column [pml_ersem_do], and re-use it here. ERSEM does not include subsidence
   !     in sink ands oruce terms. Expensive in terms of memory.
   ! (3) Rework the ERSEM code to isolate the calculation of subsidence rates from the
   !     calculation of sink and soruce terms. Bets in the long term, but requires
   !     substantial changes.

#ifdef FABM_PMLERSEM
   do n=1,I_STATE
      _SET_VERTICAL_MOVEMENT_1D_(self%id_ccc(n),-sdCCC(_DOMAIN_1D_,n)/secs_pr_day)
   end do
#endif

   ! Leave spatial loops (if any)
   _LOOP_END_1D_

   end subroutine pml_ersem_get_vertical_movement
!EOC
!-----------------------------------------------------------------------

   end module fabm_pml_ersem

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

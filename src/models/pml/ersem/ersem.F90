!$Id: pmlersem.F90 119 2010-12-27 14:23:18Z jornbr $
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_pmlersem --- PML's ERSEM biogeochemical model ,
! adapted for FABM by momm@pml.ac.uk
!
! !INTERFACE:
   module fabm_pmlersem
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
   public type_pmlersem, pmlersem_init, pmlersem_do, &
          pmlersem_get_light_extinction
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!

!  Original author(s): Momme Butenschön
!
!
! !PUBLIC DERIVED TYPES:
   type type_pmlersem
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_,allocatable      :: id_ccc(:),id_ccb(:)
      _TYPE_DEPENDENCY_ID_          :: id_EIR,id_ETW,id_x1X,id_EPW
!      _TYPE_DEPENDENCY_ID_          :: id_Water,id_SeaSurface,id_BoxDepth
!      _TYPE_DEPENDENCY_ID_          :: id_SeaFloor,id_BoxFaces,id_Bathymetry
      
!     Model parameters
!      REALTYPE :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
!      REALTYPE :: dic_per_n
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
   subroutine pmlersem_init(self,modelinfo,namlst,domainsize)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_pmlersem),  intent(out)   :: self
   type (type_model_info),intent(inout) :: modelinfo
   integer,               intent(in)    :: namlst,domainsize
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschön
!
! !LOCAL VARIABLES:
   integer :: n,ialloc
   character(len=255) :: ncdfErsemFile,ncdfErsemTitle

!   REALTYPE, parameter :: secs_pr_day = 86400.
#ifdef FABM_PMLERSEM
   namelist /fabmersem_nml/ ncdfErsemFile,ncdfErsemTitle, &
                    bioshade_feedback, &
                    nbudget,readErsemRestart, &
                    ncdfInstOut,ncdfDailyOut,ncdfWeeklyOut,ncdfMonthlyOut

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   !read(namlst,nml=fabmersem,err=99)

   N_COMP = domainsize
   call allocate_ersem()

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
   self%id_ETW = register_dependency(modelinfo, varname_temp)
   self%id_EIR = register_dependency(modelinfo, varname_swr)
   self%id_EPW = register_dependency(modelinfo, varname_pres)
   self%id_x1X = register_dependency(modelinfo, varname_salt)

   return

99 call fatal_error('pmlersem_init','Error reading namelist fabmersem')
   
   end subroutine pmlersem_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate sink and source terms (rates) of the pelagic component of ERSEM.
!
! !INTERFACE:
   subroutine pmlersem_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   type (type_pmlersem),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Momme Butenschön
!
! !LOCAL VARIABLES:
   REALTYPE, parameter        :: secs_pr_day = 86400.
   integer :: n
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_1D_

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
   !     See discussion in pmlersem_get_vertical_movement below.
   call ersem_loop()

   ! Set temporal derivatives
   do n=1,I_STATE
     _SET_ODE_1D_(self%id_ccc(n),sccc(_DOMAIN_1D_,n)/secs_pr_day)
   enddo

   ! Export diagnostic variables
!   _SET_DIAG_1D_(self%id_dPAR,par)
   call reset_ersem()

#endif
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_1D_

   end subroutine pmlersem_do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate pelagic bottom fluxes and benthic sink and source terms of ERSEM.
! Everything in units per surface area (not volume!) per time.
!
! !INTERFACE:
   subroutine pmlersem_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   type (type_pmlersem),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   REALTYPE, parameter        :: secs_pr_day = 86400.
   integer                    :: n
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_HZ_LOOP_BEGIN_1D_

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
   _FABM_HZ_LOOP_END_1D_

   end subroutine pmlersem_do_benthos
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine pmlersem_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   type (type_pmlersem), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
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
   _FABM_LOOP_BEGIN_1D_

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
   _FABM_LOOP_END_1D_
   
   end subroutine pmlersem_get_light_extinction
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical velocity of pelagic biogeochemical variables
! [not needed as long as ERSEM handles subsidence as part of its sink and source terms]
!
! !INTERFACE:
   subroutine pmlersem_get_vertical_movement(self,_FABM_ARGS_GET_VERTICAL_MOVEMENT_)
!
! !INPUT PARAMETERS:
   type (type_pmlersem), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_VERTICAL_MOVEMENT_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                    :: n
   REALTYPE, parameter        :: secs_pr_day = 86400.
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_1D_
   
   ! Here a call would be needed to ensure sdCCC contains the correct velocities
   ! for the current point in horizontal space. This may be non-trivial, though,
   ! as ERSEM seems to calculate these velocities in the same code that handles sink and source
   ! terms.
   !
   ! Solutions:
   ! (1) let ERSEM handle subsidence as part of its sinks and source terms
   !     (all vertical velocities would be set to zero for FABM). This is the fastest solution.
   ! (2) Cache the full 3D field for sdCCC while calcualting sink and source terms
   !     per column [pmlersem_do], and re-use it here. ERSEM does not include subsidence
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
   _FABM_LOOP_END_1D_
   
   end subroutine pmlersem_get_vertical_movement
!EOC
!-----------------------------------------------------------------------

   end module fabm_pmlersem

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

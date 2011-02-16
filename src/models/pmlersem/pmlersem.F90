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
   use global_declarations
   use ersem, ONLY:allocate_ersem,init_ersem,ersem_loop
   use allocationHelpers, ONLY: allocerr
   use ncdfRestartErsem, ONLY: readErsemRestart

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
   subroutine pmlersem_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_pmlersem),      intent(out)   :: self
   type (type_model_info),intent(inout) :: modelinfo
   integer,               intent(in)    :: namlst 
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
! !LOCAL VARIABLES:
   integer :: n,ialloc
   character(len=255) :: ncdfErsemFile,ncdfErsemTitle

!   REALTYPE, parameter :: secs_pr_day = 86400.

   namelist /fabmersem_nml/ ncdfErsemFile,ncdfErsemTitle, &
                    bioshade_feedback, &
                    nbudget,readErsemRestart, &
                    ncdfInstOut,ncdfDailyOut,ncdfWeeklyOut,ncdfMonthlyOut

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   !read(namlst,nml=fabmersem,err=99)

   N_COMP=_GET_NLEV_
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
        self%id_ccc(n) = register_state_variable(modelinfo,cccstr(n),'undefined','undefined',     &
                                    ccc(1,n)) ! EVENTUALLY SET A MINIMUM AT LATER STAGE
   enddo
   do n=1,I_STATEBEN
        self%id_ccb(n) = register_state_variable(modelinfo,ccbstr(n),'undefined','undefined',     &
                                    ccb(1,n)) ! EVENTUALLY SET A MINIMUM AT LATER STAGE
   enddo


   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
!   self%use_dic = dic_variable.ne.''
!   if (self%use_dic) self%id_dic = register_state_dependency(modelinfo,dic_variable)

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
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine pmlersem_do(self,_FABM_ARGS_DO_RHS_)
!
! !DESCRIPTION:
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_pmlersem),       intent(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
! !LOCAL VARIABLES:
   REALTYPE, parameter        :: secs_pr_day = 86400.
   integer :: n
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_1D_

   ! Retrieve current (local) state variable values.
   DO n=1,I_STATE
        _GET_STATE_1D_(self%id_ccc(n),ccc(_DOMAIN_1D_,n)) ! pelagic
   ENDDO
!   DO n=1,I_STATEBEN
!        _GET_STATE_HZ(self%id_ccb(:,n),ccb(:,n)) ! benthic
!   ENDDO
   
   ! Retrieve current environmental conditions.
   _GET_DEPENDENCY_1D_   (self%id_EIR,EIR)  ! local short wave radiation
   _GET_DEPENDENCY_1D_   (self%id_ETW,ETW)  ! local temperature
   _GET_DEPENDENCY_1D_   (self%id_x1X,x1X)  ! local salinity
!   _GET_DEPENDENCY_1D_   (self%id_EPW,EPW)  ! local pressure

   call ersem_loop()

   ! Set temporal derivatives
   _SET_ODE_1D_(self%id_ccc,sccc)
   !_SET_ODE_1D_(self%id_ccb,sccb)

   ! If an externally maintained DIC pool is present, change the DIC pool according to the
   ! the change in nutrients (assuming constant C:N ratio)

   ! Export diagnostic variables
!   _SET_DIAG_(self%id_dPAR,par)
!   _SET_DIAG_(self%id_GPP ,primprod)
!   _SET_DIAG_(self%id_NCP ,primprod - self%rpn*p)
!   _SET_DIAG_(self%id_PPR ,primprod*secs_pr_day)
!   _SET_DIAG_(self%id_NPR ,(primprod - self%rpn*p)*secs_pr_day)
   
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_1D_

   end subroutine pmlersem_do
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

   ! Retrieve current (local) state variable values.
   DO n=1,I_STATE
        _GET_STATE_1D_(self%id_ccc(n),ccc(_DOMAIN_1D_,n)) ! pelagic
   ENDDO
   call calculate_extinction() 
   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_1D_(xEPS(_DOMAIN_1D_))

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_1D_
   
   end subroutine pmlersem_get_light_extinction
!EOC
!-----------------------------------------------------------------------

   end module fabm_pmlersem

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

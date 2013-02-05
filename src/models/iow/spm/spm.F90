#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_iow_spm --- 1-class SPM model,
!
! !INTERFACE:
   module fabm_iow_spm
!
! !USES:
   use fabm_types
   use fabm_driver

!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_iow_spm, iow_spm_init, iow_spm_do_benthos, iow_spm_get_light_extinction
!
! !PRIVATE DATA MEMBERS:
!
! !REVISION HISTORY:!
!  Original author(s): Manuel Ruiz Villarreal & Richard Hofmeister
!
!
! !PUBLIC DERIVED TYPES:
   type type_iow_spm
!     Variable identifiers
      type (type_state_variable_id)        :: id_spm !concentrations
      type (type_bottom_state_variable_id) :: id_pmpool !sediment pool
      type (type_horizontal_dependency_id) :: id_taub    !bottom stress

!     Model parameters
      REALTYPE :: ws
      REALTYPE :: c_init
      REALTYPE :: mass_sed_init
      REALTYPE :: tauc_erosion
      REALTYPE :: tauc_sedimentation
      REALTYPE :: erosion_const
      REALTYPE :: shading

   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise the sediment model
!
! !INTERFACE:
   subroutine iow_spm_init(self,modelinfo,namlst)
!
! !DESCRIPTION:
!  Here, the spm namelist is read and te variables exported
!  by the model are registered with FABM.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_iow_spm),      intent(out)   :: self
   type (type_model_info),intent(inout) :: modelinfo
   integer,               intent(in)    :: namlst
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
! !LOCAL VARIABLES:

   REALTYPE, parameter :: secs_pr_day = 86400.d0   ! s
   REALTYPE            :: c_init=5.0d0             ! mg/l
   REALTYPE            :: mass_sed_init=100.d0     ! g/m**2
   REALTYPE            :: tauc_erosion=2.d-2       ! N/m**2
   REALTYPE            :: tauc_sedimentation=2.d-2 ! N/m**2
   REALTYPE            :: erosion_const=1.d-2
   REALTYPE            :: ws=-10.0d0               ! m/d
   REALTYPE            :: shading=1.0d0            ! 1/m per mg/l

   namelist /iow_spm/ c_init, shading, &
                         mass_sed_init, &
                         tauc_erosion, &
                         tauc_sedimentation, &
                         erosion_const, &
                         ws
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read the namelist
   read(namlst,nml=iow_spm,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%tauc_erosion = tauc_erosion
   self%tauc_sedimentation = tauc_sedimentation
   self%erosion_const = erosion_const
   self%shading = shading
   self%ws = ws/secs_pr_day

   ! Register state variables
   call register_state_variable(modelinfo,self%id_spm,'spm','mg/l','concentration of SPM',     &
                                    c_init,minimum=_ZERO_, &
                                    vertical_movement=self%ws, &
                                    no_river_dilution=.true.)

   call register_state_variable(modelinfo,self%id_pmpool,'pmpool','g/sqm','mass/sqm of PM in sediment',  &
                                    mass_sed_init)

   ! Register environmental dependencies
   call register_dependency(modelinfo, self%id_taub, varname_taub)

   return

99 call fatal_error('spm_init','Error reading namelist spm')

   end subroutine iow_spm_init
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Sedimentation/Erosion
!
! !INTERFACE:
   subroutine iow_spm_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!
! !DESCRIPTION:
! Calculating the benthic fluxes
!
   implicit none

! !INPUT PARAMETERS:
   type (type_iow_spm), intent(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_

! !LOCAL VARIABLES:
   REALTYPE                     :: taub,spm,pmpool
   REALTYPE                     :: porosity
   REALTYPE, parameter          :: rho_0=1025.d0 ! [g/l]
   REALTYPE                     :: Erosion_Flux,Sedimentation_Flux
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
!EOP
!-----------------------------------------------------------------------
!BOC
   porosity=0.d0

   _FABM_HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_spm,spm)
   _GET_HORIZONTAL_(self%id_pmpool,pmpool)
   _GET_HORIZONTAL_(self%id_taub,taub)

   ! 1-spm_porosity is the fractional bed concentration
   if(pmpool .gt. _ZERO_) then
   ! if there are sediments in the pool
      Erosion_Flux = self%erosion_const / rho_0                  &
                    *(1.d0-porosity) * (taub-self%tauc_erosion )
      Erosion_Flux = max(Erosion_Flux,_ZERO_)
   else
      Erosion_Flux = _ZERO_
   end if


   !sedimentation flux:
   Sedimentation_Flux = min(_ZERO_,self%ws * spm * (1.d0-taub / self%tauc_sedimentation))

   ! unit is g/m**3 * m/s
   _SET_BOTTOM_EXCHANGE_(self%id_spm,Sedimentation_Flux+Erosion_Flux)
   ! unit is g/m**2/s
   _SET_ODE_BEN_(self%id_pmpool,-Erosion_Flux-Sedimentation_Flux)

   _FABM_HORIZONTAL_LOOP_END_

   end subroutine iow_spm_do_benthos
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: iow_spm_get_light_extinction
!
! !INTERFACE:
   subroutine iow_spm_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)

   implicit none

! !INPUT PARAMETERS:
   type (type_iow_spm), intent(in) :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
! !LOCAL VARIABLES
   REALTYPE      :: spm
!
! !REVISION HISTORY:
!  Original author(s): Richard Hofmeister
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   _GET_(self%id_spm,spm)
   _SET_EXTINCTION_(self%shading*spm)

   _FABM_LOOP_END_

   end subroutine iow_spm_get_light_extinction
!EOC
!-----------------------------------------------------------------------

   end module fabm_iow_spm


!###############################################################################
!#                                                                             #
!# aed_viruses.F90                                                              #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created March 2012                                                          #
!#                                                                             #
!###############################################################################

#ifdef _FABM_F2003_

#include "aed.h"

!
MODULE aed_viruses
!-------------------------------------------------------------------------------
! aed_viruses --- viruses biogeochemical model
!
! The AED module viruses contains equations that describe exchange of
! soluable reactive viruses across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE fabm_types
   USE fabm_driver

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_viruses, aed_viruses_create
!
   TYPE,extends(type_base_model) :: type_aed_viruses
!     Variable identifiers
      type (type_state_variable_id)      :: id_vir
      type (type_dependency_id)          :: id_temp
      type (type_diagnostic_variable_id) :: id_sed_vir
      type (type_conserved_quantity_id)  :: id_totV

!     Model parameters
      real(rk) :: num_viruses
!     LOGICAL  :: use_oxy,use_vir

      CONTAINS      ! Model Methods
!       procedure :: initialize               => aed_viruses_init
        procedure :: do                       => aed_viruses_do
        procedure :: do_ppdd                  => aed_viruses_do_ppdd
        procedure :: do_benthos               => aed_viruses_do_benthos
        procedure :: get_conserved_quantities => aed_viruses_get_conserved_quantities
   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
FUNCTION aed_viruses_create(namlst,name,parent) RESULT(self)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and the variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                       :: namlst
   CHARACTER(len=*),INTENT(in)              :: name
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_aed_viruses),POINTER :: self

   real(rk)          :: num_viruses

   real(rk),PARAMETER :: secs_pr_day = 86400.
   real(rk),PARAMETER :: vir_initial = 0.
   NAMELIST /aed_viruses/ num_viruses
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL initialize_model_info(self,name,parent)

   print *,"WARNING! aed_viruses model is currently under development"

   ! Read the namelist
   read(namlst,nml=aed_viruses,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%num_viruses = num_viruses

   ! Register state variables
   call self%register_state_variable(self%id_vir,'vir','mmol/m**3','viruses',     &
                                    vir_initial,minimum=0.0_rk,no_river_dilution=.false.)

   ! Register external state variable dependencies

   PRINT *,'AED_VIRUSES : Note this module has not been completed. Stopping.'


   ! Register conserved quantities
   call self%register_conserved_quantity(self%id_totV,'TV','mmol/m**3','Total viruses')

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   RETURN

99 CALL fatal_error('aed_viruses_init','Error reading namelist aed_viruses')

END FUNCTION aed_viruses_create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_viruses_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_viruses model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_viruses),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   real(rk)           :: vir,diff_vir
   real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_vir,vir) ! viruses

   ! Set temporal derivatives
   diff_vir = 0.

   _SET_ODE_(self%id_vir,diff_vir)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_viruses_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_viruses_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of viruses biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_viruses),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   real(rk)                   :: vir
   real(rk)                   :: diff_vir
   real(rk), parameter        :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_vir,vir) ! viruses

   ! Set temporal derivatives
   diff_vir = 0.

   _SET_PP_(self%id_vir,self%id_vir,diff_vir)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_viruses_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_viruses_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED viruses.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_viruses),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   ! Environment
   real(rk) :: temp

   ! State
   real(rk) :: vir

   ! Temporary variables
!  real(rk) :: vir_flux

   ! Parameters
   real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_(self%id_vir,vir) ! viruses



   ! Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_

END SUBROUTINE aed_viruses_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_viruses_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities (currently only viruses)
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_viruses),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   real(rk) :: vir
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
!  _GET_(self%id_vir,vir) ! viruses

   ! Total nutrient is simply the sum of all variables.
!  _SET_CONSERVED_QUANTITY_(self%id_totV,vir)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_viruses_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_viruses
#endif

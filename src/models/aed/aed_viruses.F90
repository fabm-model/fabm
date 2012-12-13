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

#include "fabm_driver.h"

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
      _TYPE_STATE_VARIABLE_ID_      :: id_vir
      _TYPE_DEPENDENCY_ID_          :: id_temp
      _TYPE_DIAGNOSTIC_VARIABLE_ID_ :: id_sed_vir
      _TYPE_CONSERVED_QUANTITY_ID_  :: id_totV

!     Model parameters
      REALTYPE :: num_viruses
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

   REALTYPE          :: num_viruses

   REALTYPE,PARAMETER :: secs_pr_day = 86400.
   REALTYPE,PARAMETER :: vir_initial = 0.
   NAMELIST /aed_viruses/ num_viruses
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL self%initialize(name,parent)

   print *,"WARNING! aed_viruses model is currently under development"

   ! Read the namelist
   read(namlst,nml=aed_viruses,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%num_viruses = num_viruses

   ! Register state variables
   self%id_vir = self%register_state_variable('vir','mmol/m**3','viruses',     &
                                    vir_initial,minimum=_ZERO_,no_river_dilution=.false.)

   ! Register external state variable dependencies

   PRINT *,'AED_VIRUSES : Note this module has not been completed. Stopping.'


   ! Register conserved quantities
   self%id_totV = self%register_conserved_quantity('TV','mmol/m**3','Total viruses')

   ! Register environmental dependencies
   self%id_temp = self%register_dependency(varname_temp)

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
   REALTYPE           :: vir,diff_vir
   REALTYPE,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_vir,vir) ! viruses

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
   REALTYPE                   :: vir
   REALTYPE                   :: diff_vir
   REALTYPE, parameter        :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_vir,vir) ! viruses

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
   REALTYPE :: temp

   ! State
   REALTYPE :: vir

   ! Temporary variables
!  REALTYPE :: vir_flux

   ! Parameters
   REALTYPE,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_HZ_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_DEPENDENCY_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_vir,vir) ! viruses



   ! Leave spatial loops (if any)
   _FABM_HZ_LOOP_END_

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
   REALTYPE :: vir
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
!  _GET_STATE_(self%id_vir,vir) ! viruses

   ! Total nutrient is simply the sum of all variables.
!  _SET_CONSERVED_QUANTITY_(self%id_totV,vir)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_viruses_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_viruses
#endif

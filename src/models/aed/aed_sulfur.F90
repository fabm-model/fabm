!###############################################################################
!#                                                                             #
!# aed_sulfur.F90                                                              #
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
MODULE aed_sulfur
!-------------------------------------------------------------------------------
! aed_sulfur --- sulfur biogeochemical model
!
! The AED module sulfur contains equations that describe exchange of
! soluable reactive sulfur across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE fabm_types
   USE fabm_driver

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_sulfur, aed_sulfur_create
!
   TYPE,extends(type_base_model) :: type_aed_sulfur
!     Variable identifiers
      type (type_state_variable_id)      :: id_so4
      type (type_dependency_id)          :: id_temp
      type (type_diagnostic_variable_id) :: id_sed_so4
      type (type_conserved_quantity_id)  :: id_totS

!     Model parameters
      real(rk) :: Fsed_dic,Ksed_dic,theta_sed_dic
      LOGICAL  :: use_oxy,use_dic

      CONTAINS     ! Model Methods
!       procedure :: initialize               => aed_sulfur_init
        procedure :: do                       => aed_sulfur_do
        procedure :: do_ppdd                  => aed_sulfur_do_ppdd
        procedure :: do_benthos               => aed_sulfur_do_benthos
        procedure :: get_conserved_quantities => aed_sulfur_get_conserved_quantities
   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
FUNCTION aed_sulfur_create(namlst,name,parent) RESULT(self)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: namlst
   CHARACTER(len=*),INTENT(in)              :: name
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_aed_sulfur),POINTER :: self

   INTEGER           :: num_sulfurs
   real(rk)          :: decay(100)
   real(rk)          :: settling(100)
   real(rk)          :: Fsed(100)

   real(rk),PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_sulfur/ num_sulfurs,decay,settling,Fsed
!
!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL initialize_model_info(self,name,parent)

   ! Read the namelist
   read(namlst,nml=aed_sulfur,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   PRINT *,'AED_SULFUR : Note this module has not been completed. Stopping.'


!  ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   RETURN

99 CALL fatal_error('aed_sulfur_init','Error reading namelist aed_sulfur')

END FUNCTION aed_sulfur_create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sulfur_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_sulfur model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_sulfur),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
!  real(rk)           :: dic,diff_dic
!  real(rk),PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
!  _GET_(self%id_dic,dic) ! sulfur

   ! Set temporal derivatives
!  diff_dic = 0.

!  _SET_ODE_(self%id_dic,diff_dic)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_sulfur_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sulfur_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of sulfur biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_sulfur),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   real(rk)                   :: dic
   real(rk)                   :: diff_dic
   real(rk), parameter        :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
!  _GET_(self%id_dic,dic) ! sulfur

   ! Set temporal derivatives
!  diff_dic = 0.

!  _SET_PP_(self%id_dic,self%id_dic,diff_dic)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_sulfur_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sulfur_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED sulfur.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_sulfur),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   ! Environment
   real(rk) :: temp

   ! State
!  real(rk) :: dic,oxy

   ! Temporary variables
!  real(rk) :: dic_flux

   ! Parameters
   real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
!  _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
!  _GET_(self%id_dic,dic) ! sulfur

!  IF (self%use_oxy) THEN
!     ! Sediment flux dependent on oxygen and temperature
!     _GET_(self%id_oxy,oxy)
!     dic_flux = self%Fsed_dic * self%Ksed_dic/(self%Ksed_dic+oxy) * (self%theta_sed_dic**(temp-20.0))
!  ELSE
!     ! Sediment flux dependent on temperature only.
!     dic_flux = self%Fsed_dic * (self%theta_sed_dic**(temp-20.0))
!  ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to FABM.
   !_SET_BOTTOM_FLUX_(self%id_dic,dic_flux/secs_pr_day)
   !_SET_ODE_SED_FLUX_(self%id_dic,dic_flux)
!  _SET_BOTTOM_EXCHANGE_(self%id_dic,dic_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_dic,-dic_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
!  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_dic,dic_flux)

   ! Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_

END SUBROUTINE aed_sulfur_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sulfur_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities (currently only sulfur)
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_sulfur),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   real(rk) :: dic
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
!  _GET_(self%id_dic,dic) ! sulfur

   ! Total nutrient is simply the sum of all variables.
!  _SET_CONSERVED_QUANTITY_(self%id_totC,dic)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_sulfur_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_sulfur
#endif

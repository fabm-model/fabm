!###############################################################################
!#                                                                             #
!# aed_sedflux.F90                                                             #
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
!# Created May 2012                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"

#define _MAX_ZONES_ 100
!
MODULE aed_sedflux
!------------------------------------------------------------------------------+
! The AED module sediment is not truely a model in itself, rather it provides  |
! sediment flux values in a unified way to simply the interface to other models|
!------------------------------------------------------------------------------+
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_type_sedflux
!
   TYPE,extends(type_base_model) :: aed_type_sedflux
!     Variable identifiers
      TYPE (type_horizontal_diagnostic_variable_id) :: id_Fsed_oxy, id_Fsed_rsi, id_Fsed_amm, id_Fsed_nit
      TYPE (type_horizontal_diagnostic_variable_id) :: id_Fsed_frp, id_Fsed_pon, id_Fsed_don
      TYPE (type_horizontal_diagnostic_variable_id) :: id_Fsed_pop, id_Fsed_dop, id_Fsed_poc, id_Fsed_doc
      TYPE (type_horizontal_diagnostic_variable_id) :: id_Fsed_dic, id_Fsed_ch4, id_Fsed_feii

!     Model parameters
      INTEGER  :: sed_modl, n_zones
      AED_REAL :: Fsed_oxy, Fsed_rsi, Fsed_amm, Fsed_nit, Fsed_frp, &
                  Fsed_pon, Fsed_don, Fsed_pop, Fsed_dop, &
                  Fsed_poc, Fsed_doc, Fsed_dic, Fsed_ch4, Fsed_feii
      AED_REAL,DIMENSION(:),ALLOCATABLE :: &
                  Fsed_oxy_P, Fsed_rsi_P, Fsed_amm_P, Fsed_nit_P, Fsed_frp_P, &
                  Fsed_pon_P, Fsed_don_P, Fsed_pop_P, Fsed_dop_P, &
                  Fsed_poc_P, Fsed_doc_P, Fsed_dic_P, Fsed_ch4_P, Fsed_feii_P

      CONTAINS     ! Model Procedures
        PROCEDURE :: initialize  => aed_init_sedflux
        PROCEDURE :: do_benthos  => aed_sedflux_do_benthos
   END TYPE

   AED_REAL,PARAMETER :: secs_pr_day = 86400.


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_init_sedflux(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_sedflux),TARGET,INTENT(inout) :: self
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status

!  AED_REAL          :: FsedA_initial=0.01
!  AED_REAL          :: FsedN_initial=0.01
   CHARACTER(len=64) :: sedflux_model=''

   INTEGER  :: nzones = 1
   AED_REAL :: Fsed_oxy = MISVAL, Fsed_rsi = MISVAL, Fsed_amm = MISVAL, Fsed_nit = MISVAL, &
               Fsed_pon = MISVAL, Fsed_don = MISVAL, Fsed_pop = MISVAL, Fsed_dop = MISVAL, &
               Fsed_poc = MISVAL, Fsed_doc = MISVAL, Fsed_dic = MISVAL, Fsed_frp = MISVAL, &
               Fsed_ch4 = MISVAL, Fsed_feii = MISVAL

   NAMELIST /aed_sedflux/ sedflux_model
   NAMELIST /aed_sed_constant/ nzones,                                           &
                               Fsed_oxy, Fsed_rsi, Fsed_amm, Fsed_nit, Fsed_frp, &
                               Fsed_pon, Fsed_don, Fsed_pop, Fsed_dop,           &
                               Fsed_poc, Fsed_doc, Fsed_dic, Fsed_ch4, Fsed_feii

!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"aed_sedflux initialization"

   ! Read the namelist
   read(namlst,nml=aed_sedflux,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_sedflux'

   self%sed_modl = -1
   IF ( sedflux_model .EQ. "Constant" ) THEN
      read(namlst,nml=aed_sed_constant,iostat=status)
      IF (status /= 0) STOP 'Error reading namelist aed_sed_constant'

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day,
      ! and are converted here to values per second.
      self%Fsed_oxy = Fsed_oxy/secs_pr_day
      self%Fsed_rsi = Fsed_rsi/secs_pr_day
      self%Fsed_amm = Fsed_amm/secs_pr_day
      self%Fsed_nit = Fsed_nit/secs_pr_day
      self%Fsed_frp = Fsed_frp/secs_pr_day
      self%Fsed_pon = Fsed_pon/secs_pr_day
      self%Fsed_don = Fsed_don/secs_pr_day
      self%Fsed_pon = Fsed_pop/secs_pr_day
      self%Fsed_don = Fsed_dop/secs_pr_day
      self%Fsed_poc = Fsed_poc/secs_pr_day
      self%Fsed_doc = Fsed_doc/secs_pr_day
      self%Fsed_dic = Fsed_dic/secs_pr_day
      self%Fsed_ch4 = Fsed_ch4/secs_pr_day
      self%Fsed_feii = Fsed_feii/secs_pr_day
      self%sed_modl = 1
   ENDIF

   ! Register state variables
   ! NOTE the benthic=.true. argument, which specifies the variable is benthic.
   IF ( Fsed_oxy .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_oxy,'Fsed_oxy','mmol/m**2',   &
                                          'sedimentation rate of oxygen', output=output_instantaneous)
   IF ( Fsed_rsi .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_rsi,'Fsed_rsi','mmol/m**2',   &
                                          'sedimentation rate of silica', output=output_instantaneous)
   IF ( Fsed_amm .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_amm,'Fsed_amm','mmol/m**2',   &
                                          'sedimentation rate of ammonia', output=output_instantaneous)
   IF ( Fsed_nit .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_nit,'Fsed_nit','mmol/m**2',   &
                                          'sedimentation rate of nitrate', output=output_instantaneous)
   IF ( Fsed_frp .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_frp,'Fsed_frp','mmol/m**2',   &
                                          'sedimentation rate of phosphorus', output=output_instantaneous)
   IF ( Fsed_pon .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_pon,'Fsed_pon','mmol/m**2',   &
                                          'sedimentation rate of pon', output=output_instantaneous)
   IF ( Fsed_don .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_don,'Fsed_don','mmol/m**2',   &
                                          'sedimentation rate of don', output=output_instantaneous)
   IF ( Fsed_pop .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_pop,'Fsed_pop','mmol/m**2',   &
                                          'sedimentation rate of pop', output=output_instantaneous)
   IF ( Fsed_dop .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_dop,'Fsed_dop','mmol/m**2',   &
                                          'sedimentation rate of dop', output=output_instantaneous)
   IF ( Fsed_poc .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_poc,'Fsed_poc','mmol/m**2',   &
                                          'sedimentation rate of poc', output=output_instantaneous)
   IF ( Fsed_doc .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_doc,'Fsed_doc','mmol/m**2',   &
                                          'sedimentation rate of doc', output=output_instantaneous)
   IF ( Fsed_dic .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_dic,'Fsed_dic','mmol/m**2',   &
                                          'sedimentation rate of carbon', output=output_instantaneous)
   IF ( Fsed_ch4 .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_ch4,'Fsed_ch4','mmol/m**2',   &
                                          'sedimentation rate of ch4', output=output_instantaneous)
   IF ( Fsed_feii .GT. MISVAL ) &
      CALL self%register_horizontal_diagnostic_variable(self%id_Fsed_feii,'Fsed_feii','mmol/m**2', &
                                          'sedimentation rate of iron', output=output_instantaneous)
END SUBROUTINE aed_init_sedflux
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sedflux_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED sediment.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_sedflux),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   AED_REAL :: Rzone
   INTEGER  :: zone
   ! Temporary variables
   AED_REAL :: Fsed_oxy, Fsed_rsi
   AED_REAL :: Fsed_amm, Fsed_nit
   AED_REAL :: Fsed_pon, Fsed_don
   AED_REAL :: Fsed_pop, Fsed_dop
   AED_REAL :: Fsed_poc, Fsed_doc
   AED_REAL :: Fsed_dic, Fsed_frp
   AED_REAL :: Fsed_ch4, Fsed_feii
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Constant model has nothing to do
   IF ( self%sed_modl .EQ. 0 ) RETURN

   !# Do this here because for the constant model these values never change.
   IF ( self%sed_modl .EQ. 1 ) THEN
      Fsed_oxy = self%Fsed_oxy
      Fsed_rsi = self%Fsed_rsi
      Fsed_amm = self%Fsed_amm
      Fsed_nit = self%Fsed_nit
      Fsed_frp = self%Fsed_frp
      Fsed_pon = self%Fsed_pon
      Fsed_don = self%Fsed_don
      Fsed_pop = self%Fsed_pop
      Fsed_dop = self%Fsed_dop
      Fsed_poc = self%Fsed_poc
      Fsed_doc = self%Fsed_doc
      Fsed_dic = self%Fsed_dic
      Fsed_ch4 = self%Fsed_ch4
      Fsed_feii = self%Fsed_feii
!     self%sed_modl = 0 ! From now on, we don't need to do this
!     Bother! can't do this because self is intent in
   ENDIF

   !# Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   !# Also store sediment flux as diagnostic variable.
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_oxy)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_oxy, Fsed_oxy)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_rsi)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_rsi, Fsed_rsi)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_amm)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_amm, Fsed_amm)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_nit)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_nit, Fsed_nit)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_frp)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_frp, Fsed_frp)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_pon)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_pon, Fsed_pon)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_don)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_don, Fsed_don)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_pop)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_pop, Fsed_pop)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_dop)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_dop, Fsed_dop)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_poc)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_poc, Fsed_poc)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_doc)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_doc, Fsed_doc)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_dic)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_dic, Fsed_dic)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_ch4)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_ch4, Fsed_ch4)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_feii)) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Fsed_feii, Fsed_feii)

   !# Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_
END SUBROUTINE aed_sedflux_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_sedflux

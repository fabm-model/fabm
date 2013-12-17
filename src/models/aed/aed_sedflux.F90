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

#ifdef _FABM_F2003_

#include "aed.h"

#ifdef SINGLE
#  define _OXY_MIN_ -30000.0
#  define _FRP_MIN_ -30000.0
#  define _RSI_MIN_ -30000.0
#  define _AMM_MIN_ -30000.0
#  define _NIT_MIN_ -30000.0
#  define _PON_MIN_ -30000.0
#  define _DON_MIN_ -30000.0
#  define _POP_MIN_ -30000.0
#  define _DOP_MIN_ -30000.0
#  define _POC_MIN_ -30000.0
#  define _DOC_MIN_ -30000.0
#  define _DIC_MIN_ -30000.0
#  define _CH4_MIN_ -30000.0
#  define _FEIIMIN_ -30000.0
#else
#  define _OXY_MIN_ -30000.0d0
#  define _FRP_MIN_ -30000.0d0
#  define _RSI_MIN_ -30000.0d0
#  define _AMM_MIN_ -30000.0d0
#  define _NIT_MIN_ -30000.0d0
#  define _PON_MIN_ -30000.0d0
#  define _DON_MIN_ -30000.0d0
#  define _POP_MIN_ -30000.0d0
#  define _DOP_MIN_ -30000.0d0
#  define _POC_MIN_ -30000.0d0
#  define _DOC_MIN_ -30000.0d0
#  define _DIC_MIN_ -30000.0d0
#  define _CH4_MIN_ -30000.0d0
#  define _FEIIMIN_ -30000.0d0
#endif



#define _MAX_ZONES_ 100
!
MODULE aed_sedflux
!------------------------------------------------------------------------------+
! The AED module sediment is not truely a model in itself, rather it provides  |
! sediment flux values in a unified way to simply the interface to other models|
!------------------------------------------------------------------------------+
   USE fabm_types

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_sedflux
!
   TYPE,extends(type_base_model) :: type_aed_sedflux
!     Variable identifiers
      type (type_horizontal_dependency_id)     :: id_zone
      type (type_bottom_state_variable_id)     :: id_Fsed_oxy, id_Fsed_rsi, id_Fsed_amm, id_Fsed_nit
      type (type_bottom_state_variable_id)     :: id_Fsed_frp, id_Fsed_pon, id_Fsed_don
      type (type_bottom_state_variable_id)     :: id_Fsed_pop, id_Fsed_dop, id_Fsed_poc, id_Fsed_doc
      type (type_bottom_state_variable_id)     :: id_Fsed_dic, id_Fsed_ch4, id_Fsed_feii
      type (type_horizontal_dependency_id)     :: id_zones
!     type (type_conserved_quantity_id) :: id_tot_sed

!     Model parameters
      INTEGER  :: sed_modl, n_zones
      real(rk) :: Fsed_oxy, Fsed_rsi, Fsed_amm, Fsed_nit, Fsed_frp, &
                  Fsed_pon, Fsed_don, Fsed_pop, Fsed_dop, &
                  Fsed_poc, Fsed_doc, Fsed_dic, Fsed_ch4, Fsed_feii
      real(rk),DIMENSION(:),ALLOCATABLE :: &
                  Fsed_oxy_P, Fsed_rsi_P, Fsed_amm_P, Fsed_nit_P, Fsed_frp_P, &
                  Fsed_pon_P, Fsed_don_P, Fsed_pop_P, Fsed_dop_P, &
                  Fsed_poc_P, Fsed_doc_P, Fsed_dic_P, Fsed_ch4_P, Fsed_feii_P

      CONTAINS     ! Model Methods
        procedure :: initialize               => aed_sedflux_init
        procedure :: do_benthos               => aed_sedflux_do_benthos
   END TYPE

   real(rk),PARAMETER :: secs_pr_day = 86400.


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE load_sed_zone_data(self,namlst)
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_sedflux),INTENT(inout) :: self
   INTEGER,INTENT(in)                       :: namlst
!
!LOCALS
   real(rk)          :: FsedA_initial=0.01
   real(rk)          :: FsedN_initial=0.01

   INTEGER  :: n_zones
   real(rk) :: Fsed_oxy(_MAX_ZONES_) = MISVAL, Fsed_rsi(_MAX_ZONES_) = MISVAL, &
               Fsed_amm(_MAX_ZONES_) = MISVAL, Fsed_nit(_MAX_ZONES_) = MISVAL, &
               Fsed_frp(_MAX_ZONES_) = MISVAL, Fsed_pon(_MAX_ZONES_) = MISVAL, &
               Fsed_don(_MAX_ZONES_) = MISVAL, Fsed_pop(_MAX_ZONES_) = MISVAL, &
               Fsed_dop(_MAX_ZONES_) = MISVAL, Fsed_poc(_MAX_ZONES_) = MISVAL, &
               Fsed_doc(_MAX_ZONES_) = MISVAL, Fsed_dic(_MAX_ZONES_) = MISVAL, &
               Fsed_ch4(_MAX_ZONES_) = MISVAL, Fsed_feii(_MAX_ZONES_) = MISVAL

   NAMELIST /aed_sed_const2d/ n_zones, &
                              Fsed_oxy, Fsed_rsi, Fsed_amm, Fsed_nit, Fsed_frp, &
                              Fsed_pon, Fsed_don, Fsed_pop, Fsed_dop,           &
                              Fsed_poc, Fsed_doc, Fsed_dic, Fsed_ch4, Fsed_feii
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_sed_const2d,err=99)

   self%n_zones = n_zones
   IF (Fsed_oxy(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_oxy_P(n_zones)) ; self%Fsed_oxy_P(1:n_zones) = Fsed_oxy(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_rsi(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_rsi_P(n_zones)) ; self%Fsed_rsi_P(1:n_zones) = Fsed_rsi(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_amm(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_amm_P(n_zones)) ; self%Fsed_amm_P(1:n_zones) = Fsed_amm(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_nit(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_nit_P(n_zones)) ; self%Fsed_nit_P(1:n_zones) = Fsed_nit(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_frp(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_frp_P(n_zones)) ; self%Fsed_frp_P(1:n_zones) = Fsed_frp(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_pon(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_pon_P(n_zones)) ; self%Fsed_pon_P(1:n_zones) = Fsed_pon(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_don(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_don_P(n_zones)) ; self%Fsed_don_P(1:n_zones) = Fsed_don(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_pop(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_pop_P(n_zones)) ; self%Fsed_pop_P(1:n_zones) = Fsed_pop(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_dop(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_dop_P(n_zones)) ; self%Fsed_dop_P(1:n_zones) = Fsed_dop(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_poc(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_poc_P(n_zones)) ; self%Fsed_poc_P(1:n_zones) = Fsed_poc(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_doc(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_doc_P(n_zones)) ; self%Fsed_doc_P(1:n_zones) = Fsed_doc(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_dic(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_dic_P(n_zones)) ; self%Fsed_dic_P(1:n_zones) = Fsed_dic(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_ch4(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_ch4_P(n_zones)) ; self%Fsed_ch4_P(1:n_zones) = Fsed_ch4(1:n_zones)/secs_pr_day
   ENDIF
   IF (Fsed_feii(1) .NE. MISVAL ) THEN
      ALLOCATE(self%Fsed_feii_P(n_zones)) ; self%Fsed_feii_P(1:n_zones) = Fsed_feii(1:n_zones)/secs_pr_day
   ENDIF
   RETURN

99 CALL self%fatal_error('aed_sedflux_init','Error reading namelist aed_sed_const2d')
END SUBROUTINE load_sed_zone_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sedflux_init(self,configunit)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_aed_sedflux), TARGET, INTENT(INOUT) :: self
   INTEGER,INTENT(in)                              :: configunit
!
!LOCALS

   real(rk)          :: FsedA_initial=0.01
   real(rk)          :: FsedN_initial=0.01
   CHARACTER(len=64) :: sedflux_model=''

   INTEGER  :: nzones = 1
   real(rk) :: Fsed_oxy = MISVAL, Fsed_rsi = MISVAL, Fsed_amm = MISVAL, Fsed_nit = MISVAL, &
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
   ! Read the namelist
   read(configunit,nml=aed_sedflux,err=98)

   self%sed_modl = -1
   IF ( sedflux_model .EQ. "Constant" ) THEN
      read(configunit,nml=aed_sed_constant,err=99)
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
   ELSEIF ( sedflux_model .EQ. "Spatially Variable" ) THEN
      call self%register_dependency(self%id_zones,'env_sed_zone')
      self%sed_modl = 2
      CALL load_sed_zone_data(self,configunit)
      IF (ALLOCATED(self%Fsed_oxy_P)) Fsed_oxy = self%Fsed_oxy_P(1)
      IF (ALLOCATED(self%Fsed_rsi_P)) Fsed_rsi = self%Fsed_rsi_P(1)
      IF (ALLOCATED(self%Fsed_amm_P)) Fsed_amm = self%Fsed_amm_P(1)
      IF (ALLOCATED(self%Fsed_nit_P)) Fsed_nit = self%Fsed_nit_P(1)
      IF (ALLOCATED(self%Fsed_frp_P)) Fsed_frp = self%Fsed_frp_P(1)
      IF (ALLOCATED(self%Fsed_pon_P)) Fsed_pon = self%Fsed_pon_P(1)
      IF (ALLOCATED(self%Fsed_don_P)) Fsed_don = self%Fsed_don_P(1)
      IF (ALLOCATED(self%Fsed_pop_P)) Fsed_pop = self%Fsed_pop_P(1)
      IF (ALLOCATED(self%Fsed_dop_P)) Fsed_dop = self%Fsed_dop_P(1)
      IF (ALLOCATED(self%Fsed_poc_P)) Fsed_poc = self%Fsed_poc_P(1)
      IF (ALLOCATED(self%Fsed_doc_P)) Fsed_doc = self%Fsed_doc_P(1)
      IF (ALLOCATED(self%Fsed_dic_P)) Fsed_dic = self%Fsed_dic_P(1)
      IF (ALLOCATED(self%Fsed_ch4_P)) Fsed_dic = self%Fsed_ch4_P(1)
      IF (ALLOCATED(self%Fsed_feii_P)) Fsed_dic = self%Fsed_feii_P(1)
   ENDIF

   ! Register state variables
   ! NOTE the benthic=.true. argument, which specifies the variable is benthic.
   IF ( Fsed_oxy .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_oxy,'Fsed_oxy','mmol/m**2',      &
                                          'sedimentation rate of oxygen',          &
                                          Fsed_oxy,minimum=_OXY_MIN_)
   IF ( Fsed_rsi .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_rsi,'Fsed_rsi','mmol/m**2',      &
                                          'sedimentation rate of silica',          &
                                          Fsed_rsi,minimum=_RSI_MIN_)
   IF ( Fsed_amm .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_amm,'Fsed_amm','mmol/m**2',      &
                                          'sedimentation rate of ammonia',         &
                                          Fsed_amm,minimum=_AMM_MIN_)
   IF ( Fsed_nit .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_nit,'Fsed_nit','mmol/m**2',      &
                                          'sedimentation rate of nitrate',         &
                                          Fsed_nit,minimum=_NIT_MIN_)
   IF ( Fsed_frp .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_frp,'Fsed_frp','mmol/m**2',      &
                                          'sedimentation rate of phosphorus',      &
                                          Fsed_frp,minimum=_FRP_MIN_)
   IF ( Fsed_pon .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_pon,'Fsed_pon','mmol/m**2',      &
                                          'sedimentation rate of pon',             &
                                          Fsed_pon,minimum=_PON_MIN_)
   IF ( Fsed_don .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_don,'Fsed_don','mmol/m**2',      &
                                          'sedimentation rate of don',             &
                                          Fsed_don,minimum=_DON_MIN_)
   IF ( Fsed_pop .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_pop,'Fsed_pop','mmol/m**2',      &
                                          'sedimentation rate of pop',             &
                                          Fsed_pop,minimum=_POP_MIN_)
   IF ( Fsed_dop .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_dop,'Fsed_dop','mmol/m**2',      &
                                          'sedimentation rate of dop',             &
                                          Fsed_dop,minimum=_DOP_MIN_)
   IF ( Fsed_poc .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_poc,'Fsed_poc','mmol/m**2',      &
                                          'sedimentation rate of poc',             &
                                          Fsed_poc,minimum=_POC_MIN_)
   IF ( Fsed_doc .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_doc,'Fsed_doc','mmol/m**2',      &
                                          'sedimentation rate of doc',             &
                                          Fsed_doc,minimum=_DOC_MIN_)
   IF ( Fsed_dic .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_dic,'Fsed_dic','mmol/m**2',      &
                                          'sedimentation rate of carbon',          &
                                          Fsed_dic,minimum=_DIC_MIN_)

   IF ( Fsed_ch4 .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_ch4,'Fsed_ch4','mmol/m**2',      &
                                          'sedimentation rate of carbon',          &
                                          Fsed_ch4,minimum=_CH4_MIN_)

   IF ( Fsed_feii .GT. MISVAL ) &
      call self%register_state_variable(self%id_Fsed_feii,'Fsed_feii','mmol/m**2',      &
                                          'sedimentation rate of carbon',          &
                                          Fsed_feii,minimum=_FEIIMIN_)

   RETURN

98 CALL self%fatal_error('aed_sedflux_init','Error reading namelist aed_sedflux')
   STOP
99 CALL self%fatal_error('aed_sedflux_init','Error reading namelist aed_sed_constant')
   STOP

END SUBROUTINE aed_sedflux_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_sedflux_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED sediment.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_sedflux),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   real(rk) :: Rzone
   INTEGER  :: zone
   ! Temporary variables
   real(rk) :: Fsed_oxy, Fsed_rsi
   real(rk) :: Fsed_amm, Fsed_nit
   real(rk) :: Fsed_pon, Fsed_don
   real(rk) :: Fsed_pop, Fsed_dop
   real(rk) :: Fsed_poc, Fsed_doc
   real(rk) :: Fsed_dic, Fsed_frp
   real(rk) :: Fsed_ch4, Fsed_feii
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
   _FABM_HORIZONTAL_LOOP_BEGIN_

   IF ( self%sed_modl .EQ. 2) THEN
      !# Get the zone array dependency
      !# select the material zone for this cell
      !# set sediment values accordingly
      _GET_HORIZONTAL_(self%id_zones,Rzone)
      zone = Rzone
!     print *,'kk ',kk,' from zone ',zone,' Rzone ',Rzone

      IF (zone .LE. 0 .OR. zone .GT. self%n_zones ) zone = 1

      IF (_VARIABLE_REGISTERED_(self%id_Fsed_oxy)) Fsed_oxy = self%Fsed_oxy_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_rsi)) Fsed_rsi = self%Fsed_rsi_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_amm)) Fsed_amm = self%Fsed_amm_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_nit)) Fsed_nit = self%Fsed_nit_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_frp)) Fsed_frp = self%Fsed_frp_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_pon)) Fsed_pon = self%Fsed_pon_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_don)) Fsed_don = self%Fsed_don_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_pop)) Fsed_pop = self%Fsed_pop_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_dop)) Fsed_dop = self%Fsed_dop_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_poc)) Fsed_poc = self%Fsed_poc_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_doc)) Fsed_doc = self%Fsed_doc_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_dic)) Fsed_dic = self%Fsed_dic_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_ch4)) Fsed_ch4 = self%Fsed_ch4_P(zone)
      IF (_VARIABLE_REGISTERED_(self%id_Fsed_feii)) Fsed_feii = self%Fsed_feii_P(zone)
   ENDIF

   !# Also store sediment flux as diagnostic variable.
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_oxy)) _SET_STATE_BEN_(self%id_Fsed_oxy, Fsed_oxy)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_rsi)) _SET_STATE_BEN_(self%id_Fsed_rsi, Fsed_rsi)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_amm)) _SET_STATE_BEN_(self%id_Fsed_amm, Fsed_amm)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_nit)) _SET_STATE_BEN_(self%id_Fsed_nit, Fsed_nit)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_frp)) _SET_STATE_BEN_(self%id_Fsed_frp, Fsed_frp)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_pon)) _SET_STATE_BEN_(self%id_Fsed_pon, Fsed_pon)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_don)) _SET_STATE_BEN_(self%id_Fsed_don, Fsed_don)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_pop)) _SET_STATE_BEN_(self%id_Fsed_pop, Fsed_pop)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_dop)) _SET_STATE_BEN_(self%id_Fsed_dop, Fsed_dop)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_poc)) _SET_STATE_BEN_(self%id_Fsed_poc, Fsed_poc)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_doc)) _SET_STATE_BEN_(self%id_Fsed_doc, Fsed_doc)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_dic)) _SET_STATE_BEN_(self%id_Fsed_dic, Fsed_dic)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_ch4)) _SET_STATE_BEN_(self%id_Fsed_ch4, Fsed_ch4)
   IF (_VARIABLE_REGISTERED_(self%id_Fsed_feii)) _SET_STATE_BEN_(self%id_Fsed_feii, Fsed_feii)

   !# Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_
END SUBROUTINE aed_sedflux_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_sedflux
#endif

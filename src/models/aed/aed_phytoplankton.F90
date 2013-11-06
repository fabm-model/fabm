!###############################################################################
!#                                                                             #
!# aed_phytoplankton.F90                                                       #
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
!# Created August 2011                                                         #
!#                                                                             #
!###############################################################################

#ifdef _FABM_F2003_

#include "aed.h"

MODULE aed_phytoplankton
!-------------------------------------------------------------------------------
!  aed_phytoplankton --- phytoplankton biogeochemical model
!-------------------------------------------------------------------------------
   USE fabm_types
   USE fabm_driver
   USE aed_util,ONLY : find_free_lun, &
                       exp_integral,  &
                       aed_bio_temp_function, &
                       fTemp_function

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC type_aed_phytoplankton, aed_phytoplankton_create
!
   TYPE phyto_data
      ! General Attributes
      CHARACTER(64) :: p_name
      real(rk) :: p0, Xcc, kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      ! Growth rate parameters
      INTEGER  :: fT_Method
      real(rk) :: R_growth, theta_growth, T_std, T_opt, T_max, kTn, aTn, bTn
      ! Light configuration and parameters
      INTEGER  :: lightModel
      real(rk) :: I_K, I_S, KePHY
      ! Respiration parameters
      real(rk) :: f_pr, R_resp, k_fdom, k_fres, theta_resp
      ! Salinity parameters
      INTEGER  :: salTol
      real(rk) :: S_bep, S_maxsp, S_opt
      ! Nitrogen parameters
      INTEGER  :: simDINUptake, simDONUptake, simNFixation, simINDynamics
      real(rk) :: N_o, K_N, X_nmin, X_nmax, X_ncon, R_nuptake, k_nfix, R_nfix
      ! Phosphorus parameters
      INTEGER  :: simDIPUptake, simIPDynamics
      real(rk) :: P_0, K_P, X_pmin, X_pmax, X_pcon, R_puptake
      ! Silica parameters
      INTEGER  :: simSiUptake
      real(rk) :: Si_0, K_Si, X_sicon
      ! Carbon parameters
      INTEGER  :: simCUptake, dic_mode
      ! Sedimentation parameters
      real(rk) :: w_p
   END TYPE


   TYPE phyto_nml_data
      CHARACTER(64) :: p_name
      real(rk) :: p_initial
      real(rk) :: p0, w_p, Xcc, R_growth !i_min,rmax,alpha,rpn,rpdu,rpdl
      INTEGER  :: fT_Method
      real(rk) :: theta_growth, T_std, T_opt, T_max
      INTEGER  :: lightModel
      real(rk) :: I_K, I_S, KePHY
      ! Respiration parameters
      real(rk) :: f_pr, R_resp, theta_resp, k_fres, k_fdom
      ! Salinity parameters
      INTEGER  :: salTol
      real(rk) :: S_bep, S_maxsp, S_opt
      ! Nitrogen parameters
      INTEGER  :: simDINUptake, simDONUptake, simNFixation, simINDynamics
      real(rk) :: N_o, K_N, X_ncon, X_nmin, X_nmax, R_nuptake, k_nfix, R_nfix
      ! Phosphorus parameters
      INTEGER  :: simDIPUptake, simIPDynamics
      real(rk) :: P_0, K_P, X_pcon, X_pmin, X_pmax, R_puptake
      ! Silica parameters
      INTEGER  :: simSiUptake
      real(rk) :: Si_0, K_Si, X_sicon
    !  ! Carbon parameters
    !  INTEGER  :: simCUptake, dic_mode
   END TYPE


   TYPE,extends(type_base_model) :: type_aed_phytoplankton
!     Variable identifiers
      type (type_state_variable_id),ALLOCATABLE :: id_p(:)
      type (type_state_variable_id),ALLOCATABLE :: id_in(:)
      type (type_state_variable_id),ALLOCATABLE :: id_ip(:)
      type (type_state_variable_id)        :: id_Pexctarget,id_Pmorttarget,id_Pupttarget(1:2)
      type (type_state_variable_id)        :: id_Nexctarget,id_Nmorttarget,id_Nupttarget(1:4)
      type (type_state_variable_id)        :: id_Cexctarget,id_Cmorttarget,id_Cupttarget
      type (type_state_variable_id)        :: id_Siexctarget,id_Simorttarget,id_Siupttarget
      type (type_state_variable_id)        :: id_DOupttarget
      type (type_dependency_id)            :: id_par, id_tem, id_sal, id_dz, id_extc
      type (type_horizontal_dependency_id) :: id_I_0
      type (type_diagnostic_variable_id)   :: id_GPP, id_NCP, id_PPR, id_NPR, id_dPAR
      type (type_diagnostic_variable_id)   :: id_TPHY, id_TCHLA, id_TIN, id_TIP
      type (type_diagnostic_variable_id)   :: id_NUP, id_PUP, id_CUP
      type (type_diagnostic_variable_id),ALLOCATABLE :: id_NtoP(:)
      type (type_diagnostic_variable_id),ALLOCATABLE :: id_fT(:), id_fI(:), id_fNit(:), id_fPho(:), id_fSil(:), id_fSal(:)
      type (type_conserved_quantity_id)    :: id_totP

!     Model parameters
      INTEGER                                   :: num_phytos
      TYPE(phyto_data),DIMENSION(:),ALLOCATABLE :: phytos
      ! LOGICAL                                 :: do_exc,do_mort,do_upt, do_N2uptake
      LOGICAL                                   :: do_Puptake, do_Nuptake, do_Cuptake
      LOGICAL                                   :: do_Siuptake, do_DOuptake, do_N2uptake
      LOGICAL                                   :: do_Pmort, do_Nmort, do_Cmort, do_Simort
      LOGICAL                                   :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER :: nnup, npup
      real(rk) :: dic_per_n

      CONTAINS     ! Model Methods
!       procedure :: initialize               => aed_phytoplankton_init
        procedure :: do                       => aed_phytoplankton_do
        procedure :: do_ppdd                  => aed_phytoplankton_do_ppdd
        procedure :: do_benthos               => aed_phytoplankton_do_benthos
        procedure :: get_light_extinction     => aed_phytoplankton_get_light_extinction
        procedure :: get_conserved_quantities => aed_phytoplankton_get_conserved_quantities
   END TYPE

   INTEGER  :: ino3, inh4,idon, in2, ifrp, idop
   real(rk) :: dtlim = 0.9 * 3600
   LOGICAL  :: extra_debug = .false.

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_phytoplankton_load_params(self, count, list)
!-------------------------------------------------------------------------------
   _CLASS_ (type_aed_phytoplankton),INTENT(inout) :: self
   INTEGER,INTENT(in)                             :: count
   INTEGER,INTENT(in)                             :: list(*)

!  real(rk)    :: p_initial=0.
!  real(rk)    :: p0=0.0225
!  real(rk)    :: w_p=-1.157407e-05
!  real(rk)    :: i_min=25.
!  real(rk)    :: rmax=1.157407e-05
!  real(rk)    :: alpha=0.3
!  real(rk)    :: rpn=1.157407e-07
!  real(rk)    :: rpdu=2.314814e-07
!  real(rk)    :: rpdl=1.157407e-06

   INTEGER     :: i,tfil
   real(rk)    :: minNut

   real(rk), parameter :: secs_pr_day = 86400.
   TYPE(phyto_nml_data) :: pd(MAX_PHYTO_TYPES)
   NAMELIST /phyto_data/ pd
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file="aed_phyto_pars.nml", status='OLD')
    read(tfil,nml=phyto_data,err=99)
    close(tfil)

    self%num_phytos = count
    ALLOCATE(self%phytos(count))
    ALLOCATE(self%id_p(count))
    ALLOCATE(self%id_in(count))
    ALLOCATE(self%id_ip(count))
    ALLOCATE(self%id_NtoP(count))
    IF (extra_debug) THEN
       ALLOCATE(self%id_fT(count))
       ALLOCATE(self%id_fI(count))
       ALLOCATE(self%id_fNit(count))
       ALLOCATE(self%id_fPho(count))
       ALLOCATE(self%id_fSil(count))
       ALLOCATE(self%id_fSal(count))
    ENDIF

    DO i=1,count
       ! Assign parameters from database to simulated groups
       self%phytos(i)%p_name       = pd(list(i))%p_name
       self%phytos(i)%p0           = pd(list(i))%p0
       self%phytos(i)%w_p          = pd(list(i))%w_p/secs_pr_day
       self%phytos(i)%Xcc          = pd(list(i))%Xcc
       self%phytos(i)%R_growth     = pd(list(i))%R_growth/secs_pr_day
       self%phytos(i)%fT_Method    = pd(list(i))%fT_Method
       self%phytos(i)%theta_growth = pd(list(i))%theta_growth
       self%phytos(i)%T_std        = pd(list(i))%T_std
       self%phytos(i)%T_opt        = pd(list(i))%T_opt
       self%phytos(i)%T_max        = pd(list(i))%T_max
       self%phytos(i)%lightModel   = pd(list(i))%lightModel
       self%phytos(i)%I_K          = pd(list(i))%I_K
       self%phytos(i)%I_S          = pd(list(i))%I_S
       self%phytos(i)%KePHY        = pd(list(i))%KePHY
       self%phytos(i)%f_pr         = pd(list(i))%f_pr
       self%phytos(i)%R_resp       = pd(list(i))%R_resp/secs_pr_day
       self%phytos(i)%theta_resp   = pd(list(i))%theta_resp
       self%phytos(i)%k_fres       = pd(list(i))%k_fres
       self%phytos(i)%k_fdom       = pd(list(i))%k_fdom
       self%phytos(i)%salTol       = pd(list(i))%salTol
       self%phytos(i)%S_bep        = pd(list(i))%S_bep
       self%phytos(i)%S_maxsp      = pd(list(i))%S_maxsp
       self%phytos(i)%S_opt        = pd(list(i))%S_opt
       self%phytos(i)%simDINUptake = pd(list(i))%simDINUptake
       self%phytos(i)%simDONUptake = pd(list(i))%simDONUptake
       self%phytos(i)%simNFixation = pd(list(i))%simNFixation
       self%phytos(i)%simINDynamics= pd(list(i))%simINDynamics
       self%phytos(i)%N_o          = pd(list(i))%N_o
       self%phytos(i)%K_N          = pd(list(i))%K_N
       self%phytos(i)%X_ncon       = pd(list(i))%X_ncon
       self%phytos(i)%X_nmin       = pd(list(i))%X_nmin
       self%phytos(i)%X_nmax       = pd(list(i))%X_nmax
       self%phytos(i)%R_nuptake    = pd(list(i))%R_nuptake/secs_pr_day
       self%phytos(i)%k_nfix       = pd(list(i))%k_nfix
       self%phytos(i)%R_nfix       = pd(list(i))%R_nfix/secs_pr_day
       self%phytos(i)%simDIPUptake = pd(list(i))%simDIPUptake
       self%phytos(i)%simIPDynamics= pd(list(i))%simIPDynamics
       self%phytos(i)%P_0          = pd(list(i))%P_0
       self%phytos(i)%K_P          = pd(list(i))%K_P
       self%phytos(i)%X_pcon       = pd(list(i))%X_pcon
       self%phytos(i)%X_pmin       = pd(list(i))%X_pmin
       self%phytos(i)%X_pmax       = pd(list(i))%X_pmax
       self%phytos(i)%R_puptake    = pd(list(i))%R_puptake/secs_pr_day
       self%phytos(i)%simSiUptake  = pd(list(i))%simSiUptake
       self%phytos(i)%Si_0         = pd(list(i))%Si_0
       self%phytos(i)%K_Si         = pd(list(i))%K_Si
       self%phytos(i)%X_sicon      = pd(list(i))%X_sicon

       ! Register group as a state variable
       CALL self%register_state_variable(self%id_p(i),                         &
                              TRIM(self%phytos(i)%p_name),                     &
                              'mmol/m**3', 'phytoplankton',                    &
                              pd(list(i))%p_initial,                           &
                              minimum=pd(list(i))%p0,                          &
                              vertical_movement = self%phytos(i)%w_p)


       IF (self%phytos(i)%simINDynamics /= 0) THEN
          IF(self%phytos(i)%simINDynamics == 1)THEN
            minNut = self%phytos(i)%p0*self%phytos(i)%X_ncon
          ELSE
            minNut = self%phytos(i)%p0*self%phytos(i)%X_nmin
          ENDIF
          ! Register IN group as a state variable
          CALL self%register_state_variable(self%id_in(i),                     &
                              TRIM(self%phytos(i)%p_name)//'_IN',              &
                              'mmol/m**3', 'phytoplankton IN',                 &
                              pd(list(i))%p_initial*self%phytos(i)%X_ncon,     &
                              minimum=minNut,                                  &
                              vertical_movement = self%phytos(i)%w_p)

       ENDIF
       IF (self%phytos(i)%simIPDynamics /= 0) THEN
          IF(self%phytos(i)%simIPDynamics == 1)THEN
            minNut = self%phytos(i)%p0*self%phytos(i)%X_pcon
          ELSE
            minNut = self%phytos(i)%p0*self%phytos(i)%X_pmin
          ENDIF
          ! Register IP group as a state variable
          CALL self%register_state_variable(self%id_ip(i),                     &
                              TRIM(self%phytos(i)%p_name)//'_IP',              &
                              'mmol/m**3', 'phytoplankton IP',                 &
                              pd(list(i))%p_initial*self%phytos(i)%X_pcon,     &
                              minimum=minNut,                                  &
                              vertical_movement = self%phytos(i)%w_p)

       ENDIF

       CALL self%register_diagnostic_variable(self%id_NtoP(i), TRIM(self%phytos(i)%p_name)//'_NtoP','mmol/m**3', 'INi/IPi', &
                     time_treatment=time_treatment_step_integrated)

       IF (extra_debug) THEN
          CALL self%register_diagnostic_variable(self%id_fT(i), TRIM(self%phytos(i)%p_name)//'_fT','', 'fT', &
                     time_treatment=time_treatment_step_integrated)
          CALL self%register_diagnostic_variable(self%id_fI(i), TRIM(self%phytos(i)%p_name)//'_fI','', 'fI', &
                     time_treatment=time_treatment_step_integrated)
          CALL self%register_diagnostic_variable(self%id_fNit(i), TRIM(self%phytos(i)%p_name)//'_fNit','', 'fNit', &
                     time_treatment=time_treatment_step_integrated)
          CALL self%register_diagnostic_variable(self%id_fPho(i), TRIM(self%phytos(i)%p_name)//'_fPho','', 'fPho', &
                     time_treatment=time_treatment_step_integrated)
          CALL self%register_diagnostic_variable(self%id_fSil(i), TRIM(self%phytos(i)%p_name)//'_fSil','', 'fSil', &
                     time_treatment=time_treatment_step_integrated)
          CALL self%register_diagnostic_variable(self%id_fSal(i), TRIM(self%phytos(i)%p_name)//'_fSal','', 'fSal', &
                     time_treatment=time_treatment_step_integrated)
       ENDIF
    ENDDO

    RETURN

99 call fatal_error('aed_phytoplankton_load_params','Error reading namelist phyto_data')
!
END SUBROUTINE aed_phytoplankton_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION aed_phytoplankton_create(namlst,name,parent) RESULT(self)
!-------------------------------------------------------------------------------
! Initialise the phytoplankton biogeochemical model
!
!  Here, the aed_p_m namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                             :: namlst
   CHARACTER(len=*),INTENT(in)              :: name
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_aed_phytoplankton),POINTER :: self

   INTEGER            :: num_phytos
   INTEGER            :: the_phytos(MAX_PHYTO_TYPES)
   CHARACTER(len=64)  :: p_excretion_target_variable=''
   CHARACTER(len=64)  :: p_mortality_target_variable=''
   CHARACTER(len=64)  :: p1_uptake_target_variable=''
   CHARACTER(len=64)  :: p2_uptake_target_variable=''
   CHARACTER(len=64)  :: n_excretion_target_variable=''
   CHARACTER(len=64)  :: n_mortality_target_variable=''
   CHARACTER(len=64)  :: n1_uptake_target_variable=''
   CHARACTER(len=64)  :: n2_uptake_target_variable=''
   CHARACTER(len=64)  :: n3_uptake_target_variable=''
   CHARACTER(len=64)  :: n4_uptake_target_variable=''
   CHARACTER(len=64)  :: c_excretion_target_variable=''
   CHARACTER(len=64)  :: c_mortality_target_variable=''
   CHARACTER(len=64)  :: c_uptake_target_variable=''
   CHARACTER(len=64)  :: do_uptake_target_variable=''
   CHARACTER(len=64)  :: si_excretion_target_variable=''
   CHARACTER(len=64)  :: si_mortality_target_variable=''
   CHARACTER(len=64)  :: si_uptake_target_variable=''


   real(rk),PARAMETER :: secs_pr_day = 86400.
   real(rk)           :: zerolimitfudgefactor = 0.9 * 3600
   NAMELIST /aed_phytoplankton/ num_phytos, the_phytos,                        &
                    p_excretion_target_variable,p_mortality_target_variable,   &
                     p1_uptake_target_variable, p2_uptake_target_variable,     &
                    n_excretion_target_variable,n_mortality_target_variable,   &
                     n1_uptake_target_variable,n2_uptake_target_variable,      &
                     n3_uptake_target_variable,n4_uptake_target_variable,      &
                    c_excretion_target_variable,c_mortality_target_variable,   &
                      c_uptake_target_variable, do_uptake_target_variable,     &
                    si_excretion_target_variable,si_mortality_target_variable, &
                      si_uptake_target_variable,                               &
                    zerolimitfudgefactor, extra_debug
!-----------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL initialize_model_info(self,name,parent)

   ! Read the namelist
   read(namlst,nml=aed_phytoplankton,err=99)
   dtlim = zerolimitfudgefactor

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   CALL aed_phytoplankton_load_params(self, num_phytos, the_phytos)

   CALL aed_bio_temp_function(self%num_phytos,              &
                              self%phytos%theta_growth,     &
                              self%phytos%T_std,            &
                              self%phytos%T_opt,            &
                              self%phytos%T_max,            &
                              self%phytos%aTn,              &
                              self%phytos%bTn,              &
                              self%phytos%kTn,              &
                              self%phytos%p_name)


   ! Register link to nutrient pools, if variable names are provided in namelist.
   self%do_Pexc = p_excretion_target_variable .NE. ''
   IF (self%do_Pexc) THEN
     CALL self%register_state_dependency(self%id_Pexctarget, p_excretion_target_variable)
   ENDIF
   self%do_Nexc = n_excretion_target_variable .NE. ''
   IF (self%do_Pexc) THEN
     CALL self%register_state_dependency(self%id_Nexctarget, n_excretion_target_variable)
   ENDIF
   self%do_Cexc = c_excretion_target_variable .NE. ''
   IF (self%do_Pexc) THEN
     CALL self%register_state_dependency(self%id_Cexctarget, c_excretion_target_variable)
   ENDIF
   self%do_Siexc = si_excretion_target_variable .NE. ''
   IF (self%do_Siexc) THEN
     CALL self%register_state_dependency(self%id_Siexctarget, si_excretion_target_variable)
   ENDIF

   self%do_Pmort = p_mortality_target_variable .NE. ''
   IF (self%do_Pmort) THEN
     CALL self%register_state_dependency(self%id_Pmorttarget, p_mortality_target_variable)
   ENDIF
   self%do_Nmort = n_mortality_target_variable .NE. ''
   IF (self%do_Nmort) THEN
     CALL self%register_state_dependency(self%id_Nmorttarget, n_mortality_target_variable)
   ENDIF
   self%do_Cmort = c_mortality_target_variable .NE. ''
   IF (self%do_Cmort) THEN
     CALL self%register_state_dependency(self%id_Cmorttarget, c_mortality_target_variable)
   ENDIF
   self%do_Simort = si_mortality_target_variable .NE. ''
   IF (self%do_Simort) THEN
     CALL self%register_state_dependency(self%id_Simorttarget, si_mortality_target_variable)
   ENDIF

   self%npup = 0
   IF (p1_uptake_target_variable .NE. '') self%npup = 1
   IF (p2_uptake_target_variable .NE. '') self%npup = 2
   self%do_Puptake = .FALSE.
   IF (self%npup>0) self%do_Puptake=.TRUE.
   IF (self%do_Puptake) THEN
     IF (self%npup>0) CALL self%register_state_dependency(self%id_Pupttarget(1), p1_uptake_target_variable); ifrp=1
     IF (self%npup>1) CALL self%register_state_dependency(self%id_Pupttarget(2), p2_uptake_target_variable); idop=2
   ENDIF
   self%nnup = 0
   IF (n1_uptake_target_variable .NE. '') self%nnup = 1
   IF (n2_uptake_target_variable .NE. '') self%nnup = 2
   IF (n3_uptake_target_variable .NE. '') self%nnup = 3
   IF (n4_uptake_target_variable .NE. '') self%nnup = 4
   self%do_Nuptake = .false.
   IF (self%nnup>0) self%do_Nuptake=.true.
   IF (self%do_Nuptake) THEN
     IF (self%nnup>0) CALL self%register_state_dependency(self%id_Nupttarget(1), n1_uptake_target_variable); ino3=1
     IF (self%nnup>1) CALL self%register_state_dependency(self%id_Nupttarget(2), n2_uptake_target_variable); inh4=2
     IF (self%nnup>2) CALL self%register_state_dependency(self%id_Nupttarget(3), n3_uptake_target_variable); idon=3
     IF (self%nnup>3) CALL self%register_state_dependency(self%id_Nupttarget(4), n4_uptake_target_variable); in2 =4
   ENDIF
   self%do_Cuptake = c_uptake_target_variable .NE. ''
   IF (self%do_Cuptake) THEN
     CALL self%register_state_dependency(self%id_Cupttarget, c_uptake_target_variable)
   ENDIF
   self%do_DOuptake = do_uptake_target_variable .NE. ''
   IF (self%do_DOuptake) THEN
     CALL self%register_state_dependency(self%id_DOupttarget, do_uptake_target_variable)
   ENDIF
   self%do_Siuptake = si_uptake_target_variable .NE. ''
   IF (self%do_Siuptake) THEN
     CALL self%register_state_dependency(self%id_Siupttarget, si_uptake_target_variable)
   ENDIF

   ! Register diagnostic variables
   CALL self%register_diagnostic_variable(self%id_GPP, 'GPP','mmol/m**3',  'gross primary production',           &
                     time_treatment=time_treatment_step_integrated)
   CALL self%register_diagnostic_variable(self%id_NCP, 'NCP','mmol/m**3',  'net community production',           &
                     time_treatment=time_treatment_step_integrated)
   CALL self%register_diagnostic_variable(self%id_PPR, 'PPR','mmol/m**3/d','gross primary production rate',      &
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_NPR, 'NPR','mmol/m**3/d','net community production rate',      &
                     time_treatment=time_treatment_averaged)

   CALL self%register_diagnostic_variable(self%id_NUP, 'NUP','mmol/m**3/d','nitrogen uptake',    &
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_PUP, 'PUP','mmol/m**3/d','phosphorous uptake', &
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_CUP, 'CUP','mmol/m**3/d','carbon uptake',      &
                     time_treatment=time_treatment_averaged)

   CALL self%register_diagnostic_variable(self%id_dPAR, 'PAR','W/m**2',  'photosynthetically active radiation',&
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_TCHLA, 'TCHLA','ug/L', 'Total Chlorophyll-a',&
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_TPHY, 'TPHYS','ug/L',  'Total Phytoplankton',&
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_TIN, 'IN','ug/L',      'Total Chlorophyll-a',&
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_TIP, 'IP','ug/L',      'Total Chlorophyll-a',&
                     time_treatment=time_treatment_averaged)

   ! Register conserved quantities
   CALL self%register_conserved_quantity(self%id_totP, 'TPHY','mmol/m**3','phytoplankton')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_tem,  standard_variables%temperature)
   CALL self%register_dependency(self%id_sal,  standard_variables%practical_salinity)
   CALL self%register_dependency(self%id_par,  standard_variables%downwelling_photosynthetic_radiative_flux)
   CALL self%register_horizontal_dependency_sn(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   CALL self%register_dependency(self%id_dz,   standard_variables%cell_thickness)
   CALL self%register_dependency(self%id_extc, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)

   RETURN

99 call fatal_error('aed_phytoplankton_init','Error reading namelist aed_phytoplankton')

END FUNCTION aed_phytoplankton_create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of phytoplankton biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   real(rk)           :: phy, tphy, tin, tip, tchla
   real(rk)           :: INi, IPi
   real(rk)           :: pup
   real(rk)           :: no3up,nh4up
   real(rk)           :: cup, rsiup
   real(rk)           :: temp,par,Io,salinity, extc,dz
   real(rk)           :: primprod(self%num_phytos), exudation(self%num_phytos), &
                         a_nfix(self%num_phytos), respiration(self%num_phytos)
   real(rk)           :: cuptake(self%num_phytos), cexcretion(self%num_phytos), cmortality(self%num_phytos)
   real(rk)           :: nuptake(self%num_phytos,1:4), nexcretion(self%num_phytos), nmortality(self%num_phytos)
   real(rk)           :: puptake(self%num_phytos,1:2), pexcretion(self%num_phytos), pmortality(self%num_phytos)
   real(rk)           :: siuptake(self%num_phytos), siexcretion(self%num_phytos), simortality(self%num_phytos)
   real(rk)           :: fT, fNit, fPho, fSil, fI, fXl, fSal, PNf
   real(rk)           :: upTot

   INTEGER            :: phy_i,c
   real(rk)           :: flux, available

! MH to fix
!  real(rk)           :: dt = 3600. ! just for now, hard code it
   real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current environmental conditions.
   _GET_   (self%id_tem,temp)     ! local temperature
   _GET_   (self%id_sal,salinity) ! local salinity
   _GET_   (self%id_par,par)      ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,Io)       ! surface short wave radiation

   pup = 0.
   ! Retrieve current (local) state variable values.
   IF (self%do_Puptake) _GET_(self%id_Pupttarget(1), pup)

   no3up = 0.
   nh4up = 0.
   IF (self%do_Nuptake) THEN
      _GET_(self%id_Nupttarget(1), no3up)
      _GET_(self%id_Nupttarget(2), nh4up)
   ENDIF
   cup = 0.
   IF (self%do_Cuptake) _GET_(self%id_Cupttarget, cup)
   rsiup = 0.
   IF (self%do_Siuptake) _GET_(self%id_Siupttarget, rsiup)

   tphy = 0.0
   tchla = 0.0
   tin  = 0.0
   tip  = 0.0

   INi = 0.
   IPi = 0.

   DO phy_i=1,self%num_phytos

!     ! Available nutrients - this is set to the concentration of nutrients from
!     ! other aed_modules unless they are limiting
!     no3up = sum(nuptake(phy_i, 1) / dtlim
!     nh4up = nuptake(phy_i, 2) / dtlim
!     pup   = puptake(phy_i, 1) / dtlim
!     cup   = cuptake(phy_i)    / dtlim
!     rsiup = siuptake(phy_i)   / dtlim

      primprod(phy_i)    = 0.0_rk
      exudation(phy_i)   = 0.0_rk
      a_nfix(phy_i)      = 0.0_rk
      respiration(phy_i) = 0.0_rk

      cuptake(phy_i)     = 0.0_rk
      cexcretion(phy_i)  = 0.0_rk
      cmortality(phy_i)  = 0.0_rk
      nuptake(phy_i,:)   = 0.0_rk
      nexcretion(phy_i)  = 0.0_rk
      nmortality(phy_i)  = 0.0_rk
      puptake(phy_i,:)   = 0.0_rk
      pexcretion(phy_i)  = 0.0_rk
      pmortality(phy_i)  = 0.0_rk

      ! Retrieve this phytoplankton group
      _GET_(self%id_p(phy_i),phy)

      ! Get the temperature limitation function
      fT = fTemp_function(self%phytos(phy_i)%fT_Method,    &
                          self%phytos(phy_i)%T_max,        &
                          self%phytos(phy_i)%T_std,        &
                          self%phytos(phy_i)%theta_growth, &
                          self%phytos(phy_i)%aTn,          &
                          self%phytos(phy_i)%bTn,          &
                          self%phytos(phy_i)%kTn,temp)

      ! Get the light and nutrient limitation.
      ! NITROGEN.
      fNit = 0.0
      IF(self%phytos(phy_i)%simINDynamics /= 0) THEN
         ! IN variable available
         _GET_(self%id_in(phy_i),INi)
      ELSE
         ! Assumed constant IN:
         INi = phy*self%phytos(phy_i)%X_ncon
      END IF

      ! Estimate fN limitation from IN or ext N value
      IF(self%phytos(phy_i)%simINDynamics > 1) THEN
         IF (phy > self%phytos(phy_i)%p0) THEN
            fNit = INi / phy
            fNit = phyto_fN(self,phy_i,IN=fNit)
         ENDIF
         IF (phy > 0.0_rk .AND. phy <= self%phytos(phy_i)%p0) THEN
            fNit = phyto_fN(self,phy_i,din=no3up+nh4up)
         ENDIF
      ELSE
         fNit = phyto_fN(self,phy_i,din=no3up+nh4up)
      ENDIF
      IF (self%phytos(phy_i)%simNFixation /= 0) THEN
         ! Nitrogen fixer: apply no N limitation. N Fixation ability
         ! depends on DIN concentration
         a_nfix = (1.0_rk - fNit)
         fNit = 1.0_rk
      ENDIF


      ! PHOSPHOROUS.
      fPho = 0.0_rk
      IF (self%phytos(phy_i)%simIPDynamics /= 0) THEN
         ! IP variable available
         _GET_(self%id_ip(phy_i),IPi)
      ELSE
         ! Assumed constant IP:
         IPi = phy*self%phytos(phy_i)%X_pcon
      END IF

      ! Estimate fP limitation from IP or ext P value
      IF (self%phytos(phy_i)%simIPDynamics > 1) THEN
         IF (phy > self%phytos(phy_i)%p0) THEN
            fPho = IPi / phy
            fPho = phyto_fP(self,phy_i,IP=fPho)
         ENDIF
         IF (phy > 0.0_rk .AND. phy <= self%phytos(phy_i)%p0) THEN
            fPho = phyto_fP(self,phy_i,frp=pup)
         ENDIF
      ELSE
         fPho = phyto_fP(self,phy_i,frp=pup)
      ENDIF

      ! SILICA.
      fSil = phyto_fSi(self,phy_i,rsiup)


      ! LIGHT
      _GET_(self%id_extc,extc)
      ! dz = 0.5     !MH: to fix
      _GET_(self%id_dz,dz)
      fI = phyto_light(self, phy_i, par, extc, Io, dz)
      ! fI = 0.1


      ! METAL AND TOXIC EFFECTS
      fXl = 1.0

      ! Primary production rate
      primprod(phy_i) = self%phytos(phy_i)%R_growth * fT * findMin(fI,fNit,fPho,fSil) * fxl

      ! Adjust primary production rate for nitrogen fixers
      IF (self%phytos(phy_i)%simNFixation /= 0) THEN
         ! Nitrogen fixing species, and the growth rate to  must be reduced
         ! to compensate for the increased metabolic cost of this process
         primprod(phy_i) = primprod(phy_i) * (self%phytos(phy_i)%k_nfix + &
                           (1.0-a_nfix(phy_i))*(1.0-self%phytos(phy_i)%k_nfix))
      ENDIF


      ! Respiration and general metabolic loss

      respiration(phy_i) = phyto_respiration(self,phy_i,temp)

      ! Salinity stress effect on respiration
      fSal =  phyto_salinity(self,phy_i,salinity)
      respiration(phy_i) = respiration(phy_i) * fSal

      ! photo-exudation
      exudation(phy_i) = primprod(phy_i)*self%phytos(phy_i)%f_pr

      ! Limit respiration if at the min biomass to prevent
      ! leak in the C mass balance
      IF (phy <= self%phytos(phy_i)%p0) THEN
         respiration(phy_i) = 0.0_rk
         exudation(phy_i) = 0.0_rk
      ENDIF

      ! write(*,"(4X,'limitations (fT,fI,fN,fP,fSi,Io, par, mu): ',9F9.2)")fT,fI,fNit,fPho,fSil,Io,par,primprod*secs_pr_day


      ! Carbon uptake and excretion

      cuptake(phy_i)    = -primprod(phy_i) * phy
      cexcretion(phy_i) = (self%phytos(phy_i)%k_fdom*(1.0-self%phytos(phy_i)%k_fres)*respiration(phy_i)+exudation(phy_i)) * phy
      cmortality(phy_i) = ((1.0-self%phytos(phy_i)%k_fdom)*(1.0-self%phytos(phy_i)%k_fres)*respiration(phy_i)) * phy

      ! Nitrogen uptake and excretion

      CALL phyto_internal_nitrogen(self,phy_i,phy,INi,primprod(phy_i),&
                             fT,no3up,nh4up,a_nfix(phy_i),respiration(phy_i),exudation(phy_i),PNf,&
                                   nuptake(phy_i,:),nexcretion(phy_i),nmortality(phy_i))

      ! Phosphorus uptake and excretion

      CALL phyto_internal_phosphorus(self,phy_i,phy,IPi,primprod(phy_i),&
                                 fT,pup,respiration(phy_i),exudation(phy_i),&
                                         puptake(phy_i,:),pexcretion(phy_i),pmortality(phy_i))

      ! Silica uptake and excretion

      IF (self%phytos(phy_i)%simSiUptake > 0) THEN
         siuptake(phy_i)    =-self%phytos(phy_i)%X_sicon * primprod(phy_i) * phy
         siexcretion(phy_i) = self%phytos(phy_i)%X_sicon * (self%phytos(phy_i)%k_fdom*respiration(phy_i)+exudation(phy_i)) * phy
         simortality(phy_i) = self%phytos(phy_i)%X_sicon * ((1.0-self%phytos(phy_i)%k_fdom)*respiration(phy_i)) * phy
      ELSE
         siuptake(phy_i)    = 0.0_rk
         siexcretion(phy_i) = 0.0_rk
         simortality(phy_i) = 0.0_rk
      ENDIF

      ! Diagnostic info

      _SET_DIAGNOSTIC_(self%id_NtoP(phy_i), INi/IPi)

      IF (extra_debug) THEN
         _SET_DIAGNOSTIC_(self%id_fT(phy_i), fT)
         _SET_DIAGNOSTIC_(self%id_fI(phy_i), fI)
         _SET_DIAGNOSTIC_(self%id_fNit(phy_i), fNit)
         _SET_DIAGNOSTIC_(self%id_fPho(phy_i), fPho)
         _SET_DIAGNOSTIC_(self%id_fSil(phy_i), fSil)
         _SET_DIAGNOSTIC_(self%id_fSal(phy_i), fSal)
      ENDIF
   END DO


   !-----------------------------------------------------------------
   ! Check uptake values for availability to prevent -ve numbers

   ! pup   - p available
   ! no3up - no3 available
   ! nh4up - nh4 available
   ! cup   - c available
   ! rsiup - Si available

   IF (self%do_Puptake) THEN
      upTot = sum(puptake(:,1))*dtlim
      IF ( upTot >= pup ) THEN
         DO phy_i=1,self%num_phytos
            puptake(phy_i,1) = (pup*0.99/dtlim) * (puptake(phy_i,1)/upTot)
         ENDDO
      ENDIF
   ENDIF

   IF (self%do_Nuptake) THEN
      upTot = sum(nuptake(:,1))*dtlim
      IF ( upTot >= no3up ) THEN
         DO phy_i=1,self%num_phytos
            nuptake(phy_i,1) = (no3up*0.99/dtlim) * (nuptake(phy_i,1)/upTot)
         ENDDO
      ENDIF

      upTot = sum(nuptake(:,2))*dtlim
      IF ( upTot >= nh4up ) THEN
         DO phy_i=1,self%num_phytos
            nuptake(phy_i,2) = (nh4up*0.99/dtlim) * (nuptake(phy_i,2)/upTot)
         ENDDO
      ENDIF
   ENDIF
   IF (self%do_Cuptake) THEN
      upTot = sum(cuptake)*dtlim
      IF ( upTot >= cup ) THEN
         DO phy_i=1,self%num_phytos
            cuptake(phy_i) = (cup*0.99/dtlim) * (cuptake(phy_i)/upTot)
         ENDDO
      ENDIF
   ENDIF
!  IF (self%do_DOuptake) THEN
!     !
!  ENDIF
   IF (self%do_Siuptake) THEN
      upTot = sum(siuptake)*dtlim
      IF ( upTot >= rsiup ) THEN
         DO phy_i=1,self%num_phytos
            siuptake(phy_i) = (rsiup*0.99/dtlim) * (siuptake(phy_i)/upTot)
         ENDDO
      ENDIF
   ENDIF

   DO phy_i=1,self%num_phytos

      ! Retrieve this phytoplankton group
      !-----------------------------------------------------------------
      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Phytoplankton production / losses
      _GET_(self%id_p(phy_i),phy)
      flux = (primprod(phy_i) - respiration(phy_i) - exudation(phy_i)) * phy
      available = MAX(0.0_rk, phy - self%phytos(phy_i)%p0)
      IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
      _SET_ODE_(self%id_p(phy_i), flux)

      IF (self%phytos(phy_i)%simINDynamics /= 0) THEN
         ! _SET_ODE_(self%id_in(phy_i), (-sum(nuptake) - nexcretion(phy_i) - nmortality(phy_i) )*INi )
         _GET_(self%id_in(phy_i),INi)
         flux = (-sum(nuptake(phy_i,:)) - nexcretion(phy_i) - nmortality(phy_i) )
         available = MAX(0.0_rk, INi - self%phytos(phy_i)%X_nmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _SET_ODE_(self%id_in(phy_i), flux)
      ENDIF
      IF (self%phytos(phy_i)%simIPDynamics /= 0) THEN
         ! _SET_ODE_(self%id_ip(phy_i), (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) ) )
         _GET_(self%id_ip(phy_i),IPi)
         flux = (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) )
         available = MAX(0.0_rk, IPi - self%phytos(phy_i)%X_pmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _SET_ODE_(self%id_ip(phy_i), flux)
      ENDIF

      ! Now manage uptake of nutrients, CO2 and DO - these cumulative fluxes already limited above loop
      IF (self%do_Puptake) THEN
         DO c = 1,self%npup
            _SET_ODE_(self%id_Pupttarget(c), puptake(phy_i,c))
         ENDDO
      ENDIF
      IF (self%do_Nuptake) THEN
         DO c = 1,self%nnup
            _SET_ODE_(self%id_Nupttarget(c), nuptake(phy_i,c))
         ENDDO
      ENDIF
      IF (self%do_Cuptake) THEN
         _SET_ODE_(self%id_Cupttarget,  cuptake(phy_i) - respiration(phy_i)*self%phytos(phy_i)%k_fres*phy )
      ENDIF
      IF (self%do_DOuptake) THEN
         _SET_ODE_(self%id_DOupttarget, -cuptake(phy_i) + respiration(phy_i)*self%phytos(phy_i)%k_fres*phy )
      ENDIF
      IF (self%do_Siuptake) THEN
         _SET_ODE_(self%id_Siupttarget, siuptake(phy_i))
      ENDIF
      ! Now manage mortality contributions to POM
      IF (self%do_Pmort) THEN
         _SET_ODE_(self%id_Pmorttarget,pmortality(phy_i))
      ENDIF
      IF (self%do_Nmort) THEN
         _SET_ODE_(self%id_Nmorttarget,nmortality(phy_i))
      ENDIF
      IF (self%do_Cmort) THEN
         _SET_ODE_(self%id_Cmorttarget,cmortality(phy_i))
      ENDIF
      IF (self%do_Simort) THEN
         _SET_ODE_(self%id_Simorttarget,simortality(phy_i))
      ENDIF
      ! Now manage excretion/exudation contributions to DOM
      IF (self%do_Pexc) THEN
         _SET_ODE_(self%id_Pexctarget,pexcretion(phy_i))
      ENDIF
      IF (self%do_Nexc) THEN
         _SET_ODE_(self%id_Nexctarget,nexcretion(phy_i))
      ENDIF
      IF (self%do_Cexc) THEN
         _SET_ODE_(self%id_Cexctarget,cexcretion(phy_i))
      ENDIF
      IF (self%do_Siexc) THEN
         _SET_ODE_(self%id_Siexctarget,siexcretion(phy_i))
      ENDIF

      !-----------------------------------------------------------------
      ! export diagnostic variables

      ! Total phytoplankton carbon
      tphy = tphy + phy

      ! Total chlorophyll-a
      IF (self%phytos(phy_i)%Xcc > 0.1) THEN
        ! Assume Xcc (mol C/ mol chla) is a constant
        tchla = tchla + ( phy / self%phytos(phy_i)%Xcc ) * 12.0
      ELSE
        ! Use dynamic equation (Eq 13: of Baklouti, Cloern et al. 1995)
        ! theta = 1/Xcc [mg Chl (mg C)1] = 0.003 + 0.0154  e^0.050T  e^0.059E mu
        tchla = tchla + ( phy * (0.003 + 0.0154 * exp(0.050*temp) * exp(0.059*par) &
                        * primprod(phy_i)))
      ENDIF

      ! Total internal nutrients
      tin = tin + INi
      tip = tip + IPi

   ENDDO

   _SET_DIAGNOSTIC_(self%id_GPP, sum(primprod))
   _SET_DIAGNOSTIC_(self%id_NCP, sum(primprod - respiration))
   _SET_DIAGNOSTIC_(self%id_PPR, sum(primprod) * secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR, sum(primprod - respiration) * secs_pr_day)

   _SET_DIAGNOSTIC_(self%id_NUP, sum(nuptake))
   _SET_DIAGNOSTIC_(self%id_PUP, sum(puptake))
   _SET_DIAGNOSTIC_(self%id_CUP, sum(cuptake))


   _SET_DIAGNOSTIC_(self%id_dPAR, par)
   _SET_DIAGNOSTIC_(self%id_TCHLA, tchla)
   _SET_DIAGNOSTIC_(self%id_TPHY, tphy)
   _SET_DIAGNOSTIC_(self%id_TIN, tin)
   _SET_DIAGNOSTIC_(self%id_TIP, tip)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_phytoplankton_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE phyto_internal_phosphorus(self,group,phy,IP,primprod,&
                                                 fT,pup,respiration,exudation,&
                                                     uptake,excretion,mortality)
!-------------------------------------------------------------------------------
! Calculates the internal phosphorus stores and fluxes
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   INTEGER,INTENT(in)                          :: group
   real(rk),INTENT(in)                         :: phy
   real(rk),INTENT(in)                         :: IP
   real(rk),INTENT(in)                         :: primprod
   real(rk),INTENT(in)                         :: fT,pup,respiration,exudation
   real(rk),INTENT(out)                        :: uptake(:),excretion,mortality
!CONSTANTS
   real(rk),PARAMETER :: one_e_neg5 = 1e-5
!LOCALS
   real(rk) :: dumdum1,dumdum2,theX_pcon
   INTEGER  :: c
!
!-------------------------------------------------------------------------------
!BEGIN
   uptake     = 0.0_rk
   excretion  = 0.0_rk
   mortality  = 0.0_rk

   ! Uptake of phosphorus
   IF (self%phytos(group)%simIPDynamics == 0 .OR. self%phytos(group)%simIPDynamics == 1) THEN

      ! Static phosphorus uptake function
      ! uptake = X_pcon * mu * phy

      theX_pcon = self%phytos(group)%X_pcon * phy
      DO c = 1,self%npup
         ! uptake is spread over relevant sources (assumes evenly)
         uptake(c) = - (theX_pcon/self%npup) * primprod
      END DO

   ELSEIF (self%phytos(group)%simIPDynamics == 2) THEN

      ! Dynamic phosphorus uptake function
      ! uptake = R_puptake*fT*phy* (X_pmax-IP/phy)/(X_pmax-X_pmin) * (PO4/(K_P + PO4))

      theX_pcon = IP
      dumdum1  = self%phytos(group)%R_puptake * fT * phy
      dumdum2  = MAX(one_e_neg5, self%phytos(group)%X_pmax - (IP / phy))
      dumdum1  = dumdum1 * dumdum2 / (self%phytos(group)%X_pmax-self%phytos(group)%X_pmin)
      uptake(1)= -dumdum1 * phyto_fP(self,group,frp=pup)
      uptake(2)= 0.0

   ELSE

      ! Unknown phosphorus uptake function
      print *,'STOP: unknown simIPDynamics (',self%phytos(group)%simIPDynamics,') for: ',self%phytos(group)%p_name
      STOP

   ENDIF

   ! Release of phosphorus due to excretion from phytoplankton and
   ! contribution of mortality and excretion to OM

   excretion = (respiration*self%phytos(group)%k_fdom + exudation)*theX_pcon
   mortality = respiration*(1.0-self%phytos(group)%k_fdom)*theX_pcon


END SUBROUTINE phyto_internal_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE phyto_internal_nitrogen(self,group,phy,IN,primprod,fT,no3up,nh4up,   &
                                   a_nfix,respiration,exudation,PNf,            &
                                   uptake,excretion,mortality)
!-------------------------------------------------------------------------------
! Calculates the internal nitrogen stores and fluxes
!-------------------------------------------------------------------------------

!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   INTEGER,INTENT(in)                          :: group
   real(rk),INTENT(in)                         :: phy
   real(rk),INTENT(in)                         :: IN
   real(rk),INTENT(in)                         :: primprod
   real(rk),INTENT(in)                         :: fT,no3up,nh4up
   real(rk),INTENT(out)                        :: a_nfix
   real(rk),INTENT(in)                         :: respiration,exudation
   real(rk),INTENT(out)                        :: PNf
   real(rk),INTENT(out)                        :: uptake(:),excretion,mortality
!
!CONSTANTS
   real(rk),PARAMETER :: one_e_neg5 = 1e-5
!
!LOCALS
   real(rk)  :: dumdum1,dumdum2,theX_ncon
!
!-------------------------------------------------------------------------------
!BEGIN
   uptake     = 0.0_rk
   excretion  = 0.0_rk
   mortality  = 0.0_rk


   ! Uptake of nitrogen
   IF (self%phytos(group)%simINDynamics == 0 .OR. self%phytos(group)%simINDynamics == 1) THEN

      ! Static nitrogen uptake function (assuming fixed stoichiometry)
      ! uptake = X_ncon * mu * phy

      theX_ncon = self%phytos(group)%X_ncon * phy
      uptake(1)  = -theX_ncon * primprod

   ELSEIF (self%phytos(group)%simINDynamics == 2) THEN

      ! Dynamic nitrogen uptake function
      ! uptake = R_nuptake*fT*phy* (X_nmax-IN/phy)/(X_nmax-X_nmin) * (DIN/(K_N + DIN))

      theX_ncon = IN
      dumdum1  = self%phytos(group)%R_nuptake * fT * phy
      dumdum2  = MAX(self%phytos(group)%X_nmax - (IN / phy),one_e_neg5)
      dumdum1  = dumdum1 * dumdum2 / (self%phytos(group)%X_nmax-self%phytos(group)%X_nmin)
      uptake(1) = dumdum1 * phyto_fN(self,group,din=no3up+nh4up)
      uptake(1) = -uptake(1)
   ELSE
      ! Unknown nitrogen uptake function
      print *,'STOP: unknown simINDynamics (',self%phytos(group)%simINDynamics,') for: ',self%phytos(group)%p_name
      STOP
   ENDIF

   ! Now find out how much of this was due to N Fixation:
   IF (self%phytos(group)%simNFixation /= 0) THEN
      a_nfix = self%phytos(group)%R_nfix * a_nfix * phy
      IF (a_nfix > uptake(1)) THEN
         ! Extreme case:
         a_nfix = uptake(1)
         uptake(1) = 0.0_rk
      ELSE
         ! Reduce nuptake by the amount fixed:
         uptake(1) = uptake(1) * (uptake(1)-a_nfix) / uptake(1)
      ENDIF
   ENDIF

   ! Disaggregate N sources to NO3, NH4, DON and N2, based on configuraiton
   PNf = phyto_pN(self,group,nh4up,no3up)

   IF (self%phytos(group)%simDINUptake /= 0) THEN
      uptake(inh4) = uptake(1) * (1.0-PNf) !inh4 == 2
      uptake(ino3) = uptake(1) * PNf       !ino3 == 1
   ENDIF
   IF (self%phytos(group)%simDONUptake /= 0) THEN
      uptake(idon) = 0.0  !MH to fix  (idon == 3)
   ENDIF
   IF (self%phytos(group)%simNFixation /= 0 .AND. self%do_N2uptake) THEN
      uptake(iN2) = a_nfix       ! iN2 == 4
   ENDIF


   ! Release of nitrogen due to excretion from phytoplankton and
   ! contribution of mortality and excretion OM:
   ! (/day +/day)* mg N/ mg C * mgC

   excretion = (respiration*self%phytos(group)%k_fdom + exudation)*theX_ncon
   mortality = respiration*(1.0-self%phytos(group)%k_fdom)*theX_ncon

   ! should check here e or m is not exceeding X_nmin

END SUBROUTINE phyto_internal_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE aed_phytoplankton_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! !IROUTINE: Calculate pelagic sedimentation of phytoplankton.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   real(rk) :: phy        ! State
   INTEGER  :: phy_i
   real(rk) :: phy_flux

   ! Parameters
!  real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN

   _FABM_HORIZONTAL_LOOP_BEGIN_

   DO phy_i=1,self%num_phytos
      ! Retrieve current (local) state variable values.
      _GET_(self%id_p(phy_i),phy) ! phytoplankton

      phy_flux = 0.0_rk  !self%phytos(phy_i)%w_p*MAX(phy,0.0_rk)

     ! Set bottom fluxes for the pelagic (change per surface area per second)
     ! Transfer sediment flux value to FABM.
     _SET_BOTTOM_EXCHANGE_(self%id_p(phy_i),phy_flux)

   ENDDO

   ! Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_

END SUBROUTINE aed_phytoplankton_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_phytoplankton_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
!LOCALS
   real(rk) :: phy
   INTEGER  :: phy_i
!
!-----------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   DO phy_i=1,self%num_phytos
      ! Retrieve current (local) state variable values.
      _GET_(self%id_p(phy_i),phy) ! phytoplankton

      ! Self-shading with explicit contribution from background phytoplankton concentration.
      _SET_EXTINCTION_(self%phytos(phy_i)%KePHY*phy)

   ENDDO

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_phytoplankton_get_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   real(rk) :: phy,p
   INTEGER  :: phy_i
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   phy = 0.0_rk
   DO phy_i=1,self%num_phytos
      ! Retrieve current (local) state variable values.
      _GET_(self%id_p(phy_i),p) ! phytoplankton
      phy = phy + p
   ENDDO

   ! Total phytoplankton is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totP,phy)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_phytoplankton_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of phytoplankton biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   real(rk)            :: n,p,par,I_0
   INTEGER             :: phy_i
   real(rk)            :: iopt,rpd,primprod
   real(rk),PARAMETER  :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   DO phy_i=1,self%num_phytos
      ! Retrieve current (local) state variable values.
      _GET_(self%id_p(phy_i),p) ! phytoplankton
      _GET_(self%id_Pupttarget(1),n) ! nutrients

      ! Retrieve current environmental conditions.
      _GET_   (self%id_par,par)  ! local photosynthetically active radiation

      _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

      ! Light acclimation formulation based on surface light intensity.
      iopt = max(0.25*I_0,self%phytos(phy_i)%I_min)

      ! Loss rate of phytoplankton to detritus depends on local light intensity.
      IF (par .ge. self%phytos(phy_i)%I_min) THEN
         rpd = self%phytos(phy_i)%rpdu
      ELSE
         rpd = self%phytos(phy_i)%rpdl
      ENDIF

      ! NEEDS FIXING

      ! Rate of primary production will be reused multiple times - calculate it once.
      primprod = 0.0 !fnp(self,n,p,par,iopt,phy_i)


   ENDDO

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_phytoplankton_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_fN(self, group, IN, din, don) RESULT(fN)
!-------------------------------------------------------------------------------
! Nitrogen limitation of phytoplankton.
! Michaelis-Menton type formulation or droop model for species with IN
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   INTEGER,INTENT(in)                          :: group
   real(rk),INTENT(in),OPTIONAL                :: IN
   real(rk),INTENT(in),OPTIONAL                :: din
   real(rk),INTENT(in),OPTIONAL                :: don
!
!LOCALS
   real(rk) :: fN
   real(rk) :: nup

!-------------------------------------------------------------------------------
!BEGIN
   fN=1.0_rk

   IF (PRESENT(din) .OR. PRESENT(don)) THEN
     ! Calculate external nutrient limitation factor
     nup = 0.0
     IF (PRESENT(din) .AND. self%phytos(group)%simDINUptake == 1) THEN
       nup = nup + din
     ENDIF
     IF (PRESENT(don) .AND. self%phytos(group)%simDONUptake == 1) THEN
       nup = nup + don
     ENDIF
     fN = (nup-self%phytos(group)%N_o) / &
           (nup-self%phytos(group)%N_o+self%phytos(group)%K_N)
   ELSE
     ! Calculate internal nutrient limitation factor
     fN =   self%phytos(group)%X_nmax*(1.0-self%phytos(group)%X_nmin/IN) / &
            (self%phytos(group)%X_nmax-self%phytos(group)%X_nmin)
   ENDIF

   IF ( fN < 0.0_rk ) fN = 0.0_rk
END FUNCTION phyto_fN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_fP(self, group, IP, frp) RESULT(fP)
!-------------------------------------------------------------------------------
! Phosphorus limitation of phytoplankton
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   INTEGER,INTENT(in)                          :: group
   real(rk),INTENT(in), OPTIONAL               :: IP
   real(rk),INTENT(in), OPTIONAL               :: frp
!
!LOCALS
   real(rk) :: fP
!
!-------------------------------------------------------------------------------
!BEGIN
   fP=1.0_rk

   IF(PRESENT(frp)) THEN
     fP = (frp-self%phytos(group)%P_0) / &
             (self%phytos(group)%K_P + (MAX(0.0_rk, (frp-self%phytos(group)%P_0))))
   ELSE
     fP = self%phytos(group)%X_pmax * (1.0 - self%phytos(group)%X_pmin/IP) / &
                          (self%phytos(group)%X_pmax-self%phytos(group)%X_pmin)
   ENDIF

   IF( fP<0.0_rk ) fP=0.0_rk

END FUNCTION phyto_fP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
FUNCTION phyto_fSi(self, group, Si) RESULT(fSi)
!-------------------------------------------------------------------------------
! Silica limitation (eg. for diatoms)
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   INTEGER,INTENT(in)                          :: group
   real(rk),INTENT(in)                         :: Si
!
!LOCALS
   real(rk) :: fSi
!
!-------------------------------------------------------------------------------
!BEGIN
   fSi = 1.0_rk

   IF (self%phytos(group)%simSiUptake == 1) THEN
     fSi = (Si-self%phytos(group)%Si_0) / &
           (Si-self%phytos(group)%Si_0+self%phytos(group)%K_Si)
     IF ( fSi < 0.0_rk ) fSi=0.0_rk
   ENDIF

END FUNCTION phyto_fSi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
FUNCTION phyto_pN(self,group,NH4,NO3) RESULT(pN)
!-------------------------------------------------------------------------------
! Calculates the relative preference of uptake by phytoplankton of
! ammonia uptake over nitrate.
!-------------------------------------------------------------------------------
   !-- Incoming
   _CLASS_ (type_aed_phytoplankton), INTENT(in) :: self
   INTEGER,INTENT(in)                           :: group
   real(rk),INTENT(IN)                          :: NH4
   real(rk),INTENT(IN)                          :: NO3
!
!LOCALS
   real(rk) :: pN
!
!-------------------------------------------------------------------------------
!BEGIN
   pN = 0.0_rk

   IF (NH4 > 0.0) THEN
      pN = NH4*NO3 / ((NH4+self%phytos(group)%K_N)*(NO3+self%phytos(group)%K_N)) &
         + NH4*self%phytos(group)%K_N / ((NH4+NO3)*(NO3+self%phytos(group)%K_N))
   ENDIF
END FUNCTION phyto_pN
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
FUNCTION findMin(a1,a2,a3,a4) RESULT(theMin)
!-------------------------------------------------------------------------------
!ARGUMENTS
   real(rk),INTENT(in) :: a1,a2,a3,a4
!LOCALS
   real(rk)     :: theMin
!
!-------------------------------------------------------------------------------
!BEGIN
   theMin = a1
   IF(a2 < theMin)      theMin = a2
   IF(a3 < theMin)      theMin = a3
   IF(a4 < theMin)      theMin = a4

   IF( theMin<0.0_rk )  theMin=0.0_rk

END FUNCTION findMin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_respiration(self,group,temp) RESULT(respiration)
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   INTEGER,INTENT(in)                          :: group
   real(rk),INTENT(in)                         :: temp
!
!LOCALS
   real(rk) :: respiration ! Returns the phytoplankton respiration.
!
!-------------------------------------------------------------------------------
!BEGIN
   respiration = self%phytos(group)%R_resp * self%phytos(group)%theta_resp**(temp-20.0)

END FUNCTION phyto_respiration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_salinity(self,group,salinity) RESULT(fSal)
!-------------------------------------------------------------------------------
! Salinity tolerance of phytoplankton
! Implmentation based on Griffin et al 2001; Robson and Hamilton, 2004
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   INTEGER,INTENT(in)                          :: group
   real(rk),INTENT(in)                         :: salinity
!
!LOCALS
   real(rk) :: fSal ! Returns the salinity function
   real(rk) :: tmp1,tmp2,tmp3
   real(rk),PARAMETER :: wq_one = 1.0
!
!-------------------------------------------------------------------------------
!BEGIN

   IF (self%phytos(group)%salTol == 0) THEN
      fSal = 1.0
   ELSEIF (self%phytos(group)%salTol == 1) THEN
      !# f(S) = 1 at S=S_opt, f(S) = S_bep at S=S_maxsp.
      tmp1 = (self%phytos(group)%S_bep-1.0) / ((self%phytos(group)%S_maxsp - self%phytos(group)%S_opt)**2.0)
      tmp2 = (self%phytos(group)%S_bep-1.0) * 2.0*self%phytos(group)%S_opt / &
            ((self%phytos(group)%S_maxsp-self%phytos(group)%S_opt)**2.0)
      tmp3 = (self%phytos(group)%S_bep-1.0) * self%phytos(group)%S_opt*self%phytos(group)%S_opt / &
            ((self%phytos(group)%S_maxsp-self%phytos(group)%S_opt)**2.0) + 1.0
      IF (salinity>self%phytos(group)%S_opt) THEN
         fSal = tmp1*(salinity**2.0)-tmp2*salinity+tmp3
      ELSE
         fSal = 1.0
      ENDIF
   ELSEIF (self%phytos(group)%salTol == 2) THEN
      !# f(S) = 1 at S=S_opt, f(S) = S_bep at S=0.
      IF (salinity<self%phytos(group)%S_opt) THEN
         fSal = (self%phytos(group)%S_bep-1.0) * (salinity**2.0)/(self%phytos(group)%S_opt**2.0) -  &
                      2.0*(self%phytos(group)%S_bep-1.0)*salinity/self%phytos(group)%S_opt+self%phytos(group)%S_bep
      ELSE
        fSal = 1.0
      ENDIF
   ELSEIF (self%phytos(group)%salTol == 3) THEN
      ! f(S) = 1 at S=S_opt, f(S) = S_bep at S=0 and 2*S_opt.
      IF (salinity < self%phytos(group)%S_opt) THEN
      fSal = (self%phytos(group)%S_bep-1.0)*(salinity**2.0)/(self%phytos(group)%S_opt**2.0)-  &
                      2.0*(self%phytos(group)%S_bep-1.0)*salinity/self%phytos(group)%S_opt+self%phytos(group)%S_bep
      ENDIF
      IF ((salinity>self%phytos(group)%S_maxsp) .AND. (salinity<(self%phytos(group)%S_maxsp + self%phytos(group)%S_opt))) THEN
         fSal = (self%phytos(group)%S_bep - wq_one)*(self%phytos(group)%S_maxsp + self%phytos(group)%S_opt - salinity)**2  &
             / (self%phytos(group)%S_opt**2) -                                                                             &
             2 * (self%phytos(group)%S_bep - wq_one) * (self%phytos(group)%S_maxsp + self%phytos(group)%S_opt - salinity)  &
             / self%phytos(group)%S_opt + self%phytos(group)%S_bep
      ENDIF
      IF ( (salinity >= self%phytos(group)%S_opt) .AND. (salinity <= self%phytos(group)%S_maxsp) ) fSal = 1
      IF ( salinity >= (self%phytos(group)%S_maxsp + self%phytos(group)%S_opt) ) fSal = self%phytos(group)%S_bep
   ELSE
      PRINT *,'STOP: Unsupported salTol flag for group: ',group,'=', self%phytos(group)%salTol
   ENDIF

   IF( fSal < 0.0_rk ) fSal = 0.0_rk

END FUNCTION phyto_salinity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION phyto_light(self, group, par, extc, Io, dz) RESULT(fI)
!-------------------------------------------------------------------------------
! Light limitation of pytoplankton via various model approaches. Refer to
! overview presented in Table 1 of:
!
! Baklouti, M., Diaz, F., Pinazo, C., Faure, V., Quéguiner, B., 2006.
!  Investigation of mechanistic formulations depicting phytoplankton dynamics for
!    models of marine pelagic ecosystems and description of a new model.
!  Progress in Oceanography 71 (1), 1-33.
!
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_phytoplankton),INTENT(in) :: self
   INTEGER,INTENT(in)                          :: group
   real(rk),INTENT(in)                         :: par
   real(rk),INTENT(in)                         :: extc
   real(rk),INTENT(in)                         :: Io
   real(rk),INTENT(in)                         :: dz
!
!CONSTANTS
   real(rk),PARAMETER :: one_e_neg3 = 1e-3
!
!LOCALS
   real(rk) :: fI !-- Returns the light limitation
   real(rk) :: par_t,par_b,par_c
   real(rk) :: z1,z2
   real(rk) :: x
   real(rk), PARAMETER :: A = 5.0, eps = 0.5
!
!-------------------------------------------------------------------------------
!BEGIN
   fI    = 0.0

   ! MH fix this
   par_t = par
   par_b = par_t * EXP( -extc * dz )
   par_c = par_t * EXP( -extc * dz/2. )

   SELECT CASE (self%phytos(group)%lightModel)
      CASE ( 0 )
         ! Light limitation without photoinhibition.
         ! This is the Webb et al (1974) model solved using the numerical
         ! integration approach as in CAEDYM (Hipsey and Hamilton, 2008)

         IF (Io == 0.0_rk) RETURN

         z1 = -par_t / self%phytos(group)%I_K
         z2 = -par_b / self%phytos(group)%I_K

         z1 = exp_integral(z1)
         z2 = exp_integral(z2)

         fI = 1.0 + (z2 - z1) / MAX(extc * dz,one_e_neg3)

         ! A simple check
         IF (par_t < 5e-5 .OR. fI < 5e-5) fI = 0.0

      CASE ( 1 )
         ! Light limitation without photoinhibition.
         ! This is the Monod (1950) model.

         x = par_c/self%phytos(group)%I_K
         fI = x / (1.0_rk + x)

      CASE ( 2 )
         ! Light limitation with photoinhibition.
         ! This is the Steele (1962) model.

         x = par_c/self%phytos(group)%I_S
         fI = x * EXP(1.0_rk - x)
         IF (par_t < 5e-5 .OR. fI < 5e-5) fI = 0.0

      CASE ( 3 )
         ! Light limitation without photoinhibition.
         ! This is the Webb et al. (1974) model.

         x = par_c/self%phytos(group)%I_K
         fI = 1.0_rk - EXP(-x)

      CASE ( 4 )
         ! Light limitation without photoinhibition.
         ! This is the Jassby and Platt (1976) model.

         x = par_c/self%phytos(group)%I_K
         fI = TANH(x)

      CASE ( 5 )
         ! Light limitation without photoinhibition.
         ! This is the Chalker (1980) model.

         x = par_c/self%phytos(group)%I_K
         fI = (EXP(x * (1.0_rk + eps)) - 1.0_rk) / &
              (EXP(x * (1.0_rk + eps)) + eps)

      CASE ( 6 )
         ! Light limitation with photoinhibition.
         ! This is the Klepper et al. (1988) / Ebenhoh et al. (1997) model.
         x = par_c/self%phytos(group)%I_S
         fI = ((2.0 + A) * x) / ( 1.0_rk + (A * x) + (x * x) )

      CASE ( 7 )
         ! Light limitation with photoinhibition.
         ! This is an integrated form of Steele model.

         fI = ( EXP(1-par_b/self%phytos(group)%I_S) - &
                EXP(1-par_t/self%phytos(group)%I_S)   ) / (extc * dz)
   END SELECT

   IF ( fI < 0.0_rk ) fI = 0.0_rk
END FUNCTION phyto_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




END MODULE aed_phytoplankton

#endif

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

#include "aed.h"

#define _MOB_CONST_ 0
#define _MOB_TEMP_  1
#define _MOB_CALC_  2


MODULE aed_phytoplankton
!-------------------------------------------------------------------------------
!  aed_phytoplankton --- phytoplankton biogeochemical model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util,ONLY : find_free_lun, &
                       exp_integral,  &
                       aed_bio_temp_function, &
                       fTemp_function
   USE aed_phyto_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_type_phytoplankton
!

   TYPE,extends(type_base_model) :: aed_type_phytoplankton
!     Variable identifiers
      TYPE (type_state_variable_id),ALLOCATABLE :: id_p(:)
      TYPE (type_state_variable_id),ALLOCATABLE :: id_in(:)
      TYPE (type_state_variable_id),ALLOCATABLE :: id_ip(:)
      TYPE (type_state_variable_id),ALLOCATABLE :: id_rho(:)
      TYPE (type_state_variable_id)        :: id_Pexctarget,id_Pmorttarget,id_Pupttarget(1:2)
      TYPE (type_state_variable_id)        :: id_Nexctarget,id_Nmorttarget,id_Nupttarget(1:4)
      TYPE (type_state_variable_id)        :: id_Cexctarget,id_Cmorttarget,id_Cupttarget
      TYPE (type_state_variable_id)        :: id_Siexctarget,id_Simorttarget,id_Siupttarget
      TYPE (type_state_variable_id)        :: id_DOupttarget
      TYPE (type_dependency_id)            :: id_par, id_tem, id_sal, id_dz, id_extc
      TYPE (type_horizontal_dependency_id) :: id_I_0
      TYPE (type_diagnostic_variable_id)   :: id_GPP, id_NCP, id_PPR, id_NPR, id_dPAR
      TYPE (type_diagnostic_variable_id)   :: id_TPHY, id_TCHLA, id_TIN, id_TIP
      TYPE (type_diagnostic_variable_id)   :: id_NUP, id_PUP, id_CUP
      TYPE (type_diagnostic_variable_id),ALLOCATABLE :: id_NtoP(:)
      TYPE (type_diagnostic_variable_id),ALLOCATABLE :: id_fT(:), id_fI(:), id_fNit(:), id_fPho(:), id_fSil(:), id_fSal(:)

!     Model parameters
      INTEGER                                   :: num_phytos
      TYPE(phyto_data),DIMENSION(:),ALLOCATABLE :: phytos
      ! LOGICAL                                 :: do_exc,do_mort,do_upt, do_N2uptake
      LOGICAL                                   :: do_Puptake, do_Nuptake, do_Cuptake
      LOGICAL                                   :: do_Siuptake, do_DOuptake, do_N2uptake
      LOGICAL                                   :: do_Pmort, do_Nmort, do_Cmort, do_Simort
      LOGICAL                                   :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER                                   :: nnup, npup
      AED_REAL                                  :: dic_per_n

      CONTAINS     ! Model Methods
        PROCEDURE :: initialize               => aed_init_phytoplankton
        PROCEDURE :: do                       => aed_phytoplankton_do
        PROCEDURE :: do_ppdd                  => aed_phytoplankton_do_ppdd
        PROCEDURE :: do_benthos               => aed_phytoplankton_do_benthos
        PROCEDURE :: get_vertical_movement    => aed_phytoplankton_get_vertical_movement
        PROCEDURE :: get_light_extinction     => aed_phytoplankton_get_light_extinction
   END TYPE

   AED_REAL :: dtlim = 0.9 * 3600
   LOGICAL  :: extra_debug = .false.

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_phytoplankton_load_params(self, dbase, count, list, w_model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phytoplankton),INTENT(inout) :: self
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
   INTEGER,INTENT(in)          :: w_model(*)
!
!LOCALS
   INTEGER  :: status
!  AED_REAL :: p_initial=0.
!  AED_REAL :: p0=0.0225
!  AED_REAL :: w_p=-1.157407e-05
!  AED_REAL :: i_min=25.
!  AED_REAL :: rmax=1.157407e-05
!  AED_REAL :: alpha=0.3
!  AED_REAL :: rpn=1.157407e-07
!  AED_REAL :: rpdu=2.314814e-07
!  AED_REAL :: rpdl=1.157407e-06

   INTEGER  :: i,tfil
   AED_REAL :: minNut

   AED_REAL, parameter :: secs_pr_day = 86400.
   TYPE(phyto_nml_data) :: pd(MAX_PHYTO_TYPES)
   NAMELIST /phyto_data/ pd
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD', iostat=status)
    IF (status /= 0) STOP 'Cannot open phyto_data namelist file'
    read(tfil,nml=phyto_data,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist phyto_data'

    self%num_phytos = count
    ALLOCATE(self%phytos(count))
    ALLOCATE(self%id_p(count))
    ALLOCATE(self%id_in(count))
    ALLOCATE(self%id_ip(count))
    ALLOCATE(self%id_rho(count))
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
       self%phytos(i)%w_model      = w_model(i)
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
                              'mmol/m**3',                                     &
                              'phytoplankton '//TRIM(self%phytos(i)%p_name),   &
                              pd(list(i))%p_initial,                           &
                              minimum=pd(list(i))%p0,                          &
                              vertical_movement = self%phytos(i)%w_p)

       IF (self%phytos(i)%w_model == _MOB_CALC_) THEN
          ! Register rho group as a state variable
          CALL self%register_state_variable(self%id_rho(i),                    &
                              TRIM(self%phytos(i)%p_name)//'_rho',             &
                              'mmol/m**3',                                     &
                        'phytoplankton '//TRIM(self%phytos(i)%p_name)//'_rho', &
                              pd(list(i))%w_p,                                 &
                              minimum=zero_,                                   &
                              vertical_movement = self%phytos(i)%w_p)

       ENDIF
       IF (self%phytos(i)%simINDynamics /= 0) THEN
          IF(self%phytos(i)%simINDynamics == 1)THEN
            minNut = self%phytos(i)%p0*self%phytos(i)%X_ncon
          ELSE
            minNut = self%phytos(i)%p0*self%phytos(i)%X_nmin
          ENDIF
          ! Register IN group as a state variable
          CALL self%register_state_variable(self%id_in(i),                     &
                              TRIM(self%phytos(i)%p_name)//'_IN',              &
                              'mmol/m**3',                                     &
                         'phytoplankton '//TRIM(self%phytos(i)%p_name)//'_IN', &
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
                              'mmol/m**3',                                     &
                         'phytoplankton '//TRIM(self%phytos(i)%p_name)//'_IP', &
                              pd(list(i))%p_initial*self%phytos(i)%X_pcon,     &
                              minimum=minNut,                                  &
                              vertical_movement = self%phytos(i)%w_p)

       ENDIF

       CALL self%register_diagnostic_variable(self%id_NtoP(i), TRIM(self%phytos(i)%p_name)//'_NtoP','mmol/m**3', 'INi/IPi')

       IF (extra_debug) THEN
          CALL self%register_diagnostic_variable(self%id_fT(i), TRIM(self%phytos(i)%p_name)//'_fT','mmol/m**3', 'fT')
          CALL self%register_diagnostic_variable(self%id_fI(i), TRIM(self%phytos(i)%p_name)//'_fI','mmol/m**3', 'fI')
          CALL self%register_diagnostic_variable(self%id_fNit(i), TRIM(self%phytos(i)%p_name)//'_fNit','mmol/m**3', 'fNit')
          CALL self%register_diagnostic_variable(self%id_fPho(i), TRIM(self%phytos(i)%p_name)//'_fPho','mmol/m**3', 'fPho')
          CALL self%register_diagnostic_variable(self%id_fSil(i), TRIM(self%phytos(i)%p_name)//'_fSil','mmol/m**3', 'fSil')
          CALL self%register_diagnostic_variable(self%id_fSal(i), TRIM(self%phytos(i)%p_name)//'_fSal','mmol/m**3', 'fSal')
       ENDIF
    ENDDO
END SUBROUTINE aed_phytoplankton_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_init_phytoplankton(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the phytoplankton biogeochemical model
!
!  Here, the aed_p_m namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phytoplankton),TARGET,INTENT(inout) :: self
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status

   INTEGER            :: num_phytos
   INTEGER            :: the_phytos(MAX_PHYTO_TYPES)
   INTEGER            :: w_model(MAX_PHYTO_TYPES)
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
   CHARACTER(len=128) :: dbase='aed_phyto_pars.nml'


   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   AED_REAL           :: zerolimitfudgefactor = 0.9 * 3600
   NAMELIST /aed_phytoplankton/ num_phytos, the_phytos, w_model,               &
                    p_excretion_target_variable,p_mortality_target_variable,   &
                     p1_uptake_target_variable, p2_uptake_target_variable,     &
                    n_excretion_target_variable,n_mortality_target_variable,   &
                     n1_uptake_target_variable,n2_uptake_target_variable,      &
                     n3_uptake_target_variable,n4_uptake_target_variable,      &
                    c_excretion_target_variable,c_mortality_target_variable,   &
                      c_uptake_target_variable, do_uptake_target_variable,     &
                    si_excretion_target_variable,si_mortality_target_variable, &
                      si_uptake_target_variable,                               &
                    dbase, zerolimitfudgefactor, extra_debug
!-----------------------------------------------------------------------
!BEGIN
   w_model = _MOB_CONST_
   ! Read the namelist
   read(namlst,nml=aed_phytoplankton,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_phytoplankton'
   dtlim = zerolimitfudgefactor

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   CALL aed_phytoplankton_load_params(self, dbase, num_phytos, the_phytos, w_model)

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
   CALL self%register_diagnostic_variable(self%id_GPP, 'GPP','mmol/m**3',  'gross primary production')
   CALL self%register_diagnostic_variable(self%id_NCP, 'NCP','mmol/m**3',  'net community production')
   CALL self%register_diagnostic_variable(self%id_PPR, 'PPR','mmol/m**3/d','gross primary production rate')
   CALL self%register_diagnostic_variable(self%id_NPR, 'NPR','mmol/m**3/d','net community production rate')

   CALL self%register_diagnostic_variable(self%id_NUP, 'NUP','mmol/m**3/d','nitrogen uptake')
   CALL self%register_diagnostic_variable(self%id_PUP, 'PUP','mmol/m**3/d','phosphorous uptake')
   CALL self%register_diagnostic_variable(self%id_CUP, 'CUP','mmol/m**3/d','carbon uptake')

   CALL self%register_diagnostic_variable(self%id_dPAR, 'PAR','W/m**2',  'photosynthetically active radiation')
   CALL self%register_diagnostic_variable(self%id_TCHLA, 'TCHLA','ug/L', 'Total Chlorophyll-a')
   CALL self%register_diagnostic_variable(self%id_TPHY, 'TPHYS','ug/L',  'Total Phytoplankton')
   CALL self%register_diagnostic_variable(self%id_TIN, 'IN','ug/L',      'Total Chlorophyll-a')
   CALL self%register_diagnostic_variable(self%id_TIP, 'IP','ug/L',      'Total Chlorophyll-a')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_tem,  standard_variables%temperature)
   CALL self%register_dependency(self%id_sal,  standard_variables%practical_salinity)
   CALL self%register_dependency(self%id_par,  standard_variables%downwelling_photosynthetic_radiative_flux)
   CALL self%register_horizontal_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
   CALL self%register_dependency(self%id_dz,   standard_variables%cell_thickness)
   CALL self%register_dependency(self%id_extc, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
END SUBROUTINE aed_init_phytoplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of phytoplankton biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phytoplankton),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   AED_REAL :: phy, tphy, tin, tip, tchla
   AED_REAL :: INi, IPi
   AED_REAL :: pup
   AED_REAL :: no3up,nh4up
   AED_REAL :: cup, rsiup
   AED_REAL :: temp, par, Io, salinity, extc, dz
   AED_REAL :: primprod(self%num_phytos), exudation(self%num_phytos), &
               a_nfix(self%num_phytos), respiration(self%num_phytos)
   AED_REAL :: cuptake(self%num_phytos), cexcretion(self%num_phytos), cmortality(self%num_phytos)
   AED_REAL :: nuptake(self%num_phytos,1:4), nexcretion(self%num_phytos), nmortality(self%num_phytos)
   AED_REAL :: puptake(self%num_phytos,1:2), pexcretion(self%num_phytos), pmortality(self%num_phytos)
   AED_REAL :: siuptake(self%num_phytos), siexcretion(self%num_phytos), simortality(self%num_phytos)
   AED_REAL :: fT, fNit, fPho, fSil, fI, fXl, fSal, PNf
   AED_REAL :: upTot

   INTEGER  :: phy_i,c
   AED_REAL :: flux, available

! MH to fix
!  AED_REAL :: dt = 3600. ! just for now, hard code it
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

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

      primprod(phy_i)    = zero_
      exudation(phy_i)   = zero_
      a_nfix(phy_i)      = zero_
      respiration(phy_i) = zero_

      cuptake(phy_i)     = zero_
      cexcretion(phy_i)  = zero_
      cmortality(phy_i)  = zero_
      nuptake(phy_i,:)   = zero_
      nexcretion(phy_i)  = zero_
      nmortality(phy_i)  = zero_
      puptake(phy_i,:)   = zero_
      pexcretion(phy_i)  = zero_
      pmortality(phy_i)  = zero_

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
            fNit = phyto_fN(self%phytos,phy_i,IN=fNit)
         ENDIF
         IF (phy > zero_ .AND. phy <= self%phytos(phy_i)%p0) THEN
            fNit = phyto_fN(self%phytos,phy_i,din=no3up+nh4up)
         ENDIF
      ELSE
         fNit = phyto_fN(self%phytos,phy_i,din=no3up+nh4up)
      ENDIF
      IF (self%phytos(phy_i)%simNFixation /= 0) THEN
         ! Nitrogen fixer: apply no N limitation. N Fixation ability
         ! depends on DIN concentration
         a_nfix = (one_ - fNit)
         fNit = one_
      ENDIF


      ! PHOSPHOROUS.
      fPho = zero_
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
            fPho = phyto_fP(self%phytos,phy_i,IP=fPho)
         ENDIF
         IF (phy > zero_ .AND. phy <= self%phytos(phy_i)%p0) THEN
            fPho = phyto_fP(self%phytos,phy_i,frp=pup)
         ENDIF
      ELSE
         fPho = phyto_fP(self%phytos,phy_i,frp=pup)
      ENDIF

      ! SILICA.
      fSil = phyto_fSi(self%phytos,phy_i,rsiup)


      ! LIGHT
      _GET_(self%id_extc,extc)
      ! dz = 0.5     !MH: to fix
      _GET_(self%id_dz,dz)
      fI = phyto_light(self%phytos, phy_i, par, extc, Io, dz)
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

      respiration(phy_i) = phyto_respiration(self%phytos,phy_i,temp)

      ! Salinity stress effect on respiration
      fSal =  phyto_salinity(self%phytos,phy_i,salinity)
      respiration(phy_i) = respiration(phy_i) * fSal

      ! photo-exudation
      exudation(phy_i) = primprod(phy_i)*self%phytos(phy_i)%f_pr

      ! Limit respiration if at the min biomass to prevent
      ! leak in the C mass balance
      IF (phy <= self%phytos(phy_i)%p0) THEN
         respiration(phy_i) = zero_
         exudation(phy_i) = zero_
      ENDIF

      ! write(*,"(4X,'limitations (fT,fI,fN,fP,fSi,Io, par, mu): ',9F9.2)")fT,fI,fNit,fPho,fSil,Io,par,primprod*secs_pr_day


      ! Carbon uptake and excretion

      cuptake(phy_i)    = -primprod(phy_i) * phy
      cexcretion(phy_i) = (self%phytos(phy_i)%k_fdom*(1.0-self%phytos(phy_i)%k_fres)*respiration(phy_i)+exudation(phy_i)) * phy
      cmortality(phy_i) = ((1.0-self%phytos(phy_i)%k_fdom)*(1.0-self%phytos(phy_i)%k_fres)*respiration(phy_i)) * phy

      ! Nitrogen uptake and excretion

      CALL phyto_internal_nitrogen(self%phytos,phy_i,self%do_N2uptake,phy,INi,primprod(phy_i),&
                             fT,no3up,nh4up,a_nfix(phy_i),respiration(phy_i),exudation(phy_i),PNf,&
                                   nuptake(phy_i,:),nexcretion(phy_i),nmortality(phy_i))

      ! Phosphorus uptake and excretion

      CALL phyto_internal_phosphorus(self%phytos,phy_i,self%npup,phy,IPi,primprod(phy_i),&
                                 fT,pup,respiration(phy_i),exudation(phy_i),&
                                         puptake(phy_i,:),pexcretion(phy_i),pmortality(phy_i))

      ! Silica uptake and excretion

      IF (self%phytos(phy_i)%simSiUptake > 0) THEN
         siuptake(phy_i)    =-self%phytos(phy_i)%X_sicon * primprod(phy_i) * phy
         siexcretion(phy_i) = self%phytos(phy_i)%X_sicon * (self%phytos(phy_i)%k_fdom*respiration(phy_i)+exudation(phy_i)) * phy
         simortality(phy_i) = self%phytos(phy_i)%X_sicon * ((1.0-self%phytos(phy_i)%k_fdom)*respiration(phy_i)) * phy
      ELSE
         siuptake(phy_i)    = zero_
         siexcretion(phy_i) = zero_
         simortality(phy_i) = zero_
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
      IF ( abs(upTot) > 1.e-10 .and. upTot >= pup ) THEN
         DO phy_i=1,self%num_phytos
            puptake(phy_i,1) = (pup*0.99/dtlim) * (puptake(phy_i,1)/upTot)
         ENDDO
      ENDIF
   ENDIF

   IF (self%do_Nuptake) THEN
      upTot = sum(nuptake(:,1))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= no3up ) THEN
         DO phy_i=1,self%num_phytos
            nuptake(phy_i,1) = (no3up*0.99/dtlim) * (nuptake(phy_i,1)/upTot)
         ENDDO
      ENDIF

      upTot = sum(nuptake(:,2))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= nh4up ) THEN
         DO phy_i=1,self%num_phytos
            nuptake(phy_i,2) = (nh4up*0.99/dtlim) * (nuptake(phy_i,2)/upTot)
         ENDDO
      ENDIF
   ENDIF
   IF (self%do_Cuptake) THEN
      upTot = sum(cuptake)*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= cup ) THEN
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
      IF ( abs(upTot) > 1.e-10 .and. upTot >= rsiup ) THEN
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
      available = MAX(zero_, phy - self%phytos(phy_i)%p0)
      IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
      _SET_ODE_(self%id_p(phy_i), flux)

      IF (self%phytos(phy_i)%simINDynamics /= 0) THEN
         ! _SET_ODE_(self%id_in(phy_i), (-sum(nuptake) - nexcretion(phy_i) - nmortality(phy_i) )*INi )
         _GET_(self%id_in(phy_i),INi)
         flux = (-sum(nuptake(phy_i,:)) - nexcretion(phy_i) - nmortality(phy_i) )
         available = MAX(zero_, INi - self%phytos(phy_i)%X_nmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _SET_ODE_(self%id_in(phy_i), flux)
      ENDIF
      IF (self%phytos(phy_i)%simIPDynamics /= 0) THEN
         ! _SET_ODE_(self%id_ip(phy_i), (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) ) )
         _GET_(self%id_ip(phy_i),IPi)
         flux = (-sum(puptake(phy_i,:)) - pexcretion(phy_i) - pmortality(phy_i) )
         available = MAX(zero_, IPi - self%phytos(phy_i)%X_pmin*phy)
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
   _LOOP_END_

END SUBROUTINE aed_phytoplankton_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic sedimentation of phytoplankton.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phytoplankton),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   AED_REAL :: phy        ! State
   INTEGER  :: phy_i
   AED_REAL :: phy_flux

   ! Parameters
!  AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN

   _HORIZONTAL_LOOP_BEGIN_

   DO phy_i=1,self%num_phytos
      ! Retrieve current (local) state variable values.
      _GET_(self%id_p(phy_i),phy) ! phytoplankton

      phy_flux = zero_  !self%phytos(phy_i)%w_p*MAX(phy,zero_)

     ! Set bottom fluxes for the pelagic (change per surface area per second)
     ! Transfer sediment flux value to FABM.
     _SET_BOTTOM_EXCHANGE_(self%id_p(phy_i),phy_flux)

   ENDDO

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

END SUBROUTINE aed_phytoplankton_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
!-------------------------------------------------------------------------------
! Get the vertical movement values
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phytoplankton),INTENT(in) :: self
   _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
!
!LOCALS
   AED_REAL :: temp
   AED_REAL :: w_rho
   AED_REAL :: phy
   INTEGER  :: phy_i
!
!-------------------------------------------------------------------------------
!BEGIN
   _LOOP_BEGIN_

   _GET_DEPENDENCY_(self%id_tem, temp)

   DO phy_i=1,self%num_phytos
      SELECT CASE (self%phytos(phy_i)%w_model)
         CASE (_MOB_CONST_)
            w_rho = self%phytos(phy_i)%w_p
         CASE (_MOB_TEMP_)
            w_rho = self%phytos(phy_i)%w_p * temp !# MH to fix
         CASE (_MOB_CALC_)
            !# MH to complete
            _GET_(self%id_p(phy_i), phy) ! phytoplankton
             w_rho = self%phytos(phy_i)%w_p
            _GET_(self%id_rho(phy_i), w_rho)
         CASE DEFAULT
             STOP
      END SELECT

      _SET_VERTICAL_MOVEMENT_(self%id_p(phy_i), w_rho)
   ENDDO

   _LOOP_END_
END SUBROUTINE aed_phytoplankton_get_vertical_movement
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phytoplankton),INTENT(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
!LOCALS
   AED_REAL :: phy
   INTEGER  :: phy_i
!
!-----------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   DO phy_i=1,self%num_phytos
      ! Retrieve current (local) state variable values.
      _GET_(self%id_p(phy_i),phy) ! phytoplankton

      ! Self-shading with explicit contribution from background phytoplankton concentration.
      _SET_EXTINCTION_(self%phytos(phy_i)%KePHY*phy)

   ENDDO

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_phytoplankton_get_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phytoplankton_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of phytoplankton biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phytoplankton),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL            :: n,p,par,I_0
   INTEGER             :: phy_i
   AED_REAL            :: iopt,rpd,primprod
   AED_REAL,PARAMETER  :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

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
   _LOOP_END_

END SUBROUTINE aed_phytoplankton_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_phytoplankton

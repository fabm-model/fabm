!###############################################################################
!#                                                                             #
!# aed_pathogens.F90                                                           #
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
!# Created July 2012                                                           #
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_pathogens
!-------------------------------------------------------------------------------
!  aed_pathogens --- pathogen biogeochemical model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util,ONLY : find_free_lun

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_type_pathogens
!



   TYPE pathogen_nml_data
      CHARACTER(64) :: p_name
      AED_REAL      :: coef_grwth_uMAX                     !-- Max growth rate at 20C
      AED_REAL      :: coef_grwth_Tmin, coef_grwth_Tmax    !-- Tmin and Tmax, f(T)
      AED_REAL      :: coef_grwth_T1, coef_grwth_T2        !-- coef_grwth_T1  and  coef_grwth_T2
      AED_REAL      :: coef_grwth_Kdoc                     !-- Half-saturation for growth, coef_grwth_Kdoc
      AED_REAL      :: coef_grwth_ic                       !-- coef_grwth_ic
      AED_REAL      :: coef_mort_kd20                      !-- Mortality rate (Dark death rate) @ 20C and 0 psu
      AED_REAL      :: coef_mort_theta                     !-- Temperature multiplier for mortality: coef_mort_theta
      AED_REAL      :: coef_mort_c_SM, coef_mort_alpha, coef_mort_beta  !-- Salinity effect on mortality
      AED_REAL      :: coef_mort_c_PHM, coef_mort_K_PHM, coef_mort_delta_M  !-- pH effect on mortality
      AED_REAL      :: coef_mort_fdoc                      !-- Fraction of mortality back to doc
      AED_REAL      :: coef_light_kb_vis, coef_light_kb_uva, coef_light_kb_uvb !-- Light inactivation
      AED_REAL      :: coef_light_cSb_vis, coef_light_cSb_uva, coef_light_cSb_uvb !-- Salinity effect on light inactivation
      AED_REAL      :: coef_light_kDOb_vis, coef_light_kDOb_uva, coef_light_kDOb_uvb !-- DO effect on light
      AED_REAL      :: coef_light_cpHb_vis, coef_light_cpHb_uva, coef_light_cpHb_uvb !-- pH effect on light inactivation
      AED_REAL      :: coef_light_KpHb_vis, coef_light_KpHb_uva, coef_light_KpHb_uvb !-- pH effect on light inactivation
      AED_REAL      :: coef_light_delb_vis, coef_light_delb_uva, coef_light_delb_uvb !-- exponent for pH effect on light inactivation
      AED_REAL      :: coef_pred_kp20, coef_pred_theta_P   !-- Loss rate due to predation and temp multiplier
      AED_REAL      :: coef_sett_fa                        !-- Attached fraction in water column
      AED_REAL      :: coef_sett_w_path      !-- Sedimentation velocity (m/d) at 20C (-ve means down) for NON-ATTACHED orgs
   END TYPE

!  TYPE pathogen_data
!     ! General Attributes
!     TYPE(pathogen_nml_data) :: par
!  END TYPE

   TYPE,extends(type_base_model) :: aed_type_pathogens
!     Variable identifiers
      TYPE (type_state_variable_id),ALLOCATABLE :: id_p(:)
      TYPE (type_state_variable_id)      :: id_growth, id_mortality, id_sunlight, id_grazing
      TYPE (type_dependency_id)          :: id_par, id_tem, id_sal
      TYPE (type_dependency_id)          :: id_oxy, id_pH,  id_doc, id_tss
      TYPE (type_horizontal_dependency_id)  :: id_I_0

!     Model parameters
      INTEGER                                   :: num_pathogens
      TYPE(pathogen_nml_data),DIMENSION(:),ALLOCATABLE :: pathogens
      LOGICAL                                   :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER :: nnup, npup
      AED_REAL :: dic_per_n

      CONTAINS     ! Model Methods
        PROCEDURE :: initialize               => aed_init_pathogens
        PROCEDURE :: do                       => aed_pathogens_do
        PROCEDURE :: do_ppdd                  => aed_pathogens_do_ppdd
        PROCEDURE :: do_benthos               => aed_pathogens_do_benthos
!       PROCEDURE :: get_light_extinction     => aed_pathogens_get_light_extinction
   END TYPE

   AED_REAL, parameter :: secs_pr_day = 86400.

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_pathogens(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the pathogen biogeochemical model
!
!  Here, the aed_p_m namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_type_pathogens),TARGET,INTENT(inout) :: self

!
!LOCALS
   INTEGER :: status

   INTEGER :: num_pathogens
   INTEGER :: the_pathogens(MAX_PATHO_TYPES)


   NAMELIST /aed_pathogens/ num_pathogens, the_pathogens
!-----------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_pathogens,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_pathogens'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   CALL aed_pathogens_load_params(self, num_pathogens, the_pathogens)

   ! Register state dependancies
!  self%do_Pexc = p_excretion_target_variable .NE. ''
!  IF (self%do_Pexc) THEN
!    self%id_Pexctarget  = self%register_state_dependency(p_excretion_target_variable)
!  ENDIF


   ! Register environmental dependencies
   CALL self%register_dependency(self%id_tem,standard_variables%temperature)
   CALL self%register_dependency(self%id_sal,standard_variables%practical_salinity)
   CALL self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   CALL self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_photosynthetic_radiative_flux)

END SUBROUTINE aed_init_pathogens
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed_pathogens_load_params(self, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_pathogens),INTENT(inout) :: self
   INTEGER,INTENT(in) :: count
   INTEGER,INTENT(in) :: list(*)
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,tfil
   AED_REAL :: minPath

   TYPE(pathogen_nml_data) :: pd(MAX_PATHO_TYPES)
   NAMELIST /pathogen_data/ pd
!-------------------------------------------------------------------------------
!BEGIN
    minPath = 1e-10
    tfil = find_free_lun()
    open(tfil,file="aed_pathogen_pars.nml", status='OLD',iostat=status)
    IF (status /= 0) STOP 'Error opening namelist pathogen_data'
    read(tfil,nml=pathogen_data,iostat=status)
    IF (status /= 0) STOP 'Error reading namelist pathogen_data'
    close(tfil)

    self%num_pathogens = count
    ALLOCATE(self%pathogens(count))
    ALLOCATE(self%id_p(count))
    DO i=1,count
       ! Assign parameters from database to simulated groups
       !self%pathogens(i)%p_name       = pd(list(i))%p_name
       self%pathogens(i)          = pd(list(i))

       ! Register group as a state variable
       CALL self%register_state_variable(self%id_p(i),                        &
                             TRIM(self%pathogens(i)%p_name),                  &
                             'mmol/m**3', 'pathogen',                         &
                             minPath,                                         &
                            ! pd(list(i))%p_initial,                          &
                             minimum=minPath,                                 &
                             !minimum=pd(list(i))%p0,                         &
                             vertical_movement = self%pathogens(i)%coef_sett_w_path)


!      IF (self%pathogens(i)%p_name == 'crypto') THEN
!         ! Register IN group as a state variable
!         self%id_in(i) = self%register_state_variable(                       &
!                             TRIM(self%pathogens(i)%p_name)//'_dd',          &
!                             'mmol/m**3', 'pathogen dd',                     &
!                             0.0,                                            &
!                             0.0,                                            &
!                             vertical_movement = self%pathogens(i)%coef_sett_w_path)

!      ENDIF
    ENDDO
END SUBROUTINE aed_pathogens_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_pathogens_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of pathogen biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_pathogens),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   AED_REAL           :: pth
   AED_REAL           :: temp,par,Io,salinity
   AED_REAL           :: growth,light,mortality, predation
   AED_REAL           :: f_AOC,f_pH,f_DO,phi,lightBW,phstar

   INTEGER            :: pth_i,c

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current environmental conditions.
   _GET_   (self%id_tem,temp)     ! local temperature
   _GET_   (self%id_sal,salinity) ! local salinity
   _GET_   (self%id_par,par)      ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,Io)       ! surface short wave radiation
   !_GET_   (self%id_doc,doc)     ! local DOC
   !_GET_   (self%id_oxy,oxy)     ! local oxygen
   !_GET_   (self%id_ph,ph)       ! local pH
   phstar = 0.0 !abs(ph-7.)


   DO pth_i=1,self%num_pathogens

      ! Retrieve this pathogen group
      _GET_(self%id_p(pth_i),pth)

      growth    = zero_
      predation = zero_

      ! Natural mortality (as impacted by T, S, pH)
      f_AOC = 1.0 ! aoc / (K_AOC + aoc)
      f_pH  = 1.0 ! + c_PH * ( pH_star**delta / (pH_star**delta+K_PH**delta) )
      mortality = self%pathogens(pth_i)%coef_mort_kd20/86400.  &
                + (self%pathogens(pth_i)%coef_mort_c_SM*salinity**self%pathogens(pth_i)%coef_mort_alpha) &
                * ((1.0-f_AOC)**self%pathogens(pth_i)%coef_mort_beta) * f_pH
      mortality = mortality * (self%pathogens(pth_i)%coef_mort_theta**(temp-20.0))


      ! Sunlight inactivation (as impacted by S, DO and pH)
      light     = zero_
      lightBW   = zero_
      phi  = 1e-6  ! Convert J to MJ as kb is in m2/MJ)
      ! Visible
      f_DO = 1.0 !oxy / (coef_light_kDOb_vis + oxy)
      f_pH = 1.0 !(1.0 + coef_light_cpHb_vis*(pH_star**coef_light_delb_vis / (coef_light_KpHb_vis**coef_light_delb_vis+pH_star**coef_light_delb_vis)))
      lightBW = phi * (self%pathogens(pth_i)%coef_light_kb_vis + self%pathogens(pth_i)%coef_light_cSb_vis*salinity)
      lightBW = lightBW * par * f_pH * f_DO
      light     = light + lightBW
      ! UV-A
      f_DO = 1.0 !oxy / (coef_light_kDOb_uva + oxy)
      f_pH = 1.0 !(1.0 + coef_light_cpHb_uva*(pH_star**coef_light_delb_uva / (coef_light_KpHb_uva**coef_light_delb_uva+pH_star**coef_light_delb_uva)))
      lightBW = phi * (self%pathogens(pth_i)%coef_light_kb_uva + self%pathogens(pth_i)%coef_light_cSb_uva*salinity)
      lightBW = lightBW * (par*0.03) * f_pH * f_DO
      light     = light + lightBW
      ! UV-B
      f_DO = 1.0 !oxy / (coef_light_kDOb_uvb + oxy)
      f_pH = 1.0 !(1.0 + coef_light_cpHb_uvb*(pH_star**coef_light_delb_uvb / (coef_light_KpHb_uvb**coef_light_delb_uvb+pH_star**coef_light_delb_uvb)))
      lightBW = phi * (self%pathogens(pth_i)%coef_light_kb_uvb + self%pathogens(pth_i)%coef_light_cSb_uvb*salinity)
      lightBW = lightBW * (par*0.003) * f_pH * f_DO
      light     = light + lightBW

      !-----------------------------------------------------------------
      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Pathogen production / losses
      _SET_ODE_(self%id_p(pth_i), (growth - light - mortality - predation)*pth )

   ENDDO

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_pathogens_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_pathogens_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic sedimentation of pathogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_pathogens),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   AED_REAL :: pth        ! State
   INTEGER  :: pth_i
   AED_REAL :: pth_flux

   ! Parameters
!
!-------------------------------------------------------------------------------
!BEGIN

   _HORIZONTAL_LOOP_BEGIN_

   DO pth_i=1,self%num_pathogens
      ! Retrieve current (local) state variable values.
      _GET_(self%id_p(pth_i),pth) ! pathogen

      pth_flux = zero_  !self%pathogens(pth_i)%w_p*MAX(pth,zero_)

     ! Set bottom fluxes for the pelagic (change per surface area per second)
     ! Transfer sediment flux value to FABM.
     _SET_BOTTOM_EXCHANGE_(self%id_p(pth_i),pth_flux)

   ENDDO

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_
END SUBROUTINE aed_pathogens_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_pathogens_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of pathogen biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_pathogens),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL            :: n,p,par,I_0
   INTEGER             :: pth_i
   AED_REAL            :: iopt,rpd,primprod
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   DO pth_i=1,self%num_pathogens
      ! Retrieve current (local) state variable values.
!     _GET_(self%id_p(pth_i),p) ! pathogen

      ! Retrieve current environmental conditions.
!     _GET_   (self%id_par,par)  ! local photosynthetically active radiation
!     _GET_DEPENDENCY_HZ_(self%id_I_0,I_0)  ! surface short wave radiation

   ENDDO

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_pathogens_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_pathogens

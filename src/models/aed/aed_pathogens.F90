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

#ifdef _FABM_F2003_

#include "aed.h"


MODULE aed_pathogens
!-------------------------------------------------------------------------------
!  aed_pathogens --- pathogen biogeochemical model
!-------------------------------------------------------------------------------
   USE fabm_types
   USE fabm_driver
   USE aed_util,ONLY : find_free_lun

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC type_aed_pathogens, aed_pathogens_create
!



   TYPE pathogen_nml_data
      CHARACTER(64) :: p_name
      real(rk)      :: coef_grwth_uMAX                     !-- Max growth rate at 20C
      real(rk)      :: coef_grwth_Tmin, coef_grwth_Tmax    !-- Tmin and Tmax, f(T)
      real(rk)      :: coef_grwth_T1, coef_grwth_T2        !-- coef_grwth_T1  and  coef_grwth_T2
      real(rk)      :: coef_grwth_Kdoc                     !-- Half-saturation for growth, coef_grwth_Kdoc
      real(rk)      :: coef_grwth_ic                       !-- coef_grwth_ic
      real(rk)      :: coef_mort_kd20                      !-- Mortality rate (Dark death rate) @ 20C and 0 psu
      real(rk)      :: coef_mort_theta                     !-- Temperature multiplier for mortality: coef_mort_theta
      real(rk)      :: coef_mort_c_SM, coef_mort_alpha, coef_mort_beta  !-- Salinity effect on mortality
      real(rk)      :: coef_mort_c_PHM, coef_mort_K_PHM, coef_mort_delta_M  !-- pH effect on mortality
      real(rk)      :: coef_mort_fdoc                      !-- Fraction of mortality back to doc
      real(rk)      :: coef_light_kb_vis, coef_light_kb_uva, coef_light_kb_uvb !-- Light inactivation
      real(rk)      :: coef_light_cSb_vis, coef_light_cSb_uva, coef_light_cSb_uvb !-- Salinity effect on light inactivation
      real(rk)      :: coef_light_kDOb_vis, coef_light_kDOb_uva, coef_light_kDOb_uvb !-- DO effect on light
      real(rk)      :: coef_light_cpHb_vis, coef_light_cpHb_uva, coef_light_cpHb_uvb !-- pH effect on light inactivation
      real(rk)      :: coef_light_KpHb_vis, coef_light_KpHb_uva, coef_light_KpHb_uvb !-- pH effect on light inactivation
      real(rk)      :: coef_light_delb_vis, coef_light_delb_uva, coef_light_delb_uvb !-- exponent for pH effect on light inactivation
      real(rk)      :: coef_pred_kp20, coef_pred_theta_P   !-- Loss rate due to predation and temp multiplier
      real(rk)      :: coef_sett_fa                        !-- Attached fraction in water column
      real(rk)      :: coef_sett_w_path      !-- Sedimentation velocity (m/d) at 20C (-ve means down) for NON-ATTACHED orgs
   END TYPE

   TYPE pathogen_data
      ! General Attributes
      TYPE(pathogen_nml_data) :: par
   END TYPE

   TYPE,extends(type_base_model) :: type_aed_pathogens
!     Variable identifiers
      type (type_state_variable_id),ALLOCATABLE :: id_p(:)
      type (type_state_variable_id)      :: id_growth, id_mortality, id_sunlight, id_grazing
      type (type_dependency_id)          :: id_par, id_tem, id_sal
      type (type_dependency_id)          :: id_oxy, id_pH,  id_doc, id_tss
      type (type_horizontal_dependency_id)  :: id_I_0
!     type (type_diagnostic_variable_id) :: ??
!     type (type_conserved_quantity_id)  :: ??

!     Model parameters
      INTEGER                                   :: num_pathogens
      TYPE(pathogen_nml_data),DIMENSION(:),ALLOCATABLE :: pathogens
      LOGICAL                                   :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER :: nnup, npup
      real(rk) :: dic_per_n

      CONTAINS     ! Model Methods
!       procedure :: initialize               => aed_pathogens_init
        procedure :: do                       => aed_pathogens_do
        procedure :: do_ppdd                  => aed_pathogens_do_ppdd
        procedure :: do_benthos               => aed_pathogens_do_benthos
!       procedure :: get_light_extinction     => aed_pathogens_get_light_extinction
!       procedure :: get_conserved_quantities => aed_pathogens_get_conserved_quantities
   END TYPE

   real(rk), parameter :: secs_pr_day = 86400.

!===============================================================================
CONTAINS



!###############################################################################
FUNCTION aed_pathogens_create(namlst,name,parent) RESULT(self)
!-------------------------------------------------------------------------------
! Initialise the pathogen biogeochemical model
!
!  Here, the aed_p_m namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)          :: namlst
   CHARACTER(len=*),INTENT(in) :: name
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_aed_pathogens),POINTER :: self

   INTEGER            :: num_pathogens
   INTEGER            :: the_pathogens(MAX_PATHO_TYPES)


   NAMELIST /aed_pathogens/ num_pathogens, the_pathogens
!-----------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL initialize_model_info(self,name,parent)

   ! Read the namelist
   read(namlst,nml=aed_pathogens,err=99)

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
   call self%register_dependency(self%id_tem,standard_variables%temperature)
   call self%register_dependency(self%id_sal,standard_variables%practical_salinity)
   call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_photosynthetic_radiative_flux)

   RETURN

99 call fatal_error('aed_pathogens_init','Error reading namelist aed_pathogens')

END FUNCTION aed_pathogens_create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed_pathogens_load_params(self, count, list)
!-------------------------------------------------------------------------------
   _CLASS_ (type_aed_pathogens),INTENT(inout) :: self
   INTEGER,INTENT(in) :: count
   INTEGER,INTENT(in) :: list(*)

   INTEGER     :: i,tfil
   real(rk)    :: minPath

   TYPE(pathogen_nml_data) :: pd(MAX_PATHO_TYPES)
   NAMELIST /pathogen_data/ pd
!-------------------------------------------------------------------------------
!BEGIN
    minPath = 1e-10
    tfil = find_free_lun()
    open(tfil,file="aed_pathogen_pars.nml", status='OLD')
    read(tfil,nml=pathogen_data,err=99)
    close(tfil)

    self%num_pathogens = count
    ALLOCATE(self%pathogens(count))
    ALLOCATE(self%id_p(count))
    DO i=1,count
       ! Assign parameters from database to simulated groups
       !self%pathogens(i)%p_name       = pd(list(i))%p_name
       self%pathogens(i)          = pd(list(i))

       ! Register group as a state variable
       call self%register_state_variable(self%id_p(i),                           &
                             TRIM(self%pathogens(i)%p_name),                  &
                             'mmol/m**3', 'pathogen',                         &
                             minPath,                                             &
                            ! pd(list(i))%p_initial,                          &
                             minimum=minPath,                                           &
                             !minimum=pd(list(i))%p0,                         &
                             vertical_movement = self%pathogens(i)%coef_sett_w_path)


!      IF (self%pathogens(i)%p_name == 'crypto') THEN
!         ! Register IN group as a state variable
!         self%id_in(i) = self%register_state_variable(                        &
!                             TRIM(self%pathogens(i)%p_name)//'_dd',           &
!                             'mmol/m**3', 'pathogen dd',                      &
!                             0.0,   &
!                             0.0,                                  &
!                             vertical_movement = self%pathogens(i)%coef_sett_w_path)

!      ENDIF


    ENDDO

    RETURN

99 call fatal_error('aed_pathogens_load_params','Error reading namelist pathogen_data')
!
END SUBROUTINE aed_pathogens_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_pathogens_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of pathogen biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_pathogens),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   real(rk)           :: pth
   real(rk)           :: temp,par,Io,salinity
   real(rk)           :: growth,light,mortality, predation
   real(rk)           :: f_AOC,f_pH,f_DO,phi,lightBW,phstar

   INTEGER            :: pth_i,c

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

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

      growth    = 0.0_rk
      predation = 0.0_rk

      ! Natural mortality (as impacted by T, S, pH)
      f_AOC = 1.0 ! aoc / (K_AOC + aoc)
      f_pH  = 1.0 ! + c_PH * ( pH_star**delta / (pH_star**delta+K_PH**delta) )
      mortality = self%pathogens(pth_i)%coef_mort_kd20/86400.  &
                + (self%pathogens(pth_i)%coef_mort_c_SM*salinity**self%pathogens(pth_i)%coef_mort_alpha) &
                * ((1.0-f_AOC)**self%pathogens(pth_i)%coef_mort_beta) * f_pH
      mortality = mortality * (self%pathogens(pth_i)%coef_mort_theta**(temp-20.0))


      ! Sunlight inactivation (as impacted by S, DO and pH)
      light     = 0.0_rk
      lightBW   = 0.0_rk
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

   !-----------------------------------------------------------------
   ! export diagnostic variables
   ! _SET_DIAGNOSTIC_(self%id_?? ,??)


   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_pathogens_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_pathogens_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! !IROUTINE: Calculate pelagic sedimentation of pathogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_pathogens),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   real(rk) :: pth        ! State
   INTEGER  :: pth_i
   real(rk) :: pth_flux

   ! Parameters
!
!-------------------------------------------------------------------------------
!BEGIN

   _FABM_HORIZONTAL_LOOP_BEGIN_

   DO pth_i=1,self%num_pathogens
      ! Retrieve current (local) state variable values.
      _GET_(self%id_p(pth_i),pth) ! pathogen

      pth_flux = 0.0_rk  !self%pathogens(pth_i)%w_p*MAX(pth,0.0_rk)

     ! Set bottom fluxes for the pelagic (change per surface area per second)
     ! Transfer sediment flux value to FABM.
     _SET_BOTTOM_EXCHANGE_(self%id_p(pth_i),pth_flux)

   ENDDO

   ! Leave spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_END_
END SUBROUTINE aed_pathogens_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_pathogens_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of pathogen biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_pathogens),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   real(rk)            :: n,p,par,I_0
   INTEGER             :: pth_i
   real(rk)            :: iopt,rpd,primprod
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   DO pth_i=1,self%num_pathogens
      ! Retrieve current (local) state variable values.
!     _GET_(self%id_p(pth_i),p) ! pathogen

      ! Retrieve current environmental conditions.
!     _GET_   (self%id_par,par)  ! local photosynthetically active radiation
!     _GET_DEPENDENCY_HZ_(self%id_I_0,I_0)  ! surface short wave radiation

   ENDDO

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_pathogens_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_pathogens

#endif

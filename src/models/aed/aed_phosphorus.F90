!###############################################################################
!#                                                                             #
!# aed_phosphorus.F90                                                          #
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
!# Created 24 August 2011                                                      #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_phosphorus
!-------------------------------------------------------------------------------
! aed_phosphorus --- phosphorus biogeochemical model
!
! The AED module phosphorus contains equations that describe exchange of
! soluable reactive phosphorus across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_type_phosphorus
!
   TYPE,extends(type_base_model) :: aed_type_phosphorus
!     Variable identifiers
      TYPE (type_state_variable_id) :: id_frp, id_frpads, id_oxy,  id_tss, id_pH
      TYPE (type_horizontal_dependency_id) :: id_Fsed_frp
      TYPE (type_dependency_id)          :: id_temp, id_tssext
      TYPE (type_horizontal_diagnostic_variable_id) :: id_sed_frp

!     Model parameters
      AED_REAL :: Fsed_frp,Ksed_frp,theta_sed_frp      ! Benthic
      LOGICAL  :: ben_use_oxy,ben_use_aedsed
      AED_REAL :: Kpo4p,Kadsratio,Qmax, w_po4ads       ! Adsorption
      LOGICAL  :: simPO4Adsorption, ads_use_pH, ads_use_external_tss
      INTEGER  :: PO4AdsorptionModel

      CONTAINS     ! Model Procedures
        PROCEDURE :: initialize               => aed_init_phosphorus
        PROCEDURE :: do_benthos               => aed_phosphorus_do_benthos
        PROCEDURE :: check_state              => aed_phosphorus_check_state
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_phosphorus(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_type_phosphorus),TARGET,INTENT(inout) :: self

!
!LOCALS
   INTEGER  :: status

   AED_REAL          :: frp_initial   = 4.5
   ! Benthic
   AED_REAL          :: Fsed_frp      = 3.5
   AED_REAL          :: Ksed_frp      = 30.0
   AED_REAL          :: theta_sed_frp = 1.05
   CHARACTER(len=64) :: phosphorus_reactant_variable=''
   CHARACTER(len=64) :: Fsed_frp_variable=''
   ! Adsorption
   LOGICAL           :: simPO4Adsorption = .FALSE.
   LOGICAL           :: ads_use_external_tss = .FALSE.
   LOGICAL           :: ads_use_pH = .FALSE.
   AED_REAL          :: Kpo4p = 1.05
   AED_REAL          :: Kadsratio = 1.05
   AED_REAL          :: Qmax = 1.05
   AED_REAL          :: w_po4ads = 0.00
   INTEGER           :: PO4AdsorptionModel = 1
   CHARACTER(len=64) :: po4sorption_target_variable=''

   AED_REAL, parameter :: secs_pr_day = 86400.

   NAMELIST /aed_phosphorus/ frp_initial,Fsed_frp,Ksed_frp,theta_sed_frp,      &
                             phosphorus_reactant_variable,Fsed_frp_variable,   &
                             simPO4Adsorption,ads_use_external_tss,            &
                             po4sorption_target_variable, PO4AdsorptionModel,  &
                             ads_use_pH,Kpo4p,Kadsratio,Qmax,w_po4ads
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_phosphorus,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_phosphorus'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%Fsed_frp             = Fsed_frp/secs_pr_day
   self%Ksed_frp             = Ksed_frp
   self%theta_sed_frp        = theta_sed_frp
   self%simPO4Adsorption     = simPO4Adsorption
   self%ads_use_external_tss = ads_use_external_tss
   self%PO4AdsorptionModel   = PO4AdsorptionModel
   self%ads_use_pH           = ads_use_pH
   self%Kpo4p                = Kpo4p
   self%Kadsratio            = Kadsratio
   self%Qmax                 = Qmax
   self%w_po4ads             = w_po4ads/secs_pr_day


   ! Register main state variable
   CALL self%register_state_variable(self%id_frp, 'frp', 'mmol/m**3', 'phosphorus',     &
                                    frp_initial,minimum=zero_, no_river_dilution=.false.)

   ! Register external state variable dependencies (for benthic flux)
   self%ben_use_oxy = phosphorus_reactant_variable .NE. '' !This means oxygen module switched on
   IF (self%ben_use_oxy) THEN
     CALL self%register_state_dependency(self%id_oxy,phosphorus_reactant_variable)
   ENDIF

   self%ben_use_aedsed = Fsed_frp_variable .NE. '' !This means aed sediment module switched on
   IF (self%ben_use_aedsed) THEN
     CALL self%register_horizontal_dependency(self%id_Fsed_frp,Fsed_frp_variable)
   ENDIF

   ! Check if particles and PO4 adsorption are simulated
   IF (self%simPO4Adsorption) THEN
     IF (self%ads_use_external_tss) THEN
         PRINT *,'PO4 adsorption is configured to use external TSS'
         CALL self%register_dependency(self%id_tssext,standard_variables%mass_concentration_of_suspended_matter)
     ELSE
       IF (po4sorption_target_variable .NE. '' ) THEN
          CALL self%register_state_dependency(self%id_tss,po4sorption_target_variable)
       ELSE
          PRINT *,'PO4 adsorption is configured but no internal or external target variable is set'
          STOP
       ENDIF
     ENDIF

     CALL self%register_state_variable(self%id_frpads,'frp_ads','mmol/m**3','adsorbed phosphorus',     &
                      zero_,minimum=zero_,no_river_dilution=.false.,vertical_movement=self%w_po4ads)

     IF (self%ads_use_pH) THEN
       CALL self%register_state_dependency(self%id_pH,'aed_carbon_pH')
     ENDIF
   ENDIF

   ! Register diagnostic variables
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_frp,'sed_frp','mmol/m**2/d', &
                                         'Filterable reactive phosphorus')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature)
END SUBROUTINE aed_init_phosphorus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phosphorus_check_state(self,_ARGUMENTS_CHECK_STATE_)
!-------------------------------------------------------------------------------
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phosphorus),INTENT(in) :: self
   _DECLARE_ARGUMENTS_CHECK_STATE_
!
!LOCALS
   ! Environment
   AED_REAL :: temp, tss

   ! State
   AED_REAL :: frp,frpads,pH

   ! Temporary variables
   AED_REAL :: SSconc, PO4dis, PO4par, PO4tot, buffer, f_pH, K, Qm

   AED_REAL,PARAMETER :: one_e_neg_ten = 1e-10

!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT. self%simPO4Adsorption) RETURN

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current environmental conditions for the cell.
   _GET_(self%id_temp,temp)  ! local temperature
   IF(self%ads_use_external_tss) THEN
     _GET_(self%id_tssext,tss)  ! externally supplied total susp solids
   END IF

    ! Retrieve current (local) state variable values.
   _GET_(self%id_frp,frp)         ! dissolved PO4
   _GET_(self%id_frpads,frpads)   ! adsorped PO4
   IF (.NOT.self%ads_use_external_tss) THEN
      _GET_(self%id_tss,tss)       ! local total susp solids
   END IF


   PO4dis   = zero_
   PO4par   = zero_
   buffer   = zero_
   f_pH     = one_
   SSconc   = tss

   ! calculate the total possible PO4 for sorption, and solids
   PO4tot  = MAX(one_e_neg_ten, frp + frpads)  ! Co in Chao (mg)
   SSconc  = MAX(one_e_neg_ten, SSconc )       ! s in Chao  (mg = mol/L * g/mol * mg/g)


   IF(self%PO4AdsorptionModel == 1) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Ji 2008:
     !
     ! Ji, Z-G. 2008. Hydrodynamics and Water Quality. Wiley Press.
     !
     PO4par = (self%Kpo4p*SSconc) / (one_+self%Kpo4p*SSconc) * PO4tot
     PO4dis = one_ / (one_+self%Kpo4p*SSconc) * PO4tot
   ELSEIF(self%PO4AdsorptionModel == 2) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Chao et al. 2010:
     !
     ! Chao, X. et al. 2010. Three-dimensional numerical simulation of
     !   water quality and sediment associated processes with application
     !   to a Mississippi delta lake. J. Environ. Manage. 91 p1456-1466.
     !
     IF(self%ads_use_pH) THEN
       _GET_(self%id_pH,pH) ! pH

       IF(pH > 11.) pH = 11.0
       IF(pH < 3.)  pH = 3.0

       ! -0.0094x2 + 0.0428x + 0.9574
       ! (ursula.salmon@uwa.edu.au: fPH for PO4 sorption to Fe in Mine Lakes)
       f_pH = -0.0094*pH*pH + 0.0428*pH + 0.9574

     ELSE

       f_pH = one_

     END IF

     ! calculate particlate fraction based on quadratic solution
     K  = self%Kadsratio
     Qm = self%Qmax

     ! Chao Eq 16
     buffer = SQRT(((PO4tot+(1./K)-(SSconc*Qm*f_pH)))**2. + (4.*f_pH*SSconc*Qm/K))
     PO4par  = 0.5 * ((PO4tot+(1./K)+(SSconc*Qm*f_pH))  - buffer  )

     ! Check for stupid solutions
     IF(PO4par > PO4tot) PO4par = PO4tot
     IF(PO4par < zero_) PO4par = zero_

     ! Now set dissolved portion
     PO4dis = PO4tot - PO4par

   ELSE
     !-----------------------------------------------------
     ! No model is selected
     RETURN
   END IF

   _SET_STATE_(self%id_frp,PO4dis)       ! Dissolved PO4
   _SET_STATE_(self%id_frpads,PO4par)    ! Adsorped PO4

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_phosphorus_check_state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phosphorus_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED phosphorus.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_phosphorus),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   ! Environment
   AED_REAL :: temp

   ! State
   AED_REAL :: frp,oxy

   ! Temporary variables
   AED_REAL :: frp_flux, Fsed_frp

   ! Parameters
   AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_(self%id_frp,frp) ! phosphorus

   IF (self%ben_use_aedsed) THEN
      _GET_HORIZONTAL_(self%id_Fsed_frp,Fsed_frp)
   ELSE
      Fsed_frp = self%Fsed_frp
   ENDIF

   IF (self%ben_use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      _GET_(self%id_oxy,oxy)
      frp_flux = Fsed_frp * self%Ksed_frp/(self%Ksed_frp+oxy) * (self%theta_sed_frp**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
      frp_flux = Fsed_frp * (self%theta_sed_frp**(temp-20.0))
   ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to FABM.
   !_SET_BOTTOM_FLUX_(self%id_frp,frp_flux/secs_pr_day)
   !_SET_ODE_SED_FLUX_(self%id_frp,frp_flux)
   _SET_BOTTOM_EXCHANGE_(self%id_frp,frp_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_frp,-frp_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_frp,frp_flux)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

END SUBROUTINE aed_phosphorus_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_phosphorus

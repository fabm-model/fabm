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

#ifdef _FABM_F2003_

#include "aed.h"

!
MODULE aed_phosphorus
!-------------------------------------------------------------------------------
! aed_phosphorus --- phosphorus biogeochemical model
!
! The AED module phosphorus contains equations that describe exchange of
! soluable reactive phosphorus across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE fabm_types

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_phosphorus
!
   TYPE,extends(type_base_model) :: type_aed_phosphorus
!     Variable identifiers
      type (type_state_variable_id) :: id_frp, id_frpads, id_oxy,  id_tssfabm, id_pH
      type (type_bottom_state_variable_id) :: id_Fsed_frp
      type (type_dependency_id)          :: id_temp, id_tssext
      type (type_horizontal_diagnostic_variable_id) :: id_sed_frp
      type (type_conserved_quantity_id)  :: id_totP

!     Model parameters
      real(rk) :: Fsed_frp,Ksed_frp,theta_sed_frp      ! Benthic
      LOGICAL  :: ben_use_oxy,ben_use_aedsed
      real(rk) :: Kpo4p,Kadsratio,Qmax, w_po4ads       ! Adsorption
      LOGICAL  :: simPO4Adsorption, ads_use_pH, ads_use_external_tss
      INTEGER  :: PO4AdsorptionModel

      CONTAINS     ! Model Procedures
        procedure :: initialize               => aed_phosphorus_init
        procedure :: do                       => aed_phosphorus_do
        procedure :: do_ppdd                  => aed_phosphorus_do_ppdd
        procedure :: do_benthos               => aed_phosphorus_do_benthos
        procedure :: get_conserved_quantities => aed_phosphorus_get_conserved_quantities
   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_phosphorus_init(self,configunit)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_aed_phosphorus),TARGET,INTENT(INOUT) :: self
   INTEGER,INTENT(in)                               :: configunit
!
!LOCALS

   real(rk)          :: frp_initial   = 4.5
   ! Benthic
   real(rk)          :: Fsed_frp      = 3.5
   real(rk)          :: Ksed_frp      = 30.0
   real(rk)          :: theta_sed_frp = 1.05
   CHARACTER(len=64) :: phosphorus_reactant_variable=''
   CHARACTER(len=64) :: Fsed_frp_variable=''
   ! Adsorption
   LOGICAL           :: simPO4Adsorption = .FALSE.
   LOGICAL           :: ads_use_external_tss = .FALSE.
   LOGICAL           :: ads_use_pH = .FALSE.
   real(rk)          :: Kpo4p = 1.05
   real(rk)          :: Kadsratio = 1.05
   real(rk)          :: Qmax = 1.05
   real(rk)          :: w_po4ads = 0.00
   INTEGER           :: PO4AdsorptionModel = 1
   CHARACTER(len=64) :: po4sorption_target_variable=''

   real(rk), parameter :: secs_pr_day = 86400.

   NAMELIST /aed_phosphorus/ frp_initial,Fsed_frp,Ksed_frp,theta_sed_frp,      &
                             phosphorus_reactant_variable,Fsed_frp_variable,   &
                             simPO4Adsorption,ads_use_external_tss,            &
                             po4sorption_target_variable, PO4AdsorptionModel,  &
                             ads_use_pH,Kpo4p,Kadsratio,Qmax,w_po4ads
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(configunit,nml=aed_phosphorus,err=99)

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
   call self%register_state_variable(self%id_frp,'frp','mmol/m**3','phosphorus',     &
                                    frp_initial,minimum=0.0_rk,no_river_dilution=.false.)

   ! Register external state variable dependencies (for benthic flux)
   self%ben_use_oxy = phosphorus_reactant_variable .NE. '' !This means oxygen module switched on
   IF (self%ben_use_oxy) THEN
     call self%register_state_dependency(self%id_oxy,phosphorus_reactant_variable)
   ENDIF

   self%ben_use_aedsed = Fsed_frp_variable .NE. '' !This means aed sediment module switched on
   IF (self%ben_use_aedsed) THEN
     call self%register_bottom_state_dependency(self%id_Fsed_frp,Fsed_frp_variable)
   ENDIF

   ! Check if particles and PO4 adsorption are simulated
   IF (self%simPO4Adsorption) THEN

     IF (self%ads_use_external_tss) THEN

         PRINT *,'PO4 adsorption is configured to use external TSS'
         call self%register_dependency(self%id_tssext,standard_variables%mass_concentration_of_suspended_matter)

     ELSE

       IF (po4sorption_target_variable .NE. '' ) THEN
         call self%register_state_dependency(self%id_tssfabm,po4sorption_target_variable)
       ELSE
         PRINT *,'PO4 adsorption is configured but no internal or external target variable is set'
       END IF

     ENDIF

     call self%register_state_variable(self%id_frpads,'frp_ads','mmol/m**3','adsorbed phosphorus',     &
                      0.0_rk,minimum=0.0_rk,no_river_dilution=.false.,vertical_movement=self%w_po4ads )

     IF (self%ads_use_pH) THEN
       call self%register_state_dependency(self%id_pH,'aed_carbon_pH')
     ENDIF

   ENDIF

   ! Register diagnostic variables
   call self%register_horizontal_diagnostic_variable(self%id_sed_frp,'sed_frp','mmol/m**2/d', &
                                         'Filterable reactive phosphorus',           &
                     time_treatment=time_treatment_step_integrated)

   ! Register conserved quantities
   call self%register_conserved_quantity(self%id_totP,'TP','mmol/m**3','Total phosphorus')

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   RETURN

99 CALL self%fatal_error('aed_phosphorus_init','Error reading namelist aed_phosphorus')

END SUBROUTINE aed_phosphorus_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phosphorus_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_phosphorus model
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_phosphorus),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
!  real(rk)           :: frp,frpads,diff_frp
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Just stuck it here for now until generic "update_state" interface is made
   CALL aed_phosphorus_update_state(self,_FABM_ARGS_DO_RHS_)

 !  ! Retrieve current (local) state variable values.
 !  _GET_(self%id_frp,frp) ! phosphorus
 !  IF(self%simPO4Adsorption) THEN
 !    _GET_(self%id_frpads,frpads) ! phosphorus
 !  END IF
 !
 !
 !  ! Set temporal derivatives
 !  diff_frp = 0.
 !
 !  _SET_ODE_(self%id_frp,diff_frp)
 !  IF(self%simPO4Adsorption) THEN
 !    _SET_ODE_(self%id_frpads,0.0_rk)
 !  END IF
 !
   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
END SUBROUTINE aed_phosphorus_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_phosphorus_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of phosphorus biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_phosphorus),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   real(rk)           :: frp
   real(rk)           :: diff_frp
   real(rk),PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_frp,frp) ! phosphorus

   ! Set temporal derivatives
   diff_frp = 0.

   _SET_PP_(self%id_frp,self%id_frp,diff_frp)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_phosphorus_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_phosphorus_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED phosphorus.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_phosphorus),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   ! Environment
   real(rk) :: temp

   ! State
   real(rk) :: frp,oxy

   ! Temporary variables
   real(rk) :: frp_flux, Fsed_frp

   ! Parameters
   real(rk),PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

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
   _FABM_HORIZONTAL_LOOP_END_

END SUBROUTINE aed_phosphorus_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phosphorus_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities (currently only phosphorus)
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_phosphorus),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   real(rk) :: frp
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_frp,frp) ! phosphorus

   ! Total nutrient is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totP,frp)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_phosphorus_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_phosphorus_update_state(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Update partitioning of phosphate between dissolved and particulate pools
! after kinetic transformations are applied
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_phosphorus),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   ! Environment
   real(rk) :: temp, tss

   ! State
   real(rk) :: frp,frpads,pH

   ! Temporary variables
   real(rk) :: SSconc, PO4dis, PO4par, PO4tot, buffer, f_pH, K, Qm

   real(rk),PARAMETER :: one_e_neg_ten = 1e-10

!-------------------------------------------------------------------------------
!BEGIN
   IF(.NOT.self%simPO4Adsorption) RETURN

   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the cell.
   _GET_(self%id_temp,temp)  ! local temperature
   IF(self%ads_use_external_tss) THEN
     _GET_(self%id_tssext,tss)  ! externally supplied total susp solids
   END IF

    ! Retrieve current (local) state variable values.
   _GET_(self%id_frp,frp)         ! dissolved PO4
   _GET_(self%id_frpads,frpads)   ! adsorped PO4
   IF(.NOT.self%ads_use_external_tss) THEN
     _GET_(self%id_tssfabm,tss)       ! local total susp solids
   END IF


   PO4dis   = 0.0_rk
   PO4par   = 0.0_rk
   buffer   = 0.0_rk
   f_pH     = 1.0_rk
   SSconc   = tss

     ! calculate the total possible PO4 for sorption, and solids
     PO4tot  = MAX(one_e_neg_ten, frp + frpads)  ! Co in Chao (mg)
     SSconc  = MAX(one_e_neg_ten, SSconc )       ! s in Chao  (mg = mol/L * g/mol * mg/g)


   IF(self%PO4AdsorptionModel == 1) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Ji 2008:
     !
     ! Ji, Z-G. 2008. Hydrodynamics and Water Quality. Wiley Press.

     PO4par = (self%Kpo4p*SSconc) / (1.0_rk+self%Kpo4p*SSconc) * PO4tot
     PO4dis = 1.0_rk / (1.0_rk+self%Kpo4p*SSconc) * PO4tot


   ELSEIF(self%PO4AdsorptionModel == 2) THEN
     !-----------------------------------------------------
     ! This is the model for PO4 sorption from Chao et al. 2010:
     !
     ! Chao, X. et al. 2010. Three-dimensional numerical simulation of
     !   water quality and sediment associated processes with application
     !   to a Mississippi delta lake. J. Environ. Manage. 91 p1456-1466.


     IF(self%ads_use_pH) THEN
       _GET_(self%id_pH,pH) ! pH

       IF(pH > 11.) pH = 11.0
       IF(pH < 3.)  pH = 3.0

       ! -0.0094x2 + 0.0428x + 0.9574
       ! (ursula.salmon@uwa.edu.au: fPH for PO4 sorption to Fe in Mine Lakes)
       f_pH = -0.0094*pH*pH + 0.0428*pH + 0.9574

     ELSE

       f_pH = 1.0_rk

     END IF

     ! calculate particlate fraction based on quadratic solution
     K  = self%Kadsratio
     Qm = self%Qmax

     ! Chao Eq 16
     buffer = SQRT(((PO4tot+(1./K)-(SSconc*Qm*f_pH)))**2. + (4.*f_pH*SSconc*Qm/K))
     PO4par  = 0.5 * ((PO4tot+(1./K)+(SSconc*Qm*f_pH))  - buffer  )

     ! Check for stupid solutions
     IF(PO4par > PO4tot) PO4par = PO4tot
     IF(PO4par < 0.0_rk) PO4par = 0.0_rk

     ! Now set dissolved portion
     PO4dis = PO4tot - PO4par

   ELSE
     !-----------------------------------------------------
     ! No model is selected
     RETURN
   END IF

!   _SET_STATE_(self%id_frp,PO4dis)       ! Dissolved PO4
!   _SET_STATE_(self%id_frpads,PO4par)    ! Adsorped PO4

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
END SUBROUTINE aed_phosphorus_update_state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_phosphorus
#endif

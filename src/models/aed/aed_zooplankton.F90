!###############################################################################
!#                                                                             #
!# aed_zooplankton.F90                                                         #
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
!# Created October 2011                                                        #
!#                                                                             #
!###############################################################################

#include "aed.h"

MODULE aed_zooplankton
!-------------------------------------------------------------------------------
!  aed_zooplankton --- multi zooplankton biogeochemical model
!-------------------------------------------------------------------------------
   USE fabm_types
   USE aed_util,ONLY : find_free_lun,aed_bio_temp_function, fTemp_function,qsort
   USE aed_phytoplankton,ONLY : type_aed_phytoplankton

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC type_aed_zooplankton
!
   TYPE type_zoop_prey
      !State variable name for zooplankton prey
      CHARACTER(64) :: zoop_prey
      !Preference factors for zooplankton predators grazing on prey
      real(rk)      :: Pzoo_prey
   END TYPE type_zoop_prey


   TYPE type_zoop_params
      ! General Attributes
      CHARACTER(64) :: zoop_name
      real(rk) :: zoop_initial, min_zoo
      ! Growth rate parameters
      real(rk) :: Rgrz_zoo, fassim_zoo, Kgrz_zoo, theta_grz_zoo
      ! Respiration, mortaility and excretion parameters
      real(rk) :: Rresp_zoo, Rmort_zoo, ffecal_zoo, fexcr_zoo, ffecal_sed
      ! Temperature limitation on zooplankton loss terms
      real(rk) :: theta_resp_zoo, Tstd_zoo, Topt_zoo, Tmax_zoo
      ! Salinity parameters
      INTEGER  :: saltfunc_zoo
      real(rk) :: Smin_zoo, Smax_zoo, Sint_zoo
      ! Nutrient parameters
       real(rk) :: INC_zoo, IPC_zoo
      ! Dissolved oxygen parameters
      real(rk) :: DOmin_zoo
      ! Minumum prey concentration parameters
      real(rk) :: Cmin_grz_zoo
      ! Prey information
      INTEGER  :: num_prey
      TYPE(type_zoop_prey)      :: prey(MAX_ZOOP_PREY)
      INTEGER  :: simDOlim
      ! Temperature limitation derived terms
      real(rk) :: kTn, aTn, bTn
   END TYPE

   TYPE,extends(type_zoop_params) :: type_zoop_data
      type (type_state_variable_id)  :: id_prey(MAX_ZOOP_PREY)
      type (type_state_variable_id)  :: id_phyIN(MAX_ZOOP_PREY),id_phyIP(MAX_ZOOP_PREY)
   END TYPE

   TYPE,extends(type_base_model) :: type_aed_zooplankton
!     Variable identifiers
      type (type_state_variable_id)      :: id_zoo(MAX_ZOOP_TYPES)
      type (type_state_variable_id)      :: id_Nexctarget,id_Nmorttarget
      type (type_state_variable_id)      :: id_Pexctarget,id_Pmorttarget
      type (type_state_variable_id)      :: id_Cexctarget,id_Cmorttarget
      type (type_state_variable_id)      :: id_DOupttarget
      type (type_dependency_id)          :: id_tem, id_sal, id_extc
      type (type_diagnostic_variable_id) :: id_grz,id_resp,id_mort
      type (type_conserved_quantity_id)  :: id_totN,id_totP,id_totC
      type (type_conserved_quantity_id)  :: id_totZOO


!     Model parameters
      INTEGER                                   :: num_zoops
      TYPE(type_zoop_data),dimension(:),allocatable :: zoops
      LOGICAL  :: simDNexcr, simDPexcr, simDCexcr
      LOGICAL  :: simPNexcr, simPPexcr, simPCexcr

      CONTAINS     ! Model Methods
        procedure :: initialize               => aed_zooplankton_init
        procedure :: do                       => aed_zooplankton_do
        procedure :: do_ppdd                  => aed_zooplankton_do_ppdd
!       procedure :: do_benthos               => aed_zooplankton_do_benthos
        procedure :: get_conserved_quantities => aed_zooplankton_get_conserved_quantities
   END TYPE

   LOGICAL :: debug = .TRUE.
   !lcb do we need INTEGER  :: ino3, inh4,idon, in2, ifrp, idop

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE aed_zooplankton_load_params(self, count, list)
!-------------------------------------------------------------------------------
   CLASS (type_aed_zooplankton),INTENT(inout) :: self
   INTEGER,INTENT(in)                         :: count !Number of zooplankton groups
   INTEGER,INTENT(in)                         :: list(*) !List of zooplankton groups to simulate

   INTEGER  :: i,j,tfil,sort_i(MAX_ZOOP_PREY)
   real(rk) :: Pzoo_prey(MAX_ZOOP_PREY)

   real(rk),PARAMETER :: secs_pr_day = 86400.
   TYPE(type_zoop_params)  :: zoop_param(MAX_ZOOP_TYPES)
   NAMELIST /zoop_params/ zoop_param
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file="aed_zoop_pars.nml", status='OLD')
    read(tfil,nml=zoop_params,err=99)
    close(tfil)

    self%num_zoops = count
    allocate(self%zoops(count))
    DO i=1,count
       ! Assign parameters from database to simulated groups
       self%zoops(i)%zoop_name         = zoop_param(list(i))%zoop_name
       self%zoops(i)%zoop_initial      = zoop_param(list(i))%zoop_initial
       self%zoops(i)%min_zoo           = zoop_param(list(i))%min_zoo
       self%zoops(i)%Rgrz_zoo          = zoop_param(list(i))%Rgrz_zoo/secs_pr_day
       self%zoops(i)%fassim_zoo        = zoop_param(list(i))%fassim_zoo
       self%zoops(i)%Kgrz_zoo          = zoop_param(list(i))%Kgrz_zoo
       self%zoops(i)%theta_grz_zoo     = zoop_param(list(i))%theta_grz_zoo
       self%zoops(i)%Rresp_zoo         = zoop_param(list(i))%Rresp_zoo/secs_pr_day
       self%zoops(i)%Rmort_zoo         = zoop_param(list(i))%Rmort_zoo/secs_pr_day
       self%zoops(i)%ffecal_zoo        = zoop_param(list(i))%ffecal_zoo
       self%zoops(i)%fexcr_zoo         = zoop_param(list(i))%fexcr_zoo
       self%zoops(i)%ffecal_sed        = zoop_param(list(i))%ffecal_sed
       self%zoops(i)%theta_resp_zoo    = zoop_param(list(i))%theta_resp_zoo
       self%zoops(i)%Tstd_zoo          = zoop_param(list(i))%Tstd_zoo
       self%zoops(i)%Topt_zoo          = zoop_param(list(i))%Topt_zoo
       self%zoops(i)%Tmax_zoo          = zoop_param(list(i))%Tmax_zoo
       self%zoops(i)%saltfunc_zoo      = zoop_param(list(i))%saltfunc_zoo
       self%zoops(i)%Smin_zoo          = zoop_param(list(i))%Smin_zoo
       self%zoops(i)%Smax_zoo          = zoop_param(list(i))%Smax_zoo
       self%zoops(i)%Sint_zoo          = zoop_param(list(i))%Sint_zoo
       self%zoops(i)%INC_zoo           = zoop_param(list(i))%INC_zoo
       self%zoops(i)%IPC_zoo           = zoop_param(list(i))%IPC_zoo
       self%zoops(i)%simDOlim          = zoop_param(list(i))%simDOlim
       self%zoops(i)%DOmin_zoo         = zoop_param(list(i))%DOmin_zoo
       self%zoops(i)%Cmin_grz_zoo      = zoop_param(list(i))%Cmin_grz_zoo
       self%zoops(i)%num_prey          = zoop_param(list(i))%num_prey
       !Loop through prey variables assigning a target variable and preference factor
       !First sort in decending order of food preferences
       DO j=1,self%zoops(i)%num_prey
          sort_i(j) = j
          Pzoo_prey(j) = zoop_param(list(i))%prey(j)%Pzoo_prey
       ENDDO
       CALL qsort(Pzoo_prey,sort_i,1,self%zoops(i)%num_prey)
       DO j=1,self%zoops(i)%num_prey
          self%zoops(i)%prey(j)%zoop_prey = zoop_param(list(i))%prey(sort_i(self%zoops(i)%num_prey-j+1))%zoop_prey
          self%zoops(i)%prey(j)%Pzoo_prey = zoop_param(list(i))%prey(sort_i(self%zoops(i)%num_prey-j+1))%Pzoo_prey
       ENDDO

       ! Register group as a state variable
       CALL self%register_state_variable(self%id_zoo(i),           &
                              zoop_param(list(i))%zoop_name,       &
                              'mmolC/m**3', 'zooplankton',         &
                              zoop_param(list(i))%zoop_initial,    &
                              minimum=zoop_param(list(i))%min_zoo)
    ENDDO

    RETURN

99 CALL self%fatal_error('aed_zooplankton_load_params','Error reading namelist zoop_params')
!
END SUBROUTINE aed_zooplankton_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_zooplankton_init(self,configunit)
!-------------------------------------------------------------------------------
! Initialise the zooplankton biogeochemical model
!
!  Here, the aed_zooplankton namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_aed_zooplankton),TARGET,INTENT(INOUT) :: self
   INTEGER,INTENT(in)                                :: configunit
!
!LOCALS

   INTEGER            :: num_zoops
   INTEGER            :: the_zoops(MAX_ZOOP_TYPES)

   CHARACTER(len=64)  :: dn_target_variable='' !dissolved nitrogen target variable
   CHARACTER(len=64)  :: pn_target_variable='' !particulate nitrogen target variable
   CHARACTER(len=64)  :: dp_target_variable='' !dissolved phosphorus target variable
   CHARACTER(len=64)  :: pp_target_variable='' !particulate phosphorus target variable
   CHARACTER(len=64)  :: dc_target_variable='' !dissolved carbon target variable
   CHARACTER(len=64)  :: pc_target_variable='' !particulate carbon target variable

   real(rk),PARAMETER :: secs_pr_day = 86400.
   INTEGER            :: zoop_i, prey_i, phy_i

   NAMELIST /aed_zooplankton/ num_zoops, the_zoops, &
                    dn_target_variable, pn_target_variable, dp_target_variable, &
                    pp_target_variable, dc_target_variable, pc_target_variable
!-----------------------------------------------------------------------
!BEGIN
!print *,'**** Reading /aed_zooplankton/ namelist'
   ! Read the namelist
   read(configunit,nml=aed_zooplankton,err=99)

    self%num_zoops = 0
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed_zooplankton_load_params to values per second.
   CALL aed_zooplankton_load_params(self, num_zoops, the_zoops)


   CALL aed_bio_temp_function(self%num_zoops,            &
                              self%zoops%theta_resp_zoo, &
                              self%zoops%Tstd_zoo,       &
                              self%zoops%Topt_zoo,       &
                              self%zoops%Tmax_zoo,       &
                              self%zoops%aTn,            &
                              self%zoops%bTn,            &
                              self%zoops%kTn,            &
                              self%zoops%zoop_name)


   !Register link to prey state variables
   DO zoop_i = 1,num_zoops
      phy_i = 0
      DO prey_i = 1,self%zoops(zoop_i)%num_prey
          CALL self%register_state_dependency(self%zoops(zoop_i)%id_prey(prey_i), &
                                       self%zoops(zoop_i)%prey(prey_i)%zoop_prey)
          !If the zooplankton prey is phytoplankton then also register state dependency on
          !internal nitrogen and phosphorus
          IF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:17).EQ.'aed_phytoplankton') THEN
              phy_i = phy_i + 1
              CALL self%register_state_dependency(self%zoops(zoop_i)%id_phyIN(phy_i), &
                                       TRIM(self%zoops(zoop_i)%prey(prey_i)%zoop_prey)//'_IN')
              CALL self%register_state_dependency(self%zoops(zoop_i)%id_phyIP(phy_i), &
                                       TRIM(self%zoops(zoop_i)%prey(prey_i)%zoop_prey)//'_IP')

          ENDIF
      ENDDO
   ENDDO

   ! Register link to nutrient pools, if variable names are provided in namelist.
   self%simDNexcr = dn_target_variable .NE. ''
   IF (self%simDNexcr) THEN
     CALL self%register_state_dependency(self%id_Nexctarget,dn_target_variable)
   ENDIF
   self%simDPexcr = dp_target_variable .NE. ''
   IF (self%simDPexcr) THEN
     CALL self%register_state_dependency(self%id_Pexctarget,dp_target_variable)
   ENDIF
   self%simDCexcr = dc_target_variable .NE. ''
   IF (self%simDCexcr) THEN
     CALL self%register_state_dependency(self%id_Cexctarget,dc_target_variable)
   ENDIF

   self%simPNexcr = pn_target_variable .NE. ''
   IF (self%simPNexcr) THEN
     CALL self%register_state_dependency(self%id_Nmorttarget,pn_target_variable)
   ENDIF
   self%simPPexcr = pp_target_variable .NE. ''
   IF (self%simPPexcr) THEN
     CALL self%register_state_dependency(self%id_Pmorttarget,pp_target_variable)
   ENDIF
   self%simPCexcr = pc_target_variable .NE. ''
   IF (self%simPCexcr) THEN
     CALL self%register_state_dependency(self%id_Cmorttarget,pc_target_variable)
   ENDIF


   ! Register diagnostic variables
   CALL self%register_diagnostic_variable(self%id_grz,'grz','mmolC/m**3',  'net zooplankton grazing',           &
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_resp,'resp','mmolC/m**3',  'net zooplankton respiration',           &
                     time_treatment=time_treatment_averaged)
   CALL self%register_diagnostic_variable(self%id_mort,'mort','mmolC/m**3/d','net zooplankton mortality',      &
                     time_treatment=time_treatment_averaged)

   ! Register conserved quantities
   CALL self%register_conserved_quantity(self%id_totZOO,'TZOO','mmolC/m**3','Total zooplankton')
   CALL self%register_conserved_quantity(self%id_totN,'TN','mmol/m**3','Total nitrogen')
   CALL self%register_conserved_quantity(self%id_totP,'TP','mmol/m**3','Total phosphorus')
   CALL self%register_conserved_quantity(self%id_totC,'TC','mmol/m**3','Total carbon')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_tem,standard_variables%temperature)
   CALL self%register_dependency(self%id_sal,standard_variables%practical_salinity)
   CALL self%register_dependency(self%id_extc,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)

   RETURN

99 CALL self%fatal_error('aed_zooplankton_init','Error reading namelist aed_zooplankton')

END SUBROUTINE aed_zooplankton_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_zooplankton_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of zooplankton biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_aed_zooplankton),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   real(rk)           :: zoo,temp,salinity !State variables
   real(rk)           :: prey(MAX_ZOOP_PREY), grazing_prey(MAX_ZOOP_PREY) !Prey state variables
   real(rk)           :: phy_INcon(MAX_ZOOP_PREY), phy_IPcon(MAX_ZOOP_PREY) !Internal nutrients for phytoplankton
   real(rk)           :: dn_excr, dp_excr, dc_excr !Excretion state variables
   real(rk)           :: pon, pop, poc !Mortaility and fecal pellet state variables
   real(rk)           :: FGrazing_Limitation, f_T, f_Salinity
   real(rk)           :: pref_factor, Ctotal_prey !total concentration of available prey
   real(rk)           :: food, grazing, respiration, mortality !Growth & decay functions
   real(rk)           :: grazing_n, grazing_p !Grazing on nutrients
   real(rk)           :: pon_excr, pop_excr, poc_excr !POM excretion rates
   real(rk)           :: don_excr, dop_excr, doc_excr, delta_C !DOM excretion rates
   INTEGER            :: zoop_i,prey_i,prey_j,phy_i
   real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current environmental conditions.
   _GET_(self%id_tem,temp)     ! local temperature
   _GET_(self%id_sal,salinity) ! local salinity

   ! Retrieve current (local) state variable values.
   IF (self%simDNexcr) _GET_(self%id_Nexctarget, dn_excr)
   IF (self%simDPexcr) _GET_(self%id_Pexctarget, dp_excr)
   IF (self%simDCexcr) _GET_(self%id_Cexctarget, dc_excr)

   IF (self%simPNexcr) _GET_(self%id_Nmorttarget, pon)
   IF (self%simPPexcr) _GET_(self%id_Pmorttarget, pop)
   IF (self%simPCexcr) _GET_(self%id_Cmorttarget, poc)

   DO zoop_i=1,self%num_zoops

      ! Retrieve this zooplankton group
      _GET_(self%id_zoo(zoop_i),zoo)
      !Retrieve prey groups
      Ctotal_prey   = 0.0_rk
      DO prey_i=1,self%zoops(zoop_i)%num_prey
         _GET_(self%zoops(zoop_i)%id_prey(prey_i),prey(prey_i))
         Ctotal_prey = Ctotal_prey + prey(prey_i)
      ENDDO

      grazing       = 0.0_rk
      respiration   = 0.0_rk
      mortality     = 0.0_rk

      ! Get the grazing limitation function
       fGrazing_Limitation = fPrey_Limitation(self,zoop_i,Ctotal_prey)

      ! Get the temperature function
       f_T = fTemp_function(1, self%zoops(zoop_i)%Tmax_zoo,       &
                               self%zoops(zoop_i)%Tstd_zoo,       &
                               self%zoops(zoop_i)%theta_resp_zoo, &
                               self%zoops(zoop_i)%aTn,            &
                               self%zoops(zoop_i)%bTn,            &
                               self%zoops(zoop_i)%kTn, temp)

      ! Get the salinity limitation.
       f_Salinity = fSalinity_Limitation(self,zoop_i,salinity)

      ! Get the growth rate (/ s)
      ! grazing is in units of mass consumed/mass zoops/unit time
      grazing = self%zoops(zoop_i)%Rgrz_zoo * fGrazing_Limitation * f_T

      ! Now dertermine available prey and limit grazing amount to
      ! availability of prey
      ! food is total amount of food in units of mass/unit volume/unit time
      food = grazing * zoo
      IF (Ctotal_prey < self%zoops(zoop_i)%num_prey * self%zoops(zoop_i)%Cmin_grz_zoo ) THEN
          food = 0.0_rk
          grazing = food / zoo
      ELSEIF (food > Ctotal_prey - self%zoops(zoop_i)%num_prey * self%zoops(zoop_i)%Cmin_grz_zoo ) THEN
          food = Ctotal_prey - self%zoops(zoop_i)%num_prey * self%zoops(zoop_i)%Cmin_grz_zoo
          grazing = food / zoo
      ENDIF


      ! Now determine prey composition based on preference factors and
      ! availability of prey

      ! Prey has been ordered in grazing preference
      ! So take food in order of preference up to availability minus
      !value of minimum residual
      ! grazing_prey is in units of mass consumed/unit volumne/unit time

      DO prey_i = 1,self%zoops(zoop_i)%num_prey
          !Add up preferences for remaining prey
          pref_factor = 0.0_rk
          DO prey_j = prey_i,self%zoops(zoop_i)%num_prey
             pref_factor = pref_factor + self%zoops(zoop_i)%prey(prey_j)%Pzoo_prey
          ENDDO
          IF (food * self%zoops(zoop_i)%prey(prey_i)%Pzoo_prey / pref_factor <= &
                        prey(prey_i) - self%zoops(zoop_i)%Cmin_grz_zoo) THEN
             !Take fraction of left over food based on preference factor
             grazing_prey(prey_i) = food * self%zoops(zoop_i)%prey(prey_i)%Pzoo_prey / pref_factor
          ELSEIF (prey(prey_i) > self%zoops(zoop_i)%Cmin_grz_zoo) THEN
             grazing_prey(prey_i) = prey(prey_i) - self%zoops(zoop_i)%Cmin_grz_zoo
          ELSE
             grazing_prey(prey_i) = 0.0_rk
          ENDIF
          !Food remaining after grazing from current prey
          food = food - grazing_prey(prey_i)
      ENDDO

      ! Now determine nutrient composition of food based on prey type
      ! At this stage only the AED model state variables have multiple
      ! nutrients (C,N&P) so assume all others have a single nutrient
      ! and thus not need to calculate nutrient excretion as is taken
      ! care of in the respiration term.  22/12/2011
      ! grazing_n is in units of mass N consumed/unit volume/unit time
      ! grazing_p is in units of mass P consumed/unit volume/unit time

      grazing_n = 0.0_rk
      grazing_p = 0.0_rk
      phy_i = 0
      DO prey_i = 1,self%zoops(zoop_i)%num_prey
         IF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey .EQ. 'aed_organic_matter_poc') THEN
            IF (poc > 0.0_rk) THEN
                grazing_n = grazing_n + grazing_prey(prey_i) * pon/poc
                grazing_p = grazing_p + grazing_prey(prey_i) * pop/poc
            ELSE
                grazing_n = 0.0_rk
                grazing_p = 0.0_rk
            ENDIF
         ELSEIF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:17).EQ.'aed_phytoplankton') THEN
            phy_i = phy_i + 1
            _GET_(self%zoops(zoop_i)%id_phyIN(phy_i),phy_INcon(phy_i))
            _GET_(self%zoops(zoop_i)%id_phyIP(phy_i),phy_IPcon(phy_i))
            grazing_n = grazing_n + grazing_prey(prey_i) / prey(prey_i) * phy_INcon(phy_i) /14.0
            grazing_p = grazing_p + grazing_prey(prey_i) / prey(prey_i) * phy_IPcon(phy_i) /31.0
         ELSEIF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:15).EQ.'aed_zooplankton') THEN
            grazing_n = grazing_n + grazing_prey(prey_i) * self%zoops(zoop_i)%INC_zoo
            grazing_p = grazing_p + grazing_prey(prey_i) * self%zoops(zoop_i)%IPC_zoo
         ENDIF
      ENDDO


      ! Get the respiration rate (/ s)
      respiration = self%zoops(zoop_i)%Rresp_zoo * f_Salinity

      ! Get the mortality rate (/ s)
      mortality = self%zoops(zoop_i)%Rmort_zoo * f_T

      ! Don't excrete or die if we are at the min biomass otherwise we have a
      ! mass conservation leak in the C mass balance
      IF (zoo <= self%zoops(zoop_i)%min_zoo) THEN
        respiration = 0.0_rk
        mortality = 0.0_rk
      ENDIF

      ! Now we know the rates of carbon consumption and excretion,
      ! calculate rates of n & p excretion to maintain internal
      ! nutrient stores
      ! Calculate excretion of particulate organic matter - Units mmol/s
      poc_excr = ((1 - self%zoops(zoop_i)%fassim_zoo) * grazing + &
               (1 - self%zoops(zoop_i)%ffecal_sed) * self%zoops(zoop_i)%ffecal_zoo * respiration +  &
                                                        mortality) * zoo
      pon_excr = (1 - self%zoops(zoop_i)%fassim_zoo) * grazing_n  + &
               ((1 - self%zoops(zoop_i)%ffecal_sed) * self%zoops(zoop_i)%ffecal_zoo * respiration + &
                              mortality) * self%zoops(zoop_i)%INC_zoo * zoo
      pop_excr = (1 - self%zoops(zoop_i)%fassim_zoo) * grazing_p + &
               ((1 - self%zoops(zoop_i)%ffecal_sed) * self%zoops(zoop_i)%ffecal_zoo * respiration + &
                              mortality) * self%zoops(zoop_i)%IPC_zoo * zoo


      ! Calculate rate of change of zooplankton carbon (mmolC/s)          !
      delta_C = (self%zoops(zoop_i)%fassim_zoo * grazing - respiration - mortality) * zoo
      ! Calculate nutrient excretion require to balance internal nutrient store
      ! Note pon_excr includes loss due to messy feeding so no need to include assimilation fraction on grazing_n & grazing_p
      don_excr = grazing_n - pon_excr - delta_C * self%zoops(zoop_i)%INC_zoo
      dop_excr = grazing_p - pop_excr - delta_C * self%zoops(zoop_i)%IPC_zoo
      !If nutrients are limiting then must excrete doc to maintain balance
      IF ((don_excr < 0.0_rk) .AND. (dop_excr < 0.0_rk)) THEN
         !Determine which nutrient is more limiting
         IF ((self%zoops(zoop_i)%INC_zoo * (grazing_n - pon_excr) - delta_C) .GT. &
            (self%zoops(zoop_i)%IPC_zoo * (grazing_p - pop_excr) - delta_C)) THEN
             don_excr = 0.0_rk
             doc_excr =  (grazing_n - pon_excr) / self%zoops(zoop_i)%INC_zoo - delta_C
             delta_C = delta_C - doc_excr
             dop_excr = grazing_p - pop_excr - delta_C*self%zoops(zoop_i)%IPC_zoo
         ELSE
             dop_excr = 0.0_rk
             doc_excr = (grazing_p - pop_excr) / self%zoops(zoop_i)%IPC_zoo - delta_C
             delta_C = delta_C - doc_excr
             don_excr = grazing_n - pon_excr - delta_C*self%zoops(zoop_i)%INC_zoo
         ENDIF
      ELSEIF (don_excr < 0.0_rk) THEN !nitrogen limited
         don_excr = 0.0_rk
         doc_excr = (grazing_n - pon_excr) / self%zoops(zoop_i)%INC_zoo - delta_C
         delta_C = delta_C - doc_excr
         dop_excr = grazing_p - pop_excr - delta_C*self%zoops(zoop_i)%IPC_zoo
      ELSEIF (dop_excr < 0.0_rk) THEN !phosphorus limited
         dop_excr = 0.0_rk
         doc_excr = (grazing_p - pop_excr) / self%zoops(zoop_i)%IPC_zoo - delta_C
         delta_C = delta_C - doc_excr
         don_excr = grazing_n - pon_excr - delta_C*self%zoops(zoop_i)%INC_zoo
      ELSE !just excrete nutrients no need to balance c
          doc_excr = 0.0_rk
      ENDIF


      !write(*,"(4X,'limitations (f_T,f_Salinity): ',2F8.2)")f_T,f_Salinity
      !write(*,"(4X,'sources/sinks (grazing,respiration,mortaility): ',3F8.2)")grazing,excretion,mortality


      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

      ! Zooplankton production / losses in mmolC/s

      _SET_ODE_(self%id_zoo(zoop_i), (self%zoops(zoop_i)%fassim_zoo * grazing - respiration - mortality)*zoo )


      ! Now take food grazed by zooplankton from food pools in mmolC/s
      phy_i = 0
      DO prey_i = 1,self%zoops(zoop_i)%num_prey
         _SET_ODE_(self%zoops(zoop_i)%id_prey(prey_i), -1.0 * grazing_prey(prey_i))
          IF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey .EQ. 'aed_organic_matter_poc') THEN
              IF (poc > 0.0_rk) THEN
                 _SET_ODE_(self%id_Nmorttarget, -1.0 * grazing_prey(prey_i) * pon/poc)
                 _SET_ODE_(self%id_Pmorttarget, -1.0 * grazing_prey(prey_i) * pop/poc)
              ENDIF
          ELSEIF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:17).EQ.'aed_phytoplankton') THEN
            phy_i = phy_i + 1
            _SET_ODE_(self%zoops(zoop_i)%id_phyIN(phy_i), -1.0 * grazing_prey(prey_i) / prey(prey_i) * phy_INcon(phy_i))
            _SET_ODE_(self%zoops(zoop_i)%id_phyIP(phy_i), -1.0 * grazing_prey(prey_i) / prey(prey_i) * phy_IPcon(phy_i))
         ENDIF
      ENDDO


      ! Now manage excretion contributions to DOM
      IF (self%simDCexcr) THEN
         _SET_ODE_(self%id_Cexctarget,self%zoops(zoop_i)%fexcr_zoo * respiration * zoo + doc_excr)
      ENDIF
      IF (self%simDNexcr) THEN
         _SET_ODE_(self%id_Nexctarget,don_excr)
      ENDIF
      IF (self%simDPexcr) THEN
         _SET_ODE_(self%id_Pexctarget,dop_excr)
      ENDIF

      ! Now manage messy feeding, fecal pellets and mortality contributions to POM
      IF (self%simPCexcr) THEN
         _SET_ODE_(self%id_Cmorttarget, poc_excr)
      ENDIF
      IF (self%simPNexcr) THEN
         _SET_ODE_(self%id_Nmorttarget, pon_excr)
      ENDIF
      IF (self%simPPexcr) THEN
         _SET_ODE_(self%id_Pmorttarget, pop_excr)
      ENDIF

      ! Export diagnostic variables
      _SET_DIAGNOSTIC_(self%id_grz ,grazing*secs_pr_day)
      _SET_DIAGNOSTIC_(self%id_resp ,respiration*secs_pr_day)
      _SET_DIAGNOSTIC_(self%id_mort ,mortality*secs_pr_day)

   ENDDO

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
END SUBROUTINE aed_zooplankton_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE aed_zooplankton_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of zooplankton biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_zooplankton),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   real(rk)           :: zoo,temp,salinity !State variables
   real(rk), ALLOCATABLE,DIMENSION(:)  :: prey !Prey state variables
   real(rk)           :: dn_excr, dp_excr, dc_excr !Excretion state variables
   real(rk)           :: pon, pop, poc !Mortaility and fecal pellet state variables
   real(rk)           :: FGrazing_Limitation, f_T, f_Salinity
   real(rk)           :: Ctotal_prey !total concentration of available prey
   real(rk)           :: grazing, respiration, mortality !Growth & decay functions
   INTEGER            :: zoop_i,prey_i
   real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current environmental conditions.
   _GET_(self%id_tem,temp)     ! local temperature
   _GET_(self%id_sal,salinity) ! local salinity

   ! Retrieve current (local) state variable values.
   IF (self%simDNexcr) _GET_(self%id_Nexctarget, dn_excr)
   IF (self%simDPexcr) _GET_(self%id_Pexctarget, dp_excr)
   IF (self%simDCexcr) _GET_(self%id_Cexctarget, dc_excr)

   IF (self%simPNexcr) _GET_(self%id_Nmorttarget, pon)
   IF (self%simPPexcr) _GET_(self%id_Pmorttarget, pop)
   IF (self%simPCexcr) _GET_(self%id_Cmorttarget, poc)

   DO zoop_i=1,self%num_zoops

      ! Retrieve this zooplankton group
      _GET_(self%id_zoo(zoop_i),zoo)

      !Retrieve prey groups
      Ctotal_prey   = 0.0_rk
      DO prey_i=1,self%zoops(zoop_i)%num_prey
         _GET_(self%zoops(zoop_i)%id_prey(prey_i),prey(prey_i))
         Ctotal_prey = Ctotal_prey + prey(prey_i)
      ENDDO

      grazing       = 0.0_rk
      respiration   = 0.0_rk
      mortality     = 0.0_rk

      ! Get the grazing limitation function
       fGrazing_Limitation = fPrey_Limitation(self,zoop_i,Ctotal_prey)

      ! Get the temperature function
       f_T = fTemp_function(1, self%zoops(zoop_i)%Tmax_zoo,       &
                               self%zoops(zoop_i)%Tstd_zoo,       &
                               self%zoops(zoop_i)%theta_resp_zoo, &
                               self%zoops(zoop_i)%aTn,            &
                               self%zoops(zoop_i)%bTn,            &
                               self%zoops(zoop_i)%kTn,temp)

      ! Get the salinity limitation.
       f_Salinity = fSalinity_Limitation(self,zoop_i,salinity)

      ! Get the growth rate (/ s)
      grazing = self%zoops(zoop_i)%Rgrz_zoo * fGrazing_Limitation * f_T

      ! Get the respiration rate (/ s)
      respiration = self%zoops(zoop_i)%Rresp_zoo * f_Salinity

      ! Get the mortality rate (/ s)
      mortality = self%zoops(zoop_i)%Rmort_zoo * f_T

      ! Don't excrete or die if we are at the min biomass otherwise we have a
      ! leak in the C mass balance
      IF (zoo <= self%zoops(zoop_i)%min_zoo) THEN
        respiration = 0.0
        mortality = 0.0
      ENDIF

      ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER

#if 0
      ! Zooplankton production / losses
      _SET_PP_(self%id_zoo(zoop_i),self%id_zoo(zoop_i), (self%zoops(zoop_i)%fassim_zoo * grazing - respiration - mortality)*zoo )


      ! Now manage excretion contributions to DOM
      IF (self%simDNexcr.gt.0) THEN
         _SET_PP_(self%id_Nexctarget,(self%id_Nexctarget,self%zoops(zoop_i)%fexcr_zoo * respiration * self%zoops(zoop_i)%INC_zoo * zoo)
      ENDIF
      IF (self%simDPexcr.gt.0) THEN
         _SET_PP_(self%id_Pexctarget,self%id_Pexctarget,self%zoops(zoop_i)%fexcr_zoo * respiration * self%zoops(zoop_i)%IPC_zoo * zoo)
      ENDIF
      IF (self%simDCexcr.gt.0) THEN
         _SET_PP_(self%id_Cexctarget,self%id_Cexctarget,self%zoops(zoop_i)%fexcr_zoo * respiration * zoo)
      ENDIF

      ! Now manage messy feeding, fecal pellets and mortality contributions to POM
      !IF (self%simPOM .EQ. 1) THEN
         !_SET_ODE_(self%id_Nmorttarget,self%ffecal_zoo * respiration * self%INC_con * zoo)
         !_SET_ODE_(self%id_Nmorttarget,self%ffecal_zoo * respiration * self%INP_con * zoo)
         !_SET_ODE_(self%id_Nmortarget,((1 - self%fassim_zoo) * grazing + self%ffecal_zoo * respiration + mortaility) * zoo
      !ENDIF

      ! Export diagnostic variables
      _SET_DIAGNOSTIC_(self%id_grz ,grazing*secs_pr_day)
      _SET_DIAGNOSTIC_(self%id_resp ,respiration*secs_pr_day)
      _SET_DIAGNOSTIC_(self%id_mort ,mortality*secs_pr_day)

#endif
   ENDDO

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_
END SUBROUTINE aed_zooplankton_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE aed_zooplankton_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_aed_zooplankton),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   real(rk)  :: zoo
   real(rk)  :: Total_zoo, TN_zoo, TP_zoo
   INTEGER   :: zoo_i
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   Total_zoo = 0.0_rk
   TN_zoo = 0.0_rk
   TP_zoo = 0.0_rk
   DO zoo_i=1,self%num_zoops
      ! Retrieve current (local) state variable values.
      _GET_(self%id_zoo(zoo_i),zoo) ! zooplankton
      Total_zoo = Total_zoo + zoo
      TN_zoo = TN_zoo + zoo * self%zoops(zoo_i)%INC_zoo
      TP_zoo = TP_zoo + zoo * self%zoops(zoo_i)%IPC_zoo
   ENDDO

   ! Total zooplankton is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totZOO,Total_zoo)

   ! Total nutrient is simply the sum of all internal zooplankton concentrations.
   _SET_CONSERVED_QUANTITY_(self%id_totN,TN_zoo)
   _SET_CONSERVED_QUANTITY_(self%id_totP,TP_zoo)
   _SET_CONSERVED_QUANTITY_(self%id_totC,Total_zoo)


   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_zooplankton_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
FUNCTION fPrey_Limitation(self,group,C) RESULT(fPlim)
!----------------------------------------------------------------------------!
! Michaelis-Menten formulation for zooplankton grazing on available
! prey is applied.
!----------------------------------------------------------------------------!
   !-- Incoming
   CLASS (type_aed_zooplankton),INTENT(in) :: self
   INTEGER                                 :: group
   real(rk),INTENT(in)                     :: C !total concentration of available prey
!
!LOCALS
   ! Returns the M-M limitation function
   real(rk)                                 :: fPlim
!
!-------------------------------------------------------------------------------
!BEGIN
   fPlim = 1.0

   fPlim = C/(self%zoops(group)%Kgrz_zoo+C)

   IF( fPlim<0.0_rk ) fPlim=0.0_rk

 END FUNCTION fPrey_Limitation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION fSalinity_Limitation(self,group,S) RESULT(fSal)
!----------------------------------------------------------------------------!
! Salinity tolerance of zooplankton                                          !
!----------------------------------------------------------------------------!
!ARGUMENTS
   CLASS (type_aed_zooplankton),INTENT(in) :: self
   INTEGER                                 :: group
   real(rk),INTENT(in)                     :: S
!
!LOCALS
   real(rk)  :: fSal ! Returns the salinity function
   real(rk)  :: Smin,Smax,Sint
!
!-------------------------------------------------------------------------------
!BEGIN
   Smin = self%zoops(group)%Smin_zoo
   Smax = self%zoops(group)%Smax_zoo
   Sint = self%zoops(group)%Sint_zoo

   !Salinity factor represents natural mortality in response to salinity stress.
   ! f(S) = 1 at Smin<=S<=Smax, f(S) = Bep at S=0 & S=2*Smax.
   IF (S < Smin) THEN
      fSal = (Sint-1.0)/(Smin**2.0)*(S**2.0) - 2*(Sint-1.0)/Smin*S + Sint
   ELSEIF(S > Smax) THEN
      fSal = (Sint-1.0)/(Smax**2.0)*(S**2.0) - 2*(Sint-1.0)/Smax*S + Sint
   ELSE
      fSal = 1.0
   ENDIF

   IF( fSal<0.0_rk ) fSal=0.0_rk
END FUNCTION fSalinity_Limitation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_zooplankton

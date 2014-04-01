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

#define _PHYLEN_ 17
#define _PHYMOD_ 'aed_phytoplankton'
#define _OGMPOC_ 'aed_organic_matter_poc'

MODULE aed_zooplankton
!-------------------------------------------------------------------------------
!  aed_zooplankton --- multi zooplankton biogeochemical model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util,ONLY : find_free_lun,aed_bio_temp_function, fTemp_function,qsort
   USE aed_zoop_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_type_zooplankton
!

   TYPE,extends(type_base_model) :: aed_type_zooplankton
!     Variable identifiers
      TYPE (type_state_variable_id)      :: id_zoo(MAX_ZOOP_TYPES)
      TYPE (type_state_variable_id)      :: id_Nexctarget,id_Nmorttarget
      TYPE (type_state_variable_id)      :: id_Pexctarget,id_Pmorttarget
      TYPE (type_state_variable_id)      :: id_Cexctarget,id_Cmorttarget
      TYPE (type_state_variable_id)      :: id_DOupttarget
      TYPE (type_dependency_id)          :: id_tem, id_sal, id_extc
      TYPE (type_diagnostic_variable_id) :: id_grz,id_resp,id_mort


!     Model parameters
      INTEGER                                   :: num_zoops
      TYPE(type_zoop_data),DIMENSION(:),ALLOCATABLE :: zoops
      LOGICAL  :: simDNexcr, simDPexcr, simDCexcr
      LOGICAL  :: simPNexcr, simPPexcr, simPCexcr

      CONTAINS     ! Model Methods
        PROCEDURE :: initialize               => aed_init_zooplankton
        PROCEDURE :: do                       => aed_zooplankton_do
        PROCEDURE :: do_ppdd                  => aed_zooplankton_do_ppdd
!       PROCEDURE :: do_benthos               => aed_zooplankton_do_benthos
   END TYPE

   LOGICAL :: debug = .TRUE.

CONTAINS
!===============================================================================


!###############################################################################
SUBROUTINE aed_zooplankton_load_params(self, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_zooplankton),INTENT(inout) :: self
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count !Number of zooplankton groups
   INTEGER,INTENT(in)          :: list(*) !List of zooplankton groups to simulate
!
!LOCALS
   INTEGER  :: status

   INTEGER  :: i,j,tfil,sort_i(MAX_ZOOP_PREY)
   AED_REAL :: Pzoo_prey(MAX_ZOOP_PREY)

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   TYPE(type_zoop_params)  :: zoop_param(MAX_ZOOP_TYPES)
   NAMELIST /zoop_params/ zoop_param
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD',iostat=status)
    IF (status /= 0) STOP 'Error opening zoop_params namelist file'
    read(tfil,nml=zoop_params,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist zoop_params'

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
!
END SUBROUTINE aed_zooplankton_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_init_zooplankton(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the zooplankton biogeochemical model
!
!  Here, the aed_zooplankton namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_zooplankton),TARGET,INTENT(inout) :: self
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status

   INTEGER            :: num_zoops
   INTEGER            :: the_zoops(MAX_ZOOP_TYPES)

   CHARACTER(len=64)  :: dn_target_variable='' !dissolved nitrogen target variable
   CHARACTER(len=64)  :: pn_target_variable='' !particulate nitrogen target variable
   CHARACTER(len=64)  :: dp_target_variable='' !dissolved phosphorus target variable
   CHARACTER(len=64)  :: pp_target_variable='' !particulate phosphorus target variable
   CHARACTER(len=64)  :: dc_target_variable='' !dissolved carbon target variable
   CHARACTER(len=64)  :: pc_target_variable='' !particulate carbon target variable
   CHARACTER(len=128) :: dbase='aed_zoop_pars.nml'

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   INTEGER            :: zoop_i, prey_i, phy_i

   NAMELIST /aed_zooplankton/ num_zoops, the_zoops, &
                    dn_target_variable, pn_target_variable, dp_target_variable, &
                    pp_target_variable, dc_target_variable, pc_target_variable, &
                    dbase
!-----------------------------------------------------------------------
!BEGIN
!print *,'**** Reading /aed_zooplankton/ namelist'
   ! Read the namelist
   read(namlst,nml=aed_zooplankton,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_zooplankton'

    self%num_zoops = 0
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted in aed_zooplankton_load_params to values per second.
   CALL aed_zooplankton_load_params(self, dbase, num_zoops, the_zoops)


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
          IF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
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
   CALL self%register_diagnostic_variable(self%id_grz,'grz','mmolC/m**3',  'net zooplankton grazing')
   CALL self%register_diagnostic_variable(self%id_resp,'resp','mmolC/m**3',  'net zooplankton respiration')
   CALL self%register_diagnostic_variable(self%id_mort,'mort','mmolC/m**3/d','net zooplankton mortality')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_tem,standard_variables%temperature)
   CALL self%register_dependency(self%id_sal,standard_variables%practical_salinity)
   CALL self%register_dependency(self%id_extc,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
END SUBROUTINE aed_init_zooplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_zooplankton_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of zooplankton biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_zooplankton),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   AED_REAL           :: zoo,temp,salinity !State variables
   AED_REAL           :: prey(MAX_ZOOP_PREY), grazing_prey(MAX_ZOOP_PREY) !Prey state variables
   AED_REAL           :: phy_INcon(MAX_ZOOP_PREY), phy_IPcon(MAX_ZOOP_PREY) !Internal nutrients for phytoplankton
   AED_REAL           :: dn_excr, dp_excr, dc_excr !Excretion state variables
   AED_REAL           :: pon, pop, poc !Mortaility and fecal pellet state variables
   AED_REAL           :: FGrazing_Limitation, f_T, f_Salinity
   AED_REAL           :: pref_factor, Ctotal_prey !total concentration of available prey
   AED_REAL           :: food, grazing, respiration, mortality !Growth & decay functions
   AED_REAL           :: grazing_n, grazing_p !Grazing on nutrients
   AED_REAL           :: pon_excr, pop_excr, poc_excr !POM excretion rates
   AED_REAL           :: don_excr, dop_excr, doc_excr, delta_C !DOM excretion rates
   INTEGER            :: zoop_i,prey_i,prey_j,phy_i
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

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
      Ctotal_prey   = zero_
      DO prey_i=1,self%zoops(zoop_i)%num_prey
         _GET_(self%zoops(zoop_i)%id_prey(prey_i),prey(prey_i))
         Ctotal_prey = Ctotal_prey + prey(prey_i)
      ENDDO

      grazing       = zero_
      respiration   = zero_
      mortality     = zero_

      ! Get the grazing limitation function
       fGrazing_Limitation = fPrey_Limitation(self%zoops,zoop_i,Ctotal_prey)

      ! Get the temperature function
       f_T = fTemp_function(1, self%zoops(zoop_i)%Tmax_zoo,       &
                               self%zoops(zoop_i)%Tstd_zoo,       &
                               self%zoops(zoop_i)%theta_resp_zoo, &
                               self%zoops(zoop_i)%aTn,            &
                               self%zoops(zoop_i)%bTn,            &
                               self%zoops(zoop_i)%kTn, temp)

      ! Get the salinity limitation.
       f_Salinity = fSalinity_Limitation(self%zoops,zoop_i,salinity)

      ! Get the growth rate (/ s)
      ! grazing is in units of mass consumed/mass zoops/unit time
      grazing = self%zoops(zoop_i)%Rgrz_zoo * fGrazing_Limitation * f_T

      ! Now dertermine available prey and limit grazing amount to
      ! availability of prey
      ! food is total amount of food in units of mass/unit volume/unit time
      food = grazing * zoo
      IF (Ctotal_prey < self%zoops(zoop_i)%num_prey * self%zoops(zoop_i)%Cmin_grz_zoo ) THEN
          food = zero_
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
          pref_factor = zero_
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
             grazing_prey(prey_i) = zero_
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

      grazing_n = zero_
      grazing_p = zero_
      phy_i = 0
      DO prey_i = 1,self%zoops(zoop_i)%num_prey
         IF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey .EQ. _OGMPOC_) THEN
            IF (poc > zero_) THEN
                grazing_n = grazing_n + grazing_prey(prey_i) * pon/poc
                grazing_p = grazing_p + grazing_prey(prey_i) * pop/poc
            ELSE
                grazing_n = zero_
                grazing_p = zero_
            ENDIF
         ELSEIF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
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
        respiration = zero_
        mortality = zero_
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
      IF ((don_excr < zero_) .AND. (dop_excr < zero_)) THEN
         !Determine which nutrient is more limiting
         IF ((self%zoops(zoop_i)%INC_zoo * (grazing_n - pon_excr) - delta_C) .GT. &
            (self%zoops(zoop_i)%IPC_zoo * (grazing_p - pop_excr) - delta_C)) THEN
             don_excr = zero_
             doc_excr =  (grazing_n - pon_excr) / self%zoops(zoop_i)%INC_zoo - delta_C
             delta_C = delta_C - doc_excr
             dop_excr = grazing_p - pop_excr - delta_C*self%zoops(zoop_i)%IPC_zoo
         ELSE
             dop_excr = zero_
             doc_excr = (grazing_p - pop_excr) / self%zoops(zoop_i)%IPC_zoo - delta_C
             delta_C = delta_C - doc_excr
             don_excr = grazing_n - pon_excr - delta_C*self%zoops(zoop_i)%INC_zoo
         ENDIF
      ELSEIF (don_excr < zero_) THEN !nitrogen limited
         don_excr = zero_
         doc_excr = (grazing_n - pon_excr) / self%zoops(zoop_i)%INC_zoo - delta_C
         delta_C = delta_C - doc_excr
         dop_excr = grazing_p - pop_excr - delta_C*self%zoops(zoop_i)%IPC_zoo
      ELSEIF (dop_excr < zero_) THEN !phosphorus limited
         dop_excr = zero_
         doc_excr = (grazing_p - pop_excr) / self%zoops(zoop_i)%IPC_zoo - delta_C
         delta_C = delta_C - doc_excr
         don_excr = grazing_n - pon_excr - delta_C*self%zoops(zoop_i)%INC_zoo
      ELSE !just excrete nutrients no need to balance c
          doc_excr = zero_
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
          IF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey .EQ. _OGMPOC_) THEN
              IF (poc > zero_) THEN
                 _SET_ODE_(self%id_Nmorttarget, -1.0 * grazing_prey(prey_i) * pon/poc)
                 _SET_ODE_(self%id_Pmorttarget, -1.0 * grazing_prey(prey_i) * pop/poc)
              ENDIF
          ELSEIF (self%zoops(zoop_i)%prey(prey_i)%zoop_prey(1:_PHYLEN_).EQ. _PHYMOD_) THEN
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
   _LOOP_END_
END SUBROUTINE aed_zooplankton_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_zooplankton_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of zooplankton biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_zooplankton),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL           :: zoo,temp,salinity !State variables
   AED_REAL, ALLOCATABLE,DIMENSION(:)  :: prey !Prey state variables
   AED_REAL           :: dn_excr, dp_excr, dc_excr !Excretion state variables
   AED_REAL           :: pon, pop, poc !Mortaility and fecal pellet state variables
   AED_REAL           :: FGrazing_Limitation, f_T, f_Salinity
   AED_REAL           :: Ctotal_prey !total concentration of available prey
   AED_REAL           :: grazing, respiration, mortality !Growth & decay functions
   INTEGER            :: zoop_i,prey_i
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

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
      Ctotal_prey   = zero_
      DO prey_i=1,self%zoops(zoop_i)%num_prey
         _GET_(self%zoops(zoop_i)%id_prey(prey_i),prey(prey_i))
         Ctotal_prey = Ctotal_prey + prey(prey_i)
      ENDDO

      grazing       = zero_
      respiration   = zero_
      mortality     = zero_

      ! Get the grazing limitation function
       fGrazing_Limitation = fPrey_Limitation(self%zoops,zoop_i,Ctotal_prey)

      ! Get the temperature function
       f_T = fTemp_function(1, self%zoops(zoop_i)%Tmax_zoo,       &
                               self%zoops(zoop_i)%Tstd_zoo,       &
                               self%zoops(zoop_i)%theta_resp_zoo, &
                               self%zoops(zoop_i)%aTn,            &
                               self%zoops(zoop_i)%bTn,            &
                               self%zoops(zoop_i)%kTn,temp)

      ! Get the salinity limitation.
       f_Salinity = fSalinity_Limitation(self%zoops,zoop_i,salinity)

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
   _LOOP_END_
END SUBROUTINE aed_zooplankton_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_zooplankton

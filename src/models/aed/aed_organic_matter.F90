!###############################################################################
!#                                                                             #
!# aed_organic_matter.F90                                                      #
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
!# Created 6 June 2011                                                         #
!#                                                                             #
!###############################################################################

#include "aed.h"

MODULE aed_organic_matter
!-------------------------------------------------------------------------------
! aed_organic_matter --- organic matter biogeochemical model
!
! Organic Matter module contains equations for mineralisation
! of particulate and dissolved organic matter.
! In future to include adsorption/desorption to suspended solids.
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE  ! By default make everything private
!
   PUBLIC aed_type_organic_matter
!
   TYPE,extends(type_base_model) :: aed_type_organic_matter
!     Variable identifiers
      TYPE (type_state_variable_id)      :: id_pon,id_don ! particulate and dissolved organic nitrogen
      TYPE (type_state_variable_id)      :: id_pop,id_dop ! particulate and dissolved organic phosphorus
      TYPE (type_state_variable_id)      :: id_poc,id_doc ! particulate and dissolved organic carbon
      TYPE (type_state_variable_id)      :: id_oxy,id_amm,id_frp,id_dic
      TYPE (type_horizontal_dependency_id)      :: id_Fsed_pon, id_Fsed_don ! sed. rate organic nitrogen
      TYPE (type_horizontal_dependency_id)      :: id_Fsed_pop, id_Fsed_dop ! sed. rate organic phosphorus
      TYPE (type_horizontal_dependency_id)      :: id_Fsed_poc, id_Fsed_doc ! sed. rate organic carbon
      TYPE (type_horizontal_dependency_id)      :: id_Psed_poc, id_Psed_pon, id_Psed_pop ! sedimentation rates
      TYPE (type_dependency_id)          :: id_temp
      TYPE (type_diagnostic_variable_id) :: id_pon_miner, id_don_miner
      TYPE (type_diagnostic_variable_id) :: id_pop_miner, id_dop_miner
      TYPE (type_diagnostic_variable_id) :: id_poc_miner, id_doc_miner
      TYPE (type_horizontal_diagnostic_variable_id) :: id_sed_pon, id_sed_don
      TYPE (type_horizontal_diagnostic_variable_id) :: id_sed_pop, id_sed_dop
      TYPE (type_horizontal_diagnostic_variable_id) :: id_sed_poc, id_sed_doc
      TYPE (type_diagnostic_variable_id) :: id_bod

      !# Model parameters
      AED_REAL :: w_pon,Rpon_miner,Rdon_miner,Fsed_pon,Fsed_don, &
                          Kpon_miner, Kdon_miner, Ksed_don, &
                          theta_pon_miner, theta_don_miner, theta_sed_don
      AED_REAL :: w_pop,Rpop_miner,Rdop_miner,Fsed_pop,Fsed_dop, &
                          Kpop_miner, Kdop_miner, Ksed_dop, &
                          theta_pop_miner, theta_dop_miner, theta_sed_dop
      AED_REAL :: w_poc,Rpoc_miner,Rdoc_miner,Fsed_poc,Fsed_doc, &
                          Kpoc_miner, Kdoc_miner, Ksed_doc, &
                          theta_poc_miner, theta_doc_miner, theta_sed_doc, &
                          KeDOM, KePOM
      LOGICAL  :: use_oxy, use_amm, use_frp, use_dic, use_sed_model, use_sedmtn_model

      CONTAINS    ! Model Procedures
        PROCEDURE :: initialize               => aed_init_organic_matter
        PROCEDURE :: do                       => aed_organic_matter_do
        PROCEDURE :: do_ppdd                  => aed_organic_matter_do_ppdd
        PROCEDURE :: do_benthos               => aed_organic_matter_do_benthos
        PROCEDURE :: get_light_extinction     => aed_organic_matter_get_light_extinction
   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_init_organic_matter(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                              :: namlst
   CLASS (aed_type_organic_matter),TARGET,INTENT(inout) :: self
!
!LOCALS
   INTEGER  :: status

   AED_REAL                  :: pon_initial = 4.5
   AED_REAL                  :: don_initial = 4.5
   AED_REAL                  :: w_pon       = 0.0
   AED_REAL                  :: Rpon_miner  = 0.01
   AED_REAL                  :: Rdon_miner  = 0.01
   AED_REAL                  :: Fsed_pon    =  0.0
   AED_REAL                  :: Fsed_don    = 30.0
   AED_REAL                  :: Kpon_miner  = 30.0
   AED_REAL                  :: Kdon_miner  = 30.0
   AED_REAL                  :: Ksed_don    = 4.5
   AED_REAL                  :: theta_pon_miner = 1.0
   AED_REAL                  :: theta_don_miner = 1.0
   AED_REAL                  :: theta_sed_don   = 1.0
   CHARACTER(len=64)         :: don_miner_product_variable=''
   CHARACTER(len=64)         :: Fsed_pon_variable=''
   CHARACTER(len=64)         :: Fsed_don_variable=''

   AED_REAL                  :: pop_initial = 4.5
   AED_REAL                  :: dop_initial = 4.5
   AED_REAL                  :: w_pop       = 0.0
   AED_REAL                  :: Rpop_miner  = 0.01
   AED_REAL                  :: Rdop_miner  = 0.01
   AED_REAL                  :: Fsed_pop    =  0.0
   AED_REAL                  :: Fsed_dop    = 30.0
   AED_REAL                  :: Kpop_miner  = 30.0
   AED_REAL                  :: Kdop_miner  = 30.0
   AED_REAL                  :: Ksed_dop    = 4.5
   AED_REAL                  :: theta_pop_miner = 1.0
   AED_REAL                  :: theta_dop_miner = 1.0
   AED_REAL                  :: theta_sed_dop   = 1.0
   CHARACTER(len=64)         :: dop_miner_product_variable=''
   CHARACTER(len=64)         :: Fsed_pop_variable=''
   CHARACTER(len=64)         :: Fsed_dop_variable=''

   AED_REAL                  :: poc_initial = 4.5
   AED_REAL                  :: doc_initial = 4.5
   AED_REAL                  :: w_poc       = 0.0
   AED_REAL                  :: Rpoc_miner  = 0.01
   AED_REAL                  :: Rdoc_miner  = 0.01
   AED_REAL                  :: Fsed_poc    =  0.0
   AED_REAL                  :: Fsed_doc    = 30.0
   AED_REAL                  :: Kpoc_miner  = 30.0
   AED_REAL                  :: Kdoc_miner  = 30.0
   AED_REAL                  :: Ksed_doc    = 4.5
   AED_REAL                  :: theta_poc_miner = 1.0
   AED_REAL                  :: theta_doc_miner = 1.0
   AED_REAL                  :: theta_sed_doc   = 1.0
   AED_REAL                  :: KeDOM = 0.01
   AED_REAL                  :: KePOM = 0.01
   CHARACTER(len=64)         :: doc_miner_product_variable=''
   CHARACTER(len=64)         :: doc_miner_reactant_variable=''
   CHARACTER(len=64)         :: Fsed_poc_variable=''
   CHARACTER(len=64)         :: Fsed_doc_variable=''

   CHARACTER(len=64)         :: Psed_poc_variable=''
   CHARACTER(len=64)         :: Psed_pon_variable=''
   CHARACTER(len=64)         :: Psed_pop_variable=''


   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_organic_matter/ &
             pon_initial, don_initial, w_pon, Rpon_miner, Rdon_miner, Fsed_pon, Fsed_don, &
             Kpon_miner, Kdon_miner, Ksed_don,                &
             theta_pon_miner, theta_don_miner, theta_sed_don, &
             don_miner_product_variable,                      &
             Fsed_pon_variable, Fsed_don_variable,            &
             pop_initial, dop_initial, w_pop, Rpop_miner, Rdop_miner, Fsed_pop, Fsed_dop, &
             Kpop_miner, Kdop_miner, Ksed_dop,                &
             theta_pop_miner, theta_dop_miner, theta_sed_dop, &
             dop_miner_product_variable,                      &
             Fsed_pop_variable, Fsed_dop_variable,            &
             poc_initial, doc_initial, w_poc, Rpoc_miner, Rdoc_miner, Fsed_poc, Fsed_doc, &
             Kpoc_miner, Kdoc_miner, Ksed_doc,                &
             theta_poc_miner, theta_doc_miner, theta_sed_doc, KeDOM, KePOM, &
             doc_miner_reactant_variable, doc_miner_product_variable, &
             Fsed_poc_variable, Fsed_doc_variable, &
             Psed_poc_variable, Psed_pon_variable, Psed_pop_variable

!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_organic_matter,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_organic_matter'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%w_pon       = w_pon/secs_pr_day
   self%Rpon_miner  = Rpon_miner/secs_pr_day
   self%Rdon_miner  = Rdon_miner/secs_pr_day
 ! self%Fsed_pon    = Fsed_pon/secs_pr_day
   self%Fsed_pon    = 0.0
   self%Fsed_don    = Fsed_don/secs_pr_day
   self%Kpon_miner  = Kpon_miner
   self%Kdon_miner  = Kdon_miner
   self%Ksed_don  = Ksed_don
   self%theta_pon_miner = theta_pon_miner
   self%theta_don_miner = theta_don_miner
   self%theta_sed_don = theta_sed_don

   self%w_pop       = w_pop/secs_pr_day
   self%Rpop_miner  = Rpop_miner/secs_pr_day
   self%Rdop_miner  = Rdop_miner/secs_pr_day
 ! self%Fsed_pop    = Fsed_pop/secs_pr_day
   self%Fsed_pop    = 0.0
   self%Fsed_dop    = Fsed_dop/secs_pr_day
   self%Kpop_miner  = Kpop_miner
   self%Kdop_miner  = Kdop_miner
   self%Ksed_dop  = Ksed_dop
   self%theta_pop_miner = theta_pop_miner
   self%theta_dop_miner = theta_dop_miner
   self%theta_sed_dop = theta_sed_dop

   self%w_poc       = w_poc/secs_pr_day
   self%Rpoc_miner  = Rpoc_miner/secs_pr_day
   self%Rdoc_miner  = Rdoc_miner/secs_pr_day
 ! self%Fsed_poc    = Fsed_poc/secs_pr_day
   self%Fsed_poc    = 0.0
   self%Fsed_doc    = Fsed_doc/secs_pr_day
   self%Kpoc_miner  = Kpoc_miner
   self%Kdoc_miner  = Kdoc_miner
   self%Ksed_doc  = Ksed_doc
   self%theta_poc_miner = theta_poc_miner
   self%theta_doc_miner = theta_doc_miner
   self%theta_sed_doc = theta_sed_doc
   self%KeDOM       = KeDOM
   self%KePOM       = KePOM

   ! Register state variables
   CALL self%register_state_variable(self%id_don,'don','mmol/m**3','dissolved organic nitrogen',     &
                                    don_initial,minimum=zero_,no_river_dilution=.true.)

   CALL self%register_state_variable(self%id_pon,'pon','mmol/m**3','particulate organic nitrogen',   &
                                    pon_initial,minimum=zero_,no_river_dilution=.true.,vertical_movement=self%w_pon)

   CALL self%register_state_variable(self%id_dop,'dop','mmol/m**3','dissolved organic phosphorus',   &
                                    dop_initial,minimum=zero_,no_river_dilution=.true.)
   CALL self%register_state_variable(self%id_pop,'pop','mmol/m**3','particulate organic phosphorus', &
                                    pop_initial,minimum=zero_,no_river_dilution=.true.,vertical_movement=self%w_pop)

   CALL self%register_state_variable(self%id_doc,'doc','mmol/m**3','dissolved organic carbon',       &
                                    doc_initial,minimum=zero_,no_river_dilution=.true.)
   CALL self%register_state_variable(self%id_poc,'poc','mmol/m**3','particulate organic carbon',     &
                                    poc_initial,minimum=zero_,no_river_dilution=.true.,vertical_movement=self%w_poc)

   ! Register external state variable dependencies (carbon)
   self%use_oxy = doc_miner_reactant_variable .NE. '' !This means oxygen module switched on
   IF (self%use_oxy) THEN
     CALL self%register_state_dependency(self%id_oxy,doc_miner_reactant_variable)
   ENDIF

   self%use_dic = doc_miner_product_variable .NE. '' !This means carbon module switched on
   IF (self%use_dic) THEN
     CALL self%register_state_dependency(self%id_dic,doc_miner_product_variable)
   ENDIF

   ! Register external state variable dependencies (nitrogen)
   self%use_amm = don_miner_product_variable .NE. '' !This means nitrogen module switched on
   IF (self%use_amm) THEN
     CALL self%register_state_dependency(self%id_amm,don_miner_product_variable)
   ENDIF

   ! Register external state variable dependencies (phosphorous)
   self%use_frp = dop_miner_product_variable .NE. '' !This means phosphorus module switched on
   IF (self%use_frp) THEN
     CALL self%register_state_dependency(self%id_frp,dop_miner_product_variable)
   ENDIF

   self%use_sed_model = Fsed_pon_variable .NE. ''
   IF (self%use_sed_model) THEN
     CALL self%register_horizontal_dependency(self%id_Fsed_pon,Fsed_pon_variable)
     CALL self%request_coupling(self%id_Fsed_pon,Fsed_pon_variable)
     CALL self%register_horizontal_dependency(self%id_Fsed_don,Fsed_don_variable)
     CALL self%request_coupling(self%id_Fsed_don,Fsed_don_variable)
     CALL self%register_horizontal_dependency(self%id_Fsed_pop,Fsed_pop_variable)
     CALL self%request_coupling(self%id_Fsed_pop,Fsed_pop_variable)
     CALL self%register_horizontal_dependency(self%id_Fsed_dop,Fsed_dop_variable)
     CALL self%request_coupling(self%id_Fsed_dop,Fsed_dop_variable)
     CALL self%register_horizontal_dependency(self%id_Fsed_poc,Fsed_poc_variable)
     CALL self%request_coupling(self%id_Fsed_poc,Fsed_poc_variable)
     CALL self%register_horizontal_dependency(self%id_Fsed_doc,Fsed_doc_variable)
     CALL self%request_coupling(self%id_Fsed_doc,Fsed_doc_variable)
   ENDIF

   self%use_sedmtn_model = Psed_poc_variable .NE. ''
   IF (self%use_sedmtn_model) THEN
      CALL self%register_horizontal_dependency(self%id_Psed_poc,Psed_poc_variable)
      CALL self%register_horizontal_dependency(self%id_Psed_pon,Psed_pon_variable)
      CALL self%register_horizontal_dependency(self%id_Psed_pop,Psed_pop_variable)
   ENDIF

   ! Register diagnostic variables
   CALL self%register_diagnostic_variable(self%id_pon_miner,'pon_miner','mmol/m**3/d',  'PON mineralisation')
   CALL self%register_diagnostic_variable(self%id_don_miner,'don_miner','mmol/m**3/d',  'DON mineralisation')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_pon,'sed_pon','mmol/m**2/d',  'PON sediment flux')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_don,'sed_don','mmol/m**2/d',  'DON sediment flux')

   CALL self%register_diagnostic_variable(self%id_pop_miner,'pop_miner','mmol/m**3/d',  'POP mineralisation')
   CALL self%register_diagnostic_variable(self%id_dop_miner,'dop_miner','mmol/m**3/d',  'DOP mineralisation')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_pop,'sed_pop','mmol/m**2/d',  'POP sediment flux')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_dop,'sed_dop','mmol/m**2/d',  'DOP sediment flux')

   CALL self%register_diagnostic_variable(self%id_poc_miner,'poc_miner','mmol/m**3/d',  'POC mineralisation')
   CALL self%register_diagnostic_variable(self%id_doc_miner,'doc_miner','mmol/m**3/d',  'DOC mineralisation')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_poc,'sed_poc','mmol/m**2/d',  'POC sediment flux')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_doc,'sed_doc','mmol/m**2/d',  'DOC sediment flux')

   CALL self%register_diagnostic_variable(self%id_bod,'BOD','mmol/m**3',  'Biochemical Oxygen Demand (BOD)')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature)
END SUBROUTINE aed_init_organic_matter
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_organic_matter_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_organic_matter model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_organic_matter),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   AED_REAL                   :: pon,don,amm,oxy,temp !State variables
   AED_REAL                   :: pon_mineralisation, don_mineralisation
   AED_REAL                   :: pop,dop,frp !State variables
   AED_REAL                   :: pop_mineralisation, dop_mineralisation
   AED_REAL                   :: poc,doc,dic !State variables
   AED_REAL                   :: poc_mineralisation, doc_mineralisation
   AED_REAL, parameter        :: secs_pr_day = 86400.
 ! AED_REAL, parameter        :: Yoxy_don_miner = 6.625 !ratio of oxygen to nitrogen utilised during don mineralisation
 ! AED_REAL, parameter        :: Yoxy_dop_miner = 6.625 !ratio of oxygen to phosphoros utilised during dop mineralisation
   AED_REAL, parameter        :: Yoxy_doc_miner = 32./12. !ratio of oxygen to carbon utilised during doc mineralisation

!-----------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
    !call log_message('model aed_organic_matter enter do loop successfully.')

   ! Retrieve current (local) state variable values.
   _GET_(self%id_pon,pon) ! particulate organic nitrogen
   _GET_(self%id_don,don) ! dissolved organic nitrogen
   _GET_(self%id_pop,pop) ! particulate organic phosphorus
   _GET_(self%id_dop,dop) ! dissolved organic phosphorus
   _GET_(self%id_poc,poc) ! particulate organic carbon
   _GET_(self%id_doc,doc) ! dissolved organic carbon


   IF (self%use_oxy) THEN ! & use_oxy
      _GET_(self%id_oxy,oxy) ! oxygen
   ELSE
      oxy = 0.0
   ENDIF
   IF (self%use_dic) THEN ! & use_amm
      _GET_(self%id_dic,dic) ! disolved inorganic carbon
   ELSE
      dic = 0.0
   ENDIF
   IF (self%use_amm) THEN ! & use_amm
      _GET_(self%id_amm,amm) ! ammonium
   ELSE
      amm = 0.0
   ENDIF
   IF (self%use_frp) THEN ! & use_frp
      _GET_(self%id_frp,frp) ! phosphate
   ELSE
      frp = 0.0
   ENDIF

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)  ! temperature

   ! Define some intermediate quantities units mmol N/m3/day
   pon_mineralisation = fpon_miner(self%use_oxy,self%Rpon_miner,self%Kpon_miner,self%theta_pon_miner,oxy,temp)
   don_mineralisation = fdon_miner(self%use_oxy,self%Rdon_miner,self%Kdon_miner,self%theta_don_miner,oxy,temp)
   pop_mineralisation = fpop_miner(self%use_oxy,self%Rpop_miner,self%Kpop_miner,self%theta_pop_miner,oxy,temp)
   dop_mineralisation = fdop_miner(self%use_oxy,self%Rdop_miner,self%Kdop_miner,self%theta_dop_miner,oxy,temp)
   poc_mineralisation = fpoc_miner(self%use_oxy,self%Rpoc_miner,self%Kpoc_miner,self%theta_poc_miner,oxy,temp)
   doc_mineralisation = fdoc_miner(self%use_oxy,self%Rdoc_miner,self%Kdoc_miner,self%theta_doc_miner,oxy,temp)

   ! Set temporal derivatives
   _SET_ODE_(self%id_pon,-pon*pon_mineralisation)
   _SET_ODE_(self%id_don,pon*pon_mineralisation-don*don_mineralisation)
   _SET_ODE_(self%id_pop,-pop*pop_mineralisation)
   _SET_ODE_(self%id_dop,pop*pop_mineralisation-dop*dop_mineralisation)
   _SET_ODE_(self%id_poc,-poc*poc_mineralisation)
   _SET_ODE_(self%id_doc,poc*poc_mineralisation-doc*doc_mineralisation)

   ! If an externally maintained oxygen pool is present, take mineralisation from it
   IF (self%use_oxy) THEN
      _SET_ODE_(self%id_oxy,-Yoxy_doc_miner*doc*doc_mineralisation)
   ENDIF
   if (self%use_dic) THEN
      _SET_ODE_(self%id_dic,doc*doc_mineralisation)
   ENDIF
   IF (self%use_amm) THEN
      _SET_ODE_(self%id_amm,don*don_mineralisation)
   ENDIF
   IF (self%use_frp) THEN
      _SET_ODE_(self%id_frp,dop*dop_mineralisation)
   ENDIF

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_pon_miner,pon_mineralisation)
   _SET_DIAGNOSTIC_(self%id_don_miner,don_mineralisation)
   _SET_DIAGNOSTIC_(self%id_pop_miner,pop_mineralisation)
   _SET_DIAGNOSTIC_(self%id_dop_miner,dop_mineralisation)
   _SET_DIAGNOSTIC_(self%id_poc_miner,poc_mineralisation)
   _SET_DIAGNOSTIC_(self%id_doc_miner,doc_mineralisation)

   _SET_DIAGNOSTIC_(self%id_bod,poc+doc)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_organic_matter_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_organic_matter_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of oxygen biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_organic_matter),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL           :: pon,don,amm,oxy,temp !State variables
   AED_REAL           :: pon_mineralisation, don_mineralisation
   AED_REAL           :: pop,dop,frp !State variables
   AED_REAL           :: pop_mineralisation, dop_mineralisation
   AED_REAL           :: poc,doc,dic !State variables
   AED_REAL           :: poc_mineralisation, doc_mineralisation
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
 ! AED_REAL,PARAMETER :: Yoxy_don_miner = 6.625   ! ratio of oxygen to nitrogen utilised during don mineralisation
 ! AED_REAL,PARAMETER :: Yoxy_dop_miner = 6.625   ! ratio of oxygen to phosphoros utilised during dop mineralisation
   AED_REAL,PARAMETER :: Yoxy_doc_miner = 32./12. ! ratio of oxygen to carbon utilised during doc mineralisation

!-----------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_
    !call log_message('model aed_organic_matter enter do loop successfully.')

   ! Retrieve current (local) state variable values.
   _GET_(self%id_pon,pon) ! particulate organic nitrogen
   _GET_(self%id_don,don) ! dissolved organic nitrogen
   _GET_(self%id_pop,pop) ! particulate organic phosphorus
   _GET_(self%id_dop,dop) ! dissolved organic phosphorus
   _GET_(self%id_poc,poc) ! particulate organic carbon
   _GET_(self%id_doc,doc) ! dissolved organic carbon

   IF (self%use_oxy) THEN
      _GET_(self%id_oxy,oxy) ! oxygen
   ELSE
      oxy = 0.0
   ENDIF
   IF (self%use_dic) THEN
      _GET_(self%id_dic,dic) ! disolved inorganic carbon
   ELSE
      dic = 0.0
   ENDIF
   IF (self%use_amm) THEN
      _GET_(self%id_amm,amm) ! ammonium
   ELSE
      amm = 0.0
   ENDIF
   IF (self%use_frp) THEN
      _GET_(self%id_frp,frp) ! phosphate
   ELSE
      frp = 0.0
   ENDIF

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)  ! temperature

   ! Define some intermediate quantities units mmol N/m3/day
   pon_mineralisation = fpon_miner(self%use_oxy,self%Rpon_miner,self%Kpon_miner,self%theta_pon_miner,oxy,temp)
   don_mineralisation = fdon_miner(self%use_oxy,self%Rdon_miner,self%Kdon_miner,self%theta_don_miner,oxy,temp)
   pop_mineralisation = fpop_miner(self%use_oxy,self%Rpop_miner,self%Kpop_miner,self%theta_pop_miner,oxy,temp)
   dop_mineralisation = fdop_miner(self%use_oxy,self%Rdop_miner,self%Kdop_miner,self%theta_dop_miner,oxy,temp)
   poc_mineralisation = fpoc_miner(self%use_oxy,self%Rpoc_miner,self%Kpoc_miner,self%theta_poc_miner,oxy,temp)
   doc_mineralisation = fdoc_miner(self%use_oxy,self%Rdoc_miner,self%Kdoc_miner,self%theta_doc_miner,oxy,temp)

   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_(i,j,val) as opposed to _SET_DD_(i,j,val),
   ! assignments to dd(i,j) are automatically assigned to pp(j,i) as well.
   !Set for particulate organic matter mineralisation
#if 0
   _SET_DD_SYM_(self%id_pon,self%id_don,pon*pon_mineralisation)
   _SET_DD_SYM_(self%id_pop,self%id_dop,pop*pop_mineralisation)
   _SET_DD_SYM_(self%id_poc,self%id_doc,poc*poc_mineralisation)

   ! If an externally maintained oxygen pool is present, take mineralisation from it
   IF (self%use_oxy) THEN
      _SET_DD_(self%id_oxy,self%id_oxy,Yoxy_doc_miner*doc*doc_mineralisation)
   ENDIF

   !If simulating inorganic nutrients then add dissolved organic mineralisation
   !Otherwise just take from the dissolved organic matter
   IF (self%use_amm) THEN
      _SET_DD_SYM_(self%id_don,self%id_amm,don*don_mineralisation)
   ELSE
      _SET_DD_(self%id_don,self%id_don,don*don_mineralisation)
   ENDIF
   IF (self%use_frp) THEN
      _SET_DD_SYM_(self%id_dop,self%id_frp,dop*dop_mineralisation)
   ELSE
      _SET_DD_(self%id_dop,self%id_dop,dop*dop_mineralisation)
   ENDIF
   IF (self%use_dic) THEN
      _SET_DD_SYM_(self%id_doc,self%id_dic,doc*doc_mineralisation)
   ELSE
      _SET_DD_(self%id_doc,self%id_doc,doc*doc_mineralisation)
   ENDIF
#endif
   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_pon_miner,pon_mineralisation)
   _SET_DIAGNOSTIC_(self%id_don_miner,don_mineralisation)
   _SET_DIAGNOSTIC_(self%id_pop_miner,pop_mineralisation)
   _SET_DIAGNOSTIC_(self%id_dop_miner,dop_mineralisation)
   _SET_DIAGNOSTIC_(self%id_poc_miner,poc_mineralisation)
   _SET_DIAGNOSTIC_(self%id_doc_miner,doc_mineralisation)

   _SET_DIAGNOSTIC_(self%id_bod,poc+doc)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_organic_matter_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_organic_matter_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED nitrogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_organic_matter),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: pon,don
   AED_REAL :: pop,dop
   AED_REAL :: poc,doc

   ! Temporary variables
   AED_REAL :: pon_flux,don_flux
   AED_REAL :: pop_flux,dop_flux
   AED_REAL :: poc_flux,doc_flux

   AED_REAL :: Fsed_pon,Fsed_don
   AED_REAL :: Fsed_pop,Fsed_dop
   AED_REAL :: Fsed_poc,Fsed_doc
   AED_REAL :: Psed_poc, Psed_pon, Psed_pop

   ! Parameters
   AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN

   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

   ! Retrieve current (local) state variable values.
   _GET_(self%id_pon,pon) ! particulate organic matter
   _GET_(self%id_don,don) ! particulate organic matter
   _GET_(self%id_pop,pop) ! particulate organic matter
   _GET_(self%id_dop,dop) ! particulate organic matter
   _GET_(self%id_poc,poc) ! particulate organic matter
   _GET_(self%id_doc,doc) ! particulate organic matter

   IF (self%use_sed_model) THEN
      _GET_HORIZONTAL_(self%id_Fsed_pon,Fsed_pon)
      _GET_HORIZONTAL_(self%id_Fsed_don,Fsed_don)
      _GET_HORIZONTAL_(self%id_Fsed_pop,Fsed_pop)
      _GET_HORIZONTAL_(self%id_Fsed_dop,Fsed_dop)
      _GET_HORIZONTAL_(self%id_Fsed_poc,Fsed_poc)
      _GET_HORIZONTAL_(self%id_Fsed_doc,Fsed_doc)
   ELSE
      Fsed_pon = self%Fsed_pon
      Fsed_don = self%Fsed_don * self%Ksed_don/(self%Ksed_don+don) * (self%theta_sed_don**(temp-20.0))
      Fsed_pop = self%Fsed_pop
      Fsed_dop = self%Fsed_dop * self%Ksed_dop/(self%Ksed_dop+dop) * (self%theta_sed_dop**(temp-20.0))
      Fsed_poc = self%Fsed_poc
      Fsed_doc = self%Fsed_doc * self%Ksed_doc/(self%Ksed_doc+doc) * (self%theta_sed_doc**(temp-20.0))
   ENDIF

   ! Calculate sedimentation flux (mmmol/m2/s) loss from benthos.
   IF (self%use_sedmtn_model) THEN
       Psed_poc = self%w_poc * max(zero_,poc)
       Psed_pon = self%w_pon * max(zero_,pon)
       Psed_pop = self%w_pop * max(zero_,pop)
   ELSE
       Psed_poc = zero_
       Psed_pon = zero_
       Psed_pop = zero_
   ENDIF

   pon_flux = Fsed_pon + Psed_pon
   don_flux = Fsed_don
   pop_flux = Fsed_pop + Psed_pop
   dop_flux = Fsed_dop
   poc_flux = Fsed_poc + Psed_poc
   doc_flux = Fsed_doc

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to FABM.
   _SET_BOTTOM_EXCHANGE_(self%id_pon,pon_flux)
   _SET_BOTTOM_EXCHANGE_(self%id_don,don_flux)
   _SET_BOTTOM_EXCHANGE_(self%id_pop,pop_flux)
   _SET_BOTTOM_EXCHANGE_(self%id_dop,dop_flux)
   _SET_BOTTOM_EXCHANGE_(self%id_poc,poc_flux)
   _SET_BOTTOM_EXCHANGE_(self%id_doc,doc_flux)


   ! Set sedimentation flux (mmmol/m2) as calculated by organic matter.
   IF (self%use_sedmtn_model) THEN
      !_SET_STATE_BEN_(self%id_Psed_poc,Psed_poc)
      !_SET_STATE_BEN_(self%id_Psed_pon,Psed_pon)
      !_SET_STATE_BEN_(self%id_Psed_pop,Psed_pop)
   ENDIF

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_amm,-amm_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_pon,-pon_flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_don,-don_flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_pop,-pop_flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_dop,-dop_flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_poc,-poc_flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_doc,-doc_flux)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

END SUBROUTINE aed_organic_matter_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_organic_matter_get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_organic_matter),INTENT(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
!LOCALS
   AED_REAL :: doc,poc
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_doc,doc)
   _GET_(self%id_poc,poc)

   ! Self-shading with explicit contribution from background OM concentration.
   _SET_EXTINCTION_(self%KeDOM*doc +self%KePOM*poc)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_organic_matter_get_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fpon_miner(use_oxy,Rpon_miner,Kpon_miner,theta_pon_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Nitrogen
!
! Michaelis-Menten formulation for mineralisation
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpon_miner,Kpon_miner,theta_pon_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpon_miner = Rpon_miner * oxy/(Kpon_miner+oxy) * (theta_pon_miner**(temp-20.0))
   ELSE
      fpon_miner = Rpon_miner * (theta_pon_miner**(temp-20.0))
   ENDIF

END FUNCTION fpon_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdon_miner(use_oxy,Rdon_miner,Kdon_miner,theta_don_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdon_miner,Kdon_miner,theta_don_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdon_miner = Rdon_miner * oxy/(Kdon_miner+oxy) * (theta_don_miner**(temp-20.0))
   ELSE
      fdon_miner = Rdon_miner * (theta_don_miner**(temp-20.0))
   ENDIF

END FUNCTION fdon_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
PURE AED_REAL FUNCTION fpop_miner(use_oxy,Rpop_miner,Kpop_miner,theta_pop_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Phosphorus
!
! Michaelis-Menten formulation for mineralisation
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpop_miner,Kpop_miner,theta_pop_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpop_miner = Rpop_miner * oxy/(Kpop_miner+oxy) * (theta_pop_miner**(temp-20.0))
   ELSE
      fpop_miner = Rpop_miner * (theta_pop_miner**(temp-20.0))
   ENDIF

END FUNCTION fpop_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdop_miner(use_oxy,Rdop_miner,Kdop_miner,theta_dop_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdop_miner,Kdop_miner,theta_dop_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdop_miner = Rdop_miner * oxy/(Kdop_miner+oxy) * (theta_dop_miner**(temp-20.0))
   ELSE
      fdop_miner = Rdop_miner * (theta_dop_miner**(temp-20.0))
   ENDIF

END FUNCTION fdop_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
PURE AED_REAL FUNCTION fpoc_miner(use_oxy,Rpoc_miner,Kpoc_miner,theta_poc_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Carbon
!
! Michaelis-Menten formulation for mineralisation
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rpoc_miner,Kpoc_miner,theta_poc_miner,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fpoc_miner = Rpoc_miner * oxy/(Kpoc_miner+oxy) * (theta_poc_miner**(temp-20.0))
   ELSE
      fpoc_miner = Rpoc_miner * (theta_poc_miner**(temp-20.0))
   ENDIF

END FUNCTION fpoc_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdoc_miner(use_oxy,Rdoc_miner,Kdoc_miner,theta_doc_miner,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for mineralisation added 18/7/11
!
! Here, the classical Michaelis-Menten formulation for mineralisation
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in) :: use_oxy
   AED_REAL,INTENT(in) :: Rdoc_miner,Kdoc_miner,theta_doc_miner,oxy,temp
!
!-----------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdoc_miner = Rdoc_miner * oxy/(Kdoc_miner+oxy) * (theta_doc_miner**(temp-20.0))
   ELSE
      fdoc_miner = Rdoc_miner * (theta_doc_miner**(temp-20.0))
   ENDIF

END FUNCTION fdoc_miner
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_organic_matter

!###############################################################################
!#                                                                             #
!# aed_chlorophylla.F90                                                        #
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

#include "fabm_driver.h"

MODULE aed_chlorophylla
!-------------------------------------------------------------------------------
! aed_chlorphylla --- simple lumped chl-a model
!-------------------------------------------------------------------------------
   USE fabm_types
   USE fabm_driver

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC type_aed_chla, aed_chla_create
!
   TYPE,extends(type_base_model) :: type_aed_chla
!     Variable identifiers
      _TYPE_STATE_VARIABLE_ID_      :: id_p
      _TYPE_STATE_VARIABLE_ID_      :: id_exctarget,id_morttarget,id_upttarget
      _TYPE_DEPENDENCY_ID_          :: id_par,id_I_0
      _TYPE_DIAGNOSTIC_VARIABLE_ID_ :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR
      _TYPE_CONSERVED_QUANTITY_ID_  :: id_totN

!     Model parameters
      REALTYPE :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      REALTYPE :: dic_per_n
      logical  :: do_exc,do_mort,do_upt

      CONTAINS  ! Model Procedures
!       procedure :: initialize               => aed_chla_init
        procedure :: do                       => aed_chla_do
        procedure :: do_ppdd                  => aed_chla_do_ppdd
        procedure :: get_light_extinction     => aed_chla_get_light_extinction
        procedure :: get_conserved_quantities => aed_chla_get_conserved_quantities
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
FUNCTION aed_chla_create(namlst,name,parent) RESULT(self)
!-------------------------------------------------------------------------------
! Initialise the NPZD model
!
!  Here, the npzd namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                    :: namlst
   CHARACTER(len=*),INTENT(in)           :: name
   _CLASS_ (type_model_info),TARGET,INTENT(inout) :: parent
!
!LOCALS
   _CLASS_ (type_aed_chla),POINTER       :: self

   REALTYPE           :: p_initial=0.
   REALTYPE           :: p0=0.0225
   REALTYPE           :: w_p=-1.157407e-05
   REALTYPE           :: i_min=25.
   REALTYPE           :: rmax=1.157407e-05
   REALTYPE           :: alpha=0.3
   REALTYPE           :: rpn=1.157407e-07
   REALTYPE           :: rpdu=2.314814e-07
   REALTYPE           :: rpdl=1.157407e-06
   CHARACTER(len=64)  :: excretion_target_variable=''
   CHARACTER(len=64)  :: mortality_target_variable=''
   CHARACTER(len=64)  :: uptake_target_variable=''

   REALTYPE,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_chla/ p_initial,p0,w_p,i_min,rmax,alpha,rpn,rpdu,rpdl, &
                    excretion_target_variable,mortality_target_variable,uptake_target_variable

!-------------------------------------------------------------------------------
!BEGIN
   ALLOCATE(self)
   CALL self%initialize(name,parent)

   ! Read the namelist
   read(namlst,nml=aed_chla,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%p0    = p0
   self%i_min = i_min
   self%rmax  = rmax/secs_pr_day
   self%alpha = alpha
   self%rpn  = rpn /secs_pr_day
   self%rpdu = rpdu/secs_pr_day
   self%rpdl = rpdl/secs_pr_day

   ! Register state variables
   self%id_p = self%register_state_variable('phy','mmol/m**3','phytoplankton', &
                                    p_initial,minimum=_ZERO_,vertical_movement=w_p/secs_pr_day)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   self%do_exc = excretion_target_variable .NE. ''
   IF (self%do_exc) self%id_exctarget = self%register_state_dependency(excretion_target_variable)
   self%do_mort = mortality_target_variable .NE. ''
   IF (self%do_mort) self%id_morttarget = self%register_state_dependency(mortality_target_variable)
   self%do_upt = uptake_target_variable .NE. ''
   IF (self%do_upt) self%id_upttarget = self%register_state_dependency(uptake_target_variable)


   ! Register diagnostic variables
   self%id_GPP  = self%register_diagnostic_variable('GPP','mmol/m**3',  'gross primary production',           &
                     time_treatment=time_treatment_step_integrated)
   self%id_NCP  = self%register_diagnostic_variable('NCP','mmol/m**3',  'net community production',           &
                     time_treatment=time_treatment_step_integrated)
   self%id_PPR  = self%register_diagnostic_variable('PPR','mmol/m**3/d','gross primary production rate',      &
                     time_treatment=time_treatment_averaged)
   self%id_NPR  = self%register_diagnostic_variable('NPR','mmol/m**3/d','net community production rate',      &
                     time_treatment=time_treatment_averaged)
   self%id_dPAR = self%register_diagnostic_variable('PAR','W/m**2',     'photosynthetically active radiation',&
                     time_treatment=time_treatment_averaged)

   ! Register conserved quantities
   self%id_totN = self%register_conserved_quantity('N','mmol/m**3','nitrogen')

   ! Register environmental dependencies
   self%id_par = self%register_dependency(varname_par)
   self%id_I_0 = self%register_dependency(varname_par_sf, shape=shape_hz)


   PRINT *,'AED_CHLA : Note this module has not been completed. Stopping.'
   STOP

   RETURN

99 call fatal_error('aed_chla_init','Error reading namelist aed_chla')

END FUNCTION aed_chla_create
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_chla_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of chlorophylla model
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_chla),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   REALTYPE           :: n,p,par,I_0
   REALTYPE           :: iopt,rpd,primprod
   REALTYPE,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p,p) ! phytoplankton
   _GET_STATE_(self%id_upttarget,n) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_DEPENDENCY_   (self%id_par,par)  ! local photosynthetically active radiation
   _GET_DEPENDENCY_HZ_(self%id_I_0,I_0)  ! surface short wave radiation

   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,self%I_min)

   ! Loss rate of phytoplankton to detritus depends on local light intensity.
   IF (par .ge. self%I_min) THEN
      rpd = self%rpdu
   ELSE
      rpd = self%rpdl
   ENDIF

   ! Define some intermediate quantities that will be reused multiple times.
   primprod = fnp(self,n,p,par,iopt)

   ! Set temporal derivatives
   _SET_ODE_(self%id_p,primprod - self%rpn*p - rpd*p)

   ! If an externally maintained ...
   IF (self%do_upt) THEN
      _SET_ODE_(self%id_upttarget,-primprod)
   ENDIF
   IF (self%do_mort) THEN
      _SET_ODE_(self%id_morttarget,rpd*p)
   ENDIF
   IF (self%do_exc) THEN
      _SET_ODE_(self%id_exctarget,self%rpn*p)
   ENDIF

   ! Export diagnostic variables
   _SET_DIAG_(self%id_dPAR,par)
   _SET_DIAG_(self%id_GPP ,primprod)
   _SET_DIAG_(self%id_NCP ,primprod - self%rpn*p)
   _SET_DIAG_(self%id_PPR ,primprod*secs_pr_day)
   _SET_DIAG_(self%id_NPR ,(primprod - self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_chla_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_chla_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
!-------------------------------------------------------------------------------
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_chla),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_EXTINCTION_
!
!LOCALS
   REALTYPE                     :: p
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p,p) ! phytoplankton

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(0.0*p)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_chla_get_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_chla_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities (currently only nitrogen)
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_chla),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_

!LOCALS
   REALTYPE :: p
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p,p) ! phytoplankton

   ! Total nutrient is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totN,p)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_chla_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_chla_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of NPZD model exporting production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_chla),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   REALTYPE           :: n,p,par,I_0=0.
   REALTYPE           :: iopt,rpd,primprod
   REALTYPE,PARAMETER :: secs_pr_day = 86400.
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_STATE_(self%id_p,p) ! phytoplankton
   _GET_STATE_(self%id_upttarget,n) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_DEPENDENCY_   (self%id_par,par)  ! local photosynthetically active radiation

   ! Light acclimation formulation based on surface light intensity.
   iopt = max(0.25*I_0,self%I_min)

   ! Loss rate of phytoplankton to detritus depends on local light intensity.
   IF (par .ge. self%I_min) THEN
      rpd = self%rpdu
   ELSE
      rpd = self%rpdl
   ENDIF

   ! Rate of primary production will be reused multiple times - calculate it once.
   primprod = fnp(self,n,p,par,iopt)


   _SET_DD_(self%id_p,self%id_p,primprod)
   _SET_DD_(self%id_p,self%id_p,- self%rpn*p)
   _SET_DD_(self%id_p,self%id_p,- rpd*p)

   ! If an externally maintained ...
   IF (self%do_upt) THEN
      _SET_PP_(self%id_upttarget,self%id_upttarget,-primprod)
   ENDIF
   IF (self%do_mort) THEN
      _SET_PP_(self%id_morttarget,self%id_morttarget,rpd*p)
   ENDIF
   IF (self%do_exc) THEN
      _SET_PP_(self%id_exctarget,self%id_exctarget,self%rpn*p)
   ENDIF

   ! Export diagnostic variables
   _SET_DIAG_(self%id_dPAR,par)
   _SET_DIAG_(self%id_GPP,primprod)
   _SET_DIAG_(self%id_NCP,primprod-self%rpn*p)
   _SET_DIAG_(self%id_PPR,primprod*secs_pr_day)
   _SET_DIAG_(self%id_NPR,(primprod-self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_chla_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE REALTYPE FUNCTION fnp(self,n,p,par,iopt)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nutrient uptake
!-------------------------------------------------------------------------------
!ARGUMENTS
   _CLASS_ (type_aed_chla),INTENT(in) :: self
   REALTYPE,INTENT(in)                :: n,p,par,iopt
!
!-------------------------------------------------------------------------------
!BEGIN
   fnp = self%rmax*par/iopt*exp(_ONE_-par/iopt)*n/(self%alpha+n)*(p+self%p0)

END FUNCTION fnp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed_chlorophylla
#endif

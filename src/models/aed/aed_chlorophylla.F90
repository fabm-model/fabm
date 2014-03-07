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

#include "aed.h"

MODULE aed_chlorophylla
!-------------------------------------------------------------------------------
! aed_chlorphylla --- simple lumped chl-a model
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_type_chla
!
   TYPE,extends(type_base_model) :: aed_type_chla
!     Variable identifiers
      TYPE (type_state_variable_id)      :: id_p
      TYPE (type_state_variable_id)      :: id_exctarget,id_morttarget,id_upttarget
      TYPE (type_dependency_id)          :: id_par
      TYPE (type_horizontal_dependency_id)  :: id_I_0
      TYPE (type_diagnostic_variable_id) :: id_GPP,id_NCP,id_PPR,id_NPR,id_dPAR

!     Model parameters
      AED_REAL :: p0,z0,kc,i_min,rmax,gmax,iv,alpha,rpn,rzn,rdn,rpdu,rpdl,rzd
      AED_REAL :: dic_per_n
      LOGICAL  :: do_exc,do_mort,do_upt

      CONTAINS  ! Model Procedures
        PROCEDURE :: initialize               => aed_init_chla
        PROCEDURE :: do                       => aed_chla_do
        PROCEDURE :: do_ppdd                  => aed_chla_do_ppdd
        PROCEDURE :: get_light_extinction     => aed_chla_get_light_extinction
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_chla(self, namlst)
!-------------------------------------------------------------------------------
! Initialise the chlorophyl model
!
!  Here, the chla namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_chla),TARGET,INTENT(inout) :: self
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status

   AED_REAL           :: p_initial=0.
   AED_REAL           :: p0=0.0225
   AED_REAL           :: w_p=-1.157407e-05
   AED_REAL           :: i_min=25.
   AED_REAL           :: rmax=1.157407e-05
   AED_REAL           :: alpha=0.3
   AED_REAL           :: rpn=1.157407e-07
   AED_REAL           :: rpdu=2.314814e-07
   AED_REAL           :: rpdl=1.157407e-06
   CHARACTER(len=64)  :: excretion_target_variable=''
   CHARACTER(len=64)  :: mortality_target_variable=''
   CHARACTER(len=64)  :: uptake_target_variable=''

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_chla/ p_initial,p0,w_p,i_min,rmax,alpha,rpn,rpdu,rpdl, &
                    excretion_target_variable,mortality_target_variable,uptake_target_variable

!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_chla,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_chla'

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
   CALL self%register_state_variable(self%id_p,'phy','mmol/m**3','phytoplankton', &
                                    p_initial,minimum=zero_,vertical_movement=w_p/secs_pr_day)

   ! Register link to external DIC pool, if DIC variable name is provided in namelist.
   self%do_exc = excretion_target_variable .NE. ''
   IF (self%do_exc) call self%register_state_dependency(self%id_exctarget,excretion_target_variable)
   self%do_mort = mortality_target_variable .NE. ''
   IF (self%do_mort) call self%register_state_dependency(self%id_morttarget,mortality_target_variable)
   self%do_upt = uptake_target_variable .NE. ''
   IF (self%do_upt) call self%register_state_dependency(self%id_upttarget,uptake_target_variable)


   ! Register diagnostic variables
   CALL self%register_diagnostic_variable(self%id_GPP,'GPP','mmol/m**3',  'gross primary production')
   CALL self%register_diagnostic_variable(self%id_NCP,'NCP','mmol/m**3',  'net community production')
   CALL self%register_diagnostic_variable(self%id_PPR,'PPR','mmol/m**3/d','gross primary production rate')
   CALL self%register_diagnostic_variable(self%id_NPR,'NPR','mmol/m**3/d','net community production rate')
   CALL self%register_diagnostic_variable(self%id_dPAR,'PAR','W/m**2',    'photosynthetically active radiation')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
   CALL self%register_dependency(self%id_I_0,standard_variables%surface_downwelling_photosynthetic_radiative_flux)


   PRINT *,'AED_CHLA : Note this module has not been completed.'
END SUBROUTINE aed_init_chla
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_chla_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of chlorophylla model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_chla),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   AED_REAL           :: n,p,par,I_0
   AED_REAL           :: iopt,rpd,primprod
   AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_upttarget,n) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_   (self%id_par,par)  ! local photosynthetically active radiation
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

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
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP ,primprod)
   _SET_DIAGNOSTIC_(self%id_NCP ,primprod - self%rpn*p)
   _SET_DIAGNOSTIC_(self%id_PPR ,primprod*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR ,(primprod - self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_chla_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_chla_get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to biogeochemical variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_chla),INTENT(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
!LOCALS
   AED_REAL :: p
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(0.0*p)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_chla_get_light_extinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_chla_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of NPZD model exporting production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_chla),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL           :: n,p,par,I_0=0.
   AED_REAL           :: iopt,rpd,primprod
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p) ! phytoplankton
   _GET_(self%id_upttarget,n) ! nutrients

   ! Retrieve current environmental conditions.
   _GET_   (self%id_par,par)  ! local photosynthetically active radiation

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
   _SET_DIAGNOSTIC_(self%id_dPAR,par)
   _SET_DIAGNOSTIC_(self%id_GPP,primprod)
   _SET_DIAGNOSTIC_(self%id_NCP,primprod-self%rpn*p)
   _SET_DIAGNOSTIC_(self%id_PPR,primprod*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_NPR,(primprod-self%rpn*p)*secs_pr_day)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_chla_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fnp(self,n,p,par,iopt)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nutrient uptake
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_chla),INTENT(in) :: self
   AED_REAL,INTENT(in)                :: n,p,par,iopt
!
!-------------------------------------------------------------------------------
!BEGIN
   fnp = self%rmax*par/iopt*exp(one_-par/iopt)*n/(self%alpha+n)*(p+self%p0)

END FUNCTION fnp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



END MODULE aed_chlorophylla

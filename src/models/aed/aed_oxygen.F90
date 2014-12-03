!###############################################################################
!#                                                                             #
!# aed_oxygen.F90                                                              #
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
!# Created May 2011                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"

MODULE aed_oxygen
!-------------------------------------------------------------------------------
! aed_oxygen --- oxygen biogeochemical model
!
! The AED module oxygen contains equations that describe exchange of
! oxygen across the air/water interface and sediment flux.
!-------------------------------------------------------------------------------
   USE aed_core

   USE aed_util,  ONLY: aed_gas_piston_velocity, aed_oxygen_sat

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_type_oxygen
!
   TYPE,extends(type_base_model) :: aed_type_oxygen
!     Variable identifiers
      TYPE (type_state_variable_id)      :: id_oxy
      TYPE (type_dependency_id)          :: id_temp, id_salt
      TYPE (type_horizontal_dependency_id) :: id_wind
      TYPE (type_horizontal_dependency_id) :: id_Fsed_oxy
      TYPE (type_horizontal_diagnostic_variable_id) :: id_oxy_sat !, id_atm_oxy_exch3d
      TYPE (type_horizontal_diagnostic_variable_id) :: id_atm_oxy_exch
      TYPE (type_horizontal_diagnostic_variable_id) :: id_sed_oxy

!     Model parameters
      AED_REAL :: Fsed_oxy,Ksed_oxy,theta_sed_oxy
      LOGICAL  :: use_sed_model

      CONTAINS     ! Model Procedures
        PROCEDURE :: initialize               => aed_init_oxygen
        PROCEDURE :: do                       => aed_oxygen_do
        PROCEDURE :: do_ppdd                  => aed_oxygen_do_ppdd
        PROCEDURE :: do_benthos               => aed_oxygen_do_benthos
        PROCEDURE :: get_surface_exchange     => aed_oxygen_get_surface_exchange
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_oxygen(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the aed_oxygen model
!
!  Here, the oxygen namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_oxygen),TARGET,INTENT(inout) :: self
   INTEGER,INTENT(in) :: namlst

!
!LOCALS
   INTEGER  :: status
   AED_REAL :: oxy_initial=300.
   AED_REAL :: Fsed_oxy = 48.0
   AED_REAL :: Ksed_oxy = 30.0
   AED_REAL :: theta_sed_oxy = 1.0
   CHARACTER(len=64) :: Fsed_oxy_variable=''

   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   NAMELIST /aed_oxygen/ oxy_initial,Fsed_oxy,Ksed_oxy,theta_sed_oxy,  &
                         Fsed_oxy_variable
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_oxygen,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_oxygen'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.

   self%Fsed_oxy = Fsed_oxy/secs_pr_day
   self%Ksed_oxy = Ksed_oxy
   self%theta_sed_oxy = theta_sed_oxy
   self%use_sed_model = Fsed_oxy_variable .NE. ''

   ! Register state variables
   CALL self%register_state_variable(self%id_oxy,'oxy','mmol/m**3','oxygen',   &
                                    oxy_initial,minimum=zero_,no_river_dilution=.TRUE.)

   ! Register link to external pools

   IF (self%use_sed_model) THEN
      CALL self%register_horizontal_dependency(self%id_Fsed_oxy,Fsed_oxy_variable)
      CALL self%request_coupling(self%id_Fsed_oxy,Fsed_oxy_variable)
   ENDIF

   ! Register diagnostic variables
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_oxy,       &
                     'sed_oxy', 'mmol/m**2/d', 'Oxygen sediment flux')

   CALL self%register_horizontal_diagnostic_variable(self%id_atm_oxy_exch,  &
                     'atm_oxy_exch', 'mmol/m**2/d', 'Oxygen exchange across atm/water interface')

!  CALL self%register_horizontal_diagnostic_variable(self%id_atm_oxy_exch3d, &
!                    'atm_oxy_exch3d', 'mmol/m**2/d', 'Oxygen exchange across atm/water interface')

   CALL self%register_horizontal_diagnostic_variable(self%id_oxy_sat,       &
                     'sat', 'mmol/m**2/d', 'Oxygen saturation')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature) ! Temperature (degrees Celsius)
   CALL self%register_dependency(self%id_salt,standard_variables%practical_salinity) ! Salinity (psu)
!  CALL self%register_dependency(self%id_pres,standard_variables%pressure) ! Pressure (dbar = 10 kPa)
   CALL self%register_dependency(self%id_wind,standard_variables%wind_speed) ! Wind speed at 10 m above surface (m/s)

END SUBROUTINE aed_init_oxygen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_oxygen_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of oxygen biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_oxygen),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL :: oxy
   AED_REAL :: diff_oxy
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy) ! oxygen

   ! Set temporal derivatives
   diff_oxy = 0.

   _SET_PP_(self%id_oxy,self%id_oxy,diff_oxy)

   ! If an externally maintained pool is present, change the  pool according

   ! Export diagnostic variables

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_oxygen_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_oxygen_get_surface_exchange(self,_ARGUMENTS_DO_SURFACE_)
!-------------------------------------------------------------------------------
! Air-water exchange for the aed oxygen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_oxygen),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
!LOCALS
   ! Environment
   AED_REAL :: temp, salt, wind

   ! State
   AED_REAL :: oxy

   ! Temporary variables
   AED_REAL :: oxy_atm_flux = zero_
   AED_REAL :: Coxy_air = zero_ !Dissolved oxygen in the air phase
   AED_REAL :: koxy_trans = zero_
   AED_REAL :: windHt !, Tabs
   AED_REAL :: f_pres  = 1.0      ! Pressure correction function only applicable at high altitudes
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   !Get dependent state variables from physical driver
   _GET_(self%id_temp,temp)     ! Temperature (degrees Celsius)
   _GET_(self%id_salt,salt)     ! Salinity (psu)
   _GET_HORIZONTAL_(self%id_wind,wind)  ! Wind speed at 10 m above surface (m/s)
   windHt = 10.

    ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy) ! Concentration of oxygen in surface layer

   koxy_trans = aed_gas_piston_velocity(windHt,wind,temp,salt)

   ! First get the oxygen concentration in the air phase at interface
   ! Taken from Riley and Skirrow (1974)
   f_pres = 1.0
   Coxy_air = f_pres * aed_oxygen_sat(salt,temp)

   ! Get the oxygen flux
   oxy_atm_flux = koxy_trans * (Coxy_air - oxy)

   ! Transfer surface exchange value to FABM (mmmol/m2) converted by driver.
   _SET_SURFACE_EXCHANGE_(self%id_oxy,oxy_atm_flux)

   ! Also store oxygen flux across the atm/water interface as diagnostic variable (mmmol/m2).
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_atm_oxy_exch,oxy_atm_flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_oxy_sat, Coxy_air)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_
END SUBROUTINE aed_oxygen_get_surface_exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_oxygen_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_oxygen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS(aed_type_oxygen),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   AED_REAL :: oxy
   AED_REAL :: diff_oxy
!  AED_REAL,PARAMETER :: secs_pr_day = 86400.

!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy) ! oxygen

   ! Set temporal derivatives
   diff_oxy = 0.

   _SET_ODE_(self%id_oxy,diff_oxy)

   ! If an externally maintained pool is present, change the pool according

   ! Export diagnostic variables

   ! Leave spatial loops (if any)
   _LOOP_END_
END SUBROUTINE aed_oxygen_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_oxygen_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED oxygen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_oxygen),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: oxy

   ! Temporary variables
   AED_REAL :: oxy_flux, Fsed_oxy

   ! Parameters
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy) ! oxygen

   IF (self%use_sed_model) THEN
       _GET_HORIZONTAL_(self%id_Fsed_oxy,Fsed_oxy)
   ELSE
       Fsed_oxy = self%Fsed_oxy
   ENDIF

    ! Sediment flux dependent on oxygen and temperature
   oxy_flux = Fsed_oxy * oxy/(self%Ksed_oxy+oxy) * (self%theta_sed_oxy**(temp-20.0))


   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to FABM.
   _SET_BOTTOM_EXCHANGE_(self%id_oxy,oxy_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_oxy,-oxy_flux)

   ! Also store sediment flux as diagnostic variable.
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_oxy,oxy_flux)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_
END SUBROUTINE aed_oxygen_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_oxygen

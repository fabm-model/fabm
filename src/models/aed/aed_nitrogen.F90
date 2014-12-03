!###############################################################################
!#                                                                             #
!# aed_nitrogen.F90                                                            #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Earth & Environment                                           #
!# (C) The University of Western Australia                                     #
!#                                                                             #
!# Copyright by the AED-team @ UWA under the GNU Public License - www.gnu.org  #
!#                                                                             #
!#   ----------------------------------------------------------------------    #
!#                                                                             #
!# Created 9 May 2011                                                          #
!#                                                                             #
!###############################################################################

#include "aed.h"

!
MODULE aed_nitrogen
!-------------------------------------------------------------------------------
! aed_nitrogen --- nitrogen biogeochemical model
!
! Nitrogen module contains equations for nitrification and deitrification
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE    ! By default make everything private
!
   PUBLIC aed_type_nitrogen
!
   TYPE,extends(type_base_model) :: aed_type_nitrogen
!     Variable identifiers
      TYPE (type_state_variable_id)      :: id_nit, id_amm !nitrate & ammonium
      TYPE (type_state_variable_id)      :: id_oxy,id_denit_product
      TYPE (type_dependency_id)          :: id_temp
      TYPE (type_horizontal_dependency_id)      :: id_Fsed_amm,id_Fsed_nit
      TYPE (type_diagnostic_variable_id) :: id_nitrif,id_denit
      TYPE (type_horizontal_diagnostic_variable_id) :: id_sed_amm,id_sed_nit

!     Model parameters
      AED_REAL :: Rnitrif,Rdenit,Fsed_amm,Fsed_nit,Knitrif,Kdenit,Ksed_amm,Ksed_nit, &
                          theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit
      LOGICAL  :: use_oxy,use_no2,use_sed_model

      CONTAINS   ! Model Procedures
        PROCEDURE :: initialize               => aed_init_nitrogen
        PROCEDURE :: do                       => aed_nitrogen_do
        PROCEDURE :: do_ppdd                  => aed_nitrogen_do_ppdd
        PROCEDURE :: do_benthos               => aed_nitrogen_do_benthos
   END TYPE


!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_init_nitrogen(self,namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_type_nitrogen),TARGET,INTENT(inout) :: self
!
!LOCALS
   INTEGER  :: status

   AED_REAL          :: nit_initial=4.5
   AED_REAL          :: amm_initial=4.5
   AED_REAL          :: Rnitrif = 0.01
   AED_REAL          :: Rdenit = 0.01
   AED_REAL          :: Fsed_amm = 3.5
   AED_REAL          :: Fsed_nit = 3.5
   AED_REAL          :: Knitrif = 150.0
   AED_REAL          :: Kdenit = 150.0
   AED_REAL          :: Ksed_amm = 30.0
   AED_REAL          :: Ksed_nit = 30.0
   AED_REAL          :: theta_nitrif = 1.0
   AED_REAL          :: theta_denit = 1.0
   AED_REAL          :: theta_sed_amm = 1.0
   AED_REAL          :: theta_sed_nit = 1.0
   CHARACTER(len=64) :: nitrif_reactant_variable=''
   CHARACTER(len=64) :: denit_product_variable=''
   CHARACTER(len=64) :: Fsed_amm_variable=''
   CHARACTER(len=64) :: Fsed_nit_variable=''


   AED_REAL, parameter :: secs_pr_day = 86400.
   NAMELIST /aed_nitrogen/ nit_initial,amm_initial,Rnitrif,Rdenit,Fsed_amm,Fsed_nit, &
                    Knitrif,Kdenit,Ksed_amm,Ksed_nit,                     &
                    theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit, &
                    nitrif_reactant_variable,denit_product_variable,      &
                    Fsed_amm_variable, Fsed_nit_variable
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(namlst,nml=aed_nitrogen,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_nitrogen'

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   self%Rnitrif  = Rnitrif/secs_pr_day
   self%Rdenit   = Rdenit/secs_pr_day
   self%Fsed_amm = Fsed_amm/secs_pr_day
   self%Fsed_nit = Fsed_nit/secs_pr_day
   self%Knitrif  = Knitrif
   self%Kdenit   = Kdenit
   self%Ksed_amm  = Ksed_amm
   self%Ksed_nit  = Ksed_nit
   self%theta_nitrif = theta_nitrif
   self%theta_denit  = theta_denit
   self%theta_sed_amm = theta_sed_amm
   self%theta_sed_nit = theta_sed_nit

   ! Register state variables
   CALL self%register_state_variable(self%id_amm,'amm','mmol/m**3','ammonium',            &
                                    amm_initial,minimum=zero_,no_river_dilution=.true.)
   CALL self%register_state_variable(self%id_nit,'nit','mmol/m**3','nitrate',             &
                                    nit_initial,minimum=zero_,no_river_dilution=.true.)
   ! Register external state variable dependencies
   self%use_oxy = nitrif_reactant_variable .NE. '' !This means oxygen module switched on
   IF (self%use_oxy) THEN
     CALL self%register_state_dependency(self%id_oxy,nitrif_reactant_variable)
   ENDIF
   self%use_no2 = denit_product_variable .NE. '' !This means n2 module switched on
   IF (self%use_no2) call self%register_state_dependency(self%id_denit_product,denit_product_variable)

   self%use_sed_model = Fsed_amm_variable .NE. ''
   IF (self%use_sed_model) THEN
     CALL self%register_horizontal_dependency(self%id_Fsed_amm,Fsed_amm_variable)
     CALL self%request_coupling(self%id_Fsed_amm,Fsed_amm_variable)
     CALL self%register_horizontal_dependency(self%id_Fsed_nit,Fsed_nit_variable)
     CALL self%request_coupling(self%id_Fsed_nit,Fsed_nit_variable)
   ENDIF

   ! Register diagnostic variables
   CALL self%register_diagnostic_variable(self%id_nitrif,'nitrif','mmol/m**3/d',       &
                                                         'Nitrification rate')
   CALL self%register_diagnostic_variable(self%id_denit,'denit','mmol/m**3/d',         &
                                                         'De-nitrification rate')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_amm,'sed_amm','mmol/m**2/d',      &
                                                         'Ammonium sediment flux')
   CALL self%register_horizontal_diagnostic_variable(self%id_sed_nit,'sed_nit','mmol/m**2/d',      &
                                                         'Nitrate sediment flux')

   ! Register environmental dependencies
   CALL self%register_dependency(self%id_temp,standard_variables%temperature)
END SUBROUTINE aed_init_nitrogen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_nitrogen_do(self,_ARGUMENTS_DO_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_nitrogen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_nitrogen),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
!LOCALS
   AED_REAL           :: amm,nit,oxy,temp !State variables
   AED_REAL           :: nitrification,denitrification
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   AED_REAL,PARAMETER :: Yoxy_nitrif = 3. !ratio of oxygen to nitrogen utilised during nitrification
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_amm,amm) ! ammonium
   _GET_(self%id_nit,nit) ! nitrate
   IF (self%use_oxy) THEN ! & use_oxy
      _GET_(self%id_oxy,oxy) ! oxygen
   ELSE
      oxy = 0.0
   ENDIF

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)  ! temperature

   ! Define some intermediate quantities units mmol N/m3/day
   nitrification = fnitrif(self%use_oxy,self%Rnitrif,self%Knitrif,self%theta_nitrif,oxy,temp)
   denitrification = fdenit(self%use_oxy,self%Rdenit,self%Kdenit,self%theta_denit,oxy,temp)

   ! Set temporal derivatives
   _SET_ODE_(self%id_amm,-amm*nitrification)
   _SET_ODE_(self%id_nit,amm*nitrification - nit*denitrification)

   ! If an externally maintained oxygen pool is present, take nitrification from it
   IF (self%use_oxy) then ! & use_oxy
      _SET_ODE_(self%id_oxy,-Yoxy_nitrif*amm*nitrification)
   ENDIF
   !if (self%use_no2) then
   !   _SET_ODE_(self%id_denit_product,denitrification)
   !end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_nitrif, nitrification)
   _SET_DIAGNOSTIC_(self%id_denit, denitrification)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_nitrogen_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_nitrogen_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of nitrogen biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_nitrogen),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_
!
!LOCALS
   AED_REAL           :: amm,nit,oxy,temp !State variables
   AED_REAL           :: nitrification,denitrification
   AED_REAL,PARAMETER :: secs_pr_day = 86400.
   AED_REAL,PARAMETER :: Yoxy_nitrif = 3. !ratio of oxygen to nitrogen utilised during nitrification
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_


   ! Retrieve current (local) state variable values.
   _GET_(self%id_amm,amm) ! ammonium
   _GET_(self%id_nit,nit) ! nitrate
   IF (self%use_oxy) THEN ! & use_oxy
      _GET_(self%id_oxy,oxy) ! oxygen
   ELSE
      oxy = 0.0
   ENDIF

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,temp)  ! temperature

   ! Define some intermediate quantities units mmol N/m3/day
   nitrification = fnitrif(self%use_oxy,self%Rnitrif,self%Knitrif,self%theta_nitrif,oxy,temp)
   denitrification = fdenit(self%use_oxy,self%Rdenit,self%Kdenit,self%theta_denit,oxy,temp)

   ! Assign destruction rates to different elements of the destruction matrix.
   ! By assigning with _SET_DD_SYM_(i,j,val) as opposed to _SET_DD_(i,j,val),
   ! assignments to dd(i,j) are automatically assigned to pp(j,i) as well.
   _SET_DD_SYM_(self%id_amm,self%id_nit,amm*nitrification)

   _SET_DD_(self%id_nit,self%id_nit,nit*denitrification)

   ! If an externally maintained oxygen pool is present, take nitrification from it
   IF (self%use_oxy) then ! & use_oxy
      _SET_DD_(self%id_oxy,self%id_oxy,Yoxy_nitrif*amm*nitrification)
   ENDIF
   !if (self%use_no2) then
   !   _SET_PP_(self%id_denit_product,self%id_denit_product,denitrification)
   !end if

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_nitrif ,nitrification)
   _SET_DIAGNOSTIC_(self%id_denit ,denitrification)

   ! Leave spatial loops (if any)
   _LOOP_END_

END SUBROUTINE aed_nitrogen_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_nitrogen_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED nitrogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_type_nitrogen),INTENT(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
!LOCALS
   ! Environment
   AED_REAL :: temp !, layer_ht

   ! State
   AED_REAL :: amm,nit,oxy

   ! Temporary variables
   AED_REAL :: amm_flux,nit_flux
   AED_REAL :: Fsed_amm, Fsed_nit

   ! Parameters
!  AED_REAL,PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_(self%id_amm,amm) ! ammonium
   _GET_(self%id_nit,nit) ! nitrate

   IF (self%use_sed_model) THEN
      _GET_HORIZONTAL_(self%id_Fsed_amm,Fsed_amm)
      _GET_HORIZONTAL_(self%id_Fsed_nit,Fsed_nit)
   ELSE
      Fsed_amm = self%Fsed_amm
      Fsed_nit = self%Fsed_nit
   ENDIF

   IF (self%use_oxy) THEN
      ! Sediment flux dependent on oxygen and temperature
      _GET_(self%id_oxy,oxy)
      amm_flux = Fsed_amm * self%Ksed_amm/(self%Ksed_amm+oxy) * (self%theta_sed_amm**(temp-20.0))
      nit_flux = Fsed_nit * oxy/(self%Ksed_nit+oxy) * (self%theta_sed_nit**(temp-20.0))
   ELSE
      ! Sediment flux dependent on temperature only.
      oxy = 0.
      amm_flux = Fsed_amm * (self%theta_sed_amm**(temp-20.0))
      nit_flux = Fsed_nit * (self%theta_sed_nit**(temp-20.0))
   ENDIF

   ! TODO:
   ! (1) Get benthic sink and source terms (sccb?) for current environment
   ! (2) Get pelagic bttom fluxes (per surface area - division by layer height will be handled at a higher level)

   ! Set bottom fluxes for the pelagic (change per surface area per second)
   ! Transfer sediment flux value to FABM.
   !_SET_BOTTOM_FLUX_(self%id_amm,amm_flux/secs_pr_day)
   _SET_BOTTOM_EXCHANGE_(self%id_amm,amm_flux)
   _SET_BOTTOM_EXCHANGE_(self%id_nit,nit_flux)

   ! Set sink and source terms for the benthos (change per surface area per second)
   ! Note that this must include the fluxes to and from the pelagic.
   !_SET_ODE_BEN_(self%id_ben_amm,-amm_flux/secs_pr_day)

   ! Also store sediment flux as diagnostic variable.
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_amm,amm_flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sed_nit,nit_flux)

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

END SUBROUTINE aed_nitrogen_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fnitrif(use_oxy,Rnitrif,Knitrif,theta_nitrif,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nitrification
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rnitrif,Knitrif,theta_nitrif,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fnitrif = Rnitrif * oxy/(Knitrif+oxy) * (theta_nitrif**(temp-20.0))
   ELSE
      fnitrif = Rnitrif * (theta_nitrif**(temp-20.0))
   ENDIF

END FUNCTION fnitrif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE AED_REAL FUNCTION fdenit(use_oxy,Rdenit,Kdenit,theta_denit,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for denitrification
!
! Here, the classical Michaelis-Menten formulation for denitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL,INTENT(in)  :: use_oxy
   AED_REAL,INTENT(in) :: Rdenit,Kdenit,theta_denit,oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (use_oxy) THEN
      fdenit = Rdenit * Kdenit/(Kdenit+oxy) * (theta_denit**(temp-20.0))
   ELSE
      fdenit = Rdenit * (theta_denit**(temp-20.0))
   ENDIF

END FUNCTION fdenit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_nitrogen

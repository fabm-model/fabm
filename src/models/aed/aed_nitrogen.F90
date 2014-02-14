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
   USE fabm_types

   IMPLICIT NONE

   PRIVATE    ! By default make everything private
!
   PUBLIC type_aed_nitrogen
!
   TYPE,extends(type_base_model) :: type_aed_nitrogen
!     Variable identifiers
      type (type_state_variable_id)      :: id_nit, id_amm !nitrate & ammonium
      type (type_state_variable_id)      :: id_oxy,id_denit_product
      type (type_dependency_id)          :: id_temp
      type (type_bottom_state_variable_id)      :: id_Fsed_amm,id_Fsed_nit
      type (type_diagnostic_variable_id) :: id_nitrif,id_denit
      type (type_horizontal_diagnostic_variable_id) :: id_sed_amm,id_sed_nit
      type (type_conserved_quantity_id)  :: id_totN

!     Model parameters
      real(rk) :: Rnitrif,Rdenit,Fsed_amm,Fsed_nit,Knitrif,Kdenit,Ksed_amm,Ksed_nit, &
                          theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit
      LOGICAL  :: use_oxy,use_no2,use_sed_model

      CONTAINS   ! Model Procedures
        procedure :: initialize               => aed_nitrogen_init
        procedure :: do                       => aed_nitrogen_do
        procedure :: do_ppdd                  => aed_nitrogen_do_ppdd
        procedure :: do_benthos               => aed_nitrogen_do_benthos
        procedure :: get_conserved_quantities => aed_nitrogen_get_conserved_quantities
   END TYPE


!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_nitrogen_init(self,configunit)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with FABM.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (type_aed_nitrogen),TARGET,INTENT(INOUT) :: self
   INTEGER,INTENT(in)                             :: configunit

   real(rk)          :: nit_initial=4.5
   real(rk)          :: amm_initial=4.5
   real(rk)          :: Rnitrif = 0.01
   real(rk)          :: Rdenit = 0.01
   real(rk)          :: Fsed_amm = 3.5
   real(rk)          :: Fsed_nit = 3.5
   real(rk)          :: Knitrif = 150.0
   real(rk)          :: Kdenit = 150.0
   real(rk)          :: Ksed_amm = 30.0
   real(rk)          :: Ksed_nit = 30.0
   real(rk)          :: theta_nitrif = 1.0
   real(rk)          :: theta_denit = 1.0
   real(rk)          :: theta_sed_amm = 1.0
   real(rk)          :: theta_sed_nit = 1.0
   CHARACTER(len=64) :: nitrif_reactant_variable=''
   CHARACTER(len=64) :: denit_product_variable=''
   CHARACTER(len=64) :: Fsed_amm_variable=''
   CHARACTER(len=64) :: Fsed_nit_variable=''


   real(rk), parameter :: secs_pr_day = 86400.
   NAMELIST /aed_nitrogen/ nit_initial,amm_initial,Rnitrif,Rdenit,Fsed_amm,Fsed_nit, &
                    Knitrif,Kdenit,Ksed_amm,Ksed_nit,                     &
                    theta_nitrif,theta_denit,theta_sed_amm,theta_sed_nit, &
                    nitrif_reactant_variable,denit_product_variable,      &
                    Fsed_amm_variable, Fsed_nit_variable
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Read the namelist
   read(configunit,nml=aed_nitrogen,err=99)

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
   call self%register_state_variable(self%id_amm,'amm','mmol/m**3','ammonium',            &
                                    amm_initial,minimum=0.0_rk,no_river_dilution=.true.)
   call self%register_state_variable(self%id_nit,'nit','mmol/m**3','nitrate',             &
                                    nit_initial,minimum=0.0_rk,no_river_dilution=.true.)
   ! Register external state variable dependencies
   self%use_oxy = nitrif_reactant_variable .NE. '' !This means oxygen module switched on
   IF (self%use_oxy) THEN
     call self%register_state_dependency(self%id_oxy,nitrif_reactant_variable)
   ENDIF
   self%use_no2 = denit_product_variable .NE. '' !This means n2 module switched on
   IF (self%use_no2) call self%register_state_dependency(self%id_denit_product,denit_product_variable)

   self%use_sed_model = Fsed_amm_variable .NE. ''
   IF (self%use_sed_model) THEN
     call self%register_bottom_state_dependency(self%id_Fsed_amm,Fsed_amm_variable)
     call self%register_bottom_state_dependency(self%id_Fsed_nit,Fsed_nit_variable)
   ENDIF

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_nitrif,'nitrif','mmol/m**3/d',       &
                                                         'Nitrification rate',       &
                                        time_treatment=time_treatment_step_integrated)
   call self%register_diagnostic_variable(self%id_denit,'denit','mmol/m**3/d',         &
                                                         'De-nitrification rate',    &
                                        time_treatment=time_treatment_step_integrated)
   call self%register_horizontal_diagnostic_variable(self%id_sed_amm,'sed_amm','mmol/m**2/d',      &
                                                         'Ammonium sediment flux',   &
                        time_treatment=time_treatment_step_integrated)
   call self%register_horizontal_diagnostic_variable(self%id_sed_nit,'sed_nit','mmol/m**2/d',      &
                                                         'Nitrate sediment flux',    &
                        time_treatment=time_treatment_step_integrated)

   ! Register conserved quantities
   call self%register_conserved_quantity(self%id_totN,'TN','mmol/m**3','Total nitrogen')

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)

   RETURN

99 CALL self%fatal_error('aed_nitrogen_init','Error reading namelist aed_nitrogen')

END SUBROUTINE aed_nitrogen_init
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_nitrogen_do(self,_FABM_ARGS_DO_RHS_)
!-------------------------------------------------------------------------------
! Right hand sides of aed_nitrogen model
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_nitrogen),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_RHS_
!
!LOCALS
   real(rk)           :: amm,nit,oxy,temp !State variables
   real(rk)           :: nitrification,denitrification
   real(rk),PARAMETER :: secs_pr_day = 86400.
   real(rk),PARAMETER :: Yoxy_nitrif = 3. !ratio of oxygen to nitrogen utilised during nitrification
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

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
   nitrification = fnitrif(self,oxy,temp)
   denitrification = fdenit(self,oxy,temp)

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
   _FABM_LOOP_END_

END SUBROUTINE aed_nitrogen_do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_nitrogen_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
!-------------------------------------------------------------------------------
! Right hand sides of nitrogen biogeochemical model exporting
! production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_nitrogen),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_PPDD_
!
!LOCALS
   real(rk)           :: amm,nit,oxy,temp !State variables
   real(rk)           :: nitrification,denitrification
   real(rk),PARAMETER :: secs_pr_day = 86400.
   real(rk),PARAMETER :: Yoxy_nitrif = 3. !ratio of oxygen to nitrogen utilised during nitrification
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_


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
   nitrification = fnitrif(self,oxy,temp)
   denitrification = fdenit(self,oxy,temp)

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
   _FABM_LOOP_END_

END SUBROUTINE aed_nitrogen_do_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_nitrogen_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
!-------------------------------------------------------------------------------
! Calculate pelagic bottom fluxes and benthic sink and source terms of AED nitrogen.
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_nitrogen),INTENT(in) :: self
   _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
!
!LOCALS
   ! Environment
   real(rk) :: temp !, layer_ht

   ! State
   real(rk) :: amm,nit,oxy

   ! Temporary variables
   real(rk) :: amm_flux,nit_flux
   real(rk) :: Fsed_amm, Fsed_nit

   ! Parameters
   real(rk),PARAMETER :: secs_pr_day = 86400.
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current environmental conditions for the bottom pelagic layer.
   _GET_(self%id_temp,temp)  ! local temperature

    ! Retrieve current (local) state variable values.
   _GET_(self%id_amm,amm) ! ammonium
   _GET_(self%id_nit,nit) ! nitrate

   IF (self%use_sed_model) THEN
       _GET_HORIZONTAL_(self%id_Fsed_amm,Fsed_amm)
       _GET_HORIZONTAL_(self%id_Fsed_nit,Fsed_nit)
!print *,'Fsed_amm = ',Fsed_amm,' Fsed_nit = ',Fsed_nit
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
   _FABM_HORIZONTAL_LOOP_END_

END SUBROUTINE aed_nitrogen_do_benthos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_nitrogen_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
!-------------------------------------------------------------------------------
! Get the total of conserved quantities (currently only nitrogen)
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_nitrogen),INTENT(in) :: self
   _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
!
!LOCALS
   real(rk)                     :: amm,nit
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_amm,amm) ! ammonium
   _GET_(self%id_nit,nit) ! nitrate

   ! Total nutrient is simply the sum of all variables.
   _SET_CONSERVED_QUANTITY_(self%id_totN,amm+nit)

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

END SUBROUTINE aed_nitrogen_get_conserved_quantities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE real(rk) FUNCTION fnitrif(self,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for nitrification
!
! Here, the classical Michaelis-Menten formulation for nitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_nitrogen),INTENT(in) :: self
   real(rk),INTENT(in)                 :: oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (self%use_oxy) THEN
      fnitrif = self%Rnitrif * oxy/(self%Knitrif+oxy) * (self%theta_nitrif**(temp-20.0))
   ELSE
      fnitrif = self%Rnitrif * (self%theta_nitrif**(temp-20.0))
   ENDIF

END FUNCTION fnitrif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE real(rk) FUNCTION fdenit(self,oxy,temp)
!-------------------------------------------------------------------------------
! Michaelis-Menten formulation for denitrification
!
! Here, the classical Michaelis-Menten formulation for denitrification
! is formulated.
!-------------------------------------------------------------------------------
!ARGUMENTS
   class (type_aed_nitrogen),INTENT(in) :: self
   real(rk),INTENT(in)                    :: oxy,temp
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (self%use_oxy) THEN
      fdenit = self%Rdenit * self%Kdenit/(self%Kdenit+oxy) * (self%theta_denit**(temp-20.0))
   ELSE
      fdenit = self%Rdenit * (self%theta_denit**(temp-20.0))
   ENDIF

END FUNCTION fdenit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_nitrogen

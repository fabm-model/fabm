#include "fabm_driver.h"

!svn $Id: npzd_Franks.h 995 2020-01-10 04:01:28Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2020 The ROMS/TOMS Group        Craig V. Lewis   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  Nutrient-Phytoplankton-Zooplankton-Detritus Model.                  !
!                                                                      !
!  This routine computes the biological sources and sinks and adds     !
!  then the global biological fields.                                  !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Franks et al, 1986: Behavior of simple plankton model with        !
!      food-level acclimation by herbivores, Marine Biology, 91,       !
!      121-129.                                                        !
!                                                                      !
!  Adapted from code written originally by Craig V. Lewis.             !
!  Adapted for the FABM by Phil Wallhead (NIVA) 2020                   !
!                                                                      !
!***********************************************************************

module niva_roms_npzd_Franks

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_niva_roms_npzd_Franks
      ! Variable identifiers
      type (type_state_variable_id)      :: id_NO3_, id_Phyt, id_Zoop, id_SDet
      type (type_dependency_id)          :: id_z_r

      ! Model parameters
      real(rk) :: K_ext, K_NO3, K_Phy, Vm_NO3, PhyMR
      real(rk) :: ZooGR, ZooMR, ZooMD, ZooGA, ZooEC, DetRR
   contains
      procedure :: initialize
      procedure :: do
   end type
   real(rk), parameter :: eps = 1.0e-16_rk

contains

   subroutine initialize(self, configunit)
      class (type_niva_roms_npzd_Franks), intent(inout), target :: self
      integer,                intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk / 86400.0_rk
      real(rk) :: wDet

      ! Store parameter values in our own derived type
      call self%get_parameter(self%K_ext,'K_ext','1/m','light extinction coefficient',default=0.067_rk)
      call self%get_parameter(self%K_NO3,'K_NO3','mmol N/m3','half-saturation concentration for phytoplankton nitrate uptake',default=1.0_rk)
      !NB: This is erroneously listed as an inverse half-saturation concentration in npzd_Franks.h.
      call self%get_parameter(self%K_Phy,'K_Phy','mmol N/m3','phytoplankton saturation coefficient',default=0.4_rk)
      call self%get_parameter(self%Vm_NO3,'Vm_NO3','1/day','nitrate uptake rate',default=1.5_rk)
      call self%get_parameter(self%PhyMR,'PhyMR','1/day','phytoplankton senescence/mortality rate',default=0.1_rk)
      call self%get_parameter(self%ZooGR,'ZooGR','1/day','zooplankton maximum growth rate',default=0.52_rk)
      call self%get_parameter(self%ZooMR,'ZooMR','1/day','zooplankton mortality rate',default=0.145_rk)
      call self%get_parameter(self%ZooMD,'ZooMD','1/day','zooplankton death bits rate',default=0.05_rk)
      call self%get_parameter(self%ZooGA,'ZooGA','-','zooplankton grazing inefficiency',default=0.3_rk)
      call self%get_parameter(self%ZooEC,'ZooEC','-','zooplankton excreted fraction',default=0.15_rk)
      call self%get_parameter(self%DetRR,'DetRR','1/day','detritus remineralization rate',default=0.1_rk)
      call self%get_parameter(wDet,'wDet','m/day','detrital sinking rate',default=8.0_rk,scale_factor=d_per_s)

      ! Register state variables
      call self%register_state_variable(self%id_NO3_,'NO3_','mmol m-3','nutrients',initial_value=1.67_rk,minimum=eps,vertical_movement=0.0_rk)
      call self%register_state_variable(self%id_Phyt,'Phyt','mmol m-3','phytoplankton',initial_value=0.08_rk,minimum=eps,vertical_movement=0.0_rk)
      call self%register_state_variable(self%id_Zoop,'Zoop','mmol m-3','zooplankton',initial_value=0.06_rk,minimum=eps,vertical_movement=0.0_rk)
      call self%register_state_variable(self%id_SDet,'SDet','mmol m-3','detritus',initial_value=0.04_rk,minimum=eps,vertical_movement=-wDet)

      ! Register the contribution of all state variables to total nitrogen
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_NO3_)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_Phyt)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_Zoop)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_SDet)

      ! Register environmental dependencies
      call self%register_dependency(self%id_z_r, standard_variables%depth)

   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_niva_roms_npzd_Franks), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: cff, cff1, cff2, cff3, dtdays
      real(rk)            :: NO3_, Phyt, Zoop, SDet, z_r
      real(rk)            :: NO3__old, Phyt_old, Zoop_old, SDet_old 
      real(rk), parameter :: dtbio = 200.0_rk ! "Standard" time step in seconds
      ! ROMS-FABM should give very similar results to ROMS when dt(ng) = dtbio
      real(rk), dimension(4) :: v1
      integer :: itrmx(1), itrc

      dtdays = dtbio / 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_NO3_,NO3_) ! nutrient
         _GET_(self%id_Phyt,Phyt) ! phytoplankton
         _GET_(self%id_Zoop,Zoop) ! zooplankton
         _GET_(self%id_SDet,SDet) ! detritus

         _GET_(self%id_z_r,z_r) ! depth of layer midpoints (m)

         ! Store input values (Bio_old in ROMS code)
         NO3__old = NO3_
         Phyt_old = Phyt
         Zoop_old = Zoop
         SDet_old = SDet

         ! Determine Correction for negativity.
         cff1 = MAX(0.0_rk,eps-NO3_)+                                   &
     &          MAX(0.0_rk,eps-Phyt)+                                   &
     &          MAX(0.0_rk,eps-Zoop)+                                   &
     &          MAX(0.0_rk,eps-SDet)

         ! If correction needed, determine the largest pool to debit.
         if (cff1.gt.0.0) then
            v1 = (/NO3_, Phyt, Zoop, SDet/)

            !This is the ROMS way of doing it
            itrmx(1)=1
            cff = v1(1)
            DO itrc=1,4
              IF (v1(itrc).gt.cff) THEN
                itrmx(1) = itrc
                cff = v1(itrc)
              END IF
            END DO
            !Update new values.
            DO itrc=1,4
              v1(itrc) = MAX(eps,v1(itrc)) -                            &
     &                   cff1*(SIGN(0.5_rk, REAL(itrmx(1)-itrc,rk)**2)+ &
     &                         SIGN(0.5_rk,-REAL(itrmx(1)-itrc,rk)**2))
            END DO

            ! ! Alternative approach using MAXLOC
            ! ! (perhaps more elegant and not slower, tested 02/11/2020)
            ! cff = MAXVAL(v1)
            ! itrmx = MAXLOC(v1)
            ! ! Update new values
            ! v1 = MAX(eps,v1)
            ! v1(itrmx(1)) = v1(itrmx(1)) - cff1

            NO3_ = v1(1)
            Phyt = v1(2)
            Zoop = v1(3)
            SDet = v1(4)
         end if

         !Note: We could insert here the BioIter loop as in npzd_Franks.h
         !      However, we have not bothered to do this, since most applications
         !      will use BioIter = 1, and anyway the iterations will not extend
         !      here over the vertical sinking terms, as they do in npzd_Franks.h.
         !ITER_LOOP: DO Iter=1,self%BioIter

         ! Nutrient uptake by phytoplankton.
         cff1 = dtdays*self%Vm_NO3
         cff  = Phyt * cff1*EXP(self%K_ext*z_r) / (self%K_NO3+NO3_)
         NO3_ = NO3_ / (1.0_rk+cff)
         Phyt = Phyt + NO3_*cff

         ! Phytoplankton grazing by Zooplankton and mortality to Detritus
         ! (rate: PhyMR).
         cff1 = dtdays*self%ZooGR
         cff2 = dtdays*self%PhyMR
         cff3 = self%K_phy*self%K_phy
         cff  = Zoop*Phyt*cff1 / (cff3+Phyt*Phyt)
         Phyt = Phyt / (1.0_rk+cff+cff2)
         Zoop = Zoop + Phyt*cff*(1.0_rk-self%ZooGA)
         SDet = SDet + Phyt*(cff2+cff*(self%ZooGA-self%ZooEC))
         NO3_ = NO3_ + Phyt*cff*self%ZooEC

         ! Zooplankton excretion to nutrients and mortality to Detritus.
         cff1 = 1.0_rk/(1.0_rk+dtdays*(self%ZooMR+self%ZooMD))
         cff2 = dtdays*self%ZooMR
         cff3 = dtdays*self%ZooMD
         Zoop = Zoop*cff1
         NO3_ = NO3_ + Zoop*cff2
         SDet = SDet + Zoop*cff3

         ! Detritus breakdown to nutrients.
         cff1 = dtdays*self%DetRR
         cff2 = 1.0_rk/(1.0_rk+cff1)
         SDet = SDet*cff2
         NO3_ = NO3_ + SDet*cff1

         ! Set temporal derivatives (concentration / second)
         _ADD_SOURCE_(self%id_NO3_, (NO3_-NO3__old)/dtbio)
         _ADD_SOURCE_(self%id_Phyt, (Phyt-Phyt_old)/dtbio)
         _ADD_SOURCE_(self%id_Zoop, (Zoop-Zoop_old)/dtbio)
         _ADD_SOURCE_(self%id_SDet, (SDet-SDet_old)/dtbio)

      ! Leave spatial loops (if any)
      _LOOP_END_
   end subroutine do

end module niva_roms_npzd_Franks
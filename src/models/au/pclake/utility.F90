!*******************************************************************************
!*                                                                             *
!* basic.F90                                                                 *
!*                                                                             *
!*-----------------------------------------------------------------------------*
!*******************************************************************************
#include "fabm_driver.h"
!-------------------------------------------------------------------------------
! Module description
!-------------------------------------------------------------------------------
!  au_pclake_utility:

!  Original implemented by Fenjuan Hu(feh@dmu.dk)
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_pclake_foodweb_water
! !INTERFACE:
   module au_pclake_utility
!
! !DESCRIPTION:

! !USES:
   use fabm_types
   use fabm_expressions
   
   
   implicit none

!  default: all is private.
      private
      
      PUBLIC uFunTmAbio
      PUBLIC uFunTmBio
      PUBLIC uFunTmVeg

   contains
   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure real(rk) function uFunTmAbio(uTm,cTheta)
!
!    !DESCRIPTION:
!    The temperature function of abiotic processes is defined here
!    This function is applied for mineralization both in water and sediment
!    nitrification both in water and sediment, diffusion and sedimentation
!   
!    !INPUT PARAMETERS:

      real(rk),                         intent(in) :: uTm,cTheta
      !real(rk)              :: cTmRef=20.0_rk
!    !REVISION HISTORY:
!     Original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmAbio = cTheta**(uTm-20.0_rk)

   end function uFunTmAbio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure real(rk) function uFunTmBio(uTm,cSigTm,cTmOpt)
!
!    !DESCRIPTION:
!    The temperature function of biotic processes is defined here
!    This funcition is applied for all groups of phytoplankton both in water and sediment
!    zooplankton in water column, fish and zoobenthos
!   
!    !INPUT PARAMETERS:
      real(rk),                         intent(in) :: uTm,cSigTm,cTmOpt
      !real(rk)              :: cTmRef=20.0_rk
!   
!    !REVISION HISTORY:
!     Original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmBio = exp(-0.5_rk/cSigTm**2 *((uTm-cTmOpt)**2- (20.0_rk-cTmOpt)**2))

   end function uFunTmBio
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ivlev formulation for zooplankton grazing on phytoplankton
!
! !INTERFACE:
   pure real(rk) function uFunTmVeg(uTm,cQ10Veg)
!
!    !DESCRIPTION:
!    The temperature function of vegetation processes is defined here
!    This function is applied for macrophytes 
!   
!    !INPUT PARAMETERS:
      real(rk),                         intent(in) :: uTm,cQ10Veg
      !real(rk)              :: cTmRef=20.0_rk
!   
!    !REVISION HISTORY:
!     Original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmVeg = ((cQ10Veg )** (0.1_rk * (uTm - 20.0_rk)))

   end function uFunTmVeg

   end module au_pclake_utility

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

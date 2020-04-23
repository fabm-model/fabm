#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_ergom_split_utilities - ERGOM diatom compartment
!
! !INTERFACE:
   module fabm_uhh_ergom_split_utilities
!
! !DESCRIPTION:
! This module holds common functions and data types for the splitted ergom model
!
! !USES:
   use fabm_types

   implicit none
   
   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Michaelis-Menten formulation for nutrient uptake
!
! !INTERFACE:
   pure real(rk) function yy(a,x)
!
! !DESCRIPTION:
! This is a squared Michaelis-Menten type of limiter:
! \begin{equation}\label{Y}
! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
! \end{equation}
!
! !USES:
  IMPLICIT NONE
!
! !INPUT PARAMETERS:
  real(rk),intent(in)        :: a,x
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------
!BOC
  yy=x**2/(a**2+x**2)
  RETURN
  END function yy
!EOC
   
!EOC


   end module fabm_uhh_ergom_split_utilities

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

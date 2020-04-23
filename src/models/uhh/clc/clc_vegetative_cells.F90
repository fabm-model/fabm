#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_clcveg - CLC cyanobacteria vegetative cells
!
! !INTERFACE:
   module fabm_uhh_clcveg
!
! !DESCRIPTION:
!
! Model of the cyanobacteria vegetative cells lifestage within a cyanobacteria lifecycle.
!
! !USES:
   use fabm_uhh_clcbase
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_uhh_clcbase),public :: type_uhh_clcveg
      contains
      procedure :: init_lifestage
      procedure :: calculate_lifecycle_flux
   end type
!EOP
!-----------------------------------------------------------------------

   contains

   subroutine init_lifestage(self)
   class(type_uhh_clcveg) :: self

   self%lifestage_name = 'vegetative_cells'
   self%next_in_cycle  = 'heterocysts'

   end subroutine init_lifestage


   subroutine calculate_lifecycle_flux(self,e,q,c_flux,s_flux,g_flux)
     class(type_uhh_clcveg) :: self
     real(rk) :: e,q,c_flux,s_flux,g_flux
     if (self%q_min > 0.0_rk)  self%Sflux_per_Cflux=0.5_rk
     
     c_flux = self%tscale*max((self%q_min-self%trange-q)/self%trange,0.0_rk)
     s_flux = self%Sflux_per_Cflux * c_flux
     g_flux = self%Gflux_per_Cflux * c_flux
   end subroutine


   end module fabm_uhh_clcveg


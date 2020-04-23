#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_clcrec - CLC cyanobacteria recruitive cells
!
! !INTERFACE:
   module fabm_uhh_clcrec
!
! !DESCRIPTION:
!
! Model of the cyanobacteria recruitive cells lifestage within a cyanobacteria lifecycle.
!
! !USES:
   use fabm_uhh_clcbase
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_uhh_clcbase),public :: type_uhh_clcrec
      contains
      procedure :: init_lifestage
      procedure :: calculate_lifecycle_flux
   end type
!EOP
!-----------------------------------------------------------------------

   contains

   subroutine init_lifestage(self)
   class(type_uhh_clcrec) :: self

   self%lifestage_name = 'recruitive_cells'
   self%next_in_cycle  = 'vegetative_cells'

   end subroutine init_lifestage


   subroutine calculate_lifecycle_flux(self,e,q,c_flux,s_flux,g_flux)
     class(type_uhh_clcrec) :: self
     real(rk) :: e,q,c_flux,s_flux,g_flux
     if (self%e_max < 1.0_rk)     self%Gflux_per_Cflux=1.5_rk

     c_flux = self%tscale*max((min(1.0_rk,e)-self%e_max+self%trange)/self%trange,0.0_rk)
     s_flux = self%Sflux_per_Cflux * c_flux
     g_flux = self%Gflux_per_Cflux * c_flux
   end subroutine


   end module fabm_uhh_clcrec


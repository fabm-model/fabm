!$Id$

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Allocate memory for character strings
!
! !INTERFACE:
   subroutine bio_alloc_info
!
! !DESCRIPTION:
!  Allocates memory for the character strings describing the 
!  variables of the bio model
!
! !USES:
   use bio_var, only: numc
   use bio_var, only: var_ids,var_names,var_units,var_long

   IMPLICIT NONE


! !REVISION HISTORY:!
!  Original author(s): Lars Umlauf, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!
! !LOCAL VARIABLES:
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC

   allocate(var_ids(numc),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating var_ids)'

   allocate(var_names(numc),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating var_names)'

   allocate(var_units(numc),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating var_units)'

   allocate(var_long(numc),stat=rc)
   if (rc /= 0) stop 'allocate_memory(): Error allocating var_long)'


   return
   end subroutine bio_alloc_info
!EOC
!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FABM model library [Fortran 2003 models only]
!
! !INTERFACE:
   module fabm_library

!
! !USES:
   ! FABM modules
   use fabm_types
   
#ifdef _FABM_F2003_
   ! Specific biogeochemical models
   use fabm_examples_npzd_f2003
   ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
#endif

   implicit none
!   
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public fabm_library_create_model
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Function returning specific biogeochemical model given a model name.
!
! !INTERFACE:
   function fabm_library_create_model(name) result(model)
!
! !INPUT PARAMETERS:
      character(*),intent(in) :: name
      _CLASS_ (type_model_info),pointer :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
      nullify(model)
      
#ifdef _FABM_F2003_
      select case (name)
         case ('examples_npzd_f2003'); allocate(type_examples_npzd_f2003 :: model)
         ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
      end select
#endif
   
      if (associated(model)) call init_model_info(model)

   end function fabm_library_create_model
!EOC

!-----------------------------------------------------------------------

   end module fabm_library
   
!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

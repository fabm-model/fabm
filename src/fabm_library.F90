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
   use fabm_bb_passive
   use fabm_examples_npzd_f2003
   ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
   use aed_models
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
   function fabm_library_create_model(modelname,instancename,parent,configunit) result(model)
!
! !INPUT PARAMETERS:
      character(*),intent(in)           :: modelname,instancename
      integer,     intent(in)           :: configunit
      _CLASS_ (type_model_info),target  :: parent
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
      select case (modelname)
         case ('bb_passive')
            model => bb_passive_create(configunit,instancename,parent)
         case ('examples_npzd_f2003')
            model => examples_npzd_f2003_create(configunit,instancename,parent)
         ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
         case default
            if ( modelname(1:4) .eq. 'aed_' ) &
               model => aed_create_model(configunit,modelname,instancename,parent);
      end select
#endif

   end function fabm_library_create_model
!EOC

!-----------------------------------------------------------------------

   end module fabm_library

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

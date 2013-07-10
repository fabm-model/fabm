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
   use fabm_examples_npzd_nut
   use fabm_examples_npzd_phy
   use fabm_examples_npzd_zoo
   use fabm_examples_npzd_det
   use fabm_examples_duplicator
   use fabm_examples_npzd_f2003
   use aed_models
   use fabm_prey_pred
!   use pclake_models
   ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
#endif

   implicit none
!
!  default: all is private.
   private

#ifdef _FABM_F2003_
   type,extends(type_abstract_model_factory),public :: type_model_factory
      contains
      procedure,nopass :: create => fabm_library_create_model
   end type
#endif
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
         case ('bb_passive');          allocate(type_bb_passive::model)
         case ('examples_npzd_nut');   allocate(type_examples_npzd_nut::model)
         case ('examples_npzd_phy');   allocate(type_examples_npzd_phy::model)
         case ('examples_npzd_zoo');   allocate(type_examples_npzd_zoo::model)
         case ('examples_npzd_det');   allocate(type_examples_npzd_det::model)
         case ('examples_duplicator'); allocate(type_examples_duplicator::model)
         case ('examples_npzd_f2003'); allocate(type_examples_npzd_f2003::model)
         case ('fabm_prey_pred'); allocate(type_fabm_prey_pred::model)
         ! ADD_NEW_FORTRAN2003_MODEL_HERE - required
         case default
             if ( modelname(1:4) .eq. 'aed_' ) &
               model => aed_create_model(configunit,modelname,instancename,parent);
!             if ( modelname(1:7) .eq. 'pclake_' ) &
!             model => pclake_create_model(modelname);
!               model => pclake_create_model(configunit,modelname,instancename,parent);
!            if ( modelname(1:4) .eq. 'aed_' ) &
!               model => aed_create_model(modelname);
!!       keep PCLake model code update into API0.91 changes
!            if ( modelname(1:7) .eq. 'pclake_' ) &
!               model => pclake_create_model(modelname);
      end select

      if (.not.associated(model)) return

      ! If the model has not been initialized, do so now.
      ! This is the default - simulaneously creating and initializing the model is now deprecated.
      if (.not.associated(model%parent)) call parent%add_child(model,configunit,instancename)
#endif

   end function fabm_library_create_model
!EOC

!-----------------------------------------------------------------------

   end module fabm_library

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

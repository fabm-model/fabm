#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FABM model library
!
! !INTERFACE:
   module fabm_library
!
! !USES:
   use fabm_types, only: type_base_model_factory, type_base_model, factory

   use fabm_builtin_models
   use aed_models
   use bb_model_library
   use examples_model_library
   use gotm_model_library
   use iow_model_library
   use niva_model_library
   use fabm_metu_mnemiopsis
   use fabm_pml_carbonate
   use au_prey_predator
   use fabm_klimacampus_phy_feedback
   use fabm_hzg_omexdia_p
   use fabm_msi_ergom1
   ! Add additional external modules containing models or model factories here

   implicit none

   private

   public :: fabm_create_model_factory

   type,extends(type_base_model_factory) :: type_model_factory
      contains
      procedure :: create
   end type

!EOP
!-----------------------------------------------------------------------

   contains

   subroutine fabm_create_model_factory()
      if (.not.associated(factory)) then
         allocate(type_model_factory::factory)

         call factory%add(builtin_factory)
         call factory%add(aed_model_factory)
         call factory%add(bb_model_factory)
         call factory%add(examples_model_factory)
         call factory%add(gotm_model_factory)
         call factory%add(iow_model_factory)
         call factory%add(niva_model_factory)
         ! Add new additional model factories here

      end if
   end subroutine
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Subroutine returning a specific biogeochemical model given a model name.
!
! !INTERFACE:
   subroutine create(self,name,model)
!
! !INPUT PARAMETERS:
      class (type_model_factory),intent(in) :: self
      character(*),              intent(in) :: name
      class (type_base_model),pointer       :: model
!
!EOP
!-----------------------------------------------------------------------
!BOC
      nullify(model)

      select case (name)
         case ('au_prey_predator');          allocate(type_au_prey_predator::model)
         case ('metu_mnemiopsis');           allocate(type_metu_mnemiopsis::model)
         case ('pml_carbonate');             allocate(type_pml_carbonate::model)
         case ('klimacampus_phy_feedback');  allocate(type_klimacampus_phy_feedback::model)
         case ('hzg_omexdia_p');             allocate(type_hzg_omexdia_p::model)
         case ('msi_ergom1');                allocate(type_msi_ergom1::model)
         ! Add additional individual models here
         case default
            call self%type_base_model_factory%create(name,model)
      end select

      ! Store name that was used to create this model, so we can re-create it in the future.
      if (associated(model)) model%type_name = name
   end subroutine
!EOC

!-----------------------------------------------------------------------

   end module fabm_library

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

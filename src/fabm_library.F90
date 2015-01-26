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
   use fabm_types, only: type_base_model_factory, type_base_model

   use fabm_builtin_models

!  Add additional external modules containing models or model factories here
!  please in alphabetically order
   use aed_models
   use au_prey_predator
   use bb_model_library
   use examples_model_library
   use fabm_metu_mnemiopsis
   use fabm_klimacampus_phy_feedback
   use fabm_hzg_omexdia_p
   use fabm_msi_ergom1
   use gotm_model_library
   use iow_model_library
   use niva_model_library
   use pclake_model_library
   use pml_model_library
   use au_pclake_model_library

   implicit none

   private

   type,extends(type_base_model_factory),public :: type_model_factory
      contains
      procedure :: create
      procedure :: initialize
   end type

!EOP
!-----------------------------------------------------------------------

contains

   subroutine initialize(self)
      class (type_model_factory), intent(inout) :: self

      call self%add(builtin_factory)
      call self%add(aed_model_factory)
      call self%add(bb_model_factory,'bb')
      call self%add(examples_model_factory,'examples')
      call self%add(gotm_model_factory,'gotm')
      call self%add(iow_model_factory)
      call self%add(niva_model_factory)
      call self%add(pclake_model_factory)
      call self%add(pml_model_factory,'pml')
      call self%add(au_pclake_model_factory)
      ! Add additional model factories here

      ! Go through default initializaton steps. This also allows new added child model factories to initialize.
      call self%type_base_model_factory%initialize()
   end subroutine initialize

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
         case ('klimacampus_phy_feedback');  allocate(type_klimacampus_phy_feedback::model)
         case ('hzg_omexdia_p');             allocate(type_hzg_omexdia_p::model)
         case ('msi_ergom1');                allocate(type_msi_ergom1::model)
         ! Add additional individual models here
         case default
            call self%type_base_model_factory%create(name,model)
      end select

      ! Store name that was used to create this model, so we can re-create it in the future.
      if (associated(model)) model%type_name = name
   end subroutine create
!EOC

!-----------------------------------------------------------------------

   end module fabm_library

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

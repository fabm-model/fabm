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

   use aed_models
   use niva_model_library
   use fabm_bb_passive
   use fabm_bb_filter_feeder
   use fabm_examples_npzd_nut
   use fabm_examples_npzd_phy
   use fabm_examples_npzd_zoo
   use fabm_examples_npzd_det
   use fabm_examples_duplicator
   use fabm_examples_npzd_f2003
   use fabm_examples_benthic_predator
   use fabm_examples_mean
   use fabm_gotm_npzd
   use fabm_gotm_fasham
   use fabm_gotm_ergom
   use fabm_metu_mnemiopsis
   use fabm_pml_carbonate
   use au_prey_predator
   use fabm_klimacampus_phy_feedback
   use fabm_hzg_omexdia_p
   use fabm_iow_spm
   use fabm_iow_age
   use fabm_iow_ergom
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

         call factory%add(aed_model_factory)
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
         case ('bb_passive');                allocate(type_bb_passive::model)
         case ('bb_filter_feeder');          allocate(type_bb_filter_feeder::model)
         case ('examples_npzd_nut');         allocate(type_examples_npzd_nut::model)
         case ('examples_npzd_phy');         allocate(type_examples_npzd_phy::model)
         case ('examples_npzd_zoo');         allocate(type_examples_npzd_zoo::model)
         case ('examples_npzd_det');         allocate(type_examples_npzd_det::model)
         case ('examples_duplicator');       allocate(type_examples_duplicator::model)
         case ('examples_npzd_f2003');       allocate(type_examples_npzd_f2003::model)
         case ('examples_benthic_predator'); allocate(type_examples_benthic_predator::model)
         case ('examples_mean');             allocate(type_examples_mean::model)
         case ('gotm_npzd');                 allocate(type_gotm_npzd::model)
         case ('gotm_fasham');               allocate(type_gotm_fasham::model)
         case ('gotm_ergom');                allocate(type_gotm_ergom::model)
         case ('metu_mnemiopsis');           allocate(type_metu_mnemiopsis::model)
         case ('pml_carbonate');             allocate(type_pml_carbonate::model)
         case ('klimacampus_phy_feedback');  allocate(type_klimacampus_phy_feedback::model)
         case ('hzg_omexdia_p');             allocate(type_hzg_omexdia_p::model)
         case ('iow_spm');                   allocate(type_iow_spm::model)
         case ('iow_age');                   allocate(type_iow_age::model)
         case ('iow_ergom');                 allocate(type_iow_ergom::model)
         case ('msi_ergom1');                allocate(type_msi_ergom1::model)
         ! Add additional individual models here
         case default
            call self%type_base_model_factory%create(name,model)
      end select

   end subroutine
!EOC

!-----------------------------------------------------------------------

   end module fabm_library

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

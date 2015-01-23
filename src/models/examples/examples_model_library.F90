module examples_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: examples_model_factory

contains

   subroutine create(self,name,model)

      use examples_benthic_predator
      use examples_duplicator
      use examples_mean
      use examples_npzd
      use examples_npzd_f2003
      use nonlocal

      ! Add new examples models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('benthic_predator');      allocate(type_examples_benthic_predator::model)
         case ('duplicator');            allocate(type_examples_duplicator::model)
         case ('mean');                  allocate(type_examples_mean::model)
         case ('npzd_nut');              allocate(type_examples_npzd_nut::model)
         case ('npzd_phy');              allocate(type_examples_npzd_phy::model)
         case ('npzd_zoo');              allocate(type_examples_npzd_zoo::model)
         case ('npzd_det');              allocate(type_examples_npzd_det::model)
         case ('npzd_f2003');            allocate(type_examples_npzd_f2003::model)
         case ('depth_integral');        allocate(type_depth_integral::model)
         case ('vertical_distribution'); allocate(type_vertical_distribution::model)
         ! Add new examples models here
      end select

   end subroutine

end module

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

      ! Add new examples models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('examples_benthic_predator'); allocate(type_examples_benthic_predator::model)
         case ('examples_duplicator');       allocate(type_examples_duplicator::model)
         case ('examples_mean');             allocate(type_examples_mean::model)
         case ('examples_npzd_nut');         allocate(type_examples_npzd_nut::model)
         case ('examples_npzd_phy');         allocate(type_examples_npzd_phy::model)
         case ('examples_npzd_zoo');         allocate(type_examples_npzd_zoo::model)
         case ('examples_npzd_det');         allocate(type_examples_npzd_det::model)
         case ('examples_npzd_f2003');       allocate(type_examples_npzd_f2003::model)
         ! Add new examples models here
      end select

   end subroutine

end module

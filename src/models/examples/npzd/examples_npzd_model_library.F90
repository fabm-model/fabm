module examples_npzd_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use examples_npzd_det
   use examples_npzd_nut
   use examples_npzd_phy
   use examples_npzd_zoo
   ! Add new NPZD modules here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: examples_npzd_model_factory

contains

   subroutine create(self, name, model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('nut'); allocate(type_examples_npzd_nut::model)
         case ('phy'); allocate(type_examples_npzd_phy::model)
         case ('zoo'); allocate(type_examples_npzd_zoo::model)
         case ('det'); allocate(type_examples_npzd_det::model)
         ! Add new NPZD models here
      end select
   end subroutine create

end module examples_npzd_model_library

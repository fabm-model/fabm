module msi_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: msi_model_factory

contains

   subroutine create(self,name,model)

      use fabm_msi_ergom1
      ! Add new MSI modules here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('ergom1'); allocate(type_msi_ergom1::model)
         ! Add new MSI models here
      end select

   end subroutine

end module

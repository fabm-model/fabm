module pml_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: pml_model_factory

contains

   subroutine create(self,name,model)

      use fabm_pml_carbonate   
      ! Add new PML models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('carbonate'); allocate(type_pml_carbonate::model)
         ! Add new PML models here
      end select

   end subroutine

end module

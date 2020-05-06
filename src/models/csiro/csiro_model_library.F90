module csiro_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: csiro_model_factory

contains

   subroutine create(self,name,model)

      use csiro_seagrass
      ! Add new CSIRO models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('seagrass'); allocate(type_csiro_seagrass::model)
         ! Add new CSIRO models here
      end select

   end subroutine

end module

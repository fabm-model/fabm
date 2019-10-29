module iow_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: iow_model_factory

contains

   subroutine create(self,name,model)

      use fabm_iow_age
      use fabm_iow_ergom
      use fabm_iow_spm
      ! Add new IOW modules here

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('age');   allocate(type_iow_age::model)
         case ('ergom'); allocate(type_iow_ergom::model)
         case ('spm');   allocate(type_iow_spm::model)
         ! Add new IOW models here
      end select

   end subroutine

end module

module gotm_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: gotm_model_factory

contains

   subroutine create(self, name, model)

      use gotm_npzd
      use gotm_fasham
      use gotm_ergom
      use gotm_light
      ! Add new GOTM models here

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('npzd');   allocate(type_gotm_npzd::model)
         case ('fasham'); allocate(type_gotm_fasham::model)
         case ('ergom');  allocate(type_gotm_ergom::model)
         case ('light');  allocate(type_gotm_light::model)
         ! Add new GOTM models here
      end select

   end subroutine

end module

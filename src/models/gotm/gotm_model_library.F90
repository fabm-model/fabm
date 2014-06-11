module gotm_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: gotm_model_factory

contains

   subroutine create(self,name,model)

      use gotm_npzd
      use gotm_fasham
      use gotm_ergom
      ! Add new GOTM models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('gotm_npzd');    allocate(type_gotm_npzd::model)
         case ('gotm_fasham');  allocate(type_gotm_fasham::model)
         case ('gotm_ergom');   allocate(type_gotm_ergom::model)
         ! Add new GOTM models here
      end select

   end subroutine

end module

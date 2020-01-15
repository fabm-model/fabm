module su_ersem_model_library

   ! Library with biogeochemical models developed by Swansea University, UK.
   ! Copyright (C) 2019 - Suzana Leles

   use fabm_types, only: type_base_model_factory,type_base_model

   use su_ersem_mesozooplankton

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: su_ersem_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('mesozooplankton'); allocate(type_su_ersem_mesozooplankton::model)
         ! Add new models here
         case default
            call self%type_base_model_factory%create(name,model)
      end select
   end subroutine create

end module

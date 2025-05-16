module python_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use wrapped_python_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: python_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('wrapped_python_model'); allocate(type_wrapped_python_model::model)
         case default
            call self%type_base_model_factory%create(name,model)
      end select

   end subroutine

end module

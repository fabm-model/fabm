module jrc_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use jrc_bsem
   use jrc_med_ergom

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
!KB      procedure :: initialize
      procedure :: create
   end type

   type (type_factory),save,target,public :: jrc_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('bsem'); allocate(type_jrc_bsem::model)
         case ('med_ergom');    allocate(type_jrc_med_ergom::model)
         case default
            call self%type_base_model_factory%create(name,model)
      end select

   end subroutine

end module

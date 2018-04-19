module au_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use au_prey_predator
! Add new AU models here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: initialize
      procedure :: create
   end type

   type (type_factory),save,target,public :: au_model_factory

contains

   subroutine initialize(self)
      class (type_factory), intent(inout) :: self

      ! Add additional model factories here

      ! Go through default initializaton steps.
      ! This also allows new added child model factories to initialize.
      call self%type_base_model_factory%initialize()
   end subroutine initialize

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('lotka_volterra'); allocate(type_au_pp_lotka_volterra::model)
         case ('jacob_monod');    allocate(type_au_pp_jacob_monod::model)
         case default
            call self%type_base_model_factory%create(name,model)
!           Add new AU models here
      end select

   end subroutine

end module

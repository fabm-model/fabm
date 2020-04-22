module examples_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use examples_benthic_predator
   use examples_duplicator
   use examples_mean
   use examples_npzd_model_library
   use nonlocal
   use examples_light_cycle
   ! Add new examples modules here

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: initialize
      procedure :: create
   end type

   type (type_factory), save, target, public :: examples_model_factory

contains

   subroutine initialize(self)
      class (type_factory), intent(inout) :: self
      call self%add(examples_npzd_model_factory,'npzd')
      ! Add additional child model factories here

      ! Go through default initializaton steps.
      ! This also allows newly added child model factories to initialize.
      call self%type_base_model_factory%initialize()
   end subroutine initialize

   subroutine create(self, name, model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('benthic_predator');      allocate(type_examples_benthic_predator::model)
         case ('duplicator');            allocate(type_examples_duplicator::model)
         case ('mean');                  allocate(type_examples_mean::model)
         case ('depth_integral');        allocate(type_depth_integral::model)
         case ('vertical_distribution'); allocate(type_vertical_distribution::model)
         case ('light_cycle');           allocate(type_examples_light_cycle::model)
         ! Add individual example models here
         case default
            call self%type_base_model_factory%create(name, model)
      end select

   end subroutine

end module

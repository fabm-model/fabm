module examples_particle_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use examples_particle_passive_particle
   use examples_particle_depth_integrated_predator
   ! Add new example particle modules here

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: examples_particle_model_factory

contains

   subroutine create(self, name, model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
      case ('particle'); allocate(type_passive_particle::model)
      !case ('migration'); allocate(type_templates_migration::model)
      case ('depth_integrated_predator'); allocate(type_depth_integrated_predator::model)
      case default
         call self%type_base_model_factory%create(name, model)
      end select

   end subroutine

end module

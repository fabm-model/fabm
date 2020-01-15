module su_model_library

   ! Library with biogeochemical models developed by Swansea University, UK.
   ! Copyright (C) 2019 - Suzana Leles

   use fabm_types, only: type_base_model_factory,type_base_model

   use su_mixo
   use su_pools
   use su_algae
   use su_bacteria
   use su_npz
   use su_light_cycle
#ifdef SU_ERSEM
   use su_ersem_model_library
#endif

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
   contains
      procedure :: initialize
      procedure :: create
   end type

   type (type_factory), save, target, public :: su_model_factory

contains

   subroutine initialize(self)
      class (type_factory), intent(inout) :: self
#ifdef SU_ERSEM
      call self%add(su_ersem_model_factory, 'ersem')
#endif
      ! Add additional child model factories here

      ! Go through default initializaton steps.
      ! This also allows newly added child model factories to initialize.
      call self%type_base_model_factory%initialize()
   end subroutine initialize

   subroutine create(self,name,model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('mixo'); allocate(type_su_mixo::model)
         case ('pools'); allocate(type_su_pools::model)
         case ('algae'); allocate(type_su_algae::model)
         case ('bacteria'); allocate(type_su_bacteria::model)
         case ('npz'); allocate(type_su_npz::model)
         case ('light_cycle'); allocate(type_su_light_cycle::model)
         ! Add new models here
         case default
            call self%type_base_model_factory%create(name,model)
      end select
   end subroutine create

end module

module bb_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: bb_model_factory

contains

   subroutine create(self, name, model)

      use bb_filter_feeder
      use bb_lorenz63
      use bb_passive
      ! Add new BB models here

      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('passive');       allocate(type_bb_passive::model)
         case ('lorenz63');      allocate(type_bb_lorenz63::model)
         case ('filter_feeder'); allocate(type_bb_filter_feeder::model)
         ! Add new BB models here
      end select

   end subroutine

end module

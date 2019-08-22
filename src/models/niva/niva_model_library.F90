module niva_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: niva_model_factory

contains

   subroutine create(self,name,model)

      use fabm_niva_oxydep
      use fabm_niva_brom_bio
      use fabm_niva_brom_carb
      use fabm_niva_brom_eqconst
      use fabm_niva_brom_redox
      use fabm_niva_brom_salt
      ! Add new NIVA models here
      use niva_domcast

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('oxydep');       allocate(type_niva_oxydep::model)
         case ('brom_bio');     allocate(type_niva_brom_bio::model)
         case ('brom_carb');    allocate(type_niva_brom_carb::model)
         case ('brom_eqconst'); allocate(type_niva_brom_eqconst::model)
         case ('brom_redox');   allocate(type_niva_brom_redox::model)
         case ('brom_salt');    allocate(type_niva_brom_salt::model)
         case ('domcast');      allocate(type_niva_domcast::model)
         ! Add new NIVA models here
      end select

   end subroutine

end module

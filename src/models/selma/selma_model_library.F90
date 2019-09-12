module selma_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use selma
   use selma_phytoplankton
   use selma_zooplankton

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
!KB      procedure :: initialize
      procedure :: create
   end type

   type (type_factory),save,target,public :: selma_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('selma'); allocate(type_selma::model)
         case ('phytoplankton'); allocate(type_selma_phytoplankton::model)
         case ('zooplankton'); allocate(type_selma_zooplankton::model)
         case default
            call self%type_base_model_factory%create(name,model)
      end select

   end subroutine

end module

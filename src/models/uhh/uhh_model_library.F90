module uhh_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: uhh_model_factory

contains

   subroutine create(self,name,model)

      use fabm_uhh_phy_feedback
      use fabm_uhh_ergom_split_base
      use fabm_uhh_ergom_split_zoo
      !use fabm_uhh_clc
      use fabm_uhh_clchet
      use fabm_uhh_clcaki
      use fabm_uhh_clcrec
      use fabm_uhh_clcveg
      use fabm_uhh_dinoflag
      use fabm_uhh_diatoms
      use fabm_uhh_halogen
      use fabm_uhh_uv
      ! Add new UHH models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('phy_feedback');  allocate(type_uhh_phy_feedback::model)
         case ('ergom_split_base');      allocate(type_uhh_ergom_split_base::model)
         case ('ergom_split_zoo');       allocate(type_uhh_ergom_split_zoo::model)
         !case ('uhh_clc');      allocate(type_uhh_clc::model)
         case ('clchet');   allocate(type_uhh_clchet::model)
         case ('clcaki');   allocate(type_uhh_clcaki::model)
         case ('clcrec');   allocate(type_uhh_clcrec::model)
         case ('clcveg');   allocate(type_uhh_clcveg::model)
         case ('dinoflag'); allocate(type_uhh_dinoflag::model)
         case ('diatoms');  allocate(type_uhh_diatoms::model)
         case ('halogen');  allocate(type_uhh_halogen::model)
         case ('uv');       allocate(type_uhh_uv::model)
         ! Add new UHH models here
      end select

   end subroutine

end module

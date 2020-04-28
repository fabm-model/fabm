module akvaplan_model_library

   ! Library with biogeochemical models developed by Akvaplan-niva (http://www.akvaplan.niva.no).
   ! Copyright (C) 2016 - Akvaplan-niva

   use fabm_types, only: type_base_model_factory,type_base_model

   use akvaplan_tracer
   use akvaplan_plume_injection
   use akvaplan_tracer_sed
   use akvaplan_antiparasitic

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: akvaplan_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('tracer');          allocate(type_tracer::model)
         case ('plume_injection'); allocate(type_plume_injection::model)
         case ('tracer_sed');      allocate(type_tracer_sed::model)
         case ('antiparasitic');   allocate(type_antiparasitic::model)
         ! Add new models here
         case default
            call self%type_base_model_factory%create(name,model)
      end select
   end subroutine create

end module

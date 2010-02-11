#include"cppdefs.h"

module bio_driver

implicit none

public register_variable_dependency,getvar,getbiovar

private

integer,parameter :: not_found_id = -1
integer,parameter :: temp_id = 1, salt_id = 2, par_id = 3, I_0_id = 4, pres_id = 5, wind_id = 6

contains

function register_variable_dependency(variable) result(id)
   character(len=*),intent(in) :: variable
   integer :: id
   
   select case (variable)
      case ('temp')
         id = temp_id
      case ('salt')
         id = salt_id
      case ('par')
         id = par_id
      case ('I_0')
         id = I_0_id
      case ('pres')
         id = pres_id
      case ('wind')
         id = wind_id
      case default
         write (*,*) 'WARNING: unable to locate variable '//trim(variable)
         id = not_found_id
   end select
end function register_variable_dependency

function getbiovar(id,LOCATIONVARIABLE) result(value)
   use bio_var,only: cc

   integer,     intent(in) :: id
   LOCATIONTYPE,intent(in) :: LOCATIONVARIABLE
   REALTYPE                :: value
   
   value = cc(id,LOCATIONVARIABLE)
end function getbiovar

function getvar(id,LOCATIONVARIABLE) result(value)
   use meanflow, only: T,S,z
   use bio_var,  only: par
   use airsea,   only: I_0,wind=>w
   
   integer,     intent(in) :: id
   LOCATIONTYPE,intent(in) :: LOCATIONVARIABLE
   REALTYPE :: value

   select case (id)
      case (temp_id)
         value = T(LOCATIONVARIABLE)
      case (salt_id)
         value = S(LOCATIONVARIABLE)
      case (par_id)
         value = par(LOCATIONVARIABLE)
      case (pres_id)
         value = -z(LOCATIONVARIABLE)/10.d0
      case (I_0_id)
         value = I_0
      case (wind_id)
         value = wind
      case (not_found_id)
         stop 'Request for non-existing variable'
         value = _ZERO_
      case default
         
   end select   
end function getvar

!function register_state_variable_dependency(modelinfo,module,name,instance) result(id)
!   type (type_model_info),intent(inout) :: modelinfo
!   character(len=*), intent(in) :: module,name
!   integer,          intent(in) :: instance
!   integer                      :: id
!   
!      type (type_variable_dependency) :: dependencies_old(modelinfo%variable_dependency_count)
!      
!      dependencies_old = modelinfo%variable_dependencies
!      if (allocated(modelinfo%variable_dependencies)) deallocate(modelinfo%variable_dependencies)
!      allocate(modelinfo%variable_dependencies(modelinfo%variable_dependency_count+1))
!      
!      modelinfo%variable_dependencies(1:modelinfo%variable_dependency_count) = dependencies_old
!      modelinfo%variable_dependency_count = modelinfo%variable_dependency_count+1
!      modelinfo%variable_dependencies(modelinfo%variable_dependency_count)%module = module
!      modelinfo%variable_dependencies(modelinfo%variable_dependency_count)%name = name
!      modelinfo%variable_dependencies(modelinfo%variable_dependency_count)%instance = instance
!   
!   id = -1
!   
!end function register_state_variable_dependency

end module bio_driver
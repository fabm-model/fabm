module fabm_c_variable

   use iso_c_binding, only: c_double, c_int, c_char, c_f_pointer, c_loc, c_ptr

   use fabm_types
   use fabm_c_helper

   implicit none

contains

   subroutine variable_get_metadata(pvariable,length,name,units,long_name) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_metadata
      type (c_ptr),          intent(in), value             :: pvariable
      integer(c_int),        intent(in), value             :: length
      character(kind=c_char),intent(out),dimension(length) :: name,units,long_name

      type (type_internal_variable),pointer :: variable

      call c_f_pointer(pvariable, variable)
      call copy_to_c_string(variable%name,     name)
      call copy_to_c_string(variable%units,    units)
      call copy_to_c_string(variable%long_name,long_name)
   end subroutine variable_get_metadata

   subroutine variable_get_output_name(pvariable,length,name) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_output_name
      type (c_ptr),          intent(in), value  :: pvariable
      integer(c_int),        intent(in), value  :: length
      character(kind=c_char),intent(out),dimension(length) :: name

      type (type_internal_variable),pointer :: variable
      class (type_base_model),      pointer :: owner
      character(len=attribute_length)       :: name_

      call c_f_pointer(pvariable, variable)
      name_ = get_safe_name(variable%name)
      call copy_to_c_string(name_,name)
   end subroutine variable_get_output_name

   subroutine variable_get_long_path(pvariable,length,long_name) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_long_path
      type (c_ptr),          intent(in), value  :: pvariable
      integer(c_int),        intent(in), value  :: length
      character(kind=c_char),intent(out),dimension(length) ::long_name

      type (type_internal_variable),pointer :: variable
      class (type_base_model),       pointer :: owner
      character(len=attribute_length)        :: long_name_

      call c_f_pointer(pvariable, variable)
      long_name_ = variable%long_name
      owner => variable%owner
      do while (associated(owner%parent))
         long_name_ = trim(owner%long_name)//'/'//trim(long_name_)
         owner => owner%parent
      end do
      call copy_to_c_string(long_name_,long_name)
   end subroutine variable_get_long_path

   function variable_get_background_value(pvariable) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_background_value
      type (c_ptr), value, intent(in) :: pvariable
      real(kind=c_double)             :: value

      type (type_internal_variable),pointer :: variable

      call c_f_pointer(pvariable, variable)
      value = 0.0_rk
      if (size(variable%background_values%pointers)>0) value = variable%background_values%pointers(1)%p
   end function variable_get_background_value

   function variable_get_output(pvariable) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_output
      type (c_ptr), value, intent(in) :: pvariable
      integer(kind=c_int)             :: value

      type (type_internal_variable),pointer :: variable

      call c_f_pointer(pvariable, variable)
      value = logical2int(variable%output/=output_none)
   end function variable_get_output

   function variable_get_real_property(pvariable,name,default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_real_property
      type (c_ptr), value,          intent(in) :: pvariable
      character(kind=c_char),target,intent(in) :: name(*)
      real(kind=c_double), value,   intent(in) :: default
      real(kind=c_double)                      :: value

      type (type_internal_variable),  pointer :: variable
      character(len=attribute_length),pointer :: pname

      call c_f_pointer(pvariable, variable)
      call c_f_pointer(c_loc(name), pname)
      value = variable%properties%get_real(pname(:index(pname,C_NULL_CHAR)-1),default=default)
   end function variable_get_real_property

end module fabm_c_variable
module fabm_c_variable

   use iso_c_binding, only: c_int, c_char, c_f_pointer, c_loc, c_ptr, c_null_ptr

   use fabm_types
   use fabm_c_helper
   use fabm_properties, only: type_property, type_integer_property, type_real_property, type_logical_property, type_string_property

   implicit none

   type type_standard_variable_wrapper
      class (type_base_standard_variable), pointer :: p
   end type

contains

   subroutine variable_get_metadata(pvariable, length, name, units, long_name) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_metadata
      type (c_ptr),           intent(in), value              :: pvariable
      integer(c_int),         intent(in), value              :: length
      character(kind=c_char), intent(out), dimension(length) :: name, units, long_name

      type (type_internal_variable), pointer :: variable

      call c_f_pointer(pvariable, variable)
      call copy_to_c_string(variable%name,      name)
      call copy_to_c_string(variable%units,     units)
      call copy_to_c_string(variable%long_name, long_name)
   end subroutine variable_get_metadata

   subroutine variable_get_long_path(pvariable, length, long_name) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_long_path
      type (c_ptr),           intent(in), value              :: pvariable
      integer(c_int),         intent(in), value              :: length
      character(kind=c_char), intent(out), dimension(length) ::long_name

      type (type_internal_variable), pointer :: variable
      class (type_base_model),       pointer :: owner
      character(len=attribute_length)        :: long_name_

      call c_f_pointer(pvariable, variable)
      long_name_ = variable%long_name
      owner => variable%owner
      do while (associated(owner%parent))
         long_name_ = trim(owner%long_name) // '/' // trim(long_name_)
         owner => owner%parent
      end do
      call copy_to_c_string(long_name_,long_name)
   end subroutine variable_get_long_path

   function variable_get_background_value(pvariable) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_background_value
      type (c_ptr), value, intent(in) :: pvariable
      real(kind=rke)                  :: value

      type (type_internal_variable), pointer :: variable

      call c_f_pointer(pvariable, variable)
      value = 0.0_rk
      if (size(variable%background_values%pointers) > 0) value = variable%background_values%pointers(1)%p
   end function variable_get_background_value

   function variable_get_missing_value(pvariable) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_missing_value
      type (c_ptr), value, intent(in) :: pvariable
      real(kind=rke)                  :: value

      type (type_internal_variable), pointer :: variable

      call c_f_pointer(pvariable, variable)
      value = variable%missing_value
   end function variable_get_missing_value

   function variable_get_output(pvariable) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_output
      type (c_ptr), value, intent(in) :: pvariable
      integer(kind=c_int)             :: value

      type (type_internal_variable), pointer :: variable

      call c_f_pointer(pvariable, variable)
      value = logical2int(variable%output /= output_none)
   end function variable_get_output

   function variable_is_required(pvariable) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_is_required
      type (c_ptr), value, intent(in) :: pvariable
      integer(kind=c_int)             :: value

      type (type_internal_variable), pointer :: variable

      call c_f_pointer(pvariable, variable)
      value = logical2int(variable%presence /= presence_external_optional)
   end function variable_is_required

   function variable_get_no_river_dilution(pvariable) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_no_river_dilution
      type (c_ptr), value, intent(in) :: pvariable
      integer(kind=c_int)             :: value

      type (type_internal_variable), pointer :: variable

      call c_f_pointer(pvariable, variable)
      value = logical2int(variable%no_river_dilution)
   end function variable_get_no_river_dilution

   function variable_get_no_precipitation_dilution(pvariable) bind(c) result(value)
   !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_no_precipitation_dilution  ! reduced indentation to work around flang bug
      type (c_ptr), value, intent(in) :: pvariable
      integer(kind=c_int)             :: value

      type (type_internal_variable), pointer :: variable

      call c_f_pointer(pvariable, variable)
      value = logical2int(variable%no_precipitation_dilution)
   end function variable_get_no_precipitation_dilution

   function variable_get_property_type(pvariable, name) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_property_type
      type (c_ptr), value,            intent(in) :: pvariable
      character(kind=c_char), target, intent(in) :: name(*)
      integer(kind=c_int)                        :: value

      type (type_internal_variable),   pointer :: variable
      character(len=attribute_length), pointer :: pname
      class (type_property),           pointer :: property

      value = typecode_unknown
      call c_f_pointer(pvariable, variable)
      call c_f_pointer(c_loc(name), pname)
      property => variable%properties%get_property(pname(:index(pname, C_NULL_CHAR) - 1))
      if (associated(property)) then
         select type (property)
         class is (type_real_property)
            value = typecode_real
         class is (type_integer_property)
            value = typecode_integer
         class is (type_logical_property)
            value = typecode_logical
         class is (type_string_property)
            value = typecode_string
         end select
      end if
   end function variable_get_property_type

   function variable_get_real_property(pvariable, name, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_real_property
      type (c_ptr), value,            intent(in) :: pvariable
      character(kind=c_char), target, intent(in) :: name(*)
      real(kind=rke), value,          intent(in) :: default
      real(kind=rke)                             :: value

      type (type_internal_variable),   pointer :: variable
      character(len=attribute_length), pointer :: pname

      call c_f_pointer(pvariable, variable)
      call c_f_pointer(c_loc(name), pname)
      value = variable%properties%get_real(pname(:index(pname, C_NULL_CHAR) - 1), default=default)
   end function variable_get_real_property

   function variable_get_integer_property(pvariable, name, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_integer_property
      type (c_ptr), value,            intent(in) :: pvariable
      character(kind=c_char), target, intent(in) :: name(*)
      integer(kind=c_int), value,     intent(in) :: default
      integer(kind=c_int)                        :: value

      type (type_internal_variable),   pointer :: variable
      character(len=attribute_length), pointer :: pname

      call c_f_pointer(pvariable, variable)
      call c_f_pointer(c_loc(name), pname)
      value = variable%properties%get_integer(pname(:index(pname, C_NULL_CHAR) - 1), default=default)
   end function variable_get_integer_property

   function variable_get_logical_property(pvariable, name, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_logical_property
      type (c_ptr), value,            intent(in) :: pvariable
      character(kind=c_char), target, intent(in) :: name(*)
      integer(kind=c_int), value,     intent(in) :: default
      integer(kind=c_int)                        :: value

      type (type_internal_variable),   pointer :: variable
      character(len=attribute_length), pointer :: pname

      call c_f_pointer(pvariable, variable)
      call c_f_pointer(c_loc(name), pname)
      value = logical2int(variable%properties%get_logical(pname(:index(pname, C_NULL_CHAR) - 1), default=int2logical(default)))
   end function variable_get_logical_property

   function find_standard_variable(name) bind(c) result(pvariable)
      !DIR$ ATTRIBUTES DLLEXPORT :: find_standard_variable
      character(kind=c_char), target, intent(in) :: name(*)
      type (c_ptr) :: pvariable

      character(len=attribute_length), pointer :: pname
      type (type_standard_variable_wrapper), pointer :: wrapper

      call c_f_pointer(c_loc(name), pname)
      allocate(wrapper)
      wrapper%p => standard_variables%find(pname(:index(pname, C_NULL_CHAR) - 1))
      if (associated(wrapper%p)) then
         pvariable = c_loc(wrapper)
      else
         deallocate(wrapper)
         pvariable = c_null_ptr
      end if
   end function

end module fabm_c_variable

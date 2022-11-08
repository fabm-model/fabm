module fabm_c_parameter

   use iso_c_binding, only: c_ptr, c_loc, c_char, c_int, c_double, C_NULL_CHAR

   use fabm_c_helper
   use fabm_types, only: attribute_length
   use fabm_driver, only: driver
   use fabm_c
   use yaml_settings, only: type_scalar_value, format_real, format_integer
   use yaml_types, only: type_yaml_key_value_pair => type_key_value_pair, type_yaml_dictionary => type_dictionary, &
      type_yaml_error => type_error

   implicit none

contains

   subroutine reset_parameter(pmodel, index) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: reset_parameter
      type (c_ptr), value,   intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index

      type (type_model_wrapper),    pointer :: model
      class (type_scalar_value),    pointer :: scalar_value
      type (type_yaml_key_value_pair),   pointer :: pair, previous_pair
      class (type_yaml_dictionary), pointer :: parent_node

      call c_f_pointer(pmodel, model)

      scalar_value => get_dictionary_entry(model%p%settings, 'parameters', index)
      if (.not. associated(scalar_value%backing_store_node)) return

      select type (n => scalar_value%parent%backing_store_node)
      class is (type_yaml_dictionary)
         parent_node => n
      end select
      pair => parent_node%first
      previous_pair => null()
      do while (associated(pair))
         if (associated(pair%value, scalar_value%backing_store_node)) then
            if (associated(previous_pair)) then
               previous_pair%next => pair%next
            else
               parent_node%first => pair%next
            end if
            call pair%value%finalize()
            deallocate(pair%value)
            deallocate(pair)
            exit
         end if
         previous_pair => pair
         pair => pair%next
      end do
   
      ! Re-initialize the model using updated parameter values
      call reinitialize(model)
   end subroutine reset_parameter

   subroutine set_parameter(pmodel, name, value)
      type (c_ptr),                   intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: name(*)
      character(len=*),               intent(in) :: value

      type (type_model_wrapper),       pointer :: model
      character(len=attribute_length), pointer :: pname
      integer :: islash, n
      class (type_yaml_dictionary),    pointer :: instances, instance, parameters
      type(type_yaml_error),           pointer :: yaml_error

      call c_f_pointer(pmodel, model)
      call c_f_pointer(c_loc(name), pname)

      n = index(pname, C_NULL_CHAR) - 1
      islash = index(pname, '/')

      instances => model%p%settings%backing_store%get_dictionary('instances',required=.true.,error=yaml_error)
      instance => instances%get_dictionary(pname(1:islash),required=.true.,error=yaml_error)
      parameters => instances%get_dictionary('parameters',required=.true.,error=yaml_error)
      call parameters%set_string(pname(islash+1:n), value)

      ! Re-initialize the model using updated parameter values
      call reinitialize(model)
   end subroutine

   subroutine set_real_parameter(pmodel, name, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_real_parameter
      type (c_ptr),   value,          intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: name(*)
      real(c_double), value,          intent(in) :: value

      call set_parameter(pmodel, name, format_real(value))
   end subroutine set_real_parameter

   function get_real_parameter(pmodel, index, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_real_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default
      real(c_double)                    :: value

      type (type_model_wrapper), pointer :: model
      class (type_scalar_value), pointer :: scalar_value

      call c_f_pointer(pmodel, model)
      scalar_value => get_dictionary_entry(model%p%settings, 'parameters', index)
      select type (scalar_value)
      class is (type_real_setting)
         if (int2logical(default)) then
            value = scalar_value%default
         else
            value = scalar_value%pvalue
         end if
      class default
         call driver%fatal_error('get_real_parameter', 'not a real variable')
      end select
   end function get_real_parameter

   subroutine set_integer_parameter(pmodel, name, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_integer_parameter
      type (c_ptr),   value,          intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: name(*)
      integer(c_int), value,          intent(in) :: value

      call set_parameter(pmodel, name, format_integer(value))
   end subroutine set_integer_parameter

   function get_integer_parameter(pmodel, index, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_integer_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default
      integer(c_int)                    :: value

      type (type_model_wrapper), pointer :: model
      class (type_scalar_value), pointer :: scalar_value

      call c_f_pointer(pmodel, model)
      scalar_value => get_dictionary_entry(model%p%settings, 'parameters', index)
      select type (scalar_value)
      class is (type_integer_setting)
         if (int2logical(default)) then
            value = scalar_value%default
         else
            value = scalar_value%pvalue
         end if
      class default
         call driver%fatal_error('get_integer_parameter', 'not an integer variable')
      end select
   end function get_integer_parameter

   subroutine set_logical_parameter(pmodel, name, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_logical_parameter
      type (c_ptr),  value,           intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: name(*)
      integer(c_int), value,          intent(in) :: value

      if (int2logical(value)) then
         call set_parameter(pmodel, name, 'true')
      else
         call set_parameter(pmodel, name, 'false')
      end if
   end subroutine set_logical_parameter

   function get_logical_parameter(pmodel, index, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_logical_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default
      integer(c_int)                    :: value

      type (type_model_wrapper), pointer :: model
      class (type_scalar_value), pointer :: scalar_value

      call c_f_pointer(pmodel, model)
      scalar_value => get_dictionary_entry(model%p%settings, 'parameters', index)
      select type (scalar_value)
      class is (type_logical_setting)
         if (int2logical(default)) then
            value = scalar_value%default
         else
            value = scalar_value%pvalue
         end if
      class default
         call driver%fatal_error('get_logical_parameter', 'not a logical variable')
      end select
   end function get_logical_parameter

   subroutine set_string_parameter(pmodel, name, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_string_parameter
      type (c_ptr), value,            intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: name(*), value(*)

      character(len=attribute_length), pointer :: pvalue

      call c_f_pointer(c_loc(value), pvalue)
      call set_parameter(pmodel, name, pvalue(:index(pvalue, C_NULL_CHAR) - 1))
   end subroutine set_string_parameter

   subroutine get_string_parameter(pmodel, index, default, length, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_string_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default, length
      character(kind=c_char)            :: value(length)

      type (type_model_wrapper), pointer :: model
      class (type_scalar_value), pointer :: scalar_value

      call c_f_pointer(pmodel, model)
      scalar_value => get_dictionary_entry(model%p%settings, 'parameters', index)
      select type (scalar_value)
      class is (type_string_setting)
         if (int2logical(default)) then
            call copy_to_c_string(scalar_value%default, value)
         else
            call copy_to_c_string(scalar_value%pvalue, value)
         end if
      class default
         call driver%fatal_error('get_string_parameter', 'not a string variable')
      end select
   end subroutine get_string_parameter

end module

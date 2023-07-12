module fabm_c_parameter

   use iso_c_binding, only: c_ptr, c_loc, c_char, c_int, C_NULL_CHAR

   use fabm_c_helper
   use fabm_types, only: attribute_length, rke
   use fabm_driver, only: driver
   use fabm_c
   use yaml_settings, only: type_scalar_value, format_real, format_integer
   use yaml_types, only: type_yaml_key_value_pair => type_key_value_pair, type_yaml_dictionary => type_dictionary, &
      type_yaml_error => type_error, type_yaml_node => type_node

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

      scalar_value => get_parameter_by_index(model%p%root, index)
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
      integer                                  :: islash, n
      class (type_yaml_dictionary),    pointer :: instances, instance, parameters

      call c_f_pointer(pmodel, model)
      call c_f_pointer(c_loc(name), pname)

      n = index(pname, C_NULL_CHAR) - 1
      islash = index(pname, '/')

      instances => get_child_dictionary(model%p%settings%backing_store, 'instances')
      instance => get_child_dictionary(instances, pname(1:islash - 1))
      parameters => get_child_dictionary(instance, 'parameters')
      call parameters%set_string(pname(islash+1:n), value)

      ! Re-initialize the model using updated parameter values
      call reinitialize(model)

   contains

      function get_child_dictionary(parent, key) result(child)
         class (type_yaml_dictionary), intent(inout) :: parent
         character(len=*),             intent(in)    :: key
         class (type_yaml_dictionary), pointer :: child

         type(type_yaml_error), pointer :: yaml_error
         class(type_yamL_node), pointer :: child_node

         child => parent%get_dictionary(key,.false.,error=yaml_error)
         if (associated(yaml_error)) then
            call driver%fatal_error('set_parameter', trim(yaml_error%message))
            return
         end if
         if (.not. associated(child)) then
            allocate(child)
            call child%set_path(trim(parent%path) // '/' // trim(key))
            child_node => child
            call parent%set(key, child_node)
         end if
      end function

   end subroutine

   subroutine set_real_parameter(pmodel, name, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_real_parameter
      type (c_ptr),   value,          intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: name(*)
      real(rke), value,               intent(in) :: value

      call set_parameter(pmodel, name, format_real(value))
   end subroutine set_real_parameter

   function get_real_parameter(pmodel, index, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_real_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default
      real(rke)                         :: value

      type (type_model_wrapper), pointer :: model
      class (type_scalar_value), pointer :: scalar_value

      call c_f_pointer(pmodel, model)
      scalar_value => get_parameter_by_index(model%p%root, index)
      select type (scalar_value)
      class is (type_real_setting)
         if (int2logical(default)) then
            value = scalar_value%default
         else
            value = scalar_value%pvalue / scalar_value%scale_factor
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
      scalar_value => get_parameter_by_index(model%p%root, index)
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
      scalar_value => get_parameter_by_index(model%p%root, index)
      select type (scalar_value)
      class is (type_logical_setting)
         if (int2logical(default)) then
            value = logical2int(scalar_value%default)
         else
            value = logical2int(scalar_value%pvalue)
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
      scalar_value => get_parameter_by_index(model%p%root, index)
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

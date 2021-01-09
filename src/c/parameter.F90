module fabm_c_parameter

   use iso_c_binding, only: c_ptr, c_loc, c_char, c_int, c_double, C_NULL_CHAR

   use fabm_c_helper
   use fabm_properties
   use fabm_types, only: attribute_length
   use fabm_driver, only: driver
   use fabm_c

   implicit none

contains

   subroutine reset_parameter(pmodel, index) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: reset_parameter
      type (c_ptr), value,   intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index
      class (type_property), pointer    :: property

      type (type_model_wrapper), pointer :: model

      call c_f_pointer(pmodel, model)
      property => model%p%root%parameters%get_property(index)
      if (.not. associated(property)) return
      call model%forced_parameters%delete(property%name)

      ! Re-initialize the model using updated parameter values
      call reinitialize(model)
   end subroutine reset_parameter

   subroutine set_real_parameter(pmodel, name, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_real_parameter
      type (c_ptr),   value,          intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: name(*)
      real(c_double), value,          intent(in) :: value

      type (type_model_wrapper),       pointer :: model
      character(len=attribute_length), pointer :: pname

      call c_f_pointer(pmodel, model)
      call c_f_pointer(c_loc(name), pname)
      call model%forced_parameters%set_real(pname(:index(pname, C_NULL_CHAR) - 1), value)

      ! Re-initialize the model using updated parameter values
      call reinitialize(model)
   end subroutine set_real_parameter

   function get_real_parameter(pmodel, index, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_real_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default
      real(c_double)                    :: value

      type (type_model_wrapper), pointer :: model
      class (type_property),     pointer    :: property

      call c_f_pointer(pmodel, model)
      property => model%p%root%parameters%get_property(index)
      select type (property)
      class is (type_real_property)
         if (int2logical(default)) then
            value = property%default
         else
            value = property%value
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

      type (type_model_wrapper),       pointer :: model
      character(len=attribute_length), pointer :: pname

      call c_f_pointer(pmodel, model)
      call c_f_pointer(c_loc(name), pname)
      call model%forced_parameters%set_integer(pname(:index(pname, C_NULL_CHAR) - 1), value)

      ! Re-initialize the model using updated parameter values
      call reinitialize(model)
   end subroutine set_integer_parameter

   function get_integer_parameter(pmodel, index, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_integer_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default
      integer(c_int)                    :: value

      type (type_model_wrapper), pointer :: model
      class (type_property),     pointer :: property

      call c_f_pointer(pmodel, model)
      property => model%p%root%parameters%get_property(index)
      select type (property)
      class is (type_integer_property)
         if (int2logical(default)) then
            value = property%default
         else
            value = property%value
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

      type (type_model_wrapper),       pointer :: model
      character(len=attribute_length), pointer :: pname

      call c_f_pointer(pmodel, model)
      call c_f_pointer(c_loc(name), pname)
      call model%forced_parameters%set_logical(pname(:index(pname, C_NULL_CHAR) - 1), int2logical(value))

      ! Re-initialize the model using updated parameter values
      call reinitialize(model)
   end subroutine set_logical_parameter

   function get_logical_parameter(pmodel, index, default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_logical_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default
      integer(c_int)                    :: value

      type (type_model_wrapper), pointer :: model
      class (type_property),     pointer :: property

      call c_f_pointer(pmodel, model)
      property => model%p%root%parameters%get_property(index)
      select type (property)
      class is (type_logical_property)
         if (int2logical(default)) then
            value = logical2int(property%default)
         else
            value = logical2int(property%value)
         end if
      class default
         call driver%fatal_error('get_logical_parameter', 'not a logical variable')
      end select
   end function get_logical_parameter

   subroutine set_string_parameter(pmodel, name, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_string_parameter
      type (c_ptr), value,            intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: name(*), value(*)

      type (type_model_wrapper),       pointer :: model
      character(len=attribute_length), pointer :: pname, pvalue

      call c_f_pointer(pmodel, model)
      call c_f_pointer(c_loc(name), pname)
      call c_f_pointer(c_loc(value), pvalue)
      call model%forced_parameters%set_string(pname(:index(pname, C_NULL_CHAR) - 1), pvalue(:index(pname, C_NULL_CHAR) - 1))

      ! Re-initialize the model using updated parameter values
      call reinitialize(model)
   end subroutine set_string_parameter

   subroutine get_string_parameter(pmodel, index, default, length, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_string_parameter
      type (c_ptr),   value, intent(in) :: pmodel
      integer(c_int), value, intent(in) :: index, default, length
      character(kind=c_char)            :: value(length)

      type (type_model_wrapper), pointer :: model
      class (type_property),     pointer :: property

      call c_f_pointer(pmodel, model)
      property => model%p%root%parameters%get_property(index)
      select type (property)
      class is (type_string_property)
         if (int2logical(default)) then
            call copy_to_c_string(property%default, value)
         else
            call copy_to_c_string(property%value, value)
         end if
      class default
         call driver%fatal_error('get_string_parameter', 'not a string variable')
      end select
   end subroutine get_string_parameter

end module

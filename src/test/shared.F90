#include "fabm_driver.h"
#include "fabm_private.h"

module test_shared

   use fabm_driver
   use fabm_types, only: rke, get_free_unit
   use fabm

   use yaml, only: yaml_parse => parse, yaml_error_length => error_length
   use yaml_types, only: type_node, type_scalar, type_yaml_dictionary => type_dictionary, &
      type_yaml_key_value_pair => type_key_value_pair, yaml_real_kind => real_kind

   implicit none

   public type_test_driver, read_environment, deallocate_input

   private

   type, extends(type_base_driver) :: type_test_driver
   contains
      procedure :: fatal_error => test_driver_fatal_error
      procedure :: log_message => test_driver_log_message
   end type

   type type_input
      type (type_fabm_interior_variable_id)                :: interior_id
      type (type_fabm_horizontal_variable_id)              :: horizontal_id
      type (type_fabm_scalar_variable_id)                  :: scalar_id
      real(rke), allocatable _DIMENSION_GLOBAL_            :: interior_data
      real(rke), allocatable _DIMENSION_GLOBAL_HORIZONTAL_ :: horizontal_data
      real(rke)                                            :: scalar_data
      type (type_input), pointer                           :: next => null()
   end type
   type (type_input), pointer, save :: first_input => null()

contains

   subroutine read_environment(model, path _POSTARG_LOCATION_)
      class (type_fabm_model), intent(inout) :: model
      character(len=*),        intent(in)    :: path
      _DECLARE_ARGUMENTS_LOCATION_

      character(yaml_error_length)             :: yaml_error
      class (type_node),               pointer :: yaml_root
      type (type_yaml_key_value_pair), pointer :: yaml_pair
      real(rke)                                :: value
      logical                                  :: success
      type (type_input),               pointer :: input

      yaml_root => yaml_parse(path, get_free_unit(), yaml_error)
      if (yaml_error /= '') then
         call driver%log_message(yaml_error)
         stop 2
      end if

      call driver%log_message('Loading environment from ' // path)
      select type (yaml_root)
      class is (type_yaml_dictionary)
         yaml_pair => yaml_root%first
         do while (associated(yaml_pair))
            select type (node => yaml_pair%value)
            class is (type_scalar)
               value = node%to_real(0._yaml_real_kind, success)
               if (.not. success) then
                  call driver%log_message('Cannot parse ' // trim(node%string) // ' as real.')
                  stop 2
               end if
               allocate(input)
               input%interior_id = model%get_interior_variable_id(trim(yaml_pair%key))
               if (model%is_variable_used(input%interior_id)) then
                  allocate(input%interior_data _INDEX_LOCATION_)
                  input%interior_data = value
                  call model%link_interior_data(input%interior_id, input%interior_data)
               else
                  input%horizontal_id = model%get_horizontal_variable_id(trim(yaml_pair%key))
                  if (model%is_variable_used(input%horizontal_id)) then
                     allocate(input%horizontal_data _INDEX_HORIZONTAL_LOCATION_)
                     input%horizontal_data = value
                     call model%link_horizontal_data(input%horizontal_id, input%horizontal_data)
                  else
                     input%scalar_id = model%get_scalar_variable_id(trim(yaml_pair%key))
                     if (model%is_variable_used(input%scalar_id)) then
                        input%scalar_data = value
                        call model%link_scalar(input%scalar_id, input%scalar_data)
                     else
                        !call driver%log_message('WARNING: environment variable '//trim(yaml_pair%key) &
                        !   //' is not used by FABM model and will be ignored.')
                        deallocate(input)
                     end if
                  end if
               end if
               if (associated(input)) then
                  call driver%log_message('  ' // trim(yaml_pair%key) // ': ' // trim(node%string))
                  input%next => first_input
                  first_input => input
               end if
            end select
            yaml_pair => yaml_pair%next
         end do
      class default
         call driver%log_message(path // ' should contain a dictionary at root level')
         stop 2
      end select
      call yaml_root%finalize()
      deallocate(yaml_root)
   end subroutine read_environment

   subroutine deallocate_input()
      type (type_input), pointer :: next
      do while (associated(first_input))
         next => first_input%next
         deallocate(first_input)
         first_input => next
      end do
   end subroutine

   subroutine test_driver_fatal_error(self, location, message)
      class (type_test_driver), intent(inout) :: self
      character(len=*),         intent(in)    :: location, message

      write (*,'(a)') trim(location) // ': ' // trim(message)
      stop 1
   end subroutine

   subroutine test_driver_log_message(self, message)
      class (type_test_driver), intent(inout) :: self
      character(len=*),         intent(in)    :: message

      write (*,'(a)') trim(message)
   end subroutine

end module

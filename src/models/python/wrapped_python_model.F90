#include "fabm_driver.h"

module wrapped_python_model

   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_double, c_ptr, c_null_char, c_associated, c_loc, c_f_pointer
   use fabm_types

   use python_parameters

   implicit none

   private

   type, extends(type_base_model), public :: type_wrapped_python_model
      type (c_ptr) :: pobject
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

   interface
      integer(c_int) function embedded_python_initialize(program_name, home) bind(c)
       import c_char, c_int
       character(kind=c_char), dimension(*) :: program_name, home
      end function

      type(c_ptr) function embedded_python_get_model(module_name, class_name, pbase) bind(c, name='embedded_python_get_model')
       import c_char, c_ptr
       character(kind=c_char), dimension(*) :: module_name, class_name
       type(c_ptr), value :: pbase
      end function

      subroutine embedded_python_do_bottom(pbase, cache) bind(c)
       import c_ptr
       type(c_ptr), value :: pbase, cache
      end subroutine
   end interface

contains

   subroutine initialize(self, configunit)
      class(type_wrapped_python_model), intent(inout), target :: self
      integer,                          intent(in)            :: configunit

      character(len=attribute_length) :: module_name, class_name, python_home, cmd
      type(c_ptr) :: pmodule, pclass, pwhome, pwcmd
      integer(c_int) :: iresult
      type (type_base_model), pointer :: pself

      call self%get_parameter(module_name, 'module', '', 'Python module containing model class')
      call self%get_parameter(class_name, 'class', '', 'name of model class', default='Model')
      call self%get_parameter(python_home, 'home', '', 'location of the standard Python libraries', default=default_python_home)

      call get_command_argument(0, cmd)
      iresult = embedded_python_initialize(trim(cmd) // c_null_char, trim(python_home) // c_null_char)
      pself => self%type_base_model
      self%pobject = embedded_python_get_model(trim(module_name) // c_null_char, trim(class_name) // c_null_char, c_loc(pself))
      if (.not. c_associated(self%pobject)) call self%fatal_error('initialize', 'Unable to load Python model')
   end subroutine

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_wrapped_python_model), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      call embedded_python_do_bottom(self%pobject, c_loc(cache))
   end subroutine

   subroutine c_register_variable(self, name, units, long_name, domain, source, presence, initial_value, read_index, sms_index) bind(c)
      type(type_base_model),  intent(inout), target :: self
      character(kind=c_char), intent(in)            :: name(*), units(*), long_name(*)
      integer(c_int),         intent(in),    value  :: domain, source, presence
      real(c_double),         intent(in),    value  :: initial_value
      integer(c_int),         intent(inout), target :: read_index, sms_index

      type (type_link), pointer :: link, sms_link, link2
      character(len=attribute_length), pointer :: pname, punits, plong_name
      type (type_internal_variable), pointer :: variable, sms_variable

      call c_f_pointer(c_loc(name), pname)
      call c_f_pointer(c_loc(units), punits)
      call c_f_pointer(c_loc(long_name), plong_name)

      allocate(variable)
      variable%domain = domain
      link => null()
      call self%add_variable(variable, pname(:index(pname, C_NULL_CHAR) - 1), punits(:index(punits, C_NULL_CHAR) - 1), plong_name(:index(plong_name, C_NULL_CHAR) - 1), &
                             read_index=read_index, link=link, source=source, presence=presence, initial_value=initial_value)

      if (source == source_state) then
         allocate(sms_variable)
         sms_variable%domain = variable%domain
         sms_link => null()
         call self%add_variable(sms_variable, trim(link%name) // '_sms', &
            trim(variable%units) // '/s', trim(variable%long_name) // ' sources-sinks', fill_value=0.0_rk, &
            missing_value=0.0_rk, output=output_none, write_index=sms_index, link=sms_link, &
            source=source_do_bottom)
         sms_variable%write_operator = operator_add
         link2 => variable%sms_list%append(sms_variable, sms_variable%name)
         variable%sms => link2
      end if
   end subroutine

   subroutine c_unpack_horizontal_cache(cache, ni, ni_hz, nread, nread_hz, nread_scalar, nwrite_hz, read, read_hz, read_scalar, write_hz) bind(c)
      type (type_horizontal_cache), target :: cache
      integer,      intent(out) :: ni, ni_hz, nread, nread_hz, nread_scalar, nwrite_hz
      type (c_ptr), intent(out) :: read, read_hz, read_scalar, write_hz

#ifdef _INTERIOR_IS_VECTORIZED_
      ni = size(cache%read, 1)
      nread = size(cache%read, 2)
#else
      ni = 1
      nread = size(cache%read, 1)
#endif
#ifdef _HORIZONTAL_IS_VECTORIZED_
      ni_hz = size(cache%read_hz, 1)
      nread_hz = size(cache%read_hz, 2)
      nwrite_hz = size(cache%write_hz, 2)
#else
      ni_hz = 1
      nread_hz = size(cache%read_hz, 1)
      nwrite_hz = size(cache%write_hz, 1)
#endif
      nread_scalar = size(cache%read_scalar)
      read = c_loc(cache%read)
      read_hz = c_loc(cache%read_hz)
      read_scalar = c_loc(cache%read_scalar)
      write_hz= c_loc(cache%write_hz)
   end subroutine
end module

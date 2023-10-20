#include "fabm_driver.h"

! This module contains types and procedures for generic properties
! that internally can be of different types (real, integer, logical)
! It also offers a type that represents a property dictionary.

module fabm_properties

   use fabm_parameters, only: rk => rki

   implicit none

   private

   public type_property, type_property_dictionary, type_set
   public type_integer_property, type_real_property, type_logical_property, type_string_property

   integer, parameter :: metadata_string_length = 256
   integer, parameter :: value_string_length = 1024

   type, abstract :: type_property
      character(len=metadata_string_length) :: name        = ''
      character(len=metadata_string_length) :: long_name   = ''
      character(len=metadata_string_length) :: units       = ''
      character(len=metadata_string_length) :: key         = ''      ! name wih normalized case
      logical                               :: has_default = .false.
      class (type_property), pointer        :: next        => null()
   contains
      procedure :: to_real
      procedure :: to_integer
      procedure :: to_logical
      procedure :: to_string
   end type

   type, extends(type_property) :: type_integer_property
      integer :: value
      integer :: default = 0
   end type

   type, extends(type_property) :: type_real_property
      real(rk) :: value
      real(rk) :: default = 0.0_rk
   end type

   type, extends(type_property) :: type_logical_property
      logical :: value
      logical :: default = .true.
   end type

   type, extends(type_property) :: type_string_property
      character(len=value_string_length) :: value
      character(len=value_string_length) :: default= ''
   end type

   type type_property_dictionary
      class (type_property), pointer :: first => null()
   contains
      procedure :: set_property
      procedure :: set_real
      procedure :: set_integer
      procedure :: set_logical
      procedure :: set_string
      generic :: set => set_property, set_real, set_integer, set_logical, set_string

      procedure :: get_property
      procedure :: get_real
      procedure :: get_integer
      procedure :: get_logical
      procedure :: get_string

      procedure :: update

      procedure :: compare_keys

      procedure :: finalize
   end type

   type type_set_element
      character(len=metadata_string_length) :: string
      type (type_set_element), pointer :: next => null()
   end type

   type type_set
      type (type_set_element), pointer :: first => null()
   contains
      procedure :: add      => set_add
      procedure :: size     => set_size
      procedure :: to_array => set_to_array
      procedure :: finalize => set_finalize
   end type

contains

   function to_real(self, success, default) result(value)
      class (type_property), intent(in)  :: self
      logical,  optional,    intent(out) :: success
      real(rk), optional,    intent(in)  :: default
      real(rk)                           :: value

      if (present(success)) success = .true.
      select type (self)
      class is (type_real_property)
         value = self%value
         return
      class is (type_integer_property)
         value = self%value
         return
      class is (type_string_property)
         read(self%value,*,err=99,end=99) value
         return
      end select
99    if (present(success)) success = .false.
      if (present(default)) value = default
   end function

   function to_integer(self, success, default) result(value)
      class (type_property), intent(in)  :: self
      logical, optional,     intent(out) :: success
      integer, optional,     intent(in)  :: default
      integer                            :: value

      if (present(success)) success = .true.
      select type (self)
      class is (type_integer_property)
         value = self%value
         return
      class is (type_string_property)
         read(self%value,*,err=99,end=99) value
         return
      end select
99    if (present(success)) success = .false.
      if (present(default)) value = default
   end function

   function to_logical(self, success, default) result(value)
      class (type_property), intent(in)  :: self
      logical, optional,     intent(out) :: success
      logical, optional,     intent(in)  :: default
      logical                            :: value

      if (present(success)) success = .true.
      select type (self)
      class is (type_logical_property)
         value = self%value
         return
      class is (type_string_property)
         read(self%value,*,err=99,end=99) value
         return
      end select
99    if (present(success)) success = .false.
      if (present(default)) value = default
   end function

   function to_string(self, success, default) result(value)
      class (type_property),      intent(in)  :: self
      logical,          optional, intent(out) :: success
      character(len=*), optional, intent(in)  :: default
      character(value_string_length)          :: value

      if (present(success)) success = .true.
      select type (self)
      class is (type_real_property)
         write (value,'(g13.6)') self%value
      class is (type_integer_property)
         write (value,'(i0)') self%value
      class is (type_logical_property)
         if (self%value) then
            value = '.true.'
         else
            value = '.false.'
         end if
      class is (type_string_property)
         value = self%value
      class default
         if (present(success)) success = .false.
         if (present(default)) value = default
      end select
   end function

   function string_lower(string) result (lowerstring)
       character(len=*), intent(in) :: string
       character(len=len(string))   :: lowerstring

       integer :: i, k

       lowerstring = string
       do i = 1, len(string)
           k = iachar(string(i:i))
           if (k >= iachar('A') .and. k <= iachar('Z')) then
               k = k + iachar('a') - iachar('A')
               lowerstring(i:i) = achar(k)
           end if
       end do
   end function string_lower

   function compare_keys(self, key1, key2) result(equal)
      class (type_property_dictionary), intent(in) :: self
      character(len=*),                 intent(in) :: key1, key2
      logical                                      :: equal
      equal = string_lower(key1) == string_lower(key2)
   end function

   subroutine set_property(self, property, overwrite)
      class (type_property_dictionary), intent(inout) :: self
      class (type_property),            intent(in)    :: property
      logical,optional,                 intent(in)    :: overwrite

      class (type_property), pointer        :: current, previous
      logical                               :: overwrite_eff
      character(len=metadata_string_length) :: key

      overwrite_eff = .true.
      if (present(overwrite)) overwrite_eff = overwrite

      key = string_lower(property%name)

      ! First determine if a property with this name already exists (if so, delete it)
      previous => null()
      current => self%first
      do while (associated(current))
         if (current%key==key) then
            ! We found a property with the specified name - if we are not allowed to oevrwrite it, we're done.
            if (.not. overwrite_eff) return

            ! We are allowed to overwrite the existing property. Remove it from the list and deallocate it.
            if (associated(previous)) then
               ! Second or further down the list.
               previous%next => current%next
            else
               ! First in the list.
               self%first => current%next
            end if
            deallocate(current)
            exit
         end if
         previous => current
         current => previous%next
      end do

      if (.not. associated(self%first)) then
         ! First property in list
         allocate(self%first, source=property)
         current => self%first
      else
         ! Look for last element in list.
         current => self%first
         do while (associated(current%next))
            current => current%next
         end do
         allocate(current%next, source=property)
         current => current%next
      end if
      current%key = key
      current%next => null()
   end subroutine

   subroutine update(self, source, overwrite)
      class (type_property_dictionary), intent(inout) :: self
      class (type_property_dictionary), intent(in)    :: source
      logical,optional,                 intent(in)    :: overwrite
      class (type_property), pointer                  :: current

      current => source%first
      do while (associated(current))
         call self%set_property(current, overwrite)
         current => current%next
      end do
   end subroutine

   subroutine set_real(self, name, value)
      class (type_property_dictionary), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      real(rk),                         intent(in)    :: value
      call self%set_property(type_real_property(name=name, value=value))
   end subroutine

   subroutine set_integer(self, name, value)
      class (type_property_dictionary), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      integer,                          intent(in)    :: value
      call self%set_property(type_integer_property(name=name, value=value))
   end subroutine

   subroutine set_logical(self, name, value)
      class (type_property_dictionary), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      logical,                          intent(in)    :: value
      call self%set_property(type_logical_property(name=name, value=value))
   end subroutine

   subroutine set_string(self, name, value)
      class (type_property_dictionary), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      character(len=*),                 intent(in)    :: value
      call self%set_property(type_string_property(name=name, value=value))
   end subroutine

   function get_property(self, name) result(property)
      class (type_property_dictionary), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      class (type_property), pointer                  :: property

      character(len=len(name)) :: key

      key = string_lower(name)
      property => self%first
      do while (associated(property))
         if (property%key == key) return
         property => property%next
      end do
   end function

   function get_logical(self, name, default) result(value)
      class (type_property_dictionary), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      logical,                          intent(in)    :: default
      logical                                         :: value

      class (type_property), pointer :: property

      value = default
      property => self%get_property(name)
      if (.not. associated(property)) return
      value = property%to_logical(default=default)
   end function

   function get_integer(self, name, default) result(value)
      class (type_property_dictionary), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      integer,                          intent(in)    :: default
      integer                                         :: value

      class (type_property), pointer :: property

      value = default
      property => self%get_property(name)
      if (.not. associated(property)) return
      value = property%to_integer(default=default)
   end function

   function get_real(self, name, default) result(value)
      class (type_property_dictionary), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      real(rk),                         intent(in)    :: default
      real(rk)                                        :: value

      class (type_property), pointer :: property

      value = default
      property => self%get_property(name)
      if (.not. associated(property)) return
      value = property%to_real(default=default)
   end function

   function get_string(self, name, default) result(value)
      class (type_property_dictionary),  intent(inout) :: self
      character(len=*),                  intent(in)    :: name
      character(len=*),                  intent(in)    :: default
      character(len=value_string_length)               :: value

      class (type_property), pointer :: property

      value = default
      property => self%get_property(name)
      if (.not. associated(property)) return
      value = property%to_string(default=default)
   end function

   subroutine finalize(self)
      class (type_property_dictionary), intent(inout) :: self

      class (type_property), pointer :: property, next

      property => self%first
      do while (associated(property))
         next => property%next
         deallocate(property)
         property => next
      end do
      self%first => null()
   end subroutine finalize

   subroutine set_add(self, string)
      class (type_set), intent(inout) :: self
      character(len=*), intent(in)    :: string

      type (type_set_element), pointer :: element, previous

      if (.not. associated(self%first)) then
         allocate(self%first)
         element => self%first
      else
         element => self%first
         do while (associated(element))
            if (element%string == string) return
            previous => element
            element => element%next
         end do
         allocate(previous%next)
         element => previous%next
      end if
      element%string = string
   end subroutine

   function set_size(self) result(n)
      class (type_set), intent(in) :: self

      integer                          :: n
      type (type_set_element), pointer :: element

      n = 0
      element => self%first
      do while (associated(element))
         n = n + 1
         element => element%next
      end do
   end function

   subroutine set_to_array(self, array)
      class (type_set),             intent(in)  :: self
      character(len=*),allocatable, intent(out) :: array(:)

      integer                         :: n
      type (type_set_element),pointer :: element

      allocate(array(self%size()))
      n = 0
      element => self%first
      do while (associated(element))
         n = n + 1
         array(n) = element%string
         element => element%next
      end do
   end subroutine

   subroutine set_finalize(self)
      class (type_set), intent(inout) :: self

      type (type_set_element), pointer :: element, next

      element => self%first
      do while (associated(element))
         next => element%next
         deallocate(element)
         element => next
      end do
      self%first => null()
   end subroutine

end module fabm_properties

#include "fabm_driver.h"

! This module contains types and procedures for generic properties
! that internally can be of different types (real, integer, logical)
! It also offers a type that represents a property dictionary.

module fabm_properties

   implicit none

   private

   public type_property,type_property_dictionary,type_set
   public type_integer_property,type_real_property,type_logical_property,type_string_property

   integer, parameter :: rk = _FABM_REAL_KIND_
   integer, parameter :: metadata_string_length = 256
   integer, parameter :: value_string_length = 1024

   integer, parameter :: typecode_unknown = -1, typecode_real = 1, typecode_integer = 2, typecode_logical = 3, typecode_string = 4

   type,abstract :: type_property
      character(len=metadata_string_length) :: name       = ''
      character(len=metadata_string_length) :: long_name  = ''
      character(len=metadata_string_length) :: units      = ''
      logical                               :: accessed   = .false.
      class (type_property), pointer        :: next       => null()
   contains
      procedure :: typecode
      procedure :: to_real
      procedure :: to_integer
      procedure :: to_logical
      procedure :: to_string
   end type

   type,extends(type_property) :: type_integer_property
      integer :: value
   end type

   type,extends(type_property) :: type_real_property
      real(rk) :: value
   end type

   type,extends(type_property) :: type_logical_property
      logical :: value
   end type

   type,extends(type_property) :: type_string_property
      character(len=value_string_length) :: value
   end type

   type type_property_dictionary
      class (type_property),pointer :: first => null()

      contains

      procedure :: set_property
      procedure :: set_real
      procedure :: set_integer
      procedure :: set_logical
      procedure :: set_string

      procedure :: get_property
      procedure :: get_real
      procedure :: get_integer
      procedure :: get_logical
      procedure :: get_string

      procedure :: delete => delete_property

      procedure :: update

      procedure :: size => get_size
      procedure :: keys

      procedure :: compare_keys

      procedure :: reset_accessed
   end type

   type type_set_element
      character(len=metadata_string_length) :: string
      type (type_set_element), pointer :: next => null()
   end type

   type type_set
      type (type_set_element), pointer :: first => null()
   contains
      procedure :: contains => set_contains
      procedure :: add      => set_add
      procedure :: size     => set_size
      procedure :: to_array => set_to_array
      procedure :: finalize => set_finalize
   end type

contains

   integer function typecode(property)
      class (type_property),intent(in) :: property
      select type (property)
         class is (type_real_property)
            typecode = typecode_real
         class is (type_integer_property)
            typecode = typecode_integer
         class is (type_logical_property)
            typecode = typecode_logical
         class is (type_string_property)
            typecode = typecode_string
         class default
            typecode = typecode_unknown
      end select   
   end function

   function to_real(property,success,default) result(value)
      class (type_property),intent(in)  :: property
      logical, optional,    intent(out) :: success
      real(rk), optional,   intent(in)  :: default
      real(rk)                          :: value

      if (present(success)) success = .true.
      select type (property)
         class is (type_real_property)
            value = property%value
            return
         class is (type_integer_property)
            value = property%value
            return
         class is (type_string_property)
            read(property%value,*,err=99,end=99) value
            return
      end select
99    if (present(success)) success = .false.
      if (present(default)) value = default
   end function

   function to_integer(property,success,default) result(value)
      class (type_property),intent(in)  :: property
      logical, optional,    intent(out) :: success
      integer, optional,    intent(in)  :: default
      integer                           :: value

      if (present(success)) success = .true.
      select type (property)
         class is (type_integer_property)
            value = property%value
            return
         class is (type_string_property)
            read(property%value,*,err=99,end=99) value
            return
      end select
99    if (present(success)) success = .false.
      if (present(default)) value = default
   end function

   function to_logical(property,success,default) result(value)
      class (type_property),intent(in)  :: property
      logical, optional,    intent(out) :: success
      logical, optional,    intent(in)  :: default
      logical                           :: value

      if (present(success)) success = .true.
      select type (property)
         class is (type_logical_property)
            value = property%value
            return
         class is (type_string_property)
            read(property%value,*,err=99,end=99) value
            return
      end select
99    if (present(success)) success = .false.
      if (present(default)) value = default
   end function

   function to_string(property,success,default) result(value)
      class (type_property),     intent(in)  :: property
      logical, optional,         intent(out) :: success
      character(len=*), optional,intent(in)  :: default
      character(value_string_length)         :: value

      if (present(success)) success = .true.
      select type (property)
         class is (type_real_property)
            write (value,'(g13.6)') property%value
         class is (type_integer_property)
            write (value,'(i0)') property%value
         class is (type_logical_property)
            if (property%value) then
               value = '.true.'
            else
               value = '.false.'
            end if
         class is (type_string_property)
            value = property%value
         class default
            if (present(success)) success = .false.
            if (present(default)) value = default
      end select
   end function

   function string_lower(string) result (lowerstring) 
       character(len=*),intent(in) :: string 
       character(len=len(string))  :: lowerstring 

       integer                     :: i,k

       lowerstring = string 
       do i = 1,len(string) 
           k = iachar(string(i:i)) 
           if (k>=iachar('A').and.k<=iachar('Z')) then 
               k = k + iachar('a') - iachar('A') 
               lowerstring(i:i) = achar(k) 
           end if 
       end do 
   end function string_lower 

   function compare_keys(dictionary,key1,key2) result(equal)
      class (type_property_dictionary),intent(in) :: dictionary
      character(len=*),                intent(in) :: key1,key2
      logical                                     :: equal
      equal = string_lower(key1)==string_lower(key2)
   end function

   subroutine set_property(dictionary,property,overwrite)
      class (type_property_dictionary),intent(inout) :: dictionary
      class (type_property),           intent(in)    :: property
      logical,optional,                intent(in)    :: overwrite

      class (type_property),pointer                  :: current,next
      logical                                        :: overwrite_eff

      overwrite_eff = .true.
      if (present(overwrite)) overwrite_eff = overwrite

      nullify(next)
      if (.not.associated(dictionary%first)) then
         ! First property in list
         allocate(dictionary%first,source=property)
         current => dictionary%first
      else
         ! List already contains one or more properties.
         if (dictionary%compare_keys(dictionary%first%name,property%name)) then
            ! The provided property replaces the head of the list.
            if (.not.overwrite_eff) return
            next => dictionary%first%next
            deallocate(dictionary%first)
            allocate(dictionary%first,source=property)
            current => dictionary%first
         else
            ! Look for last element in list, or the one prior to the element with the same name.
            current => dictionary%first
            do while (associated(current%next))
               if (dictionary%compare_keys(current%next%name,property%name)) exit
               current => current%next
            end do
            if (associated(current%next)) then
               ! We are replacing current%next
               if (.not.overwrite_eff) return
               next => current%next%next
               deallocate(current%next)
            end if
            allocate(current%next,source=property)
            current => current%next
         end if
      end if
      current%accessed = .true.
      current%next => next
   end subroutine

   subroutine update(target,source,overwrite)
      class (type_property_dictionary),intent(inout) :: target
      class (type_property_dictionary),intent(in)    :: source
      logical,optional,                intent(in)    :: overwrite
      class (type_property),pointer                  :: current

      current => source%first
      do while (associated(current))
         call target%set_property(current,overwrite)
         current => current%next
      end do
   end subroutine

   subroutine set_real(dictionary,name,value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      real(rk),                        intent(in)    :: value
      call dictionary%set_property(type_real_property(name=name,value=value))
   end subroutine

   subroutine set_integer(dictionary,name,value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      integer,                         intent(in)    :: value
      call dictionary%set_property(type_integer_property(name=name,value=value))
   end subroutine

   subroutine set_logical(dictionary,name,value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      logical,                         intent(in)    :: value
      call dictionary%set_property(type_logical_property(name=name,value=value))
   end subroutine

   subroutine set_string(dictionary,name,value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      character(len=*),                intent(in)    :: value
      call dictionary%set_property(type_string_property(name=name,value=value))
   end subroutine

   function get_property(dictionary,name) result(property)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      class (type_property),pointer                  :: property

      property => dictionary%first
      do while (associated(property))
         if (dictionary%compare_keys(property%name,name)) then
            property%accessed = .true.
            return
         end if
         property => property%next
      end do
   end function

   function get_logical(dictionary,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      logical,                         intent(in)    :: default
      logical                                        :: value

      class (type_property),pointer :: property
      
      value = default
      property => dictionary%get_property(name)
      if (.not.associated(property)) return
      value = property%to_logical(default=default)
   end function

   function get_integer(dictionary,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      integer,                         intent(in)    :: default
      integer                                        :: value

      class (type_property),pointer :: property
      
      value = default
      property => dictionary%get_property(name)
      if (.not.associated(property)) return
      value = property%to_integer(default=default)
   end function

   function get_real(dictionary,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      real(rk),                        intent(in)    :: default
      real(rk)                                       :: value

      class (type_property),pointer :: property

      value = default
      property => dictionary%get_property(name)
      if (.not.associated(property)) return
      value = property%to_real(default=default)
   end function

   function get_string(dictionary,name,default) result(value)
      class (type_property_dictionary),  intent(inout) :: dictionary
      character(len=*),                  intent(in)    :: name
      character(len=*),                  intent(in)    :: default
      character(len=value_string_length)               :: value

      class (type_property),pointer :: property
      
      value = default
      property => dictionary%get_property(name)
      if (.not.associated(property)) return
      value = property%to_string(default=default)
   end function

   subroutine delete_property(dictionary,name)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      class (type_property),pointer                  :: property,previous

      ! First consume properties with this name at the start of the list.
      do while (dictionary%compare_keys(dictionary%first%name,name))
         property => dictionary%first
         dictionary%first => property%next
         deallocate(property)
      end do

      ! Now look internally for properties with this name.
      previous => dictionary%first
      property => property%next
      do while (associated(property))
         if (dictionary%compare_keys(property%name,name)) then
            previous%next => property%next
            deallocate(property)
         else
            previous => property
         end if
         property => previous%next
      end do
   end subroutine

   function get_size(dictionary) result(n)
      class (type_property_dictionary),intent(in) :: dictionary
      integer :: n
      class (type_property),pointer :: property

      n = 0
      property => dictionary%first
      do while (associated(property))
         n = n + 1
         property => property%next
      end do
   end function

   subroutine keys(dictionary,names)
      class (type_property_dictionary),intent(in) :: dictionary
      character(len=*),allocatable,intent(out) :: names(:)

      integer :: n
      class (type_property),pointer :: property

      allocate(names(dictionary%size()))
      n = 0
      property => dictionary%first
      do while (associated(property))
         n = n + 1
         names(n) = trim(property%name)
         property => property%next
      end do
   end subroutine

   subroutine reset_accessed(dictionary)
      class (type_property_dictionary),intent(inout) :: dictionary

      class (type_property),pointer :: property

      property => dictionary%first
      do while (associated(property))
         property%accessed = .false.
         property => property%next
      end do
   end subroutine

   logical function set_contains(self,string)
      class (type_set),intent(in) :: self
      character(len=*),intent(in) :: string

      type (type_set_element),pointer :: element

      element => self%first
      do while (associated(element))
         if (element%string==string) then
            set_contains = .true.
            return
         end if
         element => element%next
      end do
      set_contains = .false.
   end function

   subroutine set_add(self,string)
      class (type_set),intent(inout) :: self
      character(len=*),intent(in)    :: string

      type (type_set_element),pointer :: element,previous

      if (.not.associated(self%first)) then
         allocate(self%first)
         element => self%first
      else
         element => self%first
         do while (associated(element))
            if (element%string==string) return
            previous => element
            element => element%next
         end do
         allocate(previous%next)
         element => previous%next
      end if
      element%string = string
   end subroutine

   function set_size(self) result(n)
      class (type_set),intent(in) :: self

      integer                         :: n
      type (type_set_element),pointer :: element

      n = 0
      element => self%first
      do while (associated(element))
         n = n + 1
         element => element%next
      end do
   end function

   subroutine set_to_array(self,array)
      class (type_set),            intent(in)  :: self
      character(len=*),allocatable,intent(out) :: array(:)

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
      class (type_set),intent(inout) :: self

      type (type_set_element),pointer :: element,next

      element => self%first
      do while (associated(element))
         next => element%next
         deallocate(element)
         element => next
      end do
   end subroutine

end module fabm_properties

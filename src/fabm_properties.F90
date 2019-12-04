#include "fabm_driver.h"

! This module contains types and procedures for generic properties
! that internally can be of different types (real, integer, logical)
! It also offers a type that represents a property dictionary.

module fabm_properties

   use fabm_parameters, only: rk => rki

   implicit none

   private

   public type_property,type_property_dictionary,type_set,type_hierarchical_dictionary
   public type_integer_property,type_real_property,type_logical_property,type_string_property

   integer, parameter :: metadata_string_length = 256
   integer, parameter :: value_string_length = 1024

   integer, parameter :: typecode_unknown = -1, typecode_real = 1, typecode_integer = 2, typecode_logical = 3, typecode_string = 4

   type,abstract :: type_property
      character(len=metadata_string_length) :: name        = ''
      character(len=metadata_string_length) :: long_name   = ''
      character(len=metadata_string_length) :: units       = ''
      character(len=metadata_string_length) :: key         = ''      ! name wih normalized case
      logical                               :: has_default = .false.
      class (type_property), pointer        :: next        => null()
   contains
      procedure :: typecode
      procedure :: to_real
      procedure :: to_integer
      procedure :: to_logical
      procedure :: to_string
   end type

   type,extends(type_property) :: type_integer_property
      integer :: value
      integer :: default = 0
   end type

   type,extends(type_property) :: type_real_property
      real(rk) :: value
      real(rk) :: default = 0.0_rk
   end type

   type,extends(type_property) :: type_logical_property
      logical :: value
      logical :: default = .true.
   end type

   type,extends(type_property) :: type_string_property
      character(len=value_string_length) :: value
      character(len=value_string_length) :: default= ''
   end type

   type type_property_dictionary
      class (type_property),pointer :: first => null()

      contains

      procedure :: set_property
      procedure :: set_real
      procedure :: set_integer
      procedure :: set_logical
      procedure :: set_string
      generic :: set => set_property, set_real, set_integer, set_logical, set_string

      procedure :: get_property_by_name
      procedure :: get_property_by_index
      generic   :: get_property => get_property_by_name,get_property_by_index
      procedure :: get_real
      procedure :: get_integer
      procedure :: get_logical
      procedure :: get_string

      procedure :: delete_by_name
      procedure :: delete_by_index
      generic   :: delete => delete_by_name, delete_by_index

      procedure :: update

      procedure :: size => get_size
      procedure :: keys

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
      procedure :: contains => set_contains
      procedure :: add      => set_add
      procedure :: discard  => set_discard
      procedure :: size     => set_size
      procedure :: to_array => set_to_array
      procedure :: finalize => set_finalize
   end type

   type,extends(type_property_dictionary) :: type_hierarchical_dictionary
      type (type_set)                               :: retrieved
      type (type_set)                               :: missing
      character(len=metadata_string_length)         :: name = ''
      class (type_hierarchical_dictionary), pointer :: parent => null()
   contains
      procedure :: find_in_tree => hierarchical_dictionary_find_in_tree
      procedure :: set_in_tree  => hierarchical_dictionary_set_in_tree
      procedure :: add_child    => hierarchical_dictionary_add_child
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

      class (type_property),pointer                  :: current,previous
      logical                                        :: overwrite_eff
      character(len=metadata_string_length)          :: key

      overwrite_eff = .true.
      if (present(overwrite)) overwrite_eff = overwrite

      key = string_lower(property%name)

      ! First determine if a property with this name already exists (if so, delete it)
      previous => null()
      current => dictionary%first
      do while (associated(current))
         if (current%key==key) then
            ! We found a property with the specified name - if we are not allowed to oevrwrite it, we're done.
            if (.not.overwrite_eff) return

            ! We are allowed to overwrite the existing property. Remove it from the list and deallocate it.
            if (associated(previous)) then
               ! Second or further down the list.
               previous%next => current%next
            else
               ! First in the list.
               dictionary%first => current%next
            end if
            deallocate(current)
            exit
         end if
         previous => current
         current => previous%next
      end do

      if (.not.associated(dictionary%first)) then
         ! First property in list
         allocate(dictionary%first,source=property)
         current => dictionary%first
      else
         ! Look for last element in list.
         current => dictionary%first
         do while (associated(current%next))
            current => current%next
         end do
         allocate(current%next,source=property)
         current => current%next
      end if
      current%key = key
      current%next => null()
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

   function get_property_by_name(dictionary,name) result(property)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      class (type_property),pointer                  :: property

      character(len=len(name))                       :: key

      key = string_lower(name)
      property => dictionary%first
      do while (associated(property))
         if (property%key==key) return
         property => property%next
      end do
   end function

   function get_property_by_index(dictionary,index) result(property)
      class (type_property_dictionary),intent(inout) :: dictionary
      integer,                         intent(in)    :: index
      class (type_property),pointer                  :: property

      integer                                        :: i

      property => dictionary%first
      do i=2,index
         if (.not.associated(property)) return
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

   subroutine delete_by_name(dictionary,name)
      class (type_property_dictionary),intent(inout) :: dictionary
      character(len=*),                intent(in)    :: name
      class (type_property),pointer                  :: property,previous

      character(len=len(name))                       :: key

      key = string_lower(name)

      ! First consume properties with this name at the start of the list.
      do while (dictionary%first%key==key)
         property => dictionary%first
         dictionary%first => property%next
         deallocate(property)
      end do

      ! Now look internally for properties with this name.
      previous => dictionary%first
      property => previous%next
      do while (associated(property))
         if (property%key==key) then
            previous%next => property%next
            deallocate(property)
         else
            previous => property
         end if
         property => previous%next
      end do
   end subroutine

   subroutine delete_by_index(dictionary,index)
      class (type_property_dictionary),intent(inout) :: dictionary
      integer,                         intent(in)    :: index

      class (type_property),pointer                  :: property,previous
      integer                                        :: i

      if (.not.associated(dictionary%first)) return
      property => dictionary%first
      if (index==1) then
         ! Remove head
         dictionary%first => property%next
      else
         ! Remove non-head
         do i=2,index
            previous => property
            property => previous%next
            if (.not.associated(property)) return
         end do
         previous%next => property%next
      end if
      deallocate(property)
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

   subroutine finalize(dictionary)
      class (type_property_dictionary),intent(inout) :: dictionary

      class (type_property),pointer :: property, next

      property => dictionary%first
      do while (associated(property))
         next => property%next
         deallocate(property)
         property => next
      end do
      dictionary%first => null()
   end subroutine finalize

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

   subroutine set_discard(self,string)
      class (type_set),intent(inout) :: self
      character(len=*),intent(in)    :: string

      type (type_set_element),pointer :: previous,element

      previous => null()
      element => self%first
      do while (associated(element))
         if (element%string==string) exit
         previous => element
         element => element%next
      end do
      if (associated(element)) then
         if (associated(previous)) then
            previous%next => element%next
         else
            self%first => element%next
         end if
         deallocate(element)
      end if
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
      self%first => null()
   end subroutine

   function hierarchical_dictionary_find_in_tree(self,name) result(property)
      class (type_hierarchical_dictionary), intent(inout), target :: self
      character(len=*),                     intent(in)            :: name
      class (type_property), pointer :: property

      class (type_hierarchical_dictionary), pointer :: current_dictionary
      class (type_property),                pointer :: current_property
      character(len=metadata_string_length)         :: localname

      property => null()
      current_dictionary => self
      localname = name
      do while (associated(current_dictionary))
         ! Register that the value of this parameter was requested (i.e., used) by a biogeochemical model.
         call current_dictionary%retrieved%add(localname)

         current_property => current_dictionary%get_property(localname)
         if (associated(current_property)) property => current_property
         localname = trim(current_dictionary%name)//'/'//localname
         current_dictionary => current_dictionary%parent
      end do
      if (associated(property)) return

      ! Value not found. Register at all levels of the hierarchy that this parameter is missing.
      current_dictionary => self
      localname = name
      do while (associated(current_dictionary))
         call current_dictionary%missing%add(localname)
         localname = trim(current_dictionary%name)//'/'//localname
         current_dictionary => current_dictionary%parent
      end do
   end function hierarchical_dictionary_find_in_tree

   subroutine hierarchical_dictionary_set_in_tree(self,parameter)
   ! !INPUT PARAMETERS:
      class (type_hierarchical_dictionary), intent(inout), target :: self
      class (type_property),                intent(inout)         :: parameter
   !
   !EOP
      class (type_hierarchical_dictionary), pointer :: current_dictionary
      character(len=metadata_string_length)         :: oldname
   !-----------------------------------------------------------------------
   !BOC
      current_dictionary => self
      oldname = parameter%name
      do while (associated(current_dictionary))
         ! Store metadata
         call current_dictionary%set_property(parameter)
         parameter%name = trim(current_dictionary%name)//'/'//parameter%name
         current_dictionary => current_dictionary%parent
      end do
      parameter%name = oldname
   end subroutine hierarchical_dictionary_set_in_tree
   !EOC

   subroutine hierarchical_dictionary_add_child(self,child,name)
      class (type_hierarchical_dictionary), intent(in), target :: self
      class (type_hierarchical_dictionary), intent(inout)      :: child
      character(len=*),                     intent(in)         :: name

      child%parent => self
      child%name = name
   end subroutine hierarchical_dictionary_add_child

end module fabm_properties

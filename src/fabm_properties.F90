#include "fabm_driver.h"

! This module contains types and procedures for generic properties
! that internally can be of different types (real, integer, logical)
! It also offers a type that represents a property dictionary.

module fabm_properties

   implicit none

   private

   public type_property,type_property_dictionary,type_set,type_hierarchical_dictionary,type_key_property_pair,type_set_element
   public type_integer_property,type_real_property,type_logical_property,type_string_property

   integer, parameter :: rk = _FABM_REAL_KIND_
   integer, parameter :: metadata_string_length = 256
   integer, parameter :: value_string_length = 1024

   integer, parameter :: typecode_unknown = -1, typecode_real = 1, typecode_integer = 2, typecode_logical = 3, typecode_string = 4

   type,abstract :: type_property
      character(len=metadata_string_length) :: long_name   = ''
      character(len=metadata_string_length) :: units       = ''
      logical                               :: has_default = .false.
   contains
      procedure :: typecode
      procedure :: to_real
      procedure :: to_integer
      procedure :: to_logical
      procedure :: to_string
   end type

   type type_key_property_pair
      character(len=metadata_string_length)  :: key       = ''
      character(len=metadata_string_length)  :: name      = ''
      logical                                :: retrieved = .false.
      class (type_property),         pointer :: property  => null()
      type (type_key_property_pair), pointer :: next      => null()
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
      type (type_key_property_pair),pointer :: first => null()
      type (type_key_property_pair),pointer :: last  => null()

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
      procedure :: clear_retrieved

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

   type type_hierarchical_dictionary_node
      character(len=metadata_string_length)             :: name      = ''
      logical                                           :: retrieved = .false.
      type (type_hierarchical_dictionary),      pointer :: p         => null()
      type (type_hierarchical_dictionary_node), pointer :: next      => null()
   end type

   type,extends(type_property_dictionary) :: type_hierarchical_dictionary
      type (type_set)                                    :: missing
      class (type_hierarchical_dictionary_node), pointer :: first_child => null()
   contains
      procedure :: add_child            => hierarchical_dictionary_add_child
      procedure :: set_property         => hierarchical_dictionary_set_property
      procedure :: collect_unretrieved  => hierarchical_dictionary_collect_unretrieved
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

   function compare_keys(self,key1,key2) result(equal)
      class (type_property_dictionary),intent(in) :: self
      character(len=*),                intent(in) :: key1,key2
      logical                                     :: equal
      equal = string_lower(key1)==string_lower(key2)
   end function

   subroutine check(self)
      class (type_property_dictionary),intent(in) :: self
      type (type_key_property_pair),pointer :: current,previous
#ifndef NDEBUG
      previous => null()
      current => self%first
      do while (associated(current))
         if (.not.associated(current%property)) &
            stop 'property not associated'
         previous => current
         current => previous%next
      end do
      if (.not.associated(self%last,previous)) &
         stop 'last does not macth actual last node'
#endif
   end subroutine

   subroutine set_property(self,name,property,overwrite)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      class (type_property),           intent(in)    :: property
      logical,optional,                intent(in)    :: overwrite

      type (type_key_property_pair),pointer :: current,previous
      logical                               :: overwrite_eff
      character(len=len(name))              :: key

      overwrite_eff = .true.
      if (present(overwrite)) overwrite_eff = overwrite

      key = string_lower(name)

      ! First determine if a property with this name already exists
      ! If so, move it to the end of the list as we preserve the order is which keys are inserted.
      previous => null()
      current => self%first
      do while (associated(current))
         if (current%key==key) then
            ! We found a property with the specified name - if we are not allowed to overwrite it, we're done.
            if (.not.overwrite_eff) return

            ! We are allowed to overwrite the existing property. Move key/value pair to end of the list
            if (associated(self%last,current)) then
               ! Already last in the list
               exit
            elseif (associated(previous)) then
               ! Second or further down the list (but not last)
               previous%next => current%next
            elseif (associated(current%next)) then
               ! First in the list.
               self%first => current%next
            end if
            self%last%next => current
            self%last => current
            current%next => null()
            deallocate(current%property)
            exit
         end if
         previous => current
         current => current%next
      end do

      if (.not.associated(current)) then
         ! Key not found - create a new key-property pair and append to end of list
         if (.not.associated(self%first)) then
            ! First property in list
            allocate(self%first)
            self%last => self%first
         else
            ! Look for last element in list.
            allocate(self%last%next)
            self%last => self%last%next
         end if
         current => self%last
         current%key = key
      end if
      current%name = name
      current%retrieved = .true.
      allocate(current%property,source=property)
      call check(self)
   end subroutine

   subroutine clear_retrieved(self)
      class (type_property_dictionary),intent(inout) :: self

      type (type_key_property_pair),pointer :: current

      current => self%first
      do while (associated(current))
         current%retrieved = .false.
         current => current%next
      end do
   end subroutine

   subroutine update(self,source,overwrite)
      class (type_property_dictionary),intent(inout) :: self
      class (type_property_dictionary),intent(in)    :: source
      logical,optional,                intent(in)    :: overwrite

      type (type_key_property_pair),pointer :: current

      current => source%first
      do while (associated(current))
         call self%set_property(current%name,current%property,overwrite)
         current => current%next
      end do
   end subroutine

   subroutine set_real(self,name,value)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      real(rk),                        intent(in)    :: value
      call self%set_property(name,type_real_property(value=value))
   end subroutine

   subroutine set_integer(self,name,value)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      integer,                         intent(in)    :: value
      call self%set_property(name,type_integer_property(value=value))
   end subroutine

   subroutine set_logical(self,name,value)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      logical,                         intent(in)    :: value
      call self%set_property(name,type_logical_property(value=value))
   end subroutine

   subroutine set_string(self,name,value)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      character(len=*),                intent(in)    :: value
      call self%set_property(name,type_string_property(value=value))
   end subroutine

   function get_property_by_name(self,name) result(property)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      class (type_property),pointer                  :: property

      character(len=len(name))              :: key
      type (type_key_property_pair),pointer :: current

      key = string_lower(name)
      current => self%first
      do while (associated(current))
         if (current%key==key) then
            property => current%property
            current%retrieved = .true.
            return
         end if
         current => current%next
      end do
      property => null()
   end function

   function get_property_by_index(self,index) result(property)
      class (type_property_dictionary),intent(inout) :: self
      integer,                         intent(in)    :: index
      class (type_property),pointer                  :: property

      integer                               :: i
      type (type_key_property_pair),pointer :: current

      i = 1
      current => self%first
      do while (associated(current))
         if (i==index) then
            property => current%property
            current%retrieved = .true.
            return
         end if
         current => current%next
         i = i + 1
      end do
      property => null()
   end function

   function get_logical(self,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      logical,                         intent(in)    :: default
      logical                                        :: value

      class (type_property),pointer :: property

      value = default
      property => self%get_property(name)
      if (.not.associated(property)) return
      value = property%to_logical(default=default)
   end function

   function get_integer(self,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      integer,                         intent(in)    :: default
      integer                                        :: value

      class (type_property),pointer :: property

      value = default
      property => self%get_property(name)
      if (.not.associated(property)) return
      value = property%to_integer(default=default)
   end function

   function get_real(self,name,default) result(value)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name
      real(rk),                        intent(in)    :: default
      real(rk)                                       :: value

      class (type_property),pointer :: property

      value = default
      property => self%get_property(name)
      if (.not.associated(property)) return
      value = property%to_real(default=default)
   end function

   function get_string(self,name,default) result(value)
      class (type_property_dictionary),  intent(inout) :: self
      character(len=*),                  intent(in)    :: name
      character(len=*),                  intent(in)    :: default
      character(len=value_string_length)               :: value

      class (type_property),pointer :: property

      value = default
      property => self%get_property(name)
      if (.not.associated(property)) return
      value = property%to_string(default=default)
   end function

   subroutine delete_by_name(self,name)
      class (type_property_dictionary),intent(inout) :: self
      character(len=*),                intent(in)    :: name

      type (type_key_property_pair),pointer :: current,previous
      character(len=len(name))              :: key

      key = string_lower(name)
      previous => null()
      current => self%first
      do while (associated(current))
         if (current%key==key) then
            if (associated(previous)) then
               previous%next => current%next
            else
               self%first => current%next
            end if
            deallocate(current%property)
            deallocate(current)
            if (associated(self%last,current)) self%last => previous
            call check(self)
            return
         end if
         previous => current
         current => current%next
      end do
   end subroutine

   subroutine delete_by_index(self,index)
      class (type_property_dictionary),intent(inout) :: self
      integer,                         intent(in)    :: index

      type (type_key_property_pair),pointer :: current,previous
      integer                            :: i

      previous => null()
      current => self%first
      do while (associated(current))
         if (i==index) then
            if (associated(previous)) then
               previous%next => current%next
            else
               self%first => current%next
            end if
            deallocate(current%property)
            deallocate(current)
            if (associated(self%last,current)) self%last => previous
            call check(self)
            return
         end if
         previous => current
         current => current%next
         i = i + 1
      end do
   end subroutine

   function get_size(self) result(n)
      class (type_property_dictionary),intent(in) :: self
      integer :: n
      type (type_key_property_pair),pointer :: current

      n = 0
      current => self%first
      do while (associated(current))
         n = n + 1
         current => current%next
      end do
   end function

   subroutine keys(self,names)
      class (type_property_dictionary),intent(in) :: self
      character(len=*),allocatable,intent(out) :: names(:)

      integer :: n
      type (type_key_property_pair),pointer :: current

      allocate(names(self%size()))
      n = 0
      current => self%first
      do while (associated(current))
         n = n + 1
         names(n) = trim(current%name)
         current => current%next
      end do
   end subroutine

   subroutine finalize(self)
      class (type_property_dictionary),intent(inout) :: self

      type (type_key_property_pair),pointer :: current, next

      current => self%first
      do while (associated(current))
         next => current%next
         deallocate(current%property)
         deallocate(current)
         current => next
      end do
      self%first => null()
      self%last  => null()
      call check(self)
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

      nullify(previous)
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
      nullify(self%first)
   end subroutine

   function hierarchical_dictionary_get_child(self,name,create) result(node)
      class (type_hierarchical_dictionary), intent(inout), target :: self
      character(len=*),                     intent(in)            :: name
      logical,                              intent(in)            :: create
      type (type_hierarchical_dictionary_node), pointer :: node

      node => self%first_child
      do while (associated(node))
         if (node%name==name) return
         node => node%next
      end do
      allocate(node)
      if (create) allocate(node%p)
      node%name = name
      node%next => self%first_child
      self%first_child => node
   end function hierarchical_dictionary_get_child

   subroutine hierarchical_dictionary_add_child(self,child,name)
      class (type_hierarchical_dictionary), intent(inout), target :: self
      type  (type_hierarchical_dictionary), intent(inout), target :: child
      character(len=*),                     intent(in)            :: name

      type (type_hierarchical_dictionary_node), pointer :: node

      node => hierarchical_dictionary_get_child(self,name,create=.false.)
      if (associated(node%p)) then
         if (node%retrieved) stop 'child%retrieved'
         call child%update(node%p)
         call node%p%finalize()
         deallocate(node%p)
      end if
      node%p => child
      node%retrieved = .true.
   end subroutine hierarchical_dictionary_add_child

   recursive subroutine hierarchical_dictionary_set_property(self,name,property,overwrite)
      class (type_hierarchical_dictionary),intent(inout) :: self
      character(len=*),                    intent(in)    :: name
      class (type_property),               intent(in)    :: property
      logical,optional,                    intent(in)    :: overwrite

      integer :: islash
      type (type_hierarchical_dictionary_node), pointer :: child

      islash = index(name,'/')
      if (islash/=0) then
         child => hierarchical_dictionary_get_child(self,name(:islash-1),create=.true.)
         call child%p%set_property(name(islash+1:),property,overwrite)
      else
         call self%type_property_dictionary%set_property(name,property,overwrite)
      end if
   end subroutine

   recursive subroutine hierarchical_dictionary_collect_unretrieved(self,set,prefix)
      class (type_hierarchical_dictionary),intent(in)    :: self
      type (type_set),                     intent(inout) :: set
      character(len=*),                    intent(in)    :: prefix

      type (type_key_property_pair),            pointer :: pair
      type (type_hierarchical_dictionary_node), pointer :: node

      pair => self%first
      do while (associated(pair))
         if (.not.pair%retrieved) call set%add(prefix//pair%name)
         pair => pair%next
      end do

      node => self%first_child
      do while (associated(node))
         call node%p%collect_unretrieved(set,prefix=prefix//trim(node%name)//'/')
         node => node%next
      end do
   end subroutine

end module fabm_properties

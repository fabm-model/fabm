module fabm_config_types

   implicit none

   private

   public type_node,type_dictionary,type_key_value_pair,type_scalar,type_null

   integer,parameter :: string_length = 1024

   type,abstract :: type_node
      character(len=string_length) :: path = ''
   contains
      procedure (node_dump),deferred :: dump
      procedure                      :: set_path => node_set_path
   end type

   abstract interface
      subroutine node_dump(self,unit,indent)
         import type_node
         class (type_node),intent(in) :: self
         integer,intent(in) :: unit,indent
      end subroutine
   end interface

   type,extends(type_node) :: type_scalar
      character(len=string_length) :: string = ''
   contains
      procedure :: dump => value_dump
      procedure :: to_logical => scalar_to_logical
   end type

   type,extends(type_node) :: type_null
   contains
      procedure :: dump => null_dump
   end type

   type type_key_value_pair
      character(len=string_length)       :: key   = ''
      class (type_node),         pointer :: value => null()
      logical                            :: accessed = .false.
      type (type_key_value_pair),pointer :: next  => null()
   end type

   type,extends(type_node) :: type_dictionary
      type (type_key_value_pair),pointer :: first => null()
   contains
      procedure :: get            => dictionary_get
      procedure :: get_string     => dictionary_get_string
      procedure :: set            => dictionary_set
      procedure :: set_string     => dictionary_set_string
      procedure :: dump           => dictionary_dump
      procedure :: flatten        => dictionary_flatten
      procedure :: reset_accessed => dictionary_reset_accessed
      procedure :: set_path       => dictionary_set_path
   end type

contains

   function dictionary_get_string(self,key,default) result(value)
      class (type_dictionary),intent(in) :: self
      character(len=*),       intent(in) :: key
      character(len=*),       intent(in) :: default
      character(len=string_length)       :: value
      class(type_node),pointer :: node
      value = default
      node => self%get(key)
      if (associated(node)) then
         select type (node)
            class is (type_scalar)
               value = node%string
         end select
      end if
   end function

   subroutine dictionary_reset_accessed(self)
      class (type_dictionary),intent(in) :: self
      type (type_key_value_pair),pointer :: pair
      pair => self%first
      do while (associated(pair))
         pair%accessed = .false.
         pair => pair%next
      end do
   end subroutine

   function dictionary_get(self,key) result(value)
      class (type_dictionary),intent(in) :: self
      character(len=*),       intent(in) :: key
      class(type_node),pointer           :: value

      type (type_key_value_pair),pointer :: pair

      nullify(value)
      pair => self%first
      do while (associated(pair))
         if (pair%key==key) exit
         pair => pair%next
      end do
      if (associated(pair)) then
         value => pair%value
         pair%accessed = .true.
      end if
   end function

   subroutine dictionary_set(self,key,value)
      class (type_dictionary),intent(inout) :: self
      character(len=*),       intent(in)    :: key
      class(type_node),target               :: value

      type (type_key_value_pair),pointer :: pair

      if (.not.associated(self%first)) then
         ! This will be the first pair.
         allocate(self%first)
         pair => self%first
      else
         ! Try to find a pair with the same key, or failing that, the last pair.
         pair => self%first
         do while (associated(pair%next))
            if (pair%key==key) exit
            pair => pair%next
         end do
         if (.not.pair%key==key) then
            ! Key did not exist yet, which must mean we are operating on the last existing pair.
            ! Append a new pair.
            allocate(pair%next)
            pair => pair%next
         end if
      end if

      ! Store key and value.
      pair%key = key
      pair%value => value
   end subroutine

   subroutine dictionary_set_string(self,key,value)
      class (type_dictionary),intent(inout) :: self
      character(len=*),       intent(in)    :: key,value
      type (type_scalar),pointer :: node
      allocate(node)
      node%string = value
      call self%set(key,node)
   end subroutine

   subroutine value_dump(self,unit,indent)
      class (type_scalar),intent(in) :: self
      integer,            intent(in) :: unit,indent
      write (unit,*) repeat(' ',indent)//trim(self%string)
   end subroutine

   subroutine null_dump(self,unit,indent)
      class (type_null),intent(in) :: self
      integer,          intent(in) :: unit,indent
      write (unit,*) repeat(' ',indent)//'null'
   end subroutine

   recursive subroutine dictionary_dump(self,unit,indent)
      class (type_dictionary),intent(in) :: self
      integer,                intent(in) :: unit,indent
      type (type_key_value_pair),pointer :: pair
      pair => self%first
      do while (associated(pair))
         select type (value=>pair%value)
            class is (type_scalar)
               write (unit,*) repeat(' ',indent)//trim(pair%key)//': '//trim(value%string)
            class is (type_dictionary)
               write (unit,*) repeat(' ',indent)//trim(pair%key)//':'
               call value%dump(unit,indent+2)
         end select
         pair => pair%next
      end do
   end subroutine

   recursive subroutine dictionary_flatten(self,target,prefix)
      class (type_dictionary),intent(in)    :: self
      type (type_dictionary), intent(inout) :: target
      character(len=*),       intent(in)    :: prefix

      type (type_key_value_pair),pointer :: pair

      pair => self%first
      do while (associated(pair))
         select type (value=>pair%value)
            class is (type_scalar)
               call target%set_string(prefix//trim(pair%key),value%string)
            class is (type_dictionary)
               call value%flatten(target,prefix=prefix//trim(pair%key)//'/')
         end select
         pair => pair%next
      end do
   end subroutine

   function scalar_to_logical(self,default,success) result(value)
      class (type_scalar),intent(in)  :: self
      logical,            intent(in)  :: default
      logical,optional,   intent(out) :: success
      logical                         :: value

      value = default
      if (present(success)) success = .true.

      read(self%string,*,err=99,end=99) value
      return
99    if (present(success)) success = .false.
   end function

   recursive subroutine node_set_path(self,path)
      class (type_node),intent(inout) :: self
      character(len=*), intent(in)    :: path
      self%path = path
   end subroutine

   recursive subroutine dictionary_set_path(self,path)
      class (type_dictionary),intent(inout) :: self
      character(len=*),       intent(in)    :: path
      
      type (type_key_value_pair),pointer :: pair

      self%path = path
      pair => self%first
      do while (associated(pair))
         call pair%value%set_path(trim(self%path)//'/'//trim(pair%key))
         pair => pair%next
      end do
   end subroutine

end module fabm_config_types
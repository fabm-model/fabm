#include "fabm_driver.h"
#include "fabm_private.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_job: derived types that describe a single job (module subroutien to call and their interdependencies)
!
! !INTERFACE:
module fabm_job

   use fabm_types
   use fabm_driver

   implicit none

   private

   public type_job, type_call, type_call_list
   public find_dependencies, find_variable_dependencies

   type type_copy_command
      integer :: read_index
      integer :: write_index
   end type

   type type_variable_set_node
      type (type_internal_variable),pointer :: target               => null()
      type (type_call_pointer),     pointer :: first_dependent_node => null()
      logical                               :: copy_to_prefetch     = .false.
      type (type_variable_set_node),pointer :: next                 => null()
   end type

   type type_variable_set
      type (type_variable_set_node), pointer :: first => null()
   contains
      procedure :: add      => variable_set_add
      procedure :: contains => variable_set_contains
      procedure :: finalize => variable_set_finalize
   end type

   type type_call_pointer
      type (type_call),       pointer :: p    => null()
      type(type_call_pointer),pointer :: next => null()
   end type

   type type_call
      class (type_base_model), pointer      :: model => null()
      integer                               :: source = source_unknown
      type (type_copy_command), allocatable :: copy_commands_int(:) ! interior variables to copy from write to read cache after call completes
      type (type_copy_command), allocatable :: copy_commands_hz(:)  ! horizontal variables to copy from write to read cache after call completes
      type (type_variable_set)              :: outputs              ! computed variables that a later called model requires
      type (type_call), pointer             :: next => null()
   contains
      procedure :: finalize   => call_list_node_finalize
   end type type_call

   type type_call_list
      type (type_call),     pointer :: first => null()
      type (type_call_list),pointer :: prior => null()
   contains
      procedure :: find       => call_list_find
      procedure :: filter     => call_list_filter
      procedure :: append_node => call_list_append_node
      procedure :: append_call => call_list_append_call
      generic :: append => append_node,append_call
      procedure :: extend     => call_list_extend
      procedure :: print      => call_list_print
      procedure :: initialize => call_list_initialize
      procedure :: finalize   => call_list_finalize
      procedure :: computes   => call_list_computes
   end type type_call_list

   type type_job
      integer :: domain = -1
      type (type_call_list) :: calls
      integer, allocatable :: prefill_type(:)
      real(rk),allocatable :: prefill_values(:)
      integer, allocatable :: prefill_index(:)
      integer, allocatable :: save_sources(:)
      integer, allocatable :: save_sources_hz(:)
   end type

contains

subroutine call_list_node_finalize(self)
   class (type_call), intent(inout) :: self

   if (allocated(self%copy_commands_int)) deallocate(self%copy_commands_int)
   if (allocated(self%copy_commands_hz))  deallocate(self%copy_commands_hz)
   call self%outputs%finalize()
end subroutine call_list_node_finalize

function variable_set_add(self,variable) result(node)
   class (type_variable_set),intent(inout) :: self
   type (type_internal_variable),target    :: variable

   type (type_variable_set_node), pointer :: node

   ! Check if this variable already exists.
   node => self%first
   do while (associated(node))
      if (associated(node%target,variable)) return
      node => node%next
   end do

   ! Create a new variable object and prepend it to the list.
   allocate(node)
   node%target => variable
   node%next => self%first
   self%first => node
end function variable_set_add

logical function variable_set_contains(self,variable)
   class (type_variable_set),intent(inout) :: self
   type (type_internal_variable),target    :: variable

   type (type_variable_set_node), pointer :: node

   variable_set_contains = .true.
   node => self%first
   do while (associated(node))
      if (associated(node%target,variable)) return
      node => node%next
   end do
   variable_set_contains = .false.
end function variable_set_contains

subroutine variable_set_finalize(self)
   class (type_variable_set), intent(inout) :: self

   type (type_variable_set_node),pointer :: node, next

   node => self%first
   do while (associated(node))
      next => node%next
      deallocate(node)
      node => next
   end do
end subroutine variable_set_finalize

subroutine call_list_initialize(call_list)
   class (type_call_list),intent(inout) :: call_list

   type (type_call),pointer :: node

   node => call_list%first
   do while (associated(node))
      call call_list_node_initialize(node)
      node => node%next
   end do
end subroutine call_list_initialize

recursive function call_list_find(self,model,source,not_stale) result(node)
   class (type_call_list), intent(in)        :: self
   class (type_base_model),intent(in),target :: model
   integer,                intent(in)        :: source
   logical,optional,       intent(in)        :: not_stale

   type (type_call),pointer :: node
   logical                            :: not_stale_

   node => self%first
   do while (associated(node))
      if (associated(node%model,model).and.node%source==source) return
      node => node%next
   end do
   not_stale_ = .false.
   if (present(not_stale)) not_stale_ = not_stale
   if (.not.not_stale_.and.associated(self%prior)) node => self%prior%find(model,source,not_stale)
end function call_list_find

subroutine call_list_extend(self,other)
   class (type_call_list), intent(inout) :: self
   class (type_call_list), intent(in)    :: other

   type (type_call),pointer :: last, node2

   ! Find tail of the list [if any]
   last => self%first
   if (associated(last)) then
      do while (associated(last%next))
         last => last%next
      end do
   end if

   node2 => other%first
   do while (associated(node2))
      if (associated(last)) then
         ! List contains one or more items - append to tail.
         allocate(last%next)
         last => last%next
      else
         ! List is empty - create first node.
         allocate(self%first)
         last => self%first
      end if
      last%model => node2%model
      last%source = node2%source
      node2 => node2%next
   end do
end subroutine call_list_extend

subroutine call_list_append_node(self,node)
   class (type_call_list), intent(inout) :: self
   type (type_call),pointer    :: node

   type (type_call),pointer :: previous_node

   if (associated(self%first)) then
      ! List contains one or more items - append to tail.
      previous_node => self%first
      do while (associated(previous_node%next))
         previous_node => previous_node%next
      end do
      previous_node%next => node
   else
      ! List is empty - new node becomes head.
      self%first => node
   end if
   node%next => null()
end subroutine call_list_append_node

function call_list_append_call(self,model,source) result(node)
   class (type_call_list), intent(inout)     :: self
   class (type_base_model),intent(in),target :: model
   integer,                intent(in)        :: source

   type (type_call),pointer :: node

   allocate(node)
   node%model => model
   node%source = source
   call self%append(node)
end function call_list_append_call

subroutine call_list_filter(self,source)
   class (type_call_list), intent(inout) :: self
   integer,                intent(in)    :: source

   type (type_call),pointer :: node, previous, next

   previous => null()
   node => self%first
   do while (associated(node))
      next => node%next
      if (node%source==source) then
         if (.not.associated(previous)) then
            self%first => next
         else
            previous%next => next
         end if
         call node%finalize()
         deallocate(node)
      else
         previous => node
      end if
      node => next
   end do
end subroutine call_list_filter

subroutine call_list_print(self)
   class (type_call_list), intent(in) :: self

   type (type_call),   pointer :: node
   type (type_variable_set_node),pointer :: variable_list_node
   integer                               :: n
   type (type_call_pointer),pointer :: pnode

   node => self%first
   do while (associated(node))
      write (*,'(a,": ",a)') trim(node%model%get_path()),trim(source2string(node%source))
      n = 0
      variable_list_node => node%outputs%first
      do while (associated(variable_list_node))
         n = n + 1
         write (*,'("   ",a)') trim(variable_list_node%target%name)
         pnode => variable_list_node%first_dependent_node
         do while (associated(pnode))
            write (*,'("     <- ",a,": ",a)') trim(pnode%p%model%get_path()),trim(source2string(pnode%p%source))
            pnode => pnode%next
         end do
         variable_list_node => variable_list_node%next
      end do
      node => node%next
   end do

   contains

      character(len=32) function source2string(source)
         integer, intent(in) :: source
         select case (source)
         case (0); source2string = 'unknown'
         case (1); source2string = 'do'
         case (2); source2string = 'do_column'
         case (3); source2string = 'do_bottom'
         case (4); source2string = 'do_surface'
         case (5); source2string = 'none'
         case (6); source2string = 'get_vertical_movement'
         case (7); source2string = 'do_horizontal'
         case default
            write (source2string,'(i0)') node%source
         end select
      end function source2string

end subroutine call_list_print

subroutine call_list_finalize(self)
   class (type_call_list), intent(inout) :: self

   type (type_call),pointer :: node,next

   node => self%first
   do while (associated(node))
      next => node%next
      call node%finalize()
      deallocate(node)
      node => next
   end do
   nullify(self%first)
end subroutine call_list_finalize

logical function call_list_computes(self,variable)
   class (type_call_list), intent(inout) :: self
   type (type_internal_variable),target  :: variable

   type (type_call),pointer :: node

   call_list_computes = .true.
   node => self%first
   do while (associated(node))
      if (node%outputs%contains(variable)) return
      node => node%next
   end do
   call_list_computes = .false.
end function call_list_computes

subroutine call_list_node_initialize(self)
   type (type_call),intent(inout) :: self

   class (type_base_model),      pointer :: parent
   class (type_model_list_node), pointer :: model_list_node
   type (type_call),             pointer :: later_call
   type (type_call_pointer),     pointer :: dependent_call
   type (type_variable_set_node),pointer :: variable_node, dummy_variable_node

   ! Make sure the pointer to the model has the highest class (and not a base class)
   ! This is needed because model classes that use inheritance and call base class methods
   ! may end up with model pointers that are base class specific (and do not reference
   ! procedures overwritten at a higher level)
   parent => self%model%parent
   if (associated(parent)) then
      model_list_node => parent%children%first
      do while (associated(model_list_node))
         if (associated(self%model,model_list_node%model)) then
            ! Found ourselves in our parent - use the parent pointer to replace ours.
            self%model => model_list_node%model
            exit
         end if
         model_list_node => model_list_node%next
      end do
   end if

   ! For all output variables that other model are interested in, decide whether to copy their value
   ! from the write to read cache [if the other model will be called as part of the same job],
   ! of to save it to the persistent data store.
   variable_node => self%outputs%first
   do while (associated(variable_node))
      dependent_call => variable_node%first_dependent_node
      do while (associated(dependent_call).and..not.variable_node%copy_to_prefetch)
         later_call => self%next
         do while (associated(later_call).and..not.variable_node%copy_to_prefetch)
            if (associated(later_call,dependent_call%p)) variable_node%copy_to_prefetch = .true.
            later_call => later_call%next
         end do
         dependent_call => dependent_call%next
      end do
      variable_node => variable_node%next
   end do

   call process_domain(self%copy_commands_int, domain_interior)
   call process_domain(self%copy_commands_hz,  domain_horizontal)

contains

   subroutine process_domain(commands,domain)
      type (type_copy_command), allocatable, intent(inout) :: commands(:)
      integer,                               intent(in)    :: domain

      type (type_variable_set_node), pointer :: node
      integer :: i, n, maxwrite

      n = 0
      maxwrite = -1
      node => self%outputs%first
      do while (associated(node))
         if (node%copy_to_prefetch.and.iand(node%target%domain,domain)/=0) then
            if (node%target%write_indices%is_empty()) then
               call driver%fatal_error('call_list_node_initialize','BUG: target without write indices')
            end if
            n = n + 1
            maxwrite = max(maxwrite,node%target%write_indices%value)
         end if
         node => node%next
      end do

      ! Create list of copy commands, sorted by write index
      allocate(commands(n))
      n = 0
      do i=1,maxwrite
         node => self%outputs%first
         do while (associated(node))
            if (node%copy_to_prefetch.and.iand(node%target%domain,domain)/=0.and.node%target%write_indices%value==i) then
               n = n + 1
               commands(n)%read_index = node%target%read_indices%value
               commands(n)%write_index = node%target%write_indices%value
               exit
            end if
            node => node%next
         end do
      end do
   end subroutine

end subroutine

recursive subroutine find_dependencies2(self,source,list,node,not_stale,forbidden)
   class (type_base_model),     intent(in),target   :: self
   integer,                     intent(in)          :: source
   type (type_call_list),       intent(inout)       :: list
   type (type_call),  pointer             :: node
   logical,                     intent(in)          :: not_stale
   type (type_call_list),       intent(in),optional :: forbidden

   type (type_link),pointer           :: link
   type (type_call_list)              :: forbidden_with_self
   character(len=2048)                :: chain
   logical                            :: same_source

   ! Check whether we are already processed this call.
   node => list%find(self,source,not_stale)
   if (associated(node)) return

   ! Check the list of forbidden model (i.e., models that indirectly request the current model)
   ! If the current model is on this list, there is a circular dependency between models.
   if (present(forbidden)) then
      node => forbidden%find(self,source)
      if (associated(node)) then
         ! Circular dependency found - report as fatal error.
         chain = ''
         do while (associated(node))
            chain = trim(chain)//' '//trim(node%model%get_path())//' ->'
            node => node%next
         end do
         call fatal_error('find_dependencies','circular dependency found: '//trim(chain(2:))//' '//trim(self%get_path()))
      end if
      call forbidden_with_self%extend(forbidden)
   end if
   node => forbidden_with_self%append(self,source)

   ! Create object to represent this call
   allocate(node)
   node%model => self
   node%source = source

   ! Loop over all variables, and if they belong to some other model, first add that model to the dependency list.
   link => self%links%first
   do while (associated(link))
      same_source = link%target%source==source .or. (link%target%source==source_unknown.and.(source==source_do_surface.or.source==source_do_bottom))
      if (index(link%name,'/')==0 &                                        ! Our own link...
            .and..not.link%target%write_indices%is_empty() &                 ! ...to a diagnostic variable...
            .and..not.(associated(link%target%owner,self).and.same_source) & ! ...not set by ourselves...
            .and.associated(link%original%read_index)) then                  ! ...and we do depend on its value

         call find_variable_dependencies(link%target,list,force_copy_to_prefetch=.false.,forbidden=forbidden_with_self,requester=node)
      end if
      link => link%next
   end do

   ! We're happy - add ourselves to the list of processed models.
   call list%append(node)

   ! Clean up our temporary list.
   call forbidden_with_self%finalize()

end subroutine find_dependencies2

recursive subroutine find_variable_dependencies(variable,list,force_copy_to_prefetch,forbidden,not_stale,requester)
   type (type_internal_variable),intent(in)          :: variable
   type (type_call_list),        intent(inout)       :: list
   logical,                      intent(in)          :: force_copy_to_prefetch
   type (type_call_list),        intent(in),optional :: forbidden
   logical,                      intent(in),optional :: not_stale
   type (type_call),target,optional        :: requester

   type (type_call),pointer :: node

   if (variable%source==source_unknown) then
      call find_dependencies2(variable%owner,source_do_surface,list,node,not_stale,forbidden)
      call register_dependencies()
      call find_dependencies2(variable%owner,source_do_bottom,list,node,not_stale,forbidden)
      call register_dependencies()
   elseif (variable%source/=source_none) then
      call find_dependencies2(variable%owner,variable%source,list,node,not_stale,forbidden)
      call register_dependencies()
   end if

   contains

   subroutine register_dependencies()
      type (type_call_pointer),     pointer :: requester_pointer
      type (type_variable_set_node),pointer :: variable_node

      variable_node => node%outputs%add(variable)
      if (force_copy_to_prefetch) variable_node%copy_to_prefetch = .true.
      if (present(requester)) then
         requester_pointer => variable_node%first_dependent_node
         do while (associated(requester_pointer))
            if (associated(requester_pointer%p,requester)) exit
            requester_pointer => requester_pointer%next
         end do
         if (.not.associated(requester_pointer)) then
            allocate(requester_pointer)
            requester_pointer%p => requester
            requester_pointer%next => variable_node%first_dependent_node
            variable_node%first_dependent_node => requester_pointer
         end if
      end if
   end subroutine register_dependencies
end subroutine find_variable_dependencies

recursive subroutine find_dependencies(self,list,forbidden)
   class (type_base_model),intent(in),target   :: self
   type (type_model_list), intent(inout)       :: list
   type (type_model_list), intent(in),optional :: forbidden

   type (type_link),pointer            :: link
   type (type_model_list)              :: forbidden_with_self
   type (type_model_list_node),pointer :: node
   character(len=2048)                 :: chain

   if (associated(list%find(self))) return

   ! Check the list of forbidden model (i.e., models that indirectly request the current model)
   ! If the current model is on this list, there is a circular dependency between models.
   if (present(forbidden)) then
      node => forbidden%find(self)
      if (associated(node)) then
         ! Circular dependency found - report as fatal error.
         chain = ''
         do while (associated(node))
            chain = trim(chain)//' '//trim(node%model%get_path())//' ->'
            node => node%next
         end do
         call fatal_error('find_dependencies','circular dependency found: '//trim(chain(2:))//' '//trim(self%get_path()))
      end if
      call forbidden_with_self%extend(forbidden)
   end if
   call forbidden_with_self%append(self)

   ! Loop over all variables, and if they belong to some other model, first add that model to the dependency list.
   link => self%links%first
   do while (associated(link))
      if (index(link%name,'/')==0 &                        ! Our own link...
          .and..not.link%target%write_indices%is_empty() & ! ...to a diagnostic variable...
          .and..not.associated(link%target%owner,self) &   ! ...not set by ourselves...
          .and.associated(link%original%read_index)) &     ! ...and we do depend on its value.
         call find_dependencies(link%target%owner,list,forbidden_with_self)
      link => link%next
   end do

   ! We're happy - add ourselves to the list of processed models.
   call list%append(self)

   ! Clean up our temporary list.
   call forbidden_with_self%finalize()
end subroutine find_dependencies
   
end module fabm_job

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#include "fabm_driver.h"
#include "fabm_private.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_job: derived types that describe a single job (module subroutines to call and their dependencies)
!
! !INTERFACE:
module fabm_graph

   ! Objects to describe the directed acyclic graph (DAG) that characterizes the dependencies between models.
   ! Each graph node represents a call, that is, a specific combination of a model and one of its APIs ("source").
   ! Graph nodes maintain both pointers to their dependencies (a set of other nodes, maintained at the node level),
   ! and to the nodes that depend on them (maintained separately per output variable).

   use fabm_types
   use fabm_driver

   implicit none

   public type_graph, type_node, type_node_set_member, type_node_list, type_node_list_member
   public type_output_variable, type_input_variable

   private

   type type_node_set_member
      type (type_node),            pointer :: p    => null()
      type (type_node_set_member), pointer :: next => null()
   end type

   type type_node_set
      type (type_node_set_member), pointer :: first => null()
   contains
      procedure :: add      => node_set_add
      procedure :: contains => node_set_contains
      procedure :: finalize => node_set_finalize
   end type

   type type_node_list_member
      type (type_node),             pointer :: p        => null()
      type (type_node_list_member), pointer :: next     => null()
      type (type_node_list_member), pointer :: previous => null()
   end type

   type type_node_list
      type (type_node_list_member), pointer :: first => null()
      type (type_node_list_member), pointer :: last  => null()
   contains
      procedure :: append        => node_list_append
      procedure :: insert_before => node_list_insert_before
      procedure :: pop           => node_list_pop
      procedure :: finalize      => node_list_finalize
      procedure :: as_string     => node_list_as_string
      procedure :: copy          => node_list_copy
      procedure :: find          => node_list_find
      procedure :: find_node     => node_list_find_node
      procedure :: check         => node_list_check
   end type

   type type_output_variable
      type (type_internal_variable), pointer :: target          => null()
      type (type_node_set)                   :: dependent_nodes
      logical                                :: copy_to_cache   = .false.
      logical                                :: copy_to_store   = .false.
      type (type_output_variable),   pointer :: next            => null()
   end type

   type type_output_variable_set
      type (type_output_variable), pointer :: first => null()
   contains
      procedure :: add      => output_variable_set_add
      procedure :: finalize => output_variable_set_finalize
   end type

   type type_input_variable
      type (type_internal_variable), pointer :: target => null()
      type (type_input_variable),    pointer :: next   => null()
   end type

   type type_input_variable_set
      type (type_input_variable), pointer :: first => null()
   contains
      procedure :: add      => input_variable_set_add
      procedure :: finalize => input_variable_set_finalize
   end type

   type type_node
      class (type_base_model), pointer :: model => null()
      integer                          :: source = source_unknown
      type (type_input_variable_set)   :: inputs               ! input variables (irrespective of their source - can be constants, state variables, host-provided, or diagnostics computed by another node)
      type (type_node_set)             :: dependencies         ! direct dependencies
      type (type_node_set)             :: all_dependencies     ! direct and indirect dependencies
      type (type_output_variable_set)  :: outputs              ! output variables that a later called model requires
   contains
      procedure :: as_string               => node_as_string
      procedure :: get_maximum_path_length => node_get_maximum_path_length
   end type

   type,extends(type_node_list) :: type_graph
      type (type_node_set)      :: endpoints
      type (type_graph),pointer :: previous => null()
   contains
      procedure :: add_call     => graph_add_call
      procedure :: add_variable => graph_add_variable
      procedure :: print        => graph_print
      procedure :: save_as_dot  => graph_save_as_dot
      procedure :: get_longest_path => graph_get_longest_path
   end type

contains

function output_variable_set_add(self,variable) result(node)
   class (type_output_variable_set),intent(inout) :: self
   type (type_internal_variable),target           :: variable

   type (type_output_variable), pointer :: node

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
end function output_variable_set_add

subroutine output_variable_set_finalize(self)
   class (type_output_variable_set), intent(inout) :: self

   type (type_output_variable),pointer :: node, next

   node => self%first
   do while (associated(node))
      next => node%next
      call node%dependent_nodes%finalize()
      deallocate(node)
      node => next
   end do
   self%first => null()
end subroutine output_variable_set_finalize

subroutine input_variable_set_add(self,variable)
   class (type_input_variable_set),intent(inout) :: self
   type (type_internal_variable),target          :: variable

   type (type_input_variable), pointer :: node

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
end subroutine input_variable_set_add

subroutine input_variable_set_finalize(self)
   class (type_input_variable_set), intent(inout) :: self

   type (type_input_variable),pointer :: node, next

   node => self%first
   do while (associated(node))
      next => node%next
      deallocate(node)
      node => next
   end do
   self%first => null()
end subroutine input_variable_set_finalize

subroutine graph_print(self)
   class (type_graph), intent(in) :: self

   type (type_node_list_member),pointer :: node
   type (type_output_variable), pointer :: variable
   type (type_node_set_member), pointer :: pnode

   node => self%first
   do while (associated(node))
      write (*,'(a,": ",a)') trim(node%p%model%get_path()),trim(source2string(node%p%source))
      variable => node%p%outputs%first
      do while (associated(variable))
         write (*,'("   ",a,",write@",i0)',advance='no') trim(variable%target%name),variable%target%write_indices%value
         if (variable%copy_to_cache) write (*,'(",cache@",i0)',advance='no') variable%target%read_indices%value
         if (variable%copy_to_store) write (*,'(",store@",i0)',advance='no') variable%target%store_index
         if (variable%target%prefill==prefill_missing_value) write (*,'(",prefill=",g0.6)',advance='no') variable%target%missing_value
         if (variable%target%prefill==prefill_previous_value) write (*,'(",prefill=previous")',advance='no')
         write (*,*)
         pnode => variable%dependent_nodes%first
         do while (associated(pnode))
            write (*,'("     <- ",a,": ",a)') trim(pnode%p%model%get_path()),trim(source2string(pnode%p%source))
            pnode => pnode%next
         end do
         variable => variable%next
      end do
      node => node%next
   end do

end subroutine graph_print

recursive function graph_add_call(self,model,source,outer_calls,not_stale) result(node)
   class (type_graph),            intent(inout) :: self
   class (type_base_model),target,intent(in)    :: model
   integer,                       intent(in)    :: source
   type (type_node_list),  target,intent(inout) :: outer_calls
   logical,optional,              intent(in)    :: not_stale
   type (type_node), pointer :: node

   logical                               :: not_stale_
   type (type_graph),            pointer :: previous_graph
   type (type_node_list_member), pointer :: pnode
   character(len=2048)                   :: chain
   type (type_link),             pointer :: link
   logical                               :: same_source

   ! First see if this node is already in the graph. If so, we are done: return.
   node => self%find(model,source)
   if (associated(node)) return

   ! Now check any preceding graphs (provided not_stale is not active).
   ! If the call is found in a preceding graph, we are done: return.
   not_stale_ = .false.
   if (present(not_stale)) not_stale_ = not_stale
   if (.not.not_stale_) then
      previous_graph => self%previous
      do while (associated(previous_graph))
         node => previous_graph%find(model,source)
         if (associated(node)) return
         previous_graph => previous_graph%previous
      end do
   end if

   ! Add current call to the list of outer calls.
   allocate(node)
   node%model => model
   node%source = source

   ! Check the list of outer model calls, i.e., calls that [indirectly] request the current call.
   ! If the current call is already on this list, it is indirectly calling itself: there is a circular dependency.
   pnode => outer_calls%find_node(model,source)
   if (associated(pnode)) then
      ! Circular dependency found - report as fatal error.
      chain = ''
      do while (associated(pnode))
         chain = trim(chain)//' '//trim(pnode%p%as_string())//' ->'
         pnode => pnode%next
      end do
      call driver%fatal_error('graph::add_call','circular dependency found: '//trim(chain(2:))//' '//trim(node%as_string()))
   end if

   ! First add this call to the list of requesting calls [a list of all calls higher on the call stack]
   ! This forbids any indirect dependency on this call, as such would be a circular dependency.
   call outer_calls%append(node)

   link => model%links%first
   do while (associated(link))
      if (index(link%name,'/')==0.and.associated(link%original%read_index)) then
         ! This is the model's own variable (not inherited from child model) and the model itself originally requested read access to it.
         call node%inputs%add(link%target)
         same_source = link%target%source==source .or. (link%target%source==source_unknown.and.(source==source_do_surface.or.source==source_do_bottom))
         if (.not.(associated(link%target%owner,model).and.same_source)) call self%add_variable(link%target,outer_calls,caller=node)
      end if
      link => link%next
   end do

   ! Remove node from the list of outer calls and add it to the graph instead.
   node => outer_calls%pop()
   call self%append(node)
end function graph_add_call

recursive subroutine graph_add_variable(self,variable,outer_calls,copy_to_cache,copy_to_store,not_stale,caller)
   class (type_graph),                intent(inout) :: self
   type (type_internal_variable),     intent(in)    :: variable
   type (type_node_list),target,intent(inout) :: outer_calls
   logical,optional,                  intent(in)    :: copy_to_cache
   logical,optional,                  intent(in)    :: copy_to_store
   logical,optional,                  intent(in)    :: not_stale
   type (type_node),target,optional      :: caller

   ! If this variable is not an output of some model (e.g., a state variable or external dependency), no call is needed.
   if (variable%write_indices%is_empty()) return

   if (variable%source==source_unknown) then
      ! This variable is either written by do_surface or do_bottom - which one of these two APIs is unknown.
      call add_call(source_do_surface)
      call add_call(source_do_bottom)
   elseif (variable%source/=source_none) then
      ! This variable is written by a known BGC API [is is not a constant indicated by source_none]
      call add_call(variable%source)
   end if

contains
   
   recursive subroutine add_call(source)
      integer, intent(in) :: source

      type (type_node),        pointer :: node
      type (type_output_variable),   pointer :: variable_node
      type (type_node_set_member),pointer :: pnode

      node => self%add_call(variable%owner,source,outer_calls,not_stale)
      variable_node => node%outputs%add(variable)
      if (present(copy_to_cache)) variable_node%copy_to_cache = variable_node%copy_to_cache .or. copy_to_cache
      if (present(copy_to_store)) variable_node%copy_to_store = variable_node%copy_to_store .or. copy_to_store
      if (present(caller)) then
         call caller%dependencies%add(node)
         call caller%all_dependencies%add(node)
         pnode => node%all_dependencies%first
         do while (associated(pnode))
            call caller%all_dependencies%add(pnode%p)
            pnode => pnode%next
         end do
         call variable_node%dependent_nodes%add(caller)
      else
         call self%endpoints%add(node)
      end if
   end subroutine add_call

end subroutine graph_add_variable

subroutine node_set_add(self,node)
   class (type_node_set), intent(inout) :: self
   type (type_node),target              :: node

   type (type_node_set_member), pointer :: member

   ! First determine if the graph node is already part of the set. If so, we are done: return.
   member => self%first
   do while (associated(member))
      if (associated(member%p,node)) return
      member => member%next
   end do

   ! Graph node is not in set yet. Create pointer object and prepend it to the set.
   allocate(member)
   member%p => node
   member%next => self%first
   self%first => member
end subroutine node_set_add

function node_set_contains(self,node) result(found)
   class (type_node_set), intent(in) :: self
   type (type_node), target          :: node
   logical                           :: found

   type (type_node_set_member), pointer :: member

   ! First determine if the graph node is already part of the set. If so, we are done: return.
   found = .true.
   member => self%first
   do while (associated(member))
      if (associated(member%p,node)) return
      member => member%next
   end do
   found = .false.
end function node_set_contains

subroutine node_set_finalize(self)
   class (type_node_set), intent(inout) :: self

   type (type_node_set_member), pointer :: current,next

   current => self%first
   do while (associated(current))
      next => current%next
      deallocate(current)
      current => next
   end do
   self%first => null()
end subroutine node_set_finalize

subroutine node_list_append(self,node)
   class (type_node_list), intent(inout) :: self
   type (type_node),target               :: node

   type (type_node_list_member), pointer :: member

   allocate(member)
   member%p => node
   if (associated(self%last)) then
      member%previous => self%last
      self%last%next => member
   else
      self%first => member
   end if
   self%last => member
   !call self%check()
end subroutine node_list_append

subroutine node_list_insert_before(self,next,node)
   class (type_node_list), intent(inout) :: self
   type (type_node_list_member), target  :: next
   type (type_node),target               :: node

   type (type_node_list_member), pointer :: member

   allocate(member)
   member%p => node

   member%next => next
   member%previous => next%previous
   if (associated(next%previous)) then
      next%previous%next => member
   else
      self%first => member
   end if
   next%previous => member
   !call self%check()
end subroutine node_list_insert_before

subroutine node_list_check(self)
   class (type_node_list), intent(in) :: self
   type (type_node_list_member), pointer :: current,previous

   if (associated(self%last).and..not.associated(self%first)) call driver%fatal_error('node_list_check','BUG: list has tail but no head.')
   if (associated(self%first).and..not.associated(self%last)) call driver%fatal_error('node_list_check','BUG: list has head but no tail.')
   if (.not.associated(self%first)) return
   if (associated(self%first%previous)) call driver%fatal_error('node_list_check','BUG: head of list has pointer to previous node.')
   if (associated(self%last%next)) call driver%fatal_error('node_list_check','BUG: tail of list has pointer to next node.')
   previous => null()
   current => self%first
   do while (associated(current))
      if (associated(previous)) then
         if (.not.associated(current%previous,previous)) call driver%fatal_error('node_list_check','BUG: previous pointer does not match actual previous node.')
         if (.not.associated(previous%next,current)) call driver%fatal_error('node_list_check','BUG: next pointer does not match actual next node.')
      end if
      previous => current
      current => current%next
   end do
end subroutine node_list_check

subroutine node_list_copy(self,copy,member,member_copy)
   class (type_node_list), intent(in)  :: self
   class (type_node_list), intent(out) :: copy
   type (type_node_list_member), target, optional  :: member
   type (type_node_list_member), pointer, optional :: member_copy

   type (type_node_list_member), pointer :: current,current_copy,previous_copy

   if (.not.associated(self%first)) return

   if (present(member_copy)) member_copy => null()
   previous_copy => null()
   current => self%first
   do while (associated(current))
      allocate(current_copy)
      current_copy%p => current%p
      if (associated(previous_copy)) then
         previous_copy%next => current_copy
         current_copy%previous => previous_copy
      else
         copy%first => current_copy
      end if
      if (present(member).and.present(member_copy)) then
         if (associated(current,member)) member_copy => current_copy
      end if
      previous_copy => current_copy
      current => current%next
   end do
   copy%last => previous_copy
   if (present(member_copy)) then
      if (.not.associated(member_copy)) call driver%fatal_error('node_list_copy','BUG: provided member not found in original list.')
   end if
   !call copy%check()
end subroutine node_list_copy

function node_list_pop(self) result(node)
   class (type_node_list), intent(inout) :: self
   type (type_node),pointer              :: node

   type (type_node_list_member), pointer :: previous

   node => null()
   if (associated(self%last)) then
      node => self%last%p
      previous => self%last%previous
      if (associated(previous)) then
         ! More than one node in list.
         previous%next => null()
      else
         ! Pop-ed node is only one in list.
         self%first => null()
      end if
      deallocate(self%last)
      self%last => previous
   end if
   !call self%check()
end function node_list_pop

subroutine node_list_finalize(self)
   class (type_node_list), intent(inout) :: self

   type (type_node_list_member), pointer :: current,next

   current => self%first
   do while (associated(current))
      next => current%next
      deallocate(current)
      current => next
   end do
   self%first => null()
   self%last => null()
end subroutine node_list_finalize

function node_list_find_node(self,model,source) result(pnode)
   class (type_node_list),  intent(in) :: self
   class (type_base_model),target,intent(in) :: model
   integer,                       intent(in) :: source
   type (type_node_list_member), pointer :: pnode

   pnode => self%first
   do while (associated(pnode))
      if (associated(pnode%p%model,model).and.pnode%p%source==source) return
      pnode => pnode%next
   end do
end function node_list_find_node

function node_list_find(self,model,source) result(node)
   class (type_node_list),  intent(in) :: self
   class (type_base_model),target,intent(in) :: model
   integer,                       intent(in) :: source
   type (type_node), pointer :: node
   type (type_node_list_member), pointer :: pnode

   node => null()
   pnode => node_list_find_node(self,model,source)
   if (associated(pnode)) node => pnode%p
end function node_list_find

subroutine graph_finalize(self)
   class (type_graph),intent(inout) :: self

   type (type_node_list_member), pointer :: current

   call self%endpoints%finalize()
   current => self%first
   do while (associated(current))
      call current%p%inputs%finalize()
      call current%p%outputs%finalize()
      call current%p%dependencies%finalize()
      call current%p%all_dependencies%finalize()
      deallocate(current%p)
      current => current%next
   end do
   call self%type_node_list%finalize()
end subroutine graph_finalize

subroutine graph_save_as_dot(self,unit)
   class (type_graph),intent(in) :: self
   integer,           intent(in) :: unit

   type (type_node_list_member),pointer :: node
   type (type_node_set_member), pointer :: pnode

   write (unit,'(A)') 'digraph {'

   ! Edges
   node => self%first
   do while (associated(node))
      pnode => node%p%dependencies%first
      do while (associated(pnode))
         write (unit,'(A)') '  "'//trim(pnode%p%as_string())//'" -> "'//trim(node%p%as_string())//'";'
         pnode => pnode%next
      end do
      node => node%next
   end do

   ! Node attributes
   pnode => self%endpoints%first
   do while (associated(pnode))
      write (unit,'(A)') '  "'//trim(pnode%p%as_string())//'" [shape=box]'
      pnode => pnode%next
   end do

   write (unit,'(A)') '}'
end subroutine graph_save_as_dot

function node_as_string(node) result(string)
   class (type_node), intent(in) :: node
   character(len=attribute_length) :: string

   string = trim(node%model%get_path())//':'//trim(source2string(node%source))
end function node_as_string

subroutine graph_get_longest_path(self)
   class (type_graph),intent(in) :: self

   integer                                 :: maxlength
   type (type_node_set_member), pointer :: pnode
   type (type_node_list)             :: path

   maxlength = 0
   pnode => self%endpoints%first
   do while (associated(pnode))
      maxlength = max(maxlength,pnode%p%get_maximum_path_length())
      pnode => pnode%next
   end do

   pnode => self%endpoints%first
   do while (associated(pnode))
      call find_path_of_length(pnode%p,maxlength,path)
      pnode => pnode%next
   end do

contains

   recursive subroutine find_path_of_length(node,length,path)
      type (type_node),     intent(in) :: node
      integer,                    intent(in) :: length
      type (type_node_list),intent(inout) :: path

      type (type_node_set_member),pointer :: pnode
      type (type_node),        pointer :: path_node

      call path%append(node)
      if (length==0) then
         if (.not.associated(node%dependencies%first)) write (*,'(A)') trim(path%as_string())
      else
         pnode => node%dependencies%first
         do while (associated(pnode))
            call find_path_of_length(pnode%p,length-1,path)
            pnode => pnode%next
         end do
      end if
      path_node => path%pop()
   end subroutine find_path_of_length

end subroutine graph_get_longest_path

recursive function node_get_maximum_path_length(node) result(length)
   class (type_node), intent(in) :: node
   integer :: length

   type (type_node_set_member), pointer :: pnode

   length = 0
   pnode => node%dependencies%first
   do while (associated(pnode))
      length = max(length,1+pnode%p%get_maximum_path_length())
      pnode => pnode%next
   end do
end function node_get_maximum_path_length

function node_list_as_string(self) result(string)
   class (type_node_list),intent(in) :: self
   character(len=attribute_length) :: string
   
   type (type_node_list_member), pointer :: pnode

   string = ''
   pnode => self%first
   string = trim(pnode%p%as_string())
   do while (associated(pnode%next))
      write (string,'(A," -> ",A)') trim(pnode%next%p%as_string()),trim(string)
      pnode => pnode%next
   end do
end function node_list_as_string

end module fabm_graph
   
module fabm_job

   use fabm_types
   use fabm_driver
   use fabm_graph

   implicit none

   private

   public type_job, type_call, type_superjob
   public find_dependencies

   type type_cache_copy_command
      integer :: read_index
      integer :: write_index
   end type

   type type_call
      class (type_base_model), pointer            :: model => null()
      integer                                     :: source = source_unknown
      type (type_cache_copy_command), allocatable :: copy_commands_int(:) ! interior variables to copy from write to read cache after call completes
      type (type_cache_copy_command), allocatable :: copy_commands_hz(:)  ! horizontal variables to copy from write to read cache after call completes
      type (type_node), pointer             :: node => null()
      type (type_call), pointer                   :: next => null()
   contains
      procedure :: initialize => call_initialize
      procedure :: finalize   => call_finalize
   end type type_call

   type type_call_sequence
      type (type_call),         pointer :: first    => null()
      type (type_call),         pointer :: last     => null()
      integer                           :: source   = source_unknown

      type (type_call_sequence_set),pointer :: owner    => null()
      type (type_call_sequence),    pointer :: previous => null()
      type (type_call_sequence),    pointer :: next     => null()
   contains
      !procedure :: find        => call_list_find
      !procedure :: append      => call_list_append
      !procedure :: remove      => call_list_remove
      procedure :: print       => call_list_print
      procedure :: initialize  => call_list_initialize
      !procedure :: request_variable  => call_list_request_variable
      procedure :: finalize    => call_list_finalize
   end type type_call_sequence

   type type_call_sequence_set
      type (type_call_sequence),pointer :: first => null()
      type (type_call_sequence),pointer :: last  => null()
   contains
      procedure :: clone         => call_sequence_set_clone
      !procedure :: insert_after  => call_sequence_set_insert_after
      !procedure :: insert_before => call_sequence_set_insert_before
   end type

   type type_task_order_node
      integer                             :: source   =  source_unknown
      type (type_node_list)         :: graph_nodes
      type (type_task_order_node),pointer :: previous => null()
   contains
      procedure :: can_incorporate     => order_node_can_incorporate
      procedure :: incorporate         => order_node_incorporate
      procedure :: contains_dependency => order_node_contains_dependency
   end type

   type type_task_order
      type (type_task_order_node),pointer :: last_node => null()
      type (type_task_order),     pointer :: next      => null()
      logical                             :: valid     = .true.
   contains
      procedure :: clone     => task_order_clone
      procedure :: add       => task_order_add
      procedure :: as_string => task_order_as_string
   end type
   
   type type_task_order_set
      type (type_task_order), pointer :: first => null()
   end type

   type type_job
      integer :: domain = -1
      type (type_call_sequence) :: calls
      integer, allocatable :: prefill_type(:)
      real(rk),allocatable :: prefill_values(:)
      integer, allocatable :: prefill_index(:)
      integer, allocatable :: prefill_type_hz(:)
      real(rk),allocatable :: prefill_values_hz(:)
      integer, allocatable :: prefill_index_hz(:)
      integer, allocatable :: save_sources(:)
      integer, allocatable :: save_sources_hz(:)
      logical, allocatable :: load(:)
      logical, allocatable :: load_hz(:)
      logical, allocatable :: load_scalar(:)
      type (type_job), pointer :: next => null()
   contains
      procedure :: initialize => job_initialize
      procedure :: create_superjob => job_create_superjob
   end type

   type type_superjob
      type (type_job),   pointer  :: first => null()
      type (type_graph)           :: graph
   contains
      procedure :: initialize => superjob_initialize
      procedure :: request_variable => superjob_request_variable
      procedure :: request_variables => superjob_request_variables
      procedure :: set_next   => superjob_set_next
      procedure :: print      => superjob_print
   end type

   contains

   subroutine task_order_clone(self,clone,node,node_clone,graph_node,graph_node_clone)
      class (type_task_order), intent(in)   :: self
      type  (type_task_order), intent(out)  :: clone
      type  (type_task_order_node), target  :: node
      type  (type_task_order_node), pointer :: node_clone
      type  (type_node_list_member), target, optional  :: graph_node
      type  (type_node_list_member), pointer, optional :: graph_node_clone

      type  (type_task_order_node),pointer :: current,current_clone,next_clone

      next_clone => null()
      current => self%last_node
      do while (associated(current))
         allocate(current_clone)
         current_clone%source = current%source
         call current%graph_nodes%copy(current_clone%graph_nodes,graph_node,graph_node_clone)
         if (associated(next_clone)) then
            next_clone%previous => current_clone
         else
            clone%last_node => current_clone
         end if
         if (associated(current,node)) node_clone => current_clone
         next_clone => current_clone
         current => current%previous
      end do
   end subroutine task_order_clone

   subroutine job_initialize(self)
      class (type_job), intent(inout) :: self

      ! Initialize individual call objects
      call self%calls%initialize()

      ! Find all variables that must be written to persistent storage after this job completes.
      call create_persistent_store_commands(self%save_sources,   domain_interior)
      call create_persistent_store_commands(self%save_sources_hz,domain_horizontal)

      ! Create prefill instructions for all variables that will be written to.
      call create_prefill_commands(self%prefill_type,   self%prefill_index,   self%prefill_values,   domain_interior)
      call create_prefill_commands(self%prefill_type_hz,self%prefill_index_hz,self%prefill_values_hz,domain_horizontal)

      ! Create read cache load instructions for all input variables.
      call create_load_commands(self%load,       domain_interior)
      call create_load_commands(self%load_hz,    domain_horizontal)
      call create_load_commands(self%load_scalar,domain_scalar)

   contains

      subroutine create_prefill_commands(prefill,indices,values,domain)
         integer, intent(inout),allocatable :: prefill(:)
         integer, intent(inout),allocatable :: indices(:)
         real(rk),intent(inout),allocatable :: values(:)
         integer, intent(in)                :: domain

         integer                              :: ilast
         type (type_call),            pointer :: call_node
         type (type_output_variable), pointer :: variable_node

         ! Find the last write cache index
         ilast = 0
         call_node => self%calls%first
         do while (associated(call_node))
            variable_node => call_node%node%outputs%first
            do while (associated(variable_node))
               if (variable_node%target%prefill/=prefill_none.and.iand(variable_node%target%domain,domain)/=0) then
                  if (variable_node%target%write_indices%value<=0) &
                     call driver%fatal_error('create_prefill_commands','Variable '//trim(variable_node%target%name) &
                        //' has prefilling set, but it does not have a write cache index.')
                  if (variable_node%target%prefill==prefill_previous_value.and.variable_node%target%store_index<=0) &
                     call driver%fatal_error('create_prefill_commands','Variable '//trim(variable_node%target%name) &
                        //' has prefill==previous value, but it does not have a persistent storage index.')
                  ilast = max(ilast,variable_node%target%write_indices%value)
               end if
               variable_node => variable_node%next
            end do
            call_node => call_node%next
         end do

         if (ilast==0) return

         allocate(prefill(ilast))
         allocate(indices(ilast))
         allocate(values(ilast))
         prefill(:) = prefill_none
         indices(:) = 0

         call_node => self%calls%first
         do while (associated(call_node))
            variable_node => call_node%node%outputs%first
            do while (associated(variable_node))
               if (variable_node%target%prefill/=prefill_none.and.iand(variable_node%target%domain,domain)/=0) then
                  ilast = variable_node%target%write_indices%value
                  prefill(ilast) = variable_node%target%prefill
                  indices(ilast) = variable_node%target%store_index
                  values (ilast) = variable_node%target%missing_value
               end if
               variable_node => variable_node%next
            end do
            call_node => call_node%next
         end do
      end subroutine create_prefill_commands

      subroutine create_persistent_store_commands(commands,domain)
         integer,intent(inout),allocatable :: commands(:)
         integer,intent(in)                :: domain

         integer                              :: ilast
         type (type_call),            pointer :: call_node
         type (type_output_variable), pointer :: variable_node

         ! First find the last index in persistent storage that will be written to.
         ilast = 0
         call_node => self%calls%first
         do while (associated(call_node))
            variable_node => call_node%node%outputs%first
            do while (associated(variable_node))
               if (variable_node%copy_to_store.and.iand(variable_node%target%domain,domain)/=0) then
                  if (variable_node%target%write_indices%value<=0) &
                     call driver%fatal_error('create_prefill_commands','Variable '//trim(variable_node%target%name) &
                        //' has copy_to_store set, but it does not have a write cache index.')
                  if (variable_node%target%store_index<=0) &
                     call driver%fatal_error('create_persistent_store_commands','Variable '//trim(variable_node%target%name) &
                        //' has copy_to_store set, but it does not have a persistent storage index.')
                  ilast = max(ilast,variable_node%target%store_index)
               end if
               variable_node => variable_node%next
            end do
            call_node => call_node%next
         end do

         ! If no writes to persistent storage are needed, do not allocate the commands array.
         if (ilast==0) return

         ! Allocate the commands array (go up to the last written-to index in persistent storage only)
         allocate(commands(ilast))
         commands(:) = -1

         ! Associate indices in persistent storage with the index in the write cache at which the source variable will be found.
         call_node => self%calls%first
         do while (associated(call_node))
            variable_node => call_node%node%outputs%first
            do while (associated(variable_node))
               if (variable_node%copy_to_store.and.iand(variable_node%target%domain,domain)/=0) &
                  commands(variable_node%target%store_index) = variable_node%target%write_indices%value
               variable_node => variable_node%next
            end do
            call_node => call_node%next
         end do
      end subroutine create_persistent_store_commands

      subroutine create_load_commands(load,domain)
         logical,intent(inout),allocatable :: load(:)
         integer,intent(in)                :: domain

         integer                            :: ilast
         type (type_call),          pointer :: call_node
         type (type_input_variable),pointer :: input_variable

         ilast = 0
         call_node => self%calls%first
         do while (associated(call_node))
            input_variable => call_node%node%inputs%first
            do while (associated(input_variable))
               if (iand(input_variable%target%domain,domain)/=0) ilast = max(ilast,input_variable%target%read_indices%value)
               input_variable => input_variable%next
            end do
            call_node => call_node%next
         end do

         allocate(load(ilast))
         load(:) = .false.

         call_node => self%calls%first
         do while (associated(call_node))
            input_variable => call_node%node%inputs%first
            do while (associated(input_variable))
               if (iand(input_variable%target%domain,domain)/=0) load(input_variable%target%read_indices%value) = .true.
               input_variable => input_variable%next
            end do
            call_node => call_node%next
         end do
      end subroutine create_load_commands

   end subroutine job_initialize

   subroutine superjob_print(self)
      class (type_superjob), intent(in) :: self
      call self%graph%print()
   end subroutine

subroutine call_finalize(self)
   class (type_call), intent(inout) :: self

   if (allocated(self%copy_commands_int)) deallocate(self%copy_commands_int)
   if (allocated(self%copy_commands_hz))  deallocate(self%copy_commands_hz)
   self%node => null()
end subroutine call_finalize

subroutine call_list_initialize(call_list)
   class (type_call_sequence),intent(inout) :: call_list

   type (type_call),pointer :: node

   node => call_list%first
   do while (associated(node))
      call node%initialize()
      node => node%next
   end do
end subroutine call_list_initialize

!recursive function call_list_find(self,model,source,not_stale) result(node)
!   class (type_call_sequence),intent(in)        :: self
!   class (type_base_model),   intent(in),target :: model
!   integer,                   intent(in)        :: source
!   logical,optional,          intent(in)        :: not_stale
!
!   type (type_call),pointer :: node
!   logical                            :: not_stale_
!
!   node => self%first
!   do while (associated(node))
!      if (associated(node%model,model).and.node%source==source) return
!      node => node%next
!   end do
!   not_stale_ = .false.
!   if (present(not_stale)) not_stale_ = not_stale
!   if (.not.not_stale_.and.associated(self%previous)) node => self%previous%find(model,source,not_stale)
!end function call_list_find
!
!subroutine call_list_append(self,node)
!   class (type_call_sequence), intent(inout) :: self
!   type (type_call),pointer                  :: node
!
!   if (associated(self%first)) then
!      ! List contains one or more items - append to tail.
!      node%previous => self%last
!      if (associated(self%last)) self%last%next => node
!      self%last => node
!   else
!      ! List is empty - new node becomes both head and tail.
!      node%previous => null()
!      self%first => node
!      self%last => node
!   end if
!   node%next => null()
!end subroutine call_list_append
!
!subroutine call_list_remove(self,node)
!   class (type_call_sequence), intent(inout) :: self
!   type (type_call),pointer                  :: node
!
!   type (type_call),pointer :: current
!
!   ! Make sure the supplied node occurs in our call list to begin with.
!   current => self%first
!   do while (associated(current))
!      if (associated(current,node)) exit
!      current => current%next
!   end do
!   if (.not.associated(current)) call driver%fatal_error('call_list_remove','call list does not contain supplied node.')
!
!   if (associated(node%previous)) then
!      node%previous%next => node%next
!   else
!      self%first => node%next
!   end if
!
!   if (associated(node%next)) then
!      node%next%previous => node%previous
!   else
!      self%last => node%previous
!   end if
!
!   node%previous => null()
!   node%next => null()
!end subroutine call_list_remove

subroutine call_list_print(self)
   class (type_call_sequence), intent(in) :: self

   type (type_call),              pointer :: call_node
   type (type_output_variable),   pointer :: variable
   type (type_node_set_member),pointer :: pnode

   call_node => self%first
   do while (associated(call_node))
      write (*,'(a,": ",a)') trim(call_node%model%get_path()),trim(source2string(call_node%source))
      variable => call_node%node%outputs%first
      do while (associated(variable))
         write (*,'("   ",a,",write@",i0)',advance='no') trim(variable%target%name),variable%target%write_indices%value
         if (variable%copy_to_cache) write (*,'(",cache@",i0)',advance='no') variable%target%read_indices%value
         if (variable%copy_to_store) write (*,'(",store@",i0)',advance='no') variable%target%store_index
         if (variable%target%prefill==prefill_missing_value) write (*,'(",prefill=",g0.6)',advance='no') variable%target%missing_value
         if (variable%target%prefill==prefill_previous_value) write (*,'(",prefill=previous")',advance='no')
         write (*,*)
         pnode => variable%dependent_nodes%first
         do while (associated(pnode))
            write (*,'("     <- ",a,": ",a)') trim(pnode%p%model%get_path()),trim(source2string(pnode%p%source))
            pnode => pnode%next
         end do
         variable => variable%next
      end do
      call_node => call_node%next
   end do

end subroutine call_list_print

subroutine call_list_finalize(self)
   class (type_call_sequence), intent(inout) :: self

   type (type_call),pointer :: node,next

   node => self%first
   do while (associated(node))
      next => node%next
      call node%finalize()
      deallocate(node)
      node => next
   end do
   nullify(self%first)
   nullify(self%last)
end subroutine call_list_finalize

subroutine call_initialize(self)
   class (type_call),intent(inout) :: self

   class (type_base_model),       pointer :: parent
   class (type_model_list_node),  pointer :: model_list_node
   type (type_output_variable),   pointer :: variable_node
   type (type_node_set_member),pointer :: dependent_call
   type (type_call),              pointer :: later_call

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
   variable_node => self%node%outputs%first
   do while (associated(variable_node))
      ! For this output variable, loop over all calls that depend on it.
      dependent_call => variable_node%dependent_nodes%first
      do while (associated(dependent_call))
         ! For this dependent call, establish whether it is part of the same job, or of a later job.
         later_call => self%next
         do while (associated(later_call))
            if (associated(later_call,dependent_call%p)) exit
            later_call => later_call%next
         end do
         if (associated(later_call)) then
            ! The dependent call is part of the same job. Therefore the output needs to be copied to the read cache.
            variable_node%copy_to_cache = .true.
         else
            ! The dependent call is part of a later job. Therefore the output needs to be copied to the persistent store.
            variable_node%copy_to_store = .true.
         end if
         dependent_call => dependent_call%next
      end do
      variable_node => variable_node%next
   end do

   call create_cache_copy_commands(self%copy_commands_int, domain_interior)
   call create_cache_copy_commands(self%copy_commands_hz,  domain_horizontal)

contains

   subroutine create_cache_copy_commands(commands,domain)
      type (type_cache_copy_command), allocatable, intent(inout) :: commands(:)
      integer,                                     intent(in)    :: domain

      type (type_output_variable), pointer :: variable
      integer                              :: i, n, maxwrite

      n = 0
      maxwrite = -1
      variable => self%node%outputs%first
      do while (associated(variable))
         if (variable%copy_to_cache.and.iand(variable%target%domain,domain)/=0) then
            if (variable%target%write_indices%is_empty()) &
               call driver%fatal_error('call_list_node_initialize::create_cache_copy_commands', &
                  'BUG: '//trim(variable%target%name)//' cannot be copied from write to read cache because it lacks a write cache index.')
            if (variable%target%read_indices%is_empty()) &
               call driver%fatal_error('call_list_node_initialize::create_cache_copy_commands', &
                  'BUG: '//trim(variable%target%name)//' cannot be copied from write to read cache because it lacks a read cache index.')
            n = n + 1
            maxwrite = max(maxwrite,variable%target%write_indices%value)
         end if
         variable => variable%next
      end do

      ! Create list of cache copy commands (source index in write cache, target index in read cache)
      ! They are sorted by write index to make (write) access to memory more predictable.
      allocate(commands(n))
      n = 0
      do i=1,maxwrite
         variable => self%node%outputs%first
         do while (associated(variable))
            if (variable%copy_to_cache.and.iand(variable%target%domain,domain)/=0.and.variable%target%write_indices%value==i) then
               n = n + 1
               commands(n)%read_index = variable%target%read_indices%value
               commands(n)%write_index = variable%target%write_indices%value
               exit
            end if
            variable => variable%next
         end do
      end do
   end subroutine create_cache_copy_commands

end subroutine call_initialize
!
!recursive function ordering_add(self,node)
!   class (type_task_order),  target,intent(inout) :: list
!   class (type_node),target,intent(in)    :: node
!   integer,                          intent(in)    :: source
!   logical,optional,                 intent(in)    :: not_stale
!   type (type_call),  pointer                      :: node
!
!   integer :: list_source
!   type (type_link),pointer :: link
!
!   ! Check whether we have already processed this call. If so, incorporation was successful.
!   node => list%find(model,source,not_stale)
!   if (associated(node)) return
!
!   ! We need to add this call. Try it in the current call list; if that fails, try previous.
!   list_source = list%source ! copy the list source, because it may be replaced, and we do not want the replacement unless all else goes ok
!   if (.not.is_source_compatible(list_source,source)) then
!      ! This call cannot be processed as part of the current call sequence.
!      ! Try the previous sequence. If that does not exist, we failed to incorporate the
!      ! call in our existing set, so return with success=.false.
!      if (associated(list%previous)) node => incorporate_call(list%previous,model,source,not_stale)
!      return
!   end if
!
!   ! Create object to represent this call
!   node => create_call(model,source)
!
!   ! Loop over all input variables, and if they are written by another model call, make sure that model call is done first.
!   link => model%links%first
!   do while (associated(link))
!      if (is_dependency(model,source,link)) then
!         if (.not.incorporate_output(list,link%target,node)) exit
!      end if
!      link => link%next
!   end do
!
!   if (.not.associated(link)) then
!      ! We succeeded to incoporate this call and all its dependencies, so we're happy - add ourselves to the call list.
!      ! NB: call for all our dependencies have already been appended, and thus precede ours.
!      call list%append(node)
!      list%source = list_source
!   else
!      ! We failed to incoportae one of this call's dependencies - undo everything.
!      call node%finalize()
!      deallocate(node)
!   end if
!end function incorporate_call
!
!recursive function incorporate_output(self,variable,requester) result(success)
!   class (type_call_sequence),target,intent(inout) :: self
!   type (type_internal_variable),    intent(in)    :: variable
!   type (type_call),target                         :: requester
!
!   logical :: success
!   type (type_output_variable),pointer :: variable_node
!
!   type (type_call),pointer :: node
!
!   success = .false.
!   if (variable%source==source_unknown) then
!      ! This variable is either written by do_surface or do_bottom - which one of these two APIs is unknown.
!      node => incorporate_call(self,variable%owner,source_do_surface)
!      if (.not.associated(node)) return
!      variable_node => call_register_output(node,variable,requester)
!      node => incorporate_call(self,variable%owner,source_do_bottom)
!      if (.not.associated(node)) return
!      variable_node => call_register_output(node,variable,requester)
!   elseif (variable%source/=source_none) then
!      ! This variable is written by a known BGC API [is is not a constant indicated by source_none]
!      node => incorporate_call(self,variable%owner,variable%source)
!      if (.not.associated(node)) return
!      variable_node => call_register_output(node,variable,requester)
!   end if
!   success = .true.
!end function incorporate_output
!
!subroutine call_options_add(self,set)
!   class (type_call_options),intent(inout) :: self
!   class (type_call_sequence_set),target :: set
!
!   class (type_call_option),pointer :: node
!
!   node => self%first
!   do while (associated(node))
!      if (associated(node%set,set)) return
!      node => node%next
!   end do
!   allocate(node)
!   node%set => set
!   node%next => self%first
!   self%first => node
!end subroutine call_options_add

function order_node_has_graph_node(self,graph_node) result(found)
   class (type_task_order_node), target, intent(in) :: self
   type (type_node),             target, intent(in) :: graph_node
   logical                                          :: found

   class (type_task_order_node), pointer :: task
   type (type_node_list_member), pointer :: pnode

   found = .true.
   task => self
   do while (associated(task))
      pnode => self%graph_nodes%first
      do while (associated(pnode))
         if (associated(pnode%p,graph_node)) return
         pnode => pnode%next
      end do
      task => task%previous
   end do
   found = .false.
end function order_node_has_graph_node

recursive function order_node_can_incorporate(self,graph_node,next_graph_node) result(success)
   class (type_task_order_node), intent(in)      :: self
   type (type_node),target,      intent(in)      :: graph_node
   type (type_node_list_member),target, optional :: next_graph_node
   logical                                       :: success

   integer                             :: source
   type (type_node_set_member),pointer :: dependency

   source = self%source
   success = is_source_compatible(source,graph_node%source)
   success = success .and. (source==self%source.or.self%source==source_unknown)  ! If incorporation succeeded by changing the source type, we are not satisfied
   if (.not.success) then
      ! Cannot incorporate in current task - try previous
      if (associated(self%previous).and..not.self%contains_dependency(graph_node,next_graph_node)) &
         success = self%previous%can_incorporate(graph_node)
      return
   end if

   ! We could incorporate the graph node. Now try its dependencies
   dependency => graph_node%dependencies%first
   do while (associated(dependency))
      if (.not.order_node_has_graph_node(self,dependency%p)) then
         success = self%can_incorporate(dependency%p,next_graph_node)
         if (.not.success) exit
      end if
      dependency => dependency%next
   end do
end function order_node_can_incorporate

recursive subroutine order_node_incorporate(self,graph_node,next_graph_node)
   class (type_task_order_node), intent(inout)   :: self
   type (type_node),target,      intent(in)      :: graph_node
   type (type_node_list_member),target, optional :: next_graph_node

   integer                              :: source
   logical                              :: success
   type (type_node_set_member), pointer :: dependency
   type (type_node_list_member),pointer :: dep_next_graph_node

   source = self%source
   success = is_source_compatible(source,graph_node%source)
   success = success .and. (source==self%source.or.self%source==source_unknown)  ! If incorporation succeeded by changing the source type, we are not satisfied
   if (.not.success) then
      ! Cannot incorporate in current task - try previous
      if (.not.associated(self%previous)) call driver%fatal_error('order_node_incorporate','BUG: no previous node available.')
      if (self%contains_dependency(graph_node,next_graph_node)) call driver%fatal_error('order_node_incorporate','BUG: cannot pass incompatible current node because it provides dependencies.')
      call self%previous%incorporate(graph_node)
      return
   end if

   if (present(next_graph_node)) then
      call self%graph_nodes%insert_before(next_graph_node,graph_node)
      dep_next_graph_node => next_graph_node%previous
   else
      call self%graph_nodes%append(graph_node)
      dep_next_graph_node => self%graph_nodes%last
   end if

   ! We could incorporate the graph node. Now try its dependencies
   dependency => graph_node%dependencies%first
   do while (associated(dependency))
      if (.not.order_node_has_graph_node(self,dependency%p)) call self%incorporate(dependency%p,dep_next_graph_node)
      dependency => dependency%next
   end do
end subroutine order_node_incorporate

function order_node_contains_dependency(self,graph_node,next_graph_node) result(found)
   class (type_task_order_node), target             :: self
   type (type_node),            intent(in)          :: graph_node
   type (type_node_list_member),intent(in),optional :: next_graph_node
   logical :: found

   type (type_node_list_member), pointer :: pnode

   found = .true.
   if (present(next_graph_node)) then
      pnode => next_graph_node%previous
   else
      pnode => self%graph_nodes%last
   end if
   do while (associated(pnode))
      if (graph_node%all_dependencies%contains(pnode%p)) then
         write (*,*) 'task node found among graph node '//trim(graph_node%as_string())//' dependencies: '//trim(pnode%p%as_string())
         return
      end if
      pnode => pnode%previous
   end do
   found = .false.
end function order_node_contains_dependency

recursive subroutine task_order_add(self,task,graph_node,options,try_append,next_graph_node)
   class (type_task_order),     target,intent(inout) :: self
   type (type_task_order_node), target,intent(inout) :: task
   type (type_node),            target,intent(in)    :: graph_node
   type (type_task_order_set),         intent(inout) :: options
   type (type_node_list_member),target, optional     :: next_graph_node
   logical,                            intent(in)    :: try_append

   class (type_task_order),     pointer :: prepend_order, forward_order, integrate_order, append_order
   type (type_task_order_node), pointer :: prepend_task, forward_task, integrate_task, append_task, new_task
   type (type_node_list_member),pointer :: integrate_next_graph_node
   integer                              :: new_source
   logical                              :: dummy, try_prepend, try_integrate, try_forward, used_original_option

#ifndef NDEBUG
   type (type_node_list_member),pointer :: pnode
   if (.not.self%valid) call driver%fatal_error('task_order_add','BUG: called on invalid task order')
   new_task => self%last_node
   do while (associated(new_task))
      if (associated(new_task,task)) exit
      new_task => new_task%previous
   end do
   if (.not.associated(new_task)) call driver%fatal_error('task_order_add','BUG: task is not part of own task order')
   if (present(next_graph_node)) then
      pnode => task%graph_nodes%first
      do while (associated(pnode))
         if (associated(pnode,next_graph_node)) exit
         pnode => pnode%next
      end do
      if (.not.associated(pnode)) call driver%fatal_error('task_order_add','BUG: next_graph_node is not part of task''s own graph nodes.')
   end if
#endif

   if (order_node_has_graph_node(task,graph_node)) return

   ! First try to incorporate the current call and all its dependencies in the current task order.
   ! If that succeeds, we are done (no need to generate alternative task orders)
   if (task%can_incorporate(graph_node,next_graph_node)) then
      call task%incorporate(graph_node,next_graph_node)
      return
   end if

   ! Whether the original ordering option ("self") has already been used.
   ! If that is the case, any further possibilities will require cloning of the original.
   used_original_option = .false.

   ! Integrate in current task if sources of task and graph node are compatible.
   new_source = task%source
   try_integrate = is_source_compatible(new_source,graph_node%source)

   ! Prepend as new task if the current task cannot fully (i.e., with original source intact) integrate the call.
   try_prepend = .not.try_integrate .or. (task%source/=new_source .and. task%source/=source_unknown)

   ! Forward if there is a task to forward to (i.e., we are not at the head of the list).
   try_forward = associated(task%previous)

   if (try_prepend.or.try_forward) then
      ! First determine if we CAN move backward from the current task.
      ! This is not possible if the current task is responsible for computing one of the graph node's dependencies.
      if (task%contains_dependency(graph_node,next_graph_node)) then
         try_prepend = .false.
         try_forward = .false.
      end if
   end if

   ! Set up separate task order options for each supported operation [integrate, prepend, forward, append]
   if (try_integrate) call get_option(integrate_order,integrate_task,integrate_next_graph_node)
   if (try_prepend  ) call get_option(prepend_order,  prepend_task)
   if (try_forward  ) call get_option(forward_order,  forward_task)
   if (try_append   ) call get_option(append_order,   append_task)

   ! Option: integrate in current task (only if sources are compatible).
   if (try_integrate) then
      integrate_task%source = new_source
      if (present(next_graph_node)) then
         ! Insert call before the specified lastly-inserted one.
         call integrate_task%graph_nodes%insert_before(integrate_next_graph_node,graph_node)
         integrate_next_graph_node => integrate_next_graph_node%previous
      else
         ! Append call to the end.
         call integrate_task%graph_nodes%append(graph_node)
         integrate_next_graph_node => integrate_task%graph_nodes%last
      end if
      call process_dependencies(integrate_order,integrate_task,integrate_next_graph_node)
   end if

   ! Option: prepend a new node specifically for the call.
   if (try_prepend) then
      allocate(new_task)
      new_task%previous => prepend_task%previous
      prepend_task%previous => new_task
      dummy = is_source_compatible(new_task%source,graph_node%source)  ! make sure the new source is valid [not source_unknown]
      call new_task%graph_nodes%append(graph_node)
      call process_dependencies(prepend_order,new_task,new_task%graph_nodes%last)
   end if

   ! Option: forward to previous (only possible if not at head of list)
   if (try_forward) call forward_order%add(forward_task%previous,graph_node,options,.false.)

   ! Option: append a new node specifically for the call (only if explicitly requested).
   if (try_append) then
      if (.not.associated(append_order%last_node,append_task)) call driver%fatal_error('task_order_add','BUG append should be done only on tail')
      allocate(new_task)
      new_task%previous => append_task
      append_order%last_node => new_task
      dummy = is_source_compatible(new_task%source,graph_node%source)  ! make sure the new source is valid [not source_unknown]
      call new_task%graph_nodes%append(graph_node)
      call process_dependencies(append_order,new_task,new_task%graph_nodes%last)
   end if
   
   if (.not.used_original_option) self%valid = .false.

contains

   subroutine get_option(new_order,new_task,new_next_graph_node)
      class (type_task_order),     pointer :: new_order
      type (type_task_order_node), pointer :: new_task
      type (type_node_list_member),pointer, optional :: new_next_graph_node

      if (used_original_option) then
         ! We used the original option already. Clone it to create a new one.
         allocate(new_order)
         call self%clone(new_order,task,new_task,next_graph_node,new_next_graph_node)

         ! Prepend cloned option to set of existing options.
         new_order%next => options%first
         options%first => new_order
      else
         ! The original option has not been used yet, so just point to that.
         new_order => self
         new_task => task
         if (present(new_next_graph_node).and.present(next_graph_node)) new_next_graph_node => next_graph_node
         used_original_option = .true.
      end if
   end subroutine get_option

   recursive subroutine process_dependencies(order,task,next_graph_node)
      class (type_task_order), target,intent(inout) :: order
      type (type_task_order_node), pointer          :: task
      type (type_node_list_member), target          :: next_graph_node

      type (type_node_set_member), pointer :: dependency

      dependency => graph_node%dependencies%first
      do while (associated(dependency))
         if (order%valid) call order%add(task,dependency%p,options,.false.,next_graph_node)
         dependency => dependency%next
      end do
   end subroutine process_dependencies

end subroutine task_order_add

!function call_sequence_set_insert_after(self,after) result(list)
!   class (type_call_sequence_set),intent(inout),target :: self
!   class (type_call_sequence),target :: after
!
!   class (type_call_sequence),pointer :: list
!
!   allocate(list)
!   list%owner => self
!   list%previous => after
!   list%next => after%next
!   if (associated(after%next)) then
!      ! Inserting between two existing sequences
!      after%next%previous => list
!   else
!      ! Appending to the tail
!      self%last => list
!   end if
!   after%next => list
!end function call_sequence_set_insert_after
!
!function call_sequence_set_insert_before(self,before) result(list)
!   class (type_call_sequence_set),intent(inout),target :: self
!   class (type_call_sequence),target :: before
!
!   class (type_call_sequence),pointer :: list
!
!   allocate(list)
!   list%owner => self
!   list%next => before
!   list%previous => before%previous
!   if (associated(before%previous)) then
!      ! Inserting between two existing sequences
!      before%previous%next => list
!   else
!      ! Prepending to the head
!      self%first => list
!   end if
!   before%previous => list
!end function call_sequence_set_insert_before

subroutine superjob_request_variable(self,variable,copy_to_cache,copy_to_store,not_stale)
   class (type_superjob),target,     intent(inout)          :: self
   type (type_internal_variable),    intent(in)             :: variable
   logical,                          intent(in),   optional :: copy_to_cache
   logical,                          intent(in),   optional :: copy_to_store
   logical,                          intent(in),   optional :: not_stale

   type (type_node_list) :: outer_calls

   ! If this variable is not written [but a field provided by the host or a state variable], return immediately.
   if (variable%write_indices%is_empty().or.variable%source==source_none) return

   call self%graph%add_variable(variable,outer_calls,copy_to_cache,copy_to_store,not_stale)

   call outer_calls%finalize()
end subroutine superjob_request_variable

subroutine superjob_request_variables(self,link_list,copy_to_cache,copy_to_store,not_stale)
   class (type_superjob),     intent(inout)       :: self
   type (type_link_list),     intent(in)          :: link_list
   logical,                   intent(in),optional :: copy_to_cache
   logical,                   intent(in),optional :: copy_to_store
   logical,                   intent(in),optional :: not_stale

   type (type_link), pointer :: link

   link => link_list%first
   do while (associated(link))
      call self%request_variable(link%target,copy_to_cache,copy_to_store,not_stale)
      link => link%next
   end do
end subroutine superjob_request_variables

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
         call driver%fatal_error('find_dependencies','circular dependency found: '//trim(chain(2:))//' '//trim(self%get_path()))
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

function job_create_superjob(self) result(superjob)
   class (type_job),intent(inout) :: self
   type (type_superjob)           :: superjob

   !type (type_job), pointer :: current_job
   !type (type_call),pointer :: call_node, next_call_node
   !integer                  :: source
   !
   !allocate(superjob%first)
   !current_job => superjob%first
   !
   !source = source_unknown
   !call_node => self%calls%first
   !do while (associated(call_node))
   !   write (*,*) 'processing call '//trim(call_node%model%get_path())//':'//trim(source2string(call_node%source))
   !   next_call_node => call_node%next
   !
   !   if (.not.is_source_compatible(source,call_node%source)) then
   !      ! New call list needed
   !      write (*,*) '  creating new job - old list source = '//trim(source2string(source))
   !      source = call_node%source
   !      allocate(current_job%next)
   !      call current_job%set_next(current_job%next)
   !      current_job => current_job%next
   !   end if

      ! Add to current call list
      !call self%calls%remove(call_node)
      !call current_job%calls%append(call_node)
   !   call_node => next_call_node
   !end do
end function job_create_superjob

subroutine superjob_select_order(self)
   class (type_superjob), intent(inout) :: self

   type (type_task_order_set)             :: options
   type (type_task_order),        pointer :: order,best_order
   type (type_task_order_node),   pointer :: order_node
   type (type_node_set_member),pointer :: endpoint
   logical                                :: dummy
   integer                                :: n,cost,mincost

   ! Set up first task order (its only member contains one task only, specific to the first endpoint in the graph)
   allocate(order)
   options%first => order
   allocate(order_node)
   order%last_node => order_node
   endpoint => self%graph%endpoints%first
   dummy = is_source_compatible(order_node%source,endpoint%p%source)

   ! Add each endpoint to all possible task orders identified thus far.
   do while (associated(endpoint))
      order => options%first
      do while (associated(order))
         if (order%valid) call order%add(order%last_node,endpoint%p,options,.true.)
         order => order%next
      end do
      endpoint => endpoint%next
   end do

   ! Find task order with minimum cost
   order => options%first
   n = 0
   best_order => null()
   mincost = huge(mincost)
   do while (associated(order))
      if (order%valid) then
         n = n + 1
         cost = 0
         order_node => order%last_node
         do while (associated(order_node))
            cost = cost + 1
            order_node => order_node%previous
         end do
         write (*,*) trim(order%as_string())
         if (cost<mincost) then
            mincost = cost
            best_order => order
         end if
      end if
      order => order%next
   end do
   if (associated(best_order)) write (*,*) 'Selected order: '//trim(best_order%as_string())
end subroutine superjob_select_order

function task_order_as_string(self) result(string)
   class (type_task_order), intent(in)  :: self
   character(len=attribute_length)      :: string

   type (type_task_order_node), pointer :: order_node

   string = ''
   order_node => self%last_node
   write (string,'(i0)') order_node%source
   do while (associated(order_node%previous))
      write (string,'(i0,",",a)') order_node%previous%source,trim(string)
      order_node => order_node%previous
   end do
end function task_order_as_string

subroutine superjob_initialize(self)
   class (type_superjob), intent(inout) :: self

   type (type_job), pointer :: job
   integer :: unit,ios

   open(newunit=unit,file='graph.gv',action='write',status='replace',iostat=ios)
   call self%graph%save_as_dot(unit)
   close(unit)

   call self%graph%get_longest_path()

   call superjob_select_order(self)

   job => self%first
   do while (associated(job))
      call job%initialize()
      job => job%next
   end do
end subroutine superjob_initialize

logical function is_source_compatible(current_source,new_source)
   integer,intent(inout) :: current_source
   integer,intent(in)    :: new_source

   integer :: real_new_source

   real_new_source = new_source
   if (real_new_source==source_get_vertical_movement) real_new_source = source_do

   if (current_source==source_unknown) then
      is_source_compatible = .true.
      current_source = real_new_source
      return
   end if

   select case (real_new_source)
   case (source_do,source_do_column)
      is_source_compatible = current_source==real_new_source
   case (source_do_bottom)
      is_source_compatible = current_source==source_do_horizontal.or.current_source==source_do_bottom
      if (current_source==source_do_horizontal) current_source = source_do_bottom
   case (source_do_surface)
      is_source_compatible = current_source==source_do_horizontal.or.current_source==source_do_surface
      if (current_source==source_do_horizontal) current_source = source_do_surface
   case (source_do_horizontal)
      is_source_compatible = current_source==source_do_horizontal.or.current_source==source_do_surface.or.current_source==source_do_bottom
   case default
      call driver%fatal_error('is_source_compatible','unknown source')
   end select
end function is_source_compatible

subroutine superjob_set_next(self,next)
   class (type_superjob), intent(inout), target :: self
   type (type_superjob),  intent(inout), target :: next

   if (associated(next%graph%previous)) &
      call driver%fatal_error('superjob::set_next','A previous superjob has already been registered.')

   next%graph%previous => self%graph
end subroutine superjob_set_next

function call_sequence_set_clone(self) result(clone)
   class (type_call_sequence_set),intent(in) :: self
   class (type_call_sequence_set),pointer :: clone

   allocate(clone)
end function call_sequence_set_clone

end module fabm_job

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

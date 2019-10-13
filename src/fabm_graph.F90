#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_graph: derived types to describe the directed acyclic graph (DAG) that characterizes dependencies between models.
! Each graph node represents a call, that is, a specific combination of a model and one of its APIs ("source").
! Graph nodes maintain both pointers to their dependencies (a set of other nodes, maintained at the node level),
! and to the nodes that depend on them (maintained separately per output variable).
!
! !INTERFACE:
module fabm_graph

   use fabm_types
   use fabm_driver

   implicit none

   public type_graph, type_node, type_output_variable_set, type_output_variable_set_node, type_node_set_member, type_node_list, type_node_list_member
   public type_output_variable, source2operation, type_input_variable_set_node

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
      procedure :: append    => node_list_append
      procedure :: pop       => node_list_pop
      procedure :: remove    => node_list_remove
      procedure :: finalize  => node_list_finalize
      procedure :: find      => node_list_find
      procedure :: find_node => node_list_find_node
      procedure :: check     => node_list_check
   end type

   type type_output_variable
      type (type_internal_variable), pointer :: target          => null()
      type (type_node_set)                   :: dependent_nodes
      logical                                :: copy_to_cache   = .false.
      logical                                :: copy_to_store   = .false.
   end type

   type type_output_variable_set_node
      type (type_output_variable),          pointer :: p    => null()
      type (type_output_variable_set_node), pointer :: next => null()
   end type

   type type_output_variable_set
      type (type_output_variable_set_node), pointer :: first => null()
   contains
      procedure :: add                 => output_variable_set_add
      procedure :: add_output_variable => output_variable_set_add_output_variable
      procedure :: finalize            => output_variable_set_finalize
   end type

   type type_input_variable
      type (type_internal_variable), pointer :: target  => null()
      type (type_output_variable_set)        :: sources
   end type

   type type_input_variable_set_node
      type (type_input_variable),          pointer :: p    => null()
      type (type_input_variable_set_node), pointer :: next => null()
   end type

   type type_input_variable_set
      type (type_input_variable_set_node), pointer :: first => null()
   contains
      procedure :: add      => input_variable_set_add
      procedure :: finalize => input_variable_set_finalize
   end type

   type type_node
      class (type_base_model), pointer :: model => null()
      integer                          :: source = source_unknown
      type (type_input_variable_set)   :: inputs               ! input variables (irrespective of their source - can be constants, state variables, host-provided, or diagnostics computed by another node)
      type (type_node_set)             :: dependencies         ! direct dependencies
      type (type_output_variable_set)  :: outputs              ! output variables that a later called model requires
   contains
      procedure :: as_string => node_as_string
   end type

   type type_graph_set_member
      type (type_graph), pointer :: p => null()
      type (type_graph_set_member), pointer :: next => null()
   end type

   type type_graph_set
      type (type_graph_set_member), pointer :: first => null()
   end type

   type,extends(type_node_list) :: type_graph
      type (type_graph_set)     :: previous
      type (type_graph_set)     :: next
      integer                   :: operation = source_unknown
      logical                   :: frozen = .false.
   contains
      procedure :: connect      => graph_connect
      procedure :: add_call     => graph_add_call
      procedure :: add_variable => graph_add_variable
      procedure :: print        => graph_print
      procedure :: save_as_dot  => graph_save_as_dot
   end type

contains

function output_variable_set_add(self,variable) result(output_variable)
   class (type_output_variable_set),intent(inout) :: self
   type (type_internal_variable),target           :: variable

   type (type_output_variable_set_node), pointer :: node
   type (type_output_variable),          pointer :: output_variable

   ! Check if this variable already exists.
   node => self%first
   do while (associated(node))
      if (associated(node%p%target, variable)) then
         output_variable => node%p
         return
      end if
      node => node%next
   end do

   ! Create a new variable object and prepend it to the list.
   allocate(output_variable)
   output_variable%target => variable
   call self%add_output_variable(output_variable)
end function output_variable_set_add

subroutine output_variable_set_add_output_variable(self, output_variable)
   class (type_output_variable_set), intent(inout) :: self
   type (type_output_variable), target             :: output_variable

   class (type_output_variable_set_node), pointer :: node

   allocate(node)
   node%p => output_variable
   node%next => self%first
   self%first => node
end subroutine output_variable_set_add_output_variable

function input_variable_set_add(self, variable) result(input_variable)
   class (type_input_variable_set), intent(inout) :: self
   type (type_internal_variable), target          :: variable

   type (type_input_variable_set_node), pointer :: node
   type (type_input_variable),          pointer :: input_variable

   node => self%first
   do while (associated(node))
      if (associated(node%p%target, variable)) then
         input_variable => node%p
         call input_variable%sources%finalize(owner=.false.)
         return
      end if
      node => node%next
   end do

   allocate(input_variable)
   input_variable%target => variable

   allocate(node)
   node%p => input_variable
   node%next => self%first
   self%first => node
end function input_variable_set_add

subroutine output_variable_set_finalize(self, owner)
   class (type_output_variable_set), intent(inout) :: self
   logical,                          intent(in)    :: owner

   type (type_output_variable_set_node),pointer :: node, next

   node => self%first
   do while (associated(node))
      next => node%next
      if (owner) then
         call node%p%dependent_nodes%finalize()
         deallocate(node%p)
      end if
      deallocate(node)
      node => next
   end do
   self%first => null()
end subroutine output_variable_set_finalize

subroutine input_variable_set_finalize(self)
   class (type_input_variable_set), intent(inout) :: self

   type (type_input_variable_set_node), pointer :: node, next

   node => self%first
   do while (associated(node))
      next => node%next
      call node%p%sources%finalize(owner=.false.)
      deallocate(node%p)
      deallocate(node)
      node => next
   end do
   self%first => null()
end subroutine input_variable_set_finalize

subroutine graph_print(self)
   class (type_graph), intent(in) :: self

   type (type_node_list_member),pointer :: node
   type (type_output_variable_set_node), pointer :: variable
   type (type_node_set_member), pointer :: pnode

   node => self%first
   do while (associated(node))
      write (*,'(a,": ",a)') trim(node%p%model%get_path()),trim(source2string(node%p%source))
      variable => node%p%outputs%first
      do while (associated(variable))
         write (*,'("   ",a,",write@",i0)',advance='no') trim(variable%p%target%name),variable%p%target%write_indices%value
         if (variable%p%copy_to_cache) write (*,'(",cache@",i0)',advance='no') variable%p%target%read_indices%value
         if (variable%p%copy_to_store) write (*,'(",store@",i0)',advance='no') variable%p%target%store_index
         write (*,*)
         pnode => variable%p%dependent_nodes%first
         do while (associated(pnode))
            write (*,'("     <- ",a,": ",a)') trim(pnode%p%model%get_path()),trim(source2string(pnode%p%source))
            pnode => pnode%next
         end do
         variable => variable%next
      end do
      node => node%next
   end do

end subroutine graph_print

recursive subroutine find_node(self, model, source, graph, node)
   class (type_graph),      intent(inout), target :: self
   class (type_base_model), intent(in),    target :: model
   integer,                 intent(in)            :: source
   type (type_graph), pointer :: graph
   type (type_node), pointer :: node

   type (type_graph_set_member), pointer :: member

   graph => self
   node => self%find(model, source)
   if (associated(node)) return

   member => self%next%first
   do while (associated(member))
      call find_node(member%p, model, source, graph, node)
      if (associated(node)) return
      member => member%next
   end do

   graph => null()
end subroutine

recursive function graph_has_descendant(self, graph) result(has_descendant)
   type (type_graph), pointer :: self
   type (type_graph), pointer :: graph
   logical :: has_descendant

   type (type_graph_set_member), pointer :: member

   has_descendant = .true.
   if (associated(self, graph)) return
   member => self%next%first
   do while (associated(member))
      if (graph_has_descendant(member%p, graph)) return
      member => member%next
   end do
   has_descendant = .false.
end function

recursive function graph_add_call(self, model, source, outer_calls) result(node)
   class (type_graph),     target,intent(inout) :: self
   class (type_base_model),target,intent(in)    :: model
   integer,                       intent(in)    :: source
   type (type_node_list),  target,intent(inout) :: outer_calls
   type (type_node), pointer :: node

   type (type_node_list_member), pointer :: pnode
   character(len=2048)                   :: chain
   type (type_graph),            pointer :: root_graph, owner_graph, target_graph
   integer                               :: operation
   type (type_link),             pointer :: link
   logical                               :: same_source
   type (type_input_variable),   pointer :: input_variable

   ! Circular dependency check:
   ! Search the list of outer model calls, i.e., calls that [indirectly] request the current call.
   ! If the current call is already on this list, it is indirectly calling itself: there is a circular dependency.
   pnode => outer_calls%find_node(model, source)
   if (associated(pnode)) then
      ! Circular dependency found - report as fatal error.
      chain = ''
      do while (associated(pnode))
         chain = trim(chain)//' '//trim(pnode%p%as_string())//' ->'
         pnode => pnode%next
      end do
      call driver%fatal_error('graph::add_call','circular dependency found: '//trim(chain(2:))//' '//trim(node%as_string()))
   end if

   ! By default we add the call to the current graph (but if necessary we will target an ancestor instead)
   target_graph => self

   ! Check if this node is already in a graph of another job (recursive search from the root graph/very first job).
   root_graph => target_graph
   do while (associated(root_graph%previous%first))
      root_graph => root_graph%previous%first%p
   end do
   call find_node(root_graph, model, source, owner_graph, node)

   if (associated(node)) then
      ! If the graph that contains the target call is one of our ancestors, or us (i.e., it is scheduled to run before us or together with us),
      ! we are done - return.
      if (graph_has_descendant(owner_graph, target_graph)) return

      ! This node is in sibling graph (a job scheduled to run in parallel) or in a descendent graph (a job scheduled to run later).
      ! It needs to be moved to the last common ancestor of the currently targeted graph and the graph that already contains the call.
      ! First find the last common ancestor (i.e., the graph that needs to receive the node)
      do while (.not. graph_has_descendant(target_graph, owner_graph))
         _ASSERT_(.not. associated(target_graph%previous%first%next), 'graph::add_call', 'Multiple ancestors found when trying to find common')
         target_graph => target_graph%previous%first%p
      end do

      ! Now remove the node from the graph that currently contains it.
      call owner_graph%remove(node)
   else
      allocate(node)
      node%model => model
      node%source = source
   end if

   ! Find an ancestor graph (earlier scheduled job) with a compatible operation.
   operation = source2operation(source)
   do while (target_graph%operation /= source_unknown .and. target_graph%operation /= operation .and. .not. ((target_graph%operation == source_do_bottom .or. target_graph%operation == source_do_surface) .and. operation == source_do_horizontal))
      _ASSERT_(.not. associated(target_graph%previous%first%next), 'graph::add_call', 'Multiple ancestors found when trying to find ancestor with compatible operation')
      target_graph => target_graph%previous%first%p
   end do

   if (target_graph%frozen) call driver%fatal_error('graph_add_call','Target graph is frozen; no calls can be added.')

   ! First add this call to the list of requesting calls [a list of all calls higher on the call stack]
   ! This forbids any indirect dependency on this call, as such would be a circular dependency.
   call outer_calls%append(node)

   link => model%links%first
   do while (associated(link))
      if (index(link%name, '/') == 0 .and. associated(link%original%read_index)) then
         ! This is the model's own variable (not inherited from child model) and the model itself originally requested read access to it.
         _ASSERT_(.not. associated(link%target%write_owner), 'graph::add_call', 'BUG: required input variable is co-written.')
         input_variable => node%inputs%add(link%target)
         same_source = link%target%source == source .or. (link%target%source == source_unknown .and. (source == source_do_surface .or. source == source_do_bottom))
         if (.not. (associated(link%target%owner, model) .and. same_source)) call target_graph%add_variable(link%target, outer_calls, input_variable%sources, caller=node)
      end if
      link => link%next
   end do

   ! Remove node from the list of outer calls and add it to the graph instead.
   node => outer_calls%pop()
   call target_graph%append(node)
end function graph_add_call

recursive subroutine graph_add_variable(self, variable, outer_calls, variable_set, copy_to_store, caller)
   class (type_graph),              intent(inout) :: self
   type (type_internal_variable),   intent(in)    :: variable
   type (type_node_list),    target,intent(inout) :: outer_calls
   type (type_output_variable_set), intent(inout) :: variable_set
   logical,         optional,       intent(in)    :: copy_to_store
   type (type_node),optional,target,intent(inout) :: caller

   type (type_variable_node), pointer :: variable_node

   if (self%frozen) call driver%fatal_error('graph_add_variable','Graph is frozen; no variables can be added.')

   ! If this variable is not an output of some model (e.g., a state variable or external dependency), no call is needed.
   if (variable%write_indices%is_empty()) return

   if (variable%source == source_unknown) then
      ! This variable is either written by do_surface or do_bottom - which one of these two APIs is unknown.
      call add_call(source_do_surface)
      call add_call(source_do_bottom)
   elseif (variable%source /= source_constant .and. variable%source /= source_state .and. variable%source /= source_external) then
      ! This variable is written by a known BGC API [is is not constant/part of state/host- or user-provided]
      call add_call(variable%source)
   end if

   ! Automatically request additional value contributions (for reduction operators that accept in-place modification of the variable value)
   variable_node => variable%cowriters%first
   do while (associated(variable_node))
      call self%add_variable(variable_node%target, outer_calls, variable_set, copy_to_store, caller)
      variable_node => variable_node%next
   end do

contains

   recursive subroutine add_call(source)
      integer, intent(in) :: source

      type (type_node),           pointer :: node
      type (type_output_variable),pointer :: output_variable

      node => self%add_call(variable%owner, source, outer_calls)
      output_variable => node%outputs%add(variable)
      if (present(copy_to_store)) output_variable%copy_to_store = output_variable%copy_to_store .or. copy_to_store
      if (present(caller)) then
         call caller%dependencies%add(node)
         call output_variable%dependent_nodes%add(caller)
      end if
      call variable_set%add_output_variable(output_variable)
   end subroutine add_call

end subroutine graph_add_variable

function source2operation(source) result(operation)
   integer,intent(in) :: source
   integer            :: operation
   select case (source)
   case (source_do, source_do_column, source_do_bottom, source_do_surface, source_do_horizontal)
      operation = source
   case (source_get_vertical_movement, source_initialize_state, source_check_state, source_get_light_extinction)
      operation = source_do
   case (source_initialize_bottom_state, source_check_bottom_state)
      operation = source_do_bottom
   case (source_initialize_surface_state, source_check_surface_state, source_get_drag, source_get_albedo)
      operation = source_do_surface
   case default
      call driver%fatal_error('source2operation', 'unknown source value')
   end select
end function source2operation

subroutine graph_connect(self, next)
   class (type_graph), intent(inout), target :: self
   class (type_graph), intent(inout), target :: next

   type (type_graph_set_member), pointer :: member

   allocate(member)
   member%p => next
   member%next => self%next%first
   self%next%first => member
   allocate(member)
   member%p => self
   member%next => next%previous%first
   next%previous%first => member
end subroutine graph_connect

subroutine node_set_add(self, node)
   class (type_node_set), intent(inout) :: self
   type (type_node), target             :: node

   type (type_node_set_member), pointer :: member

   ! First determine if the graph node is already part of the set. If so, we are done: return.
   member => self%first
   do while (associated(member))
      if (associated(member%p, node)) return
      member => member%next
   end do

   ! Graph node is not in set yet. Create pointer object and prepend it to the set.
   allocate(member)
   member%p => node
   member%next => self%first
   self%first => member
end subroutine node_set_add

function node_set_contains(self, node) result(found)
   class (type_node_set), intent(in) :: self
   type (type_node), target          :: node
   logical                           :: found

   type (type_node_set_member), pointer :: member

   ! Determine if the graph node is already part of the set.
   found = .true.
   member => self%first
   do while (associated(member))
      if (associated(member%p, node)) return
      member => member%next
   end do
   found = .false.
end function node_set_contains

subroutine node_set_finalize(self)
   class (type_node_set), intent(inout) :: self

   type (type_node_set_member), pointer :: current, next

   current => self%first
   do while (associated(current))
      next => current%next
      deallocate(current)
      current => next
   end do
   self%first => null()
end subroutine node_set_finalize

subroutine node_list_append(self, node)
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

subroutine node_list_check(self)
   class (type_node_list), intent(in) :: self
   type (type_node_list_member), pointer :: current, previous

   if (associated(self%last) .and. .not. associated(self%first)) call driver%fatal_error('node_list_check','BUG: list has tail but no head.')
   if (associated(self%first) .and. .not. associated(self%last)) call driver%fatal_error('node_list_check','BUG: list has head but no tail.')
   if (.not. associated(self%first)) return
   if (associated(self%first%previous)) call driver%fatal_error('node_list_check','BUG: head of list has pointer to previous node.')
   if (associated(self%last%next)) call driver%fatal_error('node_list_check','BUG: tail of list has pointer to next node.')
   previous => null()
   current => self%first
   do while (associated(current))
      if (associated(previous)) then
         if (.not. associated(current%previous, previous)) call driver%fatal_error('node_list_check','BUG: previous pointer does not match actual previous node.')
         if (.not. associated(previous%next, current)) call driver%fatal_error('node_list_check','BUG: next pointer does not match actual next node.')
      end if
      previous => current
      current => current%next
   end do
end subroutine node_list_check

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

subroutine node_list_remove(self, node)
   class (type_node_list), intent(inout) :: self
   type (type_node), target              :: node

   type (type_node_list_member), pointer :: pnode

   pnode => self%first
   do while (.not. associated(pnode%p, node))
      pnode => pnode%next
   end do
   _ASSERT_(associated(pnode), 'node_list_remove', 'Node ' // trim(node%as_string()) // ' not found in graph.')
   if (associated(pnode%next)) then
      pnode%next%previous => pnode%previous
   else
      self%last => pnode%previous
   end if
   if (associated(pnode%previous)) then
      pnode%previous%next => pnode%next
   else
      self%first => pnode%next
   end if
   deallocate(pnode)
end subroutine node_list_remove

subroutine node_list_finalize(self)
   class (type_node_list), intent(inout) :: self

   type (type_node_list_member), pointer :: current, next

   current => self%first
   do while (associated(current))
      next => current%next
      deallocate(current)
      current => next
   end do
   self%first => null()
   self%last => null()
end subroutine node_list_finalize

function node_list_find_node(self, model, source) result(pnode)
   class (type_node_list),  intent(in) :: self
   class (type_base_model),target,intent(in) :: model
   integer,                       intent(in) :: source
   type (type_node_list_member), pointer :: pnode

   pnode => self%first
   do while (associated(pnode))
      if (associated(pnode%p%model, model) .and. pnode%p%source == source) return
      pnode => pnode%next
   end do
end function node_list_find_node

function node_list_find(self, model, source) result(node)
   class (type_node_list),  intent(in) :: self
   class (type_base_model),target,intent(in) :: model
   integer,                       intent(in) :: source
   type (type_node), pointer :: node
   type (type_node_list_member), pointer :: pnode

   node => null()
   pnode => node_list_find_node(self, model, source)
   if (associated(pnode)) node => pnode%p
end function node_list_find

subroutine graph_finalize(self)
   class (type_graph),intent(inout) :: self

   type (type_node_list_member), pointer :: current

   current => self%first
   do while (associated(current))
      call current%p%inputs%finalize()
      call current%p%outputs%finalize(owner=.true.)
      call current%p%dependencies%finalize()
      deallocate(current%p)
      current => current%next
   end do
   call self%type_node_list%finalize()
end subroutine graph_finalize

subroutine graph_save_as_dot(self, unit)
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
         write (unit,'(A)') '  "' // trim(pnode%p%as_string()) // '" -> "' // trim(node%p%as_string()) // '";'
         pnode => pnode%next
      end do
      node => node%next
   end do

   write (unit,'(A)') '}'
end subroutine graph_save_as_dot

function node_as_string(node) result(string)
   class (type_node), intent(in) :: node
   character(len=attribute_length) :: string

   string = trim(node%model%get_path()) // ':' // trim(source2string(node%source))
end function node_as_string

end module fabm_graph

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

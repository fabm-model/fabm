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
      procedure :: remove    => node_list_remove
      procedure :: finalize  => node_list_finalize
      procedure :: find      => node_list_find
      procedure :: find_node => node_list_find_node
      procedure :: check     => node_list_check
   end type

   type type_call_stack_node
      class (type_base_model),       pointer :: model               => null()
      integer                                :: source              = source_unknown
      type (type_internal_variable), pointer :: requested_variable  => null()
      type (type_call_stack_node),   pointer :: previous            => null()
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
      procedure :: as_dot    => node_as_dot
   end type

   type type_graph_set_member
      class (type_graph), pointer :: p => null()
      type (type_graph_set_member), pointer :: next => null()
   end type

   type type_graph_set
      type (type_graph_set_member), pointer :: first => null()
   end type

   type, extends(type_node_list) :: type_graph
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

   function output_variable_set_add(self, variable) result(output_variable)
      class (type_output_variable_set), intent(inout) :: self
      type (type_internal_variable),target            :: variable

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

      type (type_output_variable_set_node), pointer :: node

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
      class (type_graph), pointer :: graph
      type (type_node),   pointer :: node

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
      class (type_graph), pointer :: self
      class (type_graph), pointer :: graph
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

   recursive function graph_add_call(self, model, source, stack_top) result(node)
      class (type_graph),      target, intent(inout) :: self
      class (type_base_model), target, intent(in)    :: model
      integer,                         intent(in)    :: source
      type (type_call_stack_node), target, optional  :: stack_top
      type (type_node), pointer :: node

      type (type_call_stack_node),  pointer :: existing_stack_node, stack_node
      character(len=2048)                   :: chain
      class (type_graph),           pointer :: root_graph, owner_graph, target_graph
      integer                               :: operation
      type (type_call_stack_node),  target  :: own_stack_node
      type (type_link),             pointer :: link
      type (type_input_variable),   pointer :: input_variable

      ! Circular dependency check:
      ! Search the stack, i.e., the list of calls that [indirectly] request the current call.
      ! If the current call is already on the stack, it is indirectly calling itself: there is a circular dependency.
      existing_stack_node => null()
      if (present(stack_top)) existing_stack_node => stack_top
      do while (associated(existing_stack_node))
         if (associated(existing_stack_node%model, model) .and. existing_stack_node%source == source) exit
         existing_stack_node => existing_stack_node%previous
      end do
      if (associated(existing_stack_node)) then
         ! Circular dependency found - report as fatal error.
         chain = trim(model%get_path()) // ':' // trim(source2string(source))
         stack_node => stack_top
         do
            chain = trim(stack_node%model%get_path()) // ':' // trim(source2string(stack_node%source)) &
               // ' needs ' // trim(stack_node%requested_variable%name) // ' provided by ' // trim(chain)
            if (associated(stack_node, existing_stack_node)) exit
            stack_node => stack_node%previous
         end do
         call driver%fatal_error('graph::add_call', 'circular dependency found: ' // trim(chain))
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

      ! Push this call onto the stack [the list of requesting calls]. This will be used for circular dependency checking.
      own_stack_node%model => model
      own_stack_node%source = source
      if (present(stack_top)) own_stack_node%previous => stack_top

      link => model%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. associated(link%original%read_index)) then
            ! This is the model's own variable (not inherited from child model) and the model itself originally requested read access to it.
            _ASSERT_(.not. associated(link%target%write_owner), 'graph::add_call', 'BUG: required input variable is co-written.')
            input_variable => node%inputs%add(link%target)
            if (resolve(link%target)) then
               own_stack_node%requested_variable => link%target
               call target_graph%add_variable(link%target, input_variable%sources, stack_top=own_stack_node, caller=node)
            end if
         end if
         link => link%next
      end do

      ! Add node to the graph.
      call target_graph%append(node)

   contains

      logical function resolve(variable)
         type (type_internal_variable), intent(in) :: variable
         resolve = .false.
         if (associated(variable%owner, model) .and. variable%source == source) return
         if (source == source_get_light_extinction .or. source == source_get_drag .or. source == source_get_albedo) return
         resolve = .true.
      end function

   end function graph_add_call

   recursive subroutine graph_add_variable(self, variable, variable_set, stack_top, copy_to_store, caller)
      class (type_graph),                 intent(inout) :: self
      type (type_internal_variable),      intent(in)    :: variable
      type (type_output_variable_set),    intent(inout) :: variable_set
      type (type_call_stack_node), target, optional     :: stack_top
      logical,          optional,         intent(in)    :: copy_to_store
      type (type_node), optional, target, intent(inout) :: caller

      type (type_variable_node), pointer :: variable_node

      if (self%frozen) call driver%fatal_error('graph_add_variable','Graph is frozen; no variables can be added.')

      if (associated(variable%cowriters)) then
         variable_node => variable%cowriters%first
         do while (associated(variable_node))
            call add_call(variable_node%target)
            variable_node => variable_node%next
         end do
      else
         call add_call(variable)
      end if

   contains

      recursive subroutine add_call(variable)
         type (type_internal_variable), intent(in) :: variable

         type (type_node),            pointer :: node
         type (type_output_variable), pointer :: output_variable

         if (variable%source == source_constant .or. variable%source == source_state .or. variable%source == source_external .or. variable%source == source_unknown) return
         _ASSERT_ (.not. variable%write_indices%is_empty(), 'graph_add_variable::add_call', 'Variable "' // trim(variable%name) // '" with source ' // trim(source2string(variable%source)) // ' does not have a write index')

         node => self%add_call(variable%owner, variable%source, stack_top)
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

   subroutine graph_save_as_dot(self, unit, subgraph)
      class (type_graph),intent(in) :: self
      integer,           intent(in) :: unit
      character(len=*), optional, intent(in) :: subgraph

      type (type_node_list_member),pointer :: node
      type (type_node_set_member), pointer :: pnode

      ! Nodes
      if (present(subgraph)) write (unit,'(A)') 'subgraph "cluster' // trim(subgraph) // '" {'
      write (unit,'(A)') '  label="' // trim(subgraph) // '";'
      write (unit,'(A)') '  style=filled;'
      write (unit,'(A)') '  node [color=black,style=filled];'
      node => self%first
      do while (associated(node))
         write (unit,'(A)') '  ' // trim(node%p%as_dot()) // ';'
         node => node%next
      end do
      if (present(subgraph))  write (unit,'(A)') '}'

      ! Edges
      node => self%first
      do while (associated(node))
         pnode => node%p%dependencies%first
         if (associated(pnode)) then
            write (unit,'(A)',advance='no') '  {'
            do while (associated(pnode))
               write (unit,'(A)',advance='no') ' "' // trim(pnode%p%as_string()) // '"'
               pnode => pnode%next
            end do
            write (unit,'(A)') ' } -> "' // trim(node%p%as_string()) // '";'
         end if
         node => node%next
      end do
   end subroutine graph_save_as_dot

   function node_as_dot(node) result(string)
      class (type_node), intent(in)   :: node
      character(len=attribute_length) :: string

      character(len=10)               :: color
      character(len=attribute_length) :: path

      select case (node%source)
         case (source_do);            color = 'white'
         case (source_do_bottom);     color = 'yellow'
         case (source_do_horizontal); color = 'green'
         case (source_do_surface);    color = 'blue'
         case (source_do_column);     color = 'red'
         case default;                color = 'grey'
      end select
      path = node%model%get_path()
      string = '"' // trim(node%as_string()) // '" [label="'//trim(path(2:))//'",fillcolor=' // trim(color) // ']'
   end function node_as_dot

   function node_as_string(node) result(string)
      class (type_node), intent(in) :: node
      character(len=attribute_length) :: string

      string = trim(node%model%get_path()) // ':' // trim(source2string(node%source))
   end function node_as_string

end module fabm_graph

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

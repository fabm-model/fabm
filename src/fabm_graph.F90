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

   public type_graph, type_node, type_node_set_member, type_node_list, type_node_list_member
   public type_output_variable

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
      integer                                :: prefill         = prefill_none
      type (type_node_set), pointer          :: group           => null()
      type (type_output_variable),   pointer :: next            => null()
   end type

   type type_output_variable_set
      type (type_output_variable), pointer :: first => null()
   contains
      procedure :: add      => output_variable_set_add
      procedure :: finalize => output_variable_set_finalize
   end type

   type type_node
      class (type_base_model), pointer :: model => null()
      integer                          :: source = source_unknown
      type (type_variable_set)         :: inputs               ! input variables (irrespective of their source - can be constants, state variables, host-provided, or diagnostics computed by another node)
      type (type_node_set)             :: dependencies         ! direct dependencies
      type (type_output_variable_set)  :: outputs              ! output variables that a later called model requires
   contains
      procedure :: as_string => node_as_string
   end type

   type,extends(type_node_list) :: type_graph
      type (type_graph),pointer :: previous => null()
      type (type_variable_set)  :: unresolved_dependencies
      logical                   :: frozen = .false.
   contains
      procedure :: add_call     => graph_add_call
      procedure :: add_variable => graph_add_variable
      procedure :: print        => graph_print
      procedure :: save_as_dot  => graph_save_as_dot
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
   node%prefill = variable%prefill
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
         if (variable%prefill==prefill_constant) write (*,'(",prefill=",g0.6)',advance='no') variable%target%prefill_value
         if (variable%prefill==prefill_previous_value) write (*,'(",prefill=previous")',advance='no')
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

recursive function graph_add_call(self,model,source,outer_calls,ignore_dependencies) result(node)
   class (type_graph),     target,intent(inout) :: self
   class (type_base_model),target,intent(in)    :: model
   integer,                       intent(in)    :: source
   type (type_node_list),  target,intent(inout) :: outer_calls
   logical,optional,              intent(in)    :: ignore_dependencies
   type (type_node), pointer :: node

   type (type_graph),            pointer :: current_graph
   type (type_node_list_member), pointer :: pnode
   character(len=2048)                   :: chain
   type (type_link),             pointer :: link
   logical                               :: ignore_dependencies_
   logical                               :: same_source

   if (self%frozen) call driver%fatal_error('graph_add_call','Graph is frozen; no calls can be added.')

   ! Provide optional arguments with default value.
   ignore_dependencies_ = .false.
   if (present(ignore_dependencies)) ignore_dependencies_ = ignore_dependencies

   ! For some APIs we never dynamically resolve dependencies
   ignore_dependencies_ = ignore_dependencies_ .or. source==source_get_light_extinction

   ! Check if this node is already in the graph, or in any of the preceding graphs.
   ! If it is, we are done: return.
   current_graph => self
   do while (associated(current_graph))
      node => current_graph%find(model,source)
      if (associated(node)) return
      current_graph => current_graph%previous
   end do

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
         if (associated(link%target%write_owner)) call driver%fatal_error('graph::add_call','BUG: required input variable is co-written.')
         call node%inputs%add(link%target)
         same_source = link%target%source==source .or. (link%target%source==source_unknown.and.(source==source_do_surface.or.source==source_do_bottom))
         if (.not.(associated(link%target%owner,model).and.same_source)) then
            if (ignore_dependencies_) then
               call self%unresolved_dependencies%add(link%target)
            else
               call self%add_variable(link%target,outer_calls,caller=node)
            end if
         end if
      end if
      link => link%next
   end do

   ! Remove node from the list of outer calls and add it to the graph instead.
   node => outer_calls%pop()
   call self%append(node)
end function graph_add_call

recursive subroutine graph_add_variable(self,variable,outer_calls,copy_to_store,caller,group)
   class (type_graph),              intent(inout) :: self
   type (type_internal_variable),   intent(in)    :: variable
   type (type_node_list),    target,intent(inout) :: outer_calls
   logical,         optional,       intent(in)    :: copy_to_store
   type (type_node),optional,target,intent(inout) :: caller
   type (type_node_set),optional,target,intent(inout) :: group

   type (type_variable_node), pointer :: variable_node
   type (type_node_set), pointer :: group_

   if (self%frozen) call driver%fatal_error('graph_add_variable','Graph is frozen; no variables can be added.')

   ! If this variable is not an output of some model (e.g., a state variable or external dependency), no call is needed.
   if (variable%write_indices%value==-1) return

   group_ => null()
   if (present(group)) group_ => group
   if (associated(variable%cowriters%first)) allocate(group_)

   if (variable%source==source_unknown) then
      ! This variable is either written by do_surface or do_bottom - which one of these two APIs is unknown.
      call add_call(source_do_surface)
      call add_call(source_do_bottom)
   elseif (variable%source/=source_none) then
      ! This variable is written by a known BGC API [is is not a constant indicated by source_none]
      call add_call(variable%source)
   end if

   ! Automatically request additional value contributions (for reduction operators that accept in-place modification of the variable value)
   variable_node => variable%cowriters%first
   do while (associated(variable_node))
      call self%add_variable(variable_node%target,outer_calls,copy_to_store,caller,group_)
      variable_node => variable_node%next
   end do

contains
   
   recursive subroutine add_call(source)
      integer, intent(in) :: source

      type (type_node),           pointer :: node
      type (type_output_variable),pointer :: variable_node

      node => self%add_call(variable%owner,source,outer_calls)
      variable_node => node%outputs%add(variable)
      variable_node%group => group_
      if (present(copy_to_store)) variable_node%copy_to_store = variable_node%copy_to_store .or. copy_to_store
      if (present(caller)) then
         call caller%dependencies%add(node)
         call variable_node%dependent_nodes%add(caller)
      end if
      if (associated(group_)) call group_%add(node)
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

   ! Determine if the graph node is already part of the set.
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

   current => self%first
   do while (associated(current))
      call current%p%inputs%finalize()
      call current%p%outputs%finalize()
      call current%p%dependencies%finalize()
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

   write (unit,'(A)') '}'
end subroutine graph_save_as_dot

function node_as_string(node) result(string)
   class (type_node), intent(in) :: node
   character(len=attribute_length) :: string

   string = trim(node%model%get_path())//':'//trim(source2string(node%source))
end function node_as_string

end module fabm_graph

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

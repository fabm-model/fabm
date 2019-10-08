#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_job: derived types that describe a single job (module subroutines to call and their dependencies)
!
! !INTERFACE:
module fabm_job

   use fabm_types
   use fabm_schedule
   use fabm_driver
   use fabm_graph

   implicit none

   private

   public type_job_manager, type_job, type_task, type_call
   public type_global_variable_register

   type type_graph_subset_node_pointer
      type (type_graph_subset_node),         pointer :: p    => null()
      type (type_graph_subset_node_pointer), pointer :: next => null()
   end type

   type type_graph_subset_node_set
      type (type_graph_subset_node_pointer), pointer :: first => null()
   contains
      procedure :: collect            => graph_subset_node_set_collect
      procedure :: collect_and_branch => graph_subset_node_set_collect_and_branch
      procedure :: branch             => graph_subset_node_set_branch
      procedure :: finalize           => graph_subset_node_set_finalize
   end type

   type type_graph_subset_node
      type (type_node), pointer         :: graph_node => null()
      integer                           :: outgoing_count = 0
      type (type_graph_subset_node_set) :: incoming
   end type

   type type_task_tree_node
      integer                             :: operation = source_unknown
      type (type_task_tree_node), pointer :: next_sibling => null()
      type (type_task_tree_node), pointer :: first_child  => null()
      type (type_task_tree_node), pointer :: parent       => null()
   contains
      procedure :: to_string         => task_tree_node_to_string
      procedure :: get_minimum_depth => task_tree_node_get_minimum_depth
      procedure :: get_leaf_at_depth => task_tree_node_get_leaf_at_depth
   end type

   type type_variable_request
      type (type_internal_variable),pointer :: variable => null()
      type (type_output_variable_set)       :: output_variable_set
      logical                               :: store    =  .false.
      type (type_variable_request),pointer  :: next     => null()
   end type

   type type_call_request
      class (type_base_model),pointer  :: model  => null()
      integer                          :: source = source_unknown
      type (type_call_request),pointer :: next   => null()
   end type

   type type_cache_copy_command
      integer :: read_index
      integer :: write_index
   end type

   ! A single call to a specific API of a specific biogeochemical model.
   ! It is defined by the combination of a model object ("model") and one of its procedures ("source").
   type type_call
      class (type_base_model), pointer            :: model => null()
      integer                                     :: source = source_unknown
      type (type_node), pointer                   :: graph_node => null()
      logical                                     :: active = .true.
      type (type_cache_copy_command), allocatable :: copy_commands_int(:) ! interior variables to copy from write to read cache after call completes
      type (type_cache_copy_command), allocatable :: copy_commands_hz(:)  ! horizontal variables to copy from write to read cache after call completes
      type (type_call), pointer                   :: next => null()
   end type type_call

   ! A task contains one or more model calls that all use the same operation over the domain.
   ! Valid operations: interior in native direction, interior per column, surface only, bottom only, horizontal-only.
   type type_task
      integer                   :: operation = source_unknown
      type (type_call), pointer :: first_call => null()
      type (type_variable_set)  :: cache_preload
      type (type_variable_set)  :: write_cache_preload

      integer, allocatable :: prefill(:)
      integer, allocatable :: prefill_hz(:)
      integer, allocatable :: save_sources(:)
      integer, allocatable :: save_sources_hz(:)
      integer, allocatable :: load(:)
      integer, allocatable :: load_hz(:)
      integer, allocatable :: load_scalar(:)
      type (type_task), pointer :: next => null()
   contains
      procedure :: initialize => task_initialize
      procedure :: finalize   => task_finalize
      procedure :: print      => task_print
   end type

   ! Job states (used for debugging call order)
   integer, parameter :: job_state_none = 0
   integer, parameter :: job_state_created = 1
   integer, parameter :: job_state_graph_created = 2
   integer, parameter :: job_state_tasks_created = 3
   integer, parameter :: job_state_finalized_prefill_settings = 4
   integer, parameter :: job_state_initialized = 5

   type type_job_node
      class (type_job), pointer :: p => null()
      type (type_job_node), pointer :: next => null()
   end type

   type type_job_set
      type (type_job_node), pointer :: first => null()
   end type

   ! A job contains one or more tasks, each using their own specific operation over the domain.
   type type_job
      character(len=attribute_length) :: name = ''
      integer                   :: state = job_state_none
      logical                   :: outsource_tasks     = .false.
      integer                   :: operation = source_unknown

      type (type_variable_request), pointer :: first_variable_request => null()
      type (type_call_request),     pointer :: first_call_request     => null()

      type (type_variable_set)  :: store_prefills
      logical, allocatable      :: interior_store_prefill(:)
      logical, allocatable      :: horizontal_store_prefill(:)

      type (type_graph)         :: graph
      type (type_task), pointer :: first_task => null()
      type (type_job_set)       :: previous
   contains
      procedure :: request_variable  => job_request_variable
      procedure :: request_call      => job_request_call
      procedure :: connect           => job_connect
      procedure :: print             => job_print
   end type

   type,extends(type_job_set) :: type_job_manager
   contains
      procedure :: create          => job_manager_create
      procedure :: initialize      => job_manager_initialize
      procedure :: print           => job_manager_print
   end type

   type type_variable_register
      type (type_variable_list) :: interior
      type (type_variable_list) :: horizontal
      type (type_variable_list) :: scalar
      logical :: frozen = .false.
   end type
   
   type type_global_variable_register
      type (type_variable_register) :: catalog
      type (type_variable_register) :: store
      type (type_variable_register) :: read_cache
      type (type_variable_register) :: write_cache
   contains
      procedure :: add_to_store => variable_register_add_to_store
      procedure :: add_to_catalog => variable_register_add_to_catalog
      procedure :: add_to_read_cache => variable_register_add_to_read_cache
      procedure :: add_to_write_cache => variable_register_add_to_write_cache
      procedure :: print => variable_register_print
   end type

   contains

   subroutine create_graph_subset_node_set(graph, set)
      type (type_graph),                 intent(in)  :: graph
      type (type_graph_subset_node_set), intent(out) :: set

      type (type_node_list_member),         pointer :: graph_node
      type (type_node_set_member),          pointer :: dependency
      type (type_graph_subset_node_pointer),pointer :: set_node, set_node2, set_node3

      ! Create representatives for each original graph node.
      graph_node => graph%first
      do while (associated(graph_node))
         allocate(set_node)
         allocate(set_node%p)
         set_node%p%graph_node => graph_node%p
         set_node%next => set%first
         set%first => set_node
         graph_node => graph_node%next
      end do

      set_node => set%first
      do while (associated(set_node))
         ! Enumerate dependencies of the original graph node
         dependency => set_node%p%graph_node%dependencies%first
         do while (associated(dependency))
            ! Find set node that corresponds to the dependency [set_node2]
            ! Note that this node does not exist if the dependency could be fulfilled by an earlier processed job/graph!
            set_node2 => set%first
            do while (associated(set_node2))
               if (associated(set_node2%p%graph_node, dependency%p)) exit
               set_node2 => set_node2%next
            end do

            if (associated(set_node2)) then
               ! Add dependency as "incoming" node.
               allocate(set_node3)
               set_node3%p => set_node2%p
               set_node3%next => set_node%p%incoming%first
               set_node%p%incoming%first => set_node3

               ! Increment the number of outgoing edges of the dependency itself.
               set_node2%p%outgoing_count = set_node2%p%outgoing_count + 1
            end if
            dependency => dependency%next
         end do
         set_node => set_node%next
      end do
   end subroutine create_graph_subset_node_set

   subroutine graph_subset_node_set_finalize(self)
      class (type_graph_subset_node_set), intent(inout) :: self

      type (type_graph_subset_node_pointer), pointer :: node, next_node, dependency, next_dependency

      node => self%first
      do while (associated(node))
         next_node => node%next
         dependency => node%p%incoming%first
         do while (associated(dependency))
            next_dependency => dependency%next
            deallocate(dependency)
            dependency => next_dependency
         end do
         deallocate(node%p)
         deallocate(node)
         node => next_node
      end do
      self%first => null()
   end subroutine graph_subset_node_set_finalize

   subroutine graph_subset_node_set_collect(self, operation, removed)
      class (type_graph_subset_node_set), intent(inout) :: self
      integer,                            intent(in)    :: operation
      class (type_graph_subset_node_set), intent(inout) :: removed

      type (type_graph_subset_node_pointer), pointer :: pnode, pnode2, pnode_previous, pnode_next
      logical                                        :: found, source_match
      integer                                        :: new_operation

      do
         found = .false.
         pnode_previous => null()
         pnode => self%first
         do while (associated(pnode))
            ! Store pointer to next node as we may remove the node from the set (thus destroying its original next pointer)
            pnode_next => pnode%next

            ! Determine whether the node's source matches
            new_operation = operation
            source_match = is_source_compatible(new_operation,pnode%p%graph_node%source)
            source_match = source_match .and. new_operation == operation ! Source is considered compatible only if it does not change the original task

            if (pnode%p%outgoing_count == 0 .and. source_match) then
               found = .true.

               ! Remove node from active set...
               if (associated(pnode_previous)) then
                  ! This is NOT the first node in the set.
                  pnode_previous%next => pnode%next
               else
                  ! This is the first node in the set.
                  self%first => pnode%next
               end if

               ! ...disconnect it from remaining nodes
               pnode2 => pnode%p%incoming%first
               do while (associated(pnode2))
                  pnode2%p%outgoing_count = pnode2%p%outgoing_count - 1
                  pnode2 => pnode2%next
               end do

               ! ...and add it to the set of removed nodes.
               pnode%next => removed%first
               removed%first => pnode
            else
               pnode_previous => pnode
            end if
            pnode => pnode_next
         end do
         if (.not.found) exit
      end do
   end subroutine graph_subset_node_set_collect

   recursive subroutine graph_subset_node_set_collect_and_branch(self, tree_node)
      class (type_graph_subset_node_set), intent(inout) :: self
      type (type_task_tree_node),target,  intent(inout) :: tree_node

      type (type_graph_subset_node_pointer), pointer :: pnode, pnode2, pnode_next
      type (type_graph_subset_node_set)              :: removed

      call self%collect(tree_node%operation, removed)

      if (.not.associated(self%first)) write (*,*) '  - ' // trim(tree_node%to_string())

      ! We have processed all graph end points with compatible sources.
      ! Now process end points with other sources.
      call self%branch(tree_node)

      ! Restore removed nodes: add each removed node back to the set, and increment edge count for all nodes that connect to it.
      pnode => removed%first
      do while (associated(pnode))
         pnode_next => pnode%next
         pnode2 => pnode%p%incoming%first
         do while (associated(pnode2))
            pnode2%p%outgoing_count = pnode2%p%outgoing_count + 1
            pnode2 => pnode2%next
         end do
         pnode%next => self%first
         self%first => pnode
         pnode => pnode_next
      end do
   end subroutine graph_subset_node_set_collect_and_branch

   function task_tree_node_to_string(self) result(string)
      class (type_task_tree_node), intent(in) :: self
      character(len=512) :: string

      type (type_task_tree_node), pointer :: node

      string = source2string(self%operation)
      node => self%parent
      do while (associated(node))
         if (node%operation /= source_unknown) string = trim(string) // ', ' // trim(source2string(node%operation))
         node => node%parent
      end do
   end function task_tree_node_to_string

   recursive function task_tree_node_get_minimum_depth(self) result(depth)
      class (type_task_tree_node),intent(in) :: self
      integer :: depth

      type (type_task_tree_node), pointer :: child

      if (.not.associated(self%first_child)) then
         depth = 0
      else
         depth = huge(depth)
         child => self%first_child
         do while (associated(child))
            depth = min(depth, child%get_minimum_depth())
            child => child%next_sibling
         end do
      end if
      depth = depth + 1
   end function task_tree_node_get_minimum_depth

   recursive function task_tree_node_get_leaf_at_depth(self, length) result(node)
      class (type_task_tree_node), target, intent(in) :: self
      integer,                             intent(in) :: length
      class (type_task_tree_node), pointer :: node

      type (type_task_tree_node), pointer :: child

      node => null()
      if (length == 1) then
         if (.not. associated(self%first_child)) node => self
         return
      end if

      child => self%first_child
      do while (associated(child))
         node => child%get_leaf_at_depth(length - 1)
         if (associated(node)) return
         child => child%next_sibling
      end do
   end function

   recursive subroutine graph_subset_node_set_branch(self, parent)
      class (type_graph_subset_node_set), intent(inout) :: self
      type (type_task_tree_node),target,  intent(inout) :: parent

      type (type_graph_subset_node_pointer), pointer :: pnode
      type (type_task_tree_node),            pointer :: pbranch
      integer                                        :: operation

      pnode => self%first
      do while (associated(pnode))
         if (pnode%p%outgoing_count==0) then
            operation = source2operation(pnode%p%graph_node%source)
            pbranch => parent%first_child
            do while (associated(pbranch))
               if (pbranch%operation == operation) exit
               pbranch => pbranch%next_sibling
            end do
            if (.not. associated(pbranch)) then
               allocate(pbranch)
               pbranch%operation = operation
               pbranch%parent => parent
               pbranch%next_sibling => parent%first_child
               parent%first_child => pbranch
            end if
         end if
         pnode => pnode%next
      end do

      pbranch => parent%first_child
      do while (associated(pbranch))
         call self%collect_and_branch(pbranch)
         pbranch => pbranch%next_sibling
      end do
   end subroutine graph_subset_node_set_branch

   subroutine task_initialize(self, variable_register, schedules)
      class (type_task),                    intent(inout) :: self
      type (type_global_variable_register), intent(inout) :: variable_register
      type (type_schedules),                intent(inout) :: schedules

      type (type_call),          pointer :: call_node
      type (type_variable_node), pointer :: input_variable

      ! Initialize individual call objects, then collect all input variables in a task-encompassing set.
      call_node => self%first_call
      do while (associated(call_node))
         input_variable => call_node%graph_node%inputs%first
         do while (associated(input_variable))
            _ASSERT_(.not. input_variable%target%read_indices%is_empty(), 'task_initialize', 'variable without read indices among inputs')
            _ASSERT_(.not. associated(input_variable%target%write_owner), 'task_initialize', 'write contribution among inputs')

            ! Make sure the variable has an entry in the read register.
            ! (the assigned index can be used in the global read store and in the read cache)
            call variable_register%add_to_read_cache(input_variable%target)

            ! Make sure this input is loaded into the read cache before the task is started.
            if (input_variable%target%source /= source_constant) call self%cache_preload%add(input_variable%target)

            input_variable => input_variable%next
         end do         
         call_node => call_node%next
      end do

      call_node => self%first_call
      do while (associated(call_node))
         call call_initialize(call_node, variable_register)
         call schedules%attach(call_node%model, call_node%source, call_node%active)
         call_node => call_node%next
      end do
   end subroutine task_initialize

   subroutine task_process_indices(self, unfulfilled_dependencies)
      class (type_task),       intent(inout) :: self
      type (type_variable_set),intent(inout) :: unfulfilled_dependencies

      type (type_call),           pointer :: call_node
      type (type_output_variable_set_node),pointer :: output_variable

      call_node => self%first_call
      do while (associated(call_node))
         call call_process_indices(call_node)
         call_node => call_node%next
      end do

      ! For all variables that this task computes itself, there is no need to preload a value in cache.
      call_node => self%first_call
      do while (associated(call_node))
         output_variable => call_node%graph_node%outputs%first
         do while (associated(output_variable))
            call self%cache_preload%remove(output_variable%p%target, discard=.true.)
            output_variable => output_variable%next
         end do
         call_node => call_node%next
      end do

      ! Find all variables that must be written to persistent storage after this job completes.
      call create_persistent_store_commands(self%save_sources,   domain_interior)
      call create_persistent_store_commands(self%save_sources_hz,domain_horizontal)

      ! Create prefill instructions for all variables that will be written to.
      call create_prefill_commands(self%prefill,   domain_interior)
      call create_prefill_commands(self%prefill_hz,domain_horizontal)

      ! Create read cache load instructions for all input variables.
      call create_load_commands(self%load,       domain_interior)
      call create_load_commands(self%load_hz,    domain_horizontal)
      call create_load_commands(self%load_scalar,domain_scalar)

   contains

      subroutine create_prefill_commands(prefill, domain)
         integer, intent(out),allocatable :: prefill(:)
         integer, intent(in)              :: domain

         integer                              :: ilast
         type (type_call),            pointer :: call_node
         type (type_output_variable_set_node), pointer :: output_variable
         type (type_variable_node),   pointer :: variable_node

         ! Find the last write cache index
         ilast = 0
         call_node => self%first_call
         do while (associated(call_node))
            output_variable => call_node%graph_node%outputs%first
            do while (associated(output_variable))
               if (output_variable%p%target%prefill /= prefill_none .and. iand(output_variable%p%target%domain, domain) /= 0) then
                  _ASSERT_(output_variable%p%target%write_indices%value > 0, 'create_prefill_commands', 'Variable ' // trim(output_variable%p%target%name) // ' was registered for prefilling, but it does not have a write cache index.')
                  _ASSERT_(output_variable%p%target%prefill /= prefill_previous_value .or. (output_variable%p%target%has_data .or. output_variable%p%target%store_index /= store_index_none), 'create_prefill_commands','Variable ' // trim(output_variable%p%target%name) // ' has prefill==previous value, but it does not have data.')
                  ilast = max(ilast, output_variable%p%target%write_indices%value)
               end if
               output_variable => output_variable%next
            end do
            call_node => call_node%next
         end do
         variable_node => self%write_cache_preload%first
         do while (associated(variable_node))
            if (iand(variable_node%target%domain, domain) /= 0) then
               _ASSERT_(variable_node%target%write_indices%value > 0, 'create_prefill_commands', 'Variable ' // trim(variable_node%target%name) // ' is set to be preloaded to write cache, but it does not have a write cache index.')
               _ASSERT_(variable_node%target%has_data .or. variable_node%target%store_index /= store_index_none, 'create_prefill_commands','Variable ' // trim(variable_node%target%name) // ' requires preloading to write cache, but it does not have data.')
               ilast = max(ilast, variable_node%target%write_indices%value)
            end if
            variable_node => variable_node%next
         end do

         allocate(prefill(ilast))
         prefill(:) = prefill_none

         if (ilast == 0) return

         call_node => self%first_call
         do while (associated(call_node))
            output_variable => call_node%graph_node%outputs%first
            do while (associated(output_variable))
               if (output_variable%p%target%prefill /= prefill_none .and. iand(output_variable%p%target%domain, domain) /= 0) then
                  ilast = output_variable%p%target%write_indices%value
                  if (output_variable%p%target%prefill == prefill_previous_value) then
                     if (associated(output_variable%p%target%write_owner)) then
                        prefill(ilast) = output_variable%p%target%write_owner%catalog_index
                     else
                        prefill(ilast) = output_variable%p%target%catalog_index
                     end if
                  else
                     prefill(ilast) = output_variable%p%target%prefill
                  end if
               end if
               output_variable => output_variable%next
            end do
            call_node => call_node%next
         end do
         variable_node => self%write_cache_preload%first
         do while (associated(variable_node))
            if (iand(variable_node%target%domain, domain) /= 0) then
               ilast = variable_node%target%write_indices%value
               if (associated(variable_node%target%write_owner)) then
                  prefill(ilast) = variable_node%target%write_owner%catalog_index
               else
                  prefill(ilast) = variable_node%target%catalog_index
               end if
            end if
            variable_node => variable_node%next
         end do
      end subroutine create_prefill_commands

      subroutine create_persistent_store_commands(commands, domain)
         integer,allocatable,intent(out) :: commands(:)
         integer,            intent(in)  :: domain

         integer                              :: ilast
         type (type_call),            pointer :: call_node
         type (type_output_variable_set_node), pointer :: variable_node

         ! First find the last index in persistent storage that will be written to.
         ilast = 0
         call_node => self%first_call
         do while (associated(call_node))
            variable_node => call_node%graph_node%outputs%first
            do while (associated(variable_node))
               if (variable_node%p%copy_to_store .and. iand(variable_node%p%target%domain, domain) /= 0 .and. call_node%graph_node%source /= source_constant) then
                  _ASSERT_(variable_node%p%target%write_indices%value > 0, 'create_prefill_commands', 'Variable ' // trim(variable_node%p%target%name) // ' has copy_to_store set, but it does not have a write cache index.')
                  _ASSERT_(variable_node%p%target%store_index /= store_index_none, 'create_persistent_store_commands', 'Variable ' // trim(variable_node%p%target%name) // ' has copy_to_store set, but it does not have a persistent storage index.')
                  ilast = max(ilast,variable_node%p%target%store_index)
               end if
               variable_node => variable_node%next
            end do
            call_node => call_node%next
         end do

         ! Allocate the commands array (go up to the last written-to index in persistent storage only)
         allocate(commands(ilast))
         commands(:) = 0

         ! Associate indices in persistent storage with the index in the write cache at which the source variable will be found.
         call_node => self%first_call
         do while (associated(call_node))
            variable_node => call_node%graph_node%outputs%first
            do while (associated(variable_node))
               if (variable_node%p%copy_to_store .and. iand(variable_node%p%target%domain,domain) /= 0 .and. call_node%graph_node%source /= source_constant) &
                  commands(variable_node%p%target%store_index) = variable_node%p%target%write_indices%value
               variable_node => variable_node%next
            end do
            call_node => call_node%next
         end do
      end subroutine create_persistent_store_commands

      subroutine create_load_commands(load, domain)
         integer,allocatable,intent(out) :: load(:)
         integer,            intent(in)  :: domain

         integer                           :: ilast
         type (type_variable_node),pointer :: input_variable

         ! First determine the maximum index to preload from/to. That determines the length of the array with preloading instructions.
         ilast = 0
         input_variable => self%cache_preload%first
         do while (associated(input_variable))
            _ASSERT_(input_variable%target%read_indices%value /= -1, 'create_load_commands', 'variable without valid read index among variables marked for preloading.')
            _ASSERT_(.not. input_variable%target%read_indices%is_empty(), 'create_load_commands', 'variable without read indices among variables marked for preloading.')
            if (iand(input_variable%target%domain, domain) /= 0) then
               if (input_variable%target%has_data .or. input_variable%target%store_index /= store_index_none) then
                  ! This variable has had data assigned in the global catalog (or it is to be included in the persistent store).
                  ! This dependency is thus fulfilled. Update the maximum index to load.
                  ilast = max(ilast,input_variable%target%read_indices%value)
               elseif (input_variable%target%presence /= presence_external_optional) then
                  ! This variable has NOT had data assigned in the global read store, and it is also no optional. This is an UNFULFILLED DEPENDENCY.
                  call unfulfilled_dependencies%add(input_variable%target)
               end if
            end if
            input_variable => input_variable%next
         end do

         ! Allocate array with preloading instructions, and initialize these to "do not preload"
         allocate(load(ilast))
         load(:) = 0

         ! Flag variables that require preloading
         input_variable => self%cache_preload%first
         do while (associated(input_variable))
            if (iand(input_variable%target%domain, domain) /= 0 .and. (input_variable%target%has_data .or. input_variable%target%store_index /= store_index_none)) &
               load(input_variable%target%read_indices%value) = input_variable%target%catalog_index
            input_variable => input_variable%next
         end do
      end subroutine create_load_commands

   end subroutine task_process_indices

   subroutine job_print(self, unit, specific_variable)
      class (type_job),                                intent(in) :: self
      integer,                                         intent(in) :: unit
      type (type_internal_variable), target, optional, intent(in) :: specific_variable

      type (type_task), pointer :: task

      write (unit,'(a,a)') 'Job: ',trim(self%name)
      task => self%first_task
      do while (associated(task))
         write (unit,'(a,a)') '- TASK WITH OPERATION = ',trim(source2string(task%operation))
         call task%print(unit, indent=2, specific_variable=specific_variable)
         task => task%next
      end do
   end subroutine job_print

subroutine call_finalize(self)
   type (type_call), intent(inout) :: self

   if (allocated(self%copy_commands_int)) deallocate(self%copy_commands_int)
   if (allocated(self%copy_commands_hz))  deallocate(self%copy_commands_hz)
   self%graph_node => null()
end subroutine call_finalize

subroutine print_output_variable(variable, unit, indent)
   type (type_output_variable), intent(in) :: variable
   integer,                     intent(in) :: unit
   integer, optional,           intent(in) :: indent

   integer                             :: read_index
   type (type_node_set_member),pointer :: pnode

   read_index = variable%target%read_indices%value
   if (associated(variable%target%write_owner)) read_index = variable%target%write_owner%read_indices%value

   write (unit,'(a,"   - ",a,", write@",i0)',advance='no') repeat(' ',indent), trim(variable%target%name), variable%target%write_indices%value
   if (variable%copy_to_cache) write (unit,'(", cache@",i0)',advance='no') read_index
   if (variable%copy_to_store) write (unit,'(", store@",i0)',advance='no') variable%target%store_index
   write (unit,*)
   pnode => variable%dependent_nodes%first
   do while (associated(pnode))
      if (associated(pnode%p%model)) then
         write (unit,'(a,"     <- ",a,": ",a)') repeat(' ',indent), trim(pnode%p%model%get_path()), trim(source2string(pnode%p%source))
      else
         write (unit,'(a,"     <- host")') repeat(' ',indent)
      end if
      pnode => pnode%next
   end do
end subroutine print_output_variable

subroutine print_input_variable(variable, unit, indent)
   type (type_variable_node), intent(in) :: variable
   integer,                   intent(in) :: unit
   integer,optional,          intent(in) :: indent

   write (unit,'(a,"   - ",a,", read@",i0)') repeat(' ',indent), trim(variable%target%name), variable%target%read_indices%value
end subroutine print_input_variable

subroutine task_print(self, unit, indent, specific_variable)
   class (type_task),                               intent(in) :: self
   integer,                                         intent(in) :: unit
   integer,optional,                                intent(in) :: indent
   type (type_internal_variable), target, optional, intent(in) :: specific_variable

   integer                             :: indent_
   integer                             :: i
   type (type_call),           pointer :: call_node
   type (type_output_variable_set_node),pointer :: output_variable
   logical                             :: show
   type (type_variable_node),  pointer :: input_variable
   logical                             :: header_written
   logical                             :: subheader_written

   indent_ = 0
   if (present(indent)) indent_ = indent

   if (size(self%load) > 0 .or. size(self%load_hz) > 0 .or. size(self%load_scalar) > 0) write (unit,'(a,a)') repeat(' ', indent_), 'Read cache prefilling:'
   do i = 1, size(self%load)
      if (self%load(i) /= 0) write (unit,'(a,"  - interior[",i0,"]")') repeat(' ', indent_), i
   end do
   do i = 1, size(self%load_hz)
      if (self%load_hz(i) /= 0) write (unit,'(a,"  - horizontal[",i0,"]")') repeat(' ', indent_), i
   end do
   do i = 1, size(self%load_scalar)
      if (self%load_scalar(i) /= 0) write (unit,'(a,"  - scalar[",i0,"]")') repeat(' ', indent_), i
   end do

   if (size(self%prefill) > 0 .or. size(self%prefill_hz) > 0) write (unit,'(a,a)') repeat(' ', indent_), 'Write cache prefilling:'
   do i = 1, size(self%prefill)
      select case (self%prefill(i))
      case (prefill_none)
      case (prefill_constant)
         write (unit,'(a,"  - interior[",i0,"] = constant value")') repeat(' ', indent_), i
      case default
         write (unit,'(a,"  - interior[",i0,"] = previous value")') repeat(' ', indent_), i
      end select
   end do
   do i = 1, size(self%prefill_hz)
      select case (self%prefill_hz(i))
      case (prefill_none)
      case (prefill_constant)
         write (unit,'(a,"  - horizontal[",i0,"] = constant value")') repeat(' ', indent_), i
      case default
         write (unit,'(a,"  - horizontal[",i0,"] = previous value")') repeat(' ', indent_), i
      end select
   end do

   call_node => self%first_call
   do while (associated(call_node))
      header_written = .false.
      if (.not. present(specific_variable)) call write_header()

      subheader_written = .false.
      input_variable => call_node%graph_node%inputs%first
      do while (associated(input_variable))
         show = .true.
         if (present(specific_variable)) show = associated(input_variable%target, specific_variable)
         if (show) then
            call write_header()
            if (.not. subheader_written) then
               write (unit,'(a,"   ",a)') repeat(' ', indent_), 'inputs:'
               subheader_written = .true.
            end if
            call print_input_variable(input_variable, unit, indent_)
         end if
         input_variable => input_variable%next
      end do

      subheader_written = .false.
      output_variable => call_node%graph_node%outputs%first
      do while (associated(output_variable))
         show = .true.
         if (present(specific_variable)) show = associated(output_variable%p%target, specific_variable)
         if (show) then
            call write_header()
            if (.not. subheader_written) then
               write (unit,'(a,"   ",a)') repeat(' ', indent_), 'outputs:'
               subheader_written = .true.
            end if
            call print_output_variable(output_variable%p, unit, indent_)
         end if
         output_variable => output_variable%next
      end do

      call_node => call_node%next
   end do

   contains

   subroutine write_header()
      if (.not. header_written) then
         if (associated(call_node%model)) then
            write (unit,'(a,a,": ",a)') repeat(' ', indent_), trim(call_node%model%get_path()), trim(source2string(call_node%source))
         else
            write (unit,'(a,"host")') repeat(' ', indent_)
         end if
         header_written = .true.
      end if
   end subroutine

end subroutine task_print

subroutine task_finalize(self)
   class (type_task), intent(inout) :: self

   type (type_call),pointer :: node,next

   node => self%first_call
   do while (associated(node))
      next => node%next
      call call_finalize(node)
      deallocate(node)
      node => next
   end do
   nullify(self%first_call)
end subroutine task_finalize

subroutine call_initialize(self, variable_register)
   type (type_call),             intent(inout) :: self
   type (type_global_variable_register),intent(inout) :: variable_register

   class (type_base_model),     pointer :: parent
   class (type_model_list_node),pointer :: model_list_node
   type (type_output_variable_set_node), pointer :: variable_node

   ! Make sure the pointer to the model has the highest class (and not a base class)
   ! This is needed because model classes that use inheritance and call base class methods
   ! may end up with model pointers that are base class specific (and do not reference
   ! procedures overwritten at a higher level)
   if (associated(self%model)) then
      parent => self%model%parent
      if (associated(parent)) then
         model_list_node => parent%children%first
         do while (associated(model_list_node))
            if (associated(self%model, model_list_node%model)) then
               ! Found ourselves in our parent - use the parent pointer to replace ours.
               self%model => model_list_node%model
               exit
            end if
            model_list_node => model_list_node%next
         end do
      end if
   end if

   ! For all output variables that other models are interested in, decide whether to copy their value
   ! from the write to read cache [if the other model will be called as part of the same task],
   ! of to save it to the persistent data store.
   variable_node => self%graph_node%outputs%first
   do while (associated(variable_node))
      call variable_register%add_to_write_cache(variable_node%p%target)
      if (variable_node%p%copy_to_cache) call variable_register%add_to_read_cache(variable_node%p%target)
      if (variable_node%p%copy_to_store) call variable_register%add_to_store(variable_node%p%target)
      variable_node => variable_node%next
   end do
end subroutine call_initialize

subroutine call_process_indices(self)
   type (type_call), intent(inout) :: self

   call create_cache_copy_commands(self%copy_commands_int, domain_interior)
   call create_cache_copy_commands(self%copy_commands_hz, domain_horizontal)

contains

   subroutine create_cache_copy_commands(commands, domain)
      type (type_cache_copy_command), allocatable, intent(out) :: commands(:)
      integer,                                     intent(in)  :: domain

      type (type_output_variable_set_node), pointer :: variable
      integer                              :: n, max_write_index, read_index, write_index

      n = 0
      max_write_index = -1
      variable => self%graph_node%outputs%first
      do while (associated(variable))
         if (variable%p%copy_to_cache .and. iand(variable%p%target%domain, domain) /= 0) then
            _ASSERT_(variable%p%target%write_indices%value > 0, 'call_list_node_initialize::create_cache_copy_commands', 'BUG: ' // trim(variable%p%target%name) // ' cannot be copied from write to read cache because it lacks a write cache index.')
            read_index = variable%p%target%read_indices%value
            if (associated(variable%p%target%write_owner)) read_index = variable%p%target%write_owner%read_indices%value
            _ASSERT_(read_index > 0, 'call_list_node_initialize::create_cache_copy_commands', 'BUG: ' // trim(variable%p%target%name) // ' cannot be copied from write to read cache because it lacks a read cache index.')
            n = n + 1
            max_write_index = max(max_write_index, variable%p%target%write_indices%value)
         end if
         variable => variable%next
      end do

      ! Create list of cache copy commands (source index in write cache, target index in read cache)
      ! They are sorted by write index (= where the data comes from) to hopefully benefit cache prefetching.
      allocate(commands(n))
      commands%read_index = -1
      commands%write_index = -1
      n = 0
      do write_index = 1, max_write_index
         variable => self%graph_node%outputs%first
         do while (associated(variable))
            if (variable%p%copy_to_cache .and. iand(variable%p%target%domain, domain) /= 0 .and. variable%p%target%write_indices%value == write_index) then
               n = n + 1
               read_index = variable%p%target%read_indices%value
               if (associated(variable%p%target%write_owner)) read_index = variable%p%target%write_owner%read_indices%value
               commands(n)%read_index = read_index
               commands(n)%write_index = write_index
               exit
            end if
            variable => variable%next
         end do
      end do
      _ASSERT_(all(commands%read_index > 0), 'call_process_indices', 'one or more read_index values <= 0')
      _ASSERT_(all(commands%write_index > 0), 'call_process_indices', 'one or more write_index values <= 0')
   end subroutine create_cache_copy_commands

end subroutine call_process_indices

subroutine job_request_variable(self, variable, store)
   class (type_job),target,      intent(inout)        :: self
   type (type_internal_variable),intent(inout),target :: variable
   logical, optional,            intent(in)           :: store

   type (type_variable_request), pointer :: variable_request

   _ASSERT_(self%state >= job_state_created, 'job_request_variable', 'Job has not been created yet.')
   _ASSERT_(self%state <= job_state_created, 'job_request_variable', 'Job "'//trim(self%name)//'" has already begun initialization; variables can no longer be requested.')

   ! Make sure this variable will not be merged (thus variable request must be filed before starting to merge variables!)
   variable%write_operator = ior(variable%write_operator, operator_merge_forbidden)

   _ASSERT_(.not. associated(variable%write_owner), 'job_request_variable','BUG: requested variable is co-written.')

   allocate(variable_request)
   variable_request%variable => variable
   if (present(store)) variable_request%store = store
   variable_request%next => self%first_variable_request
   self%first_variable_request => variable_request
end subroutine job_request_variable

subroutine job_request_call(self, model, source)
   class (type_job),target,intent(inout) :: self
   class (type_base_model),intent(in),target :: model
   integer,                intent(in)    :: source

   type (type_call_request), pointer :: call_request

   _ASSERT_(self%state >= job_state_created, 'job_request_call', 'Job has not been created yet.')
   _ASSERT_(self%state <= job_state_created, 'job_request_call', 'Job "'//trim(self%name)//'" has already begun initialization; calls can no longer be requested.')

   allocate(call_request)
   call_request%model => model
   call_request%source = source
   call_request%next => self%first_call_request
   self%first_call_request => call_request
end subroutine job_request_call

subroutine job_create_graph(self, variable_register)
   class (type_job),target,intent(inout) :: self
   type (type_global_variable_register), intent(inout) :: variable_register

   type (type_job_node),        pointer :: job_node
   type (type_node_list)                :: outer_calls
   type (type_variable_request),pointer :: variable_request
   type (type_call_request),    pointer :: call_request, next_call_request
   type (type_node),            pointer :: graph_node

   _ASSERT_(self%state >= job_state_created, 'job_create_graph', 'This job has not been created yet.')
   _ASSERT_(self%state < job_state_graph_created, 'job_create_graps', trim(self%name)//': graph for this job has already been created.')

   ! If we are linked to an earlier called job, make sure its graph has been created already.
   ! This is essential because we can skip calls if they appear already in the previous job - we determine this by exploring its graph.
   job_node => self%previous%first
   do while (associated(job_node))
      _ASSERT_(job_node%p%state >= job_state_graph_created, 'job_create_graph', trim(self%name)//': graph for previous job ('//trim(job_node%p%name)//') has not been created yet.')
      job_node => job_node%next
   end do

   ! Construct the dependency graph by adding explicitly requested variables and calls.
   ! We clean up [deallocate] the variable and call requests at the same time.
   variable_request => self%first_variable_request
   do while (associated(variable_request))
      call self%graph%add_variable(variable_request%variable, outer_calls, copy_to_store=variable_request%store, variable_set=variable_request%output_variable_set)
      if (variable_request%store) then
         ! FABM must be able to provide data for this variable across the entire spatial domain.
         ! If this variable is a constant, explicitly request a data field for it in the persistent store.
         ! If it is not a constant, the above call to add_variable will ensure that if the variable needs explicit computation,
         ! its value will also be copied to the persistent store.
         if (variable_request%variable%source == source_constant) call variable_register%add_to_store(variable_request%variable)
      else
         if (.not. associated(variable_request%output_variable_set%first)) call variable_register%add_to_write_cache(variable_request%variable)
      end if
      variable_request => variable_request%next
   end do

   call_request => self%first_call_request
   do while (associated(call_request))
      next_call_request => call_request%next
      graph_node => self%graph%add_call(call_request%model, call_request%source, outer_calls)
      deallocate(call_request)
      call_request => next_call_request
   end do
   self%first_call_request => null()
   call outer_calls%finalize()

   !self%graph%frozen = .true.

   ! Save graph
   !open(newunit=unit,file='graph.gv',action='write',status='replace',iostat=ios)
   !call self%graph%save_as_dot(unit)
   !close(unit)

   self%state = job_state_graph_created
end subroutine job_create_graph

subroutine job_create_tasks(self)
   class (type_job),target,intent(inout) :: self

   type (type_job_node),        pointer :: job_node
   type (type_graph_subset_node_set)    :: subset
   type (type_task_tree_node)           :: root
   class (type_task_tree_node), pointer :: leaf
   integer                              :: ntasks
   type (type_task),            pointer :: task
   type (type_call),            pointer :: call_node, next_call_node

   _ASSERT_(self%state < job_state_tasks_created, 'job_create_tasks', trim(self%name)//': tasks for this job have already been created.')
   _ASSERT_(self%state >= job_state_graph_created, 'job_create_tasks', trim(self%name)//': the graph for this job have not been created yet.')

   write (*,'(a)') trim(self%name)

   ! If we are linked to an earlier called job, make sure its task list has already been created.
   ! This is essential if we will try to outsource our own calls to previous tasks/jobs.
   job_node => self%previous%first
   do while (associated(job_node))
      _ASSERT_(job_node%p%state >= job_state_tasks_created, 'job_create_graph', trim(self%name) // ': tasks for previous job (' // trim(job_node%p%name) // ') have not been created yet.')
      job_node => job_node%next
   end do

   ! Create tree that describes all possible task orders (root is at the right of the call list, i.e., at the end of the very last call)
   write (*,'(a)') '- possible task orders:'
   call create_graph_subset_node_set(self%graph, subset)
   root%operation = self%operation
   if (root%operation /= source_unknown) then
      call subset%collect_and_branch(root)
   else
      call subset%branch(root)
   end if

   ! Determine optimum task order and use it to create the task objects.
   ntasks = root%get_minimum_depth() - 1
   write (*,'(a,i0,a)') '- best task order contains ',ntasks,' tasks.'
   leaf => root%get_leaf_at_depth(ntasks + 1)
   write (*,'(a,a)') '- best task order: ',trim(leaf%to_string())
   call create_tasks(leaf)
   _ASSERT_(.not. associated(subset%first), 'job_select_order', 'BUG: graph subset should be empty after create_tasks.')

   if (self%outsource_tasks) then
      task => self%first_task
      do while (associated(task))
         call_node => task%first_call
         do while (associated(call_node))
            next_call_node => call_node%next
            !call move_call_backwards(self,task,call_node)
            call_node => next_call_node
         end do
         task => task%next
      end do
   end if

   _ASSERT_(self%operation == source_unknown .or. .not. associated(self%first_task%next), 'job_select_order', 'Multiple tasks created while only one source was acceptable.')

   self%state = job_state_tasks_created

contains

   recursive subroutine create_tasks(tree_node)
      type (type_task_tree_node),target,intent(in) :: tree_node

      type (type_graph_subset_node_set)             :: removed
      type (type_graph_subset_node_pointer),pointer :: pnode
      type (type_task),                     pointer :: task
      type (type_call),                     pointer :: call_node, previous_call_node

      ! Start at the root of the tree. This represents the last task that will be executed;
      ! as tasks will be prepended to the list, it must be processed first.
      if (associated(tree_node%parent)) call create_tasks(tree_node%parent)

      ! If we have a root that does NOT represent a task itself, we are done.
      if (tree_node%operation==source_unknown) return

      ! Collect all nodes that can be processed with the currently selected source.
      call subset%collect(tree_node%operation,removed)

      ! Create the task and preprend it to the list.
      allocate(task)
      task%next => self%first_task
      self%first_task => task
      task%operation = tree_node%operation

      ! Collect all calls for this task.
      ! Preserve the order in which calls appear in the "removed" set,
      ! as this also represents the desired call order.
      previous_call_node => null()
      pnode => removed%first
      do while (associated(pnode))
         allocate(call_node)
         call_node%graph_node => pnode%p%graph_node
         call_node%model      => pnode%p%graph_node%model
         call_node%source     =  pnode%p%graph_node%source
         if (associated(previous_call_node)) then
            previous_call_node%next => call_node
         else
            task%first_call => call_node
         end if
         previous_call_node => call_node
         _ASSERT_(associated(call_node%model), 'create_tasks', 'Call node does not have a model pointer.')
         _ASSERT_(call_node%source /= source_constant .and. call_node%source /= source_state .and. call_node%source /= source_external .and. call_node%source /= source_unknown, 'create_tasks', 'Call node has invalid source.')
         pnode => pnode%next
      end do

      ! Clean-up array with processed calls.
      call removed%finalize()
   end subroutine create_tasks

   !subroutine move_call_backwards(job,task,call_node)
   !   class (type_job),target :: job
   !   type (type_task),target :: task
   !   type (type_call),target :: call_node
   !
   !   class (type_job),pointer :: current_job
   !   type (type_task),pointer :: current_task, target_task
   !   type (type_call),pointer :: current_call
   !   integer :: operation_after_merge
   !   logical :: compatible
   !
   !   write (*,*) 'moving '//trim(call_node%graph_node%as_string())
   !   current_job => job
   !   current_task => task
   !   target_task => null()
   !   compatible = .true.
   !   do while (move_one_step_backwards(current_job,current_task,call_node))
   !      operation_after_merge = current_task%operation
   !      compatible = is_source_compatible(operation_after_merge,call_node%source)
   !      compatible = compatible .and. operation_after_merge==current_task%operation
   !      if (compatible) then
   !         write (*,*) '  new task is compatible'
   !         target_task => current_task
   !      end if
   !   end do
   !   if (associated(target_task)) then
   !      ! Remove node from original task
   !      if (associated(task%first_call,call_node)) then
   !         task%first_call => call_node%next
   !      else
   !         current_call => task%first_call
   !         do while (.not.associated(current_call%next,call_node))
   !            current_call => current_call%next
   !         end do
   !         current_call%next => call_node%next
   !      end if
   !
   !      ! Append to target task
   !      current_call => target_task%first_call
   !      do while (associated(current_call%next))
   !         current_call => current_call%next
   !      end do
   !      current_call%next => call_node
   !      call_node%next => null()
   !
   !      !call target_task%initialize()
   !   end if
   !end subroutine move_call_backwards
   !
   !function move_one_step_backwards(job, task, travelling_call) result(moved)
   !   class (type_job),pointer :: job
   !   type (type_task),pointer :: task
   !   type (type_call),intent(in) :: travelling_call
   !   logical :: moved
   !
   !   type (type_call),pointer             :: call_node
   !   type (type_node_set_member), pointer :: dependency
   !
   !   moved = .false.
   !
   !   ! First determine if we can leave the current task
   !   ! (we cannot if it also handles one or more of our dependencies)
   !   call_node => task%first_call
   !   do while (associated(call_node))
   !      dependency => travelling_call%graph_node%dependencies%first
   !      do while (associated(dependency))
   !         if (associated(dependency%p,call_node%graph_node)) then
   !            write (*,*) '  cannot move past '//trim(call_node%graph_node%as_string())
   !            return
   !         end if
   !         dependency => dependency%next
   !      end do
   !      call_node => call_node%next
   !   end do
   !
   !   ! Move to previous task (if any)
   !   call get_previous_task(job,task)
   !   if (.not.associated(task)) return
   !
   !   moved = .true.
   !end function move_one_step_backwards

end subroutine job_create_tasks

function get_last_task(job) result(task)
   class (type_job),intent(in) :: job
   type (type_task),pointer :: task

   task => job%first_task
   if (.not. associated(task)) return
   do while (associated(task%next))
      task => task%next
   end do
end function get_last_task

recursive function find_responsible_task(job, output_variable) result(task)
   class (type_job), intent(in)        :: job
   type (type_output_variable), target :: output_variable
   type (type_task), pointer :: task
   type (type_job_node), pointer :: job_node

   task => job%first_task
   do while (associated(task))
      if (task_is_responsible(task, output_variable)) return
      task => task%next
   end do

   job_node => job%previous%first
   do while (associated(job_node) .and. .not. associated(task))
      task => find_responsible_task(job_node%p, output_variable)
      job_node => job_node%next
   end do
end function

logical function task_is_responsible(task, output_variable)
   class (type_task), intent(in)       :: task
   type (type_output_variable), target :: output_variable

   type (type_call),                     pointer :: call_node
   type (type_output_variable_set_node), pointer :: output_variable_node

   task_is_responsible = .true.
   call_node => task%first_call
   do while (associated(call_node))
      ! Loop over all outputs of this call
      output_variable_node => call_node%graph_node%outputs%first
      do while (associated(output_variable_node))
         if (associated(output_variable_node%p, output_variable)) return
         output_variable_node => output_variable_node%next
      end do
      call_node => call_node%next
   end do
   task_is_responsible = .false.
end function

subroutine job_finalize_prefill_settings(self)
   class (type_job), target, intent(inout) :: self

   type (type_task),                     pointer :: task, last_task, current_task
   class (type_job),                     pointer :: first_job
   type (type_variable_request),         pointer :: variable_request
   type (type_output_variable_set_node), pointer :: output_variable 
   logical                                       :: multiple_tasks

   _ASSERT_(self%state < job_state_finalized_prefill_settings, 'job_finalize_prefill_settings', 'This job has already been initialized.')
   _ASSERT_(self%state >= job_state_tasks_created, 'job_finalize_prefill_settings', 'Tasks for this job have not been created yet.')

   ! Set copy-to-cache and copy-to-store based on dependencies between different calls/tasks.
   task => self%first_task
   do while (associated(task))
      call prepare_task(task)
      task => task%next
   end do

   last_task => get_last_task(self)

   first_job => self
   do while (associated(first_job%previous%first))
      first_job => first_job%previous%first%p
   end do

   variable_request => self%first_variable_request
   do while (associated(variable_request))
      if (.not. variable_request%store) then
         ! This variable needs to end up in the write cache
         ! Any contributions from tasks other than the last need to be saved in the store and loaded into the write cache by the last task.
         ! If the variable is not written by anyone, it needs to be preloaded inot the write cache by the last task
         output_variable => variable_request%output_variable_set%first
         if (.not. associated(output_variable)) call last_task%write_cache_preload%add(variable_request%variable)
         do while (associated(output_variable))
            if (.not. task_is_responsible(last_task, output_variable%p)) then
               output_variable%p%copy_to_store = .true.
               call last_task%write_cache_preload%add(output_variable%p%target)
            end if
            output_variable => output_variable%next
         end do
      end if

      if (associated(variable_request%output_variable_set%first)) then
         if (associated(variable_request%output_variable_set%first%next)) then
            ! This request requires multiple variables to be written (e.g., all incrementing a sum)
            ! If these variables are all written by the same task, we're fine (the task will initialize the sum to 0 and all variables will be added)
            ! If they are written by multiple tasks, the very first one should handle the initialization to 0.
            ! But since it can be difficult to find out which task is first, we let the very first job initialize the corrresponding entry in the store.
            ! All tasks will prefill the relevant output variable with the value from the store, and all but the last will write the updated value to the store again.
            task => null()
            multiple_tasks = .false.
            output_variable => variable_request%output_variable_set%first
            do while (associated(output_variable))
               current_task => find_responsible_task(self, output_variable%p)
               _ASSERT_(associated(current_task), 'job_finalize_prefill_settings', 'Task responsible for ' // trim(variable_request%output_variable_set%first%p%target%name) // ' not found.')
               if (associated(task) .and. .not. associated(task, current_task)) multiple_tasks = .true.
               task => current_task
               output_variable => output_variable%next
            end do
            if (multiple_tasks) then
               output_variable => variable_request%output_variable_set%first
               do while (associated(output_variable))
                  current_task => find_responsible_task(self, output_variable%p)
                  call first_job%store_prefills%add(output_variable%p%target)
                  call current_task%write_cache_preload%add(output_variable%p%target)
                  if (.not. associated(task, last_task)) output_variable%p%copy_to_store = .true.
                  output_variable => output_variable%next
               end do
            end if
         end if
      end if

      variable_request => variable_request%next
   end do

   self%state = job_state_finalized_prefill_settings

contains

   subroutine prepare_task(task)
      type (type_task),intent(inout) :: task

      type (type_call),pointer :: call_node
      call_node => task%first_call
      do while (associated(call_node))
         call prepare_call(call_node)
         call_node => call_node%next
      end do
   end subroutine prepare_task

   subroutine prepare_call(call_node)
      type (type_call),intent(inout) :: call_node

      type (type_output_variable_set_node), pointer :: variable_node
      type (type_node_set_member), pointer :: dependent_call
      type (type_call),            pointer :: later_call

      ! For all output variables that other models are interested in, decide whether to copy their value
      ! from the write to read cache [if the other model will be called as part of the same task],
      ! or to save it to the persistent data store.
      variable_node => call_node%graph_node%outputs%first
      do while (associated(variable_node))
         ! For this output variable, loop over all calls that depend on it.
         dependent_call => variable_node%p%dependent_nodes%first
         do while (associated(dependent_call))
            ! For this dependent call, establish whether it is part of the same task, or of another task.
            later_call => call_node%next
            do while (associated(later_call))
               if (associated(later_call%graph_node, dependent_call%p)) exit
               later_call => later_call%next
            end do
            if (associated(later_call)) then
               ! The dependent call is part of the same task. Therefore the output needs to be copied to the read cache.
               ! Note for cowritten variables - we flag all contributors. The body of job_finalize_prefill_settings
               ! will make sure that only the final contributing one is actually cached.
               variable_node%p%copy_to_cache = .true.
            else
               ! The dependent call is part of another task. Therefore the output needs to be copied to the persistent store.
               variable_node%p%copy_to_store = .true.
            end if
            dependent_call => dependent_call%next
         end do
         variable_node => variable_node%next
      end do
   end subroutine prepare_call

end subroutine job_finalize_prefill_settings

subroutine job_initialize(self, variable_register, schedules)
   class (type_job),target,             intent(inout) :: self
   type (type_global_variable_register),intent(inout) :: variable_register
   type (type_schedules),               intent(inout) :: schedules

   type (type_task), pointer :: task

   _ASSERT_(self%state < job_state_initialized, 'job_initialize', trim(self%name)//': this job has already been initialized.')
   _ASSERT_(self%state >= job_state_finalized_prefill_settings, 'job_initialize', 'Prefill settings for this job have not been finalized yet.')

   ! Initialize tasks
   task => self%first_task
   do while (associated(task))
      call task%initialize(variable_register, schedules)
      task => task%next
   end do

   self%state = job_state_initialized

end subroutine job_initialize

subroutine job_process_indices(self, unfulfilled_dependencies)
   class (type_job),target, intent(inout) :: self
   type (type_variable_set),intent(inout) :: unfulfilled_dependencies

   type (type_task),          pointer :: task

   task => self%first_task
   do while (associated(task))
      call task_process_indices(task, unfulfilled_dependencies)
      task => task%next
   end do

   call gather_prefill(domain_interior, self%interior_store_prefill)
   call gather_prefill(domain_horizontal, self%horizontal_store_prefill)

contains

   subroutine gather_prefill(domain, flags)
      integer, intent(in) :: domain
      logical, allocatable :: flags(:)

      type (type_variable_node), pointer :: variable_node
      integer :: prefill_max

      prefill_max = 0
      variable_node => self%store_prefills%first
      do while (associated(variable_node))
         if (iand(variable_node%target%domain, domain) /= 0) prefill_max = max(prefill_max, variable_node%target%store_index)
         variable_node => variable_node%next
      end do
      allocate(flags(prefill_max))
      flags(:) = .false.
      variable_node => self%store_prefills%first
      do while (associated(variable_node))
         if (iand(variable_node%target%domain, domain) /= 0) flags(variable_node%target%store_index) = .true.
         variable_node => variable_node%next
      end do
   end subroutine

end subroutine job_process_indices

subroutine job_manager_create(self, job, name, source, outsource_tasks, previous)
   class (type_job_manager), intent(inout) :: self
   type (type_job),target,   intent(inout) :: job
   character(len=*),         intent(in)    :: name
   integer, optional,        intent(in)    :: source
   logical, optional,        intent(in)    :: outsource_tasks
   type (type_job), target, optional       :: previous

   type (type_job_node),pointer :: node

   _ASSERT_(job%state < job_state_created, 'job_manager_create','This job has already been created with name '//trim(job%name)//'.')
   job%state = job_state_created

   job%name = name
   if (present(source)) then
      job%operation = source2operation(source)
      job%graph%operation = job%operation
   end if
   if (present(outsource_tasks)) job%outsource_tasks = outsource_tasks

   allocate(node)
   node%p => job
   node%next => self%first
   self%first => node

   if (present(previous)) call previous%connect(job)
end subroutine job_manager_create

subroutine check_graph_duplicates(self)
   class (type_job_manager), intent(in) :: self

   type (type_job_node), pointer :: node
   type (type_node_list_member), pointer :: graph_node, graph_node2
   type (type_node_list) :: global_call_list

   node => self%first
   do while (associated(node))
      graph_node => node%p%graph%first
      do while (associated(graph_node))
         graph_node2 => global_call_list%find_node(graph_node%p%model, graph_node%p%source)
         _ASSERT_(.not. associated(graph_node2), 'job_manager_initialize', 'Call ' // trim(graph_node%p%as_string()) // ' appears multiple times in global graph.')
         call global_call_list%append(graph_node%p)
         graph_node => graph_node%next
      end do
      node => node%next
   end do
end subroutine check_graph_duplicates

subroutine job_manager_initialize(self, variable_register, schedules, unfulfilled_dependencies)
   class (type_job_manager),             intent(inout) :: self
   type (type_global_variable_register), intent(inout) :: variable_register
   type (type_schedules),                intent(inout) :: schedules
   type (type_variable_set),             intent(out)   :: unfulfilled_dependencies

   type (type_job_node), pointer :: node, first_ordered

   ! Order jobs according to call order.
   ! This ensures that jobs that are scheduled to run earlier are also initialized earlier.
   ! During initialization, a job can therefore expect any preceding jobs to have initialized completely.
   first_ordered => null()
   do while (associated(self%first))
      node => self%first
      self%first => self%first%next
      call add_to_order(node%p)
      deallocate(node)
   end do
   self%first => first_ordered

   ! Create all graphs. This must be done across all jobs before other operations that use graphs,
   ! since a job can add to graphs owned by other jobs.
   node => self%first
   do while (associated(node))
      call job_create_graph(node%p, variable_register)
      node => node%next
   end do

#ifndef NDEBUG
   ! Ensure each call appears exactly once in the superset of all graphs
   call check_graph_duplicates(self)
#endif

   ! Create tasks. This must be done for all jobs before job_finalize_prefill_settings is called, as this API operates across all jobs.
   node => self%first
   do while (associated(node))
      call job_create_tasks(node%p)
      node => node%next
   end do

   ! Finalize prefill settings (this has cross-job implications, so must be done for all jobs before they initialize (initialization uses prefill settings)
   node => self%first
   do while (associated(node))
      call job_finalize_prefill_settings(node%p)
      node => node%next
   end do

   ! Initialize all jobs. This assigns indices in the read and write caches, and in the persistent store.
   node => self%first
   do while (associated(node))
      call job_initialize(node%p, variable_register, schedules)
      node => node%next
   end do

   variable_register%read_cache%frozen = .true.
   variable_register%write_cache%frozen = .true.
   variable_register%store%frozen = .true.

   ! Create cache preload and copy instructions per task and call, and simultaneosuly check whether all dependencies are fulfilled.
   node => self%first
   do while (associated(node))
      call job_process_indices(node%p, unfulfilled_dependencies)
      node => node%next
   end do

contains

   recursive subroutine add_to_order(job)
      type (type_job), target :: job

      type (type_job_node), pointer :: node

      ! Make sure job is not yet in list
      node => first_ordered
      do while (associated(node))
         if (associated(node%p, job)) return
         node => node%next
      end do

      ! Append any jobs that run earlier first
      node => job%previous%first
      do while (associated(node))
         call add_to_order(node%p)
         node => node%next
      end do
      
      ! Append to list
      if (associated(first_ordered)) then
         node => first_ordered
         do while (associated(node%next))
            node => node%next
         end do
         allocate(node%next)
         node%next%p => job
      else
         allocate(first_ordered)
         first_ordered%p => job
      end if
   end subroutine add_to_order

end subroutine job_manager_initialize

subroutine job_manager_print(self, unit, specific_variable)
   class (type_job_manager),                        intent(in) :: self
   integer,                                         intent(in) :: unit
   type (type_internal_variable), target, optional, intent(in) :: specific_variable

   type (type_job_node),pointer :: node

   node => self%first
   do while (associated(node))
      call node%p%print(unit, specific_variable)
      node => node%next
   end do
end subroutine job_manager_print

subroutine job_connect(self, next)
   class (type_job), intent(inout), target :: self
   class (type_job), intent(inout), target :: next

   type (type_job_node), pointer :: node

   _ASSERT_(self%state <= job_state_created, 'job_connect','This job ('//trim(self%name)//') has already started initialization; it is too late to specify its place in the call order.')
   !_ASSERT_(.not. associated(self%previous), 'job_connect','This job ('//trim(self%name)//') has already been connected to a subsequent one.')

   allocate(node)
   node%p => self
   _ASSERT_(.not. associated(node%p, next), 'job_connect', 'Attempt to connect job ' // trim(self%name) // ' to itself.')
   node%next => next%previous%first
   next%previous%first => node
   call self%graph%connect(next%graph)
end subroutine job_connect

logical function is_source_compatible(operation, new_source)
   integer,intent(inout) :: operation
   integer,intent(in)    :: new_source

   integer :: new_operation

   new_operation = source2operation(new_source)

   if (operation==source_unknown) then
      ! Current operation is still undefined. Adopt the current and return.
      is_source_compatible = .true.
      operation = new_operation
      return
   end if

   select case (new_operation)
   case (source_do_bottom)
      is_source_compatible = operation == source_do_horizontal .or. operation == source_do_bottom
      if (operation == source_do_horizontal) operation = source_do_bottom
   case (source_do_surface)
      is_source_compatible = operation == source_do_horizontal .or. operation == source_do_surface
      if (operation == source_do_horizontal) operation = source_do_surface
   case (source_do_horizontal)
      is_source_compatible = operation == source_do_horizontal .or. operation == source_do_surface .or. operation == source_do_bottom
   case default
      is_source_compatible = operation == new_operation
   end select
end function is_source_compatible

function variable_register_add(self, variable, share_constants) result(index)
   type (type_variable_register), intent(inout) :: self
   type (type_internal_variable), target        :: variable
   logical,                       intent(in)    :: share_constants
   integer :: index

   _ASSERT_(.not. self%frozen, 'variable_register_add', 'Cannot add '//trim(variable%name)//'; register has been frozen.')
   select case (variable%domain)
   case (domain_interior)
      call add(self%interior)
   case (domain_horizontal, domain_surface, domain_bottom)
      call add(self%horizontal)
   case (domain_scalar)
      call add(self%scalar)
   end select

contains

   subroutine add(list)
      type (type_variable_list), intent(inout) :: list

      type (type_variable_node), pointer :: node

      if (share_constants .and. variable%source == source_constant) then
         ! This is a constant. See if there is already another constant with the same value in the register.
         ! If there is, reuse that entry instead of creating a new one.
         index = 0
         node => list%first
         do while (associated(node))
            index = index + 1
            if (node%target%source == source_constant .and. node%target%prefill_value == variable%prefill_value) return
            node => node%next
         end do
      end if
      call list%append(variable, index=index)
   end subroutine

end function variable_register_add

recursive subroutine variable_register_add_to_store(self, variable)
   class (type_global_variable_register),intent(inout) :: self
   type (type_internal_variable), target        :: variable

   type (type_variable_node), pointer :: variable_node

   ! If this variable has already been added to the persistent store, we are done: return.
   if (variable%store_index /= store_index_none) return

   ! If this variable is contributing to a variable (e.g., by adding to a sum), that other variable
   ! takes control. That controlling variable will then propagate its store index to all contributors.
   if (associated(variable%write_owner)) then
      call self%add_to_store(variable%write_owner)
      return
   end if

   ! Add the variable to the store and obtain its index.
   variable%store_index = variable_register_add(self%store, variable, share_constants=.true.)

   ! Propagate store index to any contributing variables.
   variable_node => variable%cowriters%first
   do while (associated(variable_node))
      variable_node%target%store_index = variable%store_index
      variable_node => variable_node%next
   end do
end subroutine variable_register_add_to_store

recursive subroutine variable_register_add_to_read_cache(self, variable)
   class (type_global_variable_register),intent(inout) :: self
   type (type_internal_variable), target        :: variable

   integer :: index

   ! If this variable has already been added to the read cache, we are done: return.
   if (variable%read_indices%value /= -1) return

   ! NB line below commented out because the variables that contribute together to a "reduce" operation (e.g., summation)
   ! may be called in any order. Only the last one may have copy_to_cache set, and that last one is not necessarily the write_owner.
   !_ASSERT_(.not. associated(variable%write_owner), 'variable_register_add_read', 'called on variable with owner')

   ! Add the variable to the register and obtain its index.
   index = variable_register_add(self%read_cache, variable, share_constants=.true.)

   ! Assign the read index to the variable.
   call variable%read_indices%set_value(index)

   !variable_node => variable%cowriters%first
   !do while (associated(variable_node))
   !   call variable_node%target%read_indices%set_value(i)
   !   variable_node => variable_node%next
   !end do
end subroutine variable_register_add_to_read_cache

subroutine variable_register_add_to_catalog(self, variable)
   class (type_global_variable_register),intent(inout) :: self
   type (type_internal_variable), target        :: variable

   if (variable%catalog_index /= -1) return
   variable%catalog_index = variable_register_add(self%catalog, variable, share_constants=.false.)
end subroutine variable_register_add_to_catalog

recursive subroutine variable_register_add_to_write_cache(self, variable)
   class (type_global_variable_register),intent(inout) :: self
   type (type_internal_variable), target        :: variable

   integer :: index

   ! If this variable has already been added to the write cache, we are done: return.
   if (variable%write_indices%value /= -1) return

   ! If this variable is contributing to a variable (e.g., by adding to a sum), that other variable
   ! takes control. That controlling variable will then propagate its write index to all contributors.
   if (associated(variable%write_owner)) then
      call self%add_to_write_cache(variable%write_owner)
      return
   end if

   ! Add the variable to the register and obtain its index.
   index = variable_register_add(self%write_cache, variable, share_constants=.true.)

   ! Assign the write index to the variable.
   ! This automatically propagates to other variables that contribute to this one,
   ! because their write indices are also referenced by the write owner.
   call variable%write_indices%set_value(index)
end subroutine variable_register_add_to_write_cache

subroutine variable_register_print(self, unit)
   class (type_global_variable_register), intent(in) :: self
   integer,                               intent(in) :: unit

   call print_list('Interior catalog:', self%catalog%interior)
   call print_list('Horizontal catalog:', self%catalog%horizontal)
   call print_list('Scalar catalog:', self%catalog%scalar)
   call print_list('Interior store:', self%store%interior)
   call print_list('Horizontal store:', self%store%horizontal)
   call print_list('Interior read cache:', self%read_cache%interior)
   call print_list('Horizontal read cache:', self%read_cache%horizontal)
   call print_list('Scalar read cache:', self%read_cache%scalar)
   call print_list('Interior write cache:', self%write_cache%interior)
   call print_list('Horizontal write cache:', self%write_cache%horizontal)

contains

   subroutine print_list(title, list)
      character(len=*),         intent(in) :: title
      type (type_variable_list),intent(in) :: list

      integer                            :: i
      type (type_variable_node), pointer :: variable_node

      write (unit,'(a)') title
      i = 0
      variable_node => list%first
      do while (associated(variable_node))
         i = i + 1
         write (unit,'("  ",i0,": ",a)') i, trim(variable_node%target%name)
         variable_node => variable_node%next
      end do
   end subroutine

end subroutine variable_register_print

end module fabm_job

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

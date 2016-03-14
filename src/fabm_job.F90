#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_job: derived types that describe a single job (module subroutines to call and their dependencies)
!
! !INTERFACE:
module fabm_job

   use fabm_types
   use fabm_driver
   use fabm_graph

   implicit none

   private

   public type_job, type_task, type_call
   public find_dependencies

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

   type type_cache_copy_command
      integer :: read_index
      integer :: write_index
   end type

   ! A call is defined by the combination of a model object and one of its procedures ("source").
   type type_call
      class (type_base_model), pointer            :: model => null()
      integer                                     :: source = source_unknown
      type (type_node), pointer                   :: node => null()
      type (type_cache_copy_command), allocatable :: copy_commands_int(:) ! interior variables to copy from write to read cache after call completes
      type (type_cache_copy_command), allocatable :: copy_commands_hz(:)  ! horizontal variables to copy from write to read cache after call completes
      type (type_call), pointer                   :: next => null()
   contains
      procedure :: initialize => call_initialize
      procedure :: finalize   => call_finalize
   end type type_call

   ! A task contains one or more model calls that all use the same operation over the domain.
   ! Valid operations: interior in native direction, interior per column, surface only, bottom only, horizontal-only.
   type type_task
      integer                   :: operation  = source_unknown
      type (type_call), pointer :: first_call => null()

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
      type (type_task), pointer :: next => null()
   contains
      procedure :: initialize => task_initialize
      procedure :: finalize   => task_finalize
      procedure :: print      => task_print
   end type

   ! A job contains one or more tasks, each using their own specific operation over the domain.
   type type_job
      type (type_task), pointer :: first_task => null()
      type (type_task), pointer :: final_task  => null()
      type (type_graph)         :: graph
      type (type_input_variable_set) :: required_cache_loads
      class (type_job), pointer :: previous => null()
   contains
      procedure :: initialize        => job_initialize
      procedure :: request_variable  => job_request_variable
      procedure :: request_variables => job_request_variables
      procedure :: request_call      => job_request_call
      procedure :: set_next          => job_set_next
      procedure :: print             => job_print
   end type

   contains

   subroutine create_graph_subset_node_set(graph,set)
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
               if (associated(set_node2%p%graph_node,dependency%p)) exit
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

      type (type_graph_subset_node_pointer), pointer :: node,next_node,dependency,next_dependency

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

   subroutine graph_subset_node_set_collect(self,operation,removed)
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
            source_match = source_match .and. new_operation==operation ! Source is considered compatible only if it does not change the original task

            if (pnode%p%outgoing_count==0 .and. source_match) then
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

   recursive subroutine graph_subset_node_set_collect_and_branch(self,tree_node)
      class (type_graph_subset_node_set), intent(inout) :: self
      type (type_task_tree_node),target,  intent(inout) :: tree_node

      type (type_graph_subset_node_pointer), pointer :: pnode, pnode2, pnode_next
      type (type_graph_subset_node_set)              :: removed

      call self%collect(tree_node%operation,removed)

      if (.not.associated(self%first)) write (*,*) trim(tree_node%to_string())

      ! We have processed all graph end points with compatible sources.
      ! Now process end points with other sources.
      call self%branch(tree_node)

      ! Restore removed nodes
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
      class (type_task_tree_node),intent(in) :: self
      character(len=attribute_length) :: string

      type (type_task_tree_node),pointer :: node

      write (string,'(i0)') self%operation
      node => self%parent
      do while (associated(node))
         if (node%operation/=source_unknown) write (string,'(a,",",i0)') trim(string),node%operation
         node => node%parent
      end do
   end function task_tree_node_to_string

   recursive function task_tree_node_get_minimum_depth(self) result(depth)
      class (type_task_tree_node),intent(in) :: self
      integer :: depth

      type (type_task_tree_node),pointer :: child

      if (.not.associated(self%first_child)) then
         depth = 0
      else
         depth = huge(depth)
         child => self%first_child
         do while (associated(child))
            depth = min(depth,child%get_minimum_depth())
            child => child%next_sibling
         end do
      end if
      depth = depth + 1
   end function task_tree_node_get_minimum_depth

   recursive function task_tree_node_get_leaf_at_depth(self,length) result(node)
      class (type_task_tree_node),target,intent(in) :: self
      integer,                           intent(in) :: length
      class (type_task_tree_node),pointer :: node

      type (type_task_tree_node),pointer :: child

      node => null()
      if (length==1) then
         if (.not.associated(self%first_child)) node => self
         return
      end if

      child => self%first_child
      do while (associated(child))
         node => child%get_leaf_at_depth(length-1)
         if (associated(node)) return
         child => child%next_sibling
      end do
   end function

   recursive subroutine graph_subset_node_set_branch(self,parent)
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
               if (pbranch%operation==operation) exit
               pbranch => pbranch%next_sibling
            end do
            if (.not.associated(pbranch)) then
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

   subroutine task_initialize(self,extra_inputs)
      class (type_task),                      intent(inout) :: self
      type (type_input_variable_set),optional,intent(in)    :: extra_inputs

      type (type_call),          pointer :: call_node
      type (type_input_variable_set)     :: inputs
      type (type_input_variable),pointer :: input_variable

      if (present(extra_inputs)) then
         input_variable => extra_inputs%first
         do while (associated(input_variable))
            call inputs%add(input_variable%target)
            input_variable => input_variable%next
         end do
      end if

      ! Initialize individual call objects, then collect all input variables in a task-encompassing set.
      call_node => self%first_call
      do while (associated(call_node))
         call call_node%initialize()
         input_variable => call_node%node%inputs%first
         do while (associated(input_variable))
            call inputs%add(input_variable%target)
            input_variable => input_variable%next
         end do
         call_node => call_node%next
      end do

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

      ! Deallocate task-specific set of inputs - all necessary information has been incorporated in self%load*
      call inputs%finalize()

   contains

      subroutine create_prefill_commands(prefill,indices,values,domain)
         integer, intent(out),allocatable :: prefill(:)
         integer, intent(out),allocatable :: indices(:)
         real(rk),intent(out),allocatable :: values(:)
         integer, intent(in)              :: domain

         integer                              :: ilast
         type (type_call),            pointer :: call_node
         type (type_output_variable), pointer :: variable_node

         ! Find the last write cache index
         ilast = 0
         call_node => self%first_call
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

         call_node => self%first_call
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
         integer,allocatable,intent(out) :: commands(:)
         integer,            intent(in)  :: domain

         integer                              :: ilast
         type (type_call),            pointer :: call_node
         type (type_output_variable), pointer :: variable_node

         ! First find the last index in persistent storage that will be written to.
         ilast = 0
         call_node => self%first_call
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
         call_node => self%first_call
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
         logical,allocatable,intent(out) :: load(:)
         integer,            intent(in)  :: domain

         integer                            :: ilast
         type (type_input_variable),pointer :: input_variable

         ilast = 0
         input_variable => inputs%first
         do while (associated(input_variable))
            if (iand(input_variable%target%domain,domain)/=0 .and. input_variable%target%read_indices%value/=-1) &
               ilast = max(ilast,input_variable%target%read_indices%value)
            input_variable => input_variable%next
         end do

         allocate(load(ilast))
         load(:) = .false.

         input_variable => inputs%first
         do while (associated(input_variable))
            if (iand(input_variable%target%domain,domain)/=0 .and. input_variable%target%read_indices%value/=-1) &
               load(input_variable%target%read_indices%value) = .true.
            input_variable => input_variable%next
         end do
      end subroutine create_load_commands

   end subroutine task_initialize

   subroutine job_print(self)
      class (type_job), intent(in) :: self

      type (type_task), pointer :: task

      task => self%first_task
      do while (associated(task))
         write (*,'(a,a)') 'TASK WITH OPERATION = ',trim(source2string(task%operation))
         call task%print()
         task => task%next
      end do
      if (associated(self%final_task)) then
         write (*,'(a,a)') 'FINAL TASK WITH OPERATION = ',trim(source2string(self%final_task%operation))
         call self%final_task%print()
      end if
   end subroutine job_print

subroutine call_finalize(self)
   class (type_call), intent(inout) :: self

   if (allocated(self%copy_commands_int)) deallocate(self%copy_commands_int)
   if (allocated(self%copy_commands_hz))  deallocate(self%copy_commands_hz)
   self%node => null()
end subroutine call_finalize

subroutine task_print(self)
   class (type_task), intent(in) :: self

   type (type_call),           pointer :: call_node
   type (type_output_variable),pointer :: variable
   type (type_node_set_member),pointer :: pnode

   call_node => self%first_call
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
end subroutine task_print

subroutine task_finalize(self)
   class (type_task), intent(inout) :: self

   type (type_call),pointer :: node,next

   node => self%first_call
   do while (associated(node))
      next => node%next
      call node%finalize()
      deallocate(node)
      node => next
   end do
   nullify(self%first_call)
end subroutine task_finalize

subroutine call_initialize(self)
   class (type_call),intent(inout) :: self

   class (type_base_model),     pointer :: parent
   class (type_model_list_node),pointer :: model_list_node
   type (type_output_variable), pointer :: variable_node
   type (type_node_set_member), pointer :: dependent_call
   type (type_call),            pointer :: later_call

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

   ! For all output variables that other models are interested in, decide whether to copy their value
   ! from the write to read cache [if the other model will be called as part of the same task],
   ! of to save it to the persistent data store.
   variable_node => self%node%outputs%first
   do while (associated(variable_node))
      ! For this output variable, loop over all calls that depend on it.
      dependent_call => variable_node%dependent_nodes%first
      do while (associated(dependent_call))
         ! For this dependent call, establish whether it is part of the same task, or of another task.
         later_call => self%next
         do while (associated(later_call))
            if (associated(later_call%node,dependent_call%p)) exit
            later_call => later_call%next
         end do
         if (associated(later_call)) then
            ! The dependent call is part of the same task. Therefore the output needs to be copied to the read cache.
            variable_node%copy_to_cache = .true.
         else
            ! The dependent call is part of another task. Therefore the output needs to be copied to the persistent store.
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
      type (type_cache_copy_command), allocatable, intent(out) :: commands(:)
      integer,                                     intent(in)  :: domain

      type (type_output_variable), pointer :: variable
      integer                              :: iorder, n, maxorder

      n = 0
      maxorder = -1
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
            maxorder = max(maxorder,variable%target%read_indices%value)
         end if
         variable => variable%next
      end do

      ! Create list of cache copy commands (source index in write cache, target index in read cache)
      ! They are sorted by write index to make (write) access to memory more predictable.
      allocate(commands(n))
      n = 0
      do iorder=1,maxorder
         variable => self%node%outputs%first
         do while (associated(variable))
            if (variable%copy_to_cache.and.iand(variable%target%domain,domain)/=0.and.variable%target%read_indices%value==iorder) then
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

subroutine job_request_variable(self,variable,copy_to_cache,copy_to_store,not_stale)
   class (type_job),target,      intent(inout) :: self
   type (type_internal_variable),intent(in)    :: variable
   logical,optional,             intent(in)    :: copy_to_cache
   logical,optional,             intent(in)    :: copy_to_store
   logical,optional,             intent(in)    :: not_stale

   type (type_node_list) :: outer_calls

   ! If this variable is not written [but a field provided by the host or a state variable], return immediately.
   if (variable%write_indices%is_empty().or.variable%source==source_none) then
      if (present(copy_to_cache)) then
         if (copy_to_cache) call self%required_cache_loads%add(variable)
      end if
      return
   end if

   call self%graph%add_variable(variable,outer_calls,copy_to_cache,copy_to_store,not_stale)

   call outer_calls%finalize()
end subroutine job_request_variable

subroutine job_request_variables(self,link_list,copy_to_cache,copy_to_store,not_stale)
   class (type_job),     intent(inout) :: self
   type (type_link_list),intent(in)    :: link_list
   logical,optional,     intent(in)    :: copy_to_cache
   logical,optional,     intent(in)    :: copy_to_store
   logical,optional,     intent(in)    :: not_stale

   type (type_link), pointer :: link

   link => link_list%first
   do while (associated(link))
      call self%request_variable(link%target,copy_to_cache,copy_to_store,not_stale)
      link => link%next
   end do
end subroutine job_request_variables

subroutine job_request_call(self,model,source,ignore_dependencies)
   class (type_job),target,intent(inout) :: self
   class (type_base_model),intent(in)    :: model
   integer,                intent(in)    :: source
   logical,optional,       intent(in)    :: ignore_dependencies

   type (type_node_list)    :: outer_calls
   type (type_node),pointer :: graph_node

   graph_node => self%graph%add_call(model,source,outer_calls,ignore_dependencies=ignore_dependencies)

   call outer_calls%finalize()
end subroutine job_request_call

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

subroutine job_initialize(self,final_operation,outsource_tasks)
   class (type_job),target,intent(inout) :: self
   integer,optional,       intent(in)    :: final_operation
   logical,optional,       intent(in)    :: outsource_tasks

   type (type_graph_subset_node_set)   :: subset
   type (type_task_tree_node)          :: root
   class (type_task_tree_node),pointer :: leaf
   integer                             :: ntasks
   type (type_task),           pointer :: task
   integer                             :: unit,ios
   type (type_call),           pointer :: call_node, next_call_node, new_call_node
   logical                             :: outsource_tasks_

   outsource_tasks_ = .false.
   if (present(outsource_tasks)) outsource_tasks_ = outsource_tasks

   ! Save graph
   !open(newunit=unit,file='graph.gv',action='write',status='replace',iostat=ios)
   !call self%graph%save_as_dot(unit)
   !close(unit)

   ! Create tree that describes all possible task orders (root is at the right of the call list, i.e., at the end of the very last call)
   call create_graph_subset_node_set(self%graph,subset)
   if (present(final_operation)) then
      root%operation = final_operation
      call subset%collect_and_branch(root)
   else
      call subset%branch(root)
   end if

   ! Determine optimum task order and use it to create the task objects.
   ntasks = root%get_minimum_depth()-1
   write (*,'(a,i0,a)') 'Best task order contains ',ntasks,' tasks.'
   leaf => root%get_leaf_at_depth(ntasks+1)
   write (*,'(a,a)') 'Best task order: ',trim(leaf%to_string())
   call create_tasks(leaf)
   if (associated(subset%first)) call driver%fatal_error('job_select_order','BUG: graph subset should be empty after create_tasks.')

   if (outsource_tasks_) then
      task => self%first_task
      do while (associated(task))
         call_node => task%first_call
         do while (associated(call_node))
            next_call_node => call_node%next
            call move_call_backwards(self,task,call_node)
            call_node => next_call_node
         end do
         task => task%next
      end do
   end if

   if (present(final_operation)) then
      ! Separate the last task (remove from main task list)
      task => null()
      self%final_task => self%first_task
      do while (associated(self%final_task%next))
         task => self%final_task
         self%final_task => self%final_task%next
      end do
      if (associated(task)) then
         task%next => null()
      else
         self%first_task => null()
      end if

      ! For any calls that MUST be processed as part of the last task, but currently are allocated to one of the earlier tasks,
      ! prepend a separate call to the final task. Note that this implies the call will be made twice!
      task => self%first_task
      do while (associated(task))
         call_node => task%first_call
         do while (associated(call_node))
            if (call_node%node%not_stale.and.is_source_compatible(self%final_task%operation,call_node%node%source)) then
               allocate(new_call_node)
               new_call_node%node   => call_node%node
               new_call_node%model  => call_node%model
               new_call_node%source =  call_node%source
               new_call_node%next => self%final_task%first_call
               self%final_task%first_call => new_call_node
            end if
            call_node => call_node%next
         end do
         task => task%next
      end do
   end if

   ! Initialize tasks
   task => self%first_task
   do while (associated(task))
      call task%initialize()
      task => task%next
   end do
   if (associated(self%final_task)) call self%final_task%initialize(extra_inputs=self%required_cache_loads)

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
         call_node%node   => pnode%p%graph_node
         call_node%model  => pnode%p%graph_node%model
         call_node%source =  pnode%p%graph_node%source
         if (associated(previous_call_node)) then
            previous_call_node%next => call_node
         else
            task%first_call => call_node
         end if
         previous_call_node => call_node
         pnode => pnode%next
      end do

      ! Clean-up array with processed calls.
      call removed%finalize()
   end subroutine create_tasks

   subroutine move_call_backwards(job,task,call_node)
      class (type_job),target :: job
      type (type_task),target :: task
      type (type_call),target :: call_node

      class (type_job),pointer :: current_job
      type (type_task),pointer :: current_task, target_task
      type (type_call),pointer :: current_call
      integer :: operation_after_merge
      logical :: compatible

      write (*,*) 'moving '//trim(call_node%node%as_string())
      current_job => job
      current_task => task
      target_task => null()
      compatible = .true.
      do while (move_one_step_backwards(current_job,current_task,call_node))
         operation_after_merge = current_task%operation
         compatible = is_source_compatible(operation_after_merge,call_node%source)
         compatible = compatible .and. operation_after_merge==current_task%operation
         if (compatible) then
            write (*,*) '  new task is compatible'
            target_task => current_task
         end if
      end do
      if (associated(target_task)) then
         ! Remove node from original task
         if (associated(task%first_call,call_node)) then
            task%first_call => call_node%next
         else
            current_call => task%first_call
            do while (.not.associated(current_call%next,call_node))
               current_call => current_call%next
            end do
            current_call%next => call_node%next
         end if

         ! Append to target task
         current_call => target_task%first_call
         do while (associated(current_call%next))
            current_call => current_call%next
         end do
         current_call%next => call_node
         call_node%next => null()

         call target_task%initialize()
      end if
   end subroutine move_call_backwards

   function move_one_step_backwards(job,task,travelling_call) result(moved)
      class (type_job),pointer :: job
      type (type_task),pointer :: task
      type (type_call),intent(in) :: travelling_call
      logical :: moved

      type (type_call),pointer             :: call_node
      type (type_node_set_member), pointer :: dependency
      type (type_task),pointer             :: old_task

      moved = .false.

      ! First determine if we can leave the current task
      ! (we cannot if it also handles one or more of our dependencies)
      call_node => task%first_call
      do while (associated(call_node))
         dependency => travelling_call%node%dependencies%first
         do while (associated(dependency))
            if (associated(dependency%p,call_node%node)) then
               write (*,*) '  cannot move past '//trim(call_node%node%as_string())
               return
            end if
            dependency => dependency%next
         end do
         call_node => call_node%next
      end do

      if (associated(job%first_task,task).or.(.not.associated(job%first_task).and.associated(job%final_task,task))) then
         ! We are the first task in the job - move to previous job (and its last task)
         write (*,*) '  moving to previous job'
         job => job%previous
         do while (associated(job))
            if (associated(job%first_task).or.associated(job%final_task)) exit
            job => job%previous
         end do
         if (.not.associated(job)) return
         task => job%final_task
         if (.not.associated(task)) then
            task => job%first_task
            do while (associated(task%next))
               task => task%next
            end do
         end if
      else
         ! We are NOT the first task in the job - find the previous task
         write (*,*) '  moving to previous task'
         if (associated(job%final_task,task)) then
            ! We are the final task - move to last task in the list of preparatory tasks
            task => job%first_task
            do while (associated(task%next))
               task => task%next
            end do
         else
            ! We are one of the preparatory tasks - move to the preceding one.
            old_task => task
            task => job%first_task
            do while (.not.associated(task%next,old_task))
               task => task%next
            end do
         end if
      end if
      moved = .true.
   end function move_one_step_backwards

end subroutine job_initialize

subroutine job_set_next(self,next)
   class (type_job), intent(inout), target :: self
   type (type_job),  intent(inout), target :: next

   if (associated(next%graph%previous)) &
      call driver%fatal_error('job_set_next','This job has already been connected to a subsequent one.')

   next%previous => self
   next%graph%previous => self%graph
end subroutine job_set_next

function source2operation(source) result(operation)
   integer,intent(in) :: source
   integer            :: operation
   select case (source)
   case (source_do,source_do_column,source_do_bottom,source_do_surface,source_do_horizontal)
      operation = source
   case (source_get_vertical_movement,source_initialize_state,source_check_state,source_get_light_extinction)
      operation = source_do
   case (source_initialize_bottom_state,source_check_bottom_state)
      operation = source_do_bottom
   case (source_initialize_surface_state,source_check_surface_state,source_get_drag,source_get_albedo)
      operation = source_do_surface
   case default
      call driver%fatal_error('source2operation','unknown source value')
   end select
end function source2operation

logical function is_source_compatible(operation,new_source)
   integer,intent(inout) :: operation
   integer,intent(in)    :: new_source

   integer :: new_operation

   new_operation = source2operation(new_source)

   if (operation==source_unknown) then
      is_source_compatible = .true.
      operation = new_operation
      return
   end if

   select case (new_operation)
   case (source_do_bottom)
      is_source_compatible = operation==source_do_horizontal.or.operation==source_do_bottom
      if (operation==source_do_horizontal) operation = source_do_bottom
   case (source_do_surface)
      is_source_compatible = operation==source_do_horizontal.or.operation==source_do_surface
      if (operation==source_do_horizontal) operation = source_do_surface
   case (source_do_horizontal)
      is_source_compatible = operation==source_do_horizontal.or.operation==source_do_surface.or.operation==source_do_bottom
   case default
      is_source_compatible = operation==new_operation
   end select
end function is_source_compatible

end module fabm_job

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

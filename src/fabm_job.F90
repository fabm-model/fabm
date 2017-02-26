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

   public type_job_manager, type_job, type_task, type_call
   public type_variable_register

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
      type (type_internal_variable),pointer :: variable      => null()
      logical                               :: requirements = 0
      logical                               :: copy_to_store = .false.
      type (type_variable_request),pointer  :: next          => null()
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

   ! A call is defined by the combination of a model object and one of its procedures ("source").
   type type_call
      class (type_base_model), pointer            :: model => null()
      integer                                     :: source = source_unknown
      type (type_node), pointer                   :: graph_node => null()
      type (type_cache_copy_command), allocatable :: copy_commands_int(:) ! interior variables to copy from write to read cache after call completes
      type (type_cache_copy_command), allocatable :: copy_commands_hz(:)  ! horizontal variables to copy from write to read cache after call completes
      type (type_call), pointer                   :: next => null()
   contains
      procedure :: initialize      => call_initialize
      procedure :: process_indices => call_process_indices
      procedure :: finalize        => call_finalize
   end type type_call

   ! A task contains one or more model calls that all use the same operation over the domain.
   ! Valid operations: interior in native direction, interior per column, surface only, bottom only, horizontal-only.
   type type_task
      integer                   :: operation  = source_unknown
      type (type_call), pointer :: first_call => null()
      type (type_variable_set)  :: cache_preload

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

   integer,parameter :: job_state_none = 0, job_state_created = 1, job_state_graph_created = 2, job_state_tasks_created = 3, job_state_finalized_prefill_settings = 4, job_state_initialized = 5

   ! A job contains one or more tasks, each using their own specific operation over the domain.
   type type_job
      character(len=attribute_length) :: name = ''
      integer                   :: state = job_state_none
      logical                   :: ignore_dependencies = .false.
      logical                   :: outsource_tasks     = .false.
      type (type_node), pointer :: own_node => null()

      type (type_task), pointer :: first_task => null()
      type (type_task), pointer :: final_task => null()
      type (type_graph)         :: graph
      type (type_variable_set) :: required_in_write_cache
      type (type_variable_request), pointer :: first_variable_request => null()
      type (type_call_request),     pointer :: first_call_request     => null()
      class (type_job), pointer :: previous => null()
      class (type_job), pointer :: dependency_handler => null()
   contains
      procedure :: request_variable  => job_request_variable
      procedure :: request_call      => job_request_call
      procedure :: set_previous      => job_set_previous
      procedure :: print             => job_print
   end type

   type type_job_manager_item
      type (type_job),             pointer :: job  => null()
      type (type_job_manager_item),pointer :: next => null()
   end type

   type type_job_manager
      type (type_job_manager_item),pointer :: first => null()
      class (type_job),            pointer :: default_dependency_handler => null()
   contains
      procedure :: create          => job_manager_create
      procedure :: initialize      => job_manager_initialize
      procedure :: process_indices => job_manager_process_indices
      procedure :: print           => job_manager_print
   end type

   type type_variable_register
      type (type_variable_list) :: interior_store
      type (type_variable_list) :: horizontal_store
      type (type_variable_list) :: interior_read
      type (type_variable_list) :: horizontal_read
      type (type_variable_list) :: scalar_read
   contains
      procedure :: add_store => variable_register_add_store
      procedure :: add_read  => variable_register_add_read
      procedure :: print     => variable_register_print
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

   subroutine task_initialize(self,variable_register)
      class (type_task),                intent(inout) :: self
      type (type_variable_register),    intent(inout) :: variable_register

      type (type_call),         pointer :: call_node
      type (type_variable_node),pointer :: input_variable

      ! Initialize individual call objects, then collect all input variables in a task-encompassing set.
      call_node => self%first_call
      do while (associated(call_node))
         input_variable => call_node%graph_node%inputs%first
         do while (associated(input_variable))
            if (input_variable%target%read_indices%is_empty()) call fatal_error('task_initialize','BUG: variable without read indices among inputs')
            if (associated(input_variable%target%write_owner)) call fatal_error('task_initialize','BUG: write contribution among inputs')
            call self%cache_preload%add(input_variable%target)
            call variable_register%add_read(input_variable%target)
            if (input_variable%target%source==source_none.and..not.input_variable%target%write_indices%is_empty()) call variable_register%add_store(input_variable%target)
            input_variable => input_variable%next
         end do
         call_node => call_node%next
      end do

      call_node => self%first_call
      do while (associated(call_node))
         call call_node%initialize(variable_register)
         call_node => call_node%next
      end do
   end subroutine task_initialize

   subroutine task_process_indices(self,unfulfilled_dependencies)
      class (type_task),       intent(inout) :: self
      type (type_variable_set),intent(inout) :: unfulfilled_dependencies

      type (type_call),           pointer :: call_node
      type (type_output_variable),pointer :: output_variable

      call_node => self%first_call
      do while (associated(call_node))
         call call_node%process_indices()
         call_node => call_node%next
      end do

      ! For all variables that this task computes itself, there is no need to preload a value in cache.
      call_node => self%first_call
      do while (associated(call_node))
         output_variable => call_node%graph_node%outputs%first
         do while (associated(output_variable))
            call self%cache_preload%remove(output_variable%target,discard=.true.)
            output_variable => output_variable%next
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
            variable_node => call_node%graph_node%outputs%first
            do while (associated(variable_node))
               if (variable_node%prefill/=prefill_none.and.iand(variable_node%target%domain,domain)/=0) then
                  if (variable_node%target%write_indices%value<=0) &
                     call driver%fatal_error('create_prefill_commands','Variable '//trim(variable_node%target%name) &
                        //' has prefilling set, but it does not have a write cache index.')
                  if (variable_node%prefill==prefill_previous_value.and.variable_node%target%store_index==store_index_none) &
                     call driver%fatal_error('create_prefill_commands','Variable '//trim(variable_node%target%name) &
                        //' has prefill==previous value, but it does not have a persistent storage index.')
                  ilast = max(ilast,variable_node%target%write_indices%value)
               end if
               variable_node => variable_node%next
            end do
            call_node => call_node%next
         end do

         allocate(prefill(ilast))
         allocate(indices(ilast))
         allocate(values(ilast))
         prefill(:) = prefill_none
         indices(:) = 0

         if (ilast==0) return

         call_node => self%first_call
         do while (associated(call_node))
            variable_node => call_node%graph_node%outputs%first
            do while (associated(variable_node))
               if (variable_node%prefill/=prefill_none.and.iand(variable_node%target%domain,domain)/=0) then
                  ilast = variable_node%target%write_indices%value
                  prefill(ilast) = variable_node%prefill
                  indices(ilast) = variable_node%target%store_index
                  values (ilast) = variable_node%target%prefill_value
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
            variable_node => call_node%graph_node%outputs%first
            do while (associated(variable_node))
               if (variable_node%copy_to_store.and.iand(variable_node%target%domain,domain)/=0) then
                  if (variable_node%target%write_indices%value<=0) &
                     call driver%fatal_error('create_prefill_commands','Variable '//trim(variable_node%target%name) &
                        //' has copy_to_store set, but it does not have a write cache index.')
                  if (variable_node%target%store_index==store_index_none) &
                     call driver%fatal_error('create_persistent_store_commands','Variable '//trim(variable_node%target%name) &
                        //' has copy_to_store set, but it does not have a persistent storage index.')
                  ilast = max(ilast,variable_node%target%store_index)
               end if
               variable_node => variable_node%next
            end do
            call_node => call_node%next
         end do

         ! Allocate the commands array (go up to the last written-to index in persistent storage only)
         allocate(commands(ilast))
         commands(:) = -1

         ! Associate indices in persistent storage with the index in the write cache at which the source variable will be found.
         call_node => self%first_call
         do while (associated(call_node))
            variable_node => call_node%graph_node%outputs%first
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

         integer                           :: ilast
         type (type_variable_node),pointer :: input_variable

         ilast = 0
         input_variable => self%cache_preload%first
         do while (associated(input_variable))
            if (input_variable%target%read_indices%value==-1) call fatal_error('create_load_commands','BUG: variable without valid read index among variables marked for preloading.')
            if (input_variable%target%read_indices%is_empty()) call fatal_error('create_load_commands','BUG: variable without read indices among variables marked for preloading.')
            if (iand(input_variable%target%domain,domain)/=0) then
               if (input_variable%target%in_read_registry) then
                  ilast = max(ilast,input_variable%target%read_indices%value)
               elseif (input_variable%target%presence/=presence_external_optional) then
                  call unfulfilled_dependencies%add(input_variable%target)
               end if
            end if
            input_variable => input_variable%next
         end do

         allocate(load(ilast))
         load(:) = .false.

         input_variable => self%cache_preload%first
         do while (associated(input_variable))
            if (iand(input_variable%target%domain,domain)/=0 .and. input_variable%target%in_read_registry) &
               load(input_variable%target%read_indices%value) = .true.
            input_variable => input_variable%next
         end do
      end subroutine create_load_commands

   end subroutine task_process_indices

   subroutine job_print(self, specific_variable)
      class (type_job),                                intent(in) :: self
      type (type_internal_variable), target, optional, intent(in) :: specific_variable

      type (type_task), pointer :: task

      write (*,'(a,a)') 'Job: ',trim(self%name)
      task => self%first_task
      do while (associated(task))
         write (*,'(a,a)') '- TASK WITH OPERATION = ',trim(source2string(task%operation))
         call task%print(indent=2,specific_variable=specific_variable)
         task => task%next
      end do
      if (associated(self%final_task)) then
         write (*,'(a,a)') '- FINAL TASK WITH OPERATION = ',trim(source2string(self%final_task%operation))
         call self%final_task%print(indent=2,specific_variable=specific_variable)
      end if
   end subroutine job_print

subroutine call_finalize(self)
   class (type_call), intent(inout) :: self

   if (allocated(self%copy_commands_int)) deallocate(self%copy_commands_int)
   if (allocated(self%copy_commands_hz))  deallocate(self%copy_commands_hz)
   self%graph_node => null()
end subroutine call_finalize

subroutine print_output_variable(variable,indent)
   type (type_output_variable),intent(in) :: variable
   integer,optional,           intent(in) :: indent

   integer                             :: read_index
   type (type_node_set_member),pointer :: pnode

   read_index = variable%target%read_indices%value
   if (associated(variable%target%write_owner)) read_index = variable%target%write_owner%read_indices%value

   write (*,'(a,"   - ",a,", write@",i0)',advance='no') repeat(' ',indent), trim(variable%target%name), variable%target%write_indices%value
   if (variable%copy_to_cache) write (*,'(", cache@",i0)',advance='no') read_index
   if (variable%copy_to_store) write (*,'(", store@",i0)',advance='no') variable%target%store_index
   if (variable%prefill==prefill_constant) write (*,'(", prefill=",g0.6)',advance='no') variable%target%prefill_value
   if (variable%prefill==prefill_previous_value) write (*,'(", prefill=previous")',advance='no')
   write (*,*)
   pnode => variable%dependent_nodes%first
   do while (associated(pnode))
      if (associated(pnode%p%model)) then
         write (*,'(a,"     <- ",a,": ",a)') repeat(' ',indent), trim(pnode%p%model%get_path()), trim(source2string(pnode%p%source))
      else
         write (*,'(a,"     <- host")') repeat(' ',indent)
      end if
      pnode => pnode%next
   end do
end subroutine print_output_variable

subroutine print_input_variable(variable,indent)
   type (type_variable_node),intent(in) :: variable
   integer,optional,         intent(in) :: indent

   write (*,'(a,"   - ",a,", read@",i0)') repeat(' ',indent), trim(variable%target%name), variable%target%read_indices%value
end subroutine print_input_variable

subroutine task_print(self,indent,specific_variable)
   class (type_task),                               intent(in) :: self
   integer,optional,                                intent(in) :: indent
   type (type_internal_variable), target, optional, intent(in) :: specific_variable

   integer                             :: indent_
   integer                             :: i
   type (type_call),           pointer :: call_node
   type (type_output_variable),pointer :: output_variable
   logical                             :: show
   type (type_variable_node),  pointer :: input_variable
   logical                             :: header_written
   logical                             :: subheader_written

   indent_ = 0
   if (present(indent)) indent_ = indent

   if (allocated(self%prefill_type)) then
      do i=1,size(self%prefill_type)
         select case (self%prefill_type(i))
         case (prefill_constant)
            write (*,'(a,"prefilling write[",i0,"] = ",g0.6)') repeat(' ',indent_), i, self%prefill_values(i)
         case (prefill_previous_value)
            write (*,'(a,"prefilling write[",i0,"] = previous value")') repeat(' ',indent_), i
         end select
      end do
      do i=1,size(self%prefill_type_hz)
         select case (self%prefill_type_hz(i))
         case (prefill_constant)
            write (*,'(a,"prefilling write_hz[",i0,"] = ",g0.6)') repeat(' ',indent_), i, self%prefill_values_hz(i)
         case (prefill_previous_value)
            write (*,'(a,"prefilling write_hz[",i0,"] = previous value")') repeat(' ',indent_), i
         end select
      end do
   end if

   call_node => self%first_call
   do while (associated(call_node))
      header_written = .false.
      if (.not.present(specific_variable)) call write_header()

      subheader_written = .false.
      input_variable => call_node%graph_node%inputs%first
      do while (associated(input_variable))
         show = .true.
         if (present(specific_variable)) show = associated(input_variable%target, specific_variable)
         if (show) then
            call write_header()
            if (.not. subheader_written) then
               write (*,'(a,"   ",a)') repeat(' ',indent_),'inputs:'
               subheader_written = .true.
            end if
            call print_input_variable(input_variable,indent_)
         end if
         input_variable => input_variable%next
      end do

      subheader_written = .false.
      output_variable => call_node%graph_node%outputs%first
      do while (associated(output_variable))
         show = .true.
         if (present(specific_variable)) show = associated(output_variable%target, specific_variable)
         if (show) then
            call write_header()
            if (.not. subheader_written) then
               write (*,'(a,"   ",a)') repeat(' ',indent_),'outputs:'
               subheader_written = .true.
            end if
            call print_output_variable(output_variable,indent_)
         end if
         output_variable => output_variable%next
      end do

      call_node => call_node%next
   end do

   contains

   subroutine write_header()
      if (.not.header_written) then
         if (associated(call_node%model)) then
            write (*,'(a,a,": ",a)') repeat(' ',indent_),trim(call_node%model%get_path()),trim(source2string(call_node%source))
         else
            write (*,'(a,"host")') repeat(' ',indent_)
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
      call node%finalize()
      deallocate(node)
      node => next
   end do
   nullify(self%first_call)
end subroutine task_finalize

subroutine call_initialize(self,variable_register)
   class (type_call),            intent(inout) :: self
   type (type_variable_register),intent(inout) :: variable_register

   class (type_base_model),     pointer :: parent
   class (type_model_list_node),pointer :: model_list_node
   type (type_output_variable), pointer :: variable_node

   ! Make sure the pointer to the model has the highest class (and not a base class)
   ! This is needed because model classes that use inheritance and call base class methods
   ! may end up with model pointers that are base class specific (and do not reference
   ! procedures overwritten at a higher level)
   if (associated(self%model)) then
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
   end if

   ! For all output variables that other models are interested in, decide whether to copy their value
   ! from the write to read cache [if the other model will be called as part of the same task],
   ! of to save it to the persistent data store.
   variable_node => self%graph_node%outputs%first
   do while (associated(variable_node))
      if (variable_node%copy_to_cache) call variable_register%add_read(variable_node%target)
      if (variable_node%copy_to_store) call variable_register%add_store(variable_node%target)
      variable_node => variable_node%next
   end do
end subroutine call_initialize

subroutine call_process_indices(self)
   class (type_call),intent(inout) :: self

   call create_cache_copy_commands(self%copy_commands_int, domain_interior)
   call create_cache_copy_commands(self%copy_commands_hz,  domain_horizontal)

contains

   subroutine create_cache_copy_commands(commands,domain)
      type (type_cache_copy_command), allocatable, intent(out) :: commands(:)
      integer,                                     intent(in)  :: domain

      type (type_output_variable), pointer :: variable
      integer                              :: iorder, n, maxorder, read_index

      n = 0
      maxorder = -1
      variable => self%graph_node%outputs%first
      do while (associated(variable))
         if (variable%copy_to_cache.and.iand(variable%target%domain,domain)/=0) then
            if (variable%target%write_indices%value<=0) &
               call driver%fatal_error('call_list_node_initialize::create_cache_copy_commands', &
                  'BUG: '//trim(variable%target%name)//' cannot be copied from write to read cache because it lacks a write cache index.')
            read_index = variable%target%read_indices%value
            if (associated(variable%target%write_owner)) read_index = variable%target%write_owner%read_indices%value
            if (read_index<=0) &
               call driver%fatal_error('call_list_node_initialize::create_cache_copy_commands', &
                  'BUG: '//trim(variable%target%name)//' cannot be copied from write to read cache because it lacks a read cache index.')
            n = n + 1
            maxorder = max(maxorder,read_index)
         end if
         variable => variable%next
      end do

      ! Create list of cache copy commands (source index in write cache, target index in read cache)
      ! They are sorted by write index to make (write) access to memory more predictable.
      allocate(commands(n))
      n = 0
      do iorder=1,maxorder
         variable => self%graph_node%outputs%first
         do while (associated(variable))
            read_index = variable%target%read_indices%value
            if (associated(variable%target%write_owner)) read_index = variable%target%write_owner%read_indices%value
            if (variable%copy_to_cache.and.iand(variable%target%domain,domain)/=0.and.read_index==iorder) then
               n = n + 1
               commands(n)%read_index = read_index
               commands(n)%write_index = variable%target%write_indices%value
               exit
            end if
            variable => variable%next
         end do
      end do
   end subroutine create_cache_copy_commands

end subroutine call_process_indices

subroutine job_request_variable(self,variable,copy_to_cache,copy_to_store,must_be_in_write_cache)
   class (type_job),target,      intent(inout) :: self
   type (type_internal_variable),intent(inout),target :: variable
   logical,optional,             intent(in)    :: copy_to_cache
   logical,optional,             intent(in)    :: copy_to_store
   logical,optional,             intent(in)    :: must_be_in_write_cache

   type (type_variable_request), pointer :: variable_request

   if (self%state<job_state_created) call driver%fatal_error('job_request_variable','Job has not been created yet.')
   if (self%state>job_state_created) call driver%fatal_error('job_request_variable','Job "'//trim(self%name)//'" has already begun initialization; variables can no longer be requested.')

   ! Make sure this variable will not be merged
   variable%write_operator = operator_assign

   if (associated(variable%write_owner)) call driver%fatal_error('job_request_variable','BUG: requested variable is co-written.')

   allocate(variable_request)
   variable_request%variable => variable
   if (present(copy_to_store)) variable_request%copy_to_store = copy_to_store

   ! If the variable needs to be present in the write cache, make a note so we can copy written values across tasks/jobs if needed.
   if (present(must_be_in_write_cache)) then
      if (must_be_in_write_cache) then
         call self%required_in_write_cache%add(variable)
         variable_request%requirements = variable_request%requirements + 2
      end if
   end if

   ! If the variable needs to be present in the read cache, make sure that we load it from the store into the cache.
   if (present(copy_to_cache)) then
      if (copy_to_cache) variable_request%requirements = variable_request%requirements + 1
   end if

   variable_request%next => self%first_variable_request
   self%first_variable_request => variable_request
end subroutine job_request_variable

subroutine job_request_call(self,model,source)
   class (type_job),target,intent(inout) :: self
   class (type_base_model),intent(in),target :: model
   integer,                intent(in)    :: source

   type (type_call_request), pointer :: call_request

   if (self%state<job_state_created) call driver%fatal_error('job_request_call','Job has not been created yet.')
   if (self%state>job_state_created) call driver%fatal_error('job_request_call','Job "'//trim(self%name)//'" has already begun initialization; calls can no longer be requested.')

   allocate(call_request)
   call_request%model => model
   call_request%source = source
   call_request%next => self%first_call_request
   self%first_call_request => call_request
end subroutine job_request_call

subroutine job_create_graph(self)
   class (type_job),target,intent(inout) :: self

   type (type_node_list)                :: outer_calls
   type (type_variable_request),pointer :: variable_request, next_variable_request
   type (type_call_request),    pointer :: call_request, next_call_request
   type (type_node),            pointer :: graph_node
   type (type_variable_node),   pointer :: unresolved_dependency
   type (type_output_variable), pointer :: variable_node

   if (self%state<job_state_created) call driver%fatal_error('job_create_graph','This job has not been created yet.')
   if (self%state>=job_state_graph_created) call driver%fatal_error('job_create_graps','Graph for this job has already been created.')

   ! If we are linked to an earlier called job, make sure its graph has been created already.
   ! This is essential because we can skip calls if they appear already in the previous job - we determine this by exploring its graph.
   if (associated(self%previous)) then
      if (self%previous%state<job_state_graph_created) call driver%fatal_error('job_create_graph',trim(self%name)//': graph for previous job ('//trim(self%previous%name)//') has not been created yet.')
   end if
   if (associated(self%dependency_handler)) then
      if (self%dependency_handler%state>=job_state_graph_created) call driver%fatal_error('job_create_graph',trim(self%name)//': graph for dependency handling job ('//trim(self%dependency_handler%name)//') has already been created.')
   end if

   ! Construct the dependency graph by adding explicitly requested variables and calls.
   ! We clean up [deallocate] the variable and call requests at the same time.
   variable_request => self%first_variable_request
   do while (associated(variable_request))
      next_variable_request => variable_request%next
      if (variable_request%requirements /= 0) then
         if (iand(variable_request%requirements,1)/=0) call self%own_node%inputs%add(variable_request%variable)
         if (iand(variable_request%requirements,2)/=0) then
            variable_node => self%own_node%outputs%add(variable_request%variable)
            !variable_node%prefill = prefill_previous_value
         end if
         call self%graph%add_variable(variable_request%variable, outer_calls, copy_to_store=variable_request%copy_to_store, caller=self%own_node)
      else
         call self%graph%add_variable(variable_request%variable, outer_calls, copy_to_store=variable_request%copy_to_store)
      end if
      deallocate(variable_request)
      variable_request => next_variable_request
   end do
   self%first_variable_request => null()
   call_request => self%first_call_request
   do while (associated(call_request))
      next_call_request => call_request%next
      graph_node => self%graph%add_call(call_request%model,call_request%source,outer_calls,ignore_dependencies=self%ignore_dependencies)
      deallocate(call_request)
      call_request => next_call_request
   end do
   self%first_call_request => null()
   call outer_calls%finalize()

   self%graph%frozen = .true.

   if (associated(self%dependency_handler)) then
      unresolved_dependency => self%graph%unresolved_dependencies%first
      do while (associated(unresolved_dependency))
         call self%dependency_handler%request_variable(unresolved_dependency%target, copy_to_store=.true.)
         unresolved_dependency => unresolved_dependency%next
      end do
   end if

   ! Save graph
   !open(newunit=unit,file='graph.gv',action='write',status='replace',iostat=ios)
   !call self%graph%save_as_dot(unit)
   !close(unit)

   self%state = job_state_graph_created
end subroutine job_create_graph

subroutine job_create_tasks(self)
   class (type_job),target,intent(inout) :: self

   type (type_graph_subset_node_set)    :: subset
   type (type_task_tree_node)           :: root
   class (type_task_tree_node), pointer :: leaf
   integer                              :: ntasks
   type (type_task),            pointer :: task
   integer                              :: unit,ios
   type (type_call),            pointer :: call_node, next_call_node

   if (self%state>=job_state_tasks_created) call driver%fatal_error('job_create_tasks',trim(self%name)//': tasks for this job have already been created.')
   if (self%state<job_state_graph_created) call driver%fatal_error('job_create_tasks',trim(self%name)//': the graph for this job have not been created yet.')

   ! If we are linked to an earlier called job, make sure its task list has already been created.
   ! This is essential if we will try to outsource our own calls to previous tasks/jobs.
   if (associated(self%previous)) then
      if (self%previous%state<job_state_tasks_created) call driver%fatal_error('job_create_tasks',trim(self%name)//': tasks for previous job '//trim(self%previous%name)//' have not been created yet.')
   end if

   ! Create tree that describes all possible task orders (root is at the right of the call list, i.e., at the end of the very last call)
   call create_graph_subset_node_set(self%graph,subset)
   if (associated(self%own_node)) then
      root%operation = self%own_node%source
      call subset%collect_and_branch(root)
   else
      root%operation = source_unknown
      call subset%branch(root)
   end if

   ! Determine optimum task order and use it to create the task objects.
   ntasks = root%get_minimum_depth()-1
   write (*,'(a,i0,a)') 'Best task order contains ',ntasks,' tasks.'
   leaf => root%get_leaf_at_depth(ntasks+1)
   write (*,'(a,a)') 'Best task order: ',trim(leaf%to_string())
   call create_tasks(leaf)
   if (associated(subset%first)) call driver%fatal_error('job_select_order','BUG: graph subset should be empty after create_tasks.')

   if (self%outsource_tasks) then
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

   if (associated(self%own_node)) then
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
   end if

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
         if (.not.associated(call_node%model)) call_node%source = source_none
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

      write (*,*) 'moving '//trim(call_node%graph_node%as_string())
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

         !call target_task%initialize()
      end if
   end subroutine move_call_backwards

   function move_one_step_backwards(job,task,travelling_call) result(moved)
      class (type_job),pointer :: job
      type (type_task),pointer :: task
      type (type_call),intent(in) :: travelling_call
      logical :: moved

      type (type_call),pointer             :: call_node
      type (type_node_set_member), pointer :: dependency

      moved = .false.

      ! First determine if we can leave the current task
      ! (we cannot if it also handles one or more of our dependencies)
      call_node => task%first_call
      do while (associated(call_node))
         dependency => travelling_call%graph_node%dependencies%first
         do while (associated(dependency))
            if (associated(dependency%p,call_node%graph_node)) then
               write (*,*) '  cannot move past '//trim(call_node%graph_node%as_string())
               return
            end if
            dependency => dependency%next
         end do
         call_node => call_node%next
      end do

      ! Move to previous task (if any)
      call get_previous_task(job,task)
      if (.not.associated(task)) return

      moved = .true.
   end function move_one_step_backwards

end subroutine job_create_tasks

function get_last_task(job) result(task)
   class (type_job),intent(in) :: job
   type (type_task),pointer :: task

   task => job%final_task
   if (associated(task)) return

   task => job%first_task
   do while (associated(task%next))
      task => task%next
   end do
end function get_last_task

subroutine get_previous_task(job,task)
   class (type_job),pointer :: job
   type (type_task),pointer :: task

   type (type_task),pointer :: old_task

   if (associated(job%first_task,task).or.(.not.associated(job%first_task).and.associated(job%final_task,task))) then
      ! We are the first task in the job - move to previous job (and its last task)
      task => null()
      job => job%previous
      do while (associated(job))
         if (associated(job%first_task).or.associated(job%final_task)) exit
         job => job%previous
      end do
      if (associated(job)) task => get_last_task(job)
   else
      ! We are NOT the first task in the job - find the previous task
      old_task => task
      task => job%first_task
      do while (associated(task%next).and..not.associated(task%next,old_task))
         task => task%next
      end do
   end if
end subroutine get_previous_task

subroutine job_finalize_prefill_settings(self)
   class (type_job),target,intent(inout) :: self

   type (type_node_list_member),pointer :: graph_node
   type (type_output_variable), pointer :: output_variable, output_variable1, final_output_variable
   type (type_task),            pointer :: task
   class (type_job),            pointer :: current_job
   type (type_call),            pointer :: call_node
   type (type_variable_node),   pointer :: variable_node
   logical                              :: last_task
   logical                              :: copy_to_cache
   type (type_node_set)                 :: remaining_cowriters, found_cowriters

   if (self%state>=job_state_finalized_prefill_settings) call driver%fatal_error('job_finalize_prefill_settings','This job has already been initialized.')
   if (self%state<job_state_tasks_created) call driver%fatal_error('job_finalize_prefill_settings','Tasks for this job have not been created yet.')

   ! Set copy-to-cache and copy-to-store based on dependencies between different calls/tasks.
   task => self%first_task
   do while (associated(task))
      call prepare_task(task)
      task => task%next
   end do
   if (associated(self%final_task)) call prepare_task(self%final_task)

   ! Loop over all graph nodes, for each, enumerate outputs.
   ! For outputs that are jointly written, find all tasks that contribute to the variable value (in current and previous jobs)
   ! The very first task will can keep its default prefill setting (typically prefill to zero for summations), but all
   ! subsequent tasks need to prefill with the previous value.
   graph_node => self%graph%first
   do while (associated(graph_node))
      output_variable => graph_node%p%outputs%first
      do while (associated(output_variable))
         if (associated(output_variable%target%cowriters%first)) then
            ! This is the result of a reduction operation (e.g., sum) represented by output_variable.
            ! The result will be written by the controller (output_variable) and its co-writers (output_variable%target%cowriters).
            call remaining_cowriters%update(output_variable%group)
            call remaining_cowriters%add(output_variable%target)

            ! We have a set with all written variables that contribute to the current output.
            ! Now move from task to task, and for all but the earlier writer, prefill with the previous value from the store.
            last_task = .true.
            current_job => self
            task => get_last_task(self)
            do while (associated(task))
               ! Loop over all calls in this task
               call_node => task%first_call
               do while (associated(call_node))
                  ! Loop over all outputs of this call
                  output_variable1 => call_node%graph_node%outputs%first
                  do while (associated(output_variable1))
                     if (remaining_cowriters%contains(output_variable1%target)) then
                        ! This output is co-written to the target variable we are interested in.
                        call found_cowriters%add(output_variable1%target)
                        call remaining_cowriters%remove(output_variable1%target)
                        if (last_task) then
                           ! Keep a link to the very last output variable that contributes to the result,
                           ! and remember whether that needs to be copied to cache (copy-to-cache is reset by default below)
                           final_output_variable => output_variable1
                           copy_to_cache = final_output_variable%copy_to_cache
                        else
                           ! This is not the last task contributing. Make sure the result of this task is copied to the store,
                           ! so it is availble for subsequent tasks (which will use it as prefill).
                           output_variable1%copy_to_store = .true.
                        end if
                        output_variable1%copy_to_cache = .false.
                     end if
                     output_variable1 => output_variable1%next
                  end do
                  call_node => call_node%next
               end do
               if (associated(remaining_cowriters%first)) then
                  ! Some variable contributions are remaining and are therefore handled by earlier tasks.
                  ! Any contributions in the current task must therefore prefill with the previous value.
                  variable_node => found_cowriters%first
                  do while (associated(variable_node))
                     variable_node%target%prefill = prefill_previous_value
                     variable_node => variable_node%next
                  end do
               end if
               if (associated(found_cowriters%first)) last_task = .false.
               call found_cowriters%finalize()
               if (.not.associated(remaining_cowriters%first)) exit
               call get_previous_task(current_job,task)
            end do
            if (associated(remaining_cowriters%first)) then
               variable_node => remaining_cowriters%first
               do while (associated(variable_node))
                  write (*,*) 'not found: '//trim(variable_node%target%name)
                  variable_node => variable_node%next
               end do
               call fatal_error('job_finalize_prefill_settings','BUG: not all sources of cowritten variable '//trim(output_variable%target%name)//' were found within owning and earlier jobs.')
            end if
            final_output_variable%copy_to_cache = copy_to_cache
         end if
         output_variable => output_variable%next
      end do
      graph_node => graph_node%next
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

      type (type_output_variable), pointer :: variable_node
      type (type_node_set_member), pointer :: dependent_call
      type (type_call),            pointer :: later_call

      ! For all output variables that other models are interested in, decide whether to copy their value
      ! from the write to read cache [if the other model will be called as part of the same task],
      ! of to save it to the persistent data store.
      variable_node => call_node%graph_node%outputs%first
      do while (associated(variable_node))
         ! For this output variable, loop over all calls that depend on it.
         dependent_call => variable_node%dependent_nodes%first
         do while (associated(dependent_call))
            ! For this dependent call, establish whether it is part of the same task, or of another task.
            later_call => call_node%next
            do while (associated(later_call))
               if (associated(later_call%graph_node,dependent_call%p)) exit
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
   end subroutine prepare_call

end subroutine job_finalize_prefill_settings

subroutine job_initialize(self,variable_register)
   class (type_job),target,      intent(inout) :: self
   type (type_variable_register),intent(inout) :: variable_register

   type (type_task), pointer :: task

   if (self%state>=job_state_initialized) call driver%fatal_error('job_initialize','This job has already been initialized.')
   if (self%state<job_state_finalized_prefill_settings) call driver%fatal_error('job_initialize','Prefill settings for this job have not been finalized yet.')

   ! Initialize tasks
   task => self%first_task
   do while (associated(task))
      call task%initialize(variable_register)
      task => task%next
   end do
   if (associated(self%final_task)) call self%final_task%initialize(variable_register)

   self%state = job_state_initialized

end subroutine job_initialize

subroutine job_process_indices(self,unfulfilled_dependencies)
   class (type_job),target, intent(inout) :: self
   type (type_variable_set),intent(inout) :: unfulfilled_dependencies

   type (type_task), pointer :: task

   task => self%first_task
   do while (associated(task))
      call task_process_indices(task,unfulfilled_dependencies)
      task => task%next
   end do
   if (associated(self%final_task)) call task_process_indices(self%final_task,unfulfilled_dependencies)
end subroutine job_process_indices

subroutine job_manager_create(self,job,name,final_operation,outsource_tasks,dependency_handler,ignore_dependencies)
   class (type_job_manager), intent(inout) :: self
   type (type_job),target,   intent(inout) :: job
   character(len=*),         intent(in)    :: name
   integer,optional,         intent(in)    :: final_operation
   logical,optional,         intent(in)    :: outsource_tasks
   type (type_job),target,optional         :: dependency_handler
   logical,optional                        :: ignore_dependencies

   type (type_job_manager_item),pointer :: node

   if (job%state>=job_state_created) call driver%fatal_error('job_manager_create','This job has already been created with name '//trim(job%name)//'.')
   job%state = job_state_created

   job%name = name
   if (present(final_operation)) then
      allocate(job%own_node)
      job%own_node%source = final_operation
      call job%graph%append(job%own_node)
   end if
   if (present(outsource_tasks)) job%outsource_tasks = outsource_tasks

   allocate(node)
   node%job => job
   node%next => self%first
   self%first => node

   if (.not.associated(self%default_dependency_handler,job)) job%dependency_handler => self%default_dependency_handler
   if (present(dependency_handler)) job%dependency_handler => dependency_handler
   if (present(ignore_dependencies)) job%ignore_dependencies = ignore_dependencies
end subroutine job_manager_create

subroutine job_manager_initialize(self,variable_register)
   class (type_job_manager),     intent(inout) :: self
   type (type_variable_register),intent(inout) :: variable_register

   type (type_job_manager_item),pointer :: node, first_ordered

   ! Order jobs according to call order.
   ! This ensures that jobs that are scheduled to run earlier are also initialized earlier.
   ! During initialization, a job can therefore expect any preceding jobs to have initialized completely.
   first_ordered => null()
   do while (associated(self%first))
      node => self%first
      self%first => self%first%next
      call add_to_order(node%job)
      deallocate(node)
   end do
   self%first => first_ordered

   ! Create all graphs. This must be done across all jobs before other operations that use graphs,
   ! since a job can add to graphs owned by other jobs.
   node => self%first
   do while (associated(node))
      call job_create_graph(node%job)
      node => node%next
   end do

   ! Create tasks. This must be done for all jobs before job_finalize_prefill_settings is called, as this API operates across all jobs.
   node => self%first
   do while (associated(node))
      call job_create_tasks(node%job)
      node => node%next
   end do

   ! Finalize prefill settings (this has cross-job implications, so must be done for all jobs before they initialize (initialization uses prefill settings)
   node => self%first
   do while (associated(node))
      call job_finalize_prefill_settings(node%job)
      node => node%next
   end do

   ! Initialize all jobs
   node => self%first
   do while (associated(node))
      call job_initialize(node%job,variable_register)
      node => node%next
   end do

contains

   recursive subroutine add_to_order(job)
      type (type_job), target :: job

      type (type_job_manager_item),pointer :: node

      ! Make sure job is not yet in list
      node => first_ordered
      do while (associated(node))
         if (associated(node%job,job)) return
         node => node%next
      end do

      ! If this job has a previous job, append that first.
      if (associated(job%previous)) then
         call add_to_order(job%previous)
      end if

      ! Add any other jobs that use the current one as their dependency handler
      node => self%first
      do while (associated(node))
         if (associated(node%job%dependency_handler,job)) call add_to_order(node%job)
         node => node%next
      end do

      ! Append to list
      if (associated(first_ordered)) then
         node => first_ordered
         do while (associated(node%next))
            node => node%next
         end do
         allocate(node%next)
         node%next%job => job
      else
         allocate(first_ordered)
         first_ordered%job => job
      end if
   end subroutine add_to_order

end subroutine job_manager_initialize

subroutine job_manager_process_indices(self,unfulfilled_dependencies)
   class (type_job_manager), intent(inout) :: self
   type (type_variable_set), intent(out)   :: unfulfilled_dependencies

   type (type_job_manager_item),pointer :: node

   node => self%first
   do while (associated(node))
      call job_process_indices(node%job,unfulfilled_dependencies)
      node => node%next
   end do
end subroutine job_manager_process_indices

subroutine job_manager_print(self,specific_variable)
   class (type_job_manager),                        intent(in) :: self
   type (type_internal_variable), target, optional, intent(in) :: specific_variable

   type (type_job_manager_item),pointer :: node

   node => self%first
   do while (associated(node))
      call node%job%print(specific_variable)
      node => node%next
   end do
end subroutine job_manager_print

subroutine job_set_previous(self,previous)
   class (type_job), intent(inout), target :: self
   type (type_job),  intent(inout), target :: previous

   if (self%state>job_state_created) &
      call driver%fatal_error('job_set_previous','This job ('//trim(self%name)//') has already started initialization; it is too late to specify its place in the call order.')
   if (associated(self%previous)) &
      call driver%fatal_error('job_set_previous','This job has already been connected to a subsequent one.')

   self%previous => previous
   self%graph%previous => previous%graph
end subroutine job_set_previous

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

recursive subroutine variable_register_add_store(self,variable)
   class (type_variable_register),intent(inout) :: self
   type (type_internal_variable), target        :: variable

   integer                            :: i
   type (type_variable_node), pointer :: variable_node

   if (variable%store_index/=store_index_none) return

   if (associated(variable%write_owner)) then
      call self%add_store(variable%write_owner)
      return
   end if
   
   select case (variable%domain)
   case (domain_interior)
      call self%interior_store%append(variable,index=i)
   case (domain_horizontal,domain_surface,domain_bottom)
      call self%horizontal_store%append(variable,index=i)
   end select

   variable%store_index = i
   variable_node => variable%cowriters%first
   do while (associated(variable_node))
      variable_node%target%store_index = i
      variable_node => variable_node%next
   end do
end subroutine variable_register_add_store

recursive subroutine variable_register_add_read(self,variable)
   class (type_variable_register),intent(inout) :: self
   type (type_internal_variable), target        :: variable

   integer                            :: i
   type (type_variable_node), pointer :: variable_node

   if (variable%read_indices%value/=-1) return

   if (associated(variable%write_owner)) then
      call self%add_read(variable%write_owner)
      return
   end if
   
   select case (variable%domain)
   case (domain_interior)
      call self%interior_read%append(variable,index=i)
   case (domain_horizontal,domain_surface,domain_bottom)
      call self%horizontal_read%append(variable,index=i)
   case (domain_scalar)
      call self%scalar_read%append(variable,index=i)
   end select

   call variable%read_indices%set_value(i)
   !variable_node => variable%cowriters%first
   !do while (associated(variable_node))
   !   call variable_node%target%read_indices%set_value(i)
   !   variable_node => variable_node%next
   !end do
end subroutine variable_register_add_read

subroutine variable_register_print(self)
   class (type_variable_register),intent(in) :: self

   call print_list('Interior store:',  self%interior_store)
   call print_list('Horizontal store:',self%horizontal_store)
   call print_list('Interior read:',   self%interior_read)
   call print_list('Horizontal read:', self%horizontal_read)
   call print_list('Scalar read:',     self%scalar_read)

contains

   subroutine print_list(title,list)
      character(len=*),         intent(in) :: title
      type (type_variable_list),intent(in) :: list

      integer                            :: i
      type (type_variable_node), pointer :: variable_node

      write (*,'(a)') title
      i = 0
      variable_node => list%first
      do while (associated(variable_node))
         i = i + 1
         write (*,'("  ",i0,": ",a)') i,trim(variable_node%target%name)
         variable_node => variable_node%next
      end do
   end subroutine

end subroutine variable_register_print

end module fabm_job

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

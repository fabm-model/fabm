#include "fabm_driver.h"
#include "fabm_private.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_job: derived types that describe a single job (module subroutines to call and their dependencies)
!
! !INTERFACE:
module fabm_job

   use fabm_types
   use fabm_driver

   implicit none

   private

   public type_job, type_call
   public find_dependencies

   type type_cache_copy_command
      integer :: read_index
      integer :: write_index
   end type

   type type_output_variable
      type (type_internal_variable),pointer :: target               => null()
      type (type_call_pointer),     pointer :: first_dependent_node => null()
      logical                               :: copy_to_cache        = .false.
      logical                               :: copy_to_store        = .false.
      type (type_output_variable),pointer   :: next                 => null()
   end type

   type type_output_variable_set
      type (type_output_variable), pointer :: first => null()
   contains
      procedure :: add      => output_variable_set_add
      procedure :: finalize => output_variable_set_finalize
   end type

   type type_input_variable
      type (type_internal_variable),pointer :: target => null()
      type (type_input_variable),   pointer :: next   => null()
   end type

   type type_input_variable_set
      type (type_input_variable),pointer :: first => null()
   contains
      procedure :: add      => input_variable_set_add
      procedure :: finalize => input_variable_set_finalize
   end type

   type type_call_pointer
      type (type_call),       pointer :: p    => null()
      type(type_call_pointer),pointer :: next => null()
   end type

   type type_call
      class (type_base_model), pointer            :: model => null()
      integer                                     :: source = source_unknown
      type (type_cache_copy_command), allocatable :: copy_commands_int(:) ! interior variables to copy from write to read cache after call completes
      type (type_cache_copy_command), allocatable :: copy_commands_hz(:)  ! horizontal variables to copy from write to read cache after call completes
      type (type_input_variable_set)              :: inputs               ! input variables
      type (type_output_variable_set)             :: outputs              ! output variables that a later called model requires
      type (type_call), pointer                   :: next     => null()
      type (type_call), pointer                   :: previous => null()
   contains
      procedure :: initialize => call_initialize
      procedure :: finalize   => call_finalize
   end type type_call

   type type_call_list
      type (type_call),     pointer :: first => null()
      type (type_call),     pointer :: last  => null()
      type (type_call_list),pointer :: prior => null()
   contains
      procedure :: find        => call_list_find
      procedure :: append_node => call_list_append_node
      procedure :: append_call => call_list_append_call
      procedure :: remove      => call_list_remove
      generic :: append => append_node,append_call
      procedure :: print       => call_list_print
      procedure :: initialize  => call_list_initialize
      procedure :: request_variable  => call_list_request_variable
      procedure :: request_variables => call_list_request_variables
      procedure :: finalize    => call_list_finalize
   end type type_call_list

   type type_job
      integer :: domain = -1
      type (type_call_list) :: calls
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
   contains
      procedure :: initialize => job_initialize
      procedure :: print      => job_print
   end type

   contains

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
            variable_node => call_node%outputs%first
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
            variable_node => call_node%outputs%first
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
            variable_node => call_node%outputs%first
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
            variable_node => call_node%outputs%first
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
            input_variable => call_node%inputs%first
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
            input_variable => call_node%inputs%first
            do while (associated(input_variable))
               if (iand(input_variable%target%domain,domain)/=0) load(input_variable%target%read_indices%value) = .true.
               input_variable => input_variable%next
            end do
            call_node => call_node%next
         end do
      end subroutine create_load_commands

   end subroutine job_initialize

   subroutine job_print(self)
      class (type_job), intent(in) :: self
      call self%calls%print()
   end subroutine

subroutine call_finalize(self)
   class (type_call), intent(inout) :: self

   if (allocated(self%copy_commands_int)) deallocate(self%copy_commands_int)
   if (allocated(self%copy_commands_hz))  deallocate(self%copy_commands_hz)
   call self%inputs%finalize()
   call self%outputs%finalize()
end subroutine call_finalize

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
      deallocate(node)
      node => next
   end do
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
end subroutine input_variable_set_finalize

subroutine call_list_initialize(call_list)
   class (type_call_list),intent(inout) :: call_list

   type (type_call),pointer :: node

   node => call_list%first
   do while (associated(node))
      call node%initialize()
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

subroutine call_list_append_node(self,node)
   class (type_call_list), intent(inout) :: self
   type (type_call),pointer              :: node

   if (associated(self%first)) then
      ! List contains one or more items - append to tail.
      node%previous => self%last
      if (associated(self%last)) self%last%next => node
      self%last => node
   else
      ! List is empty - new node becomes both head and tail.
      node%previous => null()
      self%first => node
      self%last => node
   end if
   node%next => null()
end subroutine call_list_append_node

subroutine call_list_remove(self,node)
   class (type_call_list), intent(inout) :: self
   type (type_call),pointer              :: node

   type (type_call),pointer :: current

   ! Make sure the supplied node occurs in our call list to begin with.
   current => self%first
   do while (associated(current))
      if (associated(current,node)) exit
      current => current%next
   end do
   if (.not.associated(current)) call driver%fatal_error('call_list_remove','call list does not contain supplied node.')

   if (associated(node%previous)) then
      node%previous%next => node%next
      node%previous => null()
   else
      self%first => node%next
   end if

   if (associated(node%next)) then
      node%next%previous => node%previous
      node%next => null()
   else
      self%last => node%previous
   end if
end subroutine call_list_remove

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

subroutine call_list_print(self)
   class (type_call_list), intent(in) :: self

   type (type_call),           pointer :: node
   type (type_output_variable),pointer :: variable
   type (type_call_pointer),   pointer :: pnode

   node => self%first
   do while (associated(node))
      write (*,'(a,": ",a)') trim(node%model%get_path()),trim(source2string(node%source))
      variable => node%outputs%first
      do while (associated(variable))
         write (*,'("   ",a,",write@",i0)',advance='no') trim(variable%target%name),variable%target%write_indices%value
         if (variable%copy_to_cache) write (*,'(",cache@",i0)',advance='no') variable%target%read_indices%value
         if (variable%copy_to_store) write (*,'(",store@",i0)',advance='no') variable%target%store_index
         if (variable%target%prefill==prefill_missing_value) write (*,'(",prefill=",g0.6)',advance='no') variable%target%missing_value
         if (variable%target%prefill==prefill_previous_value) write (*,'(",prefill=previous")',advance='no')
         write (*,*)
         pnode => variable%first_dependent_node
         do while (associated(pnode))
            write (*,'("     <- ",a,": ",a)') trim(pnode%p%model%get_path()),trim(source2string(pnode%p%source))
            pnode => pnode%next
         end do
         variable => variable%next
      end do
      node => node%next
   end do

   contains

      character(len=32) function source2string(source)
         integer, intent(in) :: source
         select case (source)
         case (source_unknown);               source2string = 'unknown'
         case (source_do);                    source2string = 'do'
         case (source_do_column);             source2string = 'do_column'
         case (source_do_bottom);             source2string = 'do_bottom'
         case (source_do_surface);            source2string = 'do_surface'
         case (source_none);                  source2string = 'none'
         case (source_get_vertical_movement); source2string = 'get_vertical_movement'
         case (source_do_horizontal);         source2string = 'do_horizontal'
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
   nullify(self%last)
end subroutine call_list_finalize

subroutine call_initialize(self)
   class (type_call),intent(inout) :: self

   class (type_base_model),     pointer :: parent
   class (type_model_list_node),pointer :: model_list_node
   type (type_output_variable), pointer :: variable_node
   type (type_call_pointer),    pointer :: dependent_call
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

   ! For all output variables that other model are interested in, decide whether to copy their value
   ! from the write to read cache [if the other model will be called as part of the same job],
   ! of to save it to the persistent data store.
   variable_node => self%outputs%first
   do while (associated(variable_node))
      ! For this output variable, loop over all calls that depend on it.
      dependent_call => variable_node%first_dependent_node
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
      variable => self%outputs%first
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
         variable => self%outputs%first
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

recursive subroutine find_dependencies2(self,source,list,node,not_stale,forbidden)
   class (type_base_model),     intent(in),target :: self
   integer,                     intent(in)        :: source
   class (type_call_list),      intent(inout)     :: list
   type (type_call),  pointer                     :: node
   logical,                     intent(in)        :: not_stale
   type (type_call_list),target,intent(inout)     :: forbidden

   type (type_link),          pointer :: link
   type (type_input_variable),pointer :: input_variable
   character(len=2048)                :: chain
   logical                            :: same_source

   ! Check whether we are already processed this call.
   node => list%find(self,source,not_stale)
   if (associated(node)) return

   ! Check the list of forbidden model (i.e., models that indirectly request the current model)
   ! If the current model is on this list, there is a circular dependency between models.
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

   ! Create object to represent this call
   allocate(node)
   node%model => self
   node%source = source

   ! First add this call to the list of requesting calls [a list of all calls higher on the call stack]
   ! This forbids any indirect dependency on this call, as such would be a circular dependency.
   call forbidden%append(node)

   ! Collect all input variabes.
   link => self%links%first
   do while (associated(link))
      if (index(link%name,'/')==0.and.associated(link%original%read_index)) call node%inputs%add(link%target)
      link => link%next
   end do

   ! Loop over all input variables, and if they are written by another model call, make sure that model call is done first.
   input_variable => node%inputs%first
   do while (associated(input_variable))
      same_source = link%target%source==source .or. (link%target%source==source_unknown.and.(source==source_do_surface.or.source==source_do_bottom))
      if (.not.input_variable%target%write_indices%is_empty() &                       ! it is a variable that is written by a module...
          .and..not.(associated(input_variable%target%owner,self).and.same_source)) & ! ...and not ourselves
         call list%request_variable(input_variable%target,forbidden=forbidden,requester=node)
      input_variable => input_variable%next
   end do

   ! Clean up the list with forbidden calls (which would imply a circular dependency)
   call forbidden%remove(node)

   ! We're happy - add ourselves to the call list.
   call list%append(node)

end subroutine find_dependencies2

recursive subroutine call_list_request_variable(self,variable,copy_to_cache,copy_to_store,not_stale,forbidden,requester)
   class (type_call_list),       intent(inout)          :: self
   type (type_internal_variable),intent(in)             :: variable
   logical,                      intent(in),   optional :: copy_to_cache
   logical,                      intent(in),   optional :: copy_to_store
   logical,                      intent(in),   optional :: not_stale
   type (type_call_list),target, intent(inout),optional :: forbidden
   type (type_call),target,                    optional :: requester

   type (type_call),pointer :: node
   type (type_call_list),pointer :: forbidden_

   if (present(forbidden)) then
      forbidden_ => forbidden
   else
      allocate(forbidden_)
   end if

   ! If this variable is not written [but a field provided by the host or a state variable], return immediately.
   if (variable%write_indices%is_empty()) return

   if (variable%source==source_unknown) then
      ! This variable is either written by do_surface or do_bottom - which one of these two APIs is unknown.
      call find_dependencies2(variable%owner,source_do_surface,self,node,not_stale,forbidden_)
      call register_dependencies()
      call find_dependencies2(variable%owner,source_do_bottom,self,node,not_stale,forbidden_)
      call register_dependencies()
   elseif (variable%source/=source_none) then
      ! This variable is written by a known BGC API [is is not a constant indicated by source_none]
      call find_dependencies2(variable%owner,variable%source,self,node,not_stale,forbidden_)
      call register_dependencies()
   end if

   if (.not.present(forbidden)) then
      call forbidden_%finalize()
      deallocate(forbidden_)
   end if

   contains

   subroutine register_dependencies()
      type (type_call_pointer),   pointer :: requester_pointer
      type (type_output_variable),pointer :: variable_node

      variable_node => node%outputs%add(variable)
      if (present(copy_to_cache)) variable_node%copy_to_cache = variable_node%copy_to_cache .or. copy_to_cache
      if (present(copy_to_store)) variable_node%copy_to_store = variable_node%copy_to_store .or. copy_to_store
      if (present(requester)) then
         requester_pointer => variable_node%first_dependent_node
         do while (associated(requester_pointer))
            if (associated(requester_pointer%p,requester)) exit
            requester_pointer => requester_pointer%next
         end do
         if (.not.associated(requester_pointer)) then
            ! This dependent model was not included among the variable's dependent calls yet. Prepend it to the list.
            allocate(requester_pointer)
            requester_pointer%p => requester
            requester_pointer%next => variable_node%first_dependent_node
            variable_node%first_dependent_node => requester_pointer
         end if
      end if
   end subroutine register_dependencies
end subroutine call_list_request_variable

subroutine call_list_request_variables(self,link_list,copy_to_cache,copy_to_store,not_stale)
   class (type_call_list),intent(inout)       :: self
   type (type_link_list), intent(in)          :: link_list
   logical,               intent(in),optional :: copy_to_cache
   logical,               intent(in),optional :: copy_to_store
   logical,               intent(in),optional :: not_stale

   type (type_link), pointer :: link

   link => link_list%first
   do while (associated(link))
      call self%request_variable(link%target,copy_to_cache,copy_to_store,not_stale)
      link => link%next
   end do
end subroutine call_list_request_variables

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

subroutine partition_calls(self)
   class (type_call_list),intent(inout) :: self

   type (type_call_list), pointer :: current_call_list
   type (type_call),      pointer :: call_node
   integer                        :: source

   source = source_unknown
   call_node => self%first
   do while (associated(call_node))
      if (.not.is_source_compatible(source,call_node%source)) then
         ! New call list needed
         allocate(current_call_list)
      end if
      ! Add to current call list
      call_node => call_node%next
   end do
end subroutine partition_calls

logical function is_source_compatible(current_source,new_source)
   integer,intent(inout) :: current_source
   integer,intent(in)    :: new_source

   integer :: real_new_source

   real_new_source = new_source
   if (real_new_source==source_get_vertical_movement) real_new_source = source_do

   if (current_source==source_unknown) then
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
      call driver%fatal_error('source2operator','unknown source')
   end select
end function is_source_compatible

end module fabm_job

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

module fabm_coupling
   use fabm_types
   use fabm_standard_variables
   use fabm_properties
   use fabm_builtin_models
   use fabm_driver

   implicit none

   private

   public freeze_model_info, find_dependencies
   public collect_aggregate_variables, type_aggregate_variable_list, type_aggregate_variable
   public type_call_list_node, type_call_list, find_variable_dependencies

   logical,parameter :: debug_coupling = .false.

   integer,parameter :: couple_explicit                     = 0
   integer,parameter :: couple_aggregate_standard_variables = 1
   integer,parameter :: couple_flux_sums                    = 2
   integer,parameter :: couple_final                        = 3

   type type_copy_command
      integer :: read_index
      integer :: write_index
   end type

   type type_variable_set_node
      type (type_internal_variable), pointer :: target => null()
      type (type_variable_set_node), pointer :: next   => null()
   end type

   type type_variable_set
      type (type_variable_set_node), pointer :: first => null()
   contains
      procedure :: add      => variable_set_add
      procedure :: contains => variable_set_contains
      procedure :: finalize => variable_set_finalize
   end type

   type type_call_list_node
      class (type_base_model),          pointer :: model => null()
      logical                                   :: active = .true.
      integer                                   :: source = source_unknown
      type (type_copy_command), allocatable     :: copy_commands_int(:)
      type (type_copy_command), allocatable     :: copy_commands_hz(:)
      type (type_variable_set)                  :: written_variables
      type (type_variable_set)                  :: computed_variables
      type (type_call_list_node), pointer :: next => null()
   contains
      procedure :: finalize   => call_list_node_finalize
   end type

   type type_call_list
      type (type_call_list_node), pointer :: first => null()
   contains
      procedure :: find       => call_list_find
      procedure :: filter     => call_list_filter
      procedure :: append     => call_list_append
      procedure :: extend     => call_list_extend
      procedure :: print      => call_list_print
      procedure :: initialize => call_list_initialize
      procedure :: finalize   => call_list_finalize
      procedure :: computes   => call_list_computes
   end type

   type type_contributing_variable
      type (type_link),                 pointer :: link               => null()
      real(rk)                                  :: scale_factor       = 1.0_rk
      logical                                   :: include_background = .false.
      type (type_contributing_variable),pointer :: next               => null()
   end type

   type type_aggregate_variable
      type (type_bulk_standard_variable), pointer :: standard_variable           => null()
      type (type_contributing_variable),  pointer :: first_contributing_variable => null()
      type (type_aggregate_variable),     pointer :: next                        => null()
   end type

   type type_aggregate_variable_list
      type (type_aggregate_variable), pointer :: first => null()
   contains
      procedure :: get   => aggregate_variable_list_get
      procedure :: print => aggregate_variable_list_print
   end type

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Make model information read-only.
!
! !INTERFACE:
   subroutine freeze_model_info(self)
!
! !DESCRIPTION:
!  This function finalizes model initialization. It will resolve all remaining
!  internal dependencies (coupling commands) and generate final authorative lists
!  of state variables, diagnostic variables, conserved quantities and readable
!  variables ("dependencies").
!
! !INPUT/OUTPUT PARAMETER:
      class (type_base_model),intent(inout),target :: self
!
!EOP
!-----------------------------------------------------------------------
!BOC
      if (associated(self%parent)) call self%fatal_error('freeze_model_info', &
         'BUG: freeze_model_info can only operate on the root model.')

      call before_coupling(self)

      ! Coupling stage 1: implicit - couple variables based on overlapping standard identities.
      call couple_standard_variables(self)

      ! Coupling stage 2: explicit - resolve user- or model-specified links between variables.
      call process_coupling_tasks(self,couple_explicit)

      ! Coupling stage 3: try to automatically fulfil remaining dependencies on aggregate standard variables
      ! This may create new child models to handle the necessary summations.
      ! Note that these child models may also add source terms (if the sum is treated as state variable),
      ! so any source term treatment should happen after this is complete!
      call process_coupling_tasks(self,couple_aggregate_standard_variables)

      ! Create models for aggregate variables at root level, to be used to compute conserved quantities.
      ! After this step, the set of variables that contribute to aggregate quantities may not be modified.
      ! That is, no new such variables may be added, and no such variables may be coupled.
      call create_aggregate_models(self)

      ! Perform coupling for any new aggregate models.
      ! This may append items to existing lists of source terms and bottom/surface fluxes,
      ! so it has to proceed couple_flux_sums
      call process_coupling_tasks(self,couple_explicit)

      ! Now that coupling for non-rate variables is complete, contributions to aggregate quantities
      ! (including ones of slave variables) are final.
      ! Create conservation checks where needed.
      call create_conservation_checks(self)

      ! Create summations of source terms and surface/bottom fluxes where they are needed,
      ! and then process the resulting coupling tasks.
      call process_coupling_tasks(self,couple_flux_sums)
      call process_coupling_tasks(self,couple_final)

      ! Allow inheriting models to perform additional tasks after coupling.
      call after_coupling(self)

      ! Check whether units of coupled variables match
      !call check_coupling_units(self)

      call freeze(self)
   end subroutine freeze_model_info
!EOC

recursive subroutine before_coupling(self)
   class (type_base_model),intent(inout),target :: self

   type (type_model_list_node), pointer :: node

   call self%before_coupling()
   node => self%children%first
   do while (associated(node))
      call before_coupling(node%model)
      node => node%next
   end do
end subroutine

recursive subroutine after_coupling(self)
   class (type_base_model),intent(inout),target :: self

   type (type_model_list_node), pointer :: node

   call self%after_coupling()
   node => self%children%first
   do while (associated(node))
      call after_coupling(node%model)
      node => node%next
   end do
end subroutine

recursive subroutine freeze(self)
   class (type_base_model),intent(inout),target :: self

   type (type_model_list_node), pointer :: node

   self%frozen = .true.
   node => self%children%first
   do while (associated(node))
      call freeze(node%model)
      node => node%next
   end do
end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Automatically couple variables that represent the same standard variable.
!
! !INTERFACE:
   subroutine couple_standard_variables(model)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_base_model),intent(inout),target :: model
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_link),pointer :: link,first_link
      type (type_standard_variable_set) :: standard_variables
      type (type_standard_variable_node), pointer :: node
!
!-----------------------------------------------------------------------
!BOC
      ! Build a list of all unique standard variables.
      link => model%links%first
      do while (associated(link))
         call standard_variables%update(link%target%standard_variables)
         link => link%next
      end do

      ! Looop over all unique standard variable and collect and couple associated model variables.
      node => standard_variables%first
      do while (associated(node))
         first_link => null()
         link => model%links%first
         do while (associated(link))
            if (link%target%standard_variables%contains(node%p)) then
               if (associated(first_link)) then
                  if (link%target%write_indices%is_empty().and..not.(first_link%target%presence/=presence_internal.and.link%target%presence==presence_internal)) then
                     ! Default coupling: early variable (first_link) is master, later variable (link) is slave.
                     call couple_variables(model,first_link%target,link%target)
                  else
                     ! Later variable (link) is write-only and therefore can only be master. Try coupling with early variable (first_link) as slave.
                     call couple_variables(model,link%target,first_link%target)
                  end if
               else
                  first_link => link
               end if
            end if
            link => link%next
         end do
         node => node%next
      end do
      call standard_variables%finalize()
   end subroutine couple_standard_variables
!EOC

   subroutine collect_user_specified_couplings(self)
      class (type_base_model),intent(inout) :: self

      type (type_link),         pointer :: link
      class (type_property),    pointer :: master_name
      class (type_coupling_task),pointer :: task

      link => self%links%first
      do while (associated(link))
         ! Only process own links (those without slash in the name)
         if (index(link%name,'/')==0) then
            master_name => self%couplings%find_in_tree(link%name)
            if (associated(master_name)) then
               call self%coupling_task_list%add(link,.true.,task)
               task%user_specified = .true.
               select type (master_name)
                  class is (type_string_property)
                     task%master_name = master_name%value
               end select
            end if    ! Coupling provided
         end if   ! Our own link, which may be coupled
         link => link%next
      end do

      self%coupling_task_list%includes_custom = .true.
   end subroutine collect_user_specified_couplings

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process all model-specific coupling tasks.
!
! !INTERFACE:
   recursive subroutine process_coupling_tasks(self,stage)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_base_model),intent(inout),target :: self
      integer,                intent(in)           :: stage
!
!EOP
!
! !LOCAL VARIABLES:
      class (type_base_model),      pointer :: root
      type (type_model_list_node),  pointer :: child
      class (type_coupling_task),   pointer :: coupling, next_coupling
      type (type_internal_variable),pointer :: master
      type (type_link),             pointer :: link
!
!-----------------------------------------------------------------------
!BOC
      ! Find root model, which will handle the individual coupling tasks.
      root => self
      do while (associated(root%parent))
         root => root%parent
      end do

      if (.not.self%coupling_task_list%includes_custom) call collect_user_specified_couplings(self)

      ! For each variable, determine if a coupling command is provided.
      coupling => self%coupling_task_list%first
      do while (associated(coupling))

         nullify(master)
         select case (stage)
            case (couple_explicit,couple_final)
               if (associated(coupling%master_standard_variable)) then
                  ! Coupling to a standard variable - first try to find the corresponding standard variable.
                  link => root%links%first
                  do while (associated(link))
                     if (link%target%standard_variables%contains(coupling%master_standard_variable)) exit
                     link => link%next
                  end do

                  if (stage==couple_final.and..not.associated(link)) then
                     ! This is our last chance - create an appropriate variable at the root level.
                     select type (standard_variable=>coupling%master_standard_variable)
                     class is (type_bulk_standard_variable)
                        call root%add_interior_variable(standard_variable%name, standard_variable%units, standard_variable%name, &
                                                    standard_variable=standard_variable, presence=presence_external_optional, link=link)
                     class is (type_horizontal_standard_variable)
                        call root%add_horizontal_variable(standard_variable%name, standard_variable%units, standard_variable%name, &
                                                          standard_variable=standard_variable, presence=presence_external_optional, link=link)
                     class is (type_global_standard_variable)
                        call root%add_scalar_variable(standard_variable%name, standard_variable%units, standard_variable%name, &
                                                      standard_variable=standard_variable, presence=presence_external_optional, link=link)
                     end select
                  end if
                  if (associated(link)) master => link%target
               else
                  ! Try to find the master variable among the variables of the requesting model or its parents.
                  if (coupling%slave%name/=coupling%master_name) then
                     ! Master and slave name differ: start master search in current model, then move up tree.
                     master => self%find_object(coupling%master_name,recursive=.true.,exact=.false.)
                  elseif (associated(self%parent)) then
                     ! Master and slave name are identical: start master search in parent model, then move up tree.
                     master => self%parent%find_object(coupling%master_name,recursive=.true.,exact=.false.)
                  else
                     call self%fatal_error('process_coupling_tasks', &
                        'Master and slave name are identical: "'//trim(coupling%master_name)//'". This is not valid at the root of the model tree.')
                  end if
               end if
            case (couple_aggregate_standard_variables)
               if (associated(coupling%master_standard_variable)) master => generate_standard_master(root,coupling)
            case (couple_flux_sums)
               if (.not.associated(coupling%master_standard_variable)) master => generate_master(self,coupling%master_name)
         end select

         ! Save pointer to the next coupling task in advance, because current task may
         ! be deallocated from self%coupling_task_list%remove.
         next_coupling => coupling%next

         if (associated(master)) then
            ! Target variable found: perform the coupling.
            call couple_variables(root,master,coupling%slave%target)

            ! Remove coupling task from the list
            call self%coupling_task_list%remove(coupling)
         elseif (stage==couple_final) then
            call self%fatal_error('process_coupling_tasks', &
               'Coupling target "'//trim(coupling%master_name)//'" for "'//trim(coupling%slave%name)//'" was not found.')
         end if

         ! Move to next coupling task.
         coupling => next_coupling
      end do

      ! Process coupling tasks registered with child models.
      child => self%children%first
      do while (associated(child))
         call process_coupling_tasks(child%model,stage)
         child => child%next
      end do

   end subroutine process_coupling_tasks
!EOC

   subroutine create_sum(parent,link_list,name)
      class (type_base_model),intent(inout),target :: parent
      type (type_link_list),  intent(in)           :: link_list
      character(len=*),       intent(in)           :: name

      type (type_weighted_sum),pointer :: sum
      type (type_link),        pointer :: link

      allocate(sum)
      sum%result_output = output_none
      link => link_list%first
      do while (associated(link))
         call sum%add_component(link%target%name)
         link => link%next
      end do
      if (.not.sum%add_to_parent(parent,name)) deallocate(sum)
   end subroutine create_sum

   subroutine create_horizontal_sum(parent,link_list,name)
      class (type_base_model),intent(inout),target :: parent
      type (type_link_list),  intent(in)           :: link_list
      character(len=*),       intent(in)           :: name

      type (type_horizontal_weighted_sum),pointer :: sum
      type (type_link),                   pointer :: link

      allocate(sum)
      sum%result_output = output_none
      link => link_list%first
      do while (associated(link))
         call sum%add_component(link%target%name)
         link => link%next
      end do
      if (.not.sum%add_to_parent(parent,name)) deallocate(sum)
   end subroutine create_horizontal_sum

   function generate_master(self,name) result(master)
      class (type_base_model),intent(inout),target :: self
      character(len=*),       intent(in)           :: name
      type (type_internal_variable),pointer :: master

      type (type_internal_variable),pointer :: originator
      character(len=attribute_length)       :: local_name
      integer :: n

      nullify(master)
      n = len_trim(name)
      if (n<8) return
      if (name(n-7:n)=='_sms_tot'.or.name(n-7:n)=='_sfl_tot'.or.name(n-7:n)=='_bfl_tot') then
         ! This is can be a link to the sources-sinks, surface flux, or bottom flux of a state variable.
         ! First locate the related variable (same name but without the postfix), then create the necessary summation.
         originator => self%find_object(name(:n-8),recursive=.true.,exact=.false.)
         if (associated(originator)) then
            if (.not.originator%state_indices%is_empty().or.originator%fake_state_variable) then
               ! Related variable without the postfix is indeed a state variable (or acts like one)
               local_name = trim(originator%name(index(originator%name,'/',.true.)+1:))//name(n-7:n)

               ! Now we know the full path of the target variable - try to resolve it.
               ! (needed in case the derived quantity was requested for a slave)
               master => originator%owner%find_object(local_name,recursive=.false.,exact=.false.)
               if (associated(master)) return

               ! Derived quantity does not exsit yet - create it.
               select case (originator%domain)
               case (domain_interior)
                  select case (name(n-7:n))
                     case ('_sms_tot')
                        call create_sum(originator%owner,originator%sms_list,local_name)
                     case ('_sfl_tot')
                        call create_horizontal_sum(originator%owner,originator%surface_flux_list,local_name)
                     case ('_bfl_tot')
                        call create_horizontal_sum(originator%owner,originator%bottom_flux_list,local_name)
                  end select
               case (domain_horizontal,domain_surface,domain_bottom)
                  call create_horizontal_sum(originator%owner,originator%sms_list,local_name)
               end select

               ! Get pointer to newly created derived quantity
               master => originator%owner%find_object(local_name,recursive=.false.,exact=.false.)
            end if
         end if
      end if
   end function generate_master

   function generate_standard_master(self,task) result(master)
      class (type_base_model),   intent(inout),target :: self
      class (type_coupling_task),intent(inout)        :: task
      type (type_internal_variable),pointer :: master

      type (type_aggregate_variable_access), pointer :: aggregate_variable_access

      master => null()
      if (task%master_standard_variable%aggregate_variable) then
         ! Make sure that an aggregate variable will be created on the fly
         select type (aggregate_standard_variable=>task%master_standard_variable)
         class is (type_bulk_standard_variable)
            aggregate_variable_access => get_aggregate_variable_access(self,aggregate_standard_variable)
            select case (task%domain)
            case (domain_interior)
               aggregate_variable_access%interior = ior(aggregate_variable_access%interior,access_read)
               task%master_name = aggregate_standard_variable%name
            case (domain_bottom)
               aggregate_variable_access%bottom = ior(aggregate_variable_access%bottom,access_read)
               task%master_name = trim(aggregate_standard_variable%name)//'_at_bottom'
            case (domain_surface)
               aggregate_variable_access%surface = ior(aggregate_variable_access%surface,access_read)
               task%master_name = trim(aggregate_standard_variable%name)//'_at_surface'
            case default
               call self%fatal_error('generate_standard_master','BUG: unknown type of standard variable with aggregate_variable set.')
            end select
            deallocate(task%master_standard_variable)
         end select
      end if
   end function generate_standard_master

subroutine aggregate_variable_list_print(self)
   class (type_aggregate_variable_list),intent(in) :: self

   type (type_aggregate_variable),   pointer :: aggregate_variable
   type (type_contributing_variable),pointer :: contributing_variable

   aggregate_variable => self%first
   do while (associated(aggregate_variable))
      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (associated(contributing_variable%link%target,contributing_variable%link%original)) then
            call log_message('   '//trim(contributing_variable%link%name)//' (internal)')
         else
            call log_message('   '//trim(contributing_variable%link%name)//' (external)')
         end if
         contributing_variable => contributing_variable%next
      end do
      aggregate_variable => aggregate_variable%next
   end do
end subroutine aggregate_variable_list_print

function aggregate_variable_list_get(self,standard_variable) result(aggregate_variable)
   class (type_aggregate_variable_list), intent(inout)      :: self
   type (type_bulk_standard_variable),   intent(in), target :: standard_variable

   type (type_aggregate_variable), pointer :: aggregate_variable

   aggregate_variable => self%first
   do while (associated(aggregate_variable))
      if (aggregate_variable%standard_variable%compare(standard_variable)) return
      aggregate_variable => aggregate_variable%next
   end do
   allocate(aggregate_variable)
   aggregate_variable%standard_variable => standard_variable
   aggregate_variable%next => self%first
   self%first => aggregate_variable
end function

function collect_aggregate_variables(self) result(list)
   class (type_base_model),intent(in),target :: self
   type (type_aggregate_variable_list)       :: list

   type (type_link),                 pointer :: link
   type (type_contribution),         pointer :: contribution
   type (type_contributing_variable),pointer :: contributing_variable
   type (type_aggregate_variable),   pointer :: aggregate_variable

   link => self%links%first
   do while (associated(link))
      ! Enumerate the contributions of this variable to aggregate variables, and register these with
      ! aggregate variable objects on the model level.
      contribution => link%target%contributions%first
      do while (associated(contribution))
         aggregate_variable => list%get(contribution%target)
         allocate(contributing_variable)
         contributing_variable%link => link
         contributing_variable%scale_factor = contribution%scale_factor
         contributing_variable%include_background = contribution%include_background
         contributing_variable%next => aggregate_variable%first_contributing_variable
         aggregate_variable%first_contributing_variable => contributing_variable
         contribution => contribution%next
      end do
      link => link%next
   end do
end function collect_aggregate_variables

recursive subroutine create_aggregate_models(self)
   class (type_base_model),intent(inout),target :: self

   type (type_aggregate_variable_access), pointer :: aggregate_variable_access
   type (type_aggregate_variable),        pointer :: aggregate_variable
   class (type_weighted_sum),             pointer :: sum
   class (type_horizontal_weighted_sum),  pointer :: horizontal_sum,bottom_sum,surface_sum
   type (type_contributing_variable),     pointer :: contributing_variable
   type (type_aggregate_variable_list)              :: list
   type (type_model_list_node),           pointer :: child

   ! Get a list of all aggregate
   list = collect_aggregate_variables(self)

   ! Make sure that readable fields for all aggregate variables are created at the root level.
   if (.not.associated(self%parent)) then
      aggregate_variable => list%first
      do while (associated(aggregate_variable))
         aggregate_variable_access => get_aggregate_variable_access(self,aggregate_variable%standard_variable)
         aggregate_variable_access%interior = ior(aggregate_variable_access%interior,access_read)
         aggregate_variable_access%horizontal = ior(aggregate_variable_access%horizontal,access_read)
         aggregate_variable => aggregate_variable%next
      end do
   end if

   aggregate_variable_access => self%first_aggregate_variable_access
   do while (associated(aggregate_variable_access))
      nullify(sum,horizontal_sum,bottom_sum,surface_sum)
      if (aggregate_variable_access%interior  /=access_none) allocate(sum)
      if (aggregate_variable_access%horizontal/=access_none) allocate(horizontal_sum)
      if (aggregate_variable_access%bottom    /=access_none) allocate(bottom_sum)
      if (aggregate_variable_access%surface   /=access_none) allocate(surface_sum)

      aggregate_variable => list%get(aggregate_variable_access%standard_variable)
      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (associated(contributing_variable%link%target,contributing_variable%link%original) &                  ! Variable must not be coupled
             .and.(associated(self%parent).or..not.contributing_variable%link%target%fake_state_variable) &       ! Only include fake state variable for non-root models
             .and.(index(contributing_variable%link%name,'/')==0.or..not.associated(self%parent))) then           ! Variable must be owned by the model itself unless we are aggregating at root level
             if (associated(sum)) then
               select case (contributing_variable%link%target%domain)
               case (domain_interior)
                  call sum%add_component(trim(contributing_variable%link%name), &
                     weight=contributing_variable%scale_factor,include_background=contributing_variable%include_background)
               end select
            end if
            if (associated(bottom_sum)) then
               select case (contributing_variable%link%target%domain)
               case (domain_bottom)
                  call bottom_sum%add_component(trim(contributing_variable%link%name), &
                     weight=contributing_variable%scale_factor,include_background=contributing_variable%include_background)
               end select
            end if
            if (associated(surface_sum)) then
               select case (contributing_variable%link%target%domain)
               case (domain_surface)
                  call surface_sum%add_component(trim(contributing_variable%link%name), &
                     weight=contributing_variable%scale_factor,include_background=contributing_variable%include_background)
               end select
            end if
            if (associated(horizontal_sum)) then
               select case (contributing_variable%link%target%domain)
               case (domain_horizontal,domain_surface,domain_bottom)
                  call horizontal_sum%add_component(trim(contributing_variable%link%name), &
                     weight=contributing_variable%scale_factor,include_background=contributing_variable%include_background)
               end select
            end if
         end if
         contributing_variable => contributing_variable%next
      end do

      if (associated(sum)) then
         sum%units = trim(aggregate_variable%standard_variable%units)
         sum%access = aggregate_variable_access%interior
         if (associated(self%parent)) then
            sum%result_output = output_none
         else
            allocate(sum%standard_variable)
            sum%standard_variable = aggregate_variable%standard_variable
         end if
         if (.not.sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name),aggregate_variable=aggregate_variable%standard_variable)) deallocate(sum)
      end if
      if (associated(horizontal_sum)) then
         horizontal_sum%units = trim(aggregate_variable%standard_variable%units)//'*m'
         horizontal_sum%access = aggregate_variable_access%horizontal
         if (associated(self%parent)) horizontal_sum%result_output = output_none
         if (.not.horizontal_sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name)//'_at_interfaces')) deallocate(horizontal_sum)
      end if
      if (associated(bottom_sum)) then
         bottom_sum%units = trim(aggregate_variable%standard_variable%units)//'*m'
         bottom_sum%access = aggregate_variable_access%bottom
         bottom_sum%domain = domain_bottom
         if (associated(self%parent)) bottom_sum%result_output = output_none
         if (.not.bottom_sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name)//'_at_bottom',aggregate_variable=aggregate_variable%standard_variable)) deallocate(bottom_sum)
      end if
      if (associated(surface_sum)) then
         surface_sum%units = trim(aggregate_variable%standard_variable%units)//'*m'
         surface_sum%access = aggregate_variable_access%surface
         surface_sum%domain = domain_surface
         if (associated(self%parent)) surface_sum%result_output = output_none
         if (.not.surface_sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name)//'_at_surface',aggregate_variable=aggregate_variable%standard_variable)) deallocate(surface_sum)
      end if
      aggregate_variable_access => aggregate_variable_access%next
   end do

   ! Process child models
   child => self%children%first
   do while (associated(child))
      call create_aggregate_models(child%model)
      child => child%next
   end do
end subroutine create_aggregate_models

recursive subroutine create_conservation_checks(self)
   class (type_base_model),intent(inout),target :: self

   type (type_aggregate_variable_list)          :: aggregate_variable_list
   type (type_aggregate_variable),      pointer :: aggregate_variable
   class (type_weighted_sum),           pointer :: sum
   class (type_horizontal_weighted_sum),pointer :: surface_sum,bottom_sum
   type (type_contributing_variable),   pointer :: contributing_variable
   type (type_model_list_node),         pointer :: child

   if (self%check_conservation) then
      aggregate_variable_list = collect_aggregate_variables(self)

      aggregate_variable => aggregate_variable_list%first
      do while (associated(aggregate_variable))
         if (aggregate_variable%standard_variable%conserved) then
            ! Allocate objects that will do the summation across the different domains.
            allocate(sum,surface_sum,bottom_sum)

            ! Enumerate contributions to aggregate variable.
            contributing_variable => aggregate_variable%first_contributing_variable
            do while (associated(contributing_variable))
               if (index(contributing_variable%link%name,'/')==0) then
                  select case (contributing_variable%link%original%domain)
                     case (domain_interior)
                        if (associated(contributing_variable%link%original%sms))          call sum%add_component        (contributing_variable%link%original%sms%name,         contributing_variable%scale_factor)
                        if (associated(contributing_variable%link%original%surface_flux)) call surface_sum%add_component(contributing_variable%link%original%surface_flux%name,contributing_variable%scale_factor)
                        if (associated(contributing_variable%link%original%bottom_flux))  call bottom_sum%add_component (contributing_variable%link%original%bottom_flux%name, contributing_variable%scale_factor)
                     case (domain_surface)
                        if (associated(contributing_variable%link%original%sms)) call surface_sum%add_component(contributing_variable%link%original%sms%name,contributing_variable%scale_factor)
                     case (domain_bottom)
                        if (associated(contributing_variable%link%original%sms)) call bottom_sum%add_component(contributing_variable%link%original%sms%name,contributing_variable%scale_factor)
                  end select
               end if   
               contributing_variable => contributing_variable%next
            end do

            ! Process sums now that all contributing terms are known.
            sum%units = trim(aggregate_variable%standard_variable%units)//'/s'
            if (.not.sum%add_to_parent(self,'change_in_'//trim(aggregate_variable%standard_variable%name),create_for_one=.true.)) deallocate(sum)
            surface_sum%units = trim(aggregate_variable%standard_variable%units)//'*m/s'
            if (.not.surface_sum%add_to_parent(self,'change_in_'//trim(aggregate_variable%standard_variable%name)//'_at_surface',create_for_one=.true.)) deallocate(surface_sum)
            bottom_sum%units = trim(aggregate_variable%standard_variable%units)//'*m/s'
            if (.not.bottom_sum%add_to_parent(self,'change_in_'//trim(aggregate_variable%standard_variable%name)//'_at_bottom',create_for_one=.true.)) deallocate(bottom_sum)
         end if
         aggregate_variable => aggregate_variable%next
      end do
   end if

   ! Process child models
   child => self%children%first
   do while (associated(child))
      call create_conservation_checks(child%model)
      child => child%next
   end do

end subroutine create_conservation_checks

recursive subroutine couple_variables(self,master,slave)
   class (type_base_model),     intent(inout),target :: self
   type (type_internal_variable),pointer             :: master,slave

   type (type_internal_variable),pointer :: pslave
   logical                               :: slave_is_state_variable
   logical                               :: master_is_state_variable

   ! If slave and master are the same, we are done - return.
   if (associated(slave,master)) return

   slave_is_state_variable  = slave%fake_state_variable.or..not.slave%state_indices%is_empty()
   master_is_state_variable = master%fake_state_variable.or..not.master%state_indices%is_empty()

   if (associated(self%parent)) call self%fatal_error('couple_variables','BUG: must be called on root node.')
   if (.not.slave%can_be_slave) &
      call fatal_error('couple_variables','Attempt to couple write-only variable ' &
         //trim(slave%name)//' to '//trim(master%name)//'.')
   if (slave_is_state_variable) then
      ! Extra checks when coupling state variables
      if (.not.master_is_state_variable) call fatal_error('couple_variables','Attempt to couple state variable ' &
         //trim(slave%name)//' to non-state variable '//trim(master%name)//'.')
      if ((slave%domain==domain_bottom.and.master%domain==domain_surface).or.(slave%domain==domain_surface.and.master%domain==domain_bottom)) &
         call fatal_error('couple_variables', 'Cannot couple '//trim(slave%name)//' to '//trim(master%name)//', because their domains are incompatible.')
   end if
   !if (slave%domain/=master%domain.and..not.(slave%domain==domain_horizontal.and. &
   !   (master%domain==domain_surface.or.master%domain==domain_bottom))) call fatal_error('couple_variables', &
   !   'Cannot couple '//trim(slave%name)//' to '//trim(master%name)//', because their domains are incompatible.')
   if (iand(slave%domain,master%domain)==0) call fatal_error('couple_variables', &
      'Cannot couple '//trim(slave%name)//' to '//trim(master%name)//', because their domains are incompatible.')

   if (debug_coupling) call log_message(trim(slave%name)//' --> '//trim(master%name))

   ! Merge all information from the slave into the master.
   call master%state_indices%extend(slave%state_indices)
   call master%read_indices%extend(slave%read_indices)
   call master%write_indices%extend(slave%write_indices)
   call master%sms_list%extend(slave%sms_list)
   call master%background_values%extend(slave%background_values)
   call master%properties%update(slave%properties,overwrite=.false.)
   call master%surface_flux_list%extend(slave%surface_flux_list)
   call master%bottom_flux_list%extend(slave%bottom_flux_list)

   call master%standard_variables%update(slave%standard_variables)

   ! For vertical movement rates only keep the master, which all models will (over)write.
   ! NB if the slave has vertical movement but the master does not (e.g., if the master is
   ! a fake state variable, the slave variable can still be set, but won't be used).
   if (associated(slave%movement_diagnostic).and.associated(master%movement_diagnostic)) &
      call couple_variables(self,master%movement_diagnostic%target,slave%movement_diagnostic%target)
   if (master%presence==presence_external_optional.and.slave%presence/=presence_external_optional) &
      master%presence = presence_external_required

   ! Store a pointer to the slave, because the call to redirect_links will cause all pointers (from links)
   ! to the slave node to be connected to the master node. This includes the original "slave" argument.
   pslave => slave

   ! Note: in call below we provide our local copy of the pointer to the slave (pslave), not the original slave pointer (slave).
   ! Reason: ifort appears to pass the slave pointer by reference. The original pointer comes from a link, and is therefore
   ! overwritten by redirect_links. If we pass this original pointer to redirect_links for comparing object identities,
   ! the recursive redirecting will fail after the very first redirect (which destroyed the pointer to the object that we
   ! want to compare against).
   call redirect_links(self,pslave,master)
end subroutine couple_variables

recursive subroutine redirect_links(model,oldtarget,newtarget)
   class (type_base_model),     intent(inout),target :: model
   type (type_internal_variable),pointer :: oldtarget,newtarget

   type (type_link),           pointer :: link
   type (type_model_list_node),pointer :: child

   ! Process all links and if they used to refer to the specified slave,
   ! redirect them to the specified master.
   link => model%links%first
   do while (associated(link))
      if (associated(link%target,oldtarget)) link%target => newtarget
      link => link%next
   end do

   ! Allow child models to do the same.
   child => model%children%first
   do while (associated(child))
      call redirect_links(child%model,oldtarget,newtarget)
      child => child%next
   end do
end subroutine

recursive subroutine find_dependencies(self,list,forbidden)
   class (type_base_model),intent(in),target   :: self
   type (type_model_list), intent(inout)       :: list
   type (type_model_list), intent(in),optional :: forbidden

   type (type_link),pointer            :: link
   type (type_model_list)              :: forbidden_with_self
   type (type_model_list_node),pointer :: node
   character(len=2048)                 :: chain

   ! Directly below and further down: use of generic procedure list%find confuses PGI 18.10 (Jorn 2019-04-24)
   if (associated(list%find_model(self))) return

   ! Check the list of forbidden model (i.e., models that indirectly request the current model)
   ! If the current model is on this list, there is a circular dependency between models.
   if (present(forbidden)) then
      node => forbidden%find_model(self)
      if (associated(node)) then
         ! Circular dependency found - report as fatal error.
         chain = ''
         do while (associated(node))
            chain = trim(chain)//' '//trim(node%model%get_path())//' ->'
            node => node%next
         end do
         call fatal_error('find_dependencies','circular dependency found: '//trim(chain(2:))//' '//trim(self%get_path()))
         return
      end if
      call forbidden_with_self%extend(forbidden)
   end if
   call forbidden_with_self%append(self)

   ! Loop over all variables, and if they belong to some other model, first add that model to the dependency list.
   link => self%links%first
   do while (associated(link))
      if (index(link%name,'/')==0 &                        ! Our own link...
          .and..not.link%target%write_indices%is_empty() & ! ...to a diagnostic variable...
          .and.link%target%source/=source_none &           ! ...not a constant...
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

   recursive subroutine check_coupling_units(self)
      class (type_base_model),intent(in),target :: self

      type (type_link),           pointer :: link
      type (type_model_list_node),pointer :: child

      link => self%links%first
      do while (associated(link))
         if (index(link%name,'/')==0 .and. .not. associated(link%target,link%original)) then
            if (link%target%units/=''.and. link%original%units/=''.and. link%target%units/=link%original%units) &
               call log_message('WARNING: unit mismatch between master '//trim(link%target%name)//' ('//trim(link%target%units)// &
                  ') and slave '//trim(link%original%name)//' ('//trim(link%original%units)//').')
         end if
         link => link%next
      end do

      ! Process child models
      child => self%children%first
      do while (associated(child))
         call check_coupling_units(child%model)
         child => child%next
      end do
   end subroutine check_coupling_units

function call_list_find(self,model,source) result(node)
   class (type_call_list), intent(in)        :: self
   class (type_base_model),intent(in),target :: model
   integer,                intent(in)        :: source

   type (type_call_list_node),pointer :: node

   node => self%first
   do while (associated(node))
      if (associated(node%model,model).and.node%source==source) return
      node => node%next
   end do
end function call_list_find

subroutine call_list_extend(self,other)
   class (type_call_list), intent(inout) :: self
   class (type_call_list), intent(in)    :: other

   type (type_call_list_node),pointer :: last, node2

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

function call_list_append(self,model,source) result(node)
   class (type_call_list), intent(inout)     :: self
   class (type_base_model),intent(in),target :: model
   integer,                intent(in)        :: source

   type (type_call_list_node),pointer :: node

   if (associated(self%first)) then
      ! List contains one or more items - append to tail.
      node => self%first
      do while (associated(node%next))
         node => node%next
      end do
      allocate(node%next)
      node => node%next
   else
      ! List is empty - create first node.
      allocate(self%first)
      node => self%first
   end if

   ! Set new node
   node%model => model
   node%source = source
end function call_list_append

subroutine call_list_filter(self,source)
   class (type_call_list), intent(inout) :: self
   integer,                intent(in)    :: source

   type (type_call_list_node),pointer :: node, previous, next

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

   type (type_call_list_node),   pointer :: node
   type (type_variable_set_node),pointer :: variable_list_node
   integer                               :: n

   node => self%first
   do while (associated(node))
      write (*,'(a,x,i0)') trim(node%model%get_path()),node%source
      n = 0
      variable_list_node => node%written_variables%first
      do while (associated(variable_list_node))
         n = n + 1
         write (*,'(x,x,x,a,x,i0)') trim(variable_list_node%target%name)
         variable_list_node => variable_list_node%next
      end do
      node => node%next
   end do
end subroutine call_list_print

subroutine call_list_finalize(self)
   class (type_call_list), intent(inout) :: self

   type (type_call_list_node),pointer :: node,next

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

   type (type_call_list_node),pointer :: node

   call_list_computes = .true.
   node => self%first
   do while (associated(node))
      if (node%computed_variables%contains(variable)) return
      node => node%next
   end do
   call_list_computes = .false.
end function call_list_computes

subroutine call_list_node_finalize(self)
   class (type_call_list_node), intent(inout) :: self

   if (allocated(self%copy_commands_int)) deallocate(self%copy_commands_int)
   if (allocated(self%copy_commands_hz)) deallocate(self%copy_commands_hz)
   call self%written_variables%finalize()
   call self%computed_variables%finalize()
end subroutine call_list_node_finalize

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

recursive subroutine find_dependencies2(self,source,allowed_sources,list,forbidden,call_list_node)
   class (type_base_model),     intent(in),target   :: self
   integer,                     intent(in)          :: source
   integer,                     intent(in)          :: allowed_sources(:)
   type (type_call_list),       intent(inout)       :: list
   type (type_call_list),       intent(in),optional :: forbidden
   type (type_call_list_node),  pointer,   optional :: call_list_node

   type (type_link),pointer           :: link
   type (type_call_list)              :: forbidden_with_self
   type (type_call_list_node),pointer :: node
   character(len=2048)                :: chain
   logical                            :: same_source

   if (present(call_list_node)) nullify(call_list_node)

   ! Do nothing if we are looking for constants.
   if (source==source_none) return

   ! Check whether we are already processed this call.
   node => list%find(self,source)
   if (associated(node)) then
      if (present(call_list_node)) call_list_node => node
      return
   end if

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
         return
      end if
      call forbidden_with_self%extend(forbidden)
   end if
   node => forbidden_with_self%append(self,source)

   ! Loop over all variables, and if they belong to some other model, first add that model to the dependency list.
   link => self%links%first
   do while (associated(link))
      same_source = link%target%source==source .or. (link%target%source==source_unknown.and.(source==source_do_surface.or.source==source_do_bottom))
      if (index(link%name,'/')==0 &                                        ! Our own link...
          .and..not.link%target%write_indices%is_empty() &                 ! ...to a diagnostic variable...
          .and.link%target%source/=source_none &                           ! ...not a constant...
          .and..not.(associated(link%target%owner,self).and.same_source) & ! ...not set by ourselves...
          .and.associated(link%original%read_index)) then                  ! ...and we do depend on its value

         if (contains_value(allowed_sources,link%target%source)) &
            call find_variable_dependencies(link%target,allowed_sources,list,copy_to_prefetch=.true.,forbidden=forbidden_with_self)
      end if
      link => link%next
   end do

   ! We're happy - add ourselves to the list of processed models.
   node => list%append(self,source)
   if (present(call_list_node)) call_list_node => node

   ! Clean up our temporary list.
   call forbidden_with_self%finalize()
end subroutine find_dependencies2

recursive subroutine find_variable_dependencies(variable,allowed_sources,list,copy_to_prefetch,forbidden)
   type (type_internal_variable),intent(in)          :: variable
   integer,                      intent(in)          :: allowed_sources(:)
   type (type_call_list),        intent(inout)       :: list
   logical,                      intent(in)          :: copy_to_prefetch
   type (type_call_list),        intent(in),optional :: forbidden

   type (type_call_list_node),pointer :: node

   if (variable%source==source_unknown) then
      if (contains_value(allowed_sources,source_do_surface)) then
         call find_dependencies2(variable%owner,source_do_surface,allowed_sources,list,forbidden,call_list_node=node)
         if (copy_to_prefetch.and.associated(node)) call node%written_variables%add(variable)
         if (associated(node)) call node%computed_variables%add(variable)
      end if
      if (contains_value(allowed_sources,source_do_bottom)) then
         call find_dependencies2(variable%owner,source_do_bottom,allowed_sources,list,forbidden,call_list_node=node)
         if (copy_to_prefetch.and.associated(node)) call node%written_variables%add(variable)
         if (associated(node)) call node%computed_variables%add(variable)
      end if
   else
      call find_dependencies2(variable%owner,variable%source,allowed_sources,list,forbidden,call_list_node=node)
      if (copy_to_prefetch.and.associated(node)) call node%written_variables%add(variable)
      if (associated(node)) call node%computed_variables%add(variable)
   end if
end subroutine

logical function contains_value(array,value)
   integer,intent(in) :: array(:),value
   integer :: i
   contains_value = .true.
   do i=1,size(array)
      if (array(i)==value) return
   end do
   contains_value = .false.
end function

subroutine variable_set_add(self,variable)
   class (type_variable_set),intent(inout) :: self
   type (type_internal_variable),target    :: variable

   type (type_variable_set_node), pointer :: node

   node => self%first
   do while (associated(node))
      if (associated(node%target,variable)) return
      node => node%next
   end do
   allocate(node)
   node%target => variable
   node%next => self%first
   self%first => node
end subroutine variable_set_add

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

subroutine call_list_initialize(call_list)
   class (type_call_list),intent(inout) :: call_list

   type (type_call_list_node),pointer :: node

   node => call_list%first
   do while (associated(node))
      call call_list_node_initialize(node)
      node => node%next
   end do
end subroutine call_list_initialize

subroutine call_list_node_initialize(call_list_node)
   type (type_call_list_node),intent(inout) :: call_list_node

   class (type_base_model),      pointer :: parent
   class (type_model_list_node), pointer :: model_list_node

   ! Make sure the pointer to the model has the highest class (and not a base class)
   ! This is needed because model classes that use inheritance and call base class methods
   ! may end up with model pointers that are base class specific (and do not reference
   ! procedures overwritten at a higher level)
   parent => call_list_node%model%parent
   if (associated(parent)) then
      model_list_node => parent%children%first
      do while (associated(model_list_node))
         if (associated(call_list_node%model,model_list_node%model)) then
            ! Found ourselves in our parent - use the parent pointer to replace ours.
            call_list_node%model => model_list_node%model
            exit
         end if
         model_list_node => model_list_node%next
      end do
   end if

   call process_domain(call_list_node%copy_commands_int, domain_interior)
   call process_domain(call_list_node%copy_commands_hz,  domain_horizontal)

contains

   subroutine process_domain(commands,domain)
      type (type_copy_command), allocatable, intent(inout) :: commands(:)
      integer,                               intent(in)    :: domain

      type (type_variable_set_node), pointer :: node
      integer :: i, n, maxwrite

      n = 0
      maxwrite = -1
      node => call_list_node%written_variables%first
      do while (associated(node))
         if (node%target%write_indices%is_empty()) call fatal_error('call_list_node_initialize','BUG: target without write indices')
         if (iand(node%target%domain,domain)/=0) then
            n = n + 1
            maxwrite = max(maxwrite,node%target%write_indices%value)
         end if
         node => node%next
      end do

      ! Create list of copy commands, sorted by write index
      allocate(commands(n))
      n = 0
      do i=1,maxwrite
         node => call_list_node%written_variables%first
         do while (associated(node))
            if (iand(node%target%domain,domain)/=0.and.node%target%write_indices%value==i) then
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

end module

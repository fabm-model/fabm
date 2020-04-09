module fabm_coupling
   use fabm_types
   use fabm_standard_variables
   use fabm_properties
   use fabm_builtin_models
   use fabm_driver

   implicit none

   private

   public freeze_model_info
   public collect_aggregate_variables, type_aggregate_variable_list, type_aggregate_variable

   logical,parameter :: debug_coupling = .false.

   integer,parameter :: couple_explicit                     = 0
   integer,parameter :: couple_aggregate_standard_variables = 1
   integer,parameter :: couple_final                        = 3

   type type_contributing_variable
      type (type_link),                 pointer :: link               => null()
      real(rk)                                  :: scale_factor       = 1.0_rk
      logical                                   :: include_background = .false.
      type (type_contributing_variable),pointer :: next               => null()
   end type

   type type_aggregate_variable
      class (type_domain_specific_standard_variable), pointer :: standard_variable           => null()
      type (type_contributing_variable),              pointer :: first_contributing_variable => null()
      type (type_aggregate_variable),                 pointer :: next                        => null()
   end type

   type type_aggregate_variable_list
      type (type_aggregate_variable), pointer :: first => null()
   contains
      procedure :: get      => aggregate_variable_list_get
      procedure :: print    => aggregate_variable_list_print
      procedure :: finalize => aggregate_variable_list_finalize
   end type

contains

   subroutine freeze_model_info(self)
      class (type_base_model), intent(inout), target :: self

      if (associated(self%parent)) call self%fatal_error('freeze_model_info', &
         'BUG: freeze_model_info can only operate on the root model.')

      call before_coupling(self)

      ! Coupling stage 1: implicit - couple variables based on overlapping standard identities.
      call couple_standard_variables(self)

      ! Coupling stage 2: explicit - resolve user- or model-specified links between variables.
      call process_coupling_tasks(self, couple_explicit)

      ! Coupling stage 3: try to automatically fulfil remaining dependencies on aggregate standard variables
      ! This may create new child models to handle the necessary summations.
      ! Note that these child models may also add source terms (if the sum is treated as state variable),
      ! so any source term treatment should happen after this is complete!
      call process_coupling_tasks(self, couple_aggregate_standard_variables)

      ! Create models for aggregate variables at root level, to be used to compute conserved quantities.
      ! After this step, the set of variables that contribute to aggregate quantities may not be modified.
      ! That is, no new such variables may be added, and no such variables may be coupled.
      call create_aggregate_models(self)

      ! Perform coupling for any new aggregate models.
      ! This may append items to existing lists of source terms and bottom/surface fluxes,
      ! so it has to precede the call to create_flux_sums.
      call process_coupling_tasks(self, couple_explicit)

      ! Now that coupling for non-rate variables is complete, contributions to aggregate quantities
      ! (including ones of slave variables) are final.
      ! Create conservation checks where needed.
      call create_conservation_checks(self)

      ! Create summations of source terms and surface/bottom fluxes.
      call create_flux_sums(self)

      ! Process the any remaining coupling tasks.
      call process_coupling_tasks(self, couple_final)

      ! Allow inheriting models to perform additional tasks after coupling.
      call after_coupling(self)

      ! Check whether units of coupled variables match
      !call check_coupling_units(self)

      call freeze(self)
   end subroutine freeze_model_info

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
      class (type_base_model), intent(inout), target :: self

      type (type_model_list_node), pointer :: node

      call self%after_coupling()
      node => self%children%first
      do while (associated(node))
         call after_coupling(node%model)
         node => node%next
      end do
   end subroutine

   recursive subroutine freeze(self)
      class (type_base_model), intent(inout), target :: self

      type (type_model_list_node), pointer :: node

      self%frozen = .true.
      node => self%children%first
      do while (associated(node))
         call freeze(node%model)
         node => node%next
      end do
   end subroutine

   subroutine couple_standard_variables(model)
      class (type_base_model), intent(inout), target :: model

      type (type_link),                   pointer :: link, first_link
      type (type_standard_variable_set)           :: all_standard_variables
      type (type_standard_variable_node), pointer :: node

      ! Build a list of all unique standard variables.
      link => model%links%first
      do while (associated(link))
         call all_standard_variables%update(link%target%standard_variables)
         link => link%next
      end do

      ! Loop over all unique standard variable and collect and couple associated model variables.
      node => all_standard_variables%first
      do while (associated(node))
#ifndef NDEBUG
         call node%p%assert_resolved()
#endif
         first_link => null()
         link => model%links%first
         do while (associated(link))
            if (link%target%standard_variables%contains(node%p)) then
               if (associated(first_link)) then
                  if (link%target%write_indices%is_empty() .and. .not. (first_link%target%presence /= presence_internal .and. link%target%presence == presence_internal)) then
                     ! Default coupling: early variable (first_link) is master, later variable (link) is slave.
                     call couple_variables(model, first_link%target, link%target)
                  else
                     ! Later variable (link) is write-only and therefore can only be master. Try coupling with early variable (first_link) as slave.
                     call couple_variables(model, link%target, first_link%target)
                  end if
               else
                  first_link => link
               end if
            end if
            link => link%next
         end do
         node => node%next
      end do
      call all_standard_variables%finalize()
   end subroutine couple_standard_variables

   subroutine collect_user_specified_couplings(self)
      class (type_base_model), intent(inout) :: self

      type (type_link),           pointer :: link
      class (type_property),      pointer :: master_name
      class (type_coupling_task), pointer :: task

      link => self%links%first
      do while (associated(link))
         ! Only process own links (those without slash in the name)
         if (index(link%name, '/') == 0) then
            master_name => self%couplings%find_in_tree(link%name)
            if (associated(master_name)) then
               call self%coupling_task_list%add(link, .true., task)
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

   recursive subroutine process_coupling_tasks(self, stage)
      class (type_base_model), intent(inout), target :: self
      integer,                 intent(in)            :: stage

      class (type_base_model),       pointer :: root
      type (type_model_list_node),   pointer :: child
      class (type_coupling_task),    pointer :: coupling, next_coupling
      type (type_internal_variable), pointer :: master
      type (type_link),              pointer :: link

      ! Find root model, which will handle the individual coupling tasks.
      root => self
      do while (associated(root%parent))
         root => root%parent
      end do

      ! Gather all couplings that were specified by the user as part of the FABM configuration (fabm.yaml)
      ! This is done only once, as after collection the includes_custom flag will be set to .true.
      if (.not. self%coupling_task_list%includes_custom) call collect_user_specified_couplings(self)

      ! For each variable, determine if a coupling command is provided.
      coupling => self%coupling_task_list%first
      do while (associated(coupling))

         master => null()
         select case (stage)
            case (couple_explicit, couple_final)
               if (associated(coupling%master_standard_variable)) then
                  ! This is a coupling to a standard variable. First try to find the corresponding standard variable.
                  ! We search within the root model, because there all variables are found together.
                  link => root%links%first
                  do while (associated(link))
                     if (link%target%standard_variables%contains(coupling%master_standard_variable)) exit
                     link => link%next
                  end do

                  if (stage == couple_final .and. .not. associated(link)) then
                     ! Target variable was not found, but this is our last chance.
                     ! Therefore, create a placeholder variable at the root level.
                     ! This variable will still need to be provided by the host.
                     select type (standard_variable => coupling%master_standard_variable)
                     class is (type_interior_standard_variable)
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

                  ! Store link to master variable
                  if (associated(link)) master => link%target
               else
                  ! This is a coupling by variable name.
                  ! Try to find the master variable among the variables of the requesting model or its parents.
                  if (coupling%slave%name /= coupling%master_name) then
                     ! Master and slave name differ: start master search in current model, then move up tree.
                     master => self%find_object(coupling%master_name, recursive=.true., exact=.false.)
                  elseif (associated(self%parent)) then
                     ! Master and slave name are identical: start master search in parent model, then move up tree.
                     master => self%parent%find_object(coupling%master_name, recursive=.true., exact=.false.)
                  else
                     call self%fatal_error('process_coupling_tasks', &
                        'Master and slave name are identical: "' // trim(coupling%master_name) // '". This is not valid at the root of the model tree.')
                  end if
               end if
            case (couple_aggregate_standard_variables)
               if (associated(coupling%master_standard_variable)) master => generate_standard_master(root, coupling)
         end select

         ! Save pointer to the next coupling task in advance, because current task may
         ! be deallocated from self%coupling_task_list%remove.
         next_coupling => coupling%next

         if (associated(master)) then
            ! Target variable found: perform the coupling.
            call couple_variables(root, master, coupling%slave%target)

            ! Remove coupling task from the list
            call self%coupling_task_list%remove(coupling)
         elseif (stage == couple_final) then
            call self%fatal_error('process_coupling_tasks', &
               'Coupling target "' // trim(coupling%master_name) // '" for "' // trim(coupling%slave%name) // '" was not found.')
         end if

         ! Move to next coupling task.
         coupling => next_coupling
      end do

      ! Process coupling tasks registered with child models.
      child => self%children%first
      do while (associated(child))
         call process_coupling_tasks(child%model, stage)
         child => child%next
      end do

   end subroutine process_coupling_tasks

   function create_sum(parent, link_list, name) result(link)
      class (type_base_model), intent(inout), target :: parent
      type (type_link_list),   intent(in)            :: link_list
      character(len=*),        intent(in)            :: name
      type (type_link), pointer                      :: link

      type (type_weighted_sum), pointer :: sum

      allocate(sum)
      sum%missing_value = 0
      sum%result_output = output_none
      link => link_list%first
      do while (associated(link))
         call sum%add_component(link%target%name)
         link => link%next
      end do
      if (.not. sum%add_to_parent(parent, name, link=link)) deallocate(sum)
   end function create_sum

   function create_horizontal_sum(parent, link_list, name) result(link)
      class (type_base_model), intent(inout), target :: parent
      type (type_link_list),   intent(in)            :: link_list
      character(len=*),        intent(in)            :: name
      type (type_link), pointer                      :: link

      type (type_horizontal_weighted_sum), pointer :: sum

      allocate(sum)
      sum%missing_value = 0
      sum%result_output = output_none
      link => link_list%first
      do while (associated(link))
         call sum%add_component(link%target%name)
         link => link%next
      end do
      if (.not. sum%add_to_parent(parent, name, link=link)) deallocate(sum)
   end function create_horizontal_sum

   recursive subroutine create_flux_sums(self)
      class (type_base_model), intent(inout), target :: self

      type (type_link),            pointer :: link, link1
      type (type_model_list_node), pointer :: child

      link => self%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. (link%original%source == source_state .or. link%original%fake_state_variable)) then
            ! This is a state variable, or a diagnostic pretending to be one, that we have registered (it is owned by "self")
            if (associated(link%target, link%original)) then
               ! We own this variable (it has not been coupled to another). Create summations for sources-sinks and surface/bottom fluxes.
               select case (link%target%domain)
               case (domain_interior)
                  link%target%sms_sum          => create_sum(self, link%target%sms_list,                     trim(link%name) // '_sms_tot')
                  link%target%surface_flux_sum => create_horizontal_sum(self, link%target%surface_flux_list, trim(link%name) // '_sfl_tot')
                  link%target%bottom_flux_sum  => create_horizontal_sum(self, link%target%bottom_flux_list,  trim(link%name) // '_bfl_tot')
                  link%target%movement_sum     => create_sum(self, link%target%movement_list,                trim(link%name) // '_w_tot')
               case (domain_horizontal, domain_surface, domain_bottom)
                  link%target%sms_sum          => create_horizontal_sum(self, link%target%sms_list,         trim(link%name) // '_sms_tot')
               end select
            else
               ! We do not own this variable. Link to summations for sources-sinks and surface/bottom fluxes.
               select case (link%target%domain)
               case (domain_interior)
                  link1 => null()
                  call self%add_interior_variable(trim(link%name) // '_sms_tot', link=link1)
                  call self%request_coupling(link1, trim(link%target%name) // '_sms_tot')
                  link1 => null()
                  call self%add_horizontal_variable(trim(link%name) // '_sfl_tot', link=link1)
                  call self%request_coupling(link1, trim(link%target%name) // '_sfl_tot')
                  link1 => null()
                  call self%add_horizontal_variable(trim(link%name) // '_bfl_tot', link=link1)
                  call self%request_coupling(link1, trim(link%target%name) // '_bfl_tot')
                  link1 => null()
                  call self%add_interior_variable(trim(link%name) // '_w_tot', link=link1)
                  call self%request_coupling(link1, trim(link%target%name) // '_w_tot')
               case (domain_horizontal, domain_surface, domain_bottom)
                  link1 => null()
                  call self%add_horizontal_variable(trim(link%name) // '_sms_tot', link=link1)
                  call self%request_coupling(link1, trim(link%target%name) // '_sms_tot')
               end select
            end if
         end if
         link => link%next
      end do

      ! Process child models
      child => self%children%first
      do while (associated(child))
         call create_flux_sums(child%model)
         child => child%next
      end do
   end subroutine create_flux_sums

   function generate_standard_master(self, task) result(master)
      class (type_base_model),    intent(inout), target :: self
      class (type_coupling_task), intent(inout)         :: task
      type (type_internal_variable), pointer            :: master

      type (type_aggregate_variable_access), pointer :: aggregate_variable_access

      master => null()
      if (task%master_standard_variable%aggregate_variable) then
         ! Make sure that an aggregate variable will be created on the fly
         aggregate_variable_access => get_aggregate_variable_access(self, task%master_standard_variable)
         aggregate_variable_access%access = ior(aggregate_variable_access%access, access_read)
         task%master_name = trim(task%master_standard_variable%name)
         task%master_standard_variable => null()
      end if
   end function generate_standard_master

   subroutine aggregate_variable_list_print(self)
      class (type_aggregate_variable_list), intent(in) :: self

      type (type_aggregate_variable),    pointer :: aggregate_variable
      type (type_contributing_variable), pointer :: contributing_variable

      aggregate_variable => self%first
      do while (associated(aggregate_variable))
         contributing_variable => aggregate_variable%first_contributing_variable
         do while (associated(contributing_variable))
            if (associated(contributing_variable%link%target, contributing_variable%link%original)) then
               call log_message('   ' // trim(contributing_variable%link%name) // ' (internal)')
            else
               call log_message('   ' // trim(contributing_variable%link%name) // ' (external)')
            end if
            contributing_variable => contributing_variable%next
         end do
         aggregate_variable => aggregate_variable%next
      end do
   end subroutine aggregate_variable_list_print

   function aggregate_variable_list_get(self, standard_variable) result(aggregate_variable)
      class (type_aggregate_variable_list), intent(inout)    :: self
      class (type_domain_specific_standard_variable), target :: standard_variable

      type (type_aggregate_variable), pointer :: aggregate_variable

      aggregate_variable => self%first
      do while (associated(aggregate_variable))
         if (associated(aggregate_variable%standard_variable, standard_variable)) return
         aggregate_variable => aggregate_variable%next
      end do
      allocate(aggregate_variable)
      aggregate_variable%standard_variable => standard_variable
      aggregate_variable%next => self%first
      self%first => aggregate_variable
   end function

   subroutine aggregate_variable_list_finalize(self)
      class (type_aggregate_variable_list), intent(inout) :: self

      type (type_aggregate_variable),    pointer :: current, next
      type (type_contributing_variable), pointer :: contributing_variable, contributing_variable2

      current => self%first
      do while (associated(current))
         next => current%next
         contributing_variable => current%first_contributing_variable
         do while (associated(contributing_variable))
            contributing_variable2 => contributing_variable%next
            deallocate(contributing_variable)
            contributing_variable => contributing_variable2
         end do
         deallocate(current)
         current => next
      end do
      self%first => null()
   end subroutine

   function collect_aggregate_variables(self) result(list)
      class (type_base_model), intent(in), target :: self
      type (type_aggregate_variable_list)         :: list

      type (type_link),                  pointer :: link
      type (type_contribution),          pointer :: contribution
      type (type_contributing_variable), pointer :: contributing_variable
      type (type_aggregate_variable),    pointer :: aggregate_variable

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
      class (type_base_model), intent(inout), target :: self

      type (type_aggregate_variable_access), pointer :: aggregate_variable_access
      type (type_aggregate_variable),        pointer :: aggregate_variable
      class (type_weighted_sum),             pointer :: sum
      class (type_horizontal_weighted_sum),  pointer :: horizontal_sum
      type (type_contributing_variable),     pointer :: contributing_variable
      type (type_aggregate_variable_list)            :: list
      type (type_model_list_node),           pointer :: child

      ! Get a list of all aggregate variables
      list = collect_aggregate_variables(self)

      ! Make sure that readable fields for all aggregate variables are created at the root level.
      if (.not. associated(self%parent)) then
         aggregate_variable => list%first
         do while (associated(aggregate_variable))
            aggregate_variable_access => get_aggregate_variable_access(self, aggregate_variable%standard_variable)
            aggregate_variable_access%access = ior(aggregate_variable_access%access, access_read)
            aggregate_variable => aggregate_variable%next
         end do
      end if

      aggregate_variable_access => self%first_aggregate_variable_access
      do while (associated(aggregate_variable_access))
         sum => null()
         horizontal_sum => null()

         aggregate_variable => list%get(aggregate_variable_access%standard_variable)
         select type (standard_variable => aggregate_variable%standard_variable)
         class is (type_interior_standard_variable)
            allocate(sum)
         class is (type_surface_standard_variable)
            allocate(horizontal_sum)
            horizontal_sum%domain = domain_surface
         class is (type_bottom_standard_variable)
            allocate(horizontal_sum)
            horizontal_sum%domain = domain_bottom
         class is (type_horizontal_standard_variable)
            allocate(horizontal_sum)
         end select
         contributing_variable => aggregate_variable%first_contributing_variable
         do while (associated(contributing_variable))
            if (associated(contributing_variable%link%target, contributing_variable%link%original) &                  ! Variable must not be coupled
                .and. (associated(self%parent) .or. .not. contributing_variable%link%target%fake_state_variable) &    ! Only include fake state variable for non-root models
                .and. (index(contributing_variable%link%name, '/') == 0 .or. .not. associated(self%parent))) then     ! Variable must be owned by the model itself unless we are aggregating at root level
               select case (contributing_variable%link%target%domain)
               case (domain_interior)
                  call sum%add_component(trim(contributing_variable%link%name), &
                     weight=contributing_variable%scale_factor, include_background=contributing_variable%include_background)
               case (domain_bottom, domain_surface, domain_horizontal)
                  call horizontal_sum%add_component(trim(contributing_variable%link%name), &
                     weight=contributing_variable%scale_factor, include_background=contributing_variable%include_background)
               end select
            end if
            contributing_variable => contributing_variable%next
         end do

         select type (standard_variable => aggregate_variable%standard_variable)
         class is (type_interior_standard_variable)
            sum%units = trim(aggregate_variable%standard_variable%units)
            sum%access = aggregate_variable_access%access
            if (associated(self%parent)) then
               sum%result_output = output_none
            else
               sum%standard_variable => standard_variable
            end if
            if (.not. sum%add_to_parent(self, trim(aggregate_variable%standard_variable%name), aggregate_variable=standard_variable)) deallocate(sum)
         class is (type_horizontal_standard_variable)
            horizontal_sum%units = trim(standard_variable%units)
            horizontal_sum%access = aggregate_variable_access%access
            if (associated(self%parent)) then
               horizontal_sum%result_output = output_none
            else
               horizontal_sum%standard_variable => standard_variable
            end if
            if (.not. horizontal_sum%add_to_parent(self, trim(standard_variable%name), aggregate_variable=standard_variable)) deallocate(horizontal_sum)
         end select
         aggregate_variable_access => aggregate_variable_access%next
      end do

      call list%finalize()

      ! Process child models
      child => self%children%first
      do while (associated(child))
         call create_aggregate_models(child%model)
         child => child%next
      end do
   end subroutine create_aggregate_models

   recursive subroutine create_conservation_checks(self)
      class (type_base_model), intent(inout), target :: self

      type (type_aggregate_variable_list)          :: aggregate_variable_list
      type (type_aggregate_variable),      pointer :: aggregate_variable
      class (type_weighted_sum),           pointer :: sum
      class (type_horizontal_weighted_sum),pointer :: surface_sum, bottom_sum
      type (type_contributing_variable),   pointer :: contributing_variable
      type (type_model_list_node),         pointer :: child
      type (type_standard_variable_set)           :: standard_variable_set
      type (type_standard_variable_node), pointer :: standard_variable_node

      ! Process child models
      child => self%children%first
      do while (associated(child))
         call create_conservation_checks(child%model)
         child => child%next
      end do

      if (.not. self%check_conservation) return

      aggregate_variable_list = collect_aggregate_variables(self)

      ! Get list of conserved quantities (map to universal=domain-independent variables where possible)
      aggregate_variable => aggregate_variable_list%first
      do while (associated(aggregate_variable))
         if (associated(aggregate_variable%standard_variable%universal)) then
            if (aggregate_variable%standard_variable%universal%conserved) call standard_variable_set%add(aggregate_variable%standard_variable%universal)
         end if
         aggregate_variable => aggregate_variable%next
      end do

      standard_variable_node => standard_variable_set%first
      do while (associated(standard_variable_node))
         select type (standard_variable => standard_variable_node%p)
         class is (type_universal_standard_variable)
            ! Allocate objects that will sum fluxes for each of the different domains (interior, surface, bottom).
            allocate(sum, surface_sum, bottom_sum)

            ! Enumerate contributions to interior field.
            aggregate_variable => aggregate_variable_list%get(standard_variable%in_interior())
            contributing_variable => aggregate_variable%first_contributing_variable
            do while (associated(contributing_variable))
               if (index(contributing_variable%link%name, '/') == 0) then
                  if (associated(contributing_variable%link%original%sms))          call sum%add_component        (contributing_variable%link%original%sms%name,          contributing_variable%scale_factor)
                  if (associated(contributing_variable%link%original%surface_flux)) call surface_sum%add_component(contributing_variable%link%original%surface_flux%name, contributing_variable%scale_factor)
                  if (associated(contributing_variable%link%original%bottom_flux))  call bottom_sum%add_component (contributing_variable%link%original%bottom_flux%name,  contributing_variable%scale_factor)
               end if
               contributing_variable => contributing_variable%next
            end do

            ! Enumerate contributions to surface field.
            aggregate_variable => aggregate_variable_list%get(standard_variable%at_surface())
            contributing_variable => aggregate_variable%first_contributing_variable
            do while (associated(contributing_variable))
               if (index(contributing_variable%link%name, '/') == 0 .and. associated(contributing_variable%link%original%sms)) &
                  call surface_sum%add_component(contributing_variable%link%original%sms%name, contributing_variable%scale_factor)
               contributing_variable => contributing_variable%next
            end do

            ! Enumerate contributions to bottom field.
            aggregate_variable => aggregate_variable_list%get(standard_variable%at_bottom())
            contributing_variable => aggregate_variable%first_contributing_variable
            do while (associated(contributing_variable))
               if (index(contributing_variable%link%name, '/') == 0 .and. associated(contributing_variable%link%original%sms)) &
                  call bottom_sum%add_component(contributing_variable%link%original%sms%name, contributing_variable%scale_factor)
               contributing_variable => contributing_variable%next
            end do

            ! Process sums now that all contributing terms are known.
            sum%units = trim(standard_variable%units) // '/s'
            if (.not. sum%add_to_parent(self,'change_in_' // trim(standard_variable%name), create_for_one=.true.)) deallocate(sum)
            surface_sum%units = trim(standard_variable%units) // '*m/s'
            if (.not. surface_sum%add_to_parent(self,'change_in_' // trim(standard_variable%name) // '_at_surface', create_for_one=.true.)) deallocate(surface_sum)
            bottom_sum%units = trim(standard_variable%units) // '*m/s'
            if (.not. bottom_sum%add_to_parent(self,'change_in_' // trim(standard_variable%name) // '_at_bottom', create_for_one=.true.)) deallocate(bottom_sum)
         end select

         standard_variable_node => standard_variable_node%next
      end do

   end subroutine create_conservation_checks

   recursive subroutine couple_variables(self,master,slave)
      class (type_base_model), intent(inout), target :: self
      type (type_internal_variable), pointer         :: master, slave

      type (type_link_pointer), pointer :: link_pointer, next_link_pointer
      logical                           :: slave_is_state_variable
      logical                           :: master_is_state_variable

      ! If slave and master are the same, we are done - return.
      if (associated(slave, master)) return

      slave_is_state_variable  = slave%source == source_state .or. slave%fake_state_variable
      master_is_state_variable = master%source == source_state .or. master%fake_state_variable

      if (associated(self%parent)) call self%fatal_error('couple_variables', 'BUG: must be called on root node.')
      if (associated(slave%write_index)) &
         call fatal_error('couple_variables', 'Attempt to couple write-only variable ' &
            // trim(slave%name) // ' to ' // trim(master%name) // '.')
      if (slave_is_state_variable) then
         ! Extra checks when coupling state variables
         if (.not. master_is_state_variable) then
            if (slave%presence == presence_external_optional) return
            call fatal_error('couple_variables', 'Attempt to couple state variable ' &
               // trim(slave%name) // ' to non-state variable ' // trim(master%name) // '.')
         end if
         if ((slave%domain == domain_bottom .and. master%domain == domain_surface) .or. (slave%domain == domain_surface .and. master%domain == domain_bottom)) &
            call fatal_error('couple_variables', 'Cannot couple ' // trim(slave%name) // ' to '//trim(master%name) // ', because their domains are incompatible.')
      end if
      !if (slave%domain/=master%domain.and..not.(slave%domain==domain_horizontal.and. &
      !   (master%domain==domain_surface.or.master%domain==domain_bottom))) call fatal_error('couple_variables', &
      !   'Cannot couple '//trim(slave%name)//' to '//trim(master%name)//', because their domains are incompatible.')
      if (iand(slave%domain, master%domain) == 0) call fatal_error('couple_variables', &
         'Cannot couple ' // trim(slave%name) // ' to ' // trim(master%name) // ', because their domains are incompatible.')

      if (debug_coupling) call log_message(trim(slave%name) // ' --> ' // trim(master%name))

      ! Merge all information from the slave into the master.
      call master%state_indices%extend(slave%state_indices)
      call master%read_indices%extend(slave%read_indices)
      call master%write_indices%extend(slave%write_indices)
      call master%background_values%extend(slave%background_values)
      call master%properties%update(slave%properties,overwrite=.false.)
      call master%sms_list%extend(slave%sms_list)
      call master%surface_flux_list%extend(slave%surface_flux_list)
      call master%bottom_flux_list%extend(slave%bottom_flux_list)
      call master%movement_list%extend(slave%movement_list)

      call master%standard_variables%update(slave%standard_variables)

      if (master%presence == presence_external_optional .and. slave%presence /= presence_external_optional) &
         master%presence = presence_external_required

      ! Make all links that originally pointed to the slave now point to the master.
      ! Then include those links in the master's link set.
      link_pointer => slave%first_link
      slave%first_link => null()
      do while (associated(link_pointer))
         next_link_pointer => link_pointer%next
         link_pointer%p%target => master
         link_pointer%next => master%first_link
         master%first_link => link_pointer
         link_pointer => next_link_pointer
      end do
   end subroutine couple_variables

   recursive subroutine check_coupling_units(self)
      class (type_base_model), intent(in), target :: self

      type (type_link),            pointer :: link
      type (type_model_list_node), pointer :: child

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

end module

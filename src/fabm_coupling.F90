module fabm_coupling
   use fabm_types
   use fabm_builtin_sum
   use fabm_driver
   use yaml_settings, only: default_minimum_real, default_maximum_real

   implicit none

   private

   public freeze_model_info
   public collect_aggregate_variables, type_aggregate_variable_list, type_aggregate_variable

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

   subroutine freeze_model_info(self, require_initialization, coupling_log_unit)
      class (type_base_model), intent(inout), target :: self
      logical,                 intent(in)            :: require_initialization
      integer,                 intent(in)            :: coupling_log_unit

      if (associated(self%parent)) call self%fatal_error('freeze_model_info', &
         'BUG: freeze_model_info can only operate on the root model.')

      call before_coupling(self)

      call get_initial_state(self, require_initialization)

      ! Coupling stage 1: implicit - couple variables based on overlapping standard identities.
      call couple_standard_variables(self)

      ! Coupling stage 2: explicit - resolve user- or model-specified links between variables.
      if (coupling_log_unit > 0) write (coupling_log_unit, '(a)') 'process_coupling_tasks_1:'
      call process_coupling_tasks(self, final=.false., log_unit=coupling_log_unit)

      ! Create models for aggregate variables:
      ! * at root level, e.g. for totals of conserved quantities [included in default output]
      ! * at instance level if requested by coupling
      ! After this step, the set of variables that contribute to aggregate quantities may not be modified.
      ! That is, no new such variables may be added, and no such variables may be coupled.
      ! This step will typically create new child models to handle the necessary summations.
      ! These child models may also add source terms (if the sum is treated as state variable),
      ! so any source term treatment (create_conservation_checks, create_flux_sums) should happen after this is complete!
      call create_aggregate_models(self)

      ! Perform coupling for any new aggregate models.
      ! This may append items to existing lists of source terms and bottom/surface fluxes,
      ! so it has to precede the call to create_flux_sums.
      if (coupling_log_unit > 0) write (coupling_log_unit, '(a)') 'process_coupling_tasks_2:'
      call process_coupling_tasks(self, final=.false., log_unit=coupling_log_unit)

      ! Now that coupling for non-rate variables is complete, contributions to aggregate quantities
      ! (including contributions from variables that were coupled to another) are final.
      ! Create conservation checks where needed.
      call create_conservation_checks(self)

      ! Create summations of source terms and surface/bottom fluxes.
      call create_flux_sums(self)
      call request_flux_sum_coupling(self)

      ! Process the any remaining coupling tasks.
      if (coupling_log_unit > 0) write (coupling_log_unit, '(a)') 'process_coupling_tasks_3:'
      call process_coupling_tasks(self, final=.true., log_unit=coupling_log_unit)

      ! Allow inheriting models to perform additional tasks after coupling.
      call after_coupling(self)

      ! Check whether units of coupled variables match
      !call check_coupling_units(self)

      call freeze(self)

      call report_ignored_sources(self)
   end subroutine freeze_model_info

   recursive subroutine before_coupling(self)
      class (type_base_model), intent(inout), target :: self

      type (type_model_list_node), pointer :: node

      call self%before_coupling()
      node => self%children%first
      do while (associated(node))
         call before_coupling(node%model)
         node => node%next
      end do
   end subroutine

   recursive subroutine get_initial_state(self, require_initialization)
      class (type_base_model), intent(inout) :: self
      logical,                 intent(in)    :: require_initialization

      type (type_link),            pointer :: link
      real(rk)                             :: minimum
      real(rk)                             :: maximum
      type (type_model_list_node), pointer :: node

      ! Transfer user-specified initial state to the model.
      link => self%links%first
      do while (associated(link))
         minimum = default_minimum_real
         maximum = default_maximum_real
         if (link%target%minimum /= -1.e20_rk) minimum = link%target%minimum
         if (link%target%maximum /=  1.e20_rk) maximum = link%target%maximum
         if (index(link%name, '/') == 0 .and. link%target%source == source_state .and. link%target%presence == presence_internal) then
            if (require_initialization) then
               call self%initialization%get(link%target%initial_value, trim(link%name), trim(link%target%long_name), &
                  trim(link%target%units), minimum=minimum, maximum=maximum)
            else
               call self%initialization%get(link%target%initial_value, trim(link%name), trim(link%target%long_name), &
                  trim(link%target%units), minimum=minimum, maximum=maximum, default=link%target%initial_value)
            end if
         end if
         link => link%next
      end do

      node => self%children%first
      do while (associated(node))
         call get_initial_state(node%model, require_initialization)
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

   subroutine report_ignored_sources(model)
      class (type_base_model), intent(in) :: model

      type (type_link), pointer :: link

      link => model%links%first
      do while (associated(link))
         if (associated(link%original%sms) .and. .not. (link%target%source == source_state .or. link%target%fake_state_variable)) &
            call model%log_message('Warning: sources provided for ' // trim(link%name) // ' will be ignored as it has been coupled&
               & to a non-state variable (' // trim(link%target%name) // ').')
         link => link%next
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
                     ! Default coupling: later variable (link) is coupled to early variable (first_link)
                     call couple_variables(model, link%target, first_link%target)
                  else
                     ! Later variable (link) is write-only and therefore cannot be coupled to a target.
                     ! Try coupling early variable (first_link) to later variable (link).
                     call couple_variables(model, first_link%target, link%target)
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
      character(len=:), allocatable       :: target_name
      class (type_coupling_task), pointer :: task 
      integer                             :: source
      logical                             :: couplable
      integer                             :: display

      link => self%links%first
      do while (associated(link))
         ! Only process own links (those without slash in the name)
         source = link%original%source
         couplable = source == source_state .or. source == source_unknown
         if (index(link%name, '/') == 0 .and. couplable) then
            display = display_inherit
            if (link%original%presence == presence_internal .or. associated(self%coupling_task_list%find(link))) display = display_advanced
            target_name = self%couplings%get_string(trim(link%name), trim(link%original%long_name), units=trim(link%original%units), default='', display=display)
            if (len(target_name) >= 4) then
               if (target_name(len(target_name) - 3:len(target_name)) == '@old') then
                  call set_dependency_flag(link%original, source=-1, flag=dependency_flag_stale)
                  target_name = target_name(:len(target_name) - 4)
               end if
            end if
            if (target_name /= '') then
               allocate(task)
               task%link => link
               task%target_name = target_name
               call self%coupling_task_list%add(task, priority=1)
            end if    ! Coupling provided
         end if   ! Our own link, which may be coupled
         link => link%next
      end do

      self%coupling_task_list%includes_custom = .true.
   end subroutine collect_user_specified_couplings

   recursive subroutine process_coupling_tasks(self, final, log_unit)
      class (type_base_model), intent(inout), target :: self
      logical,                 intent(in)            :: final
      integer,                 intent(in)            :: log_unit

      class (type_base_model),       pointer :: root
      type (type_model_list_node),   pointer :: child
      class (type_coupling_task),    pointer :: coupling, next_coupling
      type (type_internal_variable), pointer :: target_variable
      type (type_link),              pointer :: link
      integer                                :: istart, istop

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

         ! First try if the coupling object can resolve the variable reference itself
         ! (e.g., type_coupling_from_model in fabm_particle)
         link => coupling%resolve()

         if (.not. associated(link) .and. .not. associated(coupling%target_standard_variable)) then
            ! This is a coupling by variable name.
            ! Try to find the target variable among the variables of the requesting model or its parents.
            istart = index(coupling%target_name, '(')
            if (istart /= 0) then
               ! The coupling name includes an opening parenthesis. Interpret it as a parametrized coupling (one with arguments)
               istop = len_trim(coupling%target_name)
               if (coupling%target_name(istop:istop) /= ')') call self%fatal_error('process_coupling_tasks', &
                  'Parameterized coupling ' // trim(coupling%target_name) // ' should end with closing parenthesis.')
               call resolve_parameterized_coupling(coupling%target_name(1:istart-1), coupling%target_name(istart+1:istop-1), coupling)
            elseif (coupling%link%name /= coupling%target_name) then
               ! Names of variable and its target differ: start target search in current model, then move up tree.
               link => self%find_link(coupling%target_name, recursive=.true., exact=.false.)
            elseif (associated(self%parent)) then
               ! Names of variable and its target are identical: start target search in parent model, then move up tree.
               link => self%parent%find_link(coupling%target_name, recursive=.true., exact=.false.)
            else
               call self%fatal_error('process_coupling_tasks', &
                  'Names of variable and its target are identical: "' // trim(coupling%target_name) // '". This is not valid at the root of the model tree.')
            end if
         end if

         if (.not. associated(link) .and. associated(coupling%target_standard_variable)) then
            ! This is a coupling to a standard variable. First try to find the corresponding standard variable.
            ! We search within the root model, because there all variables are found together.
            link => root%links%first
            do while (associated(link))
               if (link%target%standard_variables%contains(coupling%target_standard_variable)) exit
               link => link%next
            end do

            if (.not. associated(link) .and. coupling%target_standard_variable%aggregate_variable) &
               ! Create an aggregate variable at the level of the root model
               link => get_aggregate_variable_access(root, coupling%target_standard_variable)

            if (final .and. .not. associated(link) .and. (coupling%link%target%source /= source_state &
               .or. coupling%link%target%presence == presence_external_optional)) then
               ! Target variable was not found, but this is our last chance.
               ! Therefore, create a placeholder variable at the root level.
               ! This variable will still need to be provided by the host.
               select type (standard_variable => coupling%target_standard_variable)
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
         end if

         ! Save pointer to the next coupling task in advance, because current task may
         ! be deallocated from self%coupling_task_list%remove.
         next_coupling => coupling%next

         if (associated(link)) then
            ! Target variable found: perform the coupling.
            target_variable => link%target
            if (log_unit /= -1) write (log_unit, '(a,a,a,a)') '  ', trim(coupling%link%target%name), ': ', trim(target_variable%name)
            call couple_variables(root, coupling%link%target, target_variable)

            ! Remove coupling task from the list
            call self%coupling_task_list%remove(coupling)
         elseif (final) then
            call self%fatal_error('process_coupling_tasks', &
               'Coupling target "' // trim(coupling%target_name) // '" for "' // trim(coupling%link%name) // '" was not found.')
         end if

         ! Move to next coupling task.
         coupling => next_coupling
      end do

      ! Process coupling tasks registered with child models.
      child => self%children%first
      do while (associated(child))
         call process_coupling_tasks(child%model, final, log_unit)
         child => child%next
      end do

   contains

      subroutine resolve_parameterized_coupling(name, args, task)
         character(len=*),           intent(in)    :: name, args
         class (type_coupling_task), intent(inout) :: task

         type (type_interior_standard_variable)   :: interior_standard_variable
         type (type_bottom_standard_variable)     :: bottom_standard_variable
         type (type_surface_standard_variable)    :: surface_standard_variable
         type (type_horizontal_standard_variable) :: horizontal_standard_variable

         select case (name)
         case ('standard_variable')
            select case (task%link%target%domain)
            case (domain_interior)
               interior_standard_variable%name = args
               task%target_standard_variable => interior_standard_variable%typed_resolve()
            case (domain_bottom)
               bottom_standard_variable%name = args
               task%target_standard_variable => bottom_standard_variable%typed_resolve()
            case (domain_horizontal)
               horizontal_standard_variable%name = args
               task%target_standard_variable => horizontal_standard_variable%typed_resolve()
            case default
               call self%fatal_error('process_coupling_tasks', 'Unknown domain for ' // task%link%name // '.')
            end select
         case default
            call self%fatal_error('process_coupling_tasks', 'Unknown parameterized coupling type "' // name // '".')
         end select
      end subroutine

   end subroutine process_coupling_tasks

   recursive subroutine create_flux_sums(self)
      class (type_base_model), intent(inout), target :: self

      type (type_link),            pointer :: link
      type (type_model_list_node), pointer :: child

      link => self%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. (link%original%source == source_state .or. link%original%fake_state_variable)) then
            ! This is a state variable, or a diagnostic pretending to be one, that we have registered (it is owned by "self")

            ! Create variable links for sources and surface/bottom fluxes at the same level of the original state variable.
            ! E.g., a state variable named X will get a sibling X_sms_tot that will point to the total of its sources minus sinks.
            ! These links will be coupled to the result of the source/flux summations at the next stage.
            select case (link%target%domain)
            case (domain_interior)
               call self%add_interior_variable  (trim(link%name) // '_sms_tot', trim(link%original%units) // ' s-1',   trim(link%original%long_name) // ' total sources', link=link%original%sms_sum, output=output_none, presence=presence_external_required)
               call self%add_horizontal_variable(trim(link%name) // '_sfl_tot', trim(link%original%units) // ' m s-1', trim(link%original%long_name) // ' total surface flux', domain=domain_surface, link=link%original%surface_flux_sum, output=output_none, presence=presence_external_required)
               call self%add_horizontal_variable(trim(link%name) // '_bfl_tot', trim(link%original%units) // ' m s-1', trim(link%original%long_name) // ' total bottom flux', domain=domain_bottom, link=link%original%bottom_flux_sum, output=output_none, presence=presence_external_required)
               call self%add_interior_variable  (trim(link%name) // '_w_tot',   'm s-1',                               trim(link%original%long_name) // ' total vertical velocity', link=link%original%movement_sum, output=output_none, presence=presence_external_required)
            case (domain_horizontal, domain_surface, domain_bottom)
               call self%add_horizontal_variable(trim(link%name) // '_sms_tot', trim(link%original%units) // ' s-1',   trim(link%original%long_name) // ' total sources', domain=link%target%domain, link=link%original%sms_sum, output=output_none, presence=presence_external_required)
            end select

            if (associated(link%target, link%original)) then
               ! We own this variable (it has not been coupled to another).
               ! Create summations for sources-sinks and surface/bottom fluxes.
               select case (link%target%domain)
               case (domain_interior)
                  call create_sum(link%target%sms_sum, link%target%sms_list, domain_interior)
                  call create_sum(link%target%surface_flux_sum, link%target%surface_flux_list, domain_surface)
                  call create_sum(link%target%bottom_flux_sum, link%target%bottom_flux_list, domain_bottom)
                  call create_sum(link%target%movement_sum, link%target%movement_list, domain_interior)
               case (domain_horizontal, domain_surface, domain_bottom)
                  call create_sum(link%target%sms_sum, link%target%sms_list, link%target%domain)
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

   contains

      subroutine create_sum(link, link_list, domain)
         type (type_link),      intent(in), target :: link
         type (type_link_list), intent(in)         :: link_list
         integer,               intent(in)         :: domain

         class (type_base_sum), pointer :: sum
         type (type_link),      pointer :: component_link

         if (domain == domain_interior) then
            allocate(type_weighted_sum::sum)
         else
            allocate(type_horizontal_weighted_sum::sum)
         end if
         sum%result_output = output_always_available
         sum%missing_value = 0.0_rk
         sum%units = link%target%units
         component_link => link_list%first
         do while (associated(component_link))
            call sum%add_component(component_link)
            component_link => component_link%next
         end do
         call self%add_child(sum, trim(link%name) // '_calculator')
         call self%request_coupling(link, sum%result_link)
      end subroutine create_sum

   end subroutine create_flux_sums

   recursive subroutine request_flux_sum_coupling(self)
      class (type_base_model), intent(inout), target :: self

      type (type_link),            pointer :: link
      type (type_model_list_node), pointer :: child

      link => self%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. .not. associated(link%target, link%original) .and. (link%original%source == source_state .or. link%original%fake_state_variable) .and. (link%target%source == source_state .or. link%target%fake_state_variable)) then
            ! This is a state variable, or a diagnostic pretending to be one, that we have registered (it is owned by "self")
            ! We do not own this variable.
            ! Couple to summations for sources-sinks and surface/bottom fluxes created by the target.
            select case (link%target%domain)
            case (domain_interior)
               call self%request_coupling(link%original%sms_sum, link%target%sms_sum)
               call self%request_coupling(link%original%surface_flux_sum, link%target%surface_flux_sum)
               call self%request_coupling(link%original%bottom_flux_sum, link%target%bottom_flux_sum)
               call self%request_coupling(link%original%movement_sum, link%target%movement_sum)
            case (domain_horizontal, domain_surface, domain_bottom)
               call self%request_coupling(link%original%sms_sum, link%target%sms_sum)
            end select
         end if
         link => link%next
      end do

      ! Process child models
      child => self%children%first
      do while (associated(child))
         call request_flux_sum_coupling(child%model)
         child => child%next
      end do
   end subroutine request_flux_sum_coupling

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
      logical,                        pointer :: pmember

      pmember => standard_variable%aggregate_variable
      aggregate_variable => self%first
      do while (associated(aggregate_variable))
         ! Note: for Cray 10.0.4, the comparison below fails for class pointers! Therefore we compare type member references.
         if (associated(pmember, aggregate_variable%standard_variable%aggregate_variable)) return
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
      class (type_weighted_sum),             pointer :: interior_sum
      class (type_horizontal_weighted_sum),  pointer :: horizontal_sum
      class (type_base_sum),                 pointer :: sum
      type (type_contributing_variable),     pointer :: contributing_variable
      type (type_aggregate_variable_list)            :: list
      type (type_model_list_node),           pointer :: child
      type (type_link),                      pointer :: link

      ! Get a list of all aggregate variables
      list = collect_aggregate_variables(self)

      ! Make sure that readable fields for all aggregate variables are created at the root level.
      if (.not. associated(self%parent)) then
         aggregate_variable => list%first
         do while (associated(aggregate_variable))
            link => get_aggregate_variable_access(self, aggregate_variable%standard_variable)
            link%target%output = ior(output_instantaneous, output_always_available)
            aggregate_variable => aggregate_variable%next
         end do
      end if

      aggregate_variable_access => self%first_aggregate_variable_access
      do while (associated(aggregate_variable_access))
         aggregate_variable => list%get(aggregate_variable_access%standard_variable)

         ! Allocate and configure child model that will do the summation
         select type (standard_variable => aggregate_variable%standard_variable)
         class is (type_interior_standard_variable)
            allocate(interior_sum)
            interior_sum%aggregate_variable => standard_variable
            sum => interior_sum
         class is (type_horizontal_standard_variable)
            allocate(horizontal_sum)
            horizontal_sum%aggregate_variable => standard_variable
            horizontal_sum%domain = standard_variable2domain(standard_variable)
            sum => horizontal_sum
         end select
         sum%units = aggregate_variable_access%link%target%units
         sum%act_as_state_variable = aggregate_variable_access%link%target%fake_state_variable
         sum%result_output = output_none

         contributing_variable => aggregate_variable%first_contributing_variable
         do while (associated(contributing_variable))
            if (associated(contributing_variable%link%target, contributing_variable%link%original) &                  ! Variable must not be coupled
                .and. (associated(self%parent) .or. .not. contributing_variable%link%target%fake_state_variable) &    ! Only include fake state variable for non-root models
                .and. (index(contributing_variable%link%name, '/') == 0 .or. .not. associated(self%parent))) &        ! Variable must be owned by the model itself unless we are aggregating at root level
               call sum%add_component(contributing_variable%link, &
                  weight=contributing_variable%scale_factor, include_background=contributing_variable%include_background)
            contributing_variable => contributing_variable%next
         end do

         call self%add_child(sum, trim(aggregate_variable_access%link%name) // '_calculator')
         call self%request_coupling(aggregate_variable_access%link, sum%result_link)

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
      type (type_standard_variable_set)            :: standard_variable_set
      type (type_standard_variable_node),  pointer :: standard_variable_node
      type (type_link),                    pointer :: link

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
            sum%result_output = output_none
            surface_sum%result_output = output_none
            bottom_sum%result_output = output_none

            ! Enumerate contributions to interior field.
            aggregate_variable => aggregate_variable_list%get(standard_variable%in_interior())
            contributing_variable => aggregate_variable%first_contributing_variable
            do while (associated(contributing_variable))
               if (index(contributing_variable%link%name, '/') == 0) then
                  if (associated(contributing_variable%link%original%sms))          call sum%add_component        (contributing_variable%link%original%sms,          contributing_variable%scale_factor)
                  if (associated(contributing_variable%link%original%surface_flux)) call surface_sum%add_component(contributing_variable%link%original%surface_flux, contributing_variable%scale_factor)
                  if (associated(contributing_variable%link%original%bottom_flux))  call bottom_sum%add_component (contributing_variable%link%original%bottom_flux,  contributing_variable%scale_factor)
               end if
               contributing_variable => contributing_variable%next
            end do

            ! Enumerate contributions to surface field.
            aggregate_variable => aggregate_variable_list%get(standard_variable%at_surface())
            contributing_variable => aggregate_variable%first_contributing_variable
            do while (associated(contributing_variable))
               if (index(contributing_variable%link%name, '/') == 0 .and. associated(contributing_variable%link%original%sms)) &
                  call surface_sum%add_component(contributing_variable%link%original%sms, contributing_variable%scale_factor)
               contributing_variable => contributing_variable%next
            end do

            ! Enumerate contributions to bottom field.
            aggregate_variable => aggregate_variable_list%get(standard_variable%at_bottom())
            contributing_variable => aggregate_variable%first_contributing_variable
            do while (associated(contributing_variable))
               if (index(contributing_variable%link%name, '/') == 0 .and. associated(contributing_variable%link%original%sms)) &
                  call bottom_sum%add_component(contributing_variable%link%original%sms, contributing_variable%scale_factor)
               contributing_variable => contributing_variable%next
            end do

            ! Process sums now that all contributing terms are known.
            link => null()
            call self%add_interior_variable('change_in_' // trim(standard_variable%name), trim(standard_variable%units) // ' s-1', 'change in ' // trim(standard_variable%name), link=link, output=ior(output_instantaneous, output_always_available), presence=presence_external_required)
            call self%add_child(sum, trim(link%name) // '_calculator')
            call self%request_coupling(link, sum%result_link)
            link => null()
            call self%add_horizontal_variable('change_in_' // trim(standard_variable%name) // '_at_surface', trim(standard_variable%units) // ' m s-1', 'change in ' // trim(standard_variable%name) // ' at surface', domain=domain_surface, link=link, output=ior(output_instantaneous, output_always_available), presence=presence_external_required)
            call self%add_child(surface_sum, trim(link%name) // '_calculator')
            call self%request_coupling(link, surface_sum%result_link)
            link => null()
            call self%add_horizontal_variable('change_in_' // trim(standard_variable%name) // '_at_bottom', trim(standard_variable%units) // ' m s-1', 'change in ' // trim(standard_variable%name)// ' at bottom', domain=domain_bottom, link=link, output=ior(output_instantaneous, output_always_available), presence=presence_external_required)
            call self%add_child(bottom_sum, trim(link%name) // '_calculator')
            call self%request_coupling(link, bottom_sum%result_link)
         end select

         standard_variable_node => standard_variable_node%next
      end do

   end subroutine create_conservation_checks

   recursive subroutine couple_variables(self, variable, target_variable)
      class (type_base_model), intent(inout), target :: self
      type (type_internal_variable), pointer         :: variable, target_variable

      type (type_link_pointer), pointer :: link_pointer, next_link_pointer

      ! If the variable and its target are the same, we are done - return.
      if (associated(variable, target_variable)) return

      if (associated(self%parent)) call self%fatal_error('couple_variables', 'BUG: must be called on root node.')
      if (associated(variable%write_index)) &
         call fatal_error('couple_variables', 'Attempt to couple write-only variable ' &
            // trim(variable%name) // ' to ' // trim(target_variable%name) // '.')
      if (variable%source == source_state .or. variable%fake_state_variable) then
         ! Extra checks when coupling state variables
         if ((variable%domain == domain_bottom .and. target_variable%domain == domain_surface) .or. (variable%domain == domain_surface .and. target_variable%domain == domain_bottom)) &
            call fatal_error('couple_variables', &
               'Cannot couple ' // trim(variable%name) // ' (' // trim(domain2string(variable%domain)) // ') to ' // trim(target_variable%name) &
               // ' (' // trim(domain2string(target_variable%domain)) // '), because their domains are incompatible.')
      end if
      !if (variable%domain/=target_variable%domain .and. .not. (variable%domain==domain_horizontal.and. &
      !   (target_variable%domain==domain_surface.or.target_variable%domain==domain_bottom))) call fatal_error('couple_variables', &
      !   'Cannot couple '//trim(variable%name)//' to '//trim(target_variable%name)//', because their domains are incompatible.')
      if (iand(variable%domain, target_variable%domain) == 0) &
         call fatal_error('couple_variables', &
         'Cannot couple ' // trim(variable%name) // ' (' // trim(domain2string(variable%domain)) // ') to ' // trim(target_variable%name) &
         // ' (' // trim(domain2string(target_variable%domain)) // '), because their domains are incompatible.')

      ! Merge all information into the target variable.
      call target_variable%state_indices%extend(variable%state_indices)
      call target_variable%read_indices%extend(variable%read_indices)
      call target_variable%write_indices%extend(variable%write_indices)
      call target_variable%background_values%extend(variable%background_values)
      call target_variable%properties%update(variable%properties,overwrite=.false.)
      call target_variable%sms_list%extend(variable%sms_list)
      call target_variable%surface_flux_list%extend(variable%surface_flux_list)
      call target_variable%bottom_flux_list%extend(variable%bottom_flux_list)
      call target_variable%movement_list%extend(variable%movement_list)

      call target_variable%standard_variables%update(variable%standard_variables)

      if (target_variable%presence == presence_external_optional .and. variable%presence /= presence_external_optional) &
         target_variable%presence = presence_external_required

      ! Make all links that originally pointed to the coupled variable now point to its target.
      ! Then include those links in the target's link set.
      link_pointer => variable%first_link
      variable%first_link => null()
      do while (associated(link_pointer))
         next_link_pointer => link_pointer%next
         link_pointer%p%target => target_variable
         link_pointer%next => target_variable%first_link
         target_variable%first_link => link_pointer
         link_pointer => next_link_pointer
      end do
   end subroutine couple_variables

   recursive subroutine check_coupling_units(self)
      class (type_base_model), intent(in), target :: self

      type (type_link),            pointer :: link
      type (type_model_list_node), pointer :: child

      link => self%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. .not. associated(link%target, link%original)) then
            if (link%target%units/=''.and. link%original%units/=''.and. link%target%units /= link%original%units) &
               call log_message('WARNING: unit mismatch between ' // trim(link%original%name) &
                   // ' (' // trim(link%original%units) // ') and its coupling target ' // trim(link%target%name) &
                   // ' (' // trim(link%target%units) // ').')
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

module fabm_coupling
   use fabm_types
   use fabm_standard_variables
   use fabm_properties
   use fabm_builtin_models
   use fabm_driver

   implicit none

   private

   public freeze_model_info, find_dependencies

   logical,parameter :: debug_coupling = .false.

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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
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
      call process_coupling_tasks(self,0)

      ! Coupling stage 3: try to automatically fulfil remaining dependencies on
      ! - aggregated sources-sinks of state variables (or of diagnostics acting as such)
      ! - aggregate standard variables
      ! This may create new child models to handle the necessary summations.
      call process_coupling_tasks(self,1)

      ! Create models for aggregate variables at root level, to be used to compute conserved quantities.
      ! After this step, the set of variables that contribute to aggregate quantities may not be modified.
      ! That is, no new such variables may be added, and no such variables may be coupled.
      call build_aggregate_variables(self)
      call create_aggregate_models(self)

      ! Repeat coupling because new aggregate variables are now available.
      call process_coupling_tasks(self,2)

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
      type (type_link),pointer :: link,link2
      type (type_set)          :: processed
!
!-----------------------------------------------------------------------
!BOC
      link => model%links%first
      do while (associated(link))
         if (associated(link%target%standard_variable)) then
            if (.not.processed%contains(link%target%standard_variable%name)) then
               call processed%add(link%target%standard_variable%name)
               link2 => link%next
               do while (associated(link2))
                  if (associated(link2%target%standard_variable)) then
                     if (link%target%standard_variable%compare(link2%target%standard_variable)) then
                        if (link2%target%write_indices%is_empty().and..not.(link%target%presence/=presence_internal.and.link2%target%presence==presence_internal)) then
                           ! Default coupling: early variable is master, later variable is slave.
                           call couple_variables(model,link%target,link2%target)
                        else
                           ! Later variable is write-only and therefore can only be master. Try coupling with early variable as slave.
                           call couple_variables(model,link2%target,link%target)
                        end if
                     end if
                  end if
                  link2 => link2%next
               end do
            end if ! if link%target%standard_variable%name not processed yet
         end if ! if associated(link%target%standard_variable)
         link => link%next
      end do

      call processed%finalize()

   end subroutine couple_standard_variables
!EOC

   subroutine collect_user_specified_couplings(self)
      class (type_base_model),intent(inout) :: self

      type (type_link),         pointer :: link
      class (type_property),    pointer :: master_name
      type (type_coupling_task),pointer :: task

      link => self%links%first
      do while (associated(link))
         ! Only process own links (those without slash in the name)
         if (index(link%name,'/')==0) then
            master_name => self%couplings%find_in_tree(link%name)
            if (associated(master_name)) then
               task => self%coupling_task_list%add(link,.true.)
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
      type (type_coupling_task),    pointer :: coupling, next_coupling
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

         select case (stage)
            case (0,2)
               if (associated(coupling%master_standard_variable)) then
                  ! Coupling to a standard variable - first try to find the corresponding standard variable.
                  link => root%links%first
                  do while (associated(link))
                     if (associated(link%target%standard_variable)) then
                        if (coupling%master_standard_variable%compare(link%target%standard_variable)) exit
                     end if
                     link => link%next
                  end do

                  if (stage==2.and..not.associated(link)) then
                     ! This is our last chance - create an appropriate variable at the root level.
                     select type (standard_variable=>coupling%master_standard_variable)
                     class is (type_bulk_standard_variable)
                        call root%add_bulk_variable(standard_variable%name, standard_variable%units, standard_variable%name, &
                                                    standard_variable=standard_variable, presence=presence_external_optional, link=link)
                     class is (type_horizontal_standard_variable)
                        call root%add_horizontal_variable(standard_variable%name, standard_variable%units, standard_variable%name, &
                                                          standard_variable=standard_variable, presence=presence_external_optional, link=link)
                     class is (type_global_standard_variable)
                        call root%add_scalar_variable(standard_variable%name, standard_variable%units, standard_variable%name, &
                                                      standard_variable=standard_variable, presence=presence_external_optional, link=link)
                     end select
                  end if

                  ! If we found a link to the required master variable, save a pointer to the master.
                  nullify(master)
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
            case (1)
               if (associated(coupling%master_standard_variable)) then
                  master => generate_standard_master(root,coupling)
               else
                  master => generate_master(self,coupling%master_name)
               end if
         end select

         ! Save pointer to the next coupling task in advance, because current task may
         ! be deallocated from self%coupling_task_list%remove.
         next_coupling => coupling%next

         if (associated(master)) then
            ! Target variable found: perform the coupling.
            call couple_variables(root,master,coupling%slave%target)

            ! Remove coupling task from the list
            call self%coupling_task_list%remove(coupling)
         elseif (stage==2) then
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

   function generate_master(self,name) result(master)
      class (type_base_model),intent(inout),target :: self
      character(len=*),       intent(in)           :: name

      type (type_internal_variable),pointer :: master
      type (type_link),             pointer :: link
      integer :: istart,n
      type (type_weighted_sum),           pointer :: sum
      type (type_horizontal_weighted_sum),pointer :: sum_hz

      nullify(master)
      n = len_trim(name)
      if (n<8) return
      if (name(n-7:n)=='_sms_tot') then
         ! This is a link to the combined sources-sinks of a variable.
         ! First locate the relevant variable, then create the necessary sms summation.
         master => self%find_object(name(:n-8),recursive=.true.,exact=.false.)
         if (associated(master)) then
            if (.not.master%state_indices%is_empty().or.master%fake_state_variable) then
               ! Master variable is indeed a state variable (or acts like one)
               istart = index(master%name,'/',.true.)
               select case (master%domain)
                  case (domain_bulk)
                     ! Create sum of sink-source terms for bulk variable.
                     allocate(sum)
                     sum%result_output = output_none
                     link => master%sms_list%first
                     do while (associated(link))
                        call sum%add_component(link%target%name)
                        link => link%next
                     end do
                     if (.not.sum%add_to_parent(master%owner,trim(master%name(istart+1:))//'_sms_tot')) deallocate(sum)
                  case (domain_horizontal,domain_surface,domain_bottom)
                     ! Create sum of sink-source terms for horizontal variable.
                     allocate(sum_hz)
                     sum_hz%result_output = output_none
                     link => master%sms_list%first
                     do while (associated(link))
                        call sum_hz%add_component(link%target%name)
                        link => link%next
                     end do
                     if (.not.sum_hz%add_to_parent(master%owner,trim(master%name(istart+1:))//'_sms_tot')) deallocate(sum_hz)
               end select
               master => self%find_object(trim(master%name(istart+1:))//'_sms_tot',recursive=.true.,exact=.false.)
            end if
         end if
      end if
   end function generate_master

   function generate_standard_master(self,task) result(master)
      class (type_base_model),  intent(inout),target :: self
      type (type_coupling_task),intent(inout)        :: task
      type (type_internal_variable),pointer :: master

      type (type_aggregate_variable),   pointer :: aggregate_variable

      nullify(master)
      if (task%master_standard_variable%aggregate_variable) then
         ! Make sure that an aggregate variable will be created on the fly
         select type (aggregate_standard_variable=>task%master_standard_variable)
            class is (type_bulk_standard_variable)
               aggregate_variable => get_aggregate_variable(self,aggregate_standard_variable)
               select case (task%domain)
                  case (domain_bulk)
                     aggregate_variable%bulk_required = .true.
                     deallocate(task%master_standard_variable)
                     task%master_name = aggregate_variable%standard_variable%name
                  case (domain_bottom)
                     aggregate_variable%bottom_required = .true.
                     deallocate(task%master_standard_variable)
                     task%master_name = trim(aggregate_variable%standard_variable%name)//'_at_bottom'
                  case default
                     call self%fatal_error('generate_standard_master','BUG: unknown type of standard variable with aggregate_variable set.')
               end select
         end select
      end if
   end function generate_standard_master

subroutine print_aggregate_variable_contributions(self)
   class (type_base_model),intent(inout),target :: self

   type (type_aggregate_variable),   pointer :: aggregate_variable
   type (type_contributing_variable),pointer :: contributing_variable

   aggregate_variable => self%first_aggregate_variable
   do while (associated(aggregate_variable))
      call log_message('Model '//trim(self%name)//' contributions to '//trim(aggregate_variable%standard_variable%name))
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
end subroutine

recursive subroutine build_aggregate_variables(self)
   class (type_base_model),intent(inout),target :: self

   type (type_aggregate_variable),pointer :: aggregate_variable
   type (type_model_list_node),   pointer :: child
   type (type_link),              pointer :: link
   type (type_contribution),      pointer :: contribution
   type (type_contributing_variable),pointer :: contributing_variable

   class (type_weighted_sum),           pointer :: sum
   class (type_horizontal_weighted_sum),pointer :: surface_sum,bottom_sum

   ! This routine takes the variable->aggregate variable mappings, and creates corresponding
   ! aggregate variable->variable mappings.

   ! Enumerate all model variables, and process their contributions to aggregate variables.
   link => self%links%first
   do while (associated(link))
      ! Enumerate the contributions of this variable to aggregate variables, and register these with
      ! aggregate variable objects on the model level.
      contribution => link%target%contributions%first
      do while (associated(contribution))
         aggregate_variable => get_aggregate_variable(self,contribution%target)
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

   if (self%check_conservation) then
      aggregate_variable => self%first_aggregate_variable
      do while (associated(aggregate_variable))
         if (aggregate_variable%standard_variable%conserved) then
            ! Allocate objects that will do the summation across the different domains.
            allocate(sum,surface_sum,bottom_sum)

            ! Enumerate contributions to aggregate variable.
            contributing_variable => aggregate_variable%first_contributing_variable
            do while (associated(contributing_variable))
               if (.not.contributing_variable%link%original%state_indices%is_empty().or.contributing_variable%link%original%fake_state_variable) then
                  ! Contributing variable is a state variable
                  select case (contributing_variable%link%original%domain)
                     case (domain_bulk)
                        call sum%add_component(trim(contributing_variable%link%original%name)//'_sms',contributing_variable%scale_factor)
                        call surface_sum%add_component(trim(contributing_variable%link%original%name)//'_sfl',contributing_variable%scale_factor)
                        call bottom_sum%add_component(trim(contributing_variable%link%original%name)//'_bfl',contributing_variable%scale_factor)
                     case (domain_surface)
                        call surface_sum%add_component(trim(contributing_variable%link%original%name)//'_sms',contributing_variable%scale_factor)
                     case (domain_bottom)
                        call bottom_sum%add_component(trim(contributing_variable%link%original%name)//'_sms',contributing_variable%scale_factor)
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
      call build_aggregate_variables(child%model)
      child => child%next
   end do

   ! call print_aggregate_variable_contributions(self)

end subroutine build_aggregate_variables

recursive subroutine create_aggregate_models(self)
   class (type_base_model),       intent(inout),target :: self

   type (type_aggregate_variable),      pointer :: aggregate_variable
   class (type_weighted_sum),           pointer :: sum
   class (type_horizontal_weighted_sum),pointer :: horizontal_sum,bottom_sum
   type (type_contributing_variable),   pointer :: contributing_variable
   type (type_model_list_node),         pointer :: child

   aggregate_variable => self%first_aggregate_variable
   do while (associated(aggregate_variable))
      nullify(sum,horizontal_sum,bottom_sum)
      if (aggregate_variable%bulk_required)       allocate(sum)
      if (aggregate_variable%horizontal_required) allocate(horizontal_sum)
      if (aggregate_variable%bottom_required)     allocate(bottom_sum)

      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (associated(contributing_variable%link%target,contributing_variable%link%original) &                  ! Variable must not be coupled
             .and.(associated(self%parent).or..not.contributing_variable%link%target%fake_state_variable)) then   ! Only include fake state variable for non-root models
            if (associated(sum)) then
               select case (contributing_variable%link%target%domain)
                  case (domain_bulk)
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
         if (associated(self%parent)) sum%result_output = output_none
         if (.not.sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name))) deallocate(sum)
      end if
      if (associated(horizontal_sum)) then
         horizontal_sum%units = trim(aggregate_variable%standard_variable%units)//'*m'
         if (associated(self%parent)) horizontal_sum%result_output = output_none
         if (.not.horizontal_sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name)//'_at_interfaces')) deallocate(horizontal_sum)
      end if
      if (associated(bottom_sum)) then
         bottom_sum%units = trim(aggregate_variable%standard_variable%units)//'*m'
         if (associated(self%parent)) bottom_sum%result_output = output_none
         if (.not.bottom_sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name)//'_at_bottom')) deallocate(bottom_sum)
      end if
      aggregate_variable => aggregate_variable%next
   end do

   ! Process child models
   child => self%children%first
   do while (associated(child))
      call create_aggregate_models(child%model)
      child => child%next
   end do
end subroutine create_aggregate_models

subroutine couple_variables(self,master,slave)
   class (type_base_model),     intent(inout),target :: self
   type (type_internal_variable),pointer             :: master,slave

   type (type_internal_variable),pointer :: pslave
   type (type_contribution),     pointer :: contribution

   ! If slave and master are the same, we are done - return.
   if (associated(slave,master)) return

   if (associated(self%parent)) call self%fatal_error('couple_variables','BUG: must be called on root node.')
   if (.not.slave%write_indices%is_empty()) &
      call fatal_error('couple_variables','Attempt to couple write-only variable ' &
         //trim(slave%name)//' to '//trim(master%name)//'.')
   if (.not.slave%state_indices%is_empty().and.master%state_indices%is_empty().and..not.master%fake_state_variable) &
      call fatal_error('couple_variables','Attempt to couple state variable ' &
         //trim(slave%name)//' to non-state variable '//trim(master%name)//'.')
   !if (slave%domain/=master%domain.and..not.(slave%domain==domain_horizontal.and. &
   !   (master%domain==domain_surface.or.master%domain==domain_bottom))) call fatal_error('couple_variables', &
   !   'Cannot couple '//trim(slave%name)//' to '//trim(master%name)//', because their domains are incompatible.')
   if (iand(slave%domain,master%domain)==0) call fatal_error('couple_variables', &
      'Cannot couple '//trim(slave%name)//' to '//trim(master%name)//', because their domains are incompatible.')

   if (debug_coupling) call log_message(trim(slave%name)//' --> '//trim(master%name))

   ! Merge all information from the slave into the master.
   call master%state_indices%extend(slave%state_indices)
   call master%read_indices%extend(slave%read_indices)
   call master%sms_list%extend(slave%sms_list)
   call master%background_values%extend(slave%background_values)
   call master%properties%update(slave%properties,overwrite=.false.)
   contribution => slave%contributions%first
   do while (associated(contribution))
      call master%contributions%add(contribution%target,contribution%scale_factor)
      contribution => contribution%next
   end do
   call master%surface_flux_list%extend(slave%surface_flux_list)
   call master%bottom_flux_list%extend(slave%bottom_flux_list)
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

end module
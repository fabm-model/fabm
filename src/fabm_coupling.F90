module fabm_coupling
   use fabm_types
   use fabm_standard_variables
   use fabm_properties
   use fabm_builtin_models
   use fabm_driver

   implicit none
   
   private

   public freeze_model_info, get_aggregate_variable

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

      ! Now couple model variables
      ! Stage 1: implicit - couple variables based on overlapping standard identities.
      ! Stage 2: explicit - resolve user- or model-specified links between variables.
      call couple_standard_variables(self)
      call process_coupling_tasks(self,.false.)

      ! Create models for aggregate variables at root level, to be used to compute conserved quantities.
      ! After this step, the set of variables that contribute to aggregate quatities may not be modified.
      ! That is, no new such variables may be added, and no such variables may be coupled.
      call build_aggregate_variables(self)
      call create_aggregate_models(self)

      ! Repeat coupling because new aggregate variables are now available.
      call process_coupling_tasks(self,.true.)

      ! Allow inheriting models to perform additional tasks after coupling.
      call after_coupling(self)

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
      type (type_link),pointer                 :: link,link2
      type (type_bulk_standard_variable)       :: bulk_standard_variable
      type (type_horizontal_standard_variable) :: horizontal_standard_variable
      type (type_global_standard_variable)     :: global_standard_variable
      type (type_set)                          :: processed_bulk,processed_horizontal,processed_scalar
!
!-----------------------------------------------------------------------
!BOC
      link => model%first_link
      do while (associated(link))
         select type (object=>link%target)
            class is (type_bulk_variable)
               if (.not.object%standard_variable%is_null()) then
                  if (.not.processed_bulk%contains(object%standard_variable%name)) then
                     call processed_bulk%add(object%standard_variable%name)
                     bulk_standard_variable = object%standard_variable  ! make a copy here, because object may be deallocated later from couple_variables
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_bulk_variable)
                              if (bulk_standard_variable%compare(object2%standard_variable).and..not.object2%standard_variable%is_null()) then
                                 if (object2%write_indices%is_empty()) then
                                    ! Default coupling: early variable is master, later variable is slave.
                                    call couple_variables(model,link%target,link2%target)
                                 else
                                    ! Later variable is write-only and therefore can only be master. Try copuling with early variable as slave.
                                    call couple_variables(model,link2%target,link%target)
                                 end if
                              end if
                        end select
                        link2 => link2%next
                     end do
                  end if
               end if
            class is (type_horizontal_variable)
               if (.not.object%standard_variable%is_null()) then
                  if (.not.processed_horizontal%contains(object%standard_variable%name)) then
                     call processed_horizontal%add(object%standard_variable%name)
                     horizontal_standard_variable = object%standard_variable  ! make a copy here, because object may be deallocated later from couple_variables
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_horizontal_variable)
                              if (horizontal_standard_variable%compare(object2%standard_variable).and..not.object2%standard_variable%is_null()) then
                                 if (object2%write_indices%is_empty()) then
                                    ! Default coupling: early variable is master, later variable is slave.
                                    call couple_variables(model,link%target,link2%target)
                                 else
                                    ! Later variable is write-only and therefore can only be master. Try copuling with early variable as slave.
                                    call couple_variables(model,link2%target,link%target)
                                 end if
                              end if
                        end select
                        link2 => link2%next
                     end do
                  end if
               end if
            class is (type_scalar_variable)
               if (.not.object%standard_variable%is_null()) then
                  if (.not.processed_scalar%contains(object%standard_variable%name)) then
                     call processed_scalar%add(object%standard_variable%name)
                     global_standard_variable = object%standard_variable  ! make a copy here, because object may be deallocated later from couple_variables
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_scalar_variable)
                              if (global_standard_variable%compare(object2%standard_variable).and..not.object2%standard_variable%is_null()) then
                                 if (object2%write_indices%is_empty()) then
                                    ! Default coupling: early variable is master, later variable is slave.
                                    call couple_variables(model,link%target,link2%target)
                                 else
                                    ! Later variable is write-only and therefore can only be master. Try copuling with early variable as slave.
                                    call couple_variables(model,link2%target,link%target)
                                 end if
                              end if
                        end select
                        link2 => link2%next
                     end do
                  end if
               end if
         end select
         link => link%next
      end do

      call processed_bulk%finalize()
      call processed_horizontal%finalize()
      call processed_scalar%finalize()

   end subroutine couple_standard_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process all model-specific coupling tasks.
!
! !INTERFACE:
   recursive subroutine process_coupling_tasks(self,check)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_base_model),intent(inout),target :: self
      logical,                intent(in)           :: check
!
!EOP
!
! !LOCAL VARIABLES:
      class (type_base_model),     pointer :: root
      type (type_model_list_node), pointer :: child
      type (type_link),            pointer :: link
      class (type_property),       pointer :: master_name
      class (type_internal_object),pointer :: master
!
!-----------------------------------------------------------------------
!BOC
      ! Find root model, which will handle the individual coupling tasks.
      root => self
      do while (associated(root%parent))
         root => root%parent
      end do

      ! For each variable, determine if a coupling command is provided.
      link => self%first_link
      do while (associated(link))

         ! Only process own links (those without slash in the name)
         if (index(link%name,'/')==0) then

            ! Try to find a coupling for this variable.
            master_name => self%couplings%find_in_tree(link%name)

            if (associated(master_name)) then
               select type (master_name)
                  class is (type_string_property)
                     ! Try to find the master variable among the variables of the requesting model or its parents.
                     if (link%name/=master_name%value) then
                        ! Master and slave name differ: start master search in current model, then move up tree.
                        master => self%find_object(master_name%value,recursive=.true.,exact=.false.)
                     elseif (associated(self%parent)) then
                        ! Master and slave name are identical: start master search in parent model, then move up tree.
                        master => self%parent%find_object(master_name%value,recursive=.true.,exact=.false.)
                     end if
                     if (associated(master)) then
                        ! Target variable found: perform the coupling.
                        call couple_variables(root,master,link%target)
                     elseif (check) then
                        call self%fatal_error('process_coupling_tasks', &
                           'Coupling target "'//trim(master_name%value)//'" for "'//trim(link%name)//'" was not found.')
                     end if
               end select
            end if ! If own link (no / in name)

         end if ! If coupling was specified

         link => link%next
      end do

      ! Process coupling tasks registered with child models.
      child => self%children%first
      do while (associated(child))
         call process_coupling_tasks(child%model,check)
         child => child%next
      end do

   end subroutine process_coupling_tasks
!EOC

subroutine print_aggregate_variable_contributions(self)
   class (type_base_model),intent(inout),target :: self

   type (type_aggregate_variable),   pointer :: aggregate_variable
   type (type_contributing_variable),pointer :: contributing_variable

   aggregate_variable => self%first_aggregate_variable
   do while (associated(aggregate_variable))
      call log_message('Model '//trim(self%name)//' contributions to '//trim(aggregate_variable%standard_variable%name))
      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (contributing_variable%link%owner) then
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

   integer :: nbulk,nhz

   ! This routine takes the variable->aggregate variable mappings, and creates corresponding
   ! aggregate variable->variable mappings.

   ! Enumerate all model variables, and process their contributions to aggregate variables.
   link => self%first_link
   do while (associated(link))
      select type (variable => link%target)
         class is (type_internal_variable)
            ! This link points to a variable (rather than e.g. a model dependency).
            ! Enumerate its contributions to aggregate variables, and register these with
            ! aggregate variable objects on the model level.
            contribution => variable%contributions%first
            do while (associated(contribution))
               aggregate_variable => get_aggregate_variable(self,contribution%target)
               call add_contribution(aggregate_variable,link,contribution%scale_factor,contribution%include_background)
               contribution => contribution%next
            end do
      end select
      link => link%next
   end do

   if (self%check_conservation) then
      aggregate_variable => self%first_aggregate_variable
      do while (associated(aggregate_variable))
         ! First count number of contributing state variables (separate different physical domains).
         nbulk = 0
         nhz = 0
         contributing_variable => aggregate_variable%first_contributing_variable
         do while (associated(contributing_variable))
            select type (variable=>contributing_variable%link%target)
               class is (type_bulk_variable)
                  if (.not.variable%state_indices%is_empty()) nbulk = nbulk + 1
               class is (type_horizontal_variable)
                  if (.not.variable%state_indices%is_empty()) nhz = nhz + 1
            end select
            contributing_variable => contributing_variable%next
         end do

         ! Create diagnostic variables that hold the change in aggregate variables (specific to each physical domain).
         if (nbulk>0) call self%register_diagnostic_variable(aggregate_variable%id_rate, &
                        'change_in_'//trim(aggregate_variable%standard_variable%name), &
                        trim(aggregate_variable%standard_variable%units), &
                        'change in '//trim(aggregate_variable%standard_variable%name))
         if (nhz>0) call self%register_diagnostic_variable(aggregate_variable%id_horizontal_rate, &
                      'change_in_'//trim(aggregate_variable%standard_variable%name)//'_at_interfaces', &
                      trim(aggregate_variable%standard_variable%units)//'*m', &
                      'change in '//trim(aggregate_variable%standard_variable%name)//' at interfaces')

         ! Create arrays of indices and scale factors for all contributing state variables (specific to each physical domain).
         allocate(aggregate_variable%state_indices(nbulk))
         allocate(aggregate_variable%state_scale_factors(nbulk))

         nbulk = 0
         nhz = 0
         contributing_variable => aggregate_variable%first_contributing_variable
         do while (associated(contributing_variable))
            select type (variable=>contributing_variable%link%target)
               class is (type_bulk_variable)
                  if (.not.variable%state_indices%is_empty()) then
                     nbulk = nbulk + 1
                     call variable%state_indices%append(aggregate_variable%state_indices(nbulk))
                     aggregate_variable%state_scale_factors(nbulk) = contributing_variable%scale_factor
                  end if
               class is (type_horizontal_variable)
                  if (.not.variable%state_indices%is_empty()) then
                     nhz = nhz + 1
                  end if
            end select
            contributing_variable => contributing_variable%next
         end do

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

contains

   subroutine add_contribution(aggregate_variable,link,scale_factor,include_background)
      type (type_aggregate_variable),intent(inout) :: aggregate_variable
      type (type_link),target,       intent(inout) :: link
      real(rk),                      intent(in)    :: scale_factor
      logical,                       intent(in)    :: include_background

      type (type_contributing_variable),pointer :: contributing_variable

      if (.not.associated(aggregate_variable%first_contributing_variable)) then
         ! This aggregate variable does not have any contribution registered yet. Create the first.
         allocate(aggregate_variable%first_contributing_variable)
         contributing_variable => aggregate_variable%first_contributing_variable
      else
         ! First determine whether the contribution of this variable has already been registered.
         ! This can be the case if multiple links point to the same actual variable (i.e., if they are coupled)
         contributing_variable => aggregate_variable%first_contributing_variable
         do while (associated(contributing_variable))
            if (associated(contributing_variable%link%target,link%target)) then
               ! We already have registered a contribution from this variable.
               ! The newly provided link replaces the previous if it owns the variable (i.e., it is not coupled),
               ! so we can later check the link to determine ownership. Then we are done - return.
               if (link%owner) contributing_variable%link => link
               return
            end if
            contributing_variable => contributing_variable%next
         end do

         ! This aggregate variable has one or more contributions already. Find the last, so we can append another.
         contributing_variable => aggregate_variable%first_contributing_variable
         do while (associated(contributing_variable%next))
            contributing_variable => contributing_variable%next
         end do
         allocate(contributing_variable%next)
         contributing_variable => contributing_variable%next
      end if

      ! Store contribution properties.
      contributing_variable%link => link
      contributing_variable%scale_factor = scale_factor
      contributing_variable%include_background = include_background
   end subroutine

end subroutine build_aggregate_variables

recursive subroutine create_aggregate_models(self)
   class (type_base_model),       intent(inout),target :: self

   type (type_aggregate_variable),pointer :: aggregate_variable
   type (type_model_list_node),   pointer :: child

   aggregate_variable => self%first_aggregate_variable
   do while (associated(aggregate_variable))
      if (aggregate_variable%bulk_required)       call create_aggregate_model(self,aggregate_variable,domain_bulk)
      if (aggregate_variable%horizontal_required) call create_aggregate_model(self,aggregate_variable,domain_surface)
      aggregate_variable => aggregate_variable%next
   end do

   ! Process child models
   child => self%children%first
   do while (associated(child))
      call create_aggregate_models(child%model)
      child => child%next
   end do
end subroutine

subroutine create_aggregate_model(self,aggregate_variable,domain)
   class (type_base_model),       intent(inout),target :: self
   type (type_aggregate_variable),intent(inout)        :: aggregate_variable
   integer,                       intent(in)           :: domain

   type (type_contributing_variable),pointer :: contributing_variable
   integer                                   :: n
   logical                                   :: sumrequired
   character(len=attribute_length )          :: name
   class (type_internal_object),pointer      :: target_variable
   class (type_weighted_sum),           pointer :: sum
   class (type_horizontal_weighted_sum),pointer :: horizontal_sum

   ! This procedure takes an aggregate variable, and creates models that compute diagnostics
   ! for the total of these aggregate quantities on bulk and horizontal domains.

   select case (domain)
      case (domain_bulk)
         name = trim(aggregate_variable%standard_variable%name)
         call self%add_bulk_variable(name,aggregate_variable%standard_variable%units,name)
         call self%request_coupling(name,'zero')
      case default
         name = trim(aggregate_variable%standard_variable%name)//'_at_interfaces'
         call self%add_horizontal_variable(name,aggregate_variable%standard_variable%units,name)
         call self%request_coupling(name,'zero_hz')
   end select

   n = 0
   sumrequired = .false.
   contributing_variable => aggregate_variable%first_contributing_variable
   do while (associated(contributing_variable))
      if (contributing_variable%link%owner) then
         select type (variable=>contributing_variable%link%target)
            class is (type_bulk_variable)
               if (domain==domain_bulk) then
                  n = n + 1
                  if (contributing_variable%scale_factor/=1.0_rk.or.contributing_variable%include_background) sumrequired = .true.
                  target_variable => contributing_variable%link%target
               end if
            class is (type_horizontal_variable)
               if (domain/=domain_bulk) then
                  n = n + 1
                  if (contributing_variable%scale_factor/=1.0_rk.or.contributing_variable%include_background) sumrequired = .true.
                  target_variable => contributing_variable%link%target
               end if
         end select
      end if
      contributing_variable => contributing_variable%next
   end do

   if (n==1.and..not.sumrequired) then
      ! Always equal to another - couple aggregate variable to this other
      call self%request_coupling(name,target_variable%name)
   elseif (n>0) then
      ! Weighted sum of variables
      select case (domain)
         case (domain_bulk)
            allocate(sum)
            sum%output_units = trim(aggregate_variable%standard_variable%units)
         case default
            allocate(horizontal_sum)
            horizontal_sum%output_units = trim(aggregate_variable%standard_variable%units)//'*m'
      end select

      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (contributing_variable%link%owner) then
            select type (variable=>contributing_variable%link%target)
               class is (type_bulk_variable)
                  ! This contribution comes from a bulk variable.
                  if (domain==domain_bulk) &
                     call sum%add_component(trim(contributing_variable%link%name), &
                        weight=contributing_variable%scale_factor, include_background=contributing_variable%include_background)
               class is (type_horizontal_variable)
                  ! This contribution comes from a variable defined on a horizontal interface (top or bottom).
                  if (domain/=domain_bulk) &
                     call horizontal_sum%add_component(trim(contributing_variable%link%name), &
                        weight=contributing_variable%scale_factor,include_background=contributing_variable%include_background)
            end select
         end if
         contributing_variable => contributing_variable%next
      end do

      select case (domain)
         case (domain_bulk)
            call self%add_child(sum,trim(name)//'_calculator',configunit=-1)
         case default
            call self%add_child(horizontal_sum,trim(name)//'_calculator',configunit=-1)
      end select
      call self%request_coupling(name,trim(name)//'_calculator/result')
   end if
end subroutine create_aggregate_model

subroutine couple_variables(self,master,slave)
   class (type_base_model),     intent(inout),target :: self
   class (type_internal_object),intent(inout),target :: master
   class (type_internal_object),intent(in),  pointer :: slave

   class (type_internal_object),pointer :: pslave
      
   if (associated(self%parent)) call self%fatal_error('couple_variables','BUG: must be called on root node.')

   ! Store a pointer to the slave, because the call to redirect_links will cause all pointers (from links)
   ! to the slave node to be connected to the master node. This icludes the original "slave" argument.
   pslave => slave

   ! If slave and master are the same, we are done - return.
   if (associated(pslave,master)) return

   ! Note: in call below we provide our local copy of the pointer to the slave (pslave), not the original slave pointer (slave).
   ! Reason: ifort appears to pass the slave pointer by reference. The original pointer comes from a link, and is therefore
   ! overwritten by redirect_links. If we pass this original pointer to redirect_links for comparing object identities,
   ! the recursive redirecting will fail after the very first redirect (which destroyed the pointer to the object that we
   ! want to compare against).
   call redirect_links(self,pslave,master)

   ! Merge all information from the slave into the master.
   call merge_variables(master,pslave)

   ! Deallocate the slave, which is no longer needed.
   deallocate(pslave)
end subroutine couple_variables

recursive subroutine redirect_links(model,oldtarget,newtarget)
   class (type_base_model),     intent(inout),target :: model
   class (type_internal_object),intent(in),   target :: oldtarget,newtarget

   type (type_link),           pointer :: link
   type (type_model_list_node),pointer :: child

   ! Process all links and if they used to refer to the specified slave,
   ! redirect them to the specified master.
   link => model%first_link
   do while (associated(link))
      if (associated(link%target,oldtarget)) then
         link%target => newtarget
         link%owner = .false.
      end if
      link => link%next
   end do

   ! Allow child models to do the same.
   child => model%children%first
   do while (associated(child))
      call redirect_links(child%model,oldtarget,newtarget)
      child => child%next
   end do
end subroutine

end module
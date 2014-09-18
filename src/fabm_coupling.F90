module fabm_coupling
   use fabm_types
   use fabm_standard_variables
   use fabm_properties
   use fabm_builtin_models
   use fabm_driver

   implicit none
   
   private

   public freeze_model_info, get_aggregate_variable, find_dependencies

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
      type (type_link),pointer :: link,link2
      type (type_set)          :: processed_bulk,processed_horizontal,processed_scalar
!
!-----------------------------------------------------------------------
!BOC
      link => model%links%first
      do while (associated(link))
         select type (object=>link%target)
            class is (type_bulk_variable)
               if (.not.object%standard_variable%is_null()) then
                  if (.not.processed_bulk%contains(object%standard_variable%name)) then
                     call processed_bulk%add(object%standard_variable%name)
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_bulk_variable)
                              if (object%standard_variable%compare(object2%standard_variable).and..not.object2%standard_variable%is_null()) then
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
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_horizontal_variable)
                              if (object%standard_variable%compare(object2%standard_variable).and..not.object2%standard_variable%is_null()) then
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
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_scalar_variable)
                              if (object%standard_variable%compare(object2%standard_variable).and..not.object2%standard_variable%is_null()) then
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
      class (type_base_model),       pointer :: root
      type (type_model_list_node),   pointer :: child
      type (type_link),              pointer :: link
      class (type_property),         pointer :: master_name
      class (type_internal_variable),pointer :: master
!
!-----------------------------------------------------------------------
!BOC
      ! Find root model, which will handle the individual coupling tasks.
      root => self
      do while (associated(root%parent))
         root => root%parent
      end do

      ! For each variable, determine if a coupling command is provided.
      link => self%links%first
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
                     else
                        call self%fatal_error('process_coupling_tasks', &
                           'Master and slave name are identical: "'//trim(master_name%value)//'". This is not valid at the root of the model tree.')
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
   class (type_horizontal_weighted_sum),pointer :: horizontal_sum

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
         contribution => contribution%next
         contributing_variable%next => aggregate_variable%first_contributing_variable
         aggregate_variable%first_contributing_variable => contributing_variable
      end do
      link => link%next
   end do

   if (self%check_conservation) then
      aggregate_variable => self%first_aggregate_variable
      do while (associated(aggregate_variable))
         allocate(sum)
         allocate(horizontal_sum)

         contributing_variable => aggregate_variable%first_contributing_variable
         do while (associated(contributing_variable))
            select type (variable=>contributing_variable%link%original)
               class is (type_bulk_variable)
                  if (.not.variable%state_indices%is_empty()) then
                     call sum%add_component(trim(variable%name)//'_sms',contributing_variable%scale_factor)
                     call horizontal_sum%add_component(trim(variable%name)//'_sfl',contributing_variable%scale_factor)
                     call horizontal_sum%add_component(trim(variable%name)//'_bfl',contributing_variable%scale_factor)
                  end if
               class is (type_horizontal_variable)
                  if (.not.variable%state_indices%is_empty()) &
                     call horizontal_sum%add_component(trim(variable%name)//'_sms',contributing_variable%scale_factor)
            end select
            contributing_variable => contributing_variable%next
         end do

         sum%output_units = trim(aggregate_variable%standard_variable%units)//'/s'
         if (.not.sum%add_to_parent(self,'change_in_'//trim(aggregate_variable%standard_variable%name),create_for_one=.true.)) deallocate(sum)
         horizontal_sum%output_units = trim(aggregate_variable%standard_variable%units)//'*m/s'
         if (.not.horizontal_sum%add_to_parent(self,'change_in_'//trim(aggregate_variable%standard_variable%name)//'_at_interfaces',create_for_one=.true.)) deallocate(horizontal_sum)

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
   class (type_horizontal_weighted_sum),pointer :: horizontal_sum
   type (type_contributing_variable),   pointer :: contributing_variable
   type (type_model_list_node),         pointer :: child

   aggregate_variable => self%first_aggregate_variable
   do while (associated(aggregate_variable))
      nullify(sum,horizontal_sum)
      if (aggregate_variable%bulk_required)       allocate(sum)
      if (aggregate_variable%horizontal_required) allocate(horizontal_sum)

      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (associated(contributing_variable%link%target,contributing_variable%link%original)) then
            select type (variable=>contributing_variable%link%target)
               class is (type_bulk_variable)
                  if (associated(sum)) call sum%add_component(trim(contributing_variable%link%name), &
                     weight=contributing_variable%scale_factor, include_background=contributing_variable%include_background)
               class is (type_horizontal_variable)
                  if (associated(horizontal_sum)) call horizontal_sum%add_component(trim(contributing_variable%link%name), &
                     weight=contributing_variable%scale_factor,include_background=contributing_variable%include_background)
            end select
         end if
         contributing_variable => contributing_variable%next
      end do

      if (associated(sum)) then
         sum%output_units = trim(aggregate_variable%standard_variable%units)
         if (.not.sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name))) deallocate(sum)
      end if
      if (associated(horizontal_sum)) then
         horizontal_sum%output_units = trim(aggregate_variable%standard_variable%units)//'*m'
         if (.not.horizontal_sum%add_to_parent(self,trim(aggregate_variable%standard_variable%name)//'_at_interfaces')) deallocate(horizontal_sum)
      end if
      aggregate_variable => aggregate_variable%next
   end do

   ! Process child models
   child => self%children%first
   do while (associated(child))
      call create_aggregate_models(child%model)
      child => child%next
   end do
end subroutine

subroutine couple_variables(self,master,slave)
   class (type_base_model),     intent(inout),target :: self
   class (type_internal_variable),intent(inout),target :: master
   class (type_internal_variable),intent(in),  pointer :: slave

   class (type_internal_variable),pointer :: pslave
      
   if (associated(self%parent)) call self%fatal_error('couple_variables','BUG: must be called on root node.')

   ! Store a pointer to the slave, because the call to redirect_links will cause all pointers (from links)
   ! to the slave node to be connected to the master node. This includes the original "slave" argument.
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
end subroutine couple_variables

recursive subroutine redirect_links(model,oldtarget,newtarget)
   class (type_base_model),     intent(inout),target :: model
   class (type_internal_variable),intent(in),   target :: oldtarget,newtarget

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


subroutine merge_variables(master,slave)
   class (type_internal_variable),intent(inout) :: master
   class (type_internal_variable),intent(in)    :: slave

   type (type_contribution), pointer :: contribution

   call log_message(trim(slave%name)//' --> '//trim(master%name))

   if (.not.slave%write_indices%is_empty()) &
      call fatal_error('merge_variables','Attempt to couple write-only variable ' &
         //trim(slave%name)//' to '//trim(master%name)//'.')
   if (master%state_indices%is_empty().and..not.slave%state_indices%is_empty()) &
      call fatal_error('merge_variables','Attempt to couple state variable ' &
         //trim(slave%name)//' to non-state variable '//trim(master%name)//'.')
   if (master%presence==presence_external_optional) &
      call fatal_error('merge_variables','Attempt to couple to optional master variable "'//trim(master%name)//'".')
   if (slave%domain/=master%domain) call fatal_error('merge_variables', &
      'Cannot couple '//trim(slave%name)//' to '//trim(master%name)//', because domains do not match.')

   call master%state_indices%extend(slave%state_indices,.true.)
   call master%read_indices%extend(slave%read_indices,.true.)
   call master%sms_indices%extend(slave%sms_indices,.false.)
   call master%background_values%extend(slave%background_values)
   call master%properties%update(slave%properties,overwrite=.false.)
   contribution => slave%contributions%first
   do while (associated(contribution))
      call master%contributions%add(contribution%target,contribution%scale_factor)
      contribution => contribution%next
   end do

   select type (master)
      class is (type_bulk_variable)
         select type (slave)
            class is (type_bulk_variable)
               call master%surface_flux_indices%extend(slave%surface_flux_indices,.false.)
               call master%bottom_flux_indices%extend(slave%bottom_flux_indices,.false.)
         end select      
   end select
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
            chain = trim(chain)//trim(node%model%name)//' -> '
            node => node%next
         end do
         call fatal_error('find_dependencies','circular dependency found: '//trim(chain)//trim(self%name))
      end if
      call forbidden_with_self%extend(forbidden)
   end if
   call forbidden_with_self%append(self)

   ! Loop over all variables, and if they belong to some other model, first add that model to the dependency list.
   link => self%links%first
   do while (associated(link))
      if (index(link%name,'/')==0 &                              ! Our own link...
          .and..not.link%target%write_indices%is_empty() &       ! ...to a diagnostic variable...
          .and..not.associated(link%target%owner,self) &         ! ...not set by ourselves...
          .and..not.link%original%read_indices%is_empty()) &     ! ...and we do depend on its value.
         call find_dependencies(link%target%owner,list,forbidden_with_self)
      link => link%next
   end do

   ! We're happy - add ourselves to the list of processed models.
   call list%append(self)

   ! Now process any children that have not been processed before.
   node => self%children%first
   do while (associated(node))
      call find_dependencies(node%model,list)
      node => node%next
   end do

   call forbidden_with_self%finalize()
end subroutine

end module
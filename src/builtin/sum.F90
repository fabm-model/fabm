#include "fabm_driver.h"

module fabm_builtin_sum
   use fabm_types
   use fabm_builtin_reduction
   use fabm_builtin_scale
   use fabm_builtin_source

   implicit none

   private

   public type_weighted_sum, type_horizontal_weighted_sum

   type type_sum_term
      type (type_dependency_id) :: id
      real(rk)                  :: weight = 1._rk
   end type

   type type_horizontal_sum_term
      type (type_horizontal_dependency_id) :: id
      real(rk)                             :: weight = 1._rk
   end type

   type type_component
      character(len=attribute_length) :: name   = ''
      real(rk)                        :: weight = 1._rk
      logical                         :: include_background = .false.
      type (type_dependency_id)       :: id
      type (type_component), pointer  :: next   => null()
   end type

   type type_horizontal_component
      character(len=attribute_length)           :: name   = ''
      real(rk)                                  :: weight = 1._rk
      logical                                   :: include_background = .false.
      type (type_horizontal_dependency_id)      :: id
      type (type_horizontal_component), pointer :: next   => null()
   end type

   type, extends(type_reduction_operator) :: type_weighted_sum
      character(len=attribute_length) :: units         = ''
      integer                         :: result_output = output_instantaneous
      real(rk)                        :: offset        = 0.0_rk
      real(rk)                        :: missing_value = -2.e20_rk
      integer                         :: access        = access_read
      class (type_interior_standard_variable), pointer :: standard_variable => null()
      type (type_add_id)              :: id_output
      type (type_component), pointer  :: first => null()
      type (type_sum_term), allocatable :: sources(:)
   contains
      procedure :: initialize       => weighted_sum_initialize
      procedure :: add_component    => weighted_sum_add_component
      procedure :: do               => weighted_sum_do
      procedure :: after_coupling   => weighted_sum_after_coupling
      procedure :: add_to_parent    => weighted_sum_add_to_parent
      procedure :: merge_components => weighted_sum_merge_components
      procedure :: finalize         => weighted_sum_finalize
   end type

   type, extends(type_reduction_operator) :: type_horizontal_weighted_sum
      character(len=attribute_length) :: units         = ''
      integer                         :: result_output = output_instantaneous
      real(rk)                        :: offset        = 0.0_rk
      real(rk)                        :: missing_value = -2.e20_rk
      integer                         :: access        = access_read
      integer                         :: domain        = domain_horizontal
      class (type_horizontal_standard_variable), pointer :: standard_variable => null()
      type (type_horizontal_add_id)   :: id_output
      type (type_horizontal_component), pointer   :: first => null()
      type (type_horizontal_sum_term), allocatable :: sources(:)
   contains
      procedure :: add_component    => horizontal_weighted_sum_add_component
      procedure :: initialize       => horizontal_weighted_sum_initialize
      procedure :: do_horizontal    => horizontal_weighted_sum_do_horizontal
      procedure :: after_coupling   => horizontal_weighted_sum_after_coupling
      procedure :: add_to_parent    => horizontal_weighted_sum_add_to_parent
      procedure :: merge_components => horizontal_weighted_sum_merge_components
      procedure :: finalize         => horizontal_weighted_sum_finalize
   end type

   type, extends(type_base_model) :: type_weighted_sum_sms_distributor
      type (type_dependency_id)                  :: id_total_sms
      type (type_state_variable_id), allocatable :: id_targets(:)
      real(rk),allocatable                       :: weights(:)
   contains
      !procedure :: do => weighted_sum_sms_distributor_do
   end type

contains

   function weighted_sum_add_to_parent(self, parent, name, create_for_one, aggregate_variable, link) result(sum_used)
      class (type_weighted_sum),              intent(inout), target :: self
      class (type_base_model),                intent(inout), target :: parent
      character(len=*),                       intent(in)            :: name
      logical,optional,                       intent(in)            :: create_for_one
      type (type_interior_standard_variable), intent(in), optional  :: aggregate_variable
      type (type_link), pointer,                          optional  :: link

      logical                                        :: sum_used, create_for_one_
      type (type_link),                      pointer :: link_
      class (type_scaled_interior_variable), pointer :: scaled_variable

      create_for_one_ = .false.
      if (present(create_for_one)) create_for_one_ = create_for_one

      sum_used = .false.
      link_ => null()
      if (associated(self%standard_variable)) then
         call parent%add_interior_variable(name, self%units, name, link=link_, act_as_state_variable=iand(self%access,access_set_source) /= 0, standard_variable=self%standard_variable)
      else
         call parent%add_interior_variable(name, self%units, name, link=link_, act_as_state_variable=iand(self%access,access_set_source) /= 0)
      end if
      if (present(link)) link => link_
      if (.not. associated(self%first)) then
         ! No components - link to constant field with offset (typically 0)
         link_%target%source = source_constant
         link_%target%prefill_value = self%offset
         link_%target%missing_value = self%missing_value
         link_%target%output = self%result_output
      elseif (.not. associated(self%first%next)) then
         ! One component only.
         if (self%first%weight == 1.0_rk .and. .not. create_for_one_) then
            ! One component with scale factor 1 - directly link to the component's source variable.
            call parent%request_coupling(link_, self%first%name)
         else
            ! One component with scale factor other than 1 (or a user-specified requirement NOT to make a direct link to the source variable)
            allocate(scaled_variable)
            call parent%add_child(scaled_variable, trim(name) // '_calculator')
            call scaled_variable%register_dependency(scaled_variable%id_source, 'source', self%units, 'source variable')
            call scaled_variable%request_coupling(scaled_variable%id_source, self%first%name)
            call scaled_variable%register_diagnostic_variable(scaled_variable%id_result, 'result', self%units, 'result', &
               missing_value=self%missing_value, output=self%result_output, act_as_state_variable=iand(self%access, access_set_source) /= 0)
            scaled_variable%weight = self%first%weight
            scaled_variable%include_background = self%first%include_background
            scaled_variable%offset = self%offset
            call parent%request_coupling(link_, trim(name)//'_calculator/result')
            if (iand(self%access, access_set_source) /= 0) then
               ! This scaled variable acts as a state variable. Create a child model to distribute source terms to the original source variable.
               call copy_fluxes(scaled_variable, scaled_variable%id_result, self%first%name, scale_factor=1.0_rk / scaled_variable%weight)
               if (present(aggregate_variable)) call scaled_variable%add_to_aggregate_variable(aggregate_variable, scaled_variable%id_result)
            end if
         end if
         deallocate(self%first)
      else
         ! Multiple components. Create the sum.
         call parent%add_child(self, trim(name) // '_calculator')
         call parent%request_coupling(link_, trim(name) // '_calculator/result')
         sum_used = .true.
      end if
   end function weighted_sum_add_to_parent

   subroutine weighted_sum_initialize(self, configunit)
      class (type_weighted_sum), intent(inout), target :: self
      integer,                   intent(in)            :: configunit

      type (type_component), pointer :: component
      integer           :: ncomponents, i
      character(len=10) :: temp
      class (type_weighted_sum_sms_distributor), pointer :: sms_distributor

      call self%register_implemented_routines((/source_do/))

      call self%get_parameter(ncomponents, 'n', '', 'number of terms in summation', default=0, minimum=0)
      do i = 1, ncomponents
         call self%add_component('')
      end do
      call self%get_parameter(self%units, 'units', '', 'units', default=trim(self%units))

      ncomponents = 0
      component => self%first
      do while (associated(component))
         ncomponents = ncomponents + 1
         write (temp,'(i0)') ncomponents
         call self%get_parameter(component%weight, 'weight' // trim(temp), '-', 'weight for term ' // trim(temp), default=component%weight)
         call self%register_dependency(component%id, 'term' // trim(temp), self%units, 'term ' // trim(temp))
         if (component%name /= '') call self%request_coupling(component%id, trim(component%name))
         component => component%next
      end do

      !call self%register_diagnostic_variable(self%id_output,'result',self%units,'result',output=self%result_output) !,act_as_state_variable=iand(self%access,access_set_source)/=0)
      call self%add_interior_variable('result', self%units, 'result', fill_value=0.0_rk, missing_value=self%missing_value, &
         output=self%result_output, write_index=self%id_output%sum_index, link=self%id_output%link, source=source_do)

      if (iand(self%access, access_set_source) /= 0) then
         ! NB this does not function yet (hence the commented out act_as_state_variable above)
         ! Auto-generation of result_sms_tot fails and the do routine of type_weighted_sum_sms_distributor is not yet implemented.

         ! The sum will act as a state variable. Any source terms will have to be distributed over the individual variables that contribute to the sum.
         allocate(sms_distributor)
         call self%add_child(sms_distributor, 'sms_distributor')
         call sms_distributor%register_dependency(sms_distributor%id_total_sms, 'total_sms', trim(self%units) // '/s', 'sources-sinks of sum')
         call sms_distributor%request_coupling(sms_distributor%id_total_sms, 'result_sms_tot')
         allocate(sms_distributor%weights(ncomponents))
         allocate(sms_distributor%id_targets(ncomponents))

         ncomponents = 0
         component => self%first
         do while (associated(component))
            ncomponents = ncomponents + 1
            write (temp,'(i0)') ncomponents
            call sms_distributor%register_state_dependency(sms_distributor%id_targets(ncomponents), 'target'//trim(temp), self%units, 'target '//trim(temp))
            call sms_distributor%request_coupling(sms_distributor%id_targets(ncomponents), trim(component%name))
            sms_distributor%weights(ncomponents) = component%weight
            component => component%next
         end do
      end if
   end subroutine

   subroutine weighted_sum_add_component(self, name, weight, include_background)
      class (type_weighted_sum),intent(inout) :: self
      character(len=*),         intent(in)    :: name
      real(rk), optional,       intent(in)    :: weight
      logical,  optional,       intent(in)    :: include_background

      type (type_component), pointer :: component

      if (_VARIABLE_REGISTERED_(self%id_output)) &
         call self%fatal_error('weighted_sum_add_component', 'cannot be called after model initialization')

      if (.not. associated(self%first)) then
         allocate(self%first)
         component => self%first
      else
         component => self%first
         do while (associated(component%next))
            component => component%next
         end do
         allocate(component%next)
         component => component%next
      end if
      component%name = name
      if (present(weight)) component%weight = weight
      if (present(include_background)) component%include_background = include_background
   end subroutine

   subroutine weighted_sum_after_coupling(self)
      class (type_weighted_sum), intent(inout) :: self

      type (type_component), pointer :: component
      real(rk)                       :: background

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      background = 0
      component => self%first
      do while (associated(component))
         if (component%include_background) then
            self%offset = self%offset + component%weight * component%id%background
         else
            background = background + component%weight * component%id%background
         end if
         component => component%next
      end do
      call self%id_output%link%target%background_values%set_value(background)
   end subroutine

   subroutine weighted_sum_merge_components(self, log_unit)
      class (type_weighted_sum), intent(inout) :: self
      integer, optional,         intent(in)    :: log_unit

      integer                                :: i, n
      type (type_internal_variable), pointer :: sum_variable
      type (type_component),         pointer :: component, component_next, component_previous

      sum_variable => self%id_output%link%target
      sum_variable%prefill_value = sum_variable%prefill_value + self%offset
      allocate(sum_variable%cowriters)
      if (present(log_unit)) write (log_unit,'(a)') 'Reindexing ' // trim(sum_variable%name)
      component_previous => null()
      component => self%first
      n = 0
      do while (associated(component))
         component_next => component%next
         if (merge_component(component%id%link, component%weight, sum_variable, log_unit)) then
            ! Component was merged into target
            if (associated(component_previous)) then
               component_previous%next => component_next
            else
               self%first => component_next
            end if
            deallocate(component)
         else
            ! Component was left stand-alone (not merged)
            component_previous => component
            n = n + 1
         end if
         component => component_next
      end do

      if (n == 0) then
         if (associated(sum_variable%cowriters%first)) then
            if (present(log_unit)) write (log_unit,'(a)') '- all contributions written in place - skipping actual sum'
         else
            if (present(log_unit)) write (log_unit,'(a,g0.6)') '- no contributions left - sum becomes a constant with value ', sum_variable%prefill_value
            sum_variable%source = source_constant
         end if
         return
      end if

      call sum_variable%cowriters%add(sum_variable)

      ! Put all ids and weights in a single array to minimize memory bottlenecks when computing sum
      allocate(self%sources(n))
      component => self%first
      do i = 1, n
         self%sources(i)%weight = component%weight
         call component%id%link%target%read_indices%append(self%sources(i)%id%index)
         component => component%next
      end do
   end subroutine weighted_sum_merge_components

   subroutine weighted_sum_finalize(self)
      class (type_weighted_sum), intent(inout) :: self

      type (type_component), pointer :: component, component_next

      component => self%first
      do while (associated(component))
         component_next => component%next
         deallocate(component)
         component => component_next
      end do
      self%first => null()
      call self%type_reduction_operator%finalize()
   end subroutine weighted_sum_finalize

   logical function merge_component(component_link, weight, target_variable, log_unit)
      type (type_link),              intent(inout)         :: component_link
      type (type_internal_variable), intent(inout), target :: target_variable
      real(rk),                      intent(in)            :: weight
      integer, optional,             intent(in)            :: log_unit

      type (type_internal_variable), pointer :: component_variable

      component_variable => component_link%target

      merge_component = .false.
      if (iand(component_variable%write_operator, operator_merge_forbidden) /= 0) then
         if (present(log_unit)) write (log_unit,'(a)') '- kept ' // trim(component_variable%name) // ' because it is accessed by FABM or host'
      elseif (component_variable%write_operator /= operator_add .and. component_variable%source /= source_constant) then
         if (present(log_unit)) write (log_unit,'(a)') '- kept ' // trim(component_variable%name) // ' because it assigns directly'
      elseif (weight /= 1.0_rk) then
         if (present(log_unit)) write (log_unit,'(a,g0.6)') '- kept '// trim(component_variable%name) // ' because it uses scale factor ', weight
      elseif (size(component_variable%read_indices%pointers) > 1) then
         if (present(log_unit)) write (log_unit,'(a)') '- kept ' // trim(component_variable%name) // ' because it is accessed by other models'
      else
         ! This component does not need a separate diagnostic - it will be merged into the target
         target_variable%prefill_value = target_variable%prefill_value + component_variable%prefill_value
         if (component_variable%source == source_constant) then
            ! This component is a constant that we can just add to our offset
            if (present(log_unit)) write (log_unit,'(a,g0.6,a)') '- merged ' // trim(component_variable%name) // ' - constant ', component_variable%prefill_value, ' added to offset'
         else
            ! This component can increment the sum result directly
            if (present(log_unit)) write (log_unit,'(a,g0.6,a)') '- merged ' // trim(component_variable%name) // ' - to be written in-place'
            component_variable%write_owner => target_variable
            call target_variable%cowriters%add(component_variable)
         end if
         call component_variable%read_indices%set_value(-1)
         call component_variable%read_indices%finalize()
         component_link%original%read_index => null()
         merge_component = .true.
      end if
   end function

   subroutine horizontal_weighted_sum_merge_components(self, log_unit)
      class (type_horizontal_weighted_sum), intent(inout) :: self
      integer, optional,                    intent(in)    :: log_unit

      integer                                   :: i, n
      type (type_internal_variable),    pointer :: sum_variable
      type (type_horizontal_component), pointer :: component, component_next, component_previous

      sum_variable => self%id_output%link%target
      sum_variable%prefill_value = sum_variable%prefill_value + self%offset
      allocate(sum_variable%cowriters)
      if (present(log_unit)) write (log_unit,'(a)') 'Reindexing ' // trim(sum_variable%name)
      component_previous => null()
      component => self%first
      n = 0
      do while (associated(component))
         component_next => component%next
         if (merge_component(component%id%link, component%weight, sum_variable, log_unit)) then
            ! Component was merged into target
            if (associated(component_previous)) then
               component_previous%next => component_next
            else
               self%first => component_next
            end if
            deallocate(component)
         else
            ! Component was left stand-alone (not merged)
            component_previous => component
            n = n + 1
         end if
         component => component_next
      end do

      if (n == 0) then
         if (associated(sum_variable%cowriters%first)) then
            if (present(log_unit)) write (log_unit,'(a)') '- all contributions written in place - skipping actual sum'
         else
            if (present(log_unit)) write (log_unit,'(a,g0.6)') '- no contributions left - sum becomes a constant with value ', sum_variable%prefill_value
            sum_variable%source = source_constant
         end if
         return
      end if

      call sum_variable%cowriters%add(sum_variable)

      ! Put all ids and weights in a single array to minimize memory bottlenecks when computing sum
      allocate(self%sources(n))
      component => self%first
      do i = 1, n
         self%sources(i)%weight = component%weight
         call component%id%link%target%read_indices%append(self%sources(i)%id%horizontal_index)
         component => component%next
      end do
   end subroutine horizontal_weighted_sum_merge_components

   subroutine horizontal_weighted_sum_finalize(self)
      class (type_horizontal_weighted_sum), intent(inout) :: self

      type (type_horizontal_component), pointer :: component, component_next

      component => self%first
      do while (associated(component))
         component_next => component%next
         deallocate(component)
         component => component_next
      end do
      self%first => null()
      call self%type_reduction_operator%finalize()
   end subroutine horizontal_weighted_sum_finalize

   subroutine weighted_sum_do(self, _ARGUMENTS_DO_)
      class (type_weighted_sum), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      integer  :: i
      real(rk) :: value

      ! Enumerate components included in the sum, and add their contributions.
      do i = 1, size(self%sources)
         _CONCURRENT_LOOP_BEGIN_
            _GET_(self%sources(i)%id, value)
            _ADD_(self%id_output, self%sources(i)%weight * value)
         _LOOP_END_
      end do
   end subroutine

   function horizontal_weighted_sum_add_to_parent(self, parent, name, create_for_one, aggregate_variable, link) result(sum_used)
      class (type_horizontal_weighted_sum),      intent(inout), target :: self
      class (type_base_model),                   intent(inout), target :: parent
      character(len=*),                          intent(in)            :: name
      logical,optional,                          intent(in)            :: create_for_one
      class (type_horizontal_standard_variable), intent(in), optional  :: aggregate_variable
      type (type_link), pointer,                             optional  :: link

      logical :: sum_used, create_for_one_
      type (type_link), pointer :: link_
      class (type_scaled_horizontal_variable), pointer :: scaled_variable

      create_for_one_ = .false.
      if (present(create_for_one)) create_for_one_ = create_for_one

      sum_used = .false.
      link_ => null()
      if (associated(self%standard_variable)) then
         call parent%add_horizontal_variable(name, self%units, name, link=link_, act_as_state_variable=iand(self%access, access_set_source) /= 0, standard_variable=self%standard_variable)
      else
         call parent%add_horizontal_variable(name, self%units, name, link=link_, act_as_state_variable=iand(self%access, access_set_source) /= 0)
      end if
      if (present(link)) link => link_
      if (.not.associated(self%first)) then
         ! No components - link to constant field with offset (typically 0)
         link_%target%source = source_constant
         link_%target%prefill_value = self%offset
         link_%target%missing_value = self%missing_value
         link_%target%output = self%result_output
      elseif (.not. associated(self%first%next)) then
         ! One component only.
         if (self%first%weight == 1.0_rk .and. .not. create_for_one_) then
            ! One component with scale factor 1 - directly link to the component's source variable.
            call parent%request_coupling(link_, self%first%name)
         else
            ! One component with scale factor other than 1 (or a user-specified requirement NOT to make a direct link to the source variable)
            allocate(scaled_variable)
            call parent%add_child(scaled_variable, trim(name) // '_calculator')
            call scaled_variable%register_dependency(scaled_variable%id_source, 'source', self%units, 'source variable')
            call scaled_variable%request_coupling(scaled_variable%id_source, self%first%name)
            call scaled_variable%register_diagnostic_variable(scaled_variable%id_result, 'result', self%units, 'result', &
               missing_value=self%missing_value, output=self%result_output, act_as_state_variable= &
               iand(self%access, access_set_source) /= 0, source=source_do_horizontal, domain=self%domain)
            scaled_variable%weight = self%first%weight
            scaled_variable%include_background = self%first%include_background
            scaled_variable%offset = self%offset
            call parent%request_coupling(link_, trim(name)//'_calculator/result')
            if (iand(self%access, access_set_source) /= 0) then
               call copy_horizontal_fluxes(scaled_variable, scaled_variable%id_result, self%first%name, scale_factor=1.0_rk / scaled_variable%weight)
               if (present(aggregate_variable)) call scaled_variable%add_to_aggregate_variable(aggregate_variable, scaled_variable%id_result)
            end if
         end if
         deallocate(self%first)
      else
         ! One component with scale factor unequal to 1, or multiple components. Create the sum.
         call parent%add_child(self, trim(name) // '_calculator')
         call parent%request_coupling(link_, trim(name) // '_calculator/result')
         sum_used = .true.
      end if
   end function horizontal_weighted_sum_add_to_parent

   subroutine horizontal_weighted_sum_initialize(self, configunit)
      class (type_horizontal_weighted_sum), intent(inout), target :: self
      integer,                              intent(in)            :: configunit

      type (type_horizontal_component),pointer :: component
      integer           :: i,n
      character(len=10) :: temp

      call self%register_implemented_routines((/source_do_horizontal/))

      call self%get_parameter(n,'n','','number of terms in summation',default=0,minimum=0)
      do i=1,n
         call self%add_component('')
      end do
      call self%get_parameter(self%units,'units','','units',default=trim(self%units))

      i = 0
      component => self%first
      do while (associated(component))
         i = i + 1
         write (temp,'(i0)') i
         call self%get_parameter(component%weight,'weight'//trim(temp),'-','weight for term '//trim(temp),default=component%weight)
         call self%register_dependency(component%id,'term'//trim(temp),self%units,'term '//trim(temp))
         call self%request_coupling(component%id,trim(component%name))
         component => component%next
      end do
      call self%add_horizontal_variable('result', self%units, 'result', missing_value=self%missing_value, fill_value=0.0_rk, output=self%result_output, &
         write_index=self%id_output%horizontal_sum_index, link=self%id_output%link, source=source_do_horizontal)
   end subroutine

   subroutine horizontal_weighted_sum_after_coupling(self)
      class (type_horizontal_weighted_sum),intent(inout) :: self

      type (type_horizontal_component),pointer :: component
      real(rk) :: background

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      background = 0
      component => self%first
      do while (associated(component))
         if (component%include_background) then
            self%offset = self%offset + component%weight*component%id%background
         else
            background = background + component%weight*component%id%background
         end if
         component => component%next
      end do
      call self%id_output%link%target%background_values%set_value(background)
   end subroutine

   subroutine horizontal_weighted_sum_add_component(self,name,weight,include_background)
      class (type_horizontal_weighted_sum),intent(inout) :: self
      character(len=*),                    intent(in)    :: name
      real(rk),optional,                   intent(in)    :: weight
      logical,optional,                    intent(in)    :: include_background

      type (type_horizontal_component),pointer :: component

      if (_VARIABLE_REGISTERED_(self%id_output)) &
         call self%fatal_error('weighted_sum_add_component','cannot be called after model initialization')

      if (.not.associated(self%first)) then
         allocate(self%first)
         component => self%first
      else
         component => self%first
         do while (associated(component%next))
            component => component%next
         end do
         allocate(component%next)
         component => component%next
      end if
      component%name = name
      if (present(weight)) component%weight = weight
      if (present(include_background)) component%include_background = include_background
   end subroutine

   subroutine horizontal_weighted_sum_do_horizontal(self,_ARGUMENTS_HORIZONTAL_)
      class (type_horizontal_weighted_sum), intent(in) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_

      integer  :: i
      real(rk) :: value

      do i = 1, size(self%sources)
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            _GET_HORIZONTAL_(self%sources(i)%id, value)
            _ADD_HORIZONTAL_(self%id_output, self%sources(i)%weight * value)
         _HORIZONTAL_LOOP_END_
      end do
   end subroutine horizontal_weighted_sum_do_horizontal

end module

#include "fabm_driver.h"

module fabm_builtin_sum
   use fabm_types
   use fabm_builtin_reduction
   use fabm_builtin_scale

   implicit none

   private

   public type_base_sum, type_weighted_sum, type_horizontal_weighted_sum

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
      type (type_link), pointer       :: link => null()
      real(rk)                        :: weight = 1._rk
      logical                         :: include_background = .false.
      type (type_component), pointer  :: next   => null()
   end type

   type, extends(type_reduction_operator) :: type_base_sum
      character(len=attribute_length) :: units         = ''
      integer                         :: result_output = output_instantaneous
      real(rk)                        :: offset        = 0.0_rk
      real(rk)                        :: missing_value = -2.e20_rk
      logical                         :: act_as_state_variable = .false.
      type (type_link), pointer       :: result_link => null()

      logical, private :: components_frozen = .false.
      logical, private :: active = .false.
      type (type_component), pointer, private  :: first => null()
   contains
      procedure :: add_component_by_name
      procedure :: add_component_by_link
      generic :: add_component => add_component_by_name, add_component_by_link
      procedure :: after_coupling
      procedure :: finalize
   end type

   type, extends(type_base_sum) :: type_weighted_sum
      type (type_interior_standard_variable), pointer :: aggregate_variable => null()
      type (type_dependency_id), allocatable :: id_terms(:)
      type (type_add_id)                     :: id_result
      type (type_sum_term), private, allocatable :: sources(:)
   contains
      procedure :: initialize       => weighted_sum_initialize
      procedure :: do               => weighted_sum_do
      procedure :: merge_components => weighted_sum_merge_components
   end type

   type, extends(type_base_sum) :: type_horizontal_weighted_sum
      class (type_horizontal_standard_variable), pointer :: aggregate_variable => null()
      integer                                           :: domain = domain_horizontal
      type (type_horizontal_dependency_id), allocatable :: id_terms(:)
      type (type_horizontal_add_id)                     :: id_result
      type (type_horizontal_sum_term), private, allocatable :: sources(:)
   contains
      procedure :: initialize       => horizontal_weighted_sum_initialize
      procedure :: do_horizontal    => horizontal_weighted_sum_do_horizontal
      procedure :: merge_components => horizontal_weighted_sum_merge_components
   end type

   type, extends(type_base_model) :: type_weighted_sum_sms_distributor
      type (type_dependency_id)                  :: id_total_sms
      type (type_state_variable_id), allocatable :: id_targets(:)
      real(rk),allocatable                       :: weights(:)
   contains
      !procedure :: do => weighted_sum_sms_distributor_do
   end type

contains

   subroutine request_coupling_to_component(parent, target_link, component)
      class (type_base_model), intent(inout), target :: parent
      type (type_link), target :: target_link
      class (type_component), intent(in) :: component

      if (associated(component%link)) then
         call parent%request_coupling(target_link, component%link)
      elseif (component%name /= '') then
         call parent%request_coupling(target_link, component%name)
      end if
   end subroutine request_coupling_to_component

   function base_initialize(self) result(n)
      class (type_base_sum), intent(inout) :: self
      integer                              :: n

      integer                        :: i
      type (type_component), pointer :: component
      character(len=10)              :: temp

      call self%get_parameter(n, 'n', '', 'number of terms in summation', default=0, minimum=0)
      do i = 1, n
         call self%add_component('')
      end do
      call self%get_parameter(self%units, 'units', '', 'units', default=trim(self%units))

      n = 0
      component => self%first
      do while (associated(component))
         n = n + 1
         write (temp,'(i0)') n
         call self%get_parameter(component%weight, 'weight' // trim(temp), '-', 'weight for term ' // trim(temp), default=component%weight)
         component => component%next
      end do

      self%components_frozen = .true.
   end function base_initialize

   function append_component(self, weight, include_background) result(component)
      class (type_base_sum), intent(inout) :: self
      real(rk), optional,    intent(in)    :: weight
      logical,  optional,    intent(in)    :: include_background

      type (type_component), pointer :: component

      if (self%components_frozen) &
         call self%fatal_error('base_sum_add_component', 'cannot be called after model initialization')

      if (.not. associated(self%first)) then
         ! List is empty - add first component
         allocate(self%first)
         component => self%first
      else
         ! List is not empty - find last component and append new one
         component => self%first
         do while (associated(component%next))
            component => component%next
         end do
         allocate(component%next)
         component => component%next
      end if
      if (present(weight)) component%weight = weight
      if (present(include_background)) component%include_background = include_background
   end function append_component

   subroutine add_component_by_name(self, name, weight, include_background)
      class (type_base_sum),    intent(inout) :: self
      character(len=*),         intent(in)    :: name
      real(rk), optional,       intent(in)    :: weight
      logical,  optional,       intent(in)    :: include_background

      type (type_component), pointer :: component

      component => append_component(self, weight, include_background)
      component%name = name
   end subroutine add_component_by_name

   subroutine add_component_by_link(self, link, weight, include_background)
      class (type_base_sum),    intent(inout) :: self
      type (type_link), target                :: link
      real(rk), optional,       intent(in)    :: weight
      logical,  optional,       intent(in)    :: include_background

      type (type_component), pointer :: component

      component => append_component(self, weight, include_background)
      component%link => link
   end subroutine add_component_by_link

   subroutine after_coupling(self)
      class (type_base_sum), intent(inout) :: self

      type (type_component), pointer :: component
      real(rk)                       :: background
      integer                        :: i

      ! If we are not using the actual sum (because we have 0 or 1 components), skip this routine altogether
      if (.not. self%active) return

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      background = 0.0_rk
      i = 0
      component => self%first
      do while (associated(component))
         i = i + 1
         if (component%include_background) then
            self%offset = self%offset + component%weight * component%link%target%background_values%value
         else
            background = background + component%weight * component%link%target%background_values%value
         end if
         component => component%next
      end do
      call self%result_link%target%background_values%set_value(background)
   end subroutine after_coupling

   function base_merge_components(self, log_unit) result(n)
      class (type_base_sum),                 intent(inout) :: self
      integer, optional,                     intent(in)    :: log_unit
      integer :: n

      type (type_internal_variable), pointer :: sum_variable
      type (type_component),         pointer :: component, component_next, component_previous

      n = 0

      ! If we are not using the actual sum (because we have 0 or 1 components), skip this routine altogether
      if (.not. self%active) return

      sum_variable => self%result_link%target

      ! Include offset in the prefill value that will be used to initialize the sum
      sum_variable%prefill_value = sum_variable%prefill_value + self%offset

      allocate(sum_variable%cowriters)
      if (present(log_unit)) write (log_unit,'(a)') 'Reindexing ' // trim(sum_variable%name)
      component_previous => null()
      component => self%first
      do while (associated(component))
         component_next => component%next
         if (merge_component(component)) then
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

   contains

      logical function merge_component(component)
         type (type_component), intent(in) :: component

         type (type_internal_variable), pointer :: component_variable

         component_variable => component%link%target

         merge_component = .false.
         if (iand(component_variable%write_operator, operator_merge_forbidden) /= 0) then
            if (present(log_unit)) write (log_unit,'(a)') '- kept ' // trim(component_variable%name) // ' because it is accessed by FABM or host'
         elseif (component_variable%write_operator /= operator_add .and. component_variable%source /= source_constant) then
            if (present(log_unit)) write (log_unit,'(a)') '- kept ' // trim(component_variable%name) // ' because it assigns directly'
         elseif (component%weight /= 1.0_rk) then
            if (present(log_unit)) write (log_unit,'(a,g0.6)') '- kept '// trim(component_variable%name) // ' because it uses scale factor ', component%weight
         elseif (size(component_variable%read_indices%pointers) > 1) then
            if (present(log_unit)) write (log_unit,'(a)') '- kept ' // trim(component_variable%name) // ' because it is accessed by other models'
         else
            ! This component does not need a separate diagnostic - it will be merged into the target
            sum_variable%prefill_value = sum_variable%prefill_value + component_variable%prefill_value
            if (component_variable%source == source_constant) then
               ! This component is a constant that we can just add to our offset
               if (present(log_unit)) write (log_unit,'(a,g0.6,a)') '- merged ' // trim(component_variable%name) // ' - constant ', component_variable%prefill_value, ' added to offset'
            else
               ! This component can increment the sum result directly
               if (present(log_unit)) write (log_unit,'(a,g0.6,a)') '- merged ' // trim(component_variable%name) // ' - to be written in-place'
               component_variable%write_owner => sum_variable
               call sum_variable%cowriters%add(component_variable)
            end if
            call component_variable%read_indices%set_value(-1)
            call component_variable%read_indices%finalize()
            component%link%original%read_index => null()
            merge_component = .true.
         end if
      end function

   end function base_merge_components

   recursive subroutine finalize(self)
      class (type_base_sum), intent(inout) :: self

      type (type_component), pointer :: component, component_next

      component => self%first
      do while (associated(component))
         component_next => component%next
         deallocate(component)
         component => component_next
      end do
      self%first => null()
      call self%type_reduction_operator%finalize()
   end subroutine finalize

   subroutine weighted_sum_initialize(self, configunit)
      class (type_weighted_sum), intent(inout), target :: self
      integer,                   intent(in)            :: configunit

      type (type_component), pointer :: component
      integer           :: i, n, output
      character(len=10) :: temp
      class (type_weighted_sum_sms_distributor), pointer :: sms_distributor
      class (type_scaled_interior_variable), pointer :: scaled_variable

      n = base_initialize(self)

      allocate(self%id_terms(n))

      component => self%first
      do i = 1, n
         write (temp,'(i0)') i
         call self%register_dependency(self%id_terms(i), 'term' // trim(temp), self%units, 'term ' // trim(temp))
         call request_coupling_to_component(self, self%id_terms(i)%link, component)
         component%link => self%id_terms(i)%link
         component => component%next
      end do

      if (n == 0) then
         ! No components - create constant field with offset (typically 0)
         call self%add_interior_variable('result', self%units, 'result', fill_value=self%offset, missing_value=self%missing_value, &
            output=self%result_output, source=source_constant, link=self%result_link)
         return
      elseif (n == 1) then
         ! One component only - set up an alias that will be linked to the actual values
         output = self%result_output
         if (output /= output_none) output = ior(output, output_always_available)
         call self%add_interior_variable('result', self%units, 'result', link=self%result_link, &
            output=output, act_as_state_variable=self%act_as_state_variable, presence=presence_external_required)
         if (self%first%weight == 1.0_rk) then
            ! One component with scale factor 1 - directly link to the component's source variable.
            call request_coupling_to_component(self, self%result_link, self%first)
         else
            ! One component with scale factor other than 1 - add a child model to perform the scaling
            allocate(scaled_variable)
            scaled_variable%units = self%units
            scaled_variable%weight = self%first%weight
            scaled_variable%offset = self%offset
            scaled_variable%include_background = self%first%include_background
            scaled_variable%act_as_state_variable = self%act_as_state_variable
            scaled_variable%missing_value = self%missing_value
            scaled_variable%result_output = output_none
            call self%add_child(scaled_variable, '*')
            call request_coupling_to_component(scaled_variable, scaled_variable%id_source%link, self%first)
            call self%request_coupling(self%result_link, scaled_variable%id_result%link)
            if (self%act_as_state_variable .and. associated(self%aggregate_variable)) call scaled_variable%add_to_aggregate_variable(self%aggregate_variable, scaled_variable%id_result)
         end if
         return
      end if

      self%active = .true.

      call self%add_interior_variable('result', self%units, 'result', fill_value=0.0_rk, missing_value=self%missing_value, &
         output=self%result_output, write_index=self%id_result%sum_index, link=self%id_result%link, source=source_do, &
         act_as_state_variable=self%act_as_state_variable)
      self%result_link => self%id_result%link

      if (self%act_as_state_variable) then
         ! NB this does not function yet (hence the act_as_state_variable=.false. above)
         ! Auto-generation of result_sms_tot fails and the do routine of type_weighted_sum_sms_distributor is not yet implemented.

         ! The sum will act as a state variable. Any source terms will have to be distributed over the individual variables that contribute to the sum.
         allocate(sms_distributor)
         call self%add_child(sms_distributor, 'sms_distributor')
         call sms_distributor%register_dependency(sms_distributor%id_total_sms, 'total_sms', trim(self%units) // '/s', 'sources-sinks of sum')
         call sms_distributor%request_coupling(sms_distributor%id_total_sms, 'result_sms_tot')
         allocate(sms_distributor%weights(n))
         allocate(sms_distributor%id_targets(n))

         component => self%first
         do i = 1, n
            write (temp,'(i0)') i
            call sms_distributor%register_state_dependency(sms_distributor%id_targets(i), 'target' // trim(temp), self%units, 'target ' // trim(temp))
            call request_coupling_to_component(sms_distributor, sms_distributor%id_targets(i)%link, component)
            sms_distributor%weights(i) = component%weight
            component => component%next
         end do
      end if
   end subroutine weighted_sum_initialize

   subroutine weighted_sum_merge_components(self, log_unit)
      class (type_weighted_sum), intent(inout) :: self
      integer, optional,         intent(in)    :: log_unit

      integer                        :: i, n
      type (type_component), pointer :: component

      n = base_merge_components(self, log_unit)

      ! Put all ids and weights in a single array to minimize memory bottlenecks when computing sum
      allocate(self%sources(n))
      component => self%first
      do i = 1, n
         self%sources(i)%weight = component%weight
         call component%link%target%read_indices%append(self%sources(i)%id%index)
         component => component%next
      end do
   end subroutine weighted_sum_merge_components

   subroutine weighted_sum_do(self, _ARGUMENTS_DO_)
      class (type_weighted_sum), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      integer  :: i
      real(rk) :: value

      ! Enumerate components included in the sum, and add their contributions.
      ! Note: the offset (if any) is already included in the prefill value of id_result.
      do i = 1, size(self%sources)
         _CONCURRENT_LOOP_BEGIN_
            _GET_(self%sources(i)%id, value)
            _ADD_(self%id_result, self%sources(i)%weight * value)
         _LOOP_END_
      end do
   end subroutine weighted_sum_do

   subroutine horizontal_weighted_sum_initialize(self, configunit)
      class (type_horizontal_weighted_sum), intent(inout), target :: self
      integer,                              intent(in)            :: configunit

      type (type_component), pointer :: component
      integer           :: i, n, output
      character(len=10) :: temp
      class (type_scaled_horizontal_variable), pointer :: scaled_variable

      n = base_initialize(self)

      allocate(self%id_terms(n))

      component => self%first
      do i = 1, n
         write (temp,'(i0)') i
         call self%register_dependency(self%id_terms(i), 'term' // trim(temp), self%units, 'term ' // trim(temp))
         call request_coupling_to_component(self, self%id_terms(i)%link, component)
         component%link => self%id_terms(i)%link
         component => component%next
      end do

      if (n == 0) then
         ! No components - create constant field with offset (typically 0)
         call self%add_horizontal_variable('result', self%units, 'result', fill_value=self%offset, missing_value=self%missing_value, &
            output=self%result_output, source=source_constant, link=self%result_link, domain=self%domain)
         return
      elseif (n == 1) then
         ! One component only - set up an alias that will be linked to the actual values
         output = self%result_output
         if (output /= output_none) output = ior(output, output_always_available)
         call self%add_horizontal_variable('result', self%units, 'result', link=self%result_link, domain=self%domain, &
            output=output, act_as_state_variable=self%act_as_state_variable, presence=presence_external_required)
         if (self%first%weight == 1.0_rk) then
            ! One component with scale factor 1 - directly link to the component's source variable.
            call request_coupling_to_component(self, self%result_link, self%first)
         else
            ! One component with scale factor other than 1 - add a child model to perform the scaling
            allocate(scaled_variable)
            scaled_variable%units = self%units
            scaled_variable%weight = self%first%weight
            scaled_variable%offset = self%offset
            scaled_variable%domain = self%domain
            scaled_variable%include_background = self%first%include_background
            scaled_variable%act_as_state_variable = self%act_as_state_variable
            scaled_variable%missing_value = self%missing_value
            scaled_variable%result_output = output_none
            call self%add_child(scaled_variable, '*')
            call request_coupling_to_component(scaled_variable, scaled_variable%id_source%link, self%first)
            call self%request_coupling(self%result_link, scaled_variable%id_result%link)
            if (self%act_as_state_variable .and. associated(self%aggregate_variable)) call scaled_variable%add_to_aggregate_variable(self%aggregate_variable, scaled_variable%id_result)
         end if
         return
      end if

      self%active = .true.

      call self%add_horizontal_variable('result', self%units, 'result', missing_value=self%missing_value, fill_value=0.0_rk, output=self%result_output, &
         write_index=self%id_result%horizontal_sum_index, link=self%id_result%link, source=source_do_horizontal, domain=self%domain)
      self%result_link => self%id_result%link
   end subroutine horizontal_weighted_sum_initialize

   subroutine horizontal_weighted_sum_merge_components(self, log_unit)
      class (type_horizontal_weighted_sum), intent(inout) :: self
      integer, optional,                    intent(in)    :: log_unit

      integer                        :: i, n
      type (type_component), pointer :: component

      n = base_merge_components(self, log_unit)

      ! Put all ids and weights in a single array to minimize memory bottlenecks when computing sum
      allocate(self%sources(n))
      component => self%first
      do i = 1, n
         self%sources(i)%weight = component%weight
         call component%link%target%read_indices%append(self%sources(i)%id%horizontal_index)
         component => component%next
      end do
   end subroutine horizontal_weighted_sum_merge_components

   subroutine horizontal_weighted_sum_do_horizontal(self,_ARGUMENTS_HORIZONTAL_)
      class (type_horizontal_weighted_sum), intent(in) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_

      integer  :: i
      real(rk) :: value

      ! Enumerate components included in the sum, and add their contributions.
      ! Note: the offset (if any) is already included in the prefill value of id_result.
      do i = 1, size(self%sources)
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            _GET_HORIZONTAL_(self%sources(i)%id, value)
            _ADD_HORIZONTAL_(self%id_result, self%sources(i)%weight * value)
         _HORIZONTAL_LOOP_END_
      end do
   end subroutine horizontal_weighted_sum_do_horizontal

end module

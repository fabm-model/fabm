#include "fabm_driver.h"

module fabm_builtin_models
   use fabm_types
   use fabm_standard_variables

   implicit none

   private

   public type_reduction_operator, type_weighted_sum, type_horizontal_weighted_sum
   public type_depth_integral, type_bounded_depth_integral
   public copy_fluxes, copy_horizontal_fluxes
   public type_surface_source

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: builtin_factory

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
      type (type_component),pointer   :: next   => null()
   end type

   type type_horizontal_component
      character(len=attribute_length)          :: name   = ''
      real(rk)                                 :: weight = 1._rk
      logical                                  :: include_background = .false.
      type (type_horizontal_dependency_id)     :: id
      type (type_horizontal_component),pointer :: next   => null()
   end type

   type, extends(type_base_model) :: type_reduction_operator
   contains
      procedure :: merge_components => reduction_operator_merge_components
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
   end type

   type, extends(type_base_model) :: type_weighted_sum_sms_distributor
      type (type_dependency_id)                  :: id_total_sms
      type (type_state_variable_id), allocatable :: id_targets(:)
      real(rk),allocatable                       :: weights(:)
   contains
      !procedure :: do => weighted_sum_sms_distributor_do
   end type

   type, extends(type_base_model) :: type_scaled_interior_variable
      type (type_dependency_id)          :: id_source
      type (type_diagnostic_variable_id) :: id_result
      real(rk)                           :: weight = 1.0_rk
      real(rk)                           :: offset = 0.0_rk
      logical                            :: include_background = .false.
   contains
      procedure :: initialize     => scaled_interior_variable_initialize
      procedure :: do             => scaled_interior_variable_do
      procedure :: after_coupling => scaled_interior_variable_after_coupling
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
   end type

   type, extends(type_base_model) :: type_scaled_horizontal_variable
      type (type_horizontal_dependency_id)          :: id_source
      type (type_horizontal_diagnostic_variable_id) :: id_result
      real(rk)                                      :: weight = 1.0_rk
      real(rk)                                      :: offset = 0.0_rk
      logical                                       :: include_background = .false.
   contains
      procedure :: initialize     => scaled_horizontal_variable_initialize
      procedure :: do_horizontal  => scaled_horizontal_variable_do_horizontal
      procedure :: after_coupling => scaled_horizontal_variable_after_coupling
   end type

   type, extends(type_base_model) :: type_depth_integral
      type (type_dependency_id)                     :: id_input
      type (type_dependency_id)                     :: id_thickness
      type (type_horizontal_diagnostic_variable_id) :: id_output
      logical                                       :: average = .false.
   contains
      procedure :: initialize     => depth_integral_initialize
      procedure :: do_column      => depth_integral_do_column
      procedure :: after_coupling => depth_integral_after_coupling
   end type

   type, extends(type_depth_integral) :: type_bounded_depth_integral
      real(rk) :: minimum_depth = 0.0_rk
      real(rk) :: maximum_depth = huge(1.0_rk)
   contains
      procedure :: do_column => bounded_depth_integral_do_column
   end type

   type, extends(type_base_model) :: type_interior_constant
      type (type_diagnostic_variable_id) :: id_constant
   contains
      procedure :: initialize => interior_constant_initialize
   end type

   type, extends(type_base_model) :: type_horizontal_constant
      type (type_horizontal_diagnostic_variable_id) :: id_constant
   contains
      procedure :: initialize => horizontal_constant_initialize
   end type

   type, extends(type_base_model) :: type_horizontal_layer
      type (type_dependency_id)                     :: id_source
      type (type_horizontal_diagnostic_variable_id) :: id_result
   contains
      procedure :: after_coupling  => horizontal_layer_after_coupling
   end type

   type,extends(type_horizontal_layer) :: type_bottom_field
   contains
      procedure :: initialize => bottom_field_initialize
      procedure :: do_bottom  => bottom_field_do_bottom
   end type

   type, extends(type_horizontal_layer) :: type_surface_field
   contains
      procedure :: initialize => surface_field_initialize
      procedure :: do_surface => surface_field_do_surface
   end type

   type, extends(type_base_model) :: type_constant_surface_flux
      type (type_state_variable_id) :: id_target
      real(rk) :: flux
   contains
      procedure :: initialize => constant_surface_flux_initialize
      procedure :: do_surface => constant_surface_flux_do_surface
   end type

   type, extends(type_base_model) :: type_horizontal_flux
      type (type_link), pointer            :: target => null()
      type (type_horizontal_add_id)        :: id_target_flux
      type (type_horizontal_dependency_id) :: id_flux
      real(rk)                             :: scale_factor = 1.0_rk
   contains
      procedure :: do_horizontal => horizontal_flux_do_horizontal
   end type

   type, extends(type_horizontal_flux) :: type_external_surface_flux
   contains
      procedure :: initialize => external_surface_flux_initialize
   end type

   type, extends(type_horizontal_flux) :: type_external_bottom_flux
   contains
      procedure :: initialize => external_bottom_flux_initialize
   end type

   type, extends(type_base_model) :: type_interior_source
      type (type_link), pointer :: target => null()
      type (type_add_id)        :: id_target_sms
      type (type_dependency_id) :: id_source
      real(rk)                  :: scale_factor = 1.0_rk
   contains
      procedure :: initialize => interior_source_initialize
      procedure :: do         => interior_source_do
   end type

   type, extends(type_horizontal_flux) :: type_bottom_source
   contains
      procedure :: initialize => bottom_source_initialize
   end type

   type, extends(type_horizontal_flux) :: type_surface_source
   contains
      procedure :: initialize => surface_source_initialize
   end type

   type, extends(type_base_model) :: type_interior_relaxation
      type (type_state_variable_id) :: id_original
      type (type_dependency_id)     :: id_target
      type (type_dependency_id)     :: id_rate
      real(rk) :: rate
      logical :: rate_is_variable
   contains
      procedure :: initialize => interior_relaxation_initialize
      procedure :: do         => interior_relaxation_do
   end type

   type, extends(type_base_model) :: type_column_projection
      type (type_horizontal_dependency_id) :: id_source
      type (type_diagnostic_variable_id)   :: id_result
   contains
      procedure :: initialize => column_projection_initialize
      procedure :: do         => column_projection_do
   end type

   interface copy_fluxes
      module procedure copy_fluxes_to_id
      module procedure copy_fluxes_to_named_variable
   end interface

   contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('bulk_constant');          allocate(type_interior_constant::model)
         case ('interior_constant');      allocate(type_interior_constant::model)
         case ('horizontal_constant');    allocate(type_horizontal_constant::model)
         case ('surface_flux');           allocate(type_constant_surface_flux::model)
         case ('constant_surface_flux');  allocate(type_constant_surface_flux::model)
         case ('external_surface_flux');  allocate(type_external_surface_flux::model)
         case ('external_bottom_flux');   allocate(type_external_bottom_flux::model)
         case ('interior_source');        allocate(type_interior_source::model)
         case ('bottom_source');          allocate(type_bottom_source::model)
         case ('surface_source');         allocate(type_surface_source::model)
         case ('interior_relaxation');    allocate(type_interior_relaxation::model)
         case ('column_projection');      allocate(type_column_projection::model)
         case ('weighted_sum');           allocate(type_weighted_sum::model)
         case ('horizontal_weighted_sum');allocate(type_horizontal_weighted_sum::model)
         case ('bottom_layer');           allocate(type_bottom_field::model)
         case ('surface_layer');          allocate(type_surface_field::model)
         case ('vertical_integral');      allocate(type_depth_integral::model)
         ! Add new examples models here
      end select

   end subroutine

   subroutine reduction_operator_merge_components(self, log_unit)
      class (type_reduction_operator), intent(inout) :: self
      integer, optional,               intent(in)    :: log_unit
   end subroutine

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
            call parent%add_child(scaled_variable, trim(name) // '_calculator', configunit=-1)
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
         call parent%add_child(self, trim(name) // '_calculator', configunit=-1)
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
         call self%add_child(sms_distributor,'sms_distributor',configunit=-1)
         call sms_distributor%register_dependency(sms_distributor%id_total_sms,'total_sms',trim(self%units)//'/s','sources-sinks of sum')
         call sms_distributor%request_coupling(sms_distributor%id_total_sms,'result_sms_tot')
         allocate(sms_distributor%weights(ncomponents))
         allocate(sms_distributor%id_targets(ncomponents))

         ncomponents = 0
         component => self%first
         do while (associated(component))
            ncomponents = ncomponents + 1
            write (temp,'(i0)') ncomponents
            call sms_distributor%register_state_dependency(sms_distributor%id_targets(ncomponents),'target'//trim(temp),self%units,'target '//trim(temp))
            call sms_distributor%request_coupling(sms_distributor%id_targets(ncomponents),trim(component%name))
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

   subroutine scaled_interior_variable_initialize(self, configunit)
      class (type_scaled_interior_variable), intent(inout), target :: self
      integer,                               intent(in)            :: configunit
      call self%register_implemented_routines((/source_do/))
   end subroutine

   subroutine scaled_interior_variable_do(self, _ARGUMENTS_DO_)
      class (type_scaled_interior_variable), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value

      _CONCURRENT_LOOP_BEGIN_
         _GET_(self%id_source,value)
         _SET_DIAGNOSTIC_(self%id_result, self%offset + self%weight*value)
      _LOOP_END_
   end subroutine scaled_interior_variable_do

   subroutine scaled_interior_variable_after_coupling(self)
      class (type_scaled_interior_variable), intent(inout) :: self

      if (self%include_background) then
         self%offset = self%offset + self%weight * self%id_source%background
      else
         call self%id_result%link%target%background_values%set_value(self%weight * self%id_source%background)
      end if
   end subroutine scaled_interior_variable_after_coupling

   subroutine scaled_horizontal_variable_initialize(self, configunit)
      class (type_scaled_horizontal_variable), intent(inout), target :: self
      integer,                                 intent(in)            :: configunit
      call self%register_implemented_routines((/source_do_horizontal/))
   end subroutine

   subroutine scaled_horizontal_variable_do_horizontal(self, _ARGUMENTS_HORIZONTAL_)
      class (type_scaled_horizontal_variable), intent(in) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_

      real(rk) :: value

      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_source, value)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_result, self%offset + self%weight * value)
      _HORIZONTAL_LOOP_END_
   end subroutine scaled_horizontal_variable_do_horizontal

   subroutine scaled_horizontal_variable_after_coupling(self)
      class (type_scaled_horizontal_variable), intent(inout) :: self

      if (self%include_background) then
         self%offset = self%offset + self%weight * self%id_source%background
      else
         call self%id_result%link%target%background_values%set_value(self%weight * self%id_source%background)
      end if
   end subroutine scaled_horizontal_variable_after_coupling

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
            call parent%add_child(scaled_variable, trim(name) // '_calculator', configunit=-1)
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
         call parent%add_child(self, trim(name) // '_calculator', configunit=-1)
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

   subroutine depth_integral_initialize(self, configunit)
      class (type_depth_integral), intent(inout), target :: self
      integer,                     intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_column/))
      call self%register_dependency(self%id_input, 'source', '', 'source')
      call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
      call self%register_diagnostic_variable(self%id_output, 'result', '', 'result', source=source_do_column)
   end subroutine depth_integral_initialize

   subroutine depth_integral_after_coupling(self)
      class (type_depth_integral), intent(inout) :: self

      if (associated(self%id_output%link%target, self%id_output%link%original)) self%id_output%link%target%units = trim(self%id_input%link%target%units)//'*m'
   end subroutine depth_integral_after_coupling

   subroutine depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_depth_integral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: h,value,result,depth

      result = 0
      depth = 0
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_thickness, h)
         _GET_(self%id_input, value)
         depth = depth + h
         result = result + h * value
      _VERTICAL_LOOP_END_

      if (self%average) result = result / depth
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output, result)
   end subroutine depth_integral_do_column

   subroutine bounded_depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_bounded_depth_integral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: h,value,cum,depth
      logical :: started

      cum = 0
      depth = 0
      started = self%minimum_depth <= 0
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_thickness,h)
         _GET_(self%id_input,value)

         depth = depth + h
         if (.not.started) then
            ! Not yet at minimum depth before
            if (depth >= self%minimum_depth) then
               ! Now crossing minimum depth interface
               started = .true.
               h = depth-self%minimum_depth
            end if
         elseif (depth > self%maximum_depth) then
            ! Now crossing maximum depth interface; subtract part of layer height that is not included
            h = h - (depth - self%maximum_depth)
            _VERTICAL_LOOP_EXIT_
         end if
         cum = cum + h*value
      _VERTICAL_LOOP_END_

      if (.not.self%average) then
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,cum)
      elseif (depth>self%minimum_depth) then
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,cum/(min(self%maximum_depth,depth)-self%minimum_depth))
      endif
   end subroutine bounded_depth_integral_do_column

   subroutine interior_constant_initialize(self,configunit)
      class (type_interior_constant),intent(inout),target :: self
      integer,                       intent(in)           :: configunit

      character(len=attribute_length) :: standard_name
      real(rk)                        :: value

      call self%register_implemented_routines()
      call self%get_parameter(standard_name,'standard_name','','standard name',default='')
      call self%get_parameter(value,'value','','value')
      if (standard_name/='') then
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, standard_variable=type_interior_standard_variable(name=standard_name), source=source_constant)
      else
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, source=source_constant)
      end if
   end subroutine interior_constant_initialize

   subroutine horizontal_constant_initialize(self,configunit)
      class (type_horizontal_constant),intent(inout),target :: self
      integer,                         intent(in)           :: configunit

      character(len=attribute_length) :: standard_name
      real(rk)                        :: value

      call self%register_implemented_routines()
      call self%get_parameter(standard_name,'standard_name','','standard name',default='')
      call self%get_parameter(value,'value','','value')
      if (standard_name/='') then
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, standard_variable=type_horizontal_standard_variable(name=standard_name), source=source_constant)
      else
         call self%register_diagnostic_variable(self%id_constant,'data','','data', missing_value=value, &
            output=output_none, source=source_constant)
      end if
   end subroutine horizontal_constant_initialize

   subroutine horizontal_layer_after_coupling(self)
      class (type_horizontal_layer), intent(inout) :: self

      if (associated(self%id_result%link%target, self%id_result%link%original)) self%id_result%link%target%units = self%id_source%link%target%units
   end subroutine horizontal_layer_after_coupling

   subroutine bottom_field_initialize(self,configunit)
      class (type_bottom_field), intent(inout), target :: self
      integer,                   intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_bottom/))
      call self%register_diagnostic_variable(self%id_result, 'result', '', 'bottom values', source=source_do_bottom)
      call self%register_dependency(self%id_source, 'source', '', 'interior values')
   end subroutine bottom_field_initialize

   subroutine bottom_field_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_bottom_field), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: value

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_source, value)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_result, value)
      _HORIZONTAL_LOOP_END_
   end subroutine bottom_field_do_bottom

   subroutine surface_field_initialize(self,configunit)
      class (type_surface_field), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_surface/))
      call self%register_diagnostic_variable(self%id_result, 'result', '', 'surface values', source=source_do_surface)
      call self%register_dependency(self%id_source, 'source', '', 'interior values')
   end subroutine surface_field_initialize

   subroutine surface_field_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_surface_field), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: value

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_source, value)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_result, value)
      _HORIZONTAL_LOOP_END_
   end subroutine surface_field_do_surface

   subroutine constant_surface_flux_initialize(self,configunit)
      class (type_constant_surface_flux), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_surface/))
      call self%register_state_dependency(self%id_target,'target','UNITS m-3','target variable')
      call self%get_parameter(self%flux,'flux','UNITS m-2 s-1','flux (positive for into water)')
   end subroutine constant_surface_flux_initialize

   subroutine constant_surface_flux_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_constant_surface_flux), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      _HORIZONTAL_LOOP_BEGIN_
         _SET_SURFACE_EXCHANGE_(self%id_target,self%flux)
      _HORIZONTAL_LOOP_END_
   end subroutine constant_surface_flux_do_surface

   subroutine horizontal_flux_do_horizontal(self,_ARGUMENTS_HORIZONTAL_)
      class (type_horizontal_flux), intent(in) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_

      real(rk) :: flux

      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_flux,flux)
         _ADD_HORIZONTAL_(self%id_target_flux, self%scale_factor * flux)
      _HORIZONTAL_LOOP_END_
   end subroutine horizontal_flux_do_horizontal

   subroutine external_surface_flux_initialize(self,configunit)
      class (type_external_surface_flux),intent(inout),target :: self
      integer,                           intent(in)           :: configunit

      call self%register_implemented_routines((/source_do_horizontal/))
      call self%add_interior_variable('target', 'UNITS m-3', 'target variable', link=self%target)
      call self%register_surface_flux(self%target, self%id_target_flux, source=source_do_horizontal)
      call self%register_dependency(self%id_flux,'flux','UNITS m-2 s-1','surface flux')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine external_surface_flux_initialize

   subroutine external_bottom_flux_initialize(self,configunit)
      class (type_external_bottom_flux),intent(inout),target :: self
      integer,                          intent(in)           :: configunit

      call self%register_implemented_routines((/source_do_horizontal/))
      call self%add_interior_variable('target', 'UNITS m-3', 'target variable', link=self%target)
      call self%register_bottom_flux(self%target, self%id_target_flux, source=source_do_horizontal)
      call self%register_dependency(self%id_flux,'flux','UNITS m-2 s-1','bottom flux')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine external_bottom_flux_initialize

   subroutine interior_source_initialize(self,configunit)
      class (type_interior_source),intent(inout),target :: self
      integer,                     intent(in)           :: configunit

      call self%register_implemented_routines((/source_do/))
      call self%add_interior_variable('target', 'UNITS m-3', 'target variable', link=self%target)
      call self%register_source(self%target, self%id_target_sms)
      call self%register_dependency(self%id_source, 'source', 'UNITS m-3 s-1', 'source')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine interior_source_initialize

   subroutine interior_source_do(self,_ARGUMENTS_DO_)
      class (type_interior_source), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: source

      _LOOP_BEGIN_
         _GET_(self%id_source,source)
         _ADD_(self%id_target_sms, self%scale_factor * source)
      _LOOP_END_
   end subroutine interior_source_do

   subroutine bottom_source_initialize(self,configunit)
      class (type_bottom_source),intent(inout),target :: self
      integer,                   intent(in)           :: configunit

      call self%register_implemented_routines((/source_do_horizontal/))
      call self%add_horizontal_variable('target', 'UNITS m-2', 'target variable', domain=domain_bottom, link=self%target)
      call self%register_bottom_source(self%target, self%id_target_flux, source=source_do_horizontal)
      call self%register_dependency(self%id_flux,'source','UNITS m-2 s-1','source')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine bottom_source_initialize

   subroutine surface_source_initialize(self, configunit)
      class (type_surface_source), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_horizontal/))
      call self%add_horizontal_variable('target', 'UNITS m-2', 'target variable', domain=domain_bottom, link=self%target)
      call self%register_surface_source(self%target, self%id_target_flux, source=source_do_horizontal)
      call self%register_dependency(self%id_flux,'source','UNITS m-2 s-1','source')
      call self%get_parameter(self%scale_factor, 'scale_factor', '', 'scale factor', default=1.0_rk)
   end subroutine surface_source_initialize

   subroutine interior_relaxation_initialize(self,configunit)
      class (type_interior_relaxation),intent(inout),target :: self
      integer,                         intent(in)           :: configunit

      call self%register_implemented_routines((/source_do/))
      call self%register_state_dependency(self%id_original,'original', '', 'variable that is to be relaxed')
      call self%register_dependency(self%id_target,'target', '', 'target to relax towards')
      call self%get_parameter(self%rate_is_variable, 'rate_is_variable', '', 'use variable relaxation rate', default=.false.)
      if (self%rate_is_variable) then
         call self%register_dependency(self%id_rate, 'rate', 's-1', 'relaxation rate')
      else
         call self%get_parameter(self%rate, 'rate', 's-1', 'relaxation rate', minimum=0._rk)
      end if
   end subroutine interior_relaxation_initialize

   subroutine interior_relaxation_do(self,_ARGUMENTS_DO_)
      class (type_interior_relaxation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: original_value, target_value, relaxation_rate

      relaxation_rate = self%rate
      _LOOP_BEGIN_
         _GET_(self%id_original,original_value)
         _GET_(self%id_target,target_value)
         if (self%rate_is_variable) _GET_(self%id_rate, relaxation_rate)
         _SET_ODE_(self%id_original, relaxation_rate*(target_value - original_value))
      _LOOP_END_
   end subroutine interior_relaxation_do

   subroutine column_projection_initialize(self,configunit)
      class (type_column_projection),intent(inout),target :: self
      integer,                       intent(in)           :: configunit

      call self%register_implemented_routines((/source_do/))
      call self%register_dependency(self%id_source,'source', '', 'horizontal source')
      call self%register_diagnostic_variable(self%id_result,'result', '', 'interior result')
   end subroutine column_projection_initialize

   subroutine column_projection_do(self,_ARGUMENTS_DO_)
      class (type_column_projection), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value

      _LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_source,value)
         _SET_DIAGNOSTIC_(self%id_result,value)
      _LOOP_END_
   end subroutine column_projection_do

   subroutine copy_fluxes_to_id(source_model,source_variable,target_variable,scale_factor)
      class (type_base_model),           intent(inout), target :: source_model
      type (type_diagnostic_variable_id),intent(in)            :: source_variable
      type (type_state_variable_id),     intent(in)            :: target_variable
      real(rk),optional,                 intent(in)            :: scale_factor
      call copy_fluxes_to_named_variable(source_model,source_variable,target_variable%link%target%name,scale_factor)
   end subroutine

   subroutine copy_fluxes_to_named_variable(source_model,source_variable,target_variable,scale_factor)
      class (type_base_model),            intent(inout), target :: source_model
      type (type_diagnostic_variable_id), intent(in)            :: source_variable
      character(len=*),                   intent(in)            :: target_variable
      real(rk),optional,                  intent(in)            :: scale_factor

      class (type_interior_source),       pointer :: copier
      class (type_external_bottom_flux),  pointer :: bfl_copier
      class (type_external_surface_flux), pointer :: sfl_copier

      allocate(copier)
      call source_model%add_child(copier, '*', configunit=-1)
      if (present(scale_factor)) copier%scale_factor = scale_factor
      call copier%request_coupling(copier%target, target_variable)
      call copier%request_coupling(copier%id_source, trim(source_variable%link%target%name)//'_sms_tot')

      allocate(bfl_copier)
      call source_model%add_child(bfl_copier, '*', configunit=-1)
      if (present(scale_factor)) bfl_copier%scale_factor = scale_factor
      call bfl_copier%request_coupling(bfl_copier%target, target_variable)
      call bfl_copier%request_coupling(bfl_copier%id_flux, trim(source_variable%link%target%name)//'_bfl_tot')

      allocate(sfl_copier)
      call source_model%add_child(sfl_copier, '*', configunit=-1)
      if (present(scale_factor)) sfl_copier%scale_factor = scale_factor
      call sfl_copier%request_coupling(sfl_copier%target, target_variable)
      call sfl_copier%request_coupling(sfl_copier%id_flux, trim(source_variable%link%target%name)//'_sfl_tot')
   end subroutine

   subroutine copy_horizontal_fluxes(source_model,source_variable,target_variable,scale_factor)
      class (type_base_model),                       intent(inout), target :: source_model
      type (type_horizontal_diagnostic_variable_id), intent(in)            :: source_variable
      character(len=*),                              intent(in)            :: target_variable
      real(rk),optional,                             intent(in)            :: scale_factor

      class (type_bottom_source),  pointer :: bottom_copier
      class (type_surface_source), pointer :: surface_copier

      select case (source_variable%link%target%domain)
      case (domain_bottom)
         allocate(bottom_copier)
         call source_model%add_child(bottom_copier, '*', configunit=-1)
         if (present(scale_factor)) bottom_copier%scale_factor = scale_factor
         call bottom_copier%request_coupling(bottom_copier%target, target_variable)
         call bottom_copier%request_coupling(bottom_copier%id_flux, trim(source_variable%link%target%name)//'_sms_tot')
      case (domain_surface)
         allocate(surface_copier)
         call source_model%add_child(surface_copier, '*', configunit=-1)
         if (present(scale_factor)) surface_copier%scale_factor = scale_factor
         call surface_copier%request_coupling(surface_copier%target, target_variable)
         call surface_copier%request_coupling(surface_copier%id_flux, trim(source_variable%link%target%name)//'_sms_tot')
      case default
         call source_model%fatal_error('copy_horizontal_fluxes','source variable has unknown domain (should be either surface or bottom)')
      end select
   end subroutine copy_horizontal_fluxes

end module fabm_builtin_models

#include "fabm_driver.h"
#include "fabm_private.h"

! =============================================================================
! FABM --- Framework for Aquatic Biogeochemical Models
! -----------------------------------------------------------------------------
! This is the core module of FABM, serving as the "glue layer" between a
! physical host model (e.g., a general circulation model), and one or more
! specific biogeochemical models. A physical host model will call the
! interfaces of this module to access biogeochemistry.
!
! For more information, see the documentation at http://fabm.net/wiki.
!
! To add new biogeochemical models, add source code under src/models and
! reference your institute in src/CMakeLists.txt
! =============================================================================

module fabm

   use fabm_parameters
   use fabm_types, rki => rk, fabm_standard_variables => standard_variables
   use fabm_expressions
   use fabm_driver
   use fabm_properties
   use fabm_builtin_models
   use fabm_coupling
   use fabm_job
   use fabm_schedule
   use fabm_debug
   use fabm_work
   use fabm_config

   implicit none

   private

   ! --------------------------------------------------------------------------
   ! Public members
   ! --------------------------------------------------------------------------

   public fabm_initialize_library
   public fabm_get_version
   public fabm_create_model
   public type_fabm_model

   ! Variable identifier types by external physical drivers.
   public type_fabm_variable_id
   public type_fabm_interior_variable_id
   public type_fabm_horizontal_variable_id
   public type_fabm_scalar_variable_id
   public type_fabm_variable, type_fabm_interior_state_variable, type_fabm_horizontal_state_variable, &
      type_fabm_interior_diagnostic_variable, type_fabm_horizontal_diagnostic_variable

   ! Object with all supported standard variables as its members.
   ! Imported from fabm_types, and made available so hosts only need to "use fabm"
   public fabm_standard_variables

   integer, parameter :: status_none             = 0
   integer, parameter :: status_initialize_done  = 1
   integer, parameter, public :: status_set_domain_done  = 2
   integer, parameter, public :: status_start_done = 3

   integer, parameter, public :: data_source_none = 0
   integer, parameter, public :: data_source_host = 1
   integer, parameter, public :: data_source_fabm = 2
   integer, parameter, public :: data_source_user = 3
   integer, parameter, public :: data_source_default = data_source_host

   ! --------------------------------------------------------------------------
   ! Derived typed for variable identifiers
   ! --------------------------------------------------------------------------

   type type_fabm_variable_id
      type (type_internal_variable), pointer :: variable => null()
   end type

   type, extends(type_fabm_variable_id) :: type_fabm_interior_variable_id
   end type

   type, extends(type_fabm_variable_id) :: type_fabm_horizontal_variable_id
   end type

   type, extends(type_fabm_variable_id) :: type_fabm_scalar_variable_id
   end type

   ! --------------------------------------------------------------------------
   ! Derived types for variable metadata
   ! --------------------------------------------------------------------------

   ! Derived type for metadata of a generic variable (base type)
   type, abstract :: type_fabm_variable
      character(len=attribute_length) :: name          = ''
      character(len=attribute_length) :: long_name     = ''
      character(len=attribute_length) :: local_long_name = ''
      character(len=attribute_length) :: units         = ''
      character(len=attribute_length) :: path          = ''
      real(rke)                       :: minimum       = -1.e20_rke
      real(rke)                       :: maximum       =  1.e20_rke
      real(rke)                       :: missing_value = -2.e20_rke
      integer                         :: output        = output_instantaneous ! See output_* parameters defined in fabm_types
      type (type_property_dictionary) :: properties
      integer                         :: externalid    = 0                    ! Identifier to be used freely by host
      type (type_internal_variable), pointer :: target => null()
   end type

   ! Derived type for interior state variable metadata
   type, extends(type_fabm_variable) :: type_fabm_interior_state_variable
      class (type_interior_standard_variable), pointer :: standard_variable => null()
      real(rke)                                        :: initial_value             = 0.0_rke
      logical                                          :: no_precipitation_dilution = .false.
      logical                                          :: no_river_dilution         = .false.
   end type

   ! Derived type for horizontal (bottom/surface) state variable metadata
   type, extends(type_fabm_variable) :: type_fabm_horizontal_state_variable
      class (type_horizontal_standard_variable), pointer :: standard_variable => null()
      real(rke)                                          :: initial_value = 0.0_rke
   end type

   ! Derived type for interior diagnostic variable metadata
   type, extends(type_fabm_variable) :: type_fabm_interior_diagnostic_variable
      class (type_interior_standard_variable), pointer :: standard_variable => null()
      logical                                          :: save = .false.
      integer                                          :: source
   end type

   ! Derived type for horizontal diagnostic variable metadata
   type, extends(type_fabm_variable) :: type_fabm_horizontal_diagnostic_variable
      class (type_horizontal_standard_variable), pointer :: standard_variable => null()
      logical                                            :: save = .false.
      integer                                            :: source
   end type

   ! Derived type for conserved quantity metadata
   type, extends(type_fabm_variable) :: type_fabm_conserved_quantity
      class (type_base_standard_variable), pointer :: standard_variable => null()
      integer                                      :: index             = -1
      integer                                      :: horizontal_index  = -1
      type (type_internal_variable),       pointer :: target_hz => null()
   end type

   ! --------------------------------------------------------------------------
   ! Derived type for a biogeochemical model as seen by the host
   ! --------------------------------------------------------------------------

   type type_fabm_model
      ! Variable metadata
      type (type_fabm_interior_state_variable),        allocatable, dimension(:) :: interior_state_variables
      type (type_fabm_horizontal_state_variable),      allocatable, dimension(:) :: surface_state_variables
      type (type_fabm_horizontal_state_variable),      allocatable, dimension(:) :: bottom_state_variables
      type (type_fabm_interior_diagnostic_variable),   allocatable, dimension(:) :: interior_diagnostic_variables
      type (type_fabm_horizontal_diagnostic_variable), allocatable, dimension(:) :: horizontal_diagnostic_variables
      type (type_fabm_conserved_quantity),             allocatable, dimension(:) :: conserved_quantities

      ! Names of variables taken as input by one or more biogeochemical models.
      ! These may be accessed by the host to enumerate potential forcing variables.
      character(len=attribute_length), allocatable, dimension(:) :: dependencies
      character(len=attribute_length), allocatable, dimension(:) :: dependencies_hz
      character(len=attribute_length), allocatable, dimension(:) :: dependencies_scalar

      ! Individual jobs
      type (type_job) :: get_interior_sources_job
      type (type_job) :: get_bottom_sources_job
      type (type_job) :: get_surface_sources_job
      type (type_job) :: get_vertical_movement_job
      type (type_job) :: get_interior_conserved_quantities_job
      type (type_job) :: get_horizontal_conserved_quantities_job
      type (type_job) :: finalize_outputs_job
      type (type_job) :: prepare_inputs_job
      type (type_job) :: check_interior_state_job
      type (type_job) :: check_bottom_state_job
      type (type_job) :: check_surface_state_job
      type (type_job) :: initialize_interior_state_job
      type (type_job) :: initialize_bottom_state_job
      type (type_job) :: initialize_surface_state_job

      ! Root container of biogeochemical modules
      type (type_base_model) :: root

      integer :: status = status_none
      logical :: log = .false.

      type (type_link_list) :: links_postcoupling

      type (type_global_variable_register) :: variable_register
      type (type_job_manager)              :: job_manager
      type (type_catalog)                  :: catalog
      type (type_store)                    :: store
      type (type_schedules)                :: schedules
      type (type_domain)                   :: domain

      ! Memory caches for exchanging information with individual biogeochemical modules
      type (type_interior_cache)   :: cache_int
      type (type_horizontal_cache) :: cache_hz
      type (type_vertical_cache)   :: cache_vert

      ! Cache fill values
      type (type_cache_fill_values) :: cache_fill_values
   contains
      procedure :: initialize
      procedure :: finalize
      procedure :: set_domain
#if _FABM_DIMENSION_COUNT_>0
      procedure :: set_domain_start
      procedure :: set_domain_stop
#endif
#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_BOTTOM_INDEX_==-1
      procedure :: set_bottom_index
#endif
#ifdef _HAS_MASK_
      procedure :: set_mask
#endif
      procedure :: start

      procedure :: initialize_interior_state
      procedure :: initialize_bottom_state
      procedure :: initialize_surface_state

      procedure :: check_interior_state
      procedure :: check_bottom_state
      procedure :: check_surface_state

      procedure :: prepare_inputs1
      procedure :: prepare_inputs2
      generic :: prepare_inputs => prepare_inputs1, prepare_inputs2
      procedure :: finalize_outputs

      procedure :: get_interior_sources_rhs
      procedure :: get_interior_sources_ppdd
      generic :: get_interior_sources => get_interior_sources_rhs, get_interior_sources_ppdd
      procedure :: get_bottom_sources_rhs
      procedure :: get_bottom_sources_ppdd
      generic :: get_bottom_sources => get_bottom_sources_rhs, get_bottom_sources_ppdd
      procedure :: get_surface_sources

      procedure :: get_vertical_movement
      procedure :: get_interior_conserved_quantities
      procedure :: get_horizontal_conserved_quantities

      procedure :: link_interior_data_by_variable
      procedure :: link_interior_data_by_id
      procedure :: link_interior_data_by_sn
      procedure :: link_interior_data_by_name
      generic :: link_interior_data => link_interior_data_by_variable, link_interior_data_by_id, link_interior_data_by_sn, link_interior_data_by_name

      procedure :: link_horizontal_data_by_variable
      procedure :: link_horizontal_data_by_id
      procedure :: link_horizontal_data_by_sn
      procedure :: link_horizontal_data_by_name
      generic :: link_horizontal_data => link_horizontal_data_by_variable, link_horizontal_data_by_id, link_horizontal_data_by_sn, link_horizontal_data_by_name

      procedure :: link_scalar_by_id
      procedure :: link_scalar_by_sn
      procedure :: link_scalar_by_name
      generic :: link_scalar => link_scalar_by_id, link_scalar_by_sn, link_scalar_by_name

      procedure :: link_interior_state_data
      procedure :: link_bottom_state_data
      procedure :: link_surface_state_data
      procedure :: link_all_interior_state_data
      procedure :: link_all_bottom_state_data
      procedure :: link_all_surface_state_data

      procedure :: require_interior_data
      procedure :: require_horizontal_data
      generic :: require_data => require_interior_data, require_horizontal_data

      procedure :: get_interior_data
      procedure :: get_horizontal_data
      procedure :: get_scalar_data
      generic :: get_data => get_interior_data, get_horizontal_data, get_scalar_data

      procedure :: get_interior_diagnostic_data
      procedure :: get_horizontal_diagnostic_data

      procedure :: get_interior_variable_id_by_name
      procedure :: get_interior_variable_id_sn
      generic :: get_interior_variable_id => get_interior_variable_id_by_name, get_interior_variable_id_sn

      procedure :: get_horizontal_variable_id_by_name
      procedure :: get_horizontal_variable_id_sn
      generic :: get_horizontal_variable_id => get_horizontal_variable_id_by_name, get_horizontal_variable_id_sn

      procedure :: get_scalar_variable_id_by_name
      procedure :: get_scalar_variable_id_sn
      generic :: get_scalar_variable_id => get_scalar_variable_id_by_name, get_scalar_variable_id_sn

      procedure, nopass :: is_variable_used
      procedure :: get_variable_name

      procedure :: interior_variable_needs_values
      procedure :: interior_variable_needs_values_sn
      procedure :: horizontal_variable_needs_values
      procedure :: horizontal_variable_needs_values_sn
      procedure :: scalar_variable_needs_values
      procedure :: scalar_variable_needs_values_sn
      generic :: variable_needs_values => interior_variable_needs_values, interior_variable_needs_values_sn, &
                                          horizontal_variable_needs_values, horizontal_variable_needs_values_sn, &
                                          scalar_variable_needs_values, scalar_variable_needs_values_sn

      procedure :: process_job
      generic :: process => process_job
#if _FABM_DIMENSION_COUNT_ > 1 || (_FABM_DIMENSION_COUNT_ == 1 && !defined(_FABM_DEPTH_DIMENSION_INDEX_))
      procedure :: process_job_everywhere
      generic :: process => process_job_everywhere
#endif

   end type type_fabm_model

contains

   ! --------------------------------------------------------------------------
   ! fabm_initialize_library: initialize FABM library
   ! --------------------------------------------------------------------------
   ! This will be called automatically when creating new models.
   ! For instance, from fabm_create_model.
   ! --------------------------------------------------------------------------
   subroutine fabm_initialize_library()
      use fabm_library, only: fabm_model_factory

      ! Do nothing if already initialized.
      if (associated(factory)) return

      ! If needed, create default object for communication (e.g., logging, error reporting) with host.
      if (.not. associated(driver)) allocate(driver)

      ! Create all standard variable objects.
      call initialize_standard_variables()

      ! Create the model factory.
      factory => fabm_model_factory
      call factory%initialize()
   end subroutine fabm_initialize_library

   ! --------------------------------------------------------------------------
   ! fabm_get_version: get FABM version string
   ! --------------------------------------------------------------------------
   subroutine fabm_get_version(string)
      use fabm_version

      character(len=*), intent(out) :: string

      type (type_version), pointer :: version

      call fabm_initialize_library()
      string = git_commit_id // ' (' // git_branch_name // ' branch)'
      version => first_module_version
      do while (associated(version))
         string = trim(string) // ', ' // trim(version%module_name) // ': ' // trim(version%version_string)
         version => version%next
      end do
   end subroutine fabm_get_version

   ! --------------------------------------------------------------------------
   ! fabm_create_model: create a model from a yaml-based configuration file
   ! --------------------------------------------------------------------------
   function fabm_create_model(path, initialize, parameters, unit) result(model)
      character(len=*),                optional, intent(in) :: path
      logical,                         optional, intent(in) :: initialize
      type (type_property_dictionary), optional, intent(in) :: parameters
      integer,                         optional, intent(in) :: unit
      class (type_fabm_model), pointer                      :: model

      logical :: initialize_

      ! Make sure the library is initialized.
      call fabm_initialize_library()

      allocate(model)
      call fabm_configure_model(model%root, model%schedules, model%log, path, parameters=parameters, unit=unit)

      ! Initialize model tree
      initialize_ = .true.
      if (present(initialize)) initialize_ = initialize
      if (initialize_) call model%initialize()
   end function fabm_create_model

   ! --------------------------------------------------------------------------
   ! initialize: initialize a model object
   ! --------------------------------------------------------------------------
   ! This freezes the tree of biogeochemical modules; afterwards no new modules
   ! can be added. This routine will be called automatically from
   ! fabm_create_model unless called with initialize=.false.
   ! --------------------------------------------------------------------------
   subroutine initialize(self)
      class (type_fabm_model), target, intent(inout) :: self

      class (type_property), pointer :: property => null()
      integer                        :: islash
      integer                        :: ivar

      if (self%status >= status_initialize_done) &
         call fatal_error('initialize', 'initialize has already been called on this model object.')

      ! Create zero fields.
      call self%root%add_interior_variable('zero', act_as_state_variable=.true., source=source_constant, missing_value=0.0_rki, output=output_none)
      call self%root%add_horizontal_variable('zero_hz', act_as_state_variable=.true., source=source_constant, missing_value=0.0_rki, output=output_none)

      ! Filter out expressions that FABM can handle itself.
      ! The remainder, if any, must be handled by the host model.
      call filter_expressions(self)

      ! This will resolve all FABM dependencies and generate final authoritative lists of variables of different types.
      call freeze_model_info(self%root)

      ! Raise error for unused coupling commands.
      property => self%root%couplings%first
      do while (associated(property))
         if (.not.self%root%couplings%retrieved%contains(trim(property%name))) then
            islash = index(property%name, '/', .true.)
            call fatal_error('initialize', 'model ' // property%name(1:islash-1) // ' does not contain variable "' // trim(property%name(islash+1:)) // '" mentioned in coupling section.')
         end if
         property => property%next
      end do

      ! Build final authoritative arrays with variable metadata.
      call classify_variables(self)

      ! Create catalog for storing pointers to data per variable.
      call create_catalog(self)

      ! Create built-in jobs, which can then be chained by the host/user by calling job%set_next.
      ! (the reason for chaining is to allow later jobs to use results of earlier ones, thus reducing the number of calls needed)
      call self%job_manager%create(self%prepare_inputs_job, 'prepare_inputs')
      call self%job_manager%create(self%get_interior_sources_job, 'get_interior_sources', source=source_do, previous=self%prepare_inputs_job)
      call self%job_manager%create(self%get_surface_sources_job, 'get_surface_sources', source=source_do_surface, previous=self%prepare_inputs_job)
      call self%job_manager%create(self%get_bottom_sources_job, 'get_bottom_sources', source=source_do_bottom, previous=self%prepare_inputs_job)
      call self%job_manager%create(self%get_interior_conserved_quantities_job, 'get_interior_conserved_quantities', source=source_do, previous=self%prepare_inputs_job)
      call self%job_manager%create(self%get_horizontal_conserved_quantities_job, 'get_horizontal_conserved_quantities', source=source_do_horizontal, previous=self%prepare_inputs_job)
      call self%job_manager%create(self%finalize_outputs_job, 'finalize_outputs', outsource_tasks=.true.)
      call self%get_interior_sources_job%connect(self%finalize_outputs_job)
      call self%get_surface_sources_job%connect(self%finalize_outputs_job)
      call self%get_bottom_sources_job%connect(self%finalize_outputs_job)
      !call self%get_interior_conserved_quantities_job%connect(self%finalize_outputs_job)
      !call self%get_horizontal_conserved_quantities_job%connect(self%finalize_outputs_job)
      call self%job_manager%create(self%get_vertical_movement_job, 'get_vertical_movement', source=source_get_vertical_movement, previous=self%finalize_outputs_job)
      call self%job_manager%create(self%initialize_interior_state_job, 'initialize_interior_state', source=source_initialize_state, previous=self%finalize_outputs_job)
      call self%job_manager%create(self%initialize_bottom_state_job, 'initialize_bottom_state', source=source_initialize_bottom_state, previous=self%finalize_outputs_job)
      call self%job_manager%create(self%initialize_surface_state_job, 'initialize_surface_state', source=source_initialize_surface_state, previous=self%finalize_outputs_job)
      call self%job_manager%create(self%check_interior_state_job, 'check_interior_state', source=source_check_state, previous=self%finalize_outputs_job)
      call self%job_manager%create(self%check_bottom_state_job, 'check_bottom_state', source=source_check_bottom_state, previous=self%finalize_outputs_job)
      call self%job_manager%create(self%check_surface_state_job, 'check_surface_state', source=source_check_surface_state, previous=self%finalize_outputs_job)

      call require_flux_computation(self%get_bottom_sources_job, self%links_postcoupling, domain_bottom)
      call require_flux_computation(self%get_surface_sources_job, self%links_postcoupling, domain_surface)
      call require_flux_computation(self%get_interior_sources_job, self%links_postcoupling, domain_interior)
      call require_flux_computation(self%get_vertical_movement_job, self%links_postcoupling, domain_interior + 999)

      call require_call_all_with_state(self%initialize_interior_state_job, self%root%links, domain_interior, source_initialize_state)
      call require_call_all_with_state(self%initialize_bottom_state_job, self%root%links, domain_bottom, source_initialize_bottom_state)
      call require_call_all_with_state(self%initialize_surface_state_job, self%root%links, domain_surface, source_initialize_surface_state)
      call require_call_all_with_state(self%check_interior_state_job, self%root%links, domain_interior, source_check_state)
      call require_call_all_with_state(self%check_bottom_state_job, self%root%links, domain_bottom, source_check_bottom_state)
      call require_call_all_with_state(self%check_bottom_state_job, self%root%links, domain_interior, source_check_bottom_state)
      call require_call_all_with_state(self%check_surface_state_job, self%root%links, domain_surface, source_check_surface_state)
      call require_call_all_with_state(self%check_surface_state_job, self%root%links, domain_interior, source_check_surface_state)

      do ivar = 1, size(self%interior_state_variables)
         call self%check_interior_state_job%read_cache_loads%add(self%interior_state_variables(ivar)%target)
      end do
      do ivar = 1, size(self%bottom_state_variables)
         call self%check_bottom_state_job%read_cache_loads%add(self%bottom_state_variables(ivar)%target)
      end do
      do ivar = 1, size(self%surface_state_variables)
         call self%check_surface_state_job%read_cache_loads%add(self%surface_state_variables(ivar)%target)
      end do

      do ivar = 1, size(self%conserved_quantities)
         call self%get_interior_conserved_quantities_job%request_variable(self%conserved_quantities(ivar)%target)
         call self%get_horizontal_conserved_quantities_job%request_variable(self%conserved_quantities(ivar)%target_hz)
         call self%conserved_quantities(ivar)%target%write_indices%append(self%conserved_quantities(ivar)%index)
         call self%conserved_quantities(ivar)%target_hz%write_indices%append(self%conserved_quantities(ivar)%horizontal_index)
      end do

      self%status = status_initialize_done
   end subroutine initialize

   ! --------------------------------------------------------------------------
   ! finalize: deallocate model object
   ! --------------------------------------------------------------------------
   subroutine finalize(self)
      class (type_fabm_model),target,intent(inout) :: self
      self%status = status_none

      ! TODO: this should deallocate the memory of all biogeochemical models
   end subroutine finalize

   ! --------------------------------------------------------------------------
   ! set_domain: set extents of spatial domain and optionally time step length
   ! --------------------------------------------------------------------------
   subroutine set_domain(self _POSTARG_LOCATION_, seconds_per_time_unit)
      class (type_fabm_model), target, intent(inout) :: self
      _DECLARE_ARGUMENTS_LOCATION_
      real(rke), optional, intent(in)    :: seconds_per_time_unit

      class (type_expression), pointer :: expression

      if (self%status < status_initialize_done) call fatal_error('set_domain', 'initialize has not yet been called on this model object.')
      if (self%status >= status_set_domain_done) call fatal_error('set_domain', 'set_domain has already been called on this model object.')
      self%status = status_set_domain_done

#if _FABM_DIMENSION_COUNT_>0
      self%domain%shape(:) = (/_LOCATION_/)
      self%domain%start(:) = 1
      self%domain%stop(:) = self%domain%shape
#endif
#if _HORIZONTAL_DIMENSION_COUNT_>0
      self%domain%horizontal_shape(:) = (/_HORIZONTAL_LOCATION_/)
#endif

      if (present(seconds_per_time_unit)) then
         expression => self%root%first_expression
         do while (associated(expression))
            select type (expression)
            class is (type_interior_temporal_mean)
               expression%in = expression%link%target%catalog_index
               expression%period = expression%period / seconds_per_time_unit
               allocate(expression%history(_PREARG_LOCATION_ expression%n + 3))
               expression%history = 0.0_rke
               call self%link_interior_data(expression%output_name, &
                                            expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 3))
            class is (type_horizontal_temporal_mean)
               expression%in = expression%link%target%catalog_index
               expression%period = expression%period / seconds_per_time_unit
               allocate(expression%history(_PREARG_HORIZONTAL_LOCATION_ expression%n + 3))
               expression%history = 0.0_rke
               call self%link_horizontal_data(expression%output_name, &
                                              expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 3))
            end select
            expression => expression%next
         end do
      end if
   end subroutine set_domain

#if _FABM_DIMENSION_COUNT_>0
   ! --------------------------------------------------------------------------
   ! set_domain_start: set start index of all spatial dimensions
   ! --------------------------------------------------------------------------
   ! This is optional; by default the start index for all dimensions is 1.
   ! --------------------------------------------------------------------------
   subroutine set_domain_start(self _POSTARG_LOCATION_)
      class (type_fabm_model), target, intent(inout) :: self
      _DECLARE_ARGUMENTS_LOCATION_
      if (self%status < status_set_domain_done) &
         call fatal_error('set_domain_start', 'set_domain has not yet been called on this model object.')
      self%domain%start(:) = (/_LOCATION_/)
   end subroutine set_domain_start

   ! --------------------------------------------------------------------------
   ! set_domain_stop: set stop index of all spatial dimensions (default=domain size)
   ! --------------------------------------------------------------------------
   ! This is optional; by default the stop index for all dimensions matches
   ! the domain size provided to set_domain.
   ! --------------------------------------------------------------------------
   subroutine set_domain_stop(self _POSTARG_LOCATION_)
      class (type_fabm_model), target, intent(inout) :: self
      _DECLARE_ARGUMENTS_LOCATION_
      if (self%status < status_set_domain_done) &
         call fatal_error('set_domain_stop', 'set_domain has not yet been called on this model object.')
      self%domain%stop(:) = (/_LOCATION_/)
   end subroutine set_domain_stop
#endif

#ifdef _HAS_MASK_
   ! --------------------------------------------------------------------------
   ! set_mask: provide spatial mask
   ! --------------------------------------------------------------------------
   ! As FABM will keep a pointer to the mask, it needs to remain valid for
   ! the lifetime of the model object.
   ! --------------------------------------------------------------------------
#  ifdef _FABM_HORIZONTAL_MASK_
   subroutine set_mask(self, mask_hz)
#  else
   subroutine set_mask(self, mask, mask_hz)
#  endif
      class (type_fabm_model), target, intent(inout)                      :: self
#  ifndef _FABM_HORIZONTAL_MASK_
      _FABM_MASK_TYPE_, target, intent(in) _ATTRIBUTES_GLOBAL_            :: mask
#  endif
      _FABM_MASK_TYPE_, target, intent(in) _ATTRIBUTES_GLOBAL_HORIZONTAL_ :: mask_hz

      integer :: i

      if (self%status < status_set_domain_done) &
         call fatal_error('set_mask', 'set_domain has not yet been called on this model object.')

#  ifndef _FABM_HORIZONTAL_MASK_
#    if !defined(NDEBUG)&&_FABM_DIMENSION_COUNT_>0
      do i = 1, size(self%domain%shape)
         if (size(mask, i) /= self%domain%shape(i)) &
            call fatal_error('set_mask', 'shape of provided mask does not match domain extents provided to set_domain.')
      end do
#    endif
      self%domain%mask => mask
#  endif

#  if !defined(NDEBUG)&&_HORIZONTAL_DIMENSION_COUNT_>0
      do i = 1, size(self%domain%horizontal_shape)
         if (size(mask_hz, i) /= self%domain%horizontal_shape(i)) &
            call fatal_error('set_mask', 'shape of provided horizontal mask does not match domain extents provided to set_domain.')
      end do
#  endif
      self%domain%mask_hz => mask_hz

   end subroutine set_mask
#endif

#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_BOTTOM_INDEX_==-1
   ! --------------------------------------------------------------------------
   ! set_bottom_index: provide bottom indices for every horizontal point
   ! --------------------------------------------------------------------------
   ! As FABM will keep a pointer to the array with indices, it needs to remain
   ! valid for the lifetime of the model object.
   ! --------------------------------------------------------------------------
   subroutine set_bottom_index(self, indices)
      class (type_fabm_model), intent(inout)                             :: self
      integer, target,         intent(in) _ATTRIBUTES_GLOBAL_HORIZONTAL_ :: indices

      integer :: i

      if (self%status < status_set_domain_done) &
         call fatal_error('set_bottom_index', 'set_domain has not yet been called on this model object.')
#  if !defined(NDEBUG)&&_HORIZONTAL_DIMENSION_COUNT_>0
      do i = 1, size(self%domain%horizontal_shape)
         if (size(indices, i) /= self%domain%horizontal_shape(i)) &
            call fatal_error('set_bottom_index', 'shape of provided index array does not match domain extents provided to set_domain.')
      end do
#  endif

      self%domain%bottom_indices => indices
   end subroutine set_bottom_index
#endif

   ! --------------------------------------------------------------------------
   ! start: prepare for simulation start
   ! --------------------------------------------------------------------------
   ! This tells FABM that the user/host have finished providing (or overriding)
   ! data (link_data procedures) and have finished flagging diagnostics for
   ! output (by setting the "save" flag that is part of the variable metadata)
   ! --------------------------------------------------------------------------
   subroutine start(self)
      class (type_fabm_model), intent(inout), target :: self

      integer                           :: ivar
      logical                           :: ready
      type (type_variable_set)          :: unfulfilled_dependencies
      type (type_variable_node),pointer :: variable_node
      type (type_link), pointer         :: link
      character(len=*), parameter       :: log_prefix = 'fabm_'
      integer                           :: log_unit, ios

      if (self%status < status_set_domain_done) then
         call fatal_error('start', 'set_domain has not yet been called on this model object.')
         return
      elseif (self%status >= status_start_done) then
         ! start has been called on this model before and it must have succeeded to have this status.
         ! Reset store (e.g., all diagnostics) by setting each variable to its fill value and return.
         ! (this allows masked cells to be properly initialized if the mask changes between calls to start)
         call reset_store(self)
         return
      end if

      ready = .true.

#ifdef _HAS_MASK_
#  ifndef _FABM_HORIZONTAL_MASK_
      if (.not. associated(self%domain%mask)) then
         call log_message('spatial mask has not been set. Make sure to call set_mask.')
         ready = .false.
      end if
#  endif
      if (.not. associated(self%domain%mask_hz)) then
         call log_message('horizontal spatial mask has not been set. Make sure to call set_mask.')
         ready = .false.
      end if
#endif

#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_BOTTOM_INDEX_==-1
      if (.not. associated(self%domain%bottom_indices)) then
         call log_message('bottom indices have not been set. Make sure to call set_bottom_index.')
         ready = .false.
      end if
#endif

      ! Flag variables that have had data asssigned (by user, host or FABM).
      ! This is done only now because the user/host had till this moment to provide (or override) model fields.
      call flag_variables_with_data(self%variable_register%catalog%interior, self%catalog%interior_sources)
      call flag_variables_with_data(self%variable_register%catalog%horizontal, self%catalog%horizontal_sources)
      call flag_variables_with_data(self%variable_register%catalog%scalar, self%catalog%scalar_sources)

      ! Create job that ensures all diagnostics required by the user are computed.
      ! This is done only now because the user/host had till this moment to change the "save" flag of each diagnostic.
      do ivar = 1, size(self%interior_diagnostic_variables)
         if (self%interior_diagnostic_variables(ivar)%save) then
            select case (self%interior_diagnostic_variables(ivar)%target%source)
            case (source_check_state)
               call self%check_interior_state_job%request_variable(self%interior_diagnostic_variables(ivar)%target, store=.true.)
            case (source_get_vertical_movement)
               call self%get_vertical_movement_job%request_variable(self%interior_diagnostic_variables(ivar)%target, store=.true.)
            case default
               call self%finalize_outputs_job%request_variable(self%interior_diagnostic_variables(ivar)%target, store=.true.)
            end select
         end if
      end do
      do ivar = 1, size(self%horizontal_diagnostic_variables)
         if (self%horizontal_diagnostic_variables(ivar)%save) &
            call self%finalize_outputs_job%request_variable(self%horizontal_diagnostic_variables(ivar)%target, store=.true.)
      end do

      log_unit = -1
      if (self%log) log_unit = get_free_unit()

      ! Merge write indices when operations can be done in place
      ! This must be done after all variables are requested from the different jobs, so we know which variables
      ! will be retrieved (such variables cannot be merged)
      if (self%log) then
         open(unit=log_unit, file=log_prefix // 'merges.log', action='write', status='replace', iostat=ios)
         if (ios /= 0) call fatal_error('start', 'Unable to open ' // log_prefix // 'merges.log')
         call merge_indices(self%root, log_unit)
         close(log_unit)
      else
         call merge_indices(self%root)
      end if

      ! Initialize all jobs. This also creates registers for the read and write caches, as well as the persistent store.
      if (self%log) then
         open(unit=log_unit, file=log_prefix // 'task_order.log', action='write', status='replace', iostat=ios)
         if (ios /= 0) call fatal_error('start', 'Unable to open ' // log_prefix // 'task_order.log')
      end if
      call self%job_manager%initialize(self%variable_register, self%schedules, unfulfilled_dependencies, log_unit)
      if (self%log) then
         close(log_unit)
         open(unit=log_unit, file=log_prefix // 'graph.gv', action='write', status='replace', iostat=ios)
         if (ios /= 0) call fatal_error('start', 'Unable to open ' // log_prefix // 'graph.gv')
         call self%job_manager%write_graph(log_unit)
         close(log_unit)
      end if

      ! Create persistent store. This provides memory for all variables to be stored there.
      call create_store(self)

      ! Collect fill values and missing values for cache entries.
      self%cache_fill_values = get_cache_fill_values(self%variable_register)

      ! Create global caches for exchanging information with BGC models.
      ! This can only be done after get_cache_fill_values completes, because that determines what values to prefill the cache with.
      call cache_create(self%domain, self%cache_fill_values, self%cache_int)
      call cache_create(self%domain, self%cache_fill_values, self%cache_hz)
      call cache_create(self%domain, self%cache_fill_values, self%cache_vert)

      ! For diagnostics that are not needed, set their write index to 0 (rubbish bin)
      if (self%log) then
         open(unit=log_unit, file=log_prefix // 'discards.log', action='write', status='replace', iostat=ios)
         if (ios /= 0) call fatal_error('start', 'Unable to open ' // log_prefix // 'discards.log')
         write (log_unit,'(a)') 'Writes for the following variables are discarded:'
      end if
      link => self%links_postcoupling%first
      do while (associated(link))
         if (.not. link%target%write_indices%is_empty() .and. link%target%write_indices%value == -1) then
            call link%target%write_indices%set_value(0)
            if (self%log) write (log_unit,'("- ",a)') trim(link%target%name)
         end if
         link => link%next
      end do
      if (self%log) close(log_unit)

      if (self%log) then
         open(unit=log_unit, file=log_prefix // 'register.log', action='write', status='replace', iostat=ios)
         if (ios /= 0) call fatal_error('start', 'Unable to open ' // log_prefix // 'register.log')
         call self%variable_register%print(log_unit)
         close(log_unit)

         open(unit=log_unit, file=log_prefix // 'jobs.log', action='write', status='replace', iostat=ios)
         if (ios /= 0) call fatal_error('start', 'Unable to open ' // log_prefix // 'jobs.log')
         call self%job_manager%print(log_unit)
         close(log_unit)
      end if

      ! Report all unfulfilled dependencies.
      variable_node => unfulfilled_dependencies%first
      do while (associated(variable_node))
         call report_unfulfilled_dependency(variable_node%target)
         variable_node => variable_node%next
      end do

      if (associated(unfulfilled_dependencies%first) .or. .not. ready) &
         call fatal_error('start', 'FABM is lacking required data.')

      self%status = status_start_done

   contains

      subroutine flag_variables_with_data(variable_list, data_sources)
         type (type_variable_list), intent(inout) :: variable_list
         integer,                   intent(in)    :: data_sources(:)

         integer :: i
         type (type_variable_node), pointer :: variable_node

         variable_node => variable_list%first
         do i = 1, variable_list%count
            variable_node%target%has_data = data_sources(i) /= data_source_none
            if (data_sources(i) > data_source_fabm) then
               variable_node%target%source = source_external
               variable_node%target%write_operator = operator_assign
            elseif (data_sources(i) /= data_source_none .and. variable_node%target%source == source_unknown) then
               variable_node%target%source = source_external
            end if
            variable_node => variable_node%next
         end do
      end subroutine

      subroutine report_unfulfilled_dependency(variable)
         type (type_internal_variable), target :: variable

         type type_model_reference
            class (type_base_model),     pointer :: p    => null()
            type (type_model_reference), pointer :: next => null()
         end type

         type (type_model_reference), pointer :: first, current, next
         type (type_link),            pointer :: link
         character(len=attribute_length)      :: path

         call log_message('UNFULFILLED DEPENDENCY: ' // trim(variable%name))
         select case (variable%domain)
         case (domain_interior)
            call log_message('  This is an interior field.')
         case (domain_horizontal, domain_surface, domain_bottom)
            call log_message('  This is a horizontal field.')
         case (domain_scalar)
            call log_message('  This is a scalar field (single value valid across the entire domain).')
         end select
         if (variable%units /= '') call log_message('  It has units ' // trim(variable%units))
         call log_message('  It is needed by the following model instances:')

         first => null()
         link => self%root%links%first
         do while (associated(link))
            if (      associated(link%target, variable)    &                     ! This link points to the target variable,
                .and. associated(link%original%read_index) &                     ! the model that owns the link requests read access for it,
                .and. link%original%presence /= presence_external_optional) then ! and this access is required, not optional
               current => first
               do while (associated(current))
                  if (associated(current%p, link%original%owner)) exit
                  current => current%next
               end do
               if (.not. associated(current)) then
                  ! This model has not been reported before. Do so now and remember that we have done so.
                  allocate(current)
                  current%p => link%original%owner
                  current%next => first
                  first => current
                  path = current%p%get_path()
                  call log_message('    ' // trim(path(2:)))
               end if
            end if
            link => link%next
         end do

         ! Clean up model list
         current => first
         do while (associated(current))
            next => current%next
            deallocate(current)
            current => next
         end do
      end subroutine
   end subroutine start

   ! --------------------------------------------------------------------------
   ! get_interior_variable_id_by_name: get interior variable identifier for
   ! given variable name
   ! --------------------------------------------------------------------------
   function get_interior_variable_id_by_name(self, name) result(id)
      class (type_fabm_model), intent(in)   :: self
      character(len=*),        intent(in)   :: name
      type (type_fabm_interior_variable_id) :: id

      id%variable => get_variable_by_name(self, name, domain_interior)
   end function get_interior_variable_id_by_name

   ! --------------------------------------------------------------------------
   ! get_interior_variable_id_sn: get interior variable identifier for given
   ! standard variable
   ! --------------------------------------------------------------------------
   function get_interior_variable_id_sn(self, standard_variable) result(id)
      class (type_fabm_model),                intent(in) :: self
      type (type_interior_standard_variable), intent(in) :: standard_variable
      type (type_fabm_interior_variable_id)              :: id

      id%variable => get_variable_by_standard_variable(self, standard_variable%resolve())
   end function get_interior_variable_id_sn

   ! --------------------------------------------------------------------------
   ! get_horizontal_variable_id_by_name: get horizontal variable identifier for
   ! given variable name
   ! --------------------------------------------------------------------------
   function get_horizontal_variable_id_by_name(self, name) result(id)
      class (type_fabm_model), intent(in)     :: self
      character(len=*),        intent(in)     :: name
      type (type_fabm_horizontal_variable_id) :: id

      id%variable => get_variable_by_name(self, name, domain_horizontal)
   end function get_horizontal_variable_id_by_name

   ! --------------------------------------------------------------------------
   ! get_horizontal_variable_id_sn: get horizontal variable identifier for
   ! given standard variable
   ! --------------------------------------------------------------------------
   function get_horizontal_variable_id_sn(self, standard_variable) result(id)
      class (type_fabm_model),                   intent(in) :: self
      class (type_horizontal_standard_variable), intent(in) :: standard_variable
      type (type_fabm_horizontal_variable_id)               :: id

      id%variable => get_variable_by_standard_variable(self, standard_variable%resolve())
   end function get_horizontal_variable_id_sn

   ! --------------------------------------------------------------------------
   ! get_scalar_variable_id_by_name: get scalar variable identifier for given
   ! variable name
   ! --------------------------------------------------------------------------
   function get_scalar_variable_id_by_name(self, name) result(id)
      class (type_fabm_model), intent(in) :: self
      character(len=*),        intent(in) :: name
      type (type_fabm_scalar_variable_id) :: id

      id%variable => get_variable_by_name(self, name, domain_scalar)
   end function get_scalar_variable_id_by_name

   ! --------------------------------------------------------------------------
   ! get_scalar_variable_id_sn: get scalar variable identifier for given
   ! standard variable
   ! --------------------------------------------------------------------------
   function get_scalar_variable_id_sn(self, standard_variable) result(id)
      class (type_fabm_model),              intent(in) :: self
      type (type_global_standard_variable), intent(in) :: standard_variable
      type (type_fabm_scalar_variable_id)              :: id

      id%variable => get_variable_by_standard_variable(self, standard_variable%resolve())
   end function get_scalar_variable_id_sn

   ! --------------------------------------------------------------------------
   ! get_variable_name: get output name associated with given variable id.
   ! The name consists of alphanumeric characters and underscores only.
   ! --------------------------------------------------------------------------
   function get_variable_name(self, id) result(name)
      class (type_fabm_model),       intent(in) :: self
      class (type_fabm_variable_id), intent(in) :: id
      character(len=attribute_length)           :: name

      name = ''
      if (associated(id%variable)) name = get_safe_name(id%variable%name)
   end function get_variable_name

   ! --------------------------------------------------------------------------
   ! is_variable_used: returns whether this variable is an input for any
   ! biogeochemical module
   ! --------------------------------------------------------------------------
   function is_variable_used(id) result(used)
      class (type_fabm_variable_id), intent(in) :: id
      logical                                   :: used

      used = associated(id%variable)
      if (used) used = .not. id%variable%read_indices%is_empty()
   end function is_variable_used

   ! --------------------------------------------------------------------------
   ! interior_variable_needs_values: returns whether values still need to
   ! provided for this interior variable, identified by id.
   ! Unless these values are provided, a call to "start" will fail.
   ! --------------------------------------------------------------------------
   function interior_variable_needs_values(self, id) result(required)
      class (type_fabm_model),               intent(in) :: self
      type (type_fabm_interior_variable_id), intent(in) :: id
      logical                                           :: required

      required = associated(id%variable)
      if (required) required = .not. id%variable%read_indices%is_empty()
      if (required) required = .not. associated(self%catalog%interior(id%variable%catalog_index)%p)
   end function interior_variable_needs_values

   ! --------------------------------------------------------------------------
   ! interior_variable_needs_values_sn: returns whether values still need to
   ! provided for this interior variable, identified by standard variable.
   ! Unless these values are provided, a call to "start" will fail.
   ! --------------------------------------------------------------------------
   function interior_variable_needs_values_sn(self, standard_variable) result(required)
      class (type_fabm_model),                intent(in) :: self
      type (type_interior_standard_variable), intent(in) :: standard_variable
      logical                                            :: required

      required = interior_variable_needs_values(self, get_interior_variable_id_sn(self, standard_variable))
   end function interior_variable_needs_values_sn

   ! --------------------------------------------------------------------------
   ! horizontal_variable_needs_values: returns whether values still need to
   ! provided for this horizontal variable, identified by id.
   ! Unless these values are provided, a call to "start" will fail.
   ! --------------------------------------------------------------------------
   function horizontal_variable_needs_values(self, id) result(required)
      class (type_fabm_model),                 intent(in) :: self
      type (type_fabm_horizontal_variable_id), intent(in) :: id
      logical                                             :: required

      required = associated(id%variable)
      if (required) required = .not. id%variable%read_indices%is_empty()
      if (required) required = .not. associated(self%catalog%horizontal(id%variable%catalog_index)%p)
   end function horizontal_variable_needs_values

   ! --------------------------------------------------------------------------
   ! horizontal_variable_needs_values_sn: returns whether values still need to
   ! provided for this horizontal variable, identified by standard variable.
   ! Unless these values are provided, a call to "start" will fail.
   ! --------------------------------------------------------------------------
   function horizontal_variable_needs_values_sn(self, standard_variable) result(required)
      class (type_fabm_model),                   intent(in) :: self
      class (type_horizontal_standard_variable), intent(in) :: standard_variable
      logical                                               :: required

      required = horizontal_variable_needs_values(self, get_horizontal_variable_id_sn(self, standard_variable))
   end function horizontal_variable_needs_values_sn

   ! --------------------------------------------------------------------------
   ! scalar_variable_needs_values: returns whether a value still need to
   ! provided for this scalar variable, identified by id.
   ! Unless this value is provided, a call to "start" will fail.
   ! --------------------------------------------------------------------------
   function scalar_variable_needs_values(self, id) result(required)
      class (type_fabm_model),             intent(in) :: self
      type (type_fabm_scalar_variable_id), intent(in) :: id
      logical                                         :: required

      required = associated(id%variable)
      if (required) required = .not. id%variable%read_indices%is_empty()
      if (required) required = .not. associated(self%catalog%scalar(id%variable%catalog_index)%p)
   end function scalar_variable_needs_values

   ! --------------------------------------------------------------------------
   ! scalar_variable_needs_values: returns whether a value still need to
   ! provided for this scalar variable, identified by standard variable.
   ! Unless this value is provided, a call to "start" will fail.
   ! --------------------------------------------------------------------------
   function scalar_variable_needs_values_sn(self, standard_variable) result(required)
      class (type_fabm_model),              intent(in) :: self
      type (type_global_standard_variable), intent(in) :: standard_variable
      logical                                          :: required

      required = scalar_variable_needs_values(self, get_scalar_variable_id_sn(self, standard_variable))
   end function scalar_variable_needs_values_sn

   subroutine require_interior_data(self, standard_variable)
      class (type_fabm_model),                intent(inout) :: self
      type (type_interior_standard_variable), intent(in)    :: standard_variable

      type (type_fabm_interior_variable_id) :: id

      if (self%status < status_initialize_done) &
         call fatal_error('require_interior_data', 'This procedure can only be called after model initialization.')
      if (self%status >= status_start_done) &
         call fatal_error('require_interior_data', 'This procedure cannot be called after start is called.')

      id = self%get_interior_variable_id(standard_variable)
      if (.not. associated(id%variable)) &
         call fatal_error('require_interior_data', 'Model does not contain requested variable ' // trim(standard_variable%name))
      call self%finalize_outputs_job%request_variable(id%variable, store=.true.)
   end subroutine require_interior_data

   subroutine require_horizontal_data(self, standard_variable)
      class (type_fabm_model),                   intent(inout) :: self
      class (type_horizontal_standard_variable), intent(in)    :: standard_variable

      type (type_fabm_horizontal_variable_id) :: id

      if (self%status < status_initialize_done) &
         call fatal_error('require_horizontal_data', 'This procedure can only be called after model initialization.')
      if (self%status >= status_start_done) &
         call fatal_error('require_horizontal_data', 'This procedure cannot be called after check_ready is called.')

      id = self%get_horizontal_variable_id(standard_variable)
      if (.not. associated(id%variable)) &
         call fatal_error('require_horizontal_data', 'Model does not contain requested variable ' // trim(standard_variable%name))
      call self%finalize_outputs_job%request_variable(id%variable, store=.true.)
   end subroutine require_horizontal_data

   subroutine link_interior_data_by_variable(self, variable, dat, source)
      class (type_fabm_model),               intent(inout) :: self
      type (type_internal_variable),         intent(in)    :: variable
      real(rke) _ATTRIBUTES_GLOBAL_, target, intent(in)    :: dat
      integer, optional,                     intent(in)    :: source

      integer :: i
      integer :: source_

#if !defined(NDEBUG)&&_FABM_DIMENSION_COUNT_>0
      do i = 1, size(self%domain%shape)
         if (size(dat, i) /= self%domain%shape(i)) then
            call fatal_error('link_interior_data_by_variable', trim(variable%name) // &
               ': extents of provided array do not match domain extents.')
         end if
      end do
#endif

      i = variable%catalog_index
      if (i /= -1) then
         source_ = data_source_default
         if (present(source)) source_ = source
         if (source_ >= self%catalog%interior_sources(i)) then
            self%catalog%interior(i)%p => dat
            self%catalog%interior_sources(i) = source_
         end if
      end if
   end subroutine link_interior_data_by_variable

   subroutine link_interior_data_by_id(self, id, dat, source)
      class (type_fabm_model),               intent(inout) :: self
      type(type_fabm_interior_variable_id),  intent(in)    :: id
      real(rke) _ATTRIBUTES_GLOBAL_, target, intent(in)    :: dat
      integer,optional,                      intent(in)    :: source

      if (associated(id%variable)) call link_interior_data_by_variable(self, id%variable, dat, source)
   end subroutine link_interior_data_by_id

   subroutine link_interior_data_by_sn(self, standard_variable, dat)
      class (type_fabm_model),                intent(inout) :: self
      type (type_interior_standard_variable), intent(in)    :: standard_variable
      real(rke) _ATTRIBUTES_GLOBAL_, target,  intent(in)    :: dat

      call link_interior_data_by_id(self, get_interior_variable_id_sn(self, standard_variable), dat)
   end subroutine link_interior_data_by_sn

   subroutine link_interior_data_by_name(self, name, dat)
      class (type_fabm_model),       target, intent(inout) :: self
      character(len=*),                      intent(in)    :: name
      real(rke) _ATTRIBUTES_GLOBAL_, target, intent(in)    :: dat

      call link_interior_data_by_id(self, get_interior_variable_id_by_name(self, name), dat)
   end subroutine link_interior_data_by_name

   subroutine link_horizontal_data_by_variable(self, variable, dat, source)
      class (type_fabm_model),                          intent(inout) :: self
      type (type_internal_variable),                    intent(in)    :: variable
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, target, intent(in)    :: dat
      integer, optional,                                intent(in)    :: source

      integer :: i
      integer :: source_

#if !defined(NDEBUG)&&_HORIZONTAL_DIMENSION_COUNT_>0
      do i = 1, size(self%domain%horizontal_shape)
         if (size(dat, i) /= self%domain%horizontal_shape(i)) then
            call fatal_error('link_horizontal_data_by_variable', trim(variable%name) // &
               ': extents of provided array do not match domain extents.')
         end if
      end do
#endif

      i = variable%catalog_index
      if (i /= -1) then
         source_ = data_source_default
         if (present(source)) source_ = source
         if (source_ >= self%catalog%horizontal_sources(i)) then
            self%catalog%horizontal(i)%p => dat
            self%catalog%horizontal_sources(i) = source_
         end if
      end if
   end subroutine link_horizontal_data_by_variable

   subroutine link_horizontal_data_by_id(self, id, dat, source)
      class (type_fabm_model),                          intent(inout) :: self
      type (type_fabm_horizontal_variable_id),          intent(in)    :: id
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, target, intent(in)    :: dat
      integer, optional,                                intent(in)    :: source

      if (associated(id%variable)) call link_horizontal_data_by_variable(self, id%variable, dat, source)
   end subroutine link_horizontal_data_by_id

   subroutine link_horizontal_data_by_sn(self, standard_variable, dat)
      class (type_fabm_model),                          intent(inout) :: self
      class (type_horizontal_standard_variable),        intent(in)    :: standard_variable
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, target, intent(in)    :: dat

      call link_horizontal_data_by_id(self, get_horizontal_variable_id_sn(self, standard_variable), dat)
   end subroutine link_horizontal_data_by_sn

   subroutine link_horizontal_data_by_name(self, name, dat)
      class (type_fabm_model),                          intent(inout) :: self
      character(len=*),                                 intent(in)    :: name
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, target, intent(in)    :: dat

      call link_horizontal_data_by_id(self, get_horizontal_variable_id_by_name(self, name), dat)
   end subroutine link_horizontal_data_by_name

   subroutine link_scalar_by_id(self, id, dat, source)
      class (type_fabm_model),             intent(inout) :: self
      type (type_fabm_scalar_variable_id), intent(in)    :: id
      real(rke), target,                   intent(in)    :: dat
      integer, optional,                   intent(in)    :: source

      integer :: i
      integer :: source_

      if (.not. associated(id%variable)) return
      i = id%variable%catalog_index
      if (i /= -1) then
         source_ = data_source_default
         if (present(source)) source_ = source
         if (source_ >= self%catalog%scalar_sources(i)) then
            self%catalog%scalar(i)%p => dat
            self%catalog%scalar_sources(i) = source_
         end if
      end if
   end subroutine link_scalar_by_id

   subroutine link_scalar_by_sn(self, standard_variable, dat)
      class (type_fabm_model),             intent(inout) :: self
      type(type_global_standard_variable), intent(in)    :: standard_variable
      real(rke), target,                   intent(in)    :: dat

      call link_scalar_by_id(self, get_scalar_variable_id_sn(self, standard_variable), dat)
   end subroutine link_scalar_by_sn

   subroutine link_scalar_by_name(self, name, dat)
      class (type_fabm_model), intent(inout) :: self
      character(len=*),        intent(in)    :: name
      real(rke), target,       intent(in)    :: dat

      call link_scalar_by_id(self, get_scalar_variable_id_by_name(self, name), dat)
   end subroutine link_scalar_by_name

   subroutine link_interior_state_data(self, index, dat)
      class (type_fabm_model),               intent(inout) :: self
      integer,                               intent(in)    :: index
      real(rke) _ATTRIBUTES_GLOBAL_, target, intent(in)    :: dat

      call link_interior_data_by_variable(self, self%interior_state_variables(index)%target, dat, source=data_source_fabm)
   end subroutine link_interior_state_data

   subroutine link_bottom_state_data(self, index, dat)
      class (type_fabm_model),                          intent(inout) :: self
      integer,                                          intent(in)    :: index
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, target, intent(in)    :: dat

      call link_horizontal_data_by_variable(self, self%bottom_state_variables(index)%target, dat, source=data_source_fabm)
   end subroutine link_bottom_state_data

   subroutine link_surface_state_data(self, index, dat)
      class (type_fabm_model),                          intent(inout) :: self
      integer,                                          intent(in)    :: index
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, target, intent(in)    :: dat

      call link_horizontal_data_by_variable(self, self%surface_state_variables(index)%target, dat, source=data_source_fabm)
   end subroutine link_surface_state_data

   subroutine link_all_interior_state_data(self, dat)
      class (type_fabm_model),                     intent(inout) :: self
      real(rke) _DIMENSION_GLOBAL_PLUS_1_, target, intent(in)    :: dat

      integer :: i

#ifndef NDEBUG
      if (size(dat, _FABM_DIMENSION_COUNT_ + 1) /= size(self%interior_state_variables)) &
         call fatal_error('link_all_interior_state_data', 'size of last dimension of provided array must match number of interior state variables.')
#endif
      do i = 1, size(self%interior_state_variables)
         call link_interior_state_data(self, i, dat(_PREARG_LOCATION_DIMENSIONS_ i))
      end do
   end subroutine link_all_interior_state_data

   subroutine link_all_bottom_state_data(self, dat)
      class (type_fabm_model),                                intent(inout) :: self
      real(rke) _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_, target, intent(in)    :: dat

      integer :: i

#ifndef NDEBUG
      if (size(dat, _HORIZONTAL_DIMENSION_COUNT_ + 1) /= size(self%bottom_state_variables)) &
         call fatal_error('link_all_bottom_state_data', 'size of last dimension of provided array must match number of bottom state variables.')
#endif
      do i = 1, size(self%bottom_state_variables)
         call link_bottom_state_data(self, i, dat(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ i))
      end do
   end subroutine link_all_bottom_state_data

   subroutine link_all_surface_state_data(self, dat)
      class (type_fabm_model),                                intent(inout) :: self
      real(rke) _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_, target, intent(in)    :: dat

      integer :: i

#ifndef NDEBUG
      if (size(dat, _HORIZONTAL_DIMENSION_COUNT_ + 1) /= size(self%surface_state_variables)) &
         call fatal_error('link_all_surface_state_data', 'size of last dimension of provided array must match number of surface state variables.')
#endif
      do i = 1, size(self%surface_state_variables)
         call link_surface_state_data(self, i, dat(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ i))
      end do
   end subroutine link_all_surface_state_data

   function get_interior_diagnostic_data(self, index) result(dat)
      class (type_fabm_model), intent(in)    :: self
      integer,                 intent(in)    :: index
      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: dat

      _ASSERT_(self%status >= status_start_done, 'get_interior_diagnostic_data', 'This routine can only be called after model start.')
      dat => null()
      if (self%interior_diagnostic_variables(index)%target%catalog_index /= -1) &
         dat => self%catalog%interior(self%interior_diagnostic_variables(index)%target%catalog_index)%p
   end function get_interior_diagnostic_data

   function get_horizontal_diagnostic_data(self, index) result(dat)
      class (type_fabm_model), intent(in)               :: self
      integer,                 intent(in)               :: index
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: dat

      _ASSERT_(self%status >= status_start_done, 'get_horizontal_diagnostic_data', 'This routine can only be called after model start.')
      dat => null()
      if (self%horizontal_diagnostic_variables(index)%target%catalog_index /= -1) &
         dat => self%catalog%horizontal(self%horizontal_diagnostic_variables(index)%target%catalog_index)%p
   end function get_horizontal_diagnostic_data

   function get_interior_data(self, id) result(dat)
      class (type_fabm_model), target,       intent(in) :: self
      type (type_fabm_interior_variable_id), intent(in) :: id
      real(rke) _ATTRIBUTES_GLOBAL_, pointer            :: dat

      _ASSERT_(self%status >= status_start_done, 'get_interior_data', 'This routine can only be called after model start.')
      dat => null()
      if (.not. associated(id%variable)) return
      if (id%variable%catalog_index /= -1) dat => self%catalog%interior(id%variable%catalog_index)%p
   end function get_interior_data

   function get_horizontal_data(self, id) result(dat)
      class (type_fabm_model), target,        intent(in) :: self
      type(type_fabm_horizontal_variable_id), intent(in) :: id
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer  :: dat

      _ASSERT_(self%status >= status_start_done, 'get_horizontal_data', 'This routine can only be called after model start.')
      dat => null()
      if (.not. associated(id%variable)) return
      if (id%variable%catalog_index /= -1) dat => self%catalog%horizontal(id%variable%catalog_index)%p
   end function get_horizontal_data

   function get_scalar_data(self, id) result(dat)
      class (type_fabm_model), target,     intent(in) :: self
      type (type_fabm_scalar_variable_id), intent(in) :: id
      real(rke), pointer                              :: dat

      _ASSERT_(self%status >= status_start_done, 'get_scalar_data', 'This routine can only be called after model start.')
      dat => null()
      if (.not. associated(id%variable)) return
      if (id%variable%catalog_index /= -1) dat => self%catalog%scalar(id%variable%catalog_index)%p
   end function get_scalar_data

   subroutine initialize_interior_state(self _POSTARG_INTERIOR_IN_)
      class (type_fabm_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_

      integer :: ivar, read_index, icall
      _DECLARE_INTERIOR_INDICES_

#ifndef NDEBUG
      call check_interior_location(self%domain%start, self%domain%stop _POSTARG_INTERIOR_IN_, 'initialize_interior_state')
#endif

      call cache_pack(self%domain, self%catalog, self%cache_fill_values, self%initialize_interior_state_job%first_task, self%cache_int _POSTARG_INTERIOR_IN_)

      ! Default initialization for interior state variables
      do ivar = 1, size(self%interior_state_variables)
         read_index = self%interior_state_variables(ivar)%target%read_indices%value
         _CONCURRENT_LOOP_BEGIN_EX_(self%cache_int)
            self%cache_int%read _INDEX_SLICE_PLUS_1_(read_index) = self%interior_state_variables(ivar)%initial_value
         _LOOP_END_
      end do

      ! Allow biogeochemical models to initialize their interior state.
      do icall = 1, size(self%initialize_interior_state_job%first_task%calls)
         if (self%initialize_interior_state_job%first_task%calls(icall)%source == source_initialize_state) call self%initialize_interior_state_job%first_task%calls(icall)%model%initialize_state(self%cache_int)
      end do

      ! Copy from cache back to global data store [NB variable values have been set in the *read* cache].
      do ivar = 1, size(self%interior_state_variables)
         read_index = self%interior_state_variables(ivar)%target%read_indices%value
         if (self%catalog%interior_sources(read_index) == data_source_fabm) then
            _UNPACK_TO_GLOBAL_(self%cache_int%read, read_index, self%catalog%interior(self%interior_state_variables(ivar)%target%catalog_index)%p, self%cache_int, self%interior_state_variables(ivar)%missing_value)
         end if
      end do

      call cache_unpack(self%initialize_interior_state_job%first_task, self%cache_int, self%store _POSTARG_INTERIOR_IN_)
   end subroutine initialize_interior_state

   subroutine initialize_bottom_state(self _POSTARG_HORIZONTAL_IN_)
      class (type_fabm_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_

      integer :: icall, ivar, read_index
      _DECLARE_HORIZONTAL_INDICES_

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'initialize_bottom_state')
#endif

      call cache_pack(self%domain, self%catalog, self%cache_fill_values, self%initialize_bottom_state_job%first_task, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      ! Default initialization for bottom state variables
      do ivar = 1, size(self%bottom_state_variables)
         read_index = self%bottom_state_variables(ivar)%target%read_indices%value
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%cache_hz)
            self%cache_hz%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index) = self%bottom_state_variables(ivar)%initial_value
         _HORIZONTAL_LOOP_END_
      end do

      ! Allow biogeochemical models to initialize their bottom state.
      do icall = 1, size(self%initialize_bottom_state_job%first_task%calls)
         if (self%initialize_bottom_state_job%first_task%calls(icall)%source == source_initialize_bottom_state) call self%initialize_bottom_state_job%first_task%calls(icall)%model%initialize_bottom_state(self%cache_hz)
      end do

      ! Copy from cache back to global data store [NB variable values have been set in the *read* cache].
      do ivar = 1, size(self%bottom_state_variables)
         read_index = self%bottom_state_variables(ivar)%target%read_indices%value
         if (self%catalog%horizontal_sources(read_index) == data_source_fabm) then
            _HORIZONTAL_UNPACK_TO_GLOBAL_(self%cache_hz%read_hz, read_index, self%catalog%horizontal(self%bottom_state_variables(ivar)%target%catalog_index)%p, self%cache_hz, self%bottom_state_variables(ivar)%missing_value)
         end if
      end do

      call cache_unpack(self%initialize_bottom_state_job%first_task, self%cache_hz, self%store _POSTARG_HORIZONTAL_IN_)
   end subroutine initialize_bottom_state

   subroutine initialize_surface_state(self _POSTARG_HORIZONTAL_IN_)
      class (type_fabm_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_

      integer :: icall, ivar, read_index
      _DECLARE_HORIZONTAL_INDICES_

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'initialize_surface_state')
#endif

      call cache_pack(self%domain, self%catalog, self%cache_fill_values, self%initialize_surface_state_job%first_task, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      ! Default initialization for surface state variables
      do ivar = 1, size(self%surface_state_variables)
         read_index = self%surface_state_variables(ivar)%target%read_indices%value
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%cache_hz)
            self%cache_hz%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index) = self%surface_state_variables(ivar)%initial_value
         _HORIZONTAL_LOOP_END_
      end do

      ! Allow biogeochemical models to initialize their surface state.
      do icall = 1, size(self%initialize_surface_state_job%first_task%calls)
         if (self%initialize_surface_state_job%first_task%calls(icall)%source == source_initialize_surface_state) call self%initialize_surface_state_job%first_task%calls(icall)%model%initialize_surface_state(self%cache_hz)
      end do

      ! Copy from cache back to global data store [NB variable values have been set in the *read* cache].
      do ivar = 1, size(self%surface_state_variables)
         read_index = self%surface_state_variables(ivar)%target%read_indices%value
         if (self%catalog%horizontal_sources(read_index) == data_source_fabm) then
            _HORIZONTAL_UNPACK_TO_GLOBAL_(self%cache_hz%read_hz, read_index, self%catalog%horizontal(self%surface_state_variables(ivar)%target%catalog_index)%p, self%cache_hz, self%surface_state_variables(ivar)%missing_value)
         end if
      end do

      call cache_unpack(self%initialize_surface_state_job%first_task, self%cache_hz, self%store _POSTARG_HORIZONTAL_IN_)
   end subroutine initialize_surface_state

   subroutine get_interior_sources_rhs(self _POSTARG_INTERIOR_IN_, dy)
      class (type_fabm_model),                intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_1_, intent(inout) :: dy

      integer :: i, k
      _DECLARE_INTERIOR_INDICES_

#ifndef NDEBUG
      call check_interior_location(self%domain%start, self%domain%stop _POSTARG_INTERIOR_IN_, 'get_interior_sources_rhs')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_extents_2d(dy, _STOP_ - _START_ + 1, size(self%interior_state_variables), 'get_interior_sources_rhs', 'dy', 'stop-start+1, # interior state variables')
#  else
      call check_extents_1d(dy, size(self%interior_state_variables), 'get_interior_sources_rhs', 'dy', '# interior state variables')
#  endif
#endif

      call process_interior_slice(self%get_interior_sources_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_int _POSTARG_INTERIOR_IN_)

      ! Compose total sources-sinks for each state variable, combining model-specific contributions.
      do i = 1, size(self%get_interior_sources_job%arg1_sources)
         k = self%get_interior_sources_job%arg1_sources(i)
         _UNPACK_AND_ADD_TO_PLUS_1_(self%cache_int%write, k, dy, i, self%cache_int)
      end do
   end subroutine get_interior_sources_rhs

   subroutine get_interior_sources_ppdd(self _POSTARG_INTERIOR_IN_, pp, dd)
      class (type_fabm_model),                intent(inout) :: self
     _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_2_, intent(inout) :: pp, dd

      integer :: icall, i, j, k, ncopy
      _DECLARE_INTERIOR_INDICES_

#ifndef NDEBUG
      call check_interior_location(self%domain%start, self%domain%stop _POSTARG_INTERIOR_IN_, 'get_interior_sources_ppdd')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_extents_3d(pp, _STOP_ - _START_ + 1, size(self%interior_state_variables), size(self%interior_state_variables), 'get_interior_sources_ppdd', 'pp', 'stop-start+1, # interior state variables, # interior state variables')
      call check_extents_3d(dd, _STOP_ - _START_ + 1, size(self%interior_state_variables), size(self%interior_state_variables), 'get_interior_sources_ppdd', 'dd', 'stop-start+1, # interior state variables, # interior state variables')
#  else
      call check_extents_2d(pp, size(self%interior_state_variables), size(self%interior_state_variables), 'get_interior_sources_ppdd', 'pp', '# interior state variables, # interior state variables')
      call check_extents_2d(dd, size(self%interior_state_variables), size(self%interior_state_variables), 'get_interior_sources_ppdd', 'dd', '# interior state variables, # interior state variables')
#  endif
#endif

      call cache_pack(self%domain, self%catalog, self%cache_fill_values, self%get_interior_sources_job%first_task, self%cache_int _POSTARG_INTERIOR_IN_)

      ncopy = 0
      do icall = 1, size(self%get_interior_sources_job%first_task%calls)
         call self%get_interior_sources_job%first_task%calls(icall)%model%do_ppdd(self%cache_int, pp, dd)

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1 + ncopy,self%get_interior_sources_job%first_task%calls(icall)%ncopy_int + ncopy)
            j = self%get_interior_sources_job%first_task%copy_commands_int(i)%read_index
            k = self%get_interior_sources_job%first_task%copy_commands_int(i)%write_index
            _CONCURRENT_LOOP_BEGIN_EX_(self%cache_int)
               self%cache_int%read _INDEX_SLICE_PLUS_1_(j) = self%cache_int%write _INDEX_SLICE_PLUS_1_(k)
            _LOOP_END_
         end do
         ncopy = ncopy + self%get_interior_sources_job%first_task%calls(icall)%ncopy_int
      end do

      call cache_unpack(self%get_interior_sources_job%first_task, self%cache_int, self%store _POSTARG_INTERIOR_IN_)
   end subroutine get_interior_sources_ppdd

   subroutine check_interior_state(self _POSTARG_INTERIOR_IN_, repair, valid)
      class (type_fabm_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      logical,                 intent(in)    :: repair
      logical,                 intent(out)   :: valid

      integer            :: ivar, read_index
      real(rki)          :: value, minimum, maximum
      character(len=256) :: err
      _DECLARE_INTERIOR_INDICES_

#ifndef NDEBUG
      call check_interior_location(self%domain%start, self%domain%stop _POSTARG_INTERIOR_IN_, 'check_interior_state')
#endif

      self%cache_int%repair = repair
      self%cache_int%valid = .true.
      self%cache_int%set_interior = .false.

      call process_interior_slice(self%check_interior_state_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_int _POSTARG_INTERIOR_IN_)

      valid = self%cache_int%valid
      if (.not. (valid .or. repair)) return

      ! Finally check whether all state variable values lie within their prescribed [constant] bounds.
      ! This is always done, independently of any model-specific checks that may have been called above.

      ! Quick bounds check for the common case where all values are valid.
      do ivar = 1, size(self%interior_state_variables)
         read_index = self%interior_state_variables(ivar)%target%read_indices%value
         minimum = self%interior_state_variables(ivar)%target%minimum
         maximum = self%interior_state_variables(ivar)%target%maximum
         _LOOP_BEGIN_EX_(self%cache_int)
            value = self%cache_int%read _INDEX_SLICE_PLUS_1_(read_index)
            if (value < minimum .or. value > maximum) valid = .false.
         _LOOP_END_
      end do

      if (.not. valid) then
         ! Check boundaries for pelagic state variables specified by the models.
         ! If repair is permitted, this clips invalid values to the closest boundary.
         do ivar = 1, size(self%interior_state_variables)
            ! Shortcuts to variable information - this demonstrably helps the compiler (ifort).
            read_index = self%interior_state_variables(ivar)%target%read_indices%value
            minimum = self%interior_state_variables(ivar)%target%minimum
            maximum = self%interior_state_variables(ivar)%target%maximum

            if (repair) then
               _CONCURRENT_LOOP_BEGIN_EX_(self%cache_int)
                  value = self%cache_int%read _INDEX_SLICE_PLUS_1_(read_index)
                  self%cache_int%read _INDEX_SLICE_PLUS_1_(read_index) = max(minimum, min(maximum, value))
               _LOOP_END_
            else
               _LOOP_BEGIN_EX_(self%cache_int)
                  value = self%cache_int%read _INDEX_SLICE_PLUS_1_(read_index)
                  if (value < minimum) then
                     ! State variable value lies below prescribed minimum.
                     write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ',trim(self%interior_state_variables(ivar)%name), &
                                                                & ' below minimum value ',minimum
                     call log_message(err)
                     return
                  elseif (value > maximum) then
                     ! State variable value exceeds prescribed maximum.
                     write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ',trim(self%interior_state_variables(ivar)%name), &
                                                                & ' above maximum value ',maximum
                     call log_message(err)
                     return
                  end if
               _LOOP_END_
            end if
         end do
      end if

      if (self%cache_int%set_interior .or. .not. valid) then
         do ivar = 1, size(self%interior_state_variables)
            read_index = self%interior_state_variables(ivar)%target%read_indices%value
            if (self%catalog%interior_sources(read_index) == data_source_fabm) then
               _UNPACK_TO_GLOBAL_(self%cache_int%read, read_index, self%catalog%interior(self%interior_state_variables(ivar)%target%catalog_index)%p, self%cache_int, self%interior_state_variables(ivar)%missing_value)
            end if
         end do
      end if
   end subroutine check_interior_state

   subroutine check_bottom_state(self _POSTARG_HORIZONTAL_IN_, repair, valid)
      class (type_fabm_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      logical,                 intent(in)    :: repair
      logical,                 intent(out)   :: valid

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'check_bottom_state')
#endif

      call internal_check_horizontal_state(self, self%check_bottom_state_job _POSTARG_HORIZONTAL_IN_, 2, self%bottom_state_variables, repair, valid)
   end subroutine check_bottom_state

   subroutine check_surface_state(self _POSTARG_HORIZONTAL_IN_, repair, valid)
      class (type_fabm_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      logical,                 intent(in)    :: repair
      logical,                 intent(out)   :: valid

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'check_surface_state')
#endif

      call internal_check_horizontal_state(self, self%check_surface_state_job _POSTARG_HORIZONTAL_IN_, 1, self%surface_state_variables, repair, valid)
   end subroutine check_surface_state

   subroutine internal_check_horizontal_state(self,job _POSTARG_HORIZONTAL_IN_, flag, state_variables, repair, valid)
      class (type_fabm_model),                    intent(inout) :: self
      type (type_job),                            intent(in)    :: job
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      integer,                                    intent(in)    :: flag
      type (type_fabm_horizontal_state_variable), intent(inout) :: state_variables(:)
      logical,                                    intent(in)    :: repair
      logical,                                    intent(out)   :: valid

      integer            :: ivar, read_index
      real(rki)          :: value, minimum, maximum
      character(len=256) :: err
      _DECLARE_HORIZONTAL_INDICES_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      integer :: _VERTICAL_ITERATOR_
#endif

      self%cache_hz%repair = repair
      self%cache_hz%valid = .true.
      self%cache_hz%set_horizontal = .false.
      self%cache_hz%set_interior = .false.

      call process_horizontal_slice(job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      valid = self%cache_hz%valid
      if (.not. (valid .or. repair)) return

      ! Check boundaries for horizontal state variables, as prescribed by the owning models.
      ! If repair is permitted, this clips invalid values to the closest boundary.
      do ivar = 1, size(state_variables)
         ! Shortcuts to variable information - this demonstrably helps the compiler (ifort).
         read_index = state_variables(ivar)%target%read_indices%value
         minimum = state_variables(ivar)%target%minimum
         maximum = state_variables(ivar)%target%maximum

         _HORIZONTAL_LOOP_BEGIN_EX_(self%cache_hz)
            value = self%cache_hz%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index)
            if (value < minimum) then
               ! State variable value lies below prescribed minimum.
               valid = .false.
               if (.not. repair) then
                  write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ', &
                                                             & trim(state_variables(ivar)%name), &
                                                             & ' below minimum value ',minimum
                  call log_message(err)
                  return
               end if
               self%cache_hz%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index) = minimum
            elseif (value > maximum) then
               ! State variable value exceeds prescribed maximum.
               valid = .false.
               if (.not. repair) then
                  write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ', &
                                                             & trim(state_variables(ivar)%name), &
                                                             & ' above maximum value ',maximum
                  call log_message(err)
                  return
               end if
               self%cache_hz%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index) = maximum
            end if
         _HORIZONTAL_LOOP_END_
      end do

      if (self%cache_hz%set_horizontal .or. .not. valid) then
         do ivar = 1, size(state_variables)
            read_index = state_variables(ivar)%target%read_indices%value
            if (self%catalog%horizontal_sources(read_index)==data_source_fabm) then
               _HORIZONTAL_UNPACK_TO_GLOBAL_(self%cache_hz%read_hz, read_index, self%catalog%horizontal(state_variables(ivar)%target%catalog_index)%p, self%cache_hz, state_variables(ivar)%missing_value)
            end if
         end do
      end if

      if (self%cache_hz%set_interior) then
         ! One or more models have provided new values for an interior state variable [at the interface]

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
         if (flag == 1) then
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
            _VERTICAL_ITERATOR_ = self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)
#  else
            _VERTICAL_ITERATOR_ = self%domain%start(_FABM_DEPTH_DIMENSION_INDEX_)
#  endif
         else
#  if _FABM_BOTTOM_INDEX_==0
#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
            _VERTICAL_ITERATOR_ = self%domain%start(_FABM_DEPTH_DIMENSION_INDEX_)
#    else
            _VERTICAL_ITERATOR_ = self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)
#    endif
#  elif !defined(_HORIZONTAL_IS_VECTORIZED_)
            _VERTICAL_ITERATOR_ = self%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#  endif
         end if
#endif

         do ivar = 1, size(self%interior_state_variables)
            read_index = self%interior_state_variables(ivar)%target%read_indices%value
            if (self%catalog%interior_sources(read_index) == data_source_fabm) then
#if _FABM_BOTTOM_INDEX_==-1&&defined(_HORIZONTAL_IS_VECTORIZED_)
               if (flag == 1) then
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
#  ifdef _HAS_MASK_
                  self%catalog%interior(self%interior_state_variables(ivar)%target%catalog_index)%p _INDEX_GLOBAL_INTERIOR_(self%cache_hz%ipack(1:self%cache_hz%n)) = self%cache_hz%read(1:self%cache_hz%n, read_index)
#  else
                  _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%cache_hz)
                     self%catalog%interior(self%interior_state_variables(ivar)%target%catalog_index)%p _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1) = self%cache_hz%read _INDEX_SLICE_PLUS_1_(read_index)
                  _HORIZONTAL_LOOP_END_
#  endif
#elif defined(_INTERIOR_IS_VECTORIZED_)
                  self%catalog%interior(self%interior_state_variables(ivar)%target%catalog_index)%p _INDEX_LOCATION_ = self%cache_hz%read(1,read_index)
#else
                  self%catalog%interior(self%interior_state_variables(ivar)%target%catalog_index)%p _INDEX_LOCATION_ = self%cache_hz%read(read_index)
#endif

#if _FABM_BOTTOM_INDEX_==-1&&defined(_HORIZONTAL_IS_VECTORIZED_)
               else
                  ! Special case for bottom if vertical index of bottom point is variable.
                  _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%cache_hz)
#  ifdef _HAS_MASK_
                     _VERTICAL_ITERATOR_ = self%domain%bottom_indices _INDEX_GLOBAL_HORIZONTAL_(self%cache_hz%ipack(_J_))
                     self%catalog%interior(self%interior_state_variables(ivar)%target%catalog_index)%p _INDEX_GLOBAL_INTERIOR_(self%cache_hz%ipack(_J_)) = self%cache_hz%read _INDEX_SLICE_PLUS_1_(read_index)
#  else
                     _VERTICAL_ITERATOR_ = self%domain%bottom_indices _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1)
                     self%catalog%interior(self%interior_state_variables(ivar)%target%catalog_index)%p _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1) = self%cache_hz%read _INDEX_SLICE_PLUS_1_(read_index)
#  endif
                  _HORIZONTAL_LOOP_END_
               end if
#endif
            end if
         end do
      end if
   end subroutine internal_check_horizontal_state

   subroutine get_surface_sources(self _POSTARG_HORIZONTAL_IN_, flux_pel, flux_sf)
      class (type_fabm_model),                           intent(inout)         :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(out)           :: flux_pel
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(out), optional :: flux_sf

      integer :: i, k
      _DECLARE_HORIZONTAL_INDICES_

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'get_surface_sources')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_extents_2d(flux_pel, _STOP_ - _START_ + 1, size(self%interior_state_variables), 'get_surface_sources', 'flux_pel', 'stop-start+1, # interior state variables')
      if (present(flux_sf)) call check_extents_2d(flux_sf, _STOP_ - _START_ + 1, size(self%surface_state_variables), 'get_surface_sources', 'flux_sf', 'stop-start+1, # surface state variables')
#  else
      call check_extents_1d(flux_pel, size(self%interior_state_variables), 'get_surface_sources', 'flux_pel', '# interior state variables')
      if (present(flux_sf)) call check_extents_1d(flux_sf, size(self%surface_state_variables), 'get_surface_sources', 'flux_sf', '# surface state variables')
#  endif
#endif

      call process_horizontal_slice(self%get_surface_sources_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      ! Compose surface fluxes for each interior state variable, combining model-specific contributions.
      flux_pel = 0.0_rke
      do i = 1, size(self%get_surface_sources_job%arg1_sources)
         k = self%get_surface_sources_job%arg1_sources(i)
         _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(self%cache_hz%write_hz, k, flux_pel, i, self%cache_hz)
      end do

      ! Compose total sources-sinks for each surface-bound state variable, combining model-specific contributions.
      if (present(flux_sf)) then
         flux_sf = 0.0_rke
         do i = 1, size(self%get_surface_sources_job%arg2_sources)
            k = self%get_surface_sources_job%arg2_sources(i)
            _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(self%cache_hz%write_hz, k, flux_sf, i, self%cache_hz)
         end do
      end if
   end subroutine get_surface_sources

   subroutine get_bottom_sources_rhs(self _POSTARG_HORIZONTAL_IN_, flux_pel, flux_ben)
      class (type_fabm_model),                           intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(inout) :: flux_pel, flux_ben

      integer :: i, k
      _DECLARE_HORIZONTAL_INDICES_

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'get_bottom_sources_rhs')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_extents_2d(flux_pel, _STOP_ - _START_ + 1, size(self%interior_state_variables), 'get_bottom_sources_rhs', 'flux_pel', 'stop-start+1, # interior state variables')
      call check_extents_2d(flux_ben, _STOP_ - _START_ + 1, size(self%bottom_state_variables), 'get_bottom_sources_rhs', 'flux_ben', 'stop-start+1, # bottom state variables')
#  else
      call check_extents_1d(flux_pel, size(self%interior_state_variables), 'get_bottom_sources_rhs', 'flux_pel', '# interior state variables')
      call check_extents_1d(flux_ben, size(self%bottom_state_variables), 'get_bottom_sources_rhs', 'flux_ben', '# bottom state variables')
#  endif
#endif

      call process_horizontal_slice(self%get_bottom_sources_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      ! Compose bottom fluxes for each interior state variable, combining model-specific contributions.
      do i = 1, size(self%get_bottom_sources_job%arg1_sources)
         k = self%get_bottom_sources_job%arg1_sources(i)
         _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(self%cache_hz%write_hz, k, flux_pel, i, self%cache_hz)
      end do

      ! Compose total sources-sinks for each bottom-bound state variable, combining model-specific contributions.
      do i = 1, size(self%get_bottom_sources_job%arg2_sources)
         k = self%get_bottom_sources_job%arg2_sources(i)
         _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(self%cache_hz%write_hz, k, flux_ben, i, self%cache_hz)
      end do
   end subroutine get_bottom_sources_rhs

   subroutine get_bottom_sources_ppdd(self _POSTARG_HORIZONTAL_IN_, pp, dd, benthos_offset)
      class (type_fabm_model),                           intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      integer,                                           intent(in)    :: benthos_offset
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_2_, intent(inout) :: pp, dd

      integer                   :: icall, i, j, k, ncopy
      _DECLARE_HORIZONTAL_INDICES_

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'get_bottom_sources_ppdd')
#endif

      call cache_pack(self%domain, self%catalog, self%cache_fill_values, self%get_bottom_sources_job%first_task, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      ncopy = 0
      do icall = 1, size(self%get_bottom_sources_job%first_task%calls)
         if (self%get_bottom_sources_job%first_task%calls(icall)%source == source_do_bottom) call self%get_bottom_sources_job%first_task%calls(icall)%model%do_bottom_ppdd(self%cache_hz, pp, dd, benthos_offset)

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1 + ncopy,self%get_bottom_sources_job%first_task%calls(icall)%ncopy_hz + ncopy)
            j = self%get_bottom_sources_job%first_task%copy_commands_hz(i)%read_index
            k = self%get_bottom_sources_job%first_task%copy_commands_hz(i)%write_index
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%cache_hz)
               self%cache_hz%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = self%cache_hz%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
            _HORIZONTAL_LOOP_END_
         end do
         ncopy = ncopy + self%get_bottom_sources_job%first_task%calls(icall)%ncopy_hz
      end do

      call cache_unpack(self%get_bottom_sources_job%first_task, self%cache_hz, self%store _POSTARG_HORIZONTAL_IN_)
   end subroutine get_bottom_sources_ppdd

   subroutine get_vertical_movement(self _POSTARG_INTERIOR_IN_, velocity)
      class (type_fabm_model),                intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_1_, intent(out)   :: velocity

      integer :: i, k
      _DECLARE_INTERIOR_INDICES_

#ifndef NDEBUG
      call check_interior_location(self%domain%start, self%domain%stop _POSTARG_INTERIOR_IN_, 'get_vertical_movement')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_extents_2d(velocity, _STOP_ - _START_ + 1, size(self%interior_state_variables), 'get_vertical_movement', 'velocity', 'stop-start+1, # interior state variables')
#  else
      call check_extents_1d(velocity, size(self%interior_state_variables), 'get_vertical_movement', 'velocity', '# interior state variables')
#  endif
#endif

      call process_interior_slice(self%get_vertical_movement_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_int _POSTARG_INTERIOR_IN_)

      ! Copy vertical velocities from write cache to output array provided by host
      do i = 1, size(self%get_vertical_movement_job%arg1_sources)
         k = self%get_vertical_movement_job%arg1_sources(i)
         _UNPACK_TO_PLUS_1_(self%cache_int%write, k, velocity, i, self%cache_int, 0.0_rke)
      end do
   end subroutine get_vertical_movement

   subroutine get_interior_conserved_quantities(self _POSTARG_INTERIOR_IN_, sums)
      class (type_fabm_model),                intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_1_, intent(out)   :: sums

      integer :: i
      _DECLARE_INTERIOR_INDICES_

#ifndef NDEBUG
      call check_interior_location(self%domain%start, self%domain%stop _POSTARG_INTERIOR_IN_, 'get_interior_conserved_quantities')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_extents_2d(sums, _STOP_ - _START_ + 1, size(self%conserved_quantities), 'get_interior_conserved_quantities', 'sums', 'stop-start+1, # conserved quantities')
#  else
      call check_extents_1d(sums, size(self%conserved_quantities), 'get_interior_conserved_quantities', 'sums', '# conserved quantities')
#  endif
#endif

      call process_interior_slice(self%get_interior_conserved_quantities_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_int _POSTARG_INTERIOR_IN_)

      do i = 1, size(self%conserved_quantities)
         _UNPACK_TO_PLUS_1_(self%cache_int%write, self%conserved_quantities(i)%index, sums, i, self%cache_int, 0.0_rke)
      end do
   end subroutine get_interior_conserved_quantities

   subroutine get_horizontal_conserved_quantities(self _POSTARG_HORIZONTAL_IN_, sums)
      class (type_fabm_model),                           intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(out)   :: sums

      integer :: i
      _DECLARE_HORIZONTAL_INDICES_

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'get_horizontal_conserved_quantities')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_extents_2d(sums, _STOP_ - _START_ + 1, size(self%conserved_quantities), 'get_horizontal_conserved_quantities', 'sums', 'stop-start+1, # conserved quantities')
#  else
      call check_extents_1d(sums, size(self%conserved_quantities), 'get_horizontal_conserved_quantities', 'sums', '# conserved quantities')
#  endif
#endif

      call process_horizontal_slice(self%get_horizontal_conserved_quantities_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      do i = 1, size(self%conserved_quantities)
         _HORIZONTAL_UNPACK_TO_PLUS_1_(self%cache_hz%write_hz, self%conserved_quantities(i)%horizontal_index, sums, i, self%cache_hz, 0.0_rke)
      end do
   end subroutine get_horizontal_conserved_quantities

   subroutine process_job(self, job _ARGUMENTS_HORIZONTAL_LOCATION_RANGE_)
      class (type_fabm_model), intent(inout), target :: self
      type (type_job),         intent(in)            :: job
      _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_RANGE_

      integer                   :: i
      type (type_task), pointer :: task
      _DECLARE_LOCATION_

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      ! Jobs must be applied across the entire depth range (if any),
      ! so we set vertical start and stop indices here.
      integer :: _VERTICAL_START_, _VERTICAL_STOP_
      _VERTICAL_START_ = self%domain%start(_FABM_DEPTH_DIMENSION_INDEX_)
      _VERTICAL_STOP_ = self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)
#endif

      do i = 1, size(job%interior_store_prefill)
         if (job%interior_store_prefill(i)) then
            _BEGIN_OUTER_INTERIOR_LOOP_
               self%store%interior _INDEX_GLOBAL_INTERIOR_PLUS_1_(_START_:_STOP_, i) = self%store%interior_fill_value(i)
            _END_OUTER_INTERIOR_LOOP_
         end if
      end do
      do i = 1, size(job%horizontal_store_prefill)
         if (job%horizontal_store_prefill(i)) then
            _BEGIN_OUTER_HORIZONTAL_LOOP_
               self%store%horizontal _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_:_STOP_, i) = self%store%horizontal_fill_value(i)
            _END_OUTER_HORIZONTAL_LOOP_
         end if
      end do

      task => job%first_task
      do while (associated(task))
         select case (task%operation)
         case (source_do)
            _BEGIN_OUTER_INTERIOR_LOOP_
               call process_interior_slice(task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_int _POSTARG_INTERIOR_IN_)
            _END_OUTER_INTERIOR_LOOP_
         case (source_do_surface, source_do_bottom, source_do_horizontal)
            _BEGIN_OUTER_HORIZONTAL_LOOP_
               call process_horizontal_slice(task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_hz _POSTARG_HORIZONTAL_IN_)
            _END_OUTER_HORIZONTAL_LOOP_
         case (source_do_column)
            _BEGIN_OUTER_VERTICAL_LOOP_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==-1
#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
               _VERTICAL_START_ = self%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#    else
               _VERTICAL_STOP_ = self%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#    endif
#  endif
#endif
               if (_IS_UNMASKED_(self%domain%mask_hz _INDEX_HORIZONTAL_LOCATION_)) call process_vertical_slice(task, self%domain, &
                  self%catalog, self%cache_fill_values, self%store, self%cache_vert _POSTARG_VERTICAL_IN_)
            _END_OUTER_VERTICAL_LOOP_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
            _VERTICAL_START_ = self%domain%start(_FABM_DEPTH_DIMENSION_INDEX_)
            _VERTICAL_STOP_ = self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)
#endif
         end select
         task => task%next
      end do
   end subroutine process_job

#if _FABM_DIMENSION_COUNT_>1 || (_FABM_DIMENSION_COUNT_==1 && !defined(_FABM_DEPTH_DIMENSION_INDEX_))
   subroutine process_job_everywhere(self, job)
      class (type_fabm_model), intent(inout), target :: self
      type (type_job),         intent(in)            :: job
      integer :: _LOCATION_RANGE_
      istart__ = self%domain%start(1)
      istop__ = self%domain%stop(1)
#  if _FABM_DIMENSION_COUNT_ > 1
      jstart__ = self%domain%start(2)
      jstop__ = self%domain%stop(2)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 2
      kstart__ = self%domain%start(3)
      kstop__ = self%domain%stop(3)
#  endif
      call process_job(self, job _ARGUMENTS_HORIZONTAL_LOCATION_RANGE_)
   end subroutine process_job_everywhere
#endif

   subroutine prepare_inputs1(self, t)
      class (type_fabm_model),  intent(inout) :: self
      real(rke), optional,      intent(in)    :: t

      class (type_expression), pointer :: expression

      call self%process(self%prepare_inputs_job)

      if (.not. present(t)) return

      expression => self%root%first_expression
      do while (associated(expression))
         select type (expression)
         class is (type_interior_temporal_mean)
            call update_interior_temporal_mean(expression)
         class is (type_horizontal_temporal_mean)
            call update_horizontal_temporal_mean(expression)
         end select
         expression => expression%next
      end do

   contains

      subroutine update_interior_temporal_mean(expression)
         class (type_interior_temporal_mean), intent(inout) :: expression
         integer  :: i
         real(rke) :: weight_right, frac_outside

         if (expression%ioldest == -1) then
            ! Start of simulation
            expression%next_save_time = t + expression%period / expression%n
            expression%ioldest = 1
         end if
         do while (t >= expression%next_save_time)
            ! Weight for linear interpolation between last stored point and current point, to get at values for desired time.
            weight_right = (expression%next_save_time - expression%last_time) / (t - expression%last_time)

            ! For temporal means:
            ! - remove contribution of oldest point from historical mean (@ n + 2)
            ! - linearly interpolate to desired time (@ ioldest), by computing a weighted mean of the current value (data(expression%in)%p) and the previous value (@ n + 1)
            ! - add contribution of new point to historical mean
            expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) &
               - expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest) / expression%n
            expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest) = (1.0_rke - weight_right) * expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 1) &
               + weight_right * self%catalog%interior(expression%in)%p
            expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) &
               + expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest) / expression%n

            ! Compute next time for which we want to store output
            expression%next_save_time = expression%next_save_time + expression%period / expression%n

            ! If we just completed the first entire history, compute the running mean and record that it is now valid.
            if (expression%ioldest == expression%n .and. .not. expression%valid) then
               expression%valid = .true.
               expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) = 0.0_rke
               do i = 1, expression%n
                  expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) + expression%history(_PREARG_LOCATION_DIMENSIONS_ i)
               end do
               expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) / expression%n
            end if

            ! Increment index for oldest time point
            expression%ioldest = expression%ioldest + 1
            if (expression%ioldest > expression%n) expression%ioldest = 1
         end do

         ! Store current value to enable linear interpolation to next output time in subsequent call.
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 1) = self%catalog%interior(expression%in)%p

         if (expression%valid) then
            ! We have a full history. To compute the temporal mean:
            ! - store values at current time step (@ n + 1)
            ! - set mean (@ n + 3) to historical mean (@ n + 2) but account for change since most recent point in history.

            ! Compute extent of time period outside history
            frac_outside = (t- (expression%next_save_time - expression%period / expression%n)) / expression%period

            ! Set corrected running mean (move window by removing part of the start, and appending to the end)
            expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 3) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 2) &
               + frac_outside*(-expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest) + expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+1))
         else
            ! We do not have a full history yet; set temporal mean to missing value
            expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n + 3) = expression%missing_value
         end if

         expression%last_time = t
      end subroutine

      subroutine update_horizontal_temporal_mean(expression)
         class (type_horizontal_temporal_mean), intent(inout) :: expression
         integer  :: i
         real(rke) :: weight_right, frac_outside

         if (expression%ioldest == -1) then
            ! Start of simulation; set entire history equal to current value.
            do i = 1, expression%n + 3
               expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ i) = self%catalog%horizontal(expression%in)%p
            end do
            expression%next_save_time = t + expression%period / expression%n
            expression%ioldest = 1
         end if
         do while (t >= expression%next_save_time)
            ! Weight for linear interpolation between last stored point and current point, to get at values for desired time.
            weight_right = (expression%next_save_time - expression%last_time) / (t - expression%last_time)

            ! For temporal means:
            ! - remove contribution of oldest point from historical mean
            ! - linearly interpolate to desired time
            ! - add contribution of new point to historical mean
            expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 2) = expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 2) &
               - expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%ioldest) / expression%n
            expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%ioldest) = (1.0_rke - weight_right)*expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 1) &
               + weight_right*self%catalog%horizontal(expression%in)%p
            expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 2) = expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 2) &
               + expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%ioldest) / expression%n

            ! Compute next time for which we want to store output
            expression%next_save_time = expression%next_save_time + expression%period / expression%n

            ! Increment index for oldest time point
            expression%ioldest = expression%ioldest + 1
            if (expression%ioldest>expression%n) expression%ioldest = 1
         end do

         ! Compute extent of time period outside history
         frac_outside = (t - (expression%next_save_time - expression%period / expression%n)) / expression%period

         ! For temporal means:
         ! - store values at current time step
         ! - for current mean, use historical mean but account for change since most recent point in history.
         expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 1) = self%catalog%horizontal(expression%in)%p
         expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 3) = expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 2) &
            + frac_outside * (-expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%ioldest) + expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n + 1))

         expression%last_time = t
      end subroutine
   end subroutine prepare_inputs1

   subroutine prepare_inputs2(self, t, year, month, day, seconds)
      class (type_fabm_model), intent(inout) :: self
      real(rke),               intent(in)    :: t
      integer,                 intent(in)    :: year, month, day
      real(rke),               intent(in)    :: seconds

      call self%schedules%update(year, month, day, seconds)
      call prepare_inputs1(self, t)
   end subroutine prepare_inputs2

   subroutine finalize_outputs(self)
      class (type_fabm_model),  intent(inout) :: self

      call self%process(self%finalize_outputs_job)
   end subroutine finalize_outputs

   subroutine classify_variables(self)
      class (type_fabm_model), intent(inout), target :: self

      type (type_link),                                pointer :: link, newlink
      type (type_fabm_interior_state_variable),        pointer :: statevar
      type (type_fabm_horizontal_state_variable),      pointer :: hz_statevar
      type (type_fabm_interior_diagnostic_variable),   pointer :: diagvar
      type (type_fabm_horizontal_diagnostic_variable), pointer :: hz_diagvar
      type (type_fabm_conserved_quantity),             pointer :: consvar
      type (type_internal_variable),                   pointer :: object
      integer                                                  :: nstate, nstate_bot, nstate_surf, ndiag, ndiag_hz, ncons

      type (type_aggregate_variable_list)         :: aggregate_variable_list
      type (type_aggregate_variable),     pointer :: aggregate_variable
      type (type_set)                             :: dependencies, dependencies_hz, dependencies_scalar
      type (type_standard_variable_set)           :: standard_variable_set
      type (type_standard_variable_node), pointer :: standard_variable_node

      ! Build a list of all master variables (those that not have been coupled)
      link => self%root%links%first
      do while (associated(link))
         if (associated(link%target, link%original)) &
            newlink => self%links_postcoupling%append(link%target, link%target%name)
         link => link%next
      end do

      ! Get list of conserved quantities (map to universal=domain-independent variables where possible)
      aggregate_variable_list = collect_aggregate_variables(self%root)
      aggregate_variable => aggregate_variable_list%first
      do while (associated(aggregate_variable))
         if (associated(aggregate_variable%standard_variable%universal)) then
            if (aggregate_variable%standard_variable%universal%conserved) call standard_variable_set%add(aggregate_variable%standard_variable%universal)
         end if
         aggregate_variable => aggregate_variable%next
      end do
      call aggregate_variable_list%finalize()

      ! Count number of conserved quantities and allocate an array for them.
      ncons = 0
      standard_variable_node => standard_variable_set%first
      do while (associated(standard_variable_node))
         ncons = ncons + 1
         standard_variable_node => standard_variable_node%next
      end do
      allocate(self%conserved_quantities(ncons))

      ! Fill list of conserved quantities.
      ! This must be done before building the final authoratitive list of diagnostic variables,
      ! as the calls to append_data_pointer affect the global identifier of diagnostic variables
      ! by adding another pointer that must be set.
      ncons = 0
      standard_variable_node => standard_variable_set%first
      do while (associated(standard_variable_node))
         ncons = ncons + 1
         consvar => self%conserved_quantities(ncons)
         consvar%standard_variable => standard_variable_node%p
         consvar%name = trim(consvar%standard_variable%name)
         consvar%units = trim(consvar%standard_variable%units)
         consvar%long_name = trim(consvar%standard_variable%name)
         consvar%path = trim(consvar%standard_variable%name)
         select type (standard_variable => consvar%standard_variable)
         class is (type_universal_standard_variable)
            consvar%target => get_variable_by_standard_variable(self, standard_variable%in_interior())
            _ASSERT_(associated(consvar%target), 'classify_variables', 'Conserved quantity ' // trim(standard_variable%name) // ' not found in interior.')
            consvar%target_hz => get_variable_by_standard_variable(self, standard_variable%at_interfaces())
            _ASSERT_(associated(consvar%target_hz), 'classify_variables', 'Conserved quantity ' // trim(standard_variable%name) // ' not found at interfaces.')
         end select
         consvar%missing_value = consvar%target%missing_value
         standard_variable_node => standard_variable_node%next
      end do
      call standard_variable_set%finalize()

      ! From this point on, variables will stay as they are.
      ! Coupling is done, and the framework will not add further read indices.

      ! Count number of interior variables in various categories.
      nstate = 0
      ndiag  = 0
      nstate_bot  = 0
      nstate_surf = 0
      ndiag_hz    = 0
      link => self%links_postcoupling%first
      do while (associated(link))
         object => link%target
         _ASSERT_(object%source /= source_state .or. object%write_indices%is_empty(), 'classify_variables', 'variable ' // trim(object%name) // ' has source_state and one or more write indices.')
         select case (object%domain)
         case (domain_interior)
            select case (object%source)
            case (source_unknown)   ! Interior dependency
            case (source_state)     ! Interior state variable
               select case (object%presence)
               case (presence_internal)
                  nstate = nstate + 1
                  call object%state_indices%set_value(nstate)
               case (presence_external_required)
                  call fatal_error('classify_variables', &
                     'Variable ' // trim(link%name) // ' must be coupled to an existing state variable.')
               end select
            case default            ! Interior diagnostic variable
               ndiag = ndiag + 1
            end select
         case (domain_horizontal, domain_surface, domain_bottom)
            select case (object%source)
            case (source_unknown)   ! Horizontal dependency
            case (source_state)     ! Horizontal state variable
               select case (object%presence)
               case (presence_internal)
                  select case (object%domain)
                     case (domain_bottom)
                        nstate_bot = nstate_bot + 1
                        call object%state_indices%set_value(nstate_bot)
                     case (domain_surface)
                        nstate_surf = nstate_surf + 1
                        call object%state_indices%set_value(nstate_surf)
                  end select
               case (presence_external_required)
                  call fatal_error('classify_variables', &
                     'Variable ' // trim(link%name) // ' must be coupled to an existing state variable.')
               end select
            case default            ! Horizontal diagnostic variable
               ndiag_hz = ndiag_hz + 1
            end select
         end select
         link => link%next
      end do

      ! Allocate arrays with variable information that will be accessed by the host model.
      allocate(self%interior_state_variables       (nstate))
      allocate(self%bottom_state_variables         (nstate_bot))
      allocate(self%surface_state_variables        (nstate_surf))
      allocate(self%interior_diagnostic_variables  (ndiag))
      allocate(self%horizontal_diagnostic_variables(ndiag_hz))

      allocate(self%get_interior_sources_job%arg1_sources(nstate))
      allocate(self%get_surface_sources_job%arg1_sources(nstate), self%get_surface_sources_job%arg2_sources(nstate_surf))
      allocate(self%get_bottom_sources_job%arg1_sources(nstate), self%get_bottom_sources_job%arg2_sources(nstate_bot))
      allocate(self%get_vertical_movement_job%arg1_sources(nstate))

      ! Build lists of state variable and diagnostic variables.
      nstate = 0
      ndiag  = 0
      nstate_bot  = 0
      nstate_surf = 0
      ndiag_hz    = 0
      link => self%links_postcoupling%first
      do while (associated(link))
         object => link%target
         select case (link%target%domain)
         case (domain_interior)
            select case (object%source)
            case (source_unknown)   ! Interior dependency
            case (source_state)     ! Interior state variable
               if (object%presence == presence_internal) then
                  nstate = nstate + 1
                  statevar => self%interior_state_variables(nstate)
                  call copy_variable_metadata(object, statevar)
                  if (associated(object%standard_variables%first)) then
                     select type (standard_variable => object%standard_variables%first%p)
                     class is (type_interior_standard_variable)
                        statevar%standard_variable => standard_variable
                     end select
                  end if
                  statevar%initial_value             = object%initial_value
                  statevar%no_precipitation_dilution = object%no_precipitation_dilution
                  statevar%no_river_dilution         = object%no_river_dilution
                  call object%sms_sum%target%write_indices%append(self%get_interior_sources_job%arg1_sources(nstate))
                  call object%surface_flux_sum%target%write_indices%append(self%get_surface_sources_job%arg1_sources(nstate))
                  call object%bottom_flux_sum%target%write_indices%append(self%get_bottom_sources_job%arg1_sources(nstate))
                  call object%movement_sum%target%write_indices%append(self%get_vertical_movement_job%arg1_sources(nstate))
               end if
            case default            ! Interior diagnostic variable
               ndiag = ndiag + 1
               diagvar => self%interior_diagnostic_variables(ndiag)
               call copy_variable_metadata(object, diagvar)
               if (associated(object%standard_variables%first)) then
                  select type (standard_variable => object%standard_variables%first%p)
                  class is (type_interior_standard_variable)
                     diagvar%standard_variable => standard_variable
                  end select
               end if
               diagvar%save = diagvar%output /= output_none
               diagvar%source = object%source
            end select
         case (domain_horizontal, domain_surface, domain_bottom)
            select case (object%source)
            case (source_unknown)   ! Horizontal dependency
            case (source_state)     ! Horizontal state variable
               if (object%presence == presence_internal) then
                  select case (object%domain)
                  case (domain_bottom)
                     nstate_bot = nstate_bot + 1
                     hz_statevar => self%bottom_state_variables(nstate_bot)
                     call object%sms_sum%target%write_indices%append(self%get_bottom_sources_job%arg2_sources(nstate_bot))
                  case (domain_surface)
                     nstate_surf = nstate_surf + 1
                     hz_statevar => self%surface_state_variables(nstate_surf)
                     call object%sms_sum%target%write_indices%append(self%get_surface_sources_job%arg2_sources(nstate_surf))
                  case default
                     hz_statevar => null()
                  end select
                  call copy_variable_metadata(object, hz_statevar)
                  if (associated(object%standard_variables%first)) then
                     select type (standard_variable => object%standard_variables%first%p)
                     class is (type_horizontal_standard_variable)
                        hz_statevar%standard_variable => standard_variable
                     end select
                  end if
                  hz_statevar%initial_value = object%initial_value
               end if
            case default            ! Horizontal diagnostic variable
               ndiag_hz = ndiag_hz + 1
               hz_diagvar => self%horizontal_diagnostic_variables(ndiag_hz)
               call copy_variable_metadata(object, hz_diagvar)
               if (associated(object%standard_variables%first)) then
                  select type (standard_variable => object%standard_variables%first%p)
                  class is (type_horizontal_standard_variable)
                     hz_diagvar%standard_variable => standard_variable
                  end select
               end if
               hz_diagvar%save = hz_diagvar%output /= output_none
               hz_diagvar%source = object%source
            end select
         end select
         link => link%next
      end do

      ! Create lists of variables that may be provided by the host model.
      ! These lists include external dependencies, as well as the model's own state variables,
      ! which may be overridden by the host.
      link => self%root%links%first
      do while (associated(link))
         object =>link%target
         if (.not. object%read_indices%is_empty() .and. &
             .not. (object%presence == presence_external_optional .and. object%source == source_state)) then
            select case (object%domain)
               case (domain_interior);                                  call dependencies%add(link%name)
               case (domain_horizontal, domain_surface, domain_bottom); call dependencies_hz%add(link%name)
               case (domain_scalar);                                    call dependencies_scalar%add(link%name)
            end select

            standard_variable_node => object%standard_variables%first
            do while (associated(standard_variable_node))
               if (standard_variable_node%p%name /= '') then
                  select case (object%domain)
                     case (domain_interior);                                  call dependencies%add(standard_variable_node%p%name)
                     case (domain_horizontal, domain_surface, domain_bottom); call dependencies_hz%add(standard_variable_node%p%name)
                     case (domain_scalar);                                    call dependencies_scalar%add(standard_variable_node%p%name)
                  end select
               end if
               standard_variable_node => standard_variable_node%next
            end do
         end if
         link => link%next
      end do
      call dependencies%to_array(self%dependencies)
      call dependencies_hz%to_array(self%dependencies_hz)
      call dependencies_scalar%to_array(self%dependencies_scalar)
      call dependencies%finalize()
      call dependencies_hz%finalize()
      call dependencies_scalar%finalize()
   contains
      subroutine copy_variable_metadata(internal_variable, external_variable)
         class (type_fabm_variable),    intent(inout)      :: external_variable
         type (type_internal_variable), intent(in), target :: internal_variable

         class (type_base_model), pointer :: owner
         external_variable%name            = get_safe_name(internal_variable%name)
         external_variable%long_name       = internal_variable%long_name
         external_variable%local_long_name = internal_variable%long_name
         external_variable%units           = internal_variable%units
         external_variable%path            = internal_variable%name
         external_variable%minimum         = internal_variable%minimum
         external_variable%maximum         = internal_variable%maximum
         external_variable%missing_value   = internal_variable%missing_value
         external_variable%output          = internal_variable%output
         external_variable%target          => internal_variable

         ! Prepend long names of ancestor models to long name of variable.
         owner => internal_variable%owner
         do while (associated(owner%parent))
            external_variable%long_name = trim(owner%long_name) // ' ' // trim(external_variable%long_name)
            owner => owner%parent
         end do

         call external_variable%properties%update(internal_variable%properties)
      end subroutine
   end subroutine classify_variables

   subroutine require_flux_computation(self, link_list, domain)
      type (type_job),       intent(inout) :: self
      type (type_link_list), intent(in)    :: link_list
      integer,               intent(in)    :: domain

      type (type_link), pointer :: link

      link => link_list%first
      do while (associated(link))
         if (link%target%source == source_state) then
            ! This is a state variable
            select case (domain)
            case (domain_interior)
               if (link%target%domain == domain_interior .and. associated(link%target%sms_sum)) &
                  call self%request_variable(link%target%sms_sum%target)
            case (domain_bottom)
               if (link%target%domain == domain_bottom .and. associated(link%target%sms_sum)) then
                  call self%request_variable(link%target%sms_sum%target)
               elseif (link%target%domain == domain_interior .and. associated(link%target%bottom_flux_sum)) then
                  call self%request_variable(link%target%bottom_flux_sum%target)
               end if
            case (domain_surface)
               if (link%target%domain == domain_surface .and. associated(link%target%sms_sum)) then
                  call self%request_variable(link%target%sms_sum%target)
               elseif (link%target%domain == domain_interior .and. associated(link%target%surface_flux_sum)) then
                  call self%request_variable(link%target%surface_flux_sum%target)
               end if
            case (domain_interior + 999)
               if (link%target%domain == domain_interior .and. associated(link%target%movement_sum)) &
                  call self%request_variable(link%target%movement_sum%target)
            end select
         end if
         link => link%next
      end do
   end subroutine require_flux_computation

   subroutine require_call_all_with_state(self, link_list, domain, source)
      type (type_job),       intent(inout) :: self
      type (type_link_list), intent(in)    :: link_list
      integer,               intent(in)    :: domain
      integer,               intent(in)    :: source

      type (type_link), pointer :: link

      link => link_list%first
      do while (associated(link))
         if (link%target%domain == domain .and. link%original%source == source_state .and. link%target%source == source_state) &
            call self%request_call(link%original%owner, source)
         link => link%next
      end do
   end subroutine require_call_all_with_state

   subroutine create_catalog(self)
      class (type_fabm_model), intent(inout) :: self

      integer                   :: i
      type (type_link), pointer :: link

      ! Add all state variables to the catalog and read cache in the order the host is likely to
      ! have them in memory. This hopefully speeds up access (cache hits).
      do i = 1, size(self%interior_state_variables)
         call self%variable_register%add_to_catalog(self%interior_state_variables(i)%target)
         call self%variable_register%add_to_read_cache(self%interior_state_variables(i)%target)
      end do
      do i = 1, size(self%bottom_state_variables)
         call self%variable_register%add_to_catalog(self%bottom_state_variables(i)%target)
         call self%variable_register%add_to_read_cache(self%bottom_state_variables(i)%target)
      end do
      do i = 1, size(self%surface_state_variables)
         call self%variable_register%add_to_catalog(self%surface_state_variables(i)%target)
         call self%variable_register%add_to_read_cache(self%surface_state_variables(i)%target)
      end do

      ! Add all remaining variables to the catalog
      link => self%links_postcoupling%first
      do while (associated(link))
         call self%variable_register%add_to_catalog(link%target)
         link => link%next
      end do

      allocate(self%catalog%interior(self%variable_register%catalog%interior%count))
      allocate(self%catalog%horizontal(self%variable_register%catalog%horizontal%count))
      allocate(self%catalog%scalar(self%variable_register%catalog%scalar%count))

      ! Allocate and initialize arrays that store the source (host, fabm, user) of all data.
      allocate(self%catalog%interior_sources(size(self%catalog%interior)))
      allocate(self%catalog%horizontal_sources(size(self%catalog%horizontal)))
      allocate(self%catalog%scalar_sources(size(self%catalog%scalar)))
      self%catalog%interior_sources = data_source_none
      self%catalog%horizontal_sources = data_source_none
      self%catalog%scalar_sources = data_source_none
   end subroutine create_catalog

   function get_cache_fill_values(variable_register) result(cache_fill_values)
      type (type_global_variable_register), intent(in) :: variable_register
      type (type_cache_fill_values)                    :: cache_fill_values

      call collect(variable_register%read_cache%interior,    cache_fill_values%read,             use_missing=.false.)
      call collect(variable_register%read_cache%horizontal,  cache_fill_values%read_hz,          use_missing=.false.)
      call collect(variable_register%read_cache%scalar,      cache_fill_values%read_scalar,      use_missing=.false.)
      call collect(variable_register%write_cache%interior,   cache_fill_values%write,            use_missing=.false.)
      call collect(variable_register%write_cache%horizontal, cache_fill_values%write_hz,         use_missing=.false.)
      call collect(variable_register%write_cache%interior,   cache_fill_values%write_missing,    use_missing=.true.)
      call collect(variable_register%write_cache%horizontal, cache_fill_values%write_hz_missing, use_missing=.true.)
   contains
      subroutine collect(variable_list, values, use_missing)
         type (type_variable_list), intent(in)  :: variable_list
         real(rki), allocatable,    intent(out) :: values(:)
         logical,                   intent(in)  :: use_missing

         integer                            :: i
         type (type_variable_node), pointer :: variable_node

         allocate(values(variable_list%count))
         variable_node => variable_list%first
         do i = 1, size(values)
            if (use_missing) then
               values(i) = variable_node%target%missing_value
            else
               values(i) = variable_node%target%prefill_value
            end if
            variable_node => variable_node%next
         end do
      end subroutine
   end function

   subroutine create_store(self)
      class (type_fabm_model), intent(inout) :: self

      type (type_variable_node), pointer :: variable_node

      ! Allocate memory for persistent store
#if _FABM_DIMENSION_COUNT_==0
      call allocate_store()
#elif _FABM_DIMENSION_COUNT_==1
      call allocate_store(self%domain%shape(1))
#elif _FABM_DIMENSION_COUNT_==2
      call allocate_store(self%domain%shape(1), self%domain%shape(2))
#else
      call allocate_store(self%domain%shape(1), self%domain%shape(2), self%domain%shape(3))
#endif

      ! Collect missing values in array for faster access. These will be used to fill masked parts of outputs.
      call collect_fill_values(self%variable_register%store%interior,   self%store%interior_fill_value,      use_missing=.false.)
      call collect_fill_values(self%variable_register%store%horizontal, self%store%horizontal_fill_value,    use_missing=.false.)
      call collect_fill_values(self%variable_register%store%interior,   self%store%interior_missing_value,   use_missing=.true.)
      call collect_fill_values(self%variable_register%store%horizontal, self%store%horizontal_missing_value, use_missing=.true.)

      call reset_store(self)

      ! Register data fields from persistent store in catalog.
      variable_node => self%variable_register%catalog%interior%first
      do while (associated(variable_node))
         if (variable_node%target%store_index > 0) &
            call self%link_interior_data(variable_node%target, self%store%interior(_PREARG_LOCATION_DIMENSIONS_ variable_node%target%store_index), source=data_source_fabm)
         variable_node => variable_node%next
      end do
      variable_node => self%variable_register%catalog%horizontal%first
      do while (associated(variable_node))
         if (variable_node%target%store_index > 0) &
            call self%link_horizontal_data(variable_node%target, self%store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ variable_node%target%store_index), source=data_source_fabm)
         variable_node => variable_node%next
      end do

   contains

      subroutine allocate_store(_LOCATION_)
         _DECLARE_ARGUMENTS_LOCATION_
         allocate(self%store%interior(_PREARG_LOCATION_ 0:self%variable_register%store%interior%count))
         allocate(self%store%horizontal(_PREARG_HORIZONTAL_LOCATION_ 0:self%variable_register%store%horizontal%count))
      end subroutine

      subroutine collect_fill_values(variable_list, values, use_missing)
         type (type_variable_list), intent(in)  :: variable_list
         real(rke), allocatable,    intent(out) :: values(:)
         logical,                   intent(in)  :: use_missing

         integer                            :: i
         type (type_variable_node), pointer :: variable_node

         allocate(values(variable_list%count))
         variable_node => variable_list%first
         do i = 1, size(values)
            if (use_missing) then
               values(i) = variable_node%target%missing_value
            else
               values(i) = variable_node%target%prefill_value
            end if
            variable_node => variable_node%next
         end do
      end subroutine

   end subroutine create_store

   subroutine reset_store(self)
      class (type_fabm_model), intent(inout) :: self

      integer :: i

      ! Initialize persistent store entries to fill value.
      ! For constant outputs, their values will be set here, and never touched again.
      do i = 1, self%variable_register%store%interior%count
         self%store%interior(_PREARG_LOCATION_DIMENSIONS_ i) = self%store%interior_fill_value(i)
      end do
      do i = 1, self%variable_register%store%horizontal%count
         self%store%horizontal(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ i) = self%store%horizontal_fill_value(i)
      end do
   end subroutine

   recursive subroutine merge_indices(model, log_unit)
      class (type_base_model), intent(inout)        :: model
      integer,                 intent(in), optional :: log_unit

      type (type_model_list_node), pointer :: child

      select type (model)
      class is (type_reduction_operator)
         call model%merge_components(log_unit)
      end select

      ! Process children
      child => model%children%first
      do while (associated(child))
         call merge_indices(child%model, log_unit)
         child => child%next
      end do
   end subroutine merge_indices

   subroutine filter_expressions(self)
      class (type_fabm_model),intent(inout) :: self

      class (type_expression),             pointer :: current, previous, next
      class (type_depth_integral),         pointer :: integral
      class (type_bounded_depth_integral), pointer :: bounded_integral
      logical                                      :: filter

      previous => null()
      current => self%root%first_expression
      do while (associated(current))
         filter = .false.
         select type (current)
         class is (type_vertical_integral)
            if (current%minimum_depth == 0._rki .and. current%maximum_depth == huge(current%maximum_depth)) then
               allocate(integral)
            else
               allocate(bounded_integral)
               bounded_integral%minimum_depth = current%minimum_depth
               bounded_integral%maximum_depth = current%maximum_depth
               integral => bounded_integral
            end if
            integral%average = current%average
            call self%root%add_child(integral, trim(current%output_name) // '_calculator', configunit=-1)
            call integral%request_coupling(integral%id_input, current%input_name)
            call self%root%request_coupling(current%output_name, integral%id_output%link%target%name)
            filter = .true.
         end select

         ! If FABM handles this expression internally, remove it from the list.
         next => current%next
         if (filter) then
            if (associated(previous)) then
               previous%next => next
            else
               self%root%first_expression => next
            end if
            deallocate(current)
         else
            previous => current
         end if
         current => next
      end do
   end subroutine filter_expressions

   function get_variable_by_name(self, name, domain) result(variable)
      class (type_fabm_model), intent(in)    :: self
      character(len=*),        intent(in)    :: name
      integer,                 intent(in)    :: domain
      type (type_internal_variable), pointer :: variable

      type (type_link), pointer :: link

      variable => null()

      link => self%root%links%first
      do while (associated(link))
         if (iand(link%target%domain, domain) /= 0) then
            if (link%name == name .or. get_safe_name(link%name) == name) then
               variable => link%target
               return
            end if
         end if
         link => link%next
      end do

      ! Name not found among variable names. Now try standard names that are in use.
      link => self%root%links%first
      do while (associated(link))
         if (iand(link%target%domain, domain) /= 0 .and. link%target%standard_variables%contains(name)) then
            variable => link%target
            return
         end if
         link => link%next
      end do
   end function get_variable_by_name

   function get_variable_by_standard_variable(self, standard_variable) result(variable)
      class (type_fabm_model), intent(in)         :: self
      class (type_base_standard_variable), target :: standard_variable
      type (type_internal_variable), pointer :: variable

      type (type_link), pointer :: link

      variable => null()
      link => self%root%links%first
      do while (associated(link))
         if (link%target%standard_variables%contains(standard_variable)) then
            variable => link%target
            return
         end if
         link => link%next
      end do
      if (standard_variable%aggregate_variable) then
         select type (standard_variable)
         class is (type_interior_standard_variable)
            variable => self%root%find_object('zero')
         class is (type_horizontal_standard_variable)
            variable => self%root%find_object('zero_hz')
         end select
      end if
   end function get_variable_by_standard_variable

end module fabm

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS (GNU Public License - www.gnu.org)
!-----------------------------------------------------------------------

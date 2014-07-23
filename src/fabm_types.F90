#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_types --- Derived types and procedures for use by biogeochemical modules
!
! !INTERFACE:
   module fabm_types
!
! !DESCRIPTION:
! This module contains the derived types and procedures that are used for communication between
! biogeochemical models and FABM. This module provides types for storing model data (e.g.,
! metadata for variables and parameters), and logic for registration of model objects
! (state and diagnostic variables), retrieval of model settings (parameter values) and coupling.
!
! !USES:
   use fabm_standard_variables
   use fabm_properties
   use fabm_driver
!
   implicit none
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!
   ! Base data type for biogeochemical models.
   public type_base_model

   ! Collection with standard variables (e.g., temperature, practical_salinity)
   public standard_variables
   public type_standard_variable, type_bulk_standard_variable

   ! Variable identifier types used by biogeochemical models
   public type_variable_id
   public type_diagnostic_variable_id
   public type_horizontal_diagnostic_variable_id
   public type_state_variable_id
   public type_surface_state_variable_id
   public type_bottom_state_variable_id
   public type_dependency_id
   public type_horizontal_dependency_id
   public type_global_dependency_id
   public type_conserved_quantity_id

   ! Data types and procedures for variable management - used by FABM internally only.
   public type_link
   public type_internal_object,type_internal_variable,type_bulk_variable,type_horizontal_variable,type_scalar_variable
   public type_environment
   public freeze_model_info
   public find_dependencies

   public type_model_list,type_model_list_node

   public get_free_unit
   public get_safe_name

   public type_expression, type_bulk_expression, type_horizontal_expression

   public type_weighted_sum,type_simple_depth_integral
   public connect_bulk_state_variable_id

   public type_bulk_data_pointer,type_horizontal_data_pointer,type_scalar_data_pointer
   public type_bulk_data_pointer_pointer,type_horizontal_data_pointer_pointer,type_scalar_data_pointer_pointer
   public type_aggregate_variable, type_contributing_variable, get_aggregate_variable
   public time_treatment2output,output2time_treatment
   public append_data_pointer

! !PUBLIC DATA MEMBERS:
!
   integer, parameter, public :: attribute_length = 256

   integer, parameter, public :: rk = _FABM_REAL_KIND_

   integer, parameter, public :: domain_bulk = 0, domain_bottom = 1, domain_surface = 2

   integer, parameter, public :: presence_internal = 1, presence_external_required = 2, presence_external_optional = 6
!
! !PUBLIC TYPES:
!
   integer, parameter, public :: output_none                 = 0, &
                                 output_instantaneous        = 1, &
                                 output_time_integrated      = 2, &
                                 output_time_step_averaged   = 4, &
                                 output_time_step_integrated = 8

   ! For pre 2014-01 backward compatibility only (please use output_* instead)
   integer, parameter, public :: time_treatment_last            = 0, &
                                 time_treatment_integrated      = 1, &
                                 time_treatment_averaged        = 2, &
                                 time_treatment_step_integrated = 3

   ! ====================================================================================================
   ! Data types for pointers to variable values.
   ! ====================================================================================================

   type type_bulk_data_pointer
      real(rk),pointer _DIMENSION_GLOBAL_ :: p => null()
   end type

   type type_horizontal_data_pointer
      real(rk),pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: p => null()
   end type

   type type_scalar_data_pointer
      real(rk),pointer :: p => null()
   end type

   type type_integer_pointer
      integer,pointer :: p => null()
   end type

   type type_real_pointer
      real(rk),pointer :: p => null()
   end type

   type type_real_pointer_set
      type (type_real_pointer),allocatable :: pointers(:)
   contains
      procedure :: append    => real_pointer_set_append
      procedure :: extend    => real_pointer_set_extend
      procedure :: set_value => real_pointer_set_set_value
   end type

   type type_integer_pointer_set
      type (type_integer_pointer),allocatable :: pointers(:)
   contains
      procedure :: append    => integer_pointer_set_append
      procedure :: extend    => integer_pointer_set_extend
      procedure :: set_value => integer_pointer_set_set_value
      procedure :: is_empty  => integer_pointer_set_is_empty
   end type

   ! ====================================================================================================
   ! Data types to hold pointers to (components of) variable identifiers used by biogeochemical models.
   ! ====================================================================================================

   type type_bulk_data_pointer_pointer
      type (type_bulk_data_pointer),pointer :: p => null()
   end type

   type type_horizontal_data_pointer_pointer
      type (type_horizontal_data_pointer),pointer :: p => null()
   end type

   type type_scalar_data_pointer_pointer
      type (type_scalar_data_pointer),pointer :: p => null()
   end type

   ! ====================================================================================================
   ! Variable identifiers used by biogeochemical models.
   ! ====================================================================================================

   type,abstract :: type_variable_id
      type (type_link), pointer :: link => null()
   end type

   type,extends(type_variable_id) :: type_state_variable_id
      integer                       :: state_index = -1
      type (type_bulk_data_pointer) :: data
      real(rk)                      :: background = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_bottom_state_variable_id
      integer                             :: bottom_state_index = -1
      type (type_horizontal_data_pointer) :: horizontal_data
      real(rk)                            :: background = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_surface_state_variable_id
      integer                             :: surface_state_index = -1
      type (type_horizontal_data_pointer) :: horizontal_data
      real(rk)                            :: background = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_diagnostic_variable_id
      integer :: diag_index = -1
   end type

   type,extends(type_variable_id) :: type_horizontal_diagnostic_variable_id
      integer :: horizontal_diag_index = -1
   end type

   type,extends(type_variable_id) :: type_dependency_id
      type (type_bulk_data_pointer) :: data
      real(rk)                      :: background = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_horizontal_dependency_id
      type (type_horizontal_data_pointer) :: horizontal_data
      real(rk)                            :: background = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_global_dependency_id
      type (type_scalar_data_pointer) :: global_data
      real(rk)                        :: background = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_conserved_quantity_id
      integer :: cons_index = -1
   end type

   ! ====================================================================================================
   ! Derived types used internally to register contributions of variables to aggregate variables.
   ! ====================================================================================================

   type type_contribution
      type (type_bulk_standard_variable)  :: target
      real(rk)                            :: scale_factor = 1.0_rk
      logical                             :: include_background = .false.
      type (type_contribution),pointer    :: next => null()
   end type

   type type_contribution_list
      type (type_contribution),pointer :: first => null()
   contains
      procedure :: add => contribution_list_add
   end type

   type type_contributing_variable
      type (type_link),                 pointer :: link => null()
      real(rk)                                  :: scale_factor = 1.0_rk
      logical                                   :: include_background = .false.
      type (type_contributing_variable),pointer :: next => null()
   end type

   type type_aggregate_variable
      type (type_bulk_standard_variable)            :: standard_variable
      type (type_contributing_variable),   pointer  :: first_contributing_variable => null()
      class (type_weighted_sum),           pointer  :: sum => null()
      class (type_horizontal_weighted_sum),pointer  :: horizontal_sum => null()
      type (type_aggregate_variable),      pointer  :: next => null()
      logical                                       :: bulk_required = .false.
      logical                                       :: horizontal_required = .false.

      type (type_diagnostic_variable_id)            :: id_rate
      type (type_horizontal_diagnostic_variable_id) :: id_horizontal_rate
      integer,allocatable                           :: state_indices(:)
      real(rk),allocatable                          :: state_scale_factors(:)
   end type

   ! ====================================================================================================
   ! Derived types used internally to store information on model variables and model references.
   ! ====================================================================================================

   type,abstract :: type_internal_object
      character(len=attribute_length) :: name      = ''
      character(len=attribute_length) :: long_name = ''
      type (type_property_dictionary) :: properties
   end type

   type,extends(type_internal_object),abstract :: type_internal_variable
      character(len=attribute_length) :: units          = ''
      real(rk)                        :: minimum        = -1.e20_rk
      real(rk)                        :: maximum        =  1.e20_rk
      real(rk)                        :: missing_value  = -2.e20_rk
      real(rk)                        :: initial_value  = 0.0_rk
      integer                         :: output         = output_instantaneous
      integer                         :: presence       = presence_internal
      class (type_base_model),pointer :: source_model   => null()
      type (type_contribution_list)   :: contributions

      type (type_integer_pointer_set) :: state_indices,write_indices
      type (type_real_pointer_set)    :: background_values
   end type

   type,extends(type_internal_variable) :: type_bulk_variable
      ! Metadata
      real(rk)                                      :: vertical_movement         = 0.0_rk
      logical                                       :: no_precipitation_dilution = .false.
      logical                                       :: no_river_dilution         = .false.
      type (type_bulk_standard_variable)            :: standard_variable

      ! Arrays with all associated data and index pointers.
      type (type_bulk_data_pointer_pointer),allocatable :: alldata(:) 
   end type type_bulk_variable

   type,extends(type_internal_variable) :: type_horizontal_variable
      ! Metadata
      integer                                  :: domain         = domain_bottom
      type (type_horizontal_standard_variable) :: standard_variable
      
      ! Arrays with all associated data and index pointers.
      type (type_horizontal_data_pointer_pointer),allocatable :: alldata(:) 
   end type type_horizontal_variable

   type,extends(type_internal_variable) :: type_scalar_variable
      ! Metadata
      type (type_global_standard_variable) :: standard_variable

      ! Arrays with all associated data and index pointers.
      type (type_scalar_data_pointer_pointer),allocatable :: alldata(:) 
   end type type_scalar_variable

   type type_link
      character(len=attribute_length)       :: name    = ''
      class (type_internal_object), pointer :: target  => null()
      logical                               :: owner   = .true.
      type (type_link), pointer             :: next    => null()
   end type

   ! ====================================================================================================
   ! Types that describe custom expressions (arbitrary functions of one or more variables).
   ! ====================================================================================================

   type,abstract :: type_expression
      class (type_expression), pointer :: next        => null()
      character(len=attribute_length)  :: output_name = ''
   end type

   type,abstract,extends(type_expression) :: type_bulk_expression
      type (type_bulk_data_pointer),pointer :: out => null()
   end type

   type,abstract,extends(type_expression) :: type_horizontal_expression
      type (type_horizontal_data_pointer),pointer :: out => null()
   end type

   ! ====================================================================================================
   ! Base model type, used by biogeochemical models to inherit from, and by external host to
   ! get variable lists and metadata.
   ! ====================================================================================================

   type type_model_list_node
      class (type_base_model),    pointer :: model => null()
      type (type_model_list_node),pointer :: next  => null()
   end type

   type type_model_list
      type (type_model_list_node),pointer :: first => null()
   contains
      procedure :: append   => model_list_append
      procedure :: extend   => model_list_extend
      procedure :: find     => model_list_find
      procedure :: count    => model_list_count
      procedure :: finalize => model_list_finalize
   end type

   type type_base_model
      ! Flag determining whether the contents of the type are "frozen", i.e., they will not change anymore.
      logical :: frozen = .false.

      ! Pointers to linked models in the model tree.
      class (type_base_model),pointer :: parent => null()
      type (type_model_list)          :: children

      ! Model name and variable prefixes.
      character(len=64) :: name      = ''
      character(len=64) :: long_name = ''
      character(len=64) :: type_name = ''

      ! Models constituents: links to variables, coupling requests, parameters, expressions
      type (type_link),              pointer :: first_link               => null()
      type (type_aggregate_variable),pointer :: first_aggregate_variable => null()

      type (type_hierarchical_dictionary) :: couplings
      type (type_hierarchical_dictionary) :: parameters

      class (type_expression), pointer :: first_expression => null()

      real(rk) :: dt = 1.0_rk

      logical :: check_conservation = .false.

      ! Identifiers for variables that contain zero at all time.
      type (type_dependency_id)            :: id_zero
      type (type_horizontal_dependency_id) :: id_zero_hz
   contains

      ! Procedures that can be used to add child models during initialization.
      procedure :: add_child

      procedure :: request_coupling_for_link
      procedure :: request_coupling_for_name
      procedure :: request_coupling_for_id
      generic   :: request_coupling => request_coupling_for_link,request_coupling_for_name,request_coupling_for_id

      ! Procedures that may be used to register model variables and dependencies during initialization.
      procedure :: register_bulk_state_variable
      procedure :: register_bottom_state_variable
      procedure :: register_surface_state_variable

      procedure :: register_bulk_diagnostic_variable
      procedure :: register_horizontal_diagnostic_variable

      procedure :: register_named_bulk_dependency
      procedure :: register_standard_bulk_dependency
      procedure :: register_named_horizontal_dependency
      procedure :: register_standard_horizontal_dependency
      procedure :: register_named_global_dependency
      procedure :: register_standard_global_dependency

      generic :: register_bulk_dependency       => register_named_bulk_dependency, register_standard_bulk_dependency
      generic :: register_horizontal_dependency => register_named_horizontal_dependency, register_standard_horizontal_dependency
      generic :: register_global_dependency     => register_named_global_dependency, register_standard_global_dependency

      procedure :: register_standard_conserved_quantity
      procedure :: register_custom_conserved_quantity
      procedure :: register_bulk_state_dependency_ex
      procedure :: register_bottom_state_dependency_ex
      procedure :: register_surface_state_dependency_ex
      procedure :: register_bulk_state_dependency_old
      procedure :: register_bottom_state_dependency_old
      procedure :: register_surface_state_dependency_old

      procedure :: add_bulk_variable
      procedure :: add_horizontal_variable
      procedure :: add_scalar_variable
      procedure :: add_object
      procedure :: add_alias

      procedure :: find_link
      procedure :: find_object
      procedure :: find_model

      generic :: register_bulk_state_dependency => register_bulk_state_dependency_ex, register_bulk_state_dependency_old
      generic :: register_bottom_state_dependency => register_bottom_state_dependency_ex, register_bottom_state_dependency_old
      generic :: register_surface_state_dependency => register_surface_state_dependency_ex, register_surface_state_dependency_old

      procedure :: set_variable_property_real
      procedure :: set_variable_property_integer
      procedure :: set_variable_property_logical
      generic   :: set_variable_property => set_variable_property_real,set_variable_property_integer,set_variable_property_logical

      procedure :: add_to_aggregate_variable

      procedure :: register_bulk_expression_dependency
      procedure :: register_horizontal_expression_dependency
      generic :: register_expression_dependency => register_bulk_expression_dependency,register_horizontal_expression_dependency

      generic :: register_state_variable      => register_bulk_state_variable,register_bottom_state_variable, &
                                                 register_surface_state_variable
      generic :: register_diagnostic_variable => register_bulk_diagnostic_variable,register_horizontal_diagnostic_variable
      generic :: register_dependency          => register_named_bulk_dependency, register_standard_bulk_dependency, &
                                                 register_named_horizontal_dependency, register_standard_horizontal_dependency, &
                                                 register_named_global_dependency, register_standard_global_dependency, &
                                                 register_bulk_expression_dependency, register_horizontal_expression_dependency
      generic :: register_state_dependency    => register_bulk_state_dependency_ex,register_bottom_state_dependency_ex, &
                                                 register_surface_state_dependency_ex,register_bulk_state_dependency_old, &
                                                 register_bottom_state_dependency_old,register_surface_state_dependency_old
      generic :: register_conserved_quantity  => register_standard_conserved_quantity, register_custom_conserved_quantity

      ! Procedures that may be used to query parameter values during initialization.
      procedure :: get_real_parameter
      procedure :: get_integer_parameter
      procedure :: get_logical_parameter
      procedure :: get_string_parameter
      generic :: get_parameter => get_real_parameter,get_integer_parameter,get_logical_parameter,get_string_parameter

      ! ----------------------------------------------------------------------------------------------------
      ! Procedures below may be overridden by biogeochemical models to provide custom data or functionality.
      ! ----------------------------------------------------------------------------------------------------

      ! Model initialization.
      procedure :: initialize               => base_initialize
      procedure :: initialize_state         => base_initialize_state
      procedure :: initialize_surface_state => base_initialize_horizontal_state
      procedure :: initialize_bottom_state  => base_initialize_horizontal_state

      ! Providing process rates and diagnostics in pelagic, at surface, and at bottom.
      procedure :: do                       => base_do
      procedure :: do_check_conservation    => base_do_check_conservation
      procedure :: do_bottom                => base_do_bottom
      procedure :: do_surface               => base_do_surface
      procedure :: do_ppdd                  => base_do_ppdd
      procedure :: do_bottom_ppdd           => base_do_bottom_ppdd

      ! Advanced functionality: variable vertical movement and light attenuation, feedbacks to drag and albedo.
      procedure :: get_vertical_movement    => base_get_vertical_movement
      procedure :: get_light_extinction     => base_get_light_extinction
      procedure :: get_drag                 => base_get_drag
      procedure :: get_albedo               => base_get_albedo
      procedure :: get_light                => base_get_light

      ! Bookkeeping: calculate total of conserved quantities, check and repair model state.
      procedure :: get_conserved_quantities => base_get_conserved_quantities
      procedure :: state_to_conserved_quantities => base_state_to_conserved_quantities
      procedure :: get_horizontal_conserved_quantities => base_get_horizontal_conserved_quantities
      procedure :: check_state              => base_check_state
      procedure :: check_surface_state      => base_check_surface_state
      procedure :: check_bottom_state       => base_check_bottom_state
      procedure :: fatal_error              => base_fatal_error
      procedure :: log_message              => base_log_message

      ! For backward compatibility only - do not use these in new models!
      procedure :: set_domain               => base_set_domain
      procedure :: do_benthos               => base_do_benthos           ! superceded by do_bottom
      procedure :: do_benthos_ppdd          => base_do_benthos_ppdd      ! superceded by do_bottom_ppdd
      procedure :: get_surface_exchange     => base_get_surface_exchange ! superceded by do_surface

      ! Hooks called by FABM - usable by inheriting models
      procedure :: before_coupling => base_before_coupling
      procedure :: after_coupling  => base_after_coupling
   end type type_base_model

   ! ====================================================================================================
   ! Derived type for holding global data needed by biogeochemical model tree.
   ! ====================================================================================================

   type type_environment
      ! Declare the arrays for diagnostic variable values.
      real(rk),allocatable _DIMENSION_GLOBAL_PLUS_1_            :: diag
      real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: diag_hz
      integer                                                   :: nstate

      real(rk),allocatable _DIMENSION_GLOBAL_PLUS_1_            :: rhs

      real(rk),allocatable _DIMENSION_GLOBAL_                   :: zero
      real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_        :: zero_hz

#ifdef _FABM_MASK_
      _FABM_MASK_TYPE_,pointer _DIMENSION_GLOBAL_ :: mask => null()
#endif
   end type type_environment

   ! ====================================================================================================
   ! Base type for a model object factory (generates a model object from a model name)
   ! An implementation of this type is provided in fabm_library.F90.
   ! Institutes or groups can create inherit from this type to create their own model factories,
   ! which then need to be added to the root factory in fabm_library.F90.
   ! This makes it possible to introduce a large number of new models with only two lines added
   ! in the FABM core.
   ! ====================================================================================================

   type,public :: type_base_model_factory
      class (type_base_model_factory),pointer :: first_child => null()
      class (type_base_model_factory),pointer :: next_sibling => null()
   contains
      procedure :: create => abstract_model_factory_create
      procedure :: add    => abstract_model_factory_add
   end type

   class (type_base_model_factory),pointer,save,public :: factory => null()

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

   type,extends(type_base_model) :: type_weighted_sum
      character(len=attribute_length) :: output_name      = ''
      character(len=attribute_length) :: output_long_name = ''
      character(len=attribute_length) :: output_units     = ''
      real(rk)                        :: offset           = 0.0_rk
      type (type_diagnostic_variable_id) :: id_output
      type (type_component),pointer   :: first => null()
   contains
      procedure :: initialize     => weighted_sum_initialize
      procedure :: add_component  => weighted_sum_add_component
      procedure :: evaluate       => weighted_sum_evaluate
      procedure :: do             => weighted_sum_do
      procedure :: after_coupling => weighted_sum_after_coupling
   end type

   type,extends(type_base_model) :: type_horizontal_weighted_sum
      character(len=attribute_length) :: output_name      = ''
      character(len=attribute_length) :: output_long_name = ''
      character(len=attribute_length) :: output_units     = ''
      real(rk)                        :: offset           = 0.0_rk
      type (type_horizontal_diagnostic_variable_id) :: id_output
      type (type_horizontal_component),pointer   :: first => null()
   contains
      procedure :: add_component       => horizontal_weighted_sum_add_component
      procedure :: initialize          => horizontal_weighted_sum_initialize
      procedure :: evaluate_horizontal => horizontal_weighted_sum_evaluate_horizontal
      procedure :: do_bottom           => horizontal_weighted_sum_do_bottom
      procedure :: after_coupling      => horizontal_weighted_sum_after_coupling
   end type

   type,extends(type_base_model) :: type_simple_depth_integral
      type (type_dependency_id)                     :: id_input
      type (type_horizontal_dependency_id)          :: id_depth
      type (type_horizontal_diagnostic_variable_id) :: id_output
      real(rk)                                      :: minimum_depth = 0.0_rk
      real(rk)                                      :: maximum_depth = huge(1.0_rk)
      logical                                       :: average       = .false.
   contains
      procedure :: initialize => simple_depth_integral_initialize
      procedure :: do         => simple_depth_integral_do
   end type

   ! ====================================================================================================
   ! Interfaces
   ! ====================================================================================================

   interface append_data_pointer
      module procedure append_bulk_data_pointer
      module procedure append_horizontal_data_pointer
      module procedure append_scalar_data_pointer
   end interface

!-----------------------------------------------------------------------

   contains

   subroutine base_initialize(self,configunit)
      class (type_base_model),intent(inout),target :: self
      integer,                intent(in)           :: configunit
      call self%fatal_error('base_initialize','Model must implement the "initialize" subroutine.')
   end subroutine

   subroutine base_initialize_state(self,_ARGUMENTS_INITIALIZE_STATE_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_INITIALIZE_STATE_
   end subroutine

   subroutine base_initialize_horizontal_state(self,_ARGUMENTS_INITIALIZE_HORIZONTAL_STATE_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_INITIALIZE_HORIZONTAL_STATE_
   end subroutine

   ! Providing process rates and diagnostics
   subroutine base_do(self,_ARGUMENTS_DO_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
   end subroutine

   subroutine base_do_check_conservation(self,_ARGUMENTS_DO_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) _DIMENSION_SLICE_AUTOMATIC_ :: current_sum
      real(rk) _DIMENSION_SLICE_PLUS_1_ALLOCATABLE_,allocatable :: dy
      integer :: index,i
      real(rk) :: scale_factor
      type (type_aggregate_variable),pointer :: aggregate_variable

      ! Fist compute sink and source terms and add them to target array (rhs).
      ! NB it is relatively fast to allocate and zero the array every time, because
      ! most compilers then create a zero-on-write array that is efficient for sparse access
      ! (access here is sparse - individual models will tend to write to subsets of dy only).
      ! This appears to be more efficient than prior allocation at startup, and explicitly
      ! zeroing every time.
      allocate(dy(_SLICE_SHAPE_ environment%nstate))
      dy = 0.0_rk
      call self%do(_ARGUMENTS_ND_,dy)
      rhs = rhs + dy

      ! Now compute sink/source contributions for the different aggregate quantities.
      aggregate_variable => self%first_aggregate_variable
      do while (associated(aggregate_variable))
         if (_VARIABLE_REGISTERED_(aggregate_variable%id_rate)) then
            ! First set sum to zero across spatial domain
            current_sum = 0.0_rk

            ! Add contributions of all state variables
            do i=1,size(aggregate_variable%state_indices)
               index = aggregate_variable%state_indices(i)
               scale_factor = aggregate_variable%state_scale_factors(i)
               _LOOP_BEGIN_
                  current_sum _INDEX_SLICE_ = current_sum _INDEX_SLICE_ + scale_factor*dy _INDEX_SLICE_PLUS_1_(index)
               _LOOP_END_
            end do

            ! Store sum across spatial domain
            _LOOP_BEGIN_
               _SET_DIAGNOSTIC_(aggregate_variable%id_rate,current_sum _INDEX_SLICE_)
            _LOOP_END_
         end if
         aggregate_variable => aggregate_variable%next
      end do
   end subroutine

   subroutine base_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_PPDD_
   end subroutine

   subroutine base_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      call self%do_benthos(_ARGUMENTS_DO_BOTTOM_)
   end subroutine

   subroutine base_do_bottom_ppdd(self,_ARGUMENTS_DO_BOTTOM_PPDD_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_
      call self%do_benthos_ppdd(_ARGUMENTS_DO_BOTTOM_PPDD_)
   end subroutine

   subroutine base_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_
      call self%get_surface_exchange(_ARGUMENTS_DO_SURFACE_)
   end subroutine

   ! Vertical movement, light attenuation, feedbacks to drag and albedo
   subroutine base_get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
   end subroutine

   subroutine base_get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_EXTINCTION_
   end subroutine

   subroutine base_get_drag(self,_ARGUMENTS_GET_DRAG_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_DRAG_
   end subroutine

   subroutine base_get_albedo(self,_ARGUMENTS_GET_ALBEDO_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_ALBEDO_
   end subroutine

   subroutine base_get_light(self,_ARGUMENTS_VERT_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_VERT_
   end subroutine

   ! Bookkeeping
   subroutine base_get_conserved_quantities(self,_ARGUMENTS_GET_CONSERVED_QUANTITIES_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_CONSERVED_QUANTITIES_
   end subroutine

   subroutine base_get_horizontal_conserved_quantities(self,_ARGUMENTS_GET_HORIZONTAL_CONSERVED_QUANTITIES_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_HORIZONTAL_CONSERVED_QUANTITIES_
   end subroutine

   subroutine base_check_state(self,_ARGUMENTS_CHECK_STATE_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_CHECK_STATE_
   end subroutine

   subroutine base_check_surface_state(self,_ARGUMENTS_CHECK_SURFACE_STATE_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_CHECK_SURFACE_STATE_
   end subroutine

   subroutine base_check_bottom_state(self,_ARGUMENTS_CHECK_BOTTOM_STATE_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_CHECK_BOTTOM_STATE_
   end subroutine

   subroutine base_fatal_error(self,location,message)
      class (type_base_model), intent(in) :: self
      character(len=*),        intent(in) :: location,message
      if (self%name/='') then
         call fatal_error('model "'//trim(self%name)//'", '//trim(location),message)
      else
         call fatal_error(location,message)
      end if
   end subroutine

   subroutine base_log_message(self,message)
      class (type_base_model), intent(in) :: self
      character(len=*),        intent(in) :: message
      if (self%name/='') then
         call log_message('model "'//trim(self%name)//'": '//message)
      else
         call log_message(message)
      end if
   end subroutine

   ! For backward compatibility only
   subroutine base_set_domain(self _ARG_LOCATION_)
      class (type_base_model),intent(inout) :: self
      _DECLARE_LOCATION_ARG_
   end subroutine

   subroutine base_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
   end subroutine

   subroutine base_do_benthos_ppdd(self,_ARGUMENTS_DO_BOTTOM_PPDD_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_
   end subroutine

   subroutine base_get_surface_exchange(self,_ARGUMENTS_DO_SURFACE_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_
   end subroutine

   subroutine base_before_coupling(self)
      class (type_base_model), intent(inout) :: self
   end subroutine

   subroutine base_after_coupling(self)
      class (type_base_model), intent(inout) :: self
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes model information.
!
! !INTERFACE:
   recursive subroutine add_child(self,model,name,long_name,configunit)
!
! !DESCRIPTION:
!  This function initializes the members of a model information derived type,
!  by setting them to a reasonable default value.
!
! !INPUT/OUTPUT PARAMETER:
      class (type_base_model),target,intent(inout) :: self,model
      character(len=*),              intent(in)    :: name
      character(len=*),optional,     intent(in)    :: long_name
      integer,                       intent(in)    :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
      if (associated(model%parent)) call self%fatal_error('add_child','The provided child model has already been connected to a parent.')

      model%name = name
      if (present(long_name)) then
         model%long_name = trim(long_name)
      else
         model%long_name = trim(name)
      end if
      model%parent => self
      call self%parameters%add_child(model%parameters,name)
      call self%couplings%add_child(model%couplings,name)
      call self%children%append(model)
      call model%initialize(configunit)
   end subroutine add_child
!EOC

   subroutine set_variable_property_real(self,variable,name,value)
      class (type_base_model), intent(inout) :: self
      class (type_variable_id),intent(inout) :: variable
      character(len=*),        intent(in)    :: name
      real(rk),                intent(in)    :: value
      if (.not.associated(variable%link)) call self%fatal_error('set_variable_property_real','variable has not been registered')
      call variable%link%target%properties%set_real(name,value)
   end subroutine

   subroutine set_variable_property_integer(self,variable,name,value)
      class (type_base_model), intent(inout) :: self
      class (type_variable_id),intent(inout) :: variable
      character(len=*),        intent(in)    :: name
      integer,                 intent(in)    :: value
      if (.not.associated(variable%link)) call self%fatal_error('set_variable_property_integer','variable has not been registered')
      call variable%link%target%properties%set_integer(name,value)
   end subroutine

   subroutine set_variable_property_logical(self,variable,name,value)
      class (type_base_model), intent(inout) :: self
      class (type_variable_id),intent(inout) :: variable
      character(len=*),        intent(in)    :: name
      logical,                 intent(in)    :: value
      if (.not.associated(variable%link)) call self%fatal_error('set_variable_property_logical','variable has not been registered')
      call variable%link%target%properties%set_logical(name,value)
   end subroutine

   subroutine add_to_aggregate_variable(self,target,variable_id,scale_factor,include_background)
      class (type_base_model),           intent(inout) :: self
      type (type_bulk_standard_variable),intent(in)    :: target
      class (type_variable_id),          intent(inout) :: variable_id
      real(rk),optional,                 intent(in)    :: scale_factor
      logical,optional,                  intent(in)    :: include_background

      if (.not.associated(variable_id%link)) &
         call self%fatal_error('add_to_aggregate_variable','variable has not been registered')

      select type (variable=>variable_id%link%target)
         class is (type_internal_variable)
            call variable%contributions%add(target,scale_factor,include_background)
         class default
            call self%fatal_error('add_to_aggregate_variable','Only variables can contribute to conserved quantities at present.')
      end select

   end subroutine

   subroutine contribution_list_add(self,target,scale_factor,include_background)
      class (type_contribution_list),    intent(inout) :: self
      type (type_bulk_standard_variable),intent(in)    :: target
      real(rk),optional,                 intent(in)    :: scale_factor
      logical,optional,                  intent(in)    :: include_background

      type (type_contribution),pointer :: contribution

      if (.not.associated(self%first)) then
         ! Add first contribution
         allocate(self%first)
         contribution => self%first
      else
         ! First determine whether a contribution of this target variable already exists.
         contribution => self%first
         do while (associated(contribution))
            if (target%compare(contribution%target)) exit
            contribution => contribution%next
         end do
         if (.not.associated(contribution)) then
            ! First find tail of existing linked list.
            contribution => self%first
            do while (associated(contribution%next))
               contribution => contribution%next
            end do

            ! Append contribution to end of list.
            allocate(contribution%next)
            contribution => contribution%next
         end if
      end if

      contribution%target = target
      if (present(scale_factor)) contribution%scale_factor = scale_factor
      if (present(include_background)) contribution%include_background = include_background
   end subroutine

   subroutine model_list_append(self,model)
      class (type_model_list), intent(inout) :: self
      class (type_base_model),target         :: model

      type (type_model_list_node),pointer :: node

      if (.not.associated(self%first)) then
         allocate(self%first)
         node => self%first
      else
         node => self%first
         do while (associated(node%next))
            node => node%next
         end do
         allocate(node%next)
         node => node%next
      end if
      node%model => model
   end subroutine

   subroutine model_list_extend(self,source)
      class (type_model_list), intent(inout) :: self
      class (type_model_list), intent(in)    :: source

      type (type_model_list_node),pointer :: node

      node => source%first
      do while (associated(node))
         call self%append(node%model)
         node => node%next
      end do
   end subroutine

   function model_list_find(self,model) result(node)
      class (type_model_list),       intent(in) :: self
      class (type_base_model),target,intent(in) :: model

      type (type_model_list_node),pointer :: node

      node => self%first
      do while (associated(node))
         if (associated(node%model,model)) return
         node => node%next
      end do
      nullify(node)
   end function

   function model_list_count(self,model) result(count)
      class (type_model_list),       intent(in) :: self
      class (type_base_model),target,intent(in) :: model

      integer :: count

      type (type_model_list_node),pointer :: node

      count = 0
      node => self%first
      do while (associated(node))
         if (associated(node%model,model)) count = count + 1
         node => node%next
      end do
   end function

   subroutine model_list_finalize(self)
      class (type_model_list),       intent(in) :: self

      type (type_model_list_node),pointer :: node,next

      node => self%first
      do while (associated(node))
         next => node%next
         deallocate(node)
         node => next
      end do
   end subroutine

function create_link(self,target,name,merge,owner) result(link)
   class (type_base_model),            intent(inout) :: self
   class (type_internal_object),pointer              :: target
   character(len=*),                   intent(in)    :: name
   logical,optional,                   intent(in)    :: merge,owner

   logical                  :: merge_eff
   type (type_link),pointer :: link

   merge_eff = .false.
   if (present(merge)) merge_eff = merge

   ! First check if a link with this name exists. If so, merge new target with old target.
   link => self%find_link(name)
   if (associated(link)) then
      if (.not.merge_eff) call self%fatal_error('create_link_for_object','Link with name "'//trim(name)//'" already exists.')
      call merge_variables(link%target,target)
      deallocate(target)
      return
   end if

   ! Append a new link to the list.
   if (.not.associated(self%first_link)) then
      allocate(self%first_link)
      link => self%first_link
   else
      link => self%first_link
      do while (associated(link%next))
         link => link%next
      end do
      allocate(link%next)
      link => link%next
   end if

   ! Set link attributes.
   link%name = name
   link%target => target
   if (present(owner)) link%owner = owner
end function create_link

recursive subroutine add_alias(self,target,name)
   class (type_base_model), intent(inout) :: self
   class (type_internal_object),pointer   :: target
   character(len=*),        intent(in)    :: name

   type (type_link),pointer :: link

   link => create_link(self,target,name,owner=.false.)
   if (associated(self%parent)) call self%parent%add_alias(target,trim(self%name)//'/'//trim(name))
end subroutine

subroutine request_coupling_for_link(self,link,master)
   class (type_base_model),intent(inout)              :: self
   type (type_link),       intent(inout)              :: link
   character(len=*),       intent(in)                 :: master

   call self%request_coupling(link%name,master)
end subroutine request_coupling_for_link

subroutine request_coupling_for_name(self,slave,master)
   class (type_base_model),intent(inout)              :: self
   character(len=*),       intent(in)                 :: slave,master

   call self%couplings%set_string(slave,master)
end subroutine request_coupling_for_name

subroutine request_coupling_for_id(self,id,master)
   class (type_base_model), intent(inout) :: self
   class (type_variable_id),intent(inout) :: id
   character(len=*),        intent(in)    :: master

   call self%request_coupling(id%link,master)
end subroutine request_coupling_for_id

subroutine integer_pointer_set_append(self,value)
   class (type_integer_pointer_set),intent(inout) :: self
   integer,target :: value
   type (type_integer_pointer),allocatable :: oldarray(:)

   ! Create a new list of integer pointers, or extend it if already allocated.
   if (.not.allocated(self%pointers)) then
      allocate(self%pointers(1))
   else
      allocate(oldarray(size(self%pointers)))
      oldarray = self%pointers
      deallocate(self%pointers)
      allocate(self%pointers(size(oldarray)+1))
      self%pointers(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   ! Add pointer to provided integer to the list.
   self%pointers(size(self%pointers))%p => value
   self%pointers(size(self%pointers))%p = self%pointers(1)%p
end subroutine integer_pointer_set_append

subroutine integer_pointer_set_extend(self,other)
   class (type_integer_pointer_set),intent(inout) :: self
   class (type_integer_pointer_set),intent(in)    :: other

   integer :: i

   if (allocated(other%pointers)) then
      do i=1,size(other%pointers)
         call self%append(other%pointers(i)%p)
      end do
   end if
end subroutine integer_pointer_set_extend

subroutine integer_pointer_set_set_value(self,value)
   class (type_integer_pointer_set),intent(inout) :: self
   integer,                         intent(in)    :: value

   integer :: i

   do i=1,size(self%pointers)
      self%pointers(i)%p = value
   end do
end subroutine integer_pointer_set_set_value

logical function integer_pointer_set_is_empty(self)
   class (type_integer_pointer_set),intent(in) :: self
   integer_pointer_set_is_empty = .not.allocated(self%pointers)
end function integer_pointer_set_is_empty

subroutine real_pointer_set_append(self,value)
   class (type_real_pointer_set),intent(inout) :: self
   real(rk),target :: value
   type (type_real_pointer),allocatable :: oldarray(:)

   ! Create a new list of real pointers, or extend it if already allocated.
   if (.not.allocated(self%pointers)) then
      allocate(self%pointers(1))
   else
      allocate(oldarray(size(self%pointers)))
      oldarray = self%pointers
      deallocate(self%pointers)
      allocate(self%pointers(size(oldarray)+1))
      self%pointers(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   ! Add pointer to provided real to the list.
   self%pointers(size(self%pointers))%p => value
   self%pointers(size(self%pointers))%p = self%pointers(1)%p
end subroutine real_pointer_set_append

subroutine real_pointer_set_extend(self,other)
   class (type_real_pointer_set),intent(inout) :: self
   class (type_real_pointer_set),intent(in)    :: other

   integer :: i

   if (allocated(other%pointers)) then
      do i=1,size(other%pointers)
         call self%append(other%pointers(i)%p)
      end do
   end if
end subroutine real_pointer_set_extend

subroutine real_pointer_set_set_value(self,value)
   class (type_real_pointer_set),intent(inout) :: self
   real(rk),                     intent(in)    :: value

   integer :: i

   do i=1,size(self%pointers)
      self%pointers(i)%p = value
   end do
end subroutine real_pointer_set_set_value

subroutine append_bulk_data_pointer(array,data)
   type (type_bulk_data_pointer_pointer),allocatable :: array(:)
   type (type_bulk_data_pointer),target :: data
   type (type_bulk_data_pointer_pointer),allocatable :: oldarray(:)
   integer :: i

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not.allocated(array)) then
      allocate(array(1))
   else
      do i=1,size(array)
         if (associated(array(i)%p,data)) return
      end do

      allocate(oldarray(size(array)))
      oldarray = array
      deallocate(array)
      allocate(array(size(oldarray)+1))
      array(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   ! Add pointer to provided data to the list.
   array(size(array))%p => data
   array(size(array))%p = array(1)%p
end subroutine append_bulk_data_pointer

subroutine append_horizontal_data_pointer(array,data)
   type (type_horizontal_data_pointer_pointer),allocatable :: array(:)
   type (type_horizontal_data_pointer),target :: data
   type (type_horizontal_data_pointer_pointer),allocatable :: oldarray(:)
   integer :: i

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not.allocated(array)) then
      allocate(array(1))
   else
      do i=1,size(array)
         if (associated(array(i)%p,data)) return
      end do

      allocate(oldarray(size(array)))
      oldarray = array
      deallocate(array)
      allocate(array(size(oldarray)+1))
      array(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   ! Add pointer to provided data to the list.
   array(size(array))%p => data
   array(size(array))%p = array(1)%p
end subroutine append_horizontal_data_pointer

subroutine append_scalar_data_pointer(array,data)
   type (type_scalar_data_pointer_pointer),allocatable :: array(:)
   type (type_scalar_data_pointer),target :: data
   type (type_scalar_data_pointer_pointer),allocatable :: oldarray(:)
   integer :: i

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not.allocated(array)) then
      allocate(array(1))
   else
      do i=1,size(array)
         if (associated(array(i)%p,data)) return
      end do

      allocate(oldarray(size(array)))
      oldarray = array
      deallocate(array)
      allocate(array(size(oldarray)+1))
      array(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   ! Add pointer to provided data to the list.
   array(size(array))%p => data
   array(size(array))%p = array(1)%p
end subroutine append_scalar_data_pointer

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

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new state variable
!
! !INTERFACE:
   subroutine register_bulk_state_variable(self, id, name, units, long_name, &
                                           initial_value, vertical_movement, specific_light_extinction, &
                                           minimum, maximum, missing_value, &
                                           no_precipitation_dilution, no_river_dilution, &
                                           standard_variable,presence, background_value)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical state variable in the global model database.
!  It returns an identifier that may be used later to retrieve the value of the state variable.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),       intent(inout)        :: self
      type (type_state_variable_id), intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                   intent(in)          :: name, long_name, units
      real(rk),                           intent(in),optional :: initial_value,vertical_movement,specific_light_extinction
      real(rk),                           intent(in),optional :: minimum, maximum,missing_value,background_value
      logical,                            intent(in),optional :: no_precipitation_dilution,no_river_dilution
      type (type_bulk_standard_variable), intent(in),optional :: standard_variable
      integer,                            intent(in),optional :: presence
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable),pointer :: variable
      character(len=256)        :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_bulk_state_variable', &
         'State variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      allocate(variable)
      variable%name      = name
      variable%units     = units
      variable%long_name = long_name
      if (present(initial_value))             variable%initial_value             = initial_value
      if (present(minimum))                   variable%minimum                   = minimum
      if (present(maximum))                   variable%maximum                   = maximum
      if (present(missing_value))             variable%missing_value             = missing_value
      if (present(vertical_movement))         variable%vertical_movement         = vertical_movement
      if (present(no_precipitation_dilution)) variable%no_precipitation_dilution = no_precipitation_dilution
      if (present(no_river_dilution))         variable%no_river_dilution         = no_river_dilution
      if (present(standard_variable))         variable%standard_variable         = standard_variable
      if (present(presence))                  variable%presence                  = presence

      ! Ensure that initial value falls within prescribed valid range.
      if (variable%initial_value<variable%minimum .or. variable%initial_value>variable%maximum) then
         write (text,*) 'Initial value',variable%initial_value,'for variable "'//trim(name)//'" lies&
               &outside allowed range',variable%minimum,'to',variable%maximum
         call self%fatal_error('register_bulk_state_variable',text)
      end if

      call connect_bulk_state_variable_id(self,variable,id)
      if (present(background_value)) call variable%background_values%set_value(background_value)

      id%link => self%add_bulk_variable(variable)

      if (present(specific_light_extinction)) &
         call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,id,scale_factor=specific_light_extinction)
   end subroutine register_bulk_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new state variable
!
! !INTERFACE:
   subroutine register_bottom_state_variable(self, id, name, units, long_name, &
                                             initial_value, minimum, maximum, missing_value, &
                                             standard_variable,presence,background_value)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical state variable in the global model database.
!  It returns an identifier that may be used later to retrieve the value of the state variable.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),             intent(inout)        :: self
      type (type_bottom_state_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                         intent(in)          :: name, long_name, units
      real(rk),                                 intent(in),optional :: initial_value
      real(rk),                                 intent(in),optional :: minimum, maximum,missing_value,background_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
      integer,                                  intent(in),optional :: presence
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable),pointer :: variable
      character(len=256)                      :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_bottom_state_variable', &
         'State variables may only be registered during initialization.')

      allocate(variable)
      variable%name      = name
      variable%units     = units
      variable%long_name = long_name
      variable%domain = domain_bottom
      if (present(initial_value))     variable%initial_value     = initial_value
      if (present(minimum))           variable%minimum           = minimum
      if (present(maximum))           variable%maximum           = maximum
      if (present(missing_value))     variable%missing_value     = missing_value
      if (present(standard_variable)) variable%standard_variable = standard_variable
      if (present(presence))          variable%presence          = presence

      ! Ensure that initial value falls within prescribed valid range.
      if (variable%initial_value<variable%minimum .or. variable%initial_value>variable%maximum) then
         write (text,*) 'Initial value',variable%initial_value,'for variable "'//trim(name)//'" lies&
               &outside allowed range',variable%minimum,'to',variable%maximum
         call self%fatal_error('register_bottom_state_variable',text)
      end if

      call connect_bottom_state_variable_id(self,variable,id)
      if (present(background_value)) call variable%background_values%set_value(background_value)

      id%link => self%add_horizontal_variable(variable)

   end subroutine register_bottom_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new state variable
!
! !INTERFACE:
   subroutine register_surface_state_variable(self, id, name, units, long_name, &
                                              initial_value, minimum, maximum, missing_value, &
                                              standard_variable,presence,background_value)
!
! !DESCRIPTION:
!  This subroutine registers a new surface-bound biogeochemical state variable in the global model database.
!  The identifier "id" may be used later to retrieve the value of the state variable.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),              intent(inout)        :: self
      type (type_surface_state_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                         intent(in)          :: name, long_name, units
      real(rk),                                 intent(in),optional :: initial_value
      real(rk),                                 intent(in),optional :: minimum, maximum,missing_value,background_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
      integer,                                  intent(in),optional :: presence
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable),pointer :: variable
      character(len=256)                      :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_surface_state_variable', &
         'State variables may only be registered during initialization.')

      allocate(variable)
      variable%name      = name
      variable%units     = units
      variable%long_name = long_name
      variable%domain = domain_surface
      if (present(initial_value))     variable%initial_value     = initial_value
      if (present(minimum))           variable%minimum           = minimum
      if (present(maximum))           variable%maximum           = maximum
      if (present(missing_value))     variable%missing_value     = missing_value
      if (present(standard_variable)) variable%standard_variable = standard_variable
      if (present(presence))          variable%presence          = presence

      ! Ensure that initial value falls within prescribed valid range.
      if (variable%initial_value<variable%minimum .or. variable%initial_value>variable%maximum) then
         write (text,*) 'Initial value',variable%initial_value,'for variable "'//trim(name)//'" lies&
               &outside allowed range',variable%minimum,'to',variable%maximum
         call self%fatal_error('register_bottom_state_variable',text)
      end if

      call connect_surface_state_variable_id(self,variable,id)
      if (present(background_value)) call variable%background_values%set_value(background_value)

      id%link => self%add_horizontal_variable(variable)

   end subroutine register_surface_state_variable
!EOC

   function add_bulk_variable(self,variable) result(link)
      ! This subroutine converts the pointer to a bulk variable to a pointer
      ! to a generic object, which can then be added to the model.
      class (type_base_model),target,      intent(inout) :: self
      type (type_bulk_variable),pointer                  :: variable
      type (type_link), pointer :: link
      class (type_internal_object),pointer :: object
      if (variable%name=='' .or. variable%name/=get_safe_name(variable%name)) &
         call self%fatal_error('add_bulk_variable','Cannot register variable because its name is not valid: "'//trim(variable%name)//'".')
      object => variable
      link => add_object(self,object)
      if (.not.associated(object)) nullify(variable)
   end function

   function add_horizontal_variable(self,variable) result(link)
      ! This subroutine converts the pointer to a horizontal variable to a pointer
      ! to a generic object, which can then be added to the model.
      class (type_base_model),target,      intent(inout) :: self
      type (type_horizontal_variable),pointer            :: variable
      type (type_link), pointer :: link
      class (type_internal_object),pointer :: object
      if (variable%name=='' .or. variable%name/=get_safe_name(variable%name)) &
         call self%fatal_error('add_horizontal_variable','Cannot register variable because its name is not valid: "'//trim(variable%name)//'".')
      object => variable
      link => add_object(self,object)
      if (.not.associated(object)) nullify(variable)
   end function

   function add_scalar_variable(self,variable) result(link)
      ! This subroutine converts the pointer to a scalar variable to a pointer
      ! to a generic object, which can then be added to the model.
      class (type_base_model),target,      intent(inout) :: self
      type (type_scalar_variable),pointer                :: variable
      type (type_link), pointer :: link
      class (type_internal_object),pointer :: object
      if (variable%name=='' .or. variable%name/=get_safe_name(variable%name)) &
         call self%fatal_error('add_scalar_variable','Cannot register variable because its name is not valid: "'//trim(variable%name)//'".')
      object => variable
      link => add_object(self,object)
      if (.not.associated(object)) nullify(variable)
   end function

   recursive function add_object(self,object) result(link)
      ! This subroutine creates a link to the supplied object, then allows
      ! parent models to do the same.
      ! NB this subroutine MUST be recursive, to allow parent models to override
      ! the properties of objects added by their child models.
      class (type_base_model),target,      intent(inout) :: self
      class (type_internal_object),pointer               :: object

      type (type_link), pointer :: link,parent_link

      if (object%long_name=='') object%long_name = object%name

      link => create_link(self,object,object%name,merge=.true.)
      if (.not.associated(object)) return ! if create_link has merged the new object into an existing object

      ! Forward to parent
      if (associated(self%parent)) then
         object%name = trim(self%name)//'/'//trim(object%name)
         object%long_name = trim(self%long_name)//' '//trim(object%long_name)
         parent_link => self%parent%add_object(object)
      end if
   end function
   
   subroutine connect_bulk_state_variable_id(self,variable,id)
      class (type_base_model),       intent(in)           :: self
      type (type_bulk_variable),     intent(inout),target :: variable
      type (type_state_variable_id), intent(inout),target :: id

      if (associated(id%link)) call self%fatal_error('connect_bulk_state_variable_id', &
         'Identifier supplied for '//trim(variable%name)//' is already associated with '//trim(id%link%name)//'.')

      call variable%state_indices%append(id%state_index)
      call append_data_pointer(variable%alldata,id%data)
      call variable%background_values%append(id%background)
   end subroutine

   subroutine connect_bottom_state_variable_id(self,variable,id)
      class (type_base_model),       intent(in)           :: self
      type (type_horizontal_variable),      intent(inout),target :: variable
      type (type_bottom_state_variable_id), intent(inout),target :: id

      if (associated(id%link)) call self%fatal_error('connect_bottom_state_variable_id', &
         'Identifier supplied for '//trim(variable%name)//' is already associated with '//trim(id%link%name)//'.')

      call variable%state_indices%append(id%bottom_state_index)
      call append_data_pointer(variable%alldata,id%horizontal_data)
      call variable%background_values%append(id%background)
   end subroutine

   subroutine connect_surface_state_variable_id(self,variable,id)
      class (type_base_model),       intent(in)           :: self
      type (type_horizontal_variable),      intent(inout),target :: variable
      type (type_surface_state_variable_id),intent(inout),target :: id

      if (associated(id%link)) call self%fatal_error('connect_surface_state_variable_id', &
         'Identifier supplied for '//trim(variable%name)//' is already associated with '//trim(id%link%name)//'.')

      call variable%state_indices%append(id%surface_state_index)
      call append_data_pointer(variable%alldata,id%horizontal_data)
      call variable%background_values%append(id%background)
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   subroutine register_bulk_diagnostic_variable(self, id, name, units, long_name, &
                                                time_treatment, missing_value, standard_variable, output)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical diagnostic variable in the global model database.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),           intent(inout),target :: self
      type (type_diagnostic_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                   intent(in)          :: name, long_name, units
      integer,                            intent(in),optional :: time_treatment, output
      real(rk),                           intent(in),optional :: missing_value
      type (type_bulk_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable),pointer :: variable
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_bulk_diagnostic_variable',&
                                             'Diagnostic variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      allocate(variable)
      variable%name      = name
      variable%units     = units
      variable%long_name = long_name
      variable%source_model => self
      if (present(time_treatment)) then
         variable%output = time_treatment2output(time_treatment)
         call self%log_message('variable "'//trim(name)//'": "time_treatment" argument to register_diagnostic_variable is deprecated; &
                               &please use "output" instead (see output_* parameters near top of fabm_types.F90.')
      end if
      if (present(missing_value))     variable%missing_value     = missing_value
      if (present(standard_variable)) variable%standard_variable = standard_variable
      if (present(output))            variable%output            = output
      call variable%write_indices%append(id%diag_index)
      if (associated(id%link)) call self%fatal_error('register_bulk_diagnostic_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      id%link => self%add_bulk_variable(variable)

   end subroutine register_bulk_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   subroutine register_horizontal_diagnostic_variable(self, id, name, units, long_name, &
                                                      time_treatment, missing_value, standard_variable, output)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical diagnostic variable in the global model database.
!
! !INPUT/OUTPUT PARAMETER:
      class (type_base_model),                      intent(inout),target :: self
      type (type_horizontal_diagnostic_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                         intent(in)          :: name, long_name, units
      integer,                                  intent(in),optional :: time_treatment, output
      real(rk),                                 intent(in),optional :: missing_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable),pointer :: variable
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_horizontal_diagnostic_variable',&
                                             'Diagnostic variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      allocate(variable)
      variable%name      = name
      variable%units     = units
      variable%long_name = long_name
      variable%source_model => self
      if (present(time_treatment)) then
         variable%output = time_treatment2output(time_treatment)
         call self%log_message('variable "'//trim(name)//'": "time_treatment" argument to register_diagnostic_variable is deprecated; &
                               &please use "output" instead (see output_* parameters near top of fabm_types.F90.')
      end if
      if (present(missing_value))     variable%missing_value     = missing_value
      if (present(standard_variable)) variable%standard_variable = standard_variable
      if (present(output))            variable%output            = output

      if (associated(id%link)) call self%fatal_error('register_horizontal_diagnostic_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')
      call variable%write_indices%append(id%horizontal_diag_index)

      id%link => self%add_horizontal_variable(variable)

   end subroutine register_horizontal_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a conserved quantity taken from a standard set
!
! !INTERFACE:
   subroutine register_standard_conserved_quantity(self, id, standard_variable, name)
!
! !DESCRIPTION:
!  This function registers a new biogeochemically conserved quantity in the global
!  model database.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),          intent(inout)        :: self
      type (type_conserved_quantity_id),intent(inout),target :: id
      character(len=*),                 intent(in), optional :: name
!
! !INPUT PARAMETERS:
      type (type_bulk_standard_variable), intent(in) :: standard_variable
!
!-----------------------------------------------------------------------
!BOC
      call self%fatal_error('register_standard_conserved_quantity','register_conserved_quantity is no longer supported; please use add_to_aggregate_variable.')

   end subroutine register_standard_conserved_quantity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new conserved quantity
!
! !INTERFACE:
   subroutine register_custom_conserved_quantity(self, id, name, units, long_name)
!
! !DESCRIPTION:
!  This function registers a new biogeochemically conserved quantity in the global
!  model database.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),          intent(inout)        :: self
      type (type_conserved_quantity_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: long_name
      character(len=*), intent(in) :: units
!
!-----------------------------------------------------------------------
!BOC
      call self%fatal_error('register_standard_conserved_quantity','register_conserved_quantity is no longer supported; please use add_to_aggregate_variable.')

   end subroutine register_custom_conserved_quantity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external state variable
!
! !INTERFACE:
   subroutine register_bulk_state_dependency_ex(self,id,name,units,long_name,required)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),      intent(inout)        :: self
      type (type_state_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),intent(in)          :: name,units,long_name
      logical,         intent(in),optional :: required
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
      logical :: required_eff
      integer :: presence
!-----------------------------------------------------------------------
!BOC
      required_eff = .true.
      if (present(required)) required_eff = required

      presence = presence_external_required
      if (.not.required_eff) presence = presence_external_optional
      call register_bulk_state_variable(self, id, name, units, long_name, presence=presence)

   end subroutine register_bulk_state_dependency_ex
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external bottom state variable
!
! !INTERFACE:
   subroutine register_bottom_state_dependency_ex(model,id,name,units,long_name,required)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),             intent(inout)        :: model
      type (type_bottom_state_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),intent(in)          :: name,units,long_name
      logical,         intent(in),optional :: required
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
      logical :: required_eff
      integer :: presence
!-----------------------------------------------------------------------
!BOC
      required_eff = .true.
      if (present(required)) required_eff = required

      presence = presence_external_required
      if (.not.required_eff) presence = presence_external_optional
      call register_bottom_state_variable(model, id, name, units, long_name, presence=presence)

   end subroutine register_bottom_state_dependency_ex
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external surface-bound state variable
!
! !INTERFACE:
   subroutine register_surface_state_dependency_ex(model,id,name,units,long_name,required)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),              intent(inout)        :: model
      type (type_surface_state_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),intent(in)          :: name,units,long_name
      logical,         intent(in),optional :: required
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
      logical :: required_eff
      integer :: presence
!-----------------------------------------------------------------------
!BOC
      required_eff = .true.
      if (present(required)) required_eff = required

      presence = presence_external_required
      if (.not.required_eff) presence = presence_external_optional
      call register_surface_state_variable(model, id, name, units, long_name, presence=presence)

   end subroutine register_surface_state_dependency_ex
!EOC

   subroutine register_bulk_state_dependency_old(self,id,name)
      class (type_base_model),      intent(inout)        :: self
      type (type_state_variable_id),intent(inout),target :: id
      character(len=*),             intent(in)           :: name

      call self%register_bulk_state_dependency(id, name, '', name)
      call self%request_coupling(id,name)
   end subroutine

   subroutine register_bottom_state_dependency_old(self,id,name)
      class (type_base_model),             intent(inout)        :: self
      type (type_bottom_state_variable_id),intent(inout),target :: id
      character(len=*),                    intent(in)           :: name

      call self%register_bottom_state_dependency(id, name, '', name)
      call self%request_coupling(id,name)
   end subroutine

   subroutine register_surface_state_dependency_old(self,id,name)
      class (type_base_model),              intent(inout)        :: self
      type (type_surface_state_variable_id),intent(inout),target :: id
      character(len=*),                     intent(in)           :: name

      call self%register_surface_state_dependency(id, name, '', name)
      call self%request_coupling(id,name)
   end subroutine
   
   subroutine register_standard_bulk_dependency(model,id,standard_variable,required)
      class (type_base_model),           intent(inout)        :: model
      type (type_dependency_id),         intent(inout),target :: id
      type (type_bulk_standard_variable),intent(in)           :: standard_variable
      logical,optional,                  intent(in)           :: required

      call register_named_bulk_dependency(model,id,standard_variable%name,standard_variable%units, &
                                          standard_variable=standard_variable,required=required)
   end subroutine

   subroutine register_standard_horizontal_dependency(model,id,standard_variable,required)
      class (type_base_model),                 intent(inout)        :: model
      type (type_horizontal_dependency_id),    intent(inout),target :: id
      type (type_horizontal_standard_variable),intent(in)           :: standard_variable
      logical,optional,                        intent(in)           :: required

      call register_named_horizontal_dependency(model,id,standard_variable%name,standard_variable%units, &
                                                standard_variable=standard_variable,required=required)
   end subroutine

   subroutine register_standard_global_dependency(model,id,standard_variable,required)
      class (type_base_model),             intent(inout)        :: model
      type (type_global_dependency_id),    intent(inout),target :: id
      type (type_global_standard_variable),intent(in)           :: standard_variable
      logical,optional,                    intent(in)           :: required

      call register_named_global_dependency(model,id,standard_variable%name,standard_variable%units, &
                                            standard_variable=standard_variable,required=required)
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a variable defined on
! the full model domain.
!
! !INTERFACE:
   subroutine register_named_bulk_dependency(self,id,name,units,long_name,standard_variable,required)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),  intent(inout)        :: self
      type (type_dependency_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                  intent(in)          :: name
      character(len=*),                  intent(in),optional :: units,long_name
      type (type_bulk_standard_variable),intent(in),optional :: standard_variable
      logical,                           intent(in),optional :: required
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable),pointer :: variable
      logical                           :: required_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_bulk_dependency',&
         'Dependencies may only be registered during initialization.')

      required_eff = .true.
      if (present(required)) required_eff = required

      allocate(variable)
      variable%name = name
      if (present(units))     variable%units     = units
      if (present(long_name)) variable%long_name = long_name
      if (present(standard_variable)) variable%standard_variable = standard_variable
      variable%presence = presence_external_required
      if (.not.required_eff) variable%presence = presence_external_optional

      call append_data_pointer(variable%alldata,id%data)
      call variable%background_values%append(id%background)
      if (associated(id%link)) call self%fatal_error('register_bulk_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      id%link => self%add_bulk_variable(variable)

      !if (associated(self%parent).and..not.present(standard_variable)) call self%request_coupling(id,name,required=.false.)

   end subroutine register_named_bulk_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a variable defined on
! a horizontal slice of the model domain.
!
! !INTERFACE:
   subroutine register_named_horizontal_dependency(self,id,name,units,long_name,standard_variable,required)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),             intent(inout)        :: self
      type (type_horizontal_dependency_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                        intent(in)          :: name
      character(len=*),                        intent(in),optional :: units,long_name
      type (type_horizontal_standard_variable),intent(in),optional :: standard_variable
      logical,                                 intent(in),optional :: required
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable),pointer :: variable
      logical                                 :: required_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_horizontal_dependency',&
         'Dependencies may only be registered during initialization.')

      required_eff = .true.
      if (present(required)) required_eff = required

      allocate(variable)
      variable%name = name
      if (present(units))     variable%units     = units
      if (present(long_name)) variable%long_name = long_name
      if (present(standard_variable)) variable%standard_variable = standard_variable
      variable%presence = presence_external_required
      if (.not.required_eff) variable%presence = presence_external_optional

      call append_data_pointer(variable%alldata,id%horizontal_data)
      call variable%background_values%append(id%background)
      if (associated(id%link)) call self%fatal_error(':register_horizontal_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      id%link => self%add_horizontal_variable(variable)

      !if (associated(self%parent).and..not.present(standard_variable)) call self%request_coupling(id,name,required=.false.)

   end subroutine register_named_horizontal_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a global (space-
! independent) variable.
!
! !INTERFACE:
   subroutine register_named_global_dependency(self,id,name,units,long_name,standard_variable,required)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),         intent(inout)        :: self
      type (type_global_dependency_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                    intent(in)          :: name
      character(len=*),                    intent(in),optional :: units,long_name
      type (type_global_standard_variable),intent(in),optional :: standard_variable
      logical,                             intent(in),optional :: required
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_scalar_variable),pointer :: variable
      logical                             :: required_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error(':register_global_dependency',&
         'Dependencies may only be registered during initialization.')

      required_eff = .true.
      if (present(required)) required_eff = required

      allocate(variable)
      variable%name = name
      if (present(units))     variable%units     = units
      if (present(long_name)) variable%long_name = long_name
      if (present(standard_variable)) variable%standard_variable = standard_variable
      variable%presence = presence_external_required
      if (.not.required_eff) variable%presence = presence_external_optional

      call append_data_pointer(variable%alldata,id%global_data)
      call variable%background_values%append(id%background)
      if (associated(id%link)) call self%fatal_error(':register_global_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      id%link => self%add_scalar_variable(variable)

      !if (associated(self%parent).and..not.present(standard_variable)) call self%request_coupling(id,name,required=.false.)

   end subroutine register_named_global_dependency
!EOC

subroutine register_bulk_expression_dependency(self,id,expression)
   class (type_base_model),       intent(inout)        :: self
   type (type_dependency_id),     intent(inout),target :: id
   class (type_bulk_expression),  intent(in)           :: expression

   class (type_bulk_expression),allocatable :: copy

   allocate(copy,source=expression)
   copy%out => id%data
   call self%register_dependency(id,trim(copy%output_name))
   copy%output_name = id%link%target%name

   call register_expression(self,copy)
   deallocate(copy)
end subroutine

subroutine register_horizontal_expression_dependency(self,id,expression)
   class (type_base_model),             intent(inout)        :: self
   type (type_horizontal_dependency_id),intent(inout),target :: id
   class (type_horizontal_expression),  intent(in)           :: expression

   class (type_horizontal_expression),allocatable :: copy

   allocate(copy,source=expression)
   copy%out => id%horizontal_data
   call self%register_dependency(id,trim(copy%output_name))
   copy%output_name = id%link%target%name

   call register_expression(self,copy)
   deallocate(copy)
end subroutine

recursive subroutine register_expression(self,expression)
   class (type_base_model), intent(inout)           :: self
   class (type_expression), intent(in)              :: expression

   class (type_expression), pointer :: current

   if (.not.associated(self%first_expression)) then
      allocate(self%first_expression,source=expression)
      current => self%first_expression
   else
      current => self%first_expression
      do while (associated(current%next))
         current => current%next
      end do
      allocate(current%next,source=expression)
      current => current%next
   end if

   if (associated(self%parent)) call register_expression(self%parent,expression)
end subroutine

subroutine get_real_parameter(self,value,name,units,long_name,default,scale_factor)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target  :: self
   real(rk),                intent(inout)          :: value
   character(len=*),        intent(in)             :: name
   character(len=*),        intent(in),   optional :: units,long_name
   real(rk),                intent(in),   optional :: default,scale_factor
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),    pointer :: property
   logical                           :: success
   type (type_real_property)         :: current_parameter
!
!-----------------------------------------------------------------------
!BOC
   if (present(default)) then
      current_parameter%has_default = .true.
      current_parameter%default = default
      value = default
   end if

   ! Try to find a user-specified value for this parameter in our dictionary, and in those of our ancestors.
   property => self%parameters%find_in_tree(name)
   if (associated(property)) then
      ! Value found - try to convert to real.
      value = property%to_real(success=success)
      if (.not.success) call self%fatal_error('get_real_parameter', &
         'Value "'//trim(property%to_string())//'" for parameter "'//trim(name)//'" is not a real number.')
   elseif (.not.present(default)) then
      call self%fatal_error('get_real_parameter','No value provided for parameter "'//trim(name)//'".')
   end if

   ! Store parameter settings
   current_parameter%value = value
   call set_parameter(self,current_parameter,name,units,long_name)

   ! Apply scale factor to value provided to the model (if requested).
   if (present(scale_factor)) value = value*scale_factor
end subroutine get_real_parameter
!EOC

subroutine set_parameter(self,parameter,name,units,long_name)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target :: self
   class (type_property),   intent(inout)         :: parameter
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
!
!EOP
   class (type_base_model), pointer :: pmodel
!-----------------------------------------------------------------------
!BOC
   parameter%name = name
   if (present(units))     parameter%units     = units
   if (present(long_name)) parameter%long_name = long_name
   call self%parameters%set_in_tree(parameter)
end subroutine set_parameter
!EOC

subroutine get_integer_parameter(self,value,name,units,long_name,default)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target :: self
   integer,                 intent(inout)         :: value
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
   integer,                 intent(in),optional   :: default
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),       pointer :: property
   type (type_integer_property)         :: current_parameter
   logical                              :: success
!
!-----------------------------------------------------------------------
!BOC
   if (present(default)) then
      current_parameter%has_default = .true.
      current_parameter%default = default
      value = default
   end if

   ! Try to find a user-specified value for this parameter in our dictionary, and in those of our ancestors.
   property => self%parameters%find_in_tree(name)
   if (associated(property)) then
      ! Value found - try to convert to integer.
      value = property%to_integer(success=success)
      if (.not.success) call self%fatal_error('get_integer_parameter', &
         'Value "'//trim(property%to_string())//'" for parameter "'//trim(name)//'" is not an integer number.')
   elseif (.not.present(default)) then
      call self%fatal_error('get_integer_parameter','No value provided for parameter "'//trim(name)//'".')
   end if

   ! Store parameter settings
   current_parameter%value = value
   call set_parameter(self,current_parameter,name,units,long_name)
end subroutine get_integer_parameter
!EOC

subroutine get_logical_parameter(self,value,name,units,long_name,default)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target :: self
   logical,                 intent(inout)         :: value
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
   logical,                 intent(in),optional   :: default
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),       pointer :: property
   type (type_logical_property)         :: current_parameter
   logical                              :: success
!
!-----------------------------------------------------------------------
!BOC
   if (present(default)) then
      current_parameter%has_default = .true.
      current_parameter%default = default
      value = default
   end if

   ! Try to find a user-specified value for this parameter in our dictionary, and in those of our ancestors.
   property => self%parameters%find_in_tree(name)
   if (associated(property)) then
      ! Value found - try to convert to logical.
      value = property%to_logical(success=success)
      if (.not.success) call self%fatal_error('get_logical_parameter', &
         'Value "'//trim(property%to_string())//'" for parameter "'//trim(name)//'" is not a Boolean value.')
   elseif (.not.present(default)) then
      call self%fatal_error('get_logical_parameter','No value provided for parameter "'//trim(name)//'".')
   end if

   ! Store parameter settings
   current_parameter%value = value
   call set_parameter(self,current_parameter,name,units,long_name)
end subroutine get_logical_parameter
!EOC

recursive subroutine get_string_parameter(self,value,name,units,long_name,default)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target :: self
   character(len=*),        intent(inout)         :: value
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
   character(len=*),        intent(in),optional   :: default
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),      pointer :: property
   type (type_string_property)         :: current_parameter
   logical                             :: success
!
!-----------------------------------------------------------------------
!BOC
   if (present(default)) then
      current_parameter%has_default = .true.
      current_parameter%default = default
      value = default
   end if

   ! Try to find a user-specified value for this parameter in our dictionary, and in those of our ancestors.
   property => self%parameters%find_in_tree(name)
   if (associated(property)) then
      ! Value found - try to convert to string.
      value = property%to_string(success=success)
      if (.not.success) call self%fatal_error('get_string_parameter', &
         'Value for parameter "'//trim(name)//'" cannot be converted to string.')
   elseif (.not.present(default)) then
      call self%fatal_error('get_string_parameter','No value provided for parameter "'//trim(name)//'".')
   end if

   ! Store parameter settings
   current_parameter%value = value
   call set_parameter(self,current_parameter,name,units,long_name)
end subroutine get_string_parameter
!EOC

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
                                    call couple_variables(model,link%target,link2%target)
                                 else
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
                                    call couple_variables(model,link%target,link2%target)
                                 else
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
                                    call couple_variables(model,link%target,link2%target)
                                 else
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

subroutine merge_variables(master,slave)
   class (type_internal_object),intent(inout) :: master
   class (type_internal_object),intent(in)    :: slave

   call log_message(trim(slave%name)//' --> '//trim(master%name))
   select type (master)
      class is (type_bulk_variable)
         select type (slave)
            class is (type_bulk_variable)
               call merge_bulk_variables(master,slave)
            class is (type_internal_object) ! class default would be preferable but breaks the next line with Cray 8.1
               call fatal_error('merge_variables', &
                  'type mismatch: '//trim(master%name)//' is defined on the whole domain, '//trim(slave%name)//' is not.')
         end select      
      class is (type_horizontal_variable)
         select type (slave)
            class is (type_horizontal_variable)
               call merge_horizontal_variables(master,slave)
            class is (type_internal_object) ! class default would be preferable but breaks the next line with Cray 8.1
               call fatal_error('merge_variables', &
                  'type mismatch: '//trim(master%name)//' is defined on the horizontal domain, '//trim(slave%name)//' is not.')
         end select      
      class is (type_scalar_variable)
         select type (slave)
            class is (type_scalar_variable)
               call merge_scalar_variables(master,slave)
            class is (type_internal_object) ! class default would be preferable but breaks the next line with Cray 8.1
               call fatal_error('merge_variables', &
                  'type mismatch: '//trim(master%name)//' is defined as a scalar, '//trim(slave%name)//' is not.')
         end select      
   end select
end subroutine

subroutine merge_internal_variables(master,slave)
   class (type_internal_variable),intent(inout) :: master
   class (type_internal_variable),intent(in)    :: slave

   type (type_contribution), pointer :: contribution

   if (.not.slave%write_indices%is_empty()) &
      call fatal_error('merge_internal_variables','Attempt to couple write-only variable ' &
         //trim(slave%name)//' to '//trim(master%name)//'.')
   if (master%state_indices%is_empty().and..not.slave%state_indices%is_empty()) &
      call fatal_error('merge_internal_variables','Attempt to couple state variable ' &
         //trim(slave%name)//' to non-state variable '//trim(master%name)//'.')
   if (master%presence==presence_external_optional) &
      call fatal_error('merge_internal_variables','Attempt to couple to optional master variable "'//trim(master%name)//'".')

   call master%state_indices%extend(slave%state_indices)
   call master%background_values%extend(slave%background_values)
   call master%properties%update(slave%properties,overwrite=.false.)
   contribution => slave%contributions%first
   do while (associated(contribution))
      call master%contributions%add(contribution%target,contribution%scale_factor)
      contribution => contribution%next
   end do
end subroutine

subroutine merge_bulk_variables(master,slave)
   type (type_bulk_variable),intent(inout) :: master
   type (type_bulk_variable),intent(in)    :: slave
   integer :: i

   call merge_internal_variables(master,slave)
   if (allocated(slave%alldata)) then
      do i=1,size(slave%alldata)
         call append_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
end subroutine merge_bulk_variables

subroutine merge_horizontal_variables(master,slave)
   type (type_horizontal_variable),intent(inout) :: master
   type (type_horizontal_variable),intent(in)    :: slave
   integer :: i
   
   if (slave%domain/=master%domain) then
      call fatal_error('merge_horizontal_variables','Domains of coupled variabled ' &
         //trim(slave%name)//' to '//trim(master%name)//' do not match.')
   end if
   call merge_internal_variables(master,slave)
   if (allocated(slave%alldata)) then
      do i=1,size(slave%alldata)
         call append_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
end subroutine merge_horizontal_variables

subroutine merge_scalar_variables(master,slave)
   type (type_scalar_variable),intent(inout) :: master
   type (type_scalar_variable),intent(in)    :: slave
   integer :: i

   call merge_internal_variables(master,slave)
   if (allocated(slave%alldata)) then
      do i=1,size(slave%alldata)
         call append_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
end subroutine merge_scalar_variables

   function find_object(self,name,recursive,exact) result(object)

      class (type_base_model),  intent(in),target :: self
      character(len=*),         intent(in)        :: name
      logical,         optional,intent(in)        :: recursive,exact
      class (type_internal_object),pointer        :: object

      type (type_link), pointer :: link

      nullify(object)
      link => self%find_link(name,recursive,exact)
      if (associated(link)) object => link%target
      
   end function find_object

   function find_link(self,name,recursive,exact) result(link)

      class (type_base_model),  intent(in),target :: self
      character(len=*),         intent(in)        :: name
      logical,         optional,intent(in)        :: recursive,exact
      type (type_link),pointer                    :: link

      logical                         :: recursive_eff,exact_eff
      class (type_base_model),pointer :: current

      recursive_eff = .false.
      if (present(recursive)) recursive_eff = recursive
      nullify(link)

      ! First search self and ancestors (if allowed) based on exact name provided.
      current => self
      do while (associated(current))
         link => current%first_link
         do while (associated(link))
            if (link%name==name) return
            link => link%next
         end do
         if (.not.recursive_eff) exit
         current => current%parent
      end do

      exact_eff = .true.
      if (present(exact)) exact_eff = exact
      if (exact_eff) return

      ! Not found. Now search self and ancestors (if allowed) based on safe name (letters and underscores only).
      current => self
      do while (associated(current))
         link => current%first_link
         do while (associated(link))
            if (get_safe_name(link%name)==name) return
            link => link%next
         end do
         if (.not.recursive_eff) exit
         current => current%parent
      end do

   end function find_link

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Find a model by name.
!
! !INTERFACE:
   recursive function find_model(self,name,recursive) result(found_model)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_base_model),       intent(in),target :: self
      character(len=*),              intent(in)        :: name
      logical,optional,              intent(in)        :: recursive
      class (type_base_model),pointer                  :: found_model

      logical                             :: recursive_eff
      type (type_model_list_node),pointer :: node
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
      nullify(found_model)

      node => self%children%first
      do while (associated(node))
         if (node%model%name==name) then
            found_model => node%model
            return
         end if
         node => node%next
      end do

      ! Search among children of child
      recursive_eff = .false.
      if (present(recursive)) recursive_eff = recursive
      if (recursive_eff.and.associated(self%parent)) found_model => self%parent%find_model(name,recursive)
   end function find_model
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
            ! Enumerate its contributions to aggreagte variables, and register these with
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

function get_aggregate_variable(self,standard_variable,create) result(aggregate_variable)
   class (type_base_model),           intent(inout) :: self
   type (type_bulk_standard_variable),intent(in)    :: standard_variable
   logical,optional,                  intent(in)    :: create

   logical                                :: create_eff
   type (type_aggregate_variable),pointer :: aggregate_variable

   create_eff = .true.
   if (present(create)) create_eff = create

   nullify(aggregate_variable)
   if (.not.associated(self%first_aggregate_variable)) then
      ! No aggregate variables yet. Create one for the target variable.
      if (.not.create_eff) return
      allocate(self%first_aggregate_variable)
      aggregate_variable => self%first_aggregate_variable
   else
      ! First check whether we have already created the desired aggregate variable.
      aggregate_variable => self%first_aggregate_variable
      do while (associated(aggregate_variable))
         if (aggregate_variable%standard_variable%compare(standard_variable)) return
         aggregate_variable => aggregate_variable%next
      end do
      if (.not.create_eff) return

      ! Aggregate variable does not exist yet. Create it.
      aggregate_variable => self%first_aggregate_variable
      do while (associated(aggregate_variable%next))
         aggregate_variable => aggregate_variable%next
      end do
      allocate(aggregate_variable%next)
      aggregate_variable => aggregate_variable%next
   end if

   ! Associate the newly created aggregate variable with the requested standard variable.
   aggregate_variable%standard_variable = standard_variable

   ! Make sure that aggregate variables at the root level are computed.
   ! These are typically used by the host to check conservation.
   if (.not.associated(self%parent)) then
      aggregate_variable%bulk_required = .true.
      aggregate_variable%horizontal_required = .true.
   end if

end function

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
   class (type_base_model),pointer           :: root

   ! This procedure takes an aggregate variable, and creates models that compute diagnostics
   ! for the total of these aggregate quantities on bulk and horizontal domains.

   select case (domain)
      case (domain_bulk)
         name = trim(aggregate_variable%standard_variable%name)
      case default
         name = trim(aggregate_variable%standard_variable%name)//'_at_interfaces'
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

   if (n==0) then
      ! Always zero - create alias to zero field at root level
      root => self
      do while (associated(root%parent))
         root => root%parent
      end do
      select case (domain)
         case (domain_bulk)
            if (.not._VARIABLE_REGISTERED_(root%id_zero)) &
               call root%register_dependency(root%id_zero,'zero','-','zero')
            call self%add_alias(root%id_zero%link%target,name)
         case default
            if (.not._VARIABLE_REGISTERED_(root%id_zero_hz)) &
               call root%register_horizontal_dependency(root%id_zero_hz,'zero_hz','-','zero')
            call self%add_alias(root%id_zero_hz%link%target,name)
      end select            
   elseif (n==1.and..not.sumrequired) then
      ! Always equal to another - create alias to this other
      call self%add_alias(target_variable,name)
   else
      ! Weighted sum of variables
      select case (domain)
         case (domain_bulk)
            allocate(aggregate_variable%sum)
            aggregate_variable%sum%output_name = name
            aggregate_variable%sum%output_units = trim(aggregate_variable%standard_variable%units)
         case default
            allocate(aggregate_variable%horizontal_sum)
            aggregate_variable%horizontal_sum%output_name = name
            aggregate_variable%horizontal_sum%output_units = trim(aggregate_variable%standard_variable%units)//'*m'
      end select

      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (contributing_variable%link%owner) then
            select type (variable=>contributing_variable%link%target)
               class is (type_bulk_variable)
                  ! This contribution comes from a bulk variable.
                  if (domain==domain_bulk) &
                     call aggregate_variable%sum%add_component(trim(contributing_variable%link%name), &
                        weight=contributing_variable%scale_factor, include_background=contributing_variable%include_background)
               class is (type_horizontal_variable)
                  ! This contribution comes from a variable defined on a horizontal interface (top or bottom).
                  if (domain/=domain_bulk) &
                     call aggregate_variable%horizontal_sum%add_component(trim(contributing_variable%link%name), &
                        weight=contributing_variable%scale_factor,include_background=contributing_variable%include_background)
            end select
         end if
         contributing_variable => contributing_variable%next
      end do

      select case (domain)
         case (domain_bulk)
            call self%add_child(aggregate_variable%sum,trim(aggregate_variable%sum%output_name)//'_calculator',configunit=-1)
         case default
            call self%add_child(aggregate_variable%horizontal_sum,trim(aggregate_variable%horizontal_sum%output_name)//'_calculator',configunit=-1)
      end select
   end if
end subroutine create_aggregate_model

function get_free_unit() result(unit)
   integer :: unit
   integer, parameter :: LUN_MIN=10, LUN_MAX=1000

   logical :: opened

   do unit=LUN_MIN,LUN_MAX
      inquire(unit=unit,opened=opened)
      if (.not.opened) return
   end do
   unit = -1
end function get_free_unit

function get_safe_name(name) result(safe_name)
   character(len=*),intent(in) :: name
   character(len=len(name))    :: safe_name
   integer :: i,ch
   logical :: valid
   safe_name = name
   do i=1,len_trim(name)
      ch = iachar(name(i:i))
      valid = (ch>=iachar('a').and.ch<=iachar('z')) & ! Lower-case letter
          .or.(ch>=iachar('A').and.ch<=iachar('Z')) & ! Upper-case letter
          .or.(ch>=iachar('0').and.ch<=iachar('9')) & ! Number
          .or.(ch==iachar('_'))                       ! Underscore
      if (.not.valid) safe_name(i:i) = '_'
   end do
end function

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
   link => self%first_link
   do while (associated(link))
      select type (object => link%target)
         class is (type_internal_variable)
            if (associated(object%source_model).and..not.associated(object%source_model,self)) &
               call find_dependencies(object%source_model,list,forbidden_with_self)
      end select
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

subroutine abstract_model_factory_add(self,child)
   class (type_base_model_factory),       intent(inout) :: self
   class (type_base_model_factory),target,intent(in)    :: child

   class (type_base_model_factory),pointer :: current

   if (.not.associated(self%first_child)) then
      self%first_child => child
   else
      current => self%first_child
      do while(associated(current%next_sibling))
         current => current%next_sibling
      end do
      current%next_sibling => child
   end if
end subroutine

recursive subroutine abstract_model_factory_create(self,name,model)
   class (type_base_model_factory),intent(in) :: self
   character(len=*),               intent(in) :: name
   class (type_base_model),pointer            :: model

   class (type_base_model_factory),pointer :: child

   child => self%first_child
   do while(associated(child))
      call child%create(name,model)
      if (associated(model)) return
      child => child%next_sibling
   end do
end subroutine

   subroutine weighted_sum_initialize(self,configunit)
      class (type_weighted_sum),intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      type (type_component),pointer :: component
      integer           :: i
      character(len=10) :: temp

      i = 0
      component => self%first
      do while (associated(component))
         i = i + 1
         write (temp,'(i0)') i
         call self%register_dependency(component%id,'term'//trim(temp))
         call self%request_coupling(component%id,trim(component%name))
         component => component%next
      end do
      if (self%output_long_name=='') self%output_long_name = self%output_name
      call self%register_diagnostic_variable(self%id_output,'result',self%output_units,'result')
      call self%parent%add_alias(self%id_output%link%target,trim(self%output_name))
   end subroutine

   subroutine weighted_sum_add_component(self,name,weight,include_background)
      class (type_weighted_sum),intent(inout) :: self
      character(len=*),         intent(in)    :: name
      real(rk),optional,        intent(in)    :: weight
      logical,optional,         intent(in)    :: include_background

      type (type_component),pointer :: component

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

   subroutine weighted_sum_after_coupling(self)
      class (type_weighted_sum),intent(inout) :: self

      type (type_component),pointer :: component

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      component => self%first
      do while (associated(component))
         if (component%include_background) self%offset = self%offset + component%weight*component%id%background
         component => component%next
      end do
   end subroutine

   subroutine weighted_sum_evaluate(self,_ARGUMENTS_ND_)
      class (type_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_ND_

      type (type_component),pointer        :: component
      real(rk)                             :: value
      real(rk) _DIMENSION_SLICE_AUTOMATIC_ :: sum

      ! Initialize sum to starting value (typically zero).
      sum = self%offset

      ! Enumerate components included in the sum, and add their contributions.
      component => self%first
      do while (associated(component))
         _LOOP_BEGIN_
            _GET_(component%id,value)
            sum _INDEX_SLICE_ = sum _INDEX_SLICE_ + component%weight*value
         _LOOP_END_
         component => component%next
      end do

      ! Transfer summed values to diagnostic.
      _LOOP_BEGIN_
         _SET_DIAGNOSTIC_(self%id_output,sum _INDEX_SLICE_)
      _LOOP_END_
   end subroutine

   subroutine weighted_sum_do(self,_ARGUMENTS_DO_)
      class (type_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      call self%evaluate(_ARGUMENTS_ND_)
   end subroutine

   subroutine horizontal_weighted_sum_initialize(self,configunit)
      class (type_horizontal_weighted_sum),intent(inout),target :: self
      integer,                             intent(in)           :: configunit

      type (type_horizontal_component),pointer :: component
      integer           :: i
      character(len=10) :: temp

      i = 0
      component => self%first
      do while (associated(component))
         i = i + 1
         write (temp,'(i0)') i
         call self%register_dependency(component%id,'term'//trim(temp))
         call self%request_coupling(component%id,trim(component%name))
         component => component%next
      end do
      if (self%output_long_name=='') self%output_long_name = self%output_name
      call self%register_diagnostic_variable(self%id_output,'result',self%output_units,'result')
      call self%parent%add_alias(self%id_output%link%target,trim(self%output_name))
   end subroutine

   subroutine horizontal_weighted_sum_after_coupling(self)
      class (type_horizontal_weighted_sum),intent(inout) :: self

      type (type_horizontal_component),pointer :: component

      ! At this stage, the background values for all variables (if any) are fixed. We can therefore
      ! compute background contributions already, and add those to the space- and time-invariant offset.
      component => self%first
      do while (associated(component))
         if (component%include_background) self%offset = self%offset + component%weight*component%id%background
         component => component%next
      end do
   end subroutine

   subroutine horizontal_weighted_sum_add_component(self,name,weight,include_background)
      class (type_horizontal_weighted_sum),intent(inout) :: self
      character(len=*),                    intent(in)    :: name
      real(rk),optional,                   intent(in)    :: weight
      logical,optional,                    intent(in)    :: include_background

      type (type_horizontal_component),pointer :: component

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

   subroutine horizontal_weighted_sum_evaluate_horizontal(self,_ARGUMENTS_HZ_)
      class (type_horizontal_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_HZ_
      
      type (type_horizontal_component),pointer        :: component
      real(rk)                                        :: value
      real(rk) _DIMENSION_HORIZONTAL_SLICE_AUTOMATIC_ :: sum

      sum = self%offset

      component => self%first
      do while (associated(component))
         _HORIZONTAL_LOOP_BEGIN_
            _GET_HORIZONTAL_(component%id,value)
            sum _INDEX_HORIZONTAL_SLICE_ = sum _INDEX_HORIZONTAL_SLICE_ + component%weight*value
         _HORIZONTAL_LOOP_END_
         component => component%next
      end do

      _HORIZONTAL_LOOP_BEGIN_
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,sum _INDEX_HORIZONTAL_SLICE_)
      _HORIZONTAL_LOOP_END_
   end subroutine

   subroutine horizontal_weighted_sum_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_horizontal_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      call self%evaluate_horizontal(_ARGUMENTS_HZ_)
   end subroutine

   subroutine base_state_to_conserved_quantities(self,_ARGUMENTS_ND_,y,sums)
      class (type_base_model),          intent(in)  :: self
      _DECLARE_ARGUMENTS_ND_
      real(rk) _DIMENSION_SLICE_PLUS_1_,intent(in)  :: y
      real(rk) _DIMENSION_SLICE_PLUS_1_,intent(out) :: sums

      integer  :: i,j,index
      real(rk) :: scale_factor

      sums = 0.0_rk
      !do i=1,size(self%conserved_quantities)
      !   do j=1,size(self%conserved_quantities(i)%state_indices)
      !      index = self%conserved_quantities(i)%state_indices(j)
      !      scale_factor = self%conserved_quantities(i)%state_scale_factors(j)
      !      _LOOP_BEGIN_
      !         sums _INDEX_SLICE_PLUS_1_(i) = sums _INDEX_SLICE_PLUS_1_(i) + scale_factor*y _INDEX_SLICE_PLUS_1_(index)
      !      _LOOP_END_
      !   end do
      !end do
   end subroutine

   function time_treatment2output(time_treatment) result(output)
      integer, intent(in) :: time_treatment
      integer             :: output
      select case (time_treatment)
         case (time_treatment_last);            output = output_instantaneous
         case (time_treatment_integrated);      output = output_time_integrated
         case (time_treatment_averaged);        output = output_time_step_averaged
         case (time_treatment_step_integrated); output = output_time_step_integrated
      end select
   end function

   function output2time_treatment(output) result(time_treatment)
      integer, intent(in) :: output
      integer             :: time_treatment
      select case (output)
         case (output_none);                 time_treatment = time_treatment_last
         case (output_instantaneous);        time_treatment = time_treatment_last
         case (output_time_integrated);      time_treatment = time_treatment_integrated
         case (output_time_step_averaged);   time_treatment = time_treatment_averaged
         case (output_time_step_integrated); time_treatment = time_treatment_step_integrated
      end select
   end function

   subroutine simple_depth_integral_initialize(self,configunit)
      class (type_simple_depth_integral),intent(inout),target :: self
      integer,                           intent(in)           :: configunit
      call self%register_dependency(self%id_input,'source')
      if (.not.self%average) call self%register_dependency(self%id_depth,standard_variables%bottom_depth)
      call self%register_diagnostic_variable(self%id_output,'result','','result')
   end subroutine

   subroutine simple_depth_integral_do(self,_ARGUMENTS_DO_)
      class (type_simple_depth_integral),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value,depth

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_input,value)
         if (.not.self%average) then
            _GET_HORIZONTAL_(self%id_depth,depth)
            value = value*(min(self%maximum_depth,depth)-self%minimum_depth)
         end if
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,value)
      _HORIZONTAL_LOOP_END_
   end subroutine

   end module fabm_types

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

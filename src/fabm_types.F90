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
   public type_link, type_link_list
   public type_internal_variable
   public type_environment

   public type_model_list,type_model_list_node

   public get_free_unit
   public get_safe_name

   public type_expression, type_bulk_expression, type_horizontal_expression

   public type_bulk_data_pointer,type_horizontal_data_pointer,type_scalar_data_pointer
   public get_aggregate_variable, type_aggregate_variable, type_contributing_variable, type_contribution
   public time_treatment2output,output2time_treatment

! !PUBLIC DATA MEMBERS:
!
   integer, parameter, public :: attribute_length = 256

   integer, parameter, public :: rk = _FABM_REAL_KIND_

   integer, parameter, public :: domain_bulk = 4, domain_horizontal = 8, domain_scalar = 16, domain_bottom = 9, domain_surface = 10

   integer, parameter, public :: source_unknown = 0, source_do = 1, source_do_column = 2, source_do_bottom = 3, source_do_surface = 4

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
      integer                                 :: value = -1
   contains
      procedure :: append      => integer_pointer_set_append
      procedure :: extend      => integer_pointer_set_extend
      procedure :: set_value   => integer_pointer_set_set_value
      procedure :: is_empty    => integer_pointer_set_is_empty
      procedure :: copy_values => integer_pointer_copy_values
      procedure :: clear       => integer_pointer_clear
   end type

   ! ====================================================================================================
   ! Variable identifiers used by biogeochemical models.
   ! ====================================================================================================

   type,abstract :: type_variable_id
      type (type_link), pointer :: link => null()
   end type

   type,extends(type_variable_id) :: type_state_variable_id
      integer  :: index              = -1
      integer  :: state_index        = -1
      integer  :: sms_index          = -1
      integer  :: surface_flux_index = -1
      integer  :: bottom_flux_index  = -1
      real(rk) :: background         = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_bottom_state_variable_id
      integer  :: horizontal_index   = -1
      integer  :: bottom_sms_index   = -1
      integer  :: bottom_state_index = -1
      real(rk) :: background         = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_surface_state_variable_id
      integer  :: horizontal_index    = -1
      integer  :: surface_sms_index   = -1
      integer  :: surface_state_index = -1
      real(rk) :: background          = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_diagnostic_variable_id
      integer :: diag_index = -1
   end type

   type,extends(type_variable_id) :: type_horizontal_diagnostic_variable_id
      integer :: horizontal_diag_index = -1
   end type

   type,extends(type_variable_id) :: type_dependency_id
      integer  :: index      = -1
      real(rk) :: background = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_horizontal_dependency_id
      integer  :: horizontal_index = -1
      real(rk) :: background       = 0.0_rk
   end type

   type,extends(type_variable_id) :: type_global_dependency_id
      integer  :: global_index = -1
      real(rk) :: background   = 0.0_rk
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
      type (type_aggregate_variable),      pointer  :: next => null()
      logical                                       :: bulk_required = .false.
      logical                                       :: horizontal_required = .false.
   end type

   type type_link_list
      type (type_link),pointer :: first => null()
   contains
      procedure :: append   => link_list_append
      procedure :: find     => link_list_find
      procedure :: count    => link_list_count
      procedure :: finalize => link_list_finalize
      procedure :: extend   => link_list_extend
   end type

   ! ====================================================================================================
   ! Derived types used internally to store information on model variables and model references.
   ! ====================================================================================================

   type type_internal_variable
      character(len=attribute_length) :: name      = ''
      character(len=attribute_length) :: long_name = ''
      type (type_property_dictionary) :: properties
      character(len=attribute_length) :: units          = ''
      real(rk)                        :: minimum        = -1.e20_rk
      real(rk)                        :: maximum        =  1.e20_rk
      real(rk)                        :: missing_value  = -2.e20_rk
      real(rk)                        :: initial_value  = 0.0_rk
      integer                         :: output         = output_instantaneous
      integer                         :: presence       = presence_internal
      integer                         :: domain         = domain_bulk
      integer                         :: source         = source_unknown
      class (type_base_model),pointer :: owner          => null()
      type (type_contribution_list)   :: contributions

      class (type_standard_variable), pointer :: standard_variable => null()

      logical :: fake_state_variable = .false.

      ! Only used for bulk variables:
      real(rk) :: vertical_movement         = 0.0_rk
      logical  :: no_precipitation_dilution = .false.
      logical  :: no_river_dilution         = .false.

      integer,pointer :: read_index  => null()
      integer,pointer :: write_index => null()

      ! Collections to collect information from all coupled variables.
      type (type_integer_pointer_set) :: read_indices,state_indices,write_indices
      type (type_real_pointer_set)    :: background_values
      type (type_link_list)           :: sms_list,surface_flux_list,bottom_flux_list
   end type

   type type_link
      character(len=attribute_length)        :: name     = ''
      type (type_internal_variable), pointer :: target   => null()
      type (type_internal_variable), pointer :: original => null()
      type (type_link), pointer              :: next     => null()
   end type

   ! ====================================================================================================
   ! Types that describe custom expressions (arbitrary functions of one or more variables).
   ! ====================================================================================================

   type,abstract :: type_expression
      class (type_expression), pointer :: next        => null()
      character(len=attribute_length)  :: output_name = ''
      integer,pointer :: out => null()
   end type

   type,abstract,extends(type_expression) :: type_bulk_expression
      !type (type_bulk_data_pointer),pointer :: out => null()
   end type

   type,abstract,extends(type_expression) :: type_horizontal_expression
      !type (type_horizontal_data_pointer),pointer :: out => null()
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
      procedure :: append     => model_list_append
      procedure :: extend     => model_list_extend
      procedure :: find_name  => model_list_find_name
      procedure :: find_model => model_list_find_model
      procedure :: count      => model_list_count
      procedure :: finalize   => model_list_finalize
      procedure :: print      => model_list_print
      generic   :: find       => find_name, find_model
   end type

   type type_reused_diagnostic
      integer :: read_index
      integer :: write_index
      integer :: source
   end type

   type type_base_model
      ! Flag determining whether the contents of the type are "frozen", i.e., they will not change anymore.
      logical :: frozen = .false.

      ! Flag determining whether this model was explciitly created by the user (e.g., by requesting it in a namelist or YAML file)
      logical :: user_created = .false.

      ! Pointers to linked models in the model tree.
      class (type_base_model),pointer :: parent => null()
      type (type_model_list)          :: children

      ! Model name and variable prefixes.
      character(len=attribute_length) :: name      = ''
      character(len=attribute_length) :: long_name = ''
      character(len=attribute_length) :: type_name = ''

      ! Models constituents: links to variables, coupling requests, parameters, expressions
      type (type_link_list) :: links
      type (type_aggregate_variable),pointer :: first_aggregate_variable => null()

      type (type_hierarchical_dictionary) :: couplings
      type (type_hierarchical_dictionary) :: parameters

      class (type_expression), pointer :: first_expression => null()

      real(rk) :: dt = 1.0_rk

      logical :: check_conservation = .false.

      type (type_reused_diagnostic),dimension(:),allocatable :: reused_diag,reused_diag_hz
   contains

      ! Procedure for adding child models [during initialization only]
      procedure :: add_child

      ! Procedures for adding variables [during initialization only]
      procedure :: add_bulk_variable
      procedure :: add_horizontal_variable
      procedure :: add_scalar_variable
      procedure :: add_object

      ! Procedures for locating links, objects, models.
      procedure :: find_link
      procedure :: find_object
      procedure :: find_model

      ! Procedures for requesting coupling between variables
      procedure :: request_coupling_for_link
      procedure :: request_coupling_for_name
      procedure :: request_coupling_for_id
      generic   :: request_coupling => request_coupling_for_link,request_coupling_for_name,request_coupling_for_id

      ! Procedures that may be used to query parameter values during initialization.
      procedure :: get_real_parameter
      procedure :: get_integer_parameter
      procedure :: get_logical_parameter
      procedure :: get_string_parameter
      generic :: get_parameter => get_real_parameter,get_integer_parameter,get_logical_parameter,get_string_parameter

      procedure :: set_variable_property_real
      procedure :: set_variable_property_integer
      procedure :: set_variable_property_logical
      generic   :: set_variable_property => set_variable_property_real,set_variable_property_integer,set_variable_property_logical

      procedure :: add_to_aggregate_variable

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
      procedure :: register_named_bulk_dependency_old
      procedure :: register_named_horizontal_dependency_old
      procedure :: register_named_global_dependency_old

      generic :: register_bulk_dependency       => register_named_bulk_dependency, register_standard_bulk_dependency, &
                                                   register_named_bulk_dependency_old
      generic :: register_horizontal_dependency => register_named_horizontal_dependency, register_standard_horizontal_dependency, &
                                                   register_named_horizontal_dependency_old
      generic :: register_global_dependency     => register_named_global_dependency, register_standard_global_dependency, &
                                                   register_named_global_dependency_old

      procedure :: register_standard_conserved_quantity
      procedure :: register_custom_conserved_quantity
      procedure :: register_bulk_state_dependency_ex
      procedure :: register_bottom_state_dependency_ex
      procedure :: register_surface_state_dependency_ex
      procedure :: register_bulk_state_dependency_old
      procedure :: register_bottom_state_dependency_old
      procedure :: register_surface_state_dependency_old

      generic :: register_bulk_state_dependency => register_bulk_state_dependency_ex, register_bulk_state_dependency_old
      generic :: register_bottom_state_dependency => register_bottom_state_dependency_ex, register_bottom_state_dependency_old
      generic :: register_surface_state_dependency => register_surface_state_dependency_ex, register_surface_state_dependency_old

      procedure :: register_bulk_expression_dependency
      procedure :: register_horizontal_expression_dependency
      generic :: register_expression_dependency => register_bulk_expression_dependency,register_horizontal_expression_dependency

      generic :: register_state_variable      => register_bulk_state_variable,register_bottom_state_variable, &
                                                 register_surface_state_variable
      generic :: register_diagnostic_variable => register_bulk_diagnostic_variable,register_horizontal_diagnostic_variable
      generic :: register_dependency          => register_named_bulk_dependency, register_standard_bulk_dependency, &
                                                 register_named_bulk_dependency_old, &
                                                 register_named_horizontal_dependency, register_standard_horizontal_dependency, &
                                                 register_named_horizontal_dependency_old, &
                                                 register_named_global_dependency, register_standard_global_dependency, &
                                                 register_named_global_dependency_old, &
                                                 register_bulk_expression_dependency, register_horizontal_expression_dependency
      generic :: register_state_dependency    => register_bulk_state_dependency_ex,register_bottom_state_dependency_ex, &
                                                 register_surface_state_dependency_ex,register_bulk_state_dependency_old, &
                                                 register_bottom_state_dependency_old,register_surface_state_dependency_old
      generic :: register_conserved_quantity  => register_standard_conserved_quantity, register_custom_conserved_quantity

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
      procedure :: get_horizontal_conserved_quantities => base_get_horizontal_conserved_quantities
      procedure :: check_state              => base_check_state
      procedure :: check_surface_state      => base_check_surface_state
      procedure :: check_bottom_state       => base_check_bottom_state
      procedure :: fatal_error              => base_fatal_error
      procedure :: log_message              => base_log_message
      procedure :: get_path                 => base_get_path

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
      ! Registry with pointers to global fields of readable variables.
      ! These pointers are accessed to fill the prefetch (see below).
      type (type_bulk_data_pointer),      allocatable :: data(:)
      type (type_horizontal_data_pointer),allocatable :: data_hz(:)
      type (type_scalar_data_pointer),    allocatable :: data_scalar(:)

      ! Prefetch arrays that will hold readable data for a single domain slice.
      ! Biogeochemical models use only these to read data.
      real(rk),allocatable _DIMENSION_SLICE_PLUS_1_ALLOCATABLE_            :: prefetch
      real(rk),allocatable _DIMENSION_HORIZONTAL_SLICE_PLUS_1_ALLOCATABLE_ :: prefetch_hz
      real(rk),allocatable,dimension(:)                                    :: prefetch_scalar

      ! Scratch arrays for writing data associated with for a single domain slice.
      real(rk),allocatable _DIMENSION_SLICE_PLUS_1_ALLOCATABLE_            :: scratch
      real(rk),allocatable _DIMENSION_HORIZONTAL_SLICE_PLUS_1_ALLOCATABLE_ :: scratch_hz

#ifdef _FABM_MASK_
      _FABM_MASK_TYPE_,pointer _DIMENSION_GLOBAL_ :: mask => null()
      logical,allocatable _DIMENSION_SLICE_ALLOCATABLE_ :: prefetch_mask
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

!-----------------------------------------------------------------------

   contains

   subroutine base_initialize(self,configunit)
      class (type_base_model),intent(inout),target :: self
      integer,                intent(in)           :: configunit
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

   function base_get_path(self) result(path)
      class (type_base_model), intent(in), target :: self
      character(len=attribute_length) :: path

      class (type_base_model),pointer :: current

      path = ''
      current => self
      do while (associated(current%parent))
         path = '/'//trim(current%name)//trim(path)
         current => current%parent
      end do
   end function

   subroutine base_fatal_error(self,location,message)
      class (type_base_model), intent(in) :: self
      character(len=*),        intent(in) :: location,message
      if (self%name/='') then
         call fatal_error('model '//trim(self%get_path())//', '//trim(location),message)
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
      integer :: islash
      class (type_base_model),pointer :: parent
!EOP
!-----------------------------------------------------------------------
!BOC
      ! If a path with / is given, redirect to tentative parent model.
      islash = index(name,'/',.true.)
      if (islash/=0) then
         parent => self%find_model(name(:islash-1))
         if (.not.associated(parent)) call self%fatal_error('add_child', &
            'Proposed parent model "'//trim(name(:islash-1))//'" was not found.')
         call parent%add_child(model,name(islash+1:),long_name,configunit)
         return
      end if

      if (associated(model%parent)) call self%fatal_error('add_child', &
         'The provided child model "'//trim(name)//'" has already been assigned parent '//trim(model%parent%name)//'.')

      ! Ascertain whether the provided name is valid.
      if (name=='' .or. name/=get_safe_name(name)) call self%fatal_error('add_child', &
         'Cannot add child model "'//trim(name)//'" because its name is not valid. &
         &Model names should not be empty, and can contain letters, digits and underscores only.')
      if (len_trim(name)>len(model%name)) call fatal_error('add_child','Model name "'//trim(name)//'" exceeds maximum length.')

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

      call variable_id%link%target%contributions%add(target,scale_factor,include_background)

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

   function model_list_find_name(self,name) result(node)
      class (type_model_list),intent(in) :: self
      character(len=*),       intent(in) :: name

      type (type_model_list_node),pointer :: node

      node => self%first
      do while (associated(node))
         if (node%model%name==name) return
         node => node%next
      end do
   end function model_list_find_name

   function model_list_find_model(self,model) result(node)
      class (type_model_list),       intent(in) :: self
      class (type_base_model),target,intent(in) :: model

      type (type_model_list_node),pointer :: node

      node => self%first
      do while (associated(node))
         if (associated(node%model,model)) return
         node => node%next
      end do
   end function model_list_find_model

   subroutine model_list_print(self)
      class (type_model_list),       intent(in) :: self

      type (type_model_list_node),pointer :: node

      node => self%first
      do while (associated(node))
         call driver%log_message(node%model%get_path())
         node => node%next
      end do
   end subroutine

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

function link_list_find(self,name) result(link)
   class (type_link_list),intent(in) :: self
   character(len=*),      intent(in) :: name

   type (type_link),pointer :: link

   link => self%first
   do while (associated(link))
      if (link%name==name) return
      link => link%next
   end do
end function link_list_find

function link_list_append(self,target,name) result(link)
   class (type_link_list),             intent(inout) :: self
   type (type_internal_variable),pointer             :: target
   character(len=*),                   intent(in)    :: name

   type (type_link),pointer :: link

   ! Append a new link to the list.
   if (.not.associated(self%first)) then
      allocate(self%first)
      link => self%first
   else
      link => self%first
      do while (associated(link%next))
         link => link%next
      end do
      allocate(link%next)
      link => link%next
   end if

   ! Set link attributes.
   link%name = name
   link%target => target
   link%original => target
end function link_list_append

subroutine link_list_extend(self,source)
   class (type_link_list),intent(inout) :: self
   class (type_link_list),intent(in)    :: source

   type (type_link),pointer :: source_link,link

   source_link => source%first
   do while (associated(source_link))
      link => self%append(source_link%target,source_link%name)
      source_link => source_link%next
   end do
end subroutine link_list_extend

function link_list_count(self) result(count)
   class (type_link_list),intent(in) :: self

   type (type_link),pointer :: link
   integer                  :: count

   count = 0
   link => self%first
   do while (associated(link))
      count = count + 1
      link => link%next
   end do
end function link_list_count

subroutine link_list_finalize(self)
   class (type_link_list),intent(inout) :: self

   type (type_link),pointer :: link,next

   link => self%first
   do while (associated(link))
      next => link%next
      deallocate(link)
      link => next
   end do
   nullify(self%first)
end subroutine link_list_finalize

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
   self%pointers(size(self%pointers))%p = self%value
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

subroutine integer_pointer_clear(self)
   class (type_integer_pointer_set),intent(inout) :: self

   if (allocated(self%pointers)) deallocate(self%pointers)
   self%value = -1
end subroutine integer_pointer_clear

subroutine integer_pointer_set_set_value(self,value)
   class (type_integer_pointer_set),intent(inout) :: self
   integer,                         intent(in)    :: value

   integer :: i

   do i=1,size(self%pointers)
      self%pointers(i)%p = value
   end do
   self%value = value
end subroutine integer_pointer_set_set_value

logical function integer_pointer_set_is_empty(self)
   class (type_integer_pointer_set),intent(in) :: self
   integer_pointer_set_is_empty = .not.allocated(self%pointers)
end function integer_pointer_set_is_empty

subroutine integer_pointer_copy_values(self,target)
   class (type_integer_pointer_set),intent(in) :: self
   integer,allocatable                         :: target(:)

   integer :: i

   if (.not.allocated(self%pointers)) then
      allocate(target(0))
   else
      allocate(target(size(self%pointers)))
      do i=1,size(self%pointers)
         target(i) = self%pointers(i)%p
      end do
   end if
end subroutine

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
                                           standard_variable, presence, background_value)
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
!EOP
!-----------------------------------------------------------------------
!BOC
      if (associated(id%link)) call self%fatal_error('register_bulk_state_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      call self%add_bulk_variable(name, units, long_name, missing_value, minimum, maximum, &
                                  initial_value=initial_value, background_value=background_value, &
                                  vertical_movement=vertical_movement, specific_light_extinction=specific_light_extinction, &
                                  no_precipitation_dilution=no_precipitation_dilution, no_river_dilution=no_river_dilution, &
                                  standard_variable=standard_variable, presence=presence, &
                                  state_index=id%state_index, read_index=id%index, sms_index=id%sms_index, &
                                  surface_flux_index=id%surface_flux_index, bottom_flux_index=id%bottom_flux_index, &
                                  background=id%background, link=id%link)

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
                                             standard_variable, presence, background_value)
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
!EOP
!-----------------------------------------------------------------------
!BOC
      if (associated(id%link)) call self%fatal_error('register_bottom_state_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      call self%add_horizontal_variable(name, units, long_name, missing_value, minimum, maximum, &
                                        initial_value=initial_value, background_value=background_value, &
                                        standard_variable=standard_variable, presence=presence, domain=domain_bottom, &
                                        state_index=id%bottom_state_index, read_index=id%horizontal_index, &
                                        sms_index=id%bottom_sms_index, background=id%background, link=id%link)

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
                                              standard_variable, presence, background_value)
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
!EOP
!-----------------------------------------------------------------------
!BOC
      if (associated(id%link)) call self%fatal_error('register_surface_state_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      call self%add_horizontal_variable(name, units, long_name, missing_value, minimum, maximum, &
                                        initial_value=initial_value, background_value=background_value, &
                                        standard_variable=standard_variable, presence=presence, domain=domain_surface, &
                                        state_index=id%surface_state_index, read_index=id%horizontal_index, &
                                        sms_index=id%surface_sms_index, background=id%background, link=id%link)

   end subroutine register_surface_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adds a generic variable (bulk, horizontal or scalar) to the model and returns the link to it.
!
! !INTERFACE:
   subroutine add_variable(self, variable, name, units, long_name, missing_value, minimum, maximum, &
                           initial_value, background_value, presence, output, time_treatment, &
                           act_as_state_variable, read_index, state_index, write_index, background, link)
!
! !DESCRIPTION:
!  This function fills all generic variable fields (i.e., those independent of the variable's domain), and returns
!  a link to the variable. A pointer to the provided variable object will be kept after the function returns. Thus,
!  the object must not be deallocated (create it by allocating a pointer).
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),       target,intent(inout)       :: self
      type (type_internal_variable),pointer                     :: variable
      character(len=*),              target,intent(in)          :: name
      character(len=*),                     intent(in),optional :: long_name, units
      real(rk),                             intent(in),optional :: minimum, maximum,missing_value,initial_value,background_value
      integer,                              intent(in),optional :: presence, output, time_treatment
      logical,                              intent(in),optional :: act_as_state_variable
      integer,                       target,           optional :: read_index, state_index, write_index
      real(rk),                      target,           optional :: background
      type (type_link),              pointer,          optional :: link
!
!EOP
!
! !LOCAL VARIABLES:
      character(len=256)                     :: text
      type (type_link),              pointer :: link_
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('add_variable', &
         'Cannot register variable "'//trim(name)//'" because the model initialization phase has already completed &
         &(fabm_initialize has been called).')

      if (len_trim(name)>len(variable%name)) call self%fatal_error('add_variable', &
         'Variable name "'//trim(name)//'" exceeds maximum length.')

      ! Ascertain whether the provided name is valid.
      if (name=='' .or. name/=get_safe_name(name)) call self%fatal_error('add_variable', &
         'Cannot register variable "'//trim(name)//'" because its name is not valid. &
         &Variable names should not be empty, and can contain letters, digits and underscores only.')

      variable%name = name
      variable%owner => self
      if (present(units))         variable%units         = units
      if (present(long_name)) then
         variable%long_name = long_name
      else
         variable%long_name = name
      end if
      if (present(minimum))       variable%minimum       = minimum
      if (present(maximum))       variable%maximum       = maximum
      if (present(missing_value)) variable%missing_value = missing_value
      if (present(initial_value)) variable%initial_value = initial_value
      if (present(presence))      variable%presence      = presence
      if (present(act_as_state_variable)) variable%fake_state_variable = act_as_state_variable
      if (present(time_treatment)) then
         call self%log_message('variable "'//trim(name)//'": "time_treatment" argument is deprecated; &
                               &please use "output" instead. For possible values, see time_treatment2output in fabm_types.F90.')
         variable%output = time_treatment2output(time_treatment)
      end if
      if (present(output))        variable%output        = output

      if (present(state_index)) then
         ! Ensure that initial value falls within prescribed valid range.
         if (variable%initial_value<variable%minimum .or. variable%initial_value>variable%maximum) then
            write (text,*) 'Initial value',variable%initial_value,'for variable "'//trim(name)//'" lies&
                  &outside allowed range',variable%minimum,'to',variable%maximum
            call self%fatal_error('fill_internal_variable',text)
         end if

         ! Store a pointer to the variable that should hold the state variable index.
         call variable%state_indices%append(state_index)
      end if


      if (present(background)) then
         ! Store a pointer to the variable that should hold the background value.
         ! If the background value itself is also prescribed, use it.
         call variable%background_values%append(background)
         if (present(background_value)) call variable%background_values%set_value(background_value)
      end if

      if (present(read_index)) then
         variable%read_index => read_index
         call variable%read_indices%append(read_index)
      end if
      if (present(write_index)) then
         variable%write_index => write_index
         call variable%write_indices%append(write_index)
      end if

      ! Create a class pointer and use that to create a link.
      link_ => add_object(self,variable)
      if (present(link)) link => link_

   end subroutine add_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adds a bulk variable to the model and returns the link to it.
!
! !INTERFACE:
   recursive subroutine add_bulk_variable(self, name, units, long_name, missing_value, minimum, maximum, initial_value, &
                                          background_value, vertical_movement, specific_light_extinction, &
                                          no_precipitation_dilution, no_river_dilution, standard_variable, presence, output, &
                                          time_treatment, act_as_state_variable, source, &
                                          read_index, state_index, write_index, sms_index, surface_flux_index, bottom_flux_index, &
                                          background, link)
!
! !DESCRIPTION:
!  This function registers a new bulk variable. It is not predefined to be a state variable, diagnostic variable or dependency.
!  This is set implicitly by providing (or omitting) target variables state_index, write_index, data, background.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),target,    intent(inout)       :: self
      character(len=*),                  intent(in)          :: name
      character(len=*),                  intent(in),optional :: long_name, units
      real(rk),                          intent(in),optional :: minimum, maximum, missing_value, initial_value, background_value
      real(rk),                          intent(in),optional :: vertical_movement, specific_light_extinction
      logical,                           intent(in),optional :: no_precipitation_dilution, no_river_dilution
      type (type_bulk_standard_variable),intent(in),optional :: standard_variable
      integer,                           intent(in),optional :: presence, output, time_treatment, source
      logical,                           intent(in),optional :: act_as_state_variable

      integer,                      target,optional :: read_index, state_index, write_index
      integer,                      target,optional :: sms_index, surface_flux_index, bottom_flux_index
      real(rk),                     target,optional :: background

      type (type_link),pointer,optional :: link
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_internal_variable), pointer :: variable
      type (type_link),              pointer :: link_,link2,link_dum
!
!-----------------------------------------------------------------------
!BOC
      allocate(variable)
      variable%domain = domain_bulk

      ! Fill fields specific to bulk variables.
      variable%source = source_do
      if (present(source))                    variable%source                    = source
      if (present(vertical_movement))         variable%vertical_movement         = vertical_movement
      if (present(no_precipitation_dilution)) variable%no_precipitation_dilution = no_precipitation_dilution
      if (present(no_river_dilution))         variable%no_river_dilution         = no_river_dilution
      if (present(standard_variable))         allocate(variable%standard_variable,source=standard_variable)
      if (present(specific_light_extinction)) call variable%contributions%add( &
         standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux,scale_factor=specific_light_extinction)

      ! Process remainder of fields and creation of link generically (i.e., irrespective of variable domain).
      call add_variable(self, variable, name, units, long_name, missing_value, minimum, maximum, &
                        initial_value, background_value, presence, output, time_treatment, &
                        act_as_state_variable, read_index, state_index, write_index, background, link_)

      if (present(sms_index)) then
         call self%add_bulk_variable(trim(link_%name)//'_sms', trim(units)//'/s', trim(long_name)//' sources-sinks', &
                                     0.0_rk, output=output_none, write_index=sms_index, link=link2)
         link_dum => variable%sms_list%append(link2%target,link2%target%name)
      end if
      if (present(surface_flux_index)) then
         call self%add_horizontal_variable(trim(link_%name)//'_sfl', trim(units)//'*m/s', trim(long_name)//' surface flux', &
                                     0.0_rk, output=output_none, write_index=surface_flux_index, link=link2, &
                                     domain=domain_surface, source=source_do_surface)
         link_dum => variable%surface_flux_list%append(link2%target,link2%target%name)
      end if
      if (present(bottom_flux_index)) then
         call self%add_horizontal_variable(trim(link_%name)//'_bfl', trim(units)//'*m/s', trim(long_name)//' bottom flux', &
                                     0.0_rk, output=output_none, write_index=bottom_flux_index, link=link2, &
                                     domain=domain_bottom, source=source_do_bottom)
         link_dum => variable%bottom_flux_list%append(link2%target,link2%target%name)
      end if

      if (present(link)) link => link_
   end subroutine add_bulk_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adds a horizontal variable to the model and returns the link to it.
!
! !INTERFACE:
   recursive subroutine add_horizontal_variable(self,name,units,long_name, missing_value, minimum, maximum, initial_value, &
                                                background_value, standard_variable, presence, output, time_treatment, &
                                                act_as_state_variable, domain, source, &
                                                read_index, state_index, write_index, sms_index, background, link)
!
! !DESCRIPTION:
!  This function registers a new horizontal variable. It is not predefined to be a state variable, diagnostic variable
!  or dependency. This is set implicitly by providing (or omitting) target variables state_index, write_index, data, background.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),target,           intent(inout)       :: self
      character(len=*),                         intent(in)          :: name
      character(len=*),                         intent(in),optional :: long_name, units
      real(rk),                                 intent(in),optional :: minimum, maximum, missing_value
      real(rk),                                 intent(in),optional :: initial_value, background_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
      integer,                                  intent(in),optional :: presence, domain, output, time_treatment, source
      logical,                                  intent(in),optional :: act_as_state_variable

      integer,                            target,optional :: read_index, state_index, write_index, sms_index
      real(rk),                           target,optional :: background

      type (type_link),pointer,optional :: link
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_internal_variable),pointer :: variable
      type (type_link),             pointer :: link_,link2,link_dum
      integer                               :: sms_source
!
!-----------------------------------------------------------------------
!BOC
      allocate(variable)
      variable%domain = domain_horizontal
      variable%source = source_unknown

      ! Fill fields specific to horizontal variables.
      if (present(domain))            variable%domain = domain
      if (present(standard_variable)) allocate(variable%standard_variable,source=standard_variable)
      if (present(source))            variable%source = source

      ! Process remainder of fields and creation of link generically (i.e., irrespective of variable domain).
      call add_variable(self, variable, name, units, long_name, missing_value, minimum, maximum, &
                        initial_value, background_value, presence, output, time_treatment, &
                        act_as_state_variable, read_index, state_index, write_index, background, link_)

      if (present(sms_index)) then
         sms_source = source_do_bottom
         if (variable%domain==domain_surface) sms_source = source_do_surface
         call self%add_horizontal_variable(trim(link_%name)//'_sms', trim(units)//'/s', trim(long_name)//' sources-sinks', &
                                           0.0_rk, output=output_none, write_index=sms_index, link=link2, &
                                           domain=variable%domain, source=sms_source)
         link_dum => variable%sms_list%append(link2%target,link2%target%name)
      end if

      if (present(link)) link => link_
   end subroutine add_horizontal_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adds a scalar variable to the model and returns the link to it.
!
! !INTERFACE:
   recursive subroutine add_scalar_variable(self,name, long_name, units, missing_value, minimum, maximum, initial_value, &
                                            background_value, standard_variable, presence, output, time_treatment, &
                                            read_index, state_index, write_index, sms_index, background, link)
!
! !DESCRIPTION:
!  This function registers a new scalar variable. It is not predefined to be a state variable, diagnostic variable or dependency.
!  This is set implicitly by providing (or omitting) target variables state_index, write_index, data, background.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),target,       intent(inout)       :: self
      character(len=*),                     intent(in)          :: name
      character(len=*),                     intent(in),optional :: long_name, units
      real(rk),                             intent(in),optional :: minimum, maximum, missing_value, initial_value, background_value
      type (type_global_standard_variable), intent(in),optional :: standard_variable
      integer,                              intent(in),optional :: presence, output, time_treatment

      integer,                        target,optional :: read_index, state_index, write_index, sms_index
      real(rk),                       target,optional :: background

      type (type_link),pointer,optional :: link
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_internal_variable), pointer :: variable
!
!-----------------------------------------------------------------------
!BOC
      allocate(variable)
      variable%domain = domain_scalar

      ! Fill fields specific to scalar variables.
      if (present(standard_variable)) allocate(variable%standard_variable,source=standard_variable)

      ! Process remainder of fields and creation of link generically (i.e., irrespective of variable domain).
      call add_variable(self, variable, name, units, long_name, missing_value, minimum, maximum, &
                        initial_value, background_value, presence, output, time_treatment, &
                        .false., read_index, state_index, write_index, background, link)
   end subroutine add_scalar_variable
!EOC

   recursive function add_object(self,object) result(link)
      ! This subroutine creates a link to the supplied object, then allows
      ! parent models to do the same.
      ! NB this subroutine MUST be recursive, to allow parent models to override
      ! the properties of objects added by their child models.
      class (type_base_model),target,     intent(inout) :: self
      type (type_internal_variable),pointer             :: object

      type (type_link), pointer       :: link,parent_link
      character(len=attribute_length) :: oriname
      integer                         :: instance

      ! First check if a link with this name exists.
      link => self%links%find(object%name)

      if (associated(link)) then
         ! Link with this name exists already.
         ! Append numbers to the variable name until a unique name is found.
         oriname = object%name
         instance = 0
         do
            write (object%name,'(a,a,i0)') trim(oriname),'_',instance
            if (.not.associated(self%links%find(object%name))) exit
            instance = instance + 1
         end do
         call self%request_coupling(object%name,oriname)
      end if

      ! Create link for this object.
      link => self%links%append(object,object%name)

      ! Forward to parent
      if (associated(self%parent)) then
         if (len_trim(self%name)+1+len_trim(object%name)>len(object%name)) call fatal_error('add_object', &
            'Variable path "'//trim(self%name)//'/'//trim(object%name)//'" exceeds maximum allowed length.')
         object%name = trim(self%name)//'/'//trim(object%name)
         parent_link => self%parent%add_object(object)
      end if
   end function add_object

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   subroutine register_bulk_diagnostic_variable(self, id, name, units, long_name, &
                                                time_treatment, missing_value, standard_variable, output, source, &
                                                act_as_state_variable)
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
      integer,                            intent(in),optional :: time_treatment, output, source
      real(rk),                           intent(in),optional :: missing_value
      type (type_bulk_standard_variable), intent(in),optional :: standard_variable
      logical,                            intent(in),optional :: act_as_state_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      if (associated(id%link)) call self%fatal_error('register_bulk_diagnostic_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      call self%add_bulk_variable(name, units, long_name, missing_value, &
                                  standard_variable=standard_variable, output=output, time_treatment=time_treatment, &
                                  source=source, write_index=id%diag_index, link=id%link, act_as_state_variable=act_as_state_variable)

   end subroutine register_bulk_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   subroutine register_horizontal_diagnostic_variable(self, id, name, units, long_name, &
                                                      time_treatment, missing_value, standard_variable, output, source, &
                                                      act_as_state_variable, domain)
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
      integer,                                  intent(in),optional :: time_treatment, output, source, domain
      real(rk),                                 intent(in),optional :: missing_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
      logical,                                  intent(in),optional :: act_as_state_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
      if (associated(id%link)) call self%fatal_error('register_horizontal_diagnostic_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      call self%add_horizontal_variable(name, units, long_name, missing_value, &
                                        standard_variable=standard_variable, output=output, time_treatment=time_treatment, &
                                        source=source, write_index=id%horizontal_diag_index, link=id%link, &
                                        act_as_state_variable=act_as_state_variable, domain=domain)

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
      call self%fatal_error('register_standard_conserved_quantity', &
                            'register_conserved_quantity is no longer supported; please use add_to_aggregate_variable.')

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
      call self%fatal_error('register_standard_conserved_quantity', &
                            'register_conserved_quantity is no longer supported; please use add_to_aggregate_variable.')

   end subroutine register_custom_conserved_quantity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external state variable
!
! !INTERFACE:
   subroutine register_bulk_state_dependency_ex(self,id,name,units,long_name,required,standard_variable)
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
      character(len=*),                   intent(in)          :: name,units,long_name
      logical,                            intent(in),optional :: required
      type (type_bulk_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
      integer :: presence
!-----------------------------------------------------------------------
!BOC
      presence = presence_external_required
      if (present(required)) then
         if (.not.required) presence = presence_external_optional
      end if
      call register_bulk_state_variable(self, id, name, units, long_name, presence=presence, standard_variable=standard_variable)

   end subroutine register_bulk_state_dependency_ex
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external bottom state variable
!
! !INTERFACE:
   subroutine register_bottom_state_dependency_ex(model,id,name,units,long_name,required,standard_variable)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),              intent(inout)        :: model
      type (type_bottom_state_variable_id), intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                         intent(in)          :: name,units,long_name
      logical,                                  intent(in),optional :: required
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
      integer :: presence
!-----------------------------------------------------------------------
!BOC
      presence = presence_external_required
      if (present(required)) then
         if (.not.required) presence = presence_external_optional
      end if
      call register_bottom_state_variable(model, id, name, units, long_name, presence=presence, standard_variable=standard_variable)

   end subroutine register_bottom_state_dependency_ex
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external surface-bound state variable
!
! !INTERFACE:
   subroutine register_surface_state_dependency_ex(model,id,name,units,long_name,required,standard_variable)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model),               intent(inout)        :: model
      type (type_surface_state_variable_id), intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),                         intent(in)          :: name,units,long_name
      logical,                                  intent(in),optional :: required
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
      integer :: presence
!-----------------------------------------------------------------------
!BOC
      presence = presence_external_required
      if (present(required)) then
         if (.not.required) presence = presence_external_optional
      end if
      call register_surface_state_variable(model, id, name, units, long_name, &
                                           presence=presence, standard_variable=standard_variable)

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

      call register_named_bulk_dependency(model,id,standard_variable%name,standard_variable%units,standard_variable%name, &
                                          standard_variable=standard_variable,required=required)
   end subroutine

   subroutine register_standard_horizontal_dependency(model,id,standard_variable,required)
      class (type_base_model),                 intent(inout)        :: model
      type (type_horizontal_dependency_id),    intent(inout),target :: id
      type (type_horizontal_standard_variable),intent(in)           :: standard_variable
      logical,optional,                        intent(in)           :: required

      call register_named_horizontal_dependency(model,id,standard_variable%name,standard_variable%units,standard_variable%name, &
                                                standard_variable=standard_variable,required=required)
   end subroutine

   subroutine register_standard_global_dependency(model,id,standard_variable,required)
      class (type_base_model),             intent(inout)        :: model
      type (type_global_dependency_id),    intent(inout),target :: id
      type (type_global_standard_variable),intent(in)           :: standard_variable
      logical,optional,                    intent(in)           :: required

      call register_named_global_dependency(model,id,standard_variable%name,standard_variable%units,standard_variable%name, &
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
      character(len=*),                  intent(in)          :: units,long_name
      type (type_bulk_standard_variable),intent(in),optional :: standard_variable
      logical,                           intent(in),optional :: required
!
!EOP
!
! !LOCAL VARIABLES:
      integer :: presence
!
!-----------------------------------------------------------------------
!BOC
      if (associated(id%link)) call self%fatal_error('register_named_bulk_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      ! Dependencies MUST be fulfilled, unless explicitly specified that this is not so (required=.false.)
      presence = presence_external_required
      if (present(required)) then
         if (.not.required) presence = presence_external_optional
      end if

      call self%add_bulk_variable(name, units, long_name, &
                                  standard_variable=standard_variable, presence=presence, &
                                  read_index=id%index, background=id%background, link=id%link)

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
      character(len=*),                        intent(in)          :: units,long_name
      type (type_horizontal_standard_variable),intent(in),optional :: standard_variable
      logical,                                 intent(in),optional :: required
!
!EOP
!
! !LOCAL VARIABLES:
      integer :: presence
!
!-----------------------------------------------------------------------
!BOC
      if (associated(id%link)) call self%fatal_error('register_named_horizontal_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      ! Dependencies MUST be fulfilled, unless explicitly specified that this is not so (required=.false.)
      presence = presence_external_required
      if (present(required)) then
         if (.not.required) presence = presence_external_optional
      end if

      call self%add_horizontal_variable(name, units, long_name,&
                                        standard_variable=standard_variable, presence=presence, &
                                        read_index=id%horizontal_index, background=id%background, link=id%link)

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
      character(len=*),                    intent(in)          :: units,long_name
      type (type_global_standard_variable),intent(in),optional :: standard_variable
      logical,                             intent(in),optional :: required
!
!EOP
!
! !LOCAL VARIABLES:
      integer :: presence
!
!-----------------------------------------------------------------------
!BOC
      if (associated(id%link)) call self%fatal_error('register_named_global_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      ! Dependencies MUST be fulfilled, unless explicitly specified that this is not so (required=.false.)
      presence = presence_external_required
      if (present(required)) then
         if (.not.required) presence = presence_external_optional
      end if

      call self%add_scalar_variable(name, units, long_name, &
                                    standard_variable=standard_variable, presence=presence, &
                                    read_index=id%global_index, background=id%background, link=id%link)

   end subroutine register_named_global_dependency
!EOC

   subroutine register_named_bulk_dependency_old(self,id,name)
      class (type_base_model),  intent(inout)        :: self
      type (type_dependency_id),intent(inout),target :: id
      character(len=*),         intent(in)           :: name

      call self%log_message('Deprecated syntax for register_bulk_dependency; please call it with a local name, &
                             units, long_name. Subsequently call request_coupling if coupling to an external variable is desired.')
      call self%register_dependency(id,name, '', name)
      if (associated(self%parent)) call self%request_coupling(id,name)
   end subroutine register_named_bulk_dependency_old

   subroutine register_named_horizontal_dependency_old(self,id,name)
      class (type_base_model),             intent(inout)        :: self
      type (type_horizontal_dependency_id),intent(inout),target :: id
      character(len=*),                    intent(in)           :: name

      call self%log_message('Deprecated syntax for register_horizontal_dependency; please call it with a local name, &
                             units, long_name. Subsequently call request_coupling if coupling to an external variable is desired.')
      call self%register_dependency(id,name, '', name)
      if (associated(self%parent)) call self%request_coupling(id,name)
   end subroutine register_named_horizontal_dependency_old

   subroutine register_named_global_dependency_old(self,id,name)
      class (type_base_model),         intent(inout)        :: self
      type (type_global_dependency_id),intent(inout),target :: id
      character(len=*),                intent(in)           :: name

      call self%log_message('Deprecated syntax for register_global_dependency; please call it with a local name, &
                             units, long_name. Subsequently call request_coupling if coupling to an external variable is desired.')
      call self%register_dependency(id,name, '', name)
      if (associated(self%parent)) call self%request_coupling(id,name)
   end subroutine register_named_global_dependency_old

subroutine register_bulk_expression_dependency(self,id,expression)
   class (type_base_model),       intent(inout)        :: self
   type (type_dependency_id),     intent(inout),target :: id
   class (type_bulk_expression),  intent(in)           :: expression

   class (type_bulk_expression),allocatable :: copy

   allocate(copy,source=expression)
   copy%out => id%index
   call self%register_dependency(id,copy%output_name,'',copy%output_name)
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
   copy%out => id%horizontal_index
   call self%register_dependency(id,copy%output_name,'',copy%output_name)
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

   function find_object(self,name,recursive,exact) result(object)

      class (type_base_model),  intent(in),target :: self
      character(len=*),         intent(in)        :: name
      logical,         optional,intent(in)        :: recursive,exact
      type (type_internal_variable),pointer       :: object

      type (type_link), pointer :: link

      nullify(object)
      link => self%find_link(name,recursive,exact)
      if (associated(link)) object => link%target

   end function find_object

   recursive function find_link(self,name,recursive,exact) result(link)

      class (type_base_model),  intent(in),target :: self
      character(len=*),         intent(in)        :: name
      logical,         optional,intent(in)        :: recursive,exact
      type (type_link),pointer                    :: link

      logical                         :: recursive_eff,exact_eff
      class (type_base_model),pointer :: current

      if (name(1:1)=='/') then
         link => find_link(self,name(2:),recursive,exact)
         return
      end if

      recursive_eff = .false.
      if (present(recursive)) recursive_eff = recursive
      nullify(link)

      ! First search self and ancestors (if allowed) based on exact name provided.
      current => self
      do while (associated(current))
         link => current%links%find(name)
         if (associated(link)) return
         if (.not.recursive_eff) exit
         current => current%parent
      end do

      exact_eff = .true.
      if (present(exact)) exact_eff = exact
      if (exact_eff) return

      ! Not found. Now search self and ancestors (if allowed) based on safe name (letters and underscores only).
      current => self
      do while (associated(current))
         link => current%links%first
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
   function find_model(self,name,recursive) result(found_model)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_base_model),       intent(in),target :: self
      character(len=*),              intent(in)        :: name
      logical,optional,              intent(in)        :: recursive
      class (type_base_model),pointer                  :: found_model

      class (type_base_model),pointer     :: current_root
      logical                             :: recursive_eff
      type (type_model_list_node),pointer :: node
      integer                             :: istart,length
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
      nullify(found_model)

      ! Determine whether to also try among ancestors
      recursive_eff = .false.
      if (present(recursive)) recursive_eff = recursive

      current_root => self
      do while (associated(current_root))
         ! Process individual path components (separated by /)
         found_model => current_root
         istart = 1
         do while (associated(found_model).and.istart<=len(name))
            length = index(name(istart:),'/')-1
            if (length==-1) length = len(name) - istart + 1
            node => found_model%children%find(name(istart:istart+length-1))
            istart = istart+length+1
            nullify(found_model)
            if (associated(node)) found_model => node%model
         end do

         ! Only continue if we have not found the model and are allowed to try parent model.
         if (associated(found_model).or..not.recursive_eff) return

         current_root => current_root%parent
      end do
   end function find_model
!EOC

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
   if (.not.associated(self%parent).and.standard_variable%conserved) then
      aggregate_variable%bulk_required = .true.
      aggregate_variable%horizontal_required = .true.
   end if

end function

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

   end module fabm_types

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

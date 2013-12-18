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

   ! Variable identifier types used by biogeochemical models
   public type_diagnostic_variable_id
   public type_horizontal_diagnostic_variable_id
   public type_state_variable_id
   public type_bottom_state_variable_id
   public type_dependency_id
   public type_horizontal_dependency_id
   public type_global_dependency_id
   public type_conserved_quantity_id
   public type_model_id

   ! Variable registration procedures used by biogeochemical models.
   public register_state_variable
   public register_diagnostic_variable
   public register_conserved_quantity
   public register_state_dependency
   public register_dependency

   ! Variable identifier types by external physical drivers.
   public type_bulk_variable_id
   public type_horizontal_variable_id
   public type_scalar_variable_id
   public is_null_standard_variable

   ! Data types and procedures for variable management - used by FABM internally only.
   public type_link
   public type_internal_object,type_internal_variable,type_bulk_variable,type_horizontal_variable,type_scalar_variable
   public type_environment
   public freeze_model_info
   public find_dependencies
   public create_external_variable_id

   public type_bulk_data_pointer
   public type_horizontal_data_pointer

   public type_model_list,type_model_list_node

   public get_free_unit
   public get_safe_name

   public type_expression, type_bulk_expression, type_horizontal_expression

   public type_weighted_sum

! !PUBLIC DATA MEMBERS:
!
   integer, parameter, public :: attribute_length = 256

   integer, parameter, public :: rk = _FABM_REAL_KIND_

   integer, parameter, public :: domain_bulk = 0, domain_bottom = 1, domain_surface = 2

   integer, parameter, public :: presence_internal = 1, presence_external_required = 2, presence_external_optional = 6

   ! Alternative names for standard variables (for backward compatibility)
   type (type_bulk_standard_variable),parameter,public :: &
     varname_temp     = standard_variables%temperature,                   & ! In-situ temperature (degree_Celsius)
     varname_salt     = standard_variables%practical_salinity,                               & ! Salinity on Practical Salinity Scale (1e-3)
     varname_swr      = standard_variables%downwelling_shortwave_flux,                       & ! Shortwave [200-4000 nm] radiation (W m-2)
     varname_par      = standard_variables%downwelling_photosynthetic_radiative_flux,        & ! Photosynthetically Active [400-700 nm] Radiation (W m-2)
     varname_pres     = standard_variables%pressure,                                         & ! Pressure (dbar = 10 kPa)
     varname_dens     = standard_variables%density,                                          & ! In-situ density (kg m-3)
     varname_layer_ht = standard_variables%cell_thickness,                                   & ! Layer thickness (m)
     varname_extc     = standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, & ! Attenuation coefficient for Photosynthetically Active [400-700 nm] Radiation (m-1)
     varname_tss      = standard_variables%mass_concentration_of_suspended_matter              ! Total suspended matter or suspended solids (g m-3)

   ! Variables defined on a horizontal surface (e.g., water surface or bottom).
   type (type_horizontal_standard_variable),parameter,public :: &
     varname_lon     = standard_variables%longitude,                                         & ! Longitude (degree_East)
     varname_lat     = standard_variables%latitude,                                          & ! Latitude (degree_North)
     varname_wind_sf = standard_variables%wind_speed,                                        & ! Wind speed, defined at 10 m above water surface (m s-1)
     varname_cloud   = standard_variables%cloud_area_fraction,                               & ! Cloud cover (1), i.e., a fraction between 0 and 1
     varname_swr_sf  = standard_variables%surface_downwelling_shortwave_flux,                & ! Shortwave [200-4000 nm] radiation, defined at water surface (W m-2)
     varname_par_sf  = standard_variables%surface_downwelling_photosynthetic_radiative_flux, & ! Photosynthetically Active [400-700 nm] Radiation, defined at water surface (W m-2)
     varname_zbot    = standard_variables%bottom_depth_below_geoid,                          & ! Basin floor depth below geoid (approx. mean sea level) (m)
     varname_taub    = standard_variables%bottom_stress                                        ! Bottom stress (Pa)

   ! Non-spatial (scalar) variables.
   type (type_global_standard_variable),parameter,public :: &
     varname_yearday = standard_variables%number_of_days_since_start_of_the_year              ! Decimal day of the year (day), equal to 0.0 at 00:00 1 Jan UTC
!
! !PUBLIC TYPES:
!
   integer, parameter, public :: time_treatment_last=0,time_treatment_integrated=1, &
                                 time_treatment_averaged=2,time_treatment_step_integrated=3

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

   type,abstract :: type_id
      type (type_link), pointer :: link => null()
   end type

   type,extends(type_id) :: type_state_variable_id
      integer                         :: state_index = -1
      type (type_bulk_data_pointer)   :: data
   end type

   type,extends(type_id) :: type_bottom_state_variable_id
      integer                             :: bottom_state_index = -1
      type (type_horizontal_data_pointer) :: horizontal_data
   end type

   type,extends(type_id) :: type_surface_state_variable_id
      integer                             :: surface_state_index = -1
      type (type_horizontal_data_pointer) :: horizontal_data
   end type

   type,extends(type_id) :: type_diagnostic_variable_id
      integer :: diag_index = -1
   end type

   type,extends(type_id) :: type_horizontal_diagnostic_variable_id
      integer :: horizontal_diag_index = -1
   end type

   type,extends(type_id) :: type_dependency_id
      type (type_bulk_data_pointer) :: data
   end type

   type,extends(type_id) :: type_horizontal_dependency_id
      type (type_horizontal_data_pointer) :: horizontal_data
   end type

   type,extends(type_id) :: type_global_dependency_id
      type (type_scalar_data_pointer) :: global_data
   end type

   type,extends(type_id) :: type_conserved_quantity_id
      integer :: cons_index = -1
   end type

   ! ====================================================================================================
   ! Variable types used by FABM for both metadata and value pointers/indices.
   ! ====================================================================================================

   type type_contribution
      type (type_bulk_standard_variable)  :: target
      real(rk)                            :: scale_factor = 1.0_rk
      type (type_contribution),pointer    :: next => null()
   end type

   type type_contribution_list
      type (type_contribution),pointer :: first => null()
   contains
      procedure :: add => contribution_list_add
   end type

   type type_contributing_variable
      type (type_link),                 pointer :: link
      real(rk)                                  :: scale_factor = 1.0_rk
      integer                                   :: state_index = -1
      integer                                   :: horizontal_state_index = -1
      type (type_contributing_variable),pointer :: next => null()
   end type

   type type_aggregate_variable
      type (type_bulk_standard_variable)          :: standard_variable
      type (type_contributing_variable),  pointer :: first_contributing_variable => null()
      type (type_weighted_sum),           pointer :: sum => null()
      type (type_horizontal_weighted_sum),pointer :: horizontal_sum => null()
      type (type_aggregate_variable),     pointer :: next => null()
      type (type_diagnostic_variable_id)            :: id_rate
      type (type_horizontal_diagnostic_variable_id) :: id_horizontal_rate
      logical                                       :: has_bulk_state_component = .false.
      logical                                       :: has_horizontal_state_component = .false.
   end type

   type,abstract :: type_internal_object
      character(len=attribute_length) :: name      = ''
      character(len=attribute_length) :: long_name = ''
      type (type_property_dictionary) :: properties
   end type
   
   type,extends(type_internal_object),abstract :: type_internal_variable
      character(len=attribute_length) :: units         = ''
      real(rk)                        :: minimum       = -1.e20_rk
      real(rk)                        :: maximum       =  1.e20_rk
      real(rk)                        :: missing_value = -2.e20_rk
      real(rk)                        :: initial_value = 0.0_rk
      integer                         :: presence      = presence_internal
      class (type_base_model),pointer :: source_model  => null()
      type (type_contribution_list)   :: contributions

      type (type_integer_pointer),allocatable :: state_indices(:)
      type (type_integer_pointer),allocatable :: write_indices(:)
      type (type_integer_pointer),allocatable :: cons_indices(:)
   end type

   type,extends(type_internal_variable) :: type_bulk_variable
      ! Metadata
      real(rk)                                      :: vertical_movement         = 0.0_rk
      real(rk)                                      :: specific_light_extinction = 0.0_rk
      logical                                       :: no_precipitation_dilution = .false.
      logical                                       :: no_river_dilution         = .false.
      integer                                       :: time_treatment            = time_treatment_last
      type (type_bulk_standard_variable)            :: standard_variable

      ! Arrays with all associated data and index pointers.
      type (type_bulk_data_pointer_pointer),allocatable :: alldata(:) 
   end type type_bulk_variable

   type,extends(type_internal_variable) :: type_horizontal_variable
      ! Metadata
      integer                                  :: time_treatment = time_treatment_last
      integer                                  :: domain         = domain_bottom
      type (type_horizontal_standard_variable) :: standard_variable
      
      ! Arrays with all associated data and index pointers.
      type (type_horizontal_data_pointer_pointer),allocatable :: alldata(:) 
   end type type_horizontal_variable

   type,extends(type_internal_variable) :: type_scalar_variable
      ! Metadata
      integer                              :: time_treatment = time_treatment_last
      type (type_global_standard_variable) :: standard_variable

      ! Arrays with all associated data and index pointers.
      type (type_scalar_data_pointer_pointer),allocatable :: alldata(:) 
   end type type_scalar_variable

   type type_coupling
      character(len=attribute_length)   :: master   = ''
      type (type_model_pointer),pointer :: source   => null()
      logical                           :: required = .false.
   end type

   type type_link
      character(len=attribute_length)       :: name    = ''
      class (type_internal_object), pointer :: target  => null()
      type (type_coupling)                  :: coupling
      logical                               :: owner   = .true.
      type (type_link), pointer             :: next    => null()
   end type

   type type_model_pointer
      class (type_base_model),pointer :: p => null()
      type (type_state_variable_id),allocatable :: state(:)
   end type

   type type_model_pointer_pointer
      type (type_model_pointer),pointer :: p => null()
   end type

   type,extends(type_id) :: type_model_id
      type (type_model_pointer) :: model
   end type

   type,extends(type_internal_object) :: type_model_reference
      type (type_model_pointer_pointer),allocatable :: pointers(:)
   end type

   ! ====================================================================================================
   ! Variable identifiers used by host models.
   ! ====================================================================================================

   type type_bulk_variable_id
      class (type_bulk_variable),    pointer            :: variable => null()
      type (type_bulk_data_pointer),pointer             :: p        => null()
      type (type_bulk_data_pointer_pointer),allocatable :: alldata(:)
      integer                                           :: state_index = -1
      integer                                           :: write_index = -1
   end type

   type type_horizontal_variable_id
      class (type_horizontal_variable),    pointer            :: variable => null()
      type (type_horizontal_data_pointer),pointer             :: p        => null()
      type (type_horizontal_data_pointer_pointer),allocatable :: alldata(:)
      integer                                                 :: state_index = -1
      integer                                                 :: write_index = -1
   end type

   type type_scalar_variable_id
      class (type_scalar_variable),    pointer            :: variable => null()
      type (type_scalar_data_pointer),pointer             :: p        => null()
      type (type_scalar_data_pointer_pointer),allocatable :: alldata(:)
      integer                                             :: state_index = -1
      integer                                             :: write_index = -1
   end type

   ! ====================================================================================================
   ! Types to hold variable metadata, used by host models.
   ! ====================================================================================================

   type,abstract :: type_external_variable
      character(len=attribute_length) :: name          = ''
      character(len=attribute_length) :: long_name     = ''
      character(len=attribute_length) :: units         = ''
      real(rk)                        :: minimum       = -1.e20_rk
      real(rk)                        :: maximum       =  1.e20_rk
      real(rk)                        :: missing_value = -2.e20_rk
      type (type_property_dictionary) :: properties
   end type

!  Derived type describing a state variable
   type,extends(type_external_variable) :: type_state_variable_info
      type (type_bulk_standard_variable) :: standard_variable
      real(rk)                           :: initial_value             = 0.0_rk
      real(rk)                           :: vertical_movement         = 0.0_rk  ! Vertical movement (m/s, <0: sinking, >0: floating)
      real(rk)                           :: specific_light_extinction = 0.0_rk  ! Specific light extinction (/m/state variable unit)
      logical                            :: no_precipitation_dilution = .false.
      logical                            :: no_river_dilution         = .false.
      integer                            :: externalid                = 0       ! Identifier to be used freely by host
      type (type_bulk_variable_id)       :: globalid
   end type type_state_variable_info

   type,extends(type_external_variable) :: type_horizontal_state_variable_info
      type (type_horizontal_standard_variable) :: standard_variable
      real(rk)                                 :: initial_value = 0.0_rk
      integer                                  :: externalid    = 0             ! Identifier to be used freely by host
      type (type_horizontal_variable_id)       :: globalid
   end type type_horizontal_state_variable_info

!  Derived type describing a diagnostic variable
   type,extends(type_external_variable) :: type_diagnostic_variable_info
      type (type_bulk_standard_variable) :: standard_variable
      integer                            :: externalid     = 0                   ! Identifier to be used freely by host
      integer                            :: time_treatment = time_treatment_last ! See time_treatment_* parameters above
      type (type_bulk_variable_id)       :: globalid
   end type type_diagnostic_variable_info

   type,extends(type_external_variable) :: type_horizontal_diagnostic_variable_info
      type (type_horizontal_standard_variable) :: standard_variable
      integer                                  :: externalid     = 0                   ! Identifier to be used freely by host
      integer                                  :: time_treatment = time_treatment_last ! See time_treatment_* parameters above
      type (type_horizontal_variable_id)       :: globalid
   end type type_horizontal_diagnostic_variable_info

!  Derived type describing a conserved quantity
   type,extends(type_external_variable) :: type_conserved_quantity_info
      type (type_bulk_standard_variable)            :: standard_variable
      integer                                       :: externalid = 0       ! Identifier to be used freely by host
      type (type_weighted_sum),pointer              :: model => null()
      type (type_horizontal_weighted_sum),pointer   :: horizontal_model => null()
      integer,allocatable                           :: state_indices(:)
      real(rk),allocatable                          :: state_scale_factors(:)
      integer                                       :: rate_diag_index = -1
      integer                                       :: horizontal_rate_diag_index = -1
   end type type_conserved_quantity_info

   ! ====================================================================================================
   ! Types that describe custom expressions (arbitrary functions of one or more variables).
   ! ====================================================================================================

   type,abstract :: type_expression
      class (type_expression), pointer :: next        => null()
      character(len=attribute_length)  :: output_name = ''
   end type

   type,abstract,extends(type_expression) :: type_bulk_expression
      type (type_bulk_data_pointer),pointer :: out
   end type

   type,abstract,extends(type_expression) :: type_horizontal_expression
      type (type_horizontal_data_pointer),pointer :: out
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

      ! Arrays with variable metadata [used by hosts only]
      type (type_state_variable_info),                allocatable,dimension(:) :: state_variables                 
      type (type_horizontal_state_variable_info),     allocatable,dimension(:) :: surface_state_variables         
      type (type_horizontal_state_variable_info),     allocatable,dimension(:) :: bottom_state_variables          
      type (type_diagnostic_variable_info),           allocatable,dimension(:) :: diagnostic_variables            
      type (type_horizontal_diagnostic_variable_info),allocatable,dimension(:) :: horizontal_diagnostic_variables 
      type (type_conserved_quantity_info),            allocatable,dimension(:) :: conserved_quantities            

      ! Pointers for backward compatibility (pre 2013-06-15) [used by hosts only]
      type (type_horizontal_state_variable_info),     pointer,dimension(:) :: state_variables_ben     => null()
      type (type_horizontal_diagnostic_variable_info),pointer,dimension(:) :: diagnostic_variables_hz => null()

      ! Arrays with names of variables read by one or more biogeochemical models.
      ! These are not used within FABM, but may be accessed by the host to determine the names of
      ! potential forcing variables.
      character(len=attribute_length),allocatable,dimension(:) :: dependencies        
      character(len=attribute_length),allocatable,dimension(:) :: dependencies_hz     
      character(len=attribute_length),allocatable,dimension(:) :: dependencies_scalar 

      ! Pointers to linked models in the model tree.
      class (type_base_model),pointer :: parent => null()
      type (type_model_list)          :: children

      ! Model name and variable prefixes.
      character(len=64) :: name             = ''
      character(len=64) :: name_prefix      = ''
      character(len=64) :: long_name_prefix = ''

      ! Models constituents: links to variables, coupling requests, parameters, expressions
      type (type_link),              pointer :: first_link => null()
      type (type_aggregate_variable),pointer :: first_aggregate_variable => null()

      type (type_property_dictionary) :: parameters
      type (type_set)                 :: retrieved_parameters, missing_parameters

      class (type_expression), pointer :: first_expression => null()

      real(rk) :: dt = 1.0_rk

      logical :: check_conservation = .false.

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
      procedure :: register_bulk_dependency
      procedure :: register_bulk_dependency_sn
      procedure :: register_horizontal_dependency
      procedure :: register_horizontal_dependency_sn
      procedure :: register_global_dependency
      procedure :: register_global_dependency_sn
      procedure :: register_standard_conserved_quantity
      procedure :: register_custom_conserved_quantity
      procedure :: register_bulk_state_dependency_ex
      procedure :: register_bottom_state_dependency_ex
      procedure :: register_surface_state_dependency_ex
      procedure :: register_bulk_state_dependency_old
      procedure :: register_bottom_state_dependency_old
      procedure :: register_surface_state_dependency_old
      procedure :: register_model_dependency

      procedure :: add_object_copy
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
      generic :: register_dependency          => register_bulk_dependency, register_bulk_dependency_sn, &
                                                 register_horizontal_dependency, register_horizontal_dependency_sn, &
                                                 register_global_dependency, register_global_dependency_sn, &
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
      procedure :: do_bottom                => base_do_bottom
      procedure :: do_surface               => base_do_surface
      procedure :: do_ppdd                  => base_do_ppdd
      procedure :: do_bottom_ppdd           => base_do_bottom_ppdd

      ! Advanced functionality: variable vertical movement and light attenuation, feedbacks to drag and albedo.
      procedure :: get_vertical_movement    => base_get_vertical_movement
      procedure :: get_light_extinction     => base_get_light_extinction
      procedure :: get_drag                 => base_get_drag
      procedure :: get_albedo               => base_get_albedo

      ! Bookkeeping: calculate total of conserved quantities, check and repair model state.
      procedure :: get_conserved_quantities => base_get_conserved_quantities
      procedure :: state_to_conserved_quantities => base_state_to_conserved_quantities
      procedure :: get_horizontal_conserved_quantities => base_get_horizontal_conserved_quantities
      procedure :: check_state              => base_check_state
      procedure :: fatal_error              => base_fatal_error

      ! For backward compatibility only - do not use these in new models!
      procedure :: set_domain               => base_set_domain
      procedure :: do_benthos               => base_do_benthos           ! superceded by do_bottom
      procedure :: do_benthos_ppdd          => base_do_benthos_ppdd      ! superceded by do_bottom_ppdd
      procedure :: get_surface_exchange     => base_get_surface_exchange ! superceded by do_surface

   end type type_base_model

   ! ====================================================================================================
   ! Derived type for holding global data needed by biogeochemical model tree.
   ! ====================================================================================================

   type type_environment
      ! Declare the arrays for diagnostic variable values.
      real(rk),allocatable _DIMENSION_GLOBAL_PLUS_1_            :: diag    
      real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: diag_hz 

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

   ! ====================================================================================================
   ! Base type for a FABM host (provides routines for logging and error reporting)
   ! A host model that wants to process log message and fatal errors themselves (rather then the default
   ! behavior: log messages to stdout, fatal error to stdout followed by STOP) must create a derived
   ! type that extends type_base_driver. To use the custom type, allocate "host" in fabm_types with
   ! the custom type, e.g., with "allocate(type_custom_host::host)". This must be done before any FABM
   ! routine is called!
   ! ====================================================================================================

   type,public :: type_base_driver
   contains
      procedure :: fatal_error => base_driver_fatal_error
      procedure :: log_message => base_driver_log_message
   end type

   class (type_base_driver),pointer,save,public :: driver => null()

   type type_component
      character(len=attribute_length) :: name   = ''
      real(rk)                        :: weight = 1._rk
      type (type_dependency_id)       :: id
      type (type_component),pointer   :: next   => null()
   end type

   type type_horizontal_component
      character(len=attribute_length)          :: name   = ''
      real(rk)                                 :: weight = 1._rk
      type (type_horizontal_dependency_id)     :: id
      type (type_horizontal_component),pointer :: next   => null()
   end type

   type,extends(type_base_model) :: type_weighted_sum
      character(len=attribute_length) :: output_name      = ''
      character(len=attribute_length) :: output_long_name = ''
      character(len=attribute_length) :: output_units     = ''
      type (type_diagnostic_variable_id) :: id_output
      type (type_component),pointer   :: first => null()
   contains
      procedure :: initialize    => weighted_sum_initialize
      procedure :: add_component => weighted_sum_add_component
      procedure :: evaluate      => weighted_sum_evaluate
      procedure :: do            => weighted_sum_do
   end type

   type,extends(type_base_model) :: type_horizontal_weighted_sum
      character(len=attribute_length) :: output_name      = ''
      character(len=attribute_length) :: output_long_name = ''
      character(len=attribute_length) :: output_units     = ''
      type (type_horizontal_diagnostic_variable_id) :: id_output
      type (type_horizontal_component),pointer   :: first => null()
   contains
      procedure :: add_component       => horizontal_weighted_sum_add_component
      procedure :: initialize          => horizontal_weighted_sum_initialize
      procedure :: evaluate_horizontal => horizontal_weighted_sum_evaluate_horizontal
      procedure :: do_bottom           => horizontal_weighted_sum_do_bottom
   end type

   ! ====================================================================================================
   ! Interfaces
   ! ====================================================================================================

   interface create_external_variable_id
      module procedure create_external_bulk_id
      module procedure create_external_horizontal_id
      module procedure create_external_scalar_id
      module procedure create_external_bulk_id_for_standard_name
      module procedure create_external_horizontal_id_for_standard_name
      module procedure create_external_scalar_id_for_standard_name
   end interface

   interface register_state_variable
      module procedure register_bulk_state_variable
      module procedure register_bottom_state_variable
      module procedure register_surface_state_variable
   end interface

   interface register_state_dependency
      module procedure register_bulk_state_dependency_ex
      module procedure register_bottom_state_dependency_ex
      module procedure register_surface_state_dependency_ex
      module procedure register_bulk_state_dependency_old
      module procedure register_bottom_state_dependency_old
      module procedure register_surface_state_dependency_old
   end interface
   
   interface register_bulk_state_dependency
      module procedure register_bulk_state_dependency_ex
      module procedure register_bulk_state_dependency_old
   end interface

   interface register_bottom_state_dependency
      module procedure register_bottom_state_dependency_ex
      module procedure register_bottom_state_dependency_old
   end interface

   interface register_surface_state_dependency
      module procedure register_surface_state_dependency_ex
      module procedure register_surface_state_dependency_old
   end interface

   interface register_dependency
      module procedure register_bulk_dependency
      module procedure register_bulk_dependency_sn
      module procedure register_horizontal_dependency
      module procedure register_horizontal_dependency_sn
      module procedure register_global_dependency
      module procedure register_global_dependency_sn
   end interface

   interface register_diagnostic_variable
      module procedure register_bulk_diagnostic_variable
      module procedure register_horizontal_diagnostic_variable
   end interface

   interface register_conserved_quantity
      module procedure register_standard_conserved_quantity
      module procedure register_custom_conserved_quantity
   end interface

   interface append_data_pointer
      module procedure append_bulk_data_pointer
      module procedure append_horizontal_data_pointer
      module procedure append_scalar_data_pointer
   end interface

   interface compare_standard_variables
      module procedure compare_bulk_standard_variables
      module procedure compare_horizontal_standard_variables
      module procedure compare_global_standard_variables
   end interface
   
   interface is_null_standard_variable
      module procedure is_null_bulk_standard_variable
      module procedure is_null_horizontal_standard_variable
      module procedure is_null_global_standard_variable
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
      class (type_base_model),intent(in) ::  self
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

      integer  :: i
      real(rk) :: val,curext

      ! Use variable-specific light extinction coefficients.
      do i=1,size(self%state_variables)
         curext = self%state_variables(i)%specific_light_extinction
         if (curext/=0.0_rk) then
            _LOOP_BEGIN_
               _GET_STATE_EX_(environment,self%state_variables(i)%globalid%p,val)
               _SET_EXTINCTION_(val*curext)
            _LOOP_END_
         end if
      end do
   end subroutine

   subroutine base_get_drag(self,_ARGUMENTS_GET_DRAG_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_DRAG_
   end subroutine

   subroutine base_get_albedo(self,_ARGUMENTS_GET_ALBEDO_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_ALBEDO_
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
   end subroutine base_check_state

   subroutine base_fatal_error(self,location,message)
      class (type_base_model), intent(in) :: self
      character(len=*),        intent(in) :: location,message
      if (self%name/='') then
         call driver%fatal_error('model "'//trim(self%name)//'", '//trim(location),message)
      else
         call driver%fatal_error(location,message)
      end if
   end subroutine

   ! For backward compatibility only
   subroutine base_set_domain(self _ARG_LOCATION_)
      class (type_base_model),intent(inout) :: self
      _DECLARE_LOCATION_ARG_
   end subroutine base_set_domain

   subroutine base_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
   end subroutine base_do_benthos

   subroutine base_do_benthos_ppdd(self,_ARGUMENTS_DO_BOTTOM_PPDD_)
      class (type_base_model),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_
   end subroutine base_do_benthos_ppdd

   subroutine base_get_surface_exchange(self,_ARGUMENTS_DO_SURFACE_)
      class (type_base_model), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_
   end subroutine base_get_surface_exchange

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes model information.
!
! !INTERFACE:
   subroutine add_child(self,model,name,long_name,configunit)
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
      model%name_prefix = trim(name)//'_'
      if (present(long_name)) then
         model%long_name_prefix = trim(long_name)//' '
      else
         model%long_name_prefix = trim(name)//' '
      end if
      model%parent => self
      call self%children%append(model)
      call model%initialize(configunit)
   end subroutine add_child
!EOC

   subroutine set_variable_property_real(self,variable,name,value)
      class (type_base_model),intent(inout) :: self
      class (type_id),        intent(inout) :: variable
      character(len=*),       intent(in)    :: name
      real(rk),               intent(in)    :: value
      if (.not.associated(variable%link)) call self%fatal_error('set_variable_property_real','variable has not been registered')
      call variable%link%target%properties%set_real(name,value)
   end subroutine

   subroutine set_variable_property_integer(self,variable,name,value)
      class (type_base_model),intent(inout) :: self
      class (type_id),        intent(inout) :: variable
      character(len=*),       intent(in)    :: name
      integer,                intent(in)    :: value
      if (.not.associated(variable%link)) call self%fatal_error('set_variable_property_integer','variable has not been registered')
      call variable%link%target%properties%set_integer(name,value)
   end subroutine

   subroutine set_variable_property_logical(self,variable,name,value)
      class (type_base_model),intent(inout) :: self
      class (type_id),        intent(inout) :: variable
      character(len=*),       intent(in)    :: name
      logical,                intent(in)    :: value
      if (.not.associated(variable%link)) call self%fatal_error('set_variable_property_logical','variable has not been registered')
      call variable%link%target%properties%set_logical(name,value)
   end subroutine

   subroutine add_to_aggregate_variable(self,target,variable_id,scale_factor)
      class (type_base_model),           intent(inout) :: self
      type (type_bulk_standard_variable),intent(in)    :: target
      class (type_id),                   intent(inout) :: variable_id
      real(rk),optional,                 intent(in)    :: scale_factor

      if (.not.associated(variable_id%link)) &
         call self%fatal_error('add_to_aggregate_variable','variable has not been registered')

      select type (variable=>variable_id%link%target)
         class is (type_internal_variable)
            call variable%contributions%add(target,scale_factor)
         class default
            call self%fatal_error('add_to_aggregate_variable','Only variables can contribute to conserved quantities at present.')
      end select

   end subroutine

   subroutine contribution_list_add(self,target,scale_factor)
      class (type_contribution_list),    intent(inout) :: self
      type (type_bulk_standard_variable),intent(in)    :: target
      real(rk),optional,                 intent(in)    :: scale_factor

      class (type_contribution),pointer :: contribution

      if (.not.associated(self%first)) then
         allocate(self%first)
         contribution => self%first
      else
         contribution => self%first
         do while (associated(contribution))
            if (compare_standard_variables(target,contribution%target)) exit
            contribution => contribution%next
         end do
         if (.not.associated(contribution)) then
            contribution => self%first
            do while (associated(contribution%next))
               contribution => contribution%next
            end do
            allocate(contribution%next)
            contribution => contribution%next
         end if
      end if

      contribution%target = target
      if (present(scale_factor)) contribution%scale_factor = scale_factor
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

subroutine create_link(self,target,name,merge,owner)
   class (type_base_model),            intent(inout) :: self
   class (type_internal_object),pointer              :: target
   character(len=*),                   intent(in)    :: name
   logical,optional,                   intent(in)    :: merge,owner

   logical                  :: merge_eff
   type (type_link),pointer :: link
   class (type_internal_object),pointer :: existing_target

   merge_eff = .false.
   if (present(merge)) merge_eff = merge

   ! First check if a link with this name exists. If so, merge new target with old target.
   existing_target => self%find_object(name)
   if (associated(existing_target)) then
      if (.not.merge_eff) call self%fatal_error('create_link_for_object','Link with name "'//trim(name)//'" already exists.')
      call merge_variables(existing_target,target)
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
end subroutine create_link

recursive subroutine add_alias(self,target,name)
   class (type_base_model),            intent(inout) :: self
   class (type_id),                    intent(in)    :: target
   character(len=*),                   intent(in)    :: name

   call create_link(self,target%link%target,name,owner=.false.)
   if (associated(self%parent)) call self%parent%add_alias(target,trim(self%name_prefix)//trim(name))
end subroutine

subroutine request_coupling_for_link(self,link,master,required,source)
   class (type_base_model),intent(inout)              :: self
   type (type_link),       intent(inout)              :: link
   character(len=*),       intent(in)                 :: master
   logical,optional,       intent(in)                 :: required
   type (type_model_id),   intent(in),target,optional :: source

   link%coupling%master = master
   link%coupling%required = .true.
   nullify(link%coupling%source)
   if (present(required)) link%coupling%required = required
   if (present(source)) link%coupling%source => source%model
end subroutine request_coupling_for_link

subroutine request_coupling_for_name(self,slave,master,required,source)
   class (type_base_model),intent(inout)              :: self
   character(len=*),       intent(in)                 :: slave,master
   logical,optional,       intent(in)                 :: required
   type (type_model_id),   intent(in),target,optional :: source

   type (type_link),pointer :: link

   link => self%find_link(slave)
   if (.not.associated(link)) call self%fatal_error('request_coupling_for_name','Variable "'//trim(slave)//'" was not found.')
   call self%request_coupling(link,master,required,source)
end subroutine request_coupling_for_name

subroutine request_coupling_for_id(self,id,master,required,source)
   class (type_base_model),intent(inout)              :: self
   class (type_id),        intent(inout)              :: id
   character(len=*),       intent(in)                 :: master
   logical,optional,       intent(in)                 :: required
   type (type_model_id),   intent(in),target,optional :: source

   call self%request_coupling(id%link,master,required,source)
end subroutine request_coupling_for_id

subroutine append_index(array,index)
   type (type_integer_pointer),allocatable :: array(:)
   integer,target :: index
   type (type_integer_pointer),allocatable :: oldarray(:)

   ! Create a new list of integer pointers, or extend it if already allocated.
   if (.not.allocated(array)) then
      allocate(array(1))
   else
      allocate(oldarray(size(array)))
      oldarray = array
      deallocate(array)
      allocate(array(size(oldarray)+1))
      array(1:size(oldarray)) = oldarray
   end if

   ! Add pointer to provided integer to the list.
   array(size(array))%p => index
end subroutine append_index

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
end subroutine append_scalar_data_pointer

subroutine append_model_pointer(array,data)
   type (type_model_pointer_pointer),allocatable :: array(:)
   type (type_model_pointer),target :: data
   type (type_model_pointer_pointer),allocatable :: oldarray(:)
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
   array(size(array))%p%p => array(1)%p%p
end subroutine append_model_pointer

subroutine append_string(array,string,exists)
   character(len=attribute_length),dimension(:),allocatable :: array
   character(len=*),intent(in) :: string
   logical,intent(out),optional :: exists
   integer :: i
   character(len=attribute_length),allocatable :: oldarray(:)

   if (.not.allocated(array)) then
      allocate(array(1))
   else
      do i=1,size(array)
         if (array(i)==string) then
            if (present(exists)) exists = .true.
            return
         end if
      end do

      allocate(oldarray(size(array)))
      oldarray = array
      deallocate(array)
      allocate(array(size(oldarray)+1))
      array(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   array(size(array)) = string
   if (present(exists)) exists = .false.
end subroutine append_string

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
      class (type_base_model),intent(inout) :: self
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
   class (type_link),pointer :: link
   class (type_aggregate_variable),pointer :: aggregate_variable
!-----------------------------------------------------------------------
!BOC
      if (associated(self%parent)) call self%fatal_error('freeze_model_info', &
         'BUG: freeze_model_info can only operate on the root model.')

      ! First resolve model references, needed during variable coupling.
      ! Stage 1: resolve indirect references (references to other references)
      ! Stage 2: resolve direct references to models
      call process_coupling_tasks(self,process_model_references=.true.)
      call resolve_model_references(self)

      ! Determine whether all model references have been resolved.
      ! (resolve_model_references only processes references with explicit coupling requests,
      ! potentially leaving those without coupling request anywhere in the model tree uncoupled.)
      link => self%first_link
      do while (associated(link))
         if (link%owner) then
            select type (slave=>link%target)
               class is (type_model_reference)
                  if (.not.associated(slave%pointers(1)%p%p)) &
                     call self%fatal_error('freeze_model_info','Model reference "'//trim(link%name)//'" was not coupled to an actual model.')
            end select
         end if
         link => link%next
      end do
      
      ! Now couple model variables
      ! Stage 1: implicit - couple variables based on overlapping standard identities.
      ! Stage 2: explicit - resolve user- or model-specified links between variables.
      call couple_standard_variables(self)
      call process_coupling_tasks(self,process_model_references=.false.)

      ! Create models for aggregate variables at root level, to be used to compute conserved quantities
      call build_aggregate_variables(self)
      aggregate_variable => self%first_aggregate_variable
      do while (associated(aggregate_variable))
         call create_aggregate_model(self,aggregate_variable)
         aggregate_variable => aggregate_variable%next
      end do
      call process_coupling_tasks(self,process_model_references=.false.)

      ! Add arrays with state variable identifiers to model references (this requires variable coupling to be complete!)
      call complete_model_references(self)

      ! Build final authorative arrays with variable metadata .
      call classify_variables(self)

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
                                           standard_variable,presence)
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
      real(rk),                           intent(in),optional :: minimum, maximum,missing_value
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
      type (type_bulk_variable) :: variable
      character(len=256)        :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_bulk_state_variable', &
         'State variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      variable%name      = name
      variable%units     = units
      variable%long_name = long_name
      if (present(initial_value))             variable%initial_value             = initial_value
      if (present(minimum))                   variable%minimum                   = minimum
      if (present(maximum))                   variable%maximum                   = maximum
      if (present(missing_value))             variable%missing_value             = missing_value
      if (present(vertical_movement))         variable%vertical_movement         = vertical_movement
      if (present(specific_light_extinction)) variable%specific_light_extinction = specific_light_extinction
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

      id%link => self%add_object_copy(variable)

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
                                             standard_variable,presence)
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
      real(rk),                                 intent(in),optional :: minimum, maximum,missing_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
      integer,                                  intent(in),optional :: presence
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable) :: variable
      character(len=256)              :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_bottom_state_variable', &
         'State variables may only be registered during initialization.')

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

      id%link => self%add_object_copy(variable)

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
                                              standard_variable,presence)
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
      real(rk),                                 intent(in),optional :: minimum, maximum,missing_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
      integer,                                  intent(in),optional :: presence
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable) :: variable
      character(len=256)              :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_surface_state_variable', &
         'State variables may only be registered during initialization.')

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

      id%link => self%add_object_copy(variable)

   end subroutine register_surface_state_variable
!EOC

   function add_object_copy(self,object) result(link)
      ! This subroutine allocates and initializes a copy of the provided object,
      ! (to ensure it is persistent in memory), then creates a link to it in the list
      ! of model entities.
      class (type_base_model),target,intent(inout) :: self
      class (type_internal_object),  intent(in)    :: object
      type (type_link),pointer                     :: link

      class (type_internal_object),pointer :: copy
      character(len=attribute_length)      :: local_name

      ! Store the local object name (overwritten by self%add_object)
      local_name = object%name

      ! Create a persistent copy of the object.
      allocate(copy,source=object)

      ! Add the persistent copy of the object to the model.
      call self%add_object(copy)

      link => self%find_link(local_name)
   end function

   recursive subroutine add_object(self,object)
      ! This subroutine creates a link to the supplied object, then allows
      ! parent models to do the same.
      ! NB this subroutine MUST be recursive, to allow parent models to override
      ! the properties of objects added by their child models.
      class (type_base_model),target,      intent(inout) :: self
      class (type_internal_object),pointer               :: object

      if (object%long_name=='') object%long_name = object%name

      call create_link(self,object,object%name,merge=.true.)
      if (.not.associated(object)) return ! if create_link has merged the new object into an existing object

      ! Forward to parent
      if (associated(self%parent)) then
         object%name = trim(self%name_prefix)//trim(object%name)
         object%long_name = trim(self%long_name_prefix)//' '//trim(object%long_name)
         call self%parent%add_object(object)
      end if
   end subroutine
   
   subroutine connect_bulk_state_variable_id(self,variable,id)
      class (type_base_model),       intent(in)           :: self
      type (type_bulk_variable),     intent(inout),target :: variable
      type (type_state_variable_id), intent(inout),target :: id

      if (associated(id%link)) call self%fatal_error('connect_bulk_state_variable_id', &
         'Identifier supplied for '//trim(variable%name)//' is already associated with '//trim(id%link%name)//'.')

      call append_index(variable%state_indices,id%state_index)
      call append_data_pointer(variable%alldata,id%data)
   end subroutine

   subroutine connect_bottom_state_variable_id(self,variable,id)
      class (type_base_model),       intent(in)           :: self
      type (type_horizontal_variable),      intent(inout),target :: variable
      type (type_bottom_state_variable_id), intent(inout),target :: id

      if (associated(id%link)) call self%fatal_error('connect_bottom_state_variable_id', &
         'Identifier supplied for '//trim(variable%name)//' is already associated with '//trim(id%link%name)//'.')

      call append_index(variable%state_indices,id%bottom_state_index)
      call append_data_pointer(variable%alldata,id%horizontal_data)
   end subroutine

   subroutine connect_surface_state_variable_id(self,variable,id)
      class (type_base_model),       intent(in)           :: self
      type (type_horizontal_variable),      intent(inout),target :: variable
      type (type_surface_state_variable_id),intent(inout),target :: id

      if (associated(id%link)) call self%fatal_error('connect_surface_state_variable_id', &
         'Identifier supplied for '//trim(variable%name)//' is already associated with '//trim(id%link%name)//'.')

      call append_index(variable%state_indices,id%surface_state_index)
      call append_data_pointer(variable%alldata,id%horizontal_data)
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   subroutine register_bulk_diagnostic_variable(self, id, name, units, long_name, &
                                                time_treatment, missing_value, standard_variable)
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
      integer,                            intent(in),optional :: time_treatment
      real(rk),                           intent(in),optional :: missing_value
      type (type_bulk_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable) :: variable
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_bulk_diagnostic_variable',&
                                             'Diagnostic variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      variable%name      = name
      variable%units     = units
      variable%long_name = long_name
      variable%source_model => self
      if (present(time_treatment))    variable%time_treatment    = time_treatment
      if (present(missing_value))     variable%missing_value     = missing_value
      if (present(standard_variable)) variable%standard_variable = standard_variable
      call append_index(variable%write_indices,id%diag_index)
      if (associated(id%link)) call self%fatal_error('register_bulk_diagnostic_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      id%link => self%add_object_copy(variable)

   end subroutine register_bulk_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   subroutine register_horizontal_diagnostic_variable(self, id, name, units, long_name, &
                                                      time_treatment, missing_value, standard_variable)
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
      integer,                                  intent(in),optional :: time_treatment
      real(rk),                                 intent(in),optional :: missing_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable) :: variable
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_horizontal_diagnostic_variable',&
                                             'Diagnostic variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      variable%name      = name
      variable%units     = units
      variable%long_name = long_name
      variable%source_model => self
      if (present(time_treatment))    variable%time_treatment    = time_treatment
      if (present(missing_value))     variable%missing_value     = missing_value
      if (present(standard_variable)) variable%standard_variable = standard_variable

      if (associated(id%link)) call self%fatal_error('register_horizontal_diagnostic_variable', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')
      call append_index(variable%write_indices,id%horizontal_diag_index)

      id%link => self%add_object_copy(variable)

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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable)       :: variable
      character(len=attribute_length) :: name_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_standard_conserved_quantity',&
                                             'Conserved quantities may only be registered during initialization.')

      if (present(name)) then
         name_eff = name
      else
         name_eff = standard_variable%name
      end if

      ! Either use the provided variable object, or create a new one.
      variable%standard_variable = standard_variable
      variable%name      = standard_variable%name
      variable%units     = standard_variable%units
      variable%long_name = standard_variable%name

      if (associated(id%link)) call self%fatal_error('register_conserved_quantity', &
         'Identifier supplied for '//trim(name_eff)//' is already associated with '//trim(id%link%name)//'.')
      call append_index(variable%cons_indices,id%cons_index)

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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable) :: variable
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_custom_conserved_quantity',&
         'Conserved quantities may only be registered during initialization.')

      variable%name      = name
      variable%units     = units
      variable%long_name = long_name

      if (associated(id%link)) call self%fatal_error('register_conserved_quantity', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')
      call append_index(variable%cons_indices,id%cons_index)

      id%link => self%add_object_copy(variable)

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

   subroutine register_bulk_state_dependency_old(model,id,name)
      class (type_base_model),      intent(inout)        :: model
      type (type_state_variable_id),intent(inout),target :: id
      character(len=*),             intent(in)           :: name

      call register_bulk_state_dependency(model, id, name, '', name)
      call model%request_coupling(id,name)
   end subroutine

   subroutine register_bottom_state_dependency_old(model,id,name)
      class (type_base_model),             intent(inout)        :: model
      type (type_bottom_state_variable_id),intent(inout),target :: id
      character(len=*),                    intent(in)           :: name

      call register_bottom_state_dependency(model, id, name, '', name)
      call model%request_coupling(id,name)
   end subroutine

   subroutine register_surface_state_dependency_old(model,id,name)
      class (type_base_model),              intent(inout)        :: model
      type (type_surface_state_variable_id),intent(inout),target :: id
      character(len=*),                     intent(in)           :: name

      call register_surface_state_dependency(model, id, name, '', name)
      call model%request_coupling(id,name)
   end subroutine
   
   subroutine register_bulk_dependency_sn(model,id,standard_variable,required)
      class (type_base_model),           intent(inout)        :: model
      type (type_dependency_id),         intent(inout),target :: id
      type (type_bulk_standard_variable),intent(in)           :: standard_variable
      logical,                           intent(in),optional  :: required

      call register_bulk_dependency(model,id,standard_variable%name,standard_variable%units,standard_variable=standard_variable)
   end subroutine register_bulk_dependency_sn

   subroutine register_horizontal_dependency_sn(model,id,standard_variable,required)
      class (type_base_model),                 intent(inout)        :: model
      type (type_horizontal_dependency_id),    intent(inout),target :: id
      type (type_horizontal_standard_variable),intent(in)           :: standard_variable
      logical,                                 intent(in),optional :: required

      call register_horizontal_dependency(model,id,standard_variable%name,standard_variable%units,standard_variable=standard_variable)
   end subroutine register_horizontal_dependency_sn

   subroutine register_global_dependency_sn(model,id,standard_variable,required)
      class (type_base_model),             intent(inout)        :: model
      type (type_global_dependency_id),    intent(inout),target :: id
      type (type_global_standard_variable),intent(in)           :: standard_variable
      logical,                             intent(in),optional :: required

      call register_global_dependency(model,id,standard_variable%name,standard_variable%units,standard_variable=standard_variable)
   end subroutine register_global_dependency_sn

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a variable defined on
! the full model domain.
!
! !INTERFACE:
   subroutine register_bulk_dependency(self,id,name,units,long_name,standard_variable,required)
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
      type (type_bulk_variable) :: variable
      logical                   :: required_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_bulk_dependency',&
         'Dependencies may only be registered during initialization.')

      required_eff = .true.
      if (present(required)) required_eff = required

      variable%name = name
      if (present(units))     variable%units     = units
      if (present(long_name)) variable%long_name = long_name
      if (present(standard_variable)) variable%standard_variable = standard_variable
      variable%presence = presence_external_required
      if (.not.required_eff) variable%presence = presence_external_optional

      call append_data_pointer(variable%alldata,id%data)
      if (associated(id%link)) call self%fatal_error('register_bulk_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      id%link => self%add_object_copy(variable)

      if (associated(self%parent).and..not.present(standard_variable)) call self%request_coupling(id,name,required=.false.)

   end subroutine register_bulk_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a variable defined on
! a horizontal slice of the model domain.
!
! !INTERFACE:
   subroutine register_horizontal_dependency(self,id,name,units,long_name,standard_variable,required)
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
      type (type_horizontal_variable) :: variable
      logical                         :: required_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_horizontal_dependency',&
         'Dependencies may only be registered during initialization.')

      required_eff = .true.
      if (present(required)) required_eff = required

      variable%name = name
      if (present(units))     variable%units     = units
      if (present(long_name)) variable%long_name = long_name
      if (present(standard_variable)) variable%standard_variable = standard_variable
      variable%presence = presence_external_required
      if (.not.required_eff) variable%presence = presence_external_optional

      call append_data_pointer(variable%alldata,id%horizontal_data)
      if (associated(id%link)) call self%fatal_error(':register_horizontal_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      id%link => self%add_object_copy(variable)

      if (associated(self%parent).and..not.present(standard_variable)) call self%request_coupling(id,name,required=.false.)

   end subroutine register_horizontal_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a global (space-
! independent) variable.
!
! !INTERFACE:
   subroutine register_global_dependency(self,id,name,units,long_name,standard_variable,required)
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
      type (type_scalar_variable) :: variable
      logical                     :: required_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error(':register_global_dependency',&
         'Dependencies may only be registered during initialization.')

      required_eff = .true.
      if (present(required)) required_eff = required

      variable%name = name
      if (present(units))     variable%units     = units
      if (present(long_name)) variable%long_name = long_name
      if (present(standard_variable)) variable%standard_variable = standard_variable
      variable%presence = presence_external_required
      if (.not.required_eff) variable%presence = presence_external_optional

      call append_data_pointer(variable%alldata,id%global_data)
      if (associated(id%link)) call self%fatal_error(':register_global_dependency', &
         'Identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')

      id%link => self%add_object_copy(variable)

      if (associated(self%parent).and..not.present(standard_variable)) call self%request_coupling(id,name,required=.false.)

   end subroutine register_global_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on another model.
!
! !INTERFACE:
   subroutine register_model_dependency(self,id,name)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_base_model), intent(inout)        :: self
      type (type_model_id),    intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),        intent(in)           :: name
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_model_reference) :: reference
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error(':register_model_dependency',&
         'Dependencies may only be registered during initialization.')

      reference%name = name

      if (associated(id%link)) call self%fatal_error(':register_model_dependency', &
         'Model identifier supplied for '//trim(name)//' is already associated with '//trim(id%link%name)//'.')
      call append_model_pointer(reference%pointers,id%model)

      id%link => self%add_object_copy(reference)

   end subroutine register_model_dependency
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

recursive subroutine get_real_parameter(self,value,name,units,long_name,scale_factor,default,found)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target  :: self
   real(rk),                intent(inout)          :: value
   character(len=*),        intent(in)             :: name
   character(len=*),        intent(in),   optional :: units,long_name
   real(rk),                intent(in),   optional :: scale_factor,default
   logical,                 intent(inout),optional :: found
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),    pointer :: property
   type (type_real_property),pointer :: current_parameter
   real(rk)                          :: default_eff
   logical                           :: found_eff,success
!
!-----------------------------------------------------------------------
!BOC
   if (self%retrieved_parameters%contains(name)) call self%fatal_error('get_real_parameter', &
      'Value for parameter "'//trim(name)//'" of model "'//trim(self%name)//'" has already been retrieved once.')
   call self%retrieved_parameters%add(name)

   if (present(default)) then
      default_eff = default
   else
      default_eff = value
   end if
   if (present(found)) then
      found_eff = found
   else
      found_eff = .false.
   end if

   ! Try to find the parameter in the model's own parameter list.
   property => self%parameters%get_property(name)
   if (associated(property)) then
      found_eff = .true.
      default_eff = property%to_real(success=success)
      if (.not.success) call self%fatal_error('get_real_parameter', &
         'Value "'//trim(property%to_string())//'" for parameter "'//trim(name)//'" of model "'//trim(self%name)//'" is not a real number.')
   end if

   if (associated(self%parent)) then
      call self%parent%get_parameter(value,trim(self%name)//'/'//name,units,long_name,1.0_rk,default_eff,found_eff)
   else
      value = default_eff
   end if

   if (.not.(present(found).or.found_eff)) call self%missing_parameters%add(name)

   ! Store parameter settings
   allocate(current_parameter)
   current_parameter%value = value

   call add_parameter(self,current_parameter,name,units,long_name)

   if (present(scale_factor)) value = value*scale_factor

end subroutine get_real_parameter
!EOC

recursive subroutine get_integer_parameter(self,value,name,units,long_name,default,found)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target :: self
   integer,                 intent(inout)         :: value
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
   integer,                 intent(in),optional   :: default
   logical,                 intent(inout),optional :: found
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),       pointer :: property
   type (type_integer_property),pointer :: current_parameter
   integer                              :: default_eff
   logical                              :: found_eff,success
!
!-----------------------------------------------------------------------
!BOC
   if (self%retrieved_parameters%contains(name)) call self%fatal_error('get_integer_parameter', &
      'Value for parameter "'//trim(name)//'" of model "'//trim(self%name)//'" has already been retrieved once.')
   call self%retrieved_parameters%add(name)

   if (present(default)) then
      default_eff = default
   else
      default_eff = value
   end if
   if (present(found)) then
      found_eff = found
   else
      found_eff = .false.
   end if

   ! Try to find the parameter in the model's own parameter list.
   property => self%parameters%get_property(name)
   if (associated(property)) then
      found_eff = .true.
      default_eff = property%to_integer(success=success)
      if (.not.success) call self%fatal_error('get_integer_parameter', &
         'Value "'//trim(property%to_string())//'" for parameter "'//trim(name)//'" of model "'//trim(self%name)//'" is not an integer.')
   end if

   if (associated(self%parent)) then
      call self%parent%get_parameter(value,trim(self%name)//'/'//name,units,long_name,default_eff,found_eff)
   else
      value = default_eff
   end if

   if (.not.(present(found).or.found_eff)) call self%missing_parameters%add(name)

   ! Store parameter settings
   allocate(current_parameter)
   current_parameter%value = value

   call add_parameter(self,current_parameter,name,units,long_name)

end subroutine get_integer_parameter
!EOC

recursive subroutine get_logical_parameter(self,value,name,units,long_name,default,found)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target :: self
   logical,                 intent(inout)         :: value
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
   logical,                 intent(in),optional   :: default
   logical,                 intent(inout),optional :: found
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),       pointer :: property
   type (type_logical_property),pointer :: current_parameter
   logical                              :: default_eff
   logical                              :: found_eff,success
!
!-----------------------------------------------------------------------
!BOC
   if (self%retrieved_parameters%contains(name)) call self%fatal_error('get_logical_parameter', &
      'Value for parameter "'//trim(name)//'" of model "'//trim(self%name)//'" has already been retrieved once.')
   call self%retrieved_parameters%add(name)

   if (present(default)) then
      default_eff = default
   else
      default_eff = value
   end if
   if (present(found)) then
      found_eff = found
   else
      found_eff = .false.
   end if

   ! Try to find the parameter in the model's own parameter list.
   property => self%parameters%get_property(name)
   if (associated(property)) then
      found_eff = .true.
      default_eff = property%to_logical(success=success)
      if (.not.success) call self%fatal_error('get_logical_parameter', &
         'Value "'//trim(property%to_string())//'" for parameter "'//trim(name)//'" of model "'//trim(self%name)//'" is not an Boolean.')
   end if
   
   if (associated(self%parent)) then
      call self%parent%get_parameter(value,trim(self%name)//'/'//name,units,long_name,default_eff,found_eff)
   else
      value = default_eff
   end if

   if (.not.(present(found).or.found_eff)) call self%missing_parameters%add(name)

   ! Store parameter settings
   allocate(current_parameter)
   current_parameter%value = value

   call add_parameter(self,current_parameter,name,units,long_name)

end subroutine get_logical_parameter
!EOC

recursive subroutine get_string_parameter(self,value,name,units,long_name,default,found)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target :: self
   character(len=*),        intent(inout)         :: value
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
   character(len=*),        intent(in),optional   :: default
   logical,                 intent(inout),optional :: found
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),       pointer :: property
   type (type_string_property),pointer :: current_parameter
   character(len=1024)                 :: default_eff
   logical                              :: found_eff
!
!-----------------------------------------------------------------------
!BOC
   if (self%retrieved_parameters%contains(name)) call self%fatal_error('get_string_parameter', &
      'Value for parameter "'//trim(name)//'" of model "'//trim(self%name)//'" has already been retrieved once.')
   call self%retrieved_parameters%add(name)

   if (present(default)) then
      default_eff = default
   else
      default_eff = value
   end if
   if (present(found)) then
      found_eff = found
   else
      found_eff = .false.
   end if

   ! Try to find the parameter in the model's own parameter list.
   property => self%parameters%get_property(name)
   if (associated(property)) then
      found_eff = .true.
      default_eff = property%to_string()
   end if

   if (associated(self%parent)) then
      call self%parent%get_parameter(value,trim(self%name)//'/'//name,units,long_name,default_eff,found_eff)
   else
      value = default_eff
   end if

   if (.not.(present(found).or.found_eff)) call self%missing_parameters%add(name)

   ! Store parameter settings
   allocate(current_parameter)
   current_parameter%value = value

   call add_parameter(self,current_parameter,name,units,long_name)

end subroutine get_string_parameter
!EOC

subroutine add_parameter(model,parameter,name,units,long_name)
! !INPUT PARAMETERS:
   class (type_base_model), intent(inout), target :: model
   class (type_property),   intent(inout), target :: parameter
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store metadata
   parameter%name = name
   if (present(units))     parameter%units     = units
   if (present(long_name)) parameter%long_name = long_name
   call model%parameters%set_property(parameter)

end subroutine add_parameter
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
      type (type_link),pointer                      :: link,link2
      character(len=attribute_length),allocatable :: processed_bulk(:),processed_horizontal(:),processed_scalar(:)
      logical                                       :: exists
!
!-----------------------------------------------------------------------
!BOC
      link => model%first_link
      do while (associated(link))
         select type (object=>link%target)
            class is (type_bulk_variable)
               if (.not.is_null_standard_variable(object%standard_variable)) then
                  call append_string(processed_bulk,object%standard_variable%name,exists=exists)
                  if (.not.exists) then
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_bulk_variable)
                              if (compare_standard_variables(object%standard_variable,object2%standard_variable) .and. .not. &
                                  is_null_standard_variable(object2%standard_variable)) &
                                 call couple_variables(model,link%target,link2%target)
                        end select
                        link2 => link2%next
                     end do
                  end if
               end if
            class is (type_horizontal_variable)
               if (.not.is_null_standard_variable(object%standard_variable)) then
                  call append_string(processed_horizontal,object%standard_variable%name,exists=exists)
                  if (.not.exists) then
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_horizontal_variable)
                              if (compare_standard_variables(object%standard_variable,object2%standard_variable) .and. .not. &
                                  is_null_standard_variable(object2%standard_variable)) &
                               call couple_variables(model,link%target,link2%target)
                        end select
                        link2 => link2%next
                     end do
                  end if
               end if
            class is (type_scalar_variable)
               if (.not.is_null_standard_variable(object%standard_variable)) then
                  call append_string(processed_scalar,object%standard_variable%name,exists=exists)
                  if (.not.exists) then
                     link2 => link%next
                     do while (associated(link2))
                        select type (object2=>link2%target)
                           class is (type_scalar_variable)
                              if (compare_standard_variables(object%standard_variable,object2%standard_variable) .and. .not. &
                                  is_null_standard_variable(object2%standard_variable)) &
                                 call couple_variables(model,link%target,link2%target)
                        end select
                        link2 => link2%next
                     end do
                  end if
               end if
         end select
         link => link%next
      end do

      if (allocated(processed_bulk))       deallocate(processed_bulk)
      if (allocated(processed_horizontal)) deallocate(processed_horizontal)
      if (allocated(processed_scalar))     deallocate(processed_bulk)

   end subroutine couple_standard_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process all model-specific coupling tasks.
!
! !INTERFACE:
   recursive subroutine process_coupling_tasks(self,process_model_references)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_base_model),intent(inout),target :: self
      logical,                intent(in)           :: process_model_references
!
!EOP
!
! !LOCAL VARIABLES:
      class (type_base_model),     pointer :: root
      type (type_model_list_node), pointer :: child
      type (type_link),            pointer :: link
      logical                              :: attempt_coupling
      class (type_internal_object),pointer :: master
!
!-----------------------------------------------------------------------
!BOC
      ! Find root model, which will handle the individual coupling tasks.
      root => self
      do while (associated(root%parent))
         root => root%parent
      end do

      ! Enumerate named couplings, locate slave and target variables, and couple them.
      link => self%first_link
      do while (associated(link))
         if (link%coupling%master/='') then
            select type (slave=>link%target)
               class is (type_model_reference)
                  attempt_coupling = process_model_references
               class default
                  attempt_coupling = .not.process_model_references
            end select
            if (attempt_coupling) then
               ! Try to find master variable.
               if (associated(link%coupling%source)) then
                  ! Try to find the master variable among the links of the specified model.
                  master => link%coupling%source%p%find_object(link%coupling%master)
               else
                  ! Try to find the master variable among the variables of the requesting model or its parents.
                  if (link%name/=link%coupling%master) then
                     ! Master and slave name differ: start master search in current model, then move up tree.
                     master => self%find_object(link%coupling%master,recursive=.true.)
                  elseif (associated(self%parent)) then
                     ! Master and slave name are identical: start master search in parent model, then move up tree.
                     master => self%parent%find_object(link%coupling%master,recursive=.true.)
                  end if
               end if
               if (associated(master)) then
                  ! Target variable found: perform the coupling.
                  call couple_variables(root,master,link%target)
               elseif (link%coupling%required.and..not.process_model_references) then
                  call self%fatal_error('process_coupling_tasks','Coupling target "'//trim(link%coupling%master)//'" for "'//trim(link%name)//'" was not found.')
               end if
            end if
         end if   ! if coupling
         link => link%next
      end do

      ! Process coupling tasks registered with child models.
      child => self%children%first
      do while (associated(child))
         call process_coupling_tasks(child%model,process_model_references)
         child => child%next
      end do

   end subroutine process_coupling_tasks
!EOC

recursive subroutine resolve_model_references(self)
   class (type_base_model),intent(inout),target :: self

   type (type_link),            pointer :: link
   class (type_base_model),     pointer :: model
   type (type_model_list_node), pointer :: child
   integer                              :: i

   ! Note: this routine MUST be called recursively on each model,
   ! because links can only be resolved locally (based on their local name).
   ! The root-level link name is not the same as the name of the model that must be coupled to.

   ! Check our own coupling commands.
   link => self%first_link
   do while (associated(link))
      if (link%coupling%master/=''.and.link%owner) then
         select type (slave=>link%target)
            class is (type_model_reference)
               model => self%find_model(link%coupling%master,recursive=.true.)
               if (link%coupling%required.and..not.associated(model)) &
                  call self%fatal_error('resolve_model_references','Target model "'//trim(link%coupling%master)//'" for "'//trim(link%name)//'" was not found.')
               do i=1,size(slave%pointers)
                  slave%pointers(i)%p%p => model
               end do
         end select
      end if
      link => link%next
   end do

   ! Now process child models
   child => self%children%first
   do while (associated(child))
      call resolve_model_references(child%model)
      child => child%next
   end do
end subroutine resolve_model_references

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

   call driver%log_message(trim(slave%name)//' --> '//trim(master%name))
   select type (master)
      class is (type_bulk_variable)
         select type (slave)
            class is (type_bulk_variable)
               call merge_bulk_variables(master,slave)
            class default
               call driver%fatal_error('merge_variables', &
                  'type mismatch: '//trim(master%name)//' is defined on the whole domain, '//trim(slave%name)//' is not.')
         end select      
      class is (type_horizontal_variable)
         select type (slave)
            class is (type_horizontal_variable)
               call merge_horizontal_variables(master,slave)
            class default
               call driver%fatal_error('merge_variables', &
                  'type mismatch: '//trim(master%name)//' is defined on the horizontal domain, '//trim(slave%name)//' is not.')
         end select      
      class is (type_scalar_variable)
         select type (slave)
            class is (type_scalar_variable)
               call merge_scalar_variables(master,slave)
            class default
               call driver%fatal_error('merge_variables', &
                  'type mismatch: '//trim(master%name)//' is defined as a scalar, '//trim(slave%name)//' is not.')
         end select      
      class is (type_model_reference)
         select type (slave)
            class is (type_model_reference)
               call merge_model_references(master,slave)
            class default
               call driver%fatal_error('merge_variables', &
                  'type mismatch: '//trim(master%name)//' is as model reference, '//trim(slave%name)//' is not.')
         end select      
   end select
end subroutine

subroutine merge_internal_variables(master,slave)
   class (type_internal_variable),intent(inout) :: master
   class (type_internal_variable),intent(in)    :: slave

   integer :: i
   type (type_contribution), pointer :: contribution

   if (allocated(slave%write_indices)) &
      call driver%fatal_error('merge_internal_variables','Attempt to couple write-only variable ' &
         //trim(slave%name)//' to '//trim(master%name)//'.')
   if (allocated(slave%state_indices).and..not.allocated(master%state_indices)) &
      call driver%fatal_error('merge_internal_variables','Attempt to couple state variable ' &
         //trim(slave%name)//' to non-state variable '//trim(master%name)//'.')
   if (allocated(slave%cons_indices).and..not.allocated(master%cons_indices)) &
      call driver%fatal_error('merge_internal_variables','Attempt to couple conserved quantity ' &
         //trim(slave%name)//' with non-conserved quantity '//trim(master%name)//'.')
   if (allocated(master%cons_indices).and..not.allocated(slave%cons_indices)) &
      call driver%fatal_error('merge_internal_variables','Attempt to couple non-conserved quantity ' &
         //trim(slave%name)//' with conserved quantity '//trim(master%name)//'.')
   if (master%presence==presence_external_optional) &
      call driver%fatal_error('merge_internal_variables','Attempt to couple to optional master variable "'//trim(master%name)//'".')

   if (allocated(slave%state_indices)) then
      do i=1,size(slave%state_indices)
         call append_index(master%state_indices,slave%state_indices(i)%p)
      end do
   end if
   if (allocated(slave%cons_indices)) then
      do i=1,size(slave%cons_indices)
         call append_index(master%cons_indices,slave%cons_indices(i)%p)
      end do
   end if
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
      call driver%fatal_error('merge_horizontal_variables','Domains of coupled variabled ' &
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

subroutine merge_model_references(master,slave)
   type (type_model_reference),intent(inout) :: master
   type (type_model_reference),intent(in)    :: slave
   integer :: i

   do i=1,size(slave%pointers)
      call append_model_pointer(master%pointers,slave%pointers(i)%p)
   end do
end subroutine merge_model_references

   function find_object(self,name,recursive) result(object)

      class (type_base_model),  intent(in),target :: self
      character(len=*),         intent(in)        :: name
      logical,         optional,intent(in)        :: recursive
      class (type_internal_object),pointer        :: object

      type (type_link), pointer :: link

      nullify(object)
      link => self%find_link(name,recursive)
      if (associated(link)) object => link%target
      
   end function find_object

   recursive function find_link(self,name,recursive) result(link)

      class (type_base_model),  intent(in),target :: self
      character(len=*),         intent(in)        :: name
      logical,         optional,intent(in)        :: recursive
      type (type_link),       pointer             :: link

      logical :: recursive_eff

      link => self%first_link
      do while (associated(link))
         if (link%name==name) return
         link => link%next
      end do

      ! Object not found - try parent if allowed.
      recursive_eff = .false.
      if (present(recursive)) recursive_eff = recursive
      if (recursive_eff.and.associated(self%parent)) link => self%parent%find_link(name,recursive)

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

function create_external_bulk_id(model,variable) result(id)
   class (type_base_model),intent(in),   target :: model
   type (type_bulk_variable),intent(inout),target :: variable
   type (type_bulk_variable_id) :: id
   id%variable => variable
   if (allocated(variable%alldata)) then
      allocate(id%alldata(size(variable%alldata)))
      id%alldata = variable%alldata
      id%p => id%alldata(1)%p
   else
      allocate(id%alldata(0))
   end if
   if (allocated(variable%state_indices)) id%state_index = variable%state_indices(1)%p
   if (allocated(variable%write_indices)) id%write_index = variable%write_indices(1)%p
end function create_external_bulk_id

function create_external_horizontal_id(model,variable) result(id)
   class (type_base_model),      intent(in),   target :: model
   type (type_horizontal_variable),intent(inout),target :: variable
   type (type_horizontal_variable_id) :: id
   id%variable => variable
   if (allocated(variable%alldata)) then
      allocate(id%alldata(size(variable%alldata)))
      id%alldata = variable%alldata
      id%p => id%alldata(1)%p
   else
      allocate(id%alldata(0))
   end if
   if (allocated(variable%state_indices)) id%state_index = variable%state_indices(1)%p
   if (allocated(variable%write_indices)) id%write_index = variable%write_indices(1)%p
end function create_external_horizontal_id

function create_external_scalar_id(model,variable) result(id)
   class (type_base_model),  intent(in),   target :: model
   type (type_scalar_variable),intent(inout),target :: variable
   type (type_scalar_variable_id) :: id
   id%variable => variable
   if (allocated(variable%alldata)) then
      allocate(id%alldata(size(variable%alldata)))
      id%alldata = variable%alldata
      id%p => id%alldata(1)%p
   else
      allocate(id%alldata(0))
   end if
   if (allocated(variable%state_indices)) id%state_index = variable%state_indices(1)%p
   if (allocated(variable%write_indices)) id%write_index = variable%write_indices(1)%p
end function create_external_scalar_id

function compare_bulk_standard_variables(variable1,variable2) result(equal)
   type (type_bulk_standard_variable),intent(in) :: variable1,variable2
   logical :: equal

   equal = .false.
   if (variable1%name /='' .and. variable2%name /='' .and. variable1%name /=variable2%name ) return
   if (variable1%units/='' .and. variable2%units/='' .and. variable1%units/=variable2%units) return
   equal = .true.
end function compare_bulk_standard_variables

function compare_horizontal_standard_variables(variable1,variable2) result(equal)
   type (type_horizontal_standard_variable),intent(in) :: variable1,variable2
   logical :: equal

   equal = .false.
   if (variable1%name /='' .and. variable2%name /='' .and. variable1%name /=variable2%name ) return
   if (variable1%units/='' .and. variable2%units/='' .and. variable1%units/=variable2%units) return
   equal = .true.
end function compare_horizontal_standard_variables

function compare_global_standard_variables(variable1,variable2) result(equal)
   type (type_global_standard_variable),intent(in) :: variable1,variable2
   logical :: equal

   equal = .false.
   if (variable1%name /='' .and. variable2%name /='' .and. variable1%name /=variable2%name ) return
   if (variable1%units/='' .and. variable2%units/='' .and. variable1%units/=variable2%units) return
   equal = .true.
end function compare_global_standard_variables

function is_null_bulk_standard_variable(variable) result(isnull)
   type (type_bulk_standard_variable),intent(in) :: variable
   logical :: isnull

   isnull = (variable%name==''.and. variable%units=='')
end function

function is_null_horizontal_standard_variable(variable) result(isnull)
   type (type_horizontal_standard_variable),intent(in) :: variable
   logical :: isnull

   isnull = (variable%name==''.and. variable%units=='')
end function

function is_null_global_standard_variable(variable) result(isnull)
   type (type_global_standard_variable),intent(in) :: variable
   logical :: isnull

   isnull = (variable%name==''.and. variable%units=='')
end function

function create_external_bulk_id_for_standard_name(model,standard_variable) result(id)
   class (type_base_model),           intent(in) :: model
   type (type_bulk_standard_variable),intent(in) :: standard_variable
   type (type_bulk_variable_id)                  :: id

   type (type_link), pointer :: link
   integer                   :: i

   allocate(id%alldata(0))
   link => model%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_bulk_variable)
            if (compare_standard_variables(object%standard_variable,standard_variable).and. &
                .not.is_null_standard_variable(object%standard_variable)) then
               if (allocated(object%alldata)) then
                  do i=1,size(object%alldata)
                     call append_data_pointer(id%alldata,object%alldata(i)%p)
                  end do
               end if
               id%variable => object
            end if
      end select    
      link => link%next
   end do
   if (size(id%alldata)>0) id%p => id%alldata(1)%p
end function create_external_bulk_id_for_standard_name

function create_external_horizontal_id_for_standard_name(model,standard_variable) result(id)
   class (type_base_model),                 intent(in) :: model
   type (type_horizontal_standard_variable),intent(in) :: standard_variable
   type (type_horizontal_variable_id)                  :: id

   type (type_link), pointer :: link
   integer                   :: i

   allocate(id%alldata(0))
   link => model%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_horizontal_variable)
            if (compare_standard_variables(object%standard_variable,standard_variable).and. &
                .not.is_null_standard_variable(object%standard_variable)) then
               if (allocated(object%alldata)) then
                  do i=1,size(object%alldata)
                     call append_data_pointer(id%alldata,object%alldata(i)%p)
                  end do
               end if
               id%variable => object
            end if
      end select
      link => link%next
   end do
   if (size(id%alldata)>0) id%p => id%alldata(1)%p
end function create_external_horizontal_id_for_standard_name

function create_external_scalar_id_for_standard_name(model,standard_variable) result(id)
   class (type_base_model),             intent(in) :: model
   type (type_global_standard_variable),intent(in) :: standard_variable
   type (type_scalar_variable_id)                  :: id

   type (type_link), pointer :: link
   integer                   :: i

   allocate(id%alldata(0))
   link => model%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_scalar_variable)
            if (compare_standard_variables(object%standard_variable,standard_variable).and. &
                .not.is_null_standard_variable(object%standard_variable)) then
               if (allocated(object%alldata)) then
                  do i=1,size(object%alldata)
                     call append_data_pointer(id%alldata,object%alldata(i)%p)
                  end do
               end if
               id%variable => object
            end if
      end select
      link => link%next
   end do
   if (size(id%alldata)>0) id%p => id%alldata(1)%p
end function create_external_scalar_id_for_standard_name

recursive subroutine build_aggregate_variables(self)
   class (type_base_model),intent(inout),target :: self

   type (type_aggregate_variable),pointer :: aggregate_variable
   type (type_model_list_node),   pointer :: child
   type (type_link),              pointer :: link
   type (type_contribution),      pointer :: contribution

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
               call add_contribution(aggregate_variable,link,contribution%scale_factor)
               contribution => contribution%next
            end do
      end select
      link => link%next
   end do

   if (self%check_conservation) then
      aggregate_variable => self%first_aggregate_variable
      do while (associated(aggregate_variable))
         if (aggregate_variable%has_bulk_state_component)       call self%register_diagnostic_variable(aggregate_variable%id_rate,'change_in_'//trim(aggregate_variable%standard_variable%name),trim(aggregate_variable%standard_variable%units),'change in '//trim(aggregate_variable%standard_variable%name))
         if (aggregate_variable%has_horizontal_state_component) call self%register_diagnostic_variable(aggregate_variable%id_horizontal_rate,'change_in_'//trim(aggregate_variable%standard_variable%name)//'_at_interfaces',trim(aggregate_variable%standard_variable%units),'change in '//trim(aggregate_variable%standard_variable%name)//' at interfaces')
         aggregate_variable => aggregate_variable%next
      end do
   end if
   
   ! Process child models
   child => self%children%first
   do while (associated(child))
      call build_aggregate_variables(child%model)
      child => child%next
   end do

contains

   subroutine add_contribution(aggregate_variable,link,scale_factor)
      type (type_aggregate_variable),intent(inout) :: aggregate_variable
      type (type_link),target,       intent(in)    :: link
      real(rk),                      intent(in)    :: scale_factor

      type (type_contributing_variable),pointer :: contributing_variable

      if (.not.associated(aggregate_variable%first_contributing_variable)) then
         ! This aggregate variable does not have any contribution registered yet. Create the first.
         allocate(aggregate_variable%first_contributing_variable)
         contributing_variable => aggregate_variable%first_contributing_variable
      else
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

      ! If the contributing variable is a state variable, keep a reference to its state variable index.
      select type (variable=>link%target)
         class is (type_bulk_variable)
            if (allocated(variable%state_indices)) then
               call append_index(variable%state_indices,contributing_variable%state_index)
               aggregate_variable%has_bulk_state_component = .true.
            end if
         class is (type_horizontal_variable)
            if (allocated(variable%state_indices)) then
               call append_index(variable%state_indices,contributing_variable%horizontal_state_index)
               aggregate_variable%has_horizontal_state_component = .true.
            end if
      end select
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
         if (compare_standard_variables(aggregate_variable%standard_variable,standard_variable)) return
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
end function

subroutine create_aggregate_model(self,aggregate_variable)
   class (type_base_model),       intent(inout),target :: self
   type (type_aggregate_variable),intent(inout)        :: aggregate_variable

   type (type_contributing_variable),pointer :: contributing_variable

   ! This procedure takes an aggregate variable, and creates models that compute diagnostics
   ! for the total of these aggregate quantities on bulk and horizontal domains.
   
   contributing_variable => aggregate_variable%first_contributing_variable
   do while (associated(contributing_variable))
      if (contributing_variable%link%owner) then
         select type (variable=>contributing_variable%link%target)
            class is (type_bulk_variable)
               ! This contribution comes from a bulk variable.
               if (.not.associated(aggregate_variable%sum)) then
                  allocate(aggregate_variable%sum)
                  aggregate_variable%sum%output_name = trim(aggregate_variable%standard_variable%name)
                  aggregate_variable%sum%output_units = trim(aggregate_variable%standard_variable%units)
               end if
               call aggregate_variable%sum%add_component(trim(contributing_variable%link%name),weight=contributing_variable%scale_factor)
            class is (type_horizontal_variable)
               ! This contribution comes from a variable defined on a horizontal interface (top or bottom).
               if (.not.associated(aggregate_variable%horizontal_sum)) then
                  allocate(aggregate_variable%horizontal_sum)
                  aggregate_variable%horizontal_sum%output_name = trim(aggregate_variable%standard_variable%name)//'_at_interfaces'
                  aggregate_variable%horizontal_sum%output_units = trim(aggregate_variable%standard_variable%units)
               end if
               call aggregate_variable%horizontal_sum%add_component(trim(contributing_variable%link%name),weight=contributing_variable%scale_factor)
         end select
      end if
      contributing_variable => contributing_variable%next
   end do
   if (associated(aggregate_variable%sum))            call self%add_child(aggregate_variable%sum,trim(aggregate_variable%sum%output_name)//'_calculator',configunit=-1)
   if (associated(aggregate_variable%horizontal_sum)) call self%add_child(aggregate_variable%horizontal_sum,trim(aggregate_variable%horizontal_sum%output_name)//'_calculator',configunit=-1)
end subroutine create_aggregate_model

recursive subroutine classify_variables(self)
   class (type_base_model),intent(inout),target :: self

   type (type_link),pointer :: link

   type (type_state_variable_info),                pointer :: statevar
   type (type_horizontal_state_variable_info),     pointer :: hz_statevar
   type (type_diagnostic_variable_info),           pointer :: diagvar
   type (type_horizontal_diagnostic_variable_info),pointer :: hz_diagvar
   type (type_conserved_quantity_info),            pointer :: consvar
   integer                                                 :: nstate,nstate_bot,nstate_surf,ndiag,ndiag_hz,ncons

   integer :: i
   logical :: set_indices

   type (type_model_list_node),       pointer :: child
   type (type_aggregate_variable),    pointer :: aggregate_variable
   type (type_contributing_variable), pointer :: contributing_variable

   if (self%frozen) call self%fatal_error('classify_variables', &
      'classify_variables may only be called once (during initialization).')

   ! Make sure that no one will register any new variables from this moment onward.
   self%frozen = .true.

   ! Count number of bulk variables in various categories.
   set_indices = .not.associated(self%parent)
   nstate = 0
   ndiag  = 0
   nstate_bot  = 0
   nstate_surf = 0
   ndiag_hz    = 0
   link => self%first_link
   do while (associated(link))
      if (link%owner) then
         select type (object => link%target)
            class is (type_bulk_variable)
               ! This is a variable owned by the model.
               if (allocated(object%write_indices)) then
                  ! This is a diagnostic variable.
                  ndiag = ndiag+1
                  do i=1,size(object%write_indices)
                     if (set_indices) object%write_indices(i)%p = ndiag
                  end do
               end if
               if (allocated(object%state_indices)) then
                  ! This is a state variable.
                  select case (object%presence)
                     case (presence_internal)
                        nstate = nstate+1
                        do i=1,size(object%state_indices)
                           if (set_indices) object%state_indices(i)%p = nstate
                        end do
                     case (presence_external_required)
                        call self%fatal_error('classify_variables','Variable "'//trim(link%name)//'" must be coupled to an existing state variable.')
                     case default
                        continue
                  end select
               end if
            class is (type_horizontal_variable)
               if (allocated(object%write_indices)) then
                  ! This is a diagnostic variable.
                  ndiag_hz = ndiag_hz+1
                  do i=1,size(object%write_indices)
                     if (set_indices) object%write_indices(i)%p = ndiag_hz
                  end do
               end if
               if (allocated(object%state_indices)) then
                  ! This is a state variable.
                  select case (object%presence)
                     case (presence_internal)
                        select case (object%domain)
                           case (domain_bottom)
                              ! Bottom state variable
                              nstate_bot = nstate_bot+1
                              do i=1,size(object%state_indices)
                                 if (set_indices) object%state_indices(i)%p = nstate_bot
                              end do
                           case (domain_surface)
                              ! Surface state variable
                              nstate_surf = nstate_surf+1
                              do i=1,size(object%state_indices)
                                 if (set_indices) object%state_indices(i)%p = nstate_surf
                              end do
                        end select
                     case (presence_external_required)
                        call self%fatal_error('classify_variables','Variable "'//trim(link%name)//'" must be coupled to an existing state variable.')
                  end select
               end if
         end select
      end if
      link => link%next
   end do

   ncons = 0
   aggregate_variable => self%first_aggregate_variable
   do while (associated(aggregate_variable))
      ncons = ncons + 1
      aggregate_variable => aggregate_variable%next
   end do

   ! Allocate arrays with variable information that will be accessed by the host model.
   allocate(self%state_variables                (nstate))
   allocate(self%bottom_state_variables         (nstate_bot))
   allocate(self%surface_state_variables        (nstate_surf))
   allocate(self%diagnostic_variables           (ndiag))
   allocate(self%horizontal_diagnostic_variables(ndiag_hz))
   allocate(self%conserved_quantities           (ncons))

   ! Create empty lists that will hold the names of all readable variables.
   allocate(self%dependencies(0))
   allocate(self%dependencies_hz(0))
   allocate(self%dependencies_scalar(0))

   ! Set pointers for backward compatibility (pre 2013-06-15)
   ! Note: this MUST be doen after allocation of the target arrays, above!
   self%state_variables_ben => self%bottom_state_variables
   self%diagnostic_variables_hz => self%horizontal_diagnostic_variables

   ! Copy metadata of bulk variables
   nstate = 0
   ndiag  = 0
   ncons  = 0
   nstate_bot  = 0
   nstate_surf = 0
   ndiag_hz    = 0
   link => self%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_bulk_variable)
            if (link%owner) then
               ! The model owns this variable (no external master variable has been assigned)
               ! Transfer variable information to the array that will be accessed by the host model.
               if (allocated(object%write_indices)) then
                  ndiag = ndiag + 1
                  diagvar => self%diagnostic_variables(ndiag)
                  call copy_variable_metadata(object,diagvar)
                  diagvar%globalid          = create_external_variable_id(self,object)
                  diagvar%standard_variable = object%standard_variable
                  diagvar%time_treatment    = object%time_treatment
                  call diagvar%properties%update(object%properties)
               end if
         
               if (allocated(object%state_indices).and.object%presence==presence_internal) then
                  nstate = nstate + 1
                  statevar => self%state_variables(nstate)
                  call copy_variable_metadata(object,statevar)
                  statevar%globalid                  = create_external_variable_id(self,object)
                  statevar%standard_variable         = object%standard_variable
                  statevar%initial_value             = object%initial_value
                  statevar%vertical_movement         = object%vertical_movement
                  statevar%specific_light_extinction = object%specific_light_extinction
                  statevar%no_precipitation_dilution = object%no_precipitation_dilution
                  statevar%no_river_dilution         = object%no_river_dilution
               end if
            end if
            if (allocated(object%alldata).and..not.(allocated(object%state_indices).and.object%presence==presence_external_optional)) then
               call append_string(self%dependencies,link%name)
               if (object%standard_variable%name/='') call append_string(self%dependencies,object%standard_variable%name)
            end if
         class is (type_horizontal_variable)
            if (link%owner) then
               ! The model owns this variable (no external master variable has been assigned)
               ! Transfer variable information to the array that will be accessed by the host model.
               if (allocated(object%write_indices)) then
                  ndiag_hz = ndiag_hz + 1
                  hz_diagvar => self%horizontal_diagnostic_variables(ndiag_hz)
                  call copy_variable_metadata(object,hz_diagvar)
                  hz_diagvar%globalid          = create_external_variable_id(self,object)
                  hz_diagvar%standard_variable = object%standard_variable
                  hz_diagvar%time_treatment    = object%time_treatment
               end if
               if (allocated(object%state_indices).and.object%presence==presence_internal) then
                  select case (object%domain)
                     case (domain_bottom)
                        nstate_bot = nstate_bot + 1
                        hz_statevar => self%bottom_state_variables(nstate_bot)
                     case (domain_surface)
                        nstate_surf = nstate_surf + 1
                        hz_statevar => self%surface_state_variables(nstate_surf)
                     case default
                        nullify(hz_statevar)
                  end select
                  call copy_variable_metadata(object,hz_statevar)
                  hz_statevar%globalid          = create_external_variable_id(self,object)
                  hz_statevar%standard_variable = object%standard_variable
                  hz_statevar%initial_value     = object%initial_value
               end if
            end if
            if (allocated(object%alldata).and..not.(allocated(object%state_indices).and.object%presence==presence_external_optional)) then
               call append_string(self%dependencies_hz,link%name)
               if (object%standard_variable%name/='') call append_string(self%dependencies_hz,object%standard_variable%name)
            end if
         class is (type_scalar_variable)
            if (allocated(object%alldata)) then
               call append_string(self%dependencies_scalar,link%name)
               if (object%standard_variable%name/='') call append_string(self%dependencies_scalar,object%standard_variable%name)
            end if
      end select
      link => link%next
   end do

   ncons = 0
   aggregate_variable => self%first_aggregate_variable
   do while (associated(aggregate_variable))
      ncons = ncons + 1
      consvar => self%conserved_quantities(ncons)
      consvar%standard_variable = aggregate_variable%standard_variable
      consvar%name = trim(consvar%standard_variable%name)
      consvar%units = trim(consvar%standard_variable%units)
      consvar%long_name = trim(consvar%standard_variable%name)
      consvar%model => aggregate_variable%sum
      consvar%horizontal_model => aggregate_variable%horizontal_sum

      ! Loop over all contributing quantities, and if they are bulk state variables,
      ! store their state variable index and scale factor, so we can include them later
      ! when computing the (change in) conserved quantities from the (change in) state alone.

      nstate = 0
      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (contributing_variable%state_index/=-1) nstate = nstate+1
         contributing_variable => contributing_variable%next
      end do

      allocate(consvar%state_indices(nstate))
      allocate(consvar%state_scale_factors(nstate))
      nstate = 0
      contributing_variable => aggregate_variable%first_contributing_variable
      do while (associated(contributing_variable))
         if (contributing_variable%state_index/=-1) then
            nstate = nstate+1
            consvar%state_indices(nstate) = contributing_variable%state_index
            consvar%state_scale_factors(nstate) = contributing_variable%scale_factor
         end if
         contributing_variable => contributing_variable%next
      end do
      call remove_duplicates(consvar%state_indices,consvar%state_scale_factors)

      if (self%check_conservation) then
         consvar%rate_diag_index = aggregate_variable%id_rate%diag_index
         consvar%horizontal_rate_diag_index = aggregate_variable%id_horizontal_rate%horizontal_diag_index
      end if
      
      aggregate_variable => aggregate_variable%next
   end do

   ! Process child models
   child => self%children%first
   do while (associated(child))
      call classify_variables(child%model)
      child => child%next
   end do

end subroutine classify_variables

subroutine remove_duplicates(array1,array2)
   integer, allocatable,intent(inout) :: array1(:)
   real(rk),allocatable,intent(inout) :: array2(:)

   integer  :: i,j,inext
   logical  :: add
   integer  :: array1_tmp(size(array1))
   real(rk) :: array2_tmp(size(array2))

   inext = 0
   do i=1,size(array1)
      add = array1(i)/=-1
      do j=1,i-1
         if (array1(j)==array1(i)) add = .false.
      end do
      if (add) then
         inext = inext + 1
         array1_tmp(inext) = array1(i)
         array2_tmp(inext) = array2(i)
      end if
   end do
   deallocate(array1)
   deallocate(array2)
   allocate(array1(inext))
   allocate(array2(inext))
   array1 = array1_tmp(:inext)
   array2 = array2_tmp(:inext)
end subroutine

subroutine complete_model_references(model)
   class (type_base_model),intent(inout),target :: model

   type (type_link),pointer :: link
   integer                  :: i

   link => model%first_link
   do while (associated(link))
      if (link%owner) then
         select type (object => link%target)
            class is (type_model_reference)
               do i=1,size(object%pointers)
                  call process_model_pointer(object%pointers(i)%p)
               end do
         end select
      end if
      link => link%next
   end do

contains

   subroutine process_model_pointer(p)
      type (type_model_pointer),intent(inout) :: p

      integer                  :: n
      type (type_link),pointer :: link

      ! Count number of state variables in target model.
      n = 0
      link => p%p%first_link
      do while (associated(link))
         select type (object=>link%target)
            class is (type_bulk_variable)
               if (allocated(object%state_indices).and.object%presence==presence_internal.and.link%owner) n = n + 1
         end select
         link => link%next
      end do

      ! Allocate array to hold state variable identifiers.
      allocate(p%state(n))

      ! Connect target state variables to identifiers in model reference.
      n = 0
      link => p%p%first_link
      do while (associated(link))
         select type (object=>link%target)
            class is (type_bulk_variable)
               if (allocated(object%state_indices).and.object%presence==presence_internal.and.link%owner) then
                  n = n + 1
                  call connect_bulk_state_variable_id(model,object,p%state(n))
               end if
         end select
         link => link%next
      end do
   end subroutine process_model_pointer

end subroutine complete_model_references

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

subroutine copy_variable_metadata(internal_variable,external_variable)
   class (type_external_variable),intent(inout) :: external_variable
   class (type_internal_variable),intent(in)    :: internal_variable

   external_variable%name          = internal_variable%name
   external_variable%long_name     = internal_variable%long_name
   external_variable%units         = internal_variable%units
   external_variable%minimum       = internal_variable%minimum
   external_variable%maximum       = internal_variable%maximum
   external_variable%missing_value = internal_variable%missing_value
   call external_variable%properties%update(internal_variable%properties)
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
         call driver%fatal_error('find_dependencies','circular dependency found: '//trim(chain)//trim(self%name))
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

   subroutine base_driver_fatal_error(self,location,message)
      class (type_base_driver), intent(inout) :: self
      character(len=*),         intent(in)    :: location,message

      write (*,*) trim(location)//': '//trim(message)
      stop 1
   end subroutine

   subroutine base_driver_log_message(self,message)
      class (type_base_driver), intent(inout) :: self
      character(len=*),         intent(in)    :: message

      write (*,*) trim(message)
   end subroutine

   subroutine weighted_sum_initialize(self,configunit)
      class (type_weighted_sum),intent(inout),target :: self
      integer,                  intent(in)           :: configunit

      type (type_component),pointer :: component
      component => self%first
      do while (associated(component))
         call self%register_dependency(component%id,trim(component%name))
         component => component%next
      end do
      if (self%output_long_name=='') self%output_long_name = self%output_name
      call self%register_diagnostic_variable(self%id_output,self%output_name,self%output_units,self%output_long_name)
      call self%parent%add_alias(self%id_output,trim(self%output_name))
   end subroutine

   subroutine weighted_sum_add_component(self,name,weight)
      class (type_weighted_sum),intent(inout) :: self
      character(len=*),         intent(in)    :: name
      real(rk),optional,        intent(in)    :: weight

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
   end subroutine

   subroutine weighted_sum_evaluate(self,_ARGUMENTS_ND_)
      class (type_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_ND_

      type (type_component),pointer        :: component
      real(rk)                             :: value
      real(rk) _DIMENSION_SLICE_AUTOMATIC_ :: sum

      ! Initialize sum to zero.
      sum = 0.0_rk

      ! Enumerate components included in the sum, and add their contributions.
      component => self%first
      do while (associated(component))
         _LOOP_BEGIN_
            _GET_(component%id,value)
            sum _INDEX_OUTPUT_ = sum _INDEX_OUTPUT_ + component%weight*value
         _LOOP_END_
         component => component%next
      end do

      ! Transfer summed values to diagnostic.
      _LOOP_BEGIN_
         _SET_DIAGNOSTIC_(self%id_output,sum _INDEX_OUTPUT_)
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
      component => self%first
      do while (associated(component))
         call self%register_dependency(component%id,trim(component%name))
         component => component%next
      end do
      if (self%output_long_name=='') self%output_long_name = self%output_name
      call self%register_diagnostic_variable(self%id_output,self%output_name,self%output_units,self%output_long_name)
      call self%parent%add_alias(self%id_output,trim(self%output_name))
   end subroutine

   subroutine horizontal_weighted_sum_add_component(self,name,weight)
      class (type_horizontal_weighted_sum),intent(inout) :: self
      character(len=*),                    intent(in)    :: name
      real(rk),optional,                   intent(in)    :: weight

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
   end subroutine

   subroutine horizontal_weighted_sum_evaluate_horizontal(self,_ARGUMENTS_HZ_)
      class (type_horizontal_weighted_sum),intent(in) :: self
      _DECLARE_ARGUMENTS_HZ_
      
      type (type_horizontal_component),pointer :: component
      real(rk)                      :: sum,value

      _HORIZONTAL_LOOP_BEGIN_
         sum = 0._rk
         component => self%first
         do while (associated(component))
            _GET_HORIZONTAL_(component%id,value)
            sum = sum + component%weight*value
            component => component%next
         end do
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_output,sum)
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
      do i=1,size(self%conserved_quantities)
         do j=1,size(self%conserved_quantities(i)%state_indices)
            index = self%conserved_quantities(i)%state_indices(j)
            scale_factor = self%conserved_quantities(i)%state_scale_factors(j)
            _LOOP_BEGIN_
               sums _INDEX_OUTPUT_1D_(i) = sums _INDEX_OUTPUT_1D_(i) + scale_factor*y _INDEX_OUTPUT_1D_(index)
            _LOOP_END_
         end do
      end do
   end subroutine

   end module fabm_types

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

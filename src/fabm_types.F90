#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_types --- Derived types used by 0D biogeochemical modules
!
! !INTERFACE:
   module fabm_types
!
! !DESCRIPTION:
! This module contains the derived types that are used for communication between
! the 0D biogeochemical models and a hosting physical environment.
! Types are used to describe the 0d model, its state variables, and the local
! environment.
!
! Subroutines for intialization of the derived types are also provided. It is
! recommended that you always call these subroutines to initialize such types
! before use. This ensures that if new members are added to the types, they will
! be set to a reasonable default even if your program is not aware the member
! has been added.
!
! !USES:
   use fabm_driver, only: fatal_error, log_message
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
   public type_model_info

   ! Collection with standard variables (e.g., temeprature, practical_salinity)
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

   ! Data types and procedures for variable management - used by FABM internally only.
   public type_bulk_variable_link
   public type_horizontal_variable_link
   public type_scalar_variable_link
   public type_environment
   public initialize_model_info
   public freeze_model_info
   public create_external_variable_id

   public type_bulk_data_pointer
   public type_horizontal_data_pointer

   public type_conserved_quantity_component

#ifdef _FABM_F2003_
   public type_expression, type_bulk_expression, type_horizontal_expression
#endif

! !PUBLIC DATA MEMBERS:
!
   integer, parameter, public :: attribute_length = 256

   integer, parameter, public :: rk = _FABM_REAL_KIND_

   integer, parameter, public :: domain_bulk = 0, domain_bottom = 1, domain_surface = 2

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
      real(rk),pointer _ATTR_LOCATION_DIMENSIONS_ :: p => null()
   end type

   type type_horizontal_data_pointer
      real(rk),pointer _ATTR_LOCATION_DIMENSIONS_HZ_ :: p => null()
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
   end type type_bulk_data_pointer_pointer

   type type_horizontal_data_pointer_pointer
      type (type_horizontal_data_pointer),pointer :: p => null()
   end type type_horizontal_data_pointer_pointer

   type type_scalar_data_pointer_pointer
      type (type_scalar_data_pointer),pointer :: p => null()
   end type type_scalar_data_pointer_pointer

   ! ====================================================================================================
   ! Variable types used by FABM for both metadata and value pointers/indices.
   ! ====================================================================================================

   type type_conserved_quantity_component
      type (type_bulk_data_pointer)         :: state
      real(rk)                              :: scale_factor = 1.0_rk
      type (type_conserved_quantity_component),pointer :: next => null()
   end type

   type type_conserved_quantity_component_list
      type (type_conserved_quantity_component),pointer :: first => null()
   end type

#ifdef _FABM_F2003_
   type type_internal_variable
      type (type_property_dictionary)               :: properties
      type (type_conserved_quantity_component_list) :: components
   end type
#define _EXTENDS_INTERNAL_VARIABLE_ ,extends(type_internal_variable) ::
#else
#define _EXTENDS_INTERNAL_VARIABLE_
#endif

   type _EXTENDS_INTERNAL_VARIABLE_ type_bulk_variable
      ! Metadata
      character(len=attribute_length)    :: name                      = ''
      character(len=attribute_length)    :: long_name                 = ''
      character(len=attribute_length)    :: units                     = ''
      real(rk)                           :: minimum                   = -1.e20_rk
      real(rk)                           :: maximum                   =  1.e20_rk
      real(rk)                           :: missing_value             = -2.e20_rk
      real(rk)                           :: initial_value             = 0.0_rk
      real(rk)                           :: vertical_movement         = 0.0_rk
      real(rk)                           :: specific_light_extinction = 0.0_rk
      logical                            :: no_precipitation_dilution = .false.
      logical                            :: no_river_dilution         = .false.
      integer                            :: time_treatment            = time_treatment_last
      type (type_bulk_standard_variable) :: standard_variable

      ! Arrays with all associated data and index pointers.
      type (type_bulk_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata       _NULL_
      type (type_integer_pointer),          dimension(:),_ALLOCATABLE_ :: state_indices _NULL_
      type (type_integer_pointer),          dimension(:),_ALLOCATABLE_ :: write_indices _NULL_
      type (type_integer_pointer),          dimension(:),_ALLOCATABLE_ :: cons_indices  _NULL_
   end type type_bulk_variable

   type _EXTENDS_INTERNAL_VARIABLE_ type_horizontal_variable
      ! Metadata
      character(len=attribute_length) :: name           = ''
      character(len=attribute_length) :: long_name      = ''
      character(len=attribute_length) :: units          = ''
      real(rk)                        :: minimum        = -1.e20_rk
      real(rk)                        :: maximum        =  1.e20_rk
      real(rk)                        :: missing_value  = -2.e20_rk
      real(rk)                        :: initial_value  = 0.0_rk
      integer                         :: time_treatment = time_treatment_last
      integer                         :: domain         = domain_bottom
      type (type_horizontal_standard_variable) :: standard_variable
      
      ! Arrays with all associated data and index pointers.
      type (type_horizontal_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata       _NULL_
      type (type_integer_pointer),                dimension(:),_ALLOCATABLE_ :: state_indices _NULL_
      type (type_integer_pointer),                dimension(:),_ALLOCATABLE_ :: write_indices _NULL_
      type (type_integer_pointer),                dimension(:),_ALLOCATABLE_ :: cons_indices  _NULL_
   end type type_horizontal_variable

   type _EXTENDS_INTERNAL_VARIABLE_ type_scalar_variable
      ! Metadata
      character(len=attribute_length) :: name           = ''
      character(len=attribute_length) :: long_name      = ''
      character(len=attribute_length) :: units          = ''
      real(rk)                        :: minimum        = -1.e20_rk
      real(rk)                        :: maximum        =  1.e20_rk
      real(rk)                        :: missing_value  = -2.e20_rk
      integer                         :: time_treatment = time_treatment_last
      type (type_global_standard_variable) :: standard_variable

      ! Arrays with all associated data and index pointers.
      type (type_scalar_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata       _NULL_
      type (type_integer_pointer),            dimension(:),_ALLOCATABLE_ :: state_indices _NULL_
      type (type_integer_pointer),            dimension(:),_ALLOCATABLE_ :: write_indices _NULL_
      type (type_integer_pointer),            dimension(:),_ALLOCATABLE_ :: cons_indices  _NULL_
   end type type_scalar_variable

   type type_bulk_variable_link
      character(len=attribute_length)         :: name    = ''
      type (type_bulk_variable), pointer      :: target  => null()
      logical                                 :: coupled = .false.
      type (type_bulk_variable_link), pointer :: next    => null()
   end type type_bulk_variable_link

   type type_horizontal_variable_link
      character(len=attribute_length)               :: name    = ''
      type (type_horizontal_variable),      pointer :: target  => null()
      logical                                       :: coupled = .false.
      type (type_horizontal_variable_link), pointer :: next    => null()
   end type type_horizontal_variable_link

   type type_scalar_variable_link
      character(len=attribute_length)           :: name    = ''
      type (type_scalar_variable),      pointer :: target  => null()
      logical                                   :: coupled = .false.
      type (type_scalar_variable_link), pointer :: next    => null()
   end type type_scalar_variable_link

   type type_named_coupling
      character(len=attribute_length)     :: master = ''
      character(len=attribute_length)     :: slave  = ''
      type (type_named_coupling), pointer :: next   => null()
   end type

   type type_bulk_variable_id
      type (type_bulk_variable),    pointer                            :: variable => null()
      type (type_bulk_data_pointer),pointer                            :: p        => null()
      type (type_bulk_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata _NULL_
      integer                                                          :: state_index = -1
      integer                                                          :: write_index = -1
   end type

   type type_horizontal_variable_id
      type (type_horizontal_variable),    pointer                            :: variable => null()
      type (type_horizontal_data_pointer),pointer                            :: p        => null()
      type (type_horizontal_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata _NULL_
      integer                                                                :: state_index = -1
      integer                                                                :: write_index = -1
   end type

   type type_scalar_variable_id
      type (type_scalar_variable),    pointer                            :: variable => null()
      type (type_scalar_data_pointer),pointer                            :: p        => null()
      type (type_scalar_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata _NULL_
      integer                                                            :: state_index = -1
      integer                                                            :: write_index = -1
   end type

   ! ====================================================================================================
   ! Variable identifiers used by biogeochemical models.
   ! ====================================================================================================

#ifdef _FABM_F2003_
   type type_id
      class (type_internal_variable), pointer :: metadata => null()
   end type
#define _EXTENDS_ID_ ,extends(type_id) ::
#else
#define _EXTENDS_ID_
#endif

   type _EXTENDS_ID_ type_state_variable_id
      character(len=attribute_length) :: name        = ''
      integer                         :: state_index = -1
      type (type_bulk_data_pointer)   :: data
   end type

   type _EXTENDS_ID_ type_bottom_state_variable_id
      character(len=attribute_length)     :: name               = ''
      integer                             :: bottom_state_index = -1
      type (type_horizontal_data_pointer) :: horizontal_data
   end type

   type _EXTENDS_ID_ type_surface_state_variable_id
      character(len=attribute_length)     :: name               = ''
      integer                             :: surface_state_index = -1
      type (type_horizontal_data_pointer) :: horizontal_data
   end type

   type _EXTENDS_ID_ type_diagnostic_variable_id
      character(len=attribute_length) :: name       = ''
      integer                         :: diag_index = -1
   end type

   type _EXTENDS_ID_ type_horizontal_diagnostic_variable_id
      character(len=attribute_length) :: name                  = ''
      integer                         :: horizontal_diag_index = -1
   end type

   type type_dependency_id
      character(len=attribute_length) :: name = ''
      type (type_bulk_data_pointer)   :: data
   end type

   type type_horizontal_dependency_id
      character(len=attribute_length)     :: name = ''
      type (type_horizontal_data_pointer) :: horizontal_data
   end type

   type type_global_dependency_id
      character(len=attribute_length) :: name = ''
      type (type_scalar_data_pointer) :: global_data
   end type

   type _EXTENDS_ID_ type_conserved_quantity_id
      character(len=attribute_length) :: name       = ''
      integer                         :: cons_index = -1
   end type

   ! ====================================================================================================
   ! Types to hold variable metadata, used by the external host.
   ! ====================================================================================================

#ifdef _FABM_F2003_
   type type_external_variable
      type (type_property_dictionary) :: properties
   end type
#define _EXTENDS_EXTERNAL_VARIABLE_ ,extends(type_external_variable) ::
#else
#define _EXTENDS_EXTERNAL_VARIABLE_
#endif

!  Derived type describing a state variable
   type _EXTENDS_EXTERNAL_VARIABLE_ type_state_variable_info
      character(len=attribute_length)   :: name                      = ''
      character(len=attribute_length)   :: long_name                 = ''
      character(len=attribute_length)   :: units                     = ''
      type (type_bulk_standard_variable) :: standard_variable
      real(rk)                          :: initial_value             = 0.0_rk
      real(rk)                          :: minimum                   = -1.e20_rk
      real(rk)                          :: maximum                   =  1.e20_rk
      real(rk)                          :: missing_value             = -2.e20_rk
      real(rk)                          :: vertical_movement         = 0.0_rk  ! Vertical movement (m/s) due to sinking (<0), floating (>0).
      real(rk)                          :: specific_light_extinction = 0.0_rk  ! Specific light extinction (/m/state variable unit)
      logical                           :: no_precipitation_dilution = .false.
      logical                           :: no_river_dilution         = .false.
      integer                           :: externalid                = 0       ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      type (type_bulk_variable_id)      :: globalid
   end type type_state_variable_info

   type _EXTENDS_EXTERNAL_VARIABLE_ type_horizontal_state_variable_info
      character(len=attribute_length)         :: name          = ''
      character(len=attribute_length)         :: long_name     = ''
      character(len=attribute_length)         :: units         = ''
      type (type_horizontal_standard_variable) :: standard_variable
      real(rk)                                :: initial_value = 0.0_rk
      real(rk)                                :: minimum       = -1.e20_rk
      real(rk)                                :: maximum       =  1.e20_rk
      real(rk)                                :: missing_value = -2.e20_rk
      integer                                 :: externalid    = 0      ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      type (type_horizontal_variable_id)      :: globalid
   end type type_horizontal_state_variable_info

!  Derived type describing a diagnostic variable
   type _EXTENDS_EXTERNAL_VARIABLE_ type_diagnostic_variable_info
      character(len=attribute_length)   :: name           = ''
      character(len=attribute_length)   :: long_name      = ''
      character(len=attribute_length)   :: units          = ''
      type (type_bulk_standard_variable) :: standard_variable
      real(rk)                          :: minimum        = -1.e20_rk
      real(rk)                          :: maximum        =  1.e20_rk
      real(rk)                          :: missing_value  = -2.e20_rk
      integer                           :: externalid     = 0                    ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      integer                           :: time_treatment = time_treatment_last ! Time treatment: 0=last value, 1=time-integrated, 2=time step-averaged, 3=time step-integrated
      type (type_bulk_variable_id)      :: globalid
   end type type_diagnostic_variable_info

   type _EXTENDS_EXTERNAL_VARIABLE_ type_horizontal_diagnostic_variable_info
      character(len=attribute_length)         :: name           = ''
      character(len=attribute_length)         :: long_name      = ''
      character(len=attribute_length)         :: units          = ''
      type (type_horizontal_standard_variable) :: standard_variable
      real(rk)                                :: minimum        = -1.e20_rk
      real(rk)                                :: maximum        =  1.e20_rk
      real(rk)                                :: missing_value  = -2.e20_rk
      integer                                 :: externalid     = 0                   ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      integer                                 :: time_treatment = time_treatment_last ! Time treatment: 0=last value, 1=time-integrated, 2=time step-averaged, 3=time step-integrated
      type (type_horizontal_variable_id)      :: globalid
   end type type_horizontal_diagnostic_variable_info

!  Derived type describing a conserved quantity
   type _EXTENDS_EXTERNAL_VARIABLE_ type_conserved_quantity_info
      character(len=attribute_length)   :: name       = ''
      character(len=attribute_length)   :: long_name  = ''
      character(len=attribute_length)   :: units      = ''
      type (type_bulk_standard_variable) :: standard_variable
      integer                           :: externalid = 0       ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      type (type_bulk_variable_id)      :: globalid
      type (type_conserved_quantity_component_list) :: components
   end type type_conserved_quantity_info

#ifdef _FABM_F2003_
   type type_expression
      class (type_expression), pointer :: next        => null()
      character(len=attribute_length)  :: output_name = ''
   end type

   type,extends(type_expression) :: type_bulk_expression
      type (type_bulk_data_pointer),pointer :: out
   end type

   type,extends(type_expression) :: type_horizontal_expression
      type (type_horizontal_data_pointer),pointer :: out
   end type
#endif

   ! ====================================================================================================
   ! Base model type, used by biogeochemical models to inherit from, and by external host to
   ! get variable lists and metadata.
   ! ====================================================================================================

   type type_model_info
      ! Flag determining whether the contents of the type are "frozen", i.e., they will not change anymore.
      logical :: frozen = .false.

      ! Arrays with metadata on model variables.
      type (type_state_variable_info),                _ALLOCATABLE_,dimension(:) :: state_variables                 _NULL_
      type (type_horizontal_state_variable_info),     _ALLOCATABLE_,dimension(:) :: surface_state_variables         _NULL_
      type (type_horizontal_state_variable_info),     _ALLOCATABLE_,dimension(:) :: bottom_state_variables          _NULL_
      type (type_diagnostic_variable_info),           _ALLOCATABLE_,dimension(:) :: diagnostic_variables            _NULL_
      type (type_horizontal_diagnostic_variable_info),_ALLOCATABLE_,dimension(:) :: horizontal_diagnostic_variables _NULL_
      type (type_conserved_quantity_info),            _ALLOCATABLE_,dimension(:) :: conserved_quantities            _NULL_

      ! Pointers for backward compatibility (pre 2013-06-15)
      type (type_horizontal_state_variable_info),     pointer,dimension(:) :: state_variables_ben     => null()
      type (type_horizontal_diagnostic_variable_info),pointer,dimension(:) :: diagnostic_variables_hz => null()

      character(len=attribute_length),_ALLOCATABLE_,dimension(:) :: dependencies        _NULL_
      character(len=attribute_length),_ALLOCATABLE_,dimension(:) :: dependencies_hz     _NULL_
      character(len=attribute_length),_ALLOCATABLE_,dimension(:) :: dependencies_scalar _NULL_

      ! Pointers to linked models in the model tree.
      _CLASS_ (type_model_info),pointer :: parent       => null()
      _CLASS_ (type_model_info),pointer :: first_child  => null()
      _CLASS_ (type_model_info),pointer :: next_sibling => null()

      ! Model name and variable prefixes.
      character(len=64) :: name             = ''
      character(len=64) :: name_prefix      = ''
      character(len=64) :: long_name_prefix = ''

      type (type_bulk_variable_link),      pointer :: first_link            => null()
      type (type_horizontal_variable_link),pointer :: first_horizontal_link => null()
      type (type_scalar_variable_link),    pointer :: first_scalar_link     => null()
      type (type_named_coupling),          pointer :: first_coupling        => null()

#ifdef _FABM_F2003_
      class (type_property), pointer :: first_parameter => null()

      class (type_expression), pointer :: first_expression => null()

   contains

      ! Procedures that can be used to add child models during initialization.
      procedure :: add_child_model_object
      procedure :: add_named_child_model
      generic   :: add_child => add_child_model_object,add_named_child_model

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
      procedure :: register_bulk_state_dependency
      procedure :: register_bottom_state_dependency
      procedure :: register_surface_state_dependency

      procedure :: set_variable_property_real
      procedure :: set_variable_property_integer
      procedure :: set_variable_property_logical
      generic   :: set_variable_property => set_variable_property_real,set_variable_property_integer,set_variable_property_logical

      procedure :: add_conserved_quantity_component

      procedure :: register_bulk_expression_dependency
      procedure :: register_horizontal_expression_dependency

      generic :: register_state_variable      => register_bulk_state_variable,register_bottom_state_variable,register_surface_state_variable
      generic :: register_diagnostic_variable => register_bulk_diagnostic_variable,register_horizontal_diagnostic_variable
      generic :: register_dependency          => register_bulk_dependency, register_bulk_dependency_sn, &
                                                 register_horizontal_dependency, register_horizontal_dependency_sn, &
                                                 register_global_dependency, register_global_dependency_sn, &
                                                 register_bulk_expression_dependency, register_horizontal_expression_dependency
      generic :: register_state_dependency    => register_bulk_state_dependency,register_bottom_state_dependency,register_surface_state_dependency
      generic :: register_conserved_quantity  => register_standard_conserved_quantity, register_custom_conserved_quantity

      ! Procedures that may be used to query parameter values during initialization.
      procedure :: get_real_parameter
      procedure :: get_integer_parameter
      procedure :: get_logical_parameter
      generic :: get_parameter => get_real_parameter,get_integer_parameter,get_logical_parameter

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
      procedure :: check_state              => base_check_state

      ! For backward compatibility only - do not use these in new models!
      procedure :: set_domain               => base_set_domain
      procedure :: do_benthos               => base_do_benthos           ! superceded by do_bottom
      procedure :: do_benthos_ppdd          => base_do_benthos_ppdd      ! superceded by do_bottom_ppdd
      procedure :: get_surface_exchange     => base_get_surface_exchange ! superceded by do_surface
#endif
   end type type_model_info

   ! ====================================================================================================
   ! Derived type for holding global data needed by biogeochemical model tree.
   ! ====================================================================================================

   type type_environment
      ! Declare the arrays for diagnostic variable values.
      real(rk),_ALLOCATABLE_ _ATTR_LOCATION_DIMENSIONS_PLUS_ONE_    :: diag    _NULL_
      real(rk),_ALLOCATABLE_ _ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE_ :: diag_hz _NULL_

#ifdef _FABM_MASK_
      _FABM_MASK_TYPE_,pointer _ATTR_LOCATION_DIMENSIONS_ :: mask => null()
#endif
   end type type_environment

#ifdef _FABM_F2003_
   ! ====================================================================================================
   ! Abstract derived type for a model object factory (generates a model object from a model name)
   ! The only actual implementation of this type is fabm_library.F90.
   ! ====================================================================================================

   abstract interface
      function factory_create_model(modelname,instancename,parent,configunit) result(model)
         import
         character(*),intent(in)           :: modelname,instancename
         integer,     intent(in)           :: configunit
         _CLASS_ (type_model_info),target  :: parent
         _CLASS_ (type_model_info),pointer :: model
      end function
   end interface

   type,abstract,public :: type_abstract_model_factory
      contains
      procedure (factory_create_model),deferred,nopass :: create
   end type

   class (type_abstract_model_factory), pointer,save,public :: factory => null()
#endif

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
      module procedure register_bulk_state_dependency
      module procedure register_bottom_state_dependency
      module procedure register_surface_state_dependency
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

   interface merge_variables
      module procedure merge_bulk_variables
      module procedure merge_horizontal_variables
      module procedure merge_scalar_variables
   end interface

   interface couple_variables
      module procedure couple_bulk_variables
      module procedure couple_horizontal_variables
      module procedure couple_scalar_variables
   end interface

   interface new_link
      module procedure new_bulk_link
      module procedure new_horizontal_link
      module procedure new_scalar_link
   end interface

!-----------------------------------------------------------------------

   contains

#ifdef _FABM_F2003_
   ! Model initialization
   subroutine base_initialize(self,configunit)
      class (type_model_info),intent(inout),target :: self
      integer,                intent(in)           :: configunit
      call fatal_error('base_initialize','derived model '//trim(self%name)//' must implement the "initialize" subroutine.')
   end subroutine base_initialize   
   subroutine base_initialize_state(self,_ARGUMENTS_INITIALIZE_STATE_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_INITIALIZE_STATE_
   end subroutine base_initialize_state
   subroutine base_initialize_horizontal_state(self,_ARGUMENTS_INITIALIZE_HORIZONTAL_STATE_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_INITIALIZE_HORIZONTAL_STATE_
   end subroutine base_initialize_horizontal_state

   ! Providing process rates and diagnostics
   subroutine base_do(self,_ARGUMENTS_DO_)
      class (type_model_info),intent(in) ::  self
      _DECLARE_ARGUMENTS_DO_
   end subroutine base_do
   subroutine base_do_ppdd(self,_ARGUMENTS_DO_PPDD_)
      class (type_model_info),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_PPDD_
   end subroutine base_do_ppdd
   subroutine base_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_model_info),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      call self%do_benthos(_ARGUMENTS_DO_BOTTOM_)
   end subroutine base_do_bottom
   subroutine base_do_bottom_ppdd(self,_ARGUMENTS_DO_BOTTOM_PPDD_)
      class (type_model_info),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_
      call self%do_benthos_ppdd(_ARGUMENTS_DO_BOTTOM_PPDD_)
   end subroutine base_do_bottom_ppdd
   subroutine base_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_
      call self%get_surface_exchange(_ARGUMENTS_DO_SURFACE_)
   end subroutine base_do_surface

   ! Vertical movement, ligth attenuation, feedbacks to drag and albedo
   subroutine base_get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
   end subroutine base_get_vertical_movement
   subroutine base_get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_EXTINCTION_
   end subroutine base_get_light_extinction
   subroutine base_get_drag(self,_ARGUMENTS_GET_DRAG_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_DRAG_
   end subroutine base_get_drag
   subroutine base_get_albedo(self,_ARGUMENTS_GET_ALBEDO_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_ALBEDO_
   end subroutine base_get_albedo

   ! Bookkeeping
   subroutine base_get_conserved_quantities(self,_ARGUMENTS_GET_CONSERVED_QUANTITIES_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_CONSERVED_QUANTITIES_
   end subroutine base_get_conserved_quantities
   subroutine base_check_state(self,_ARGUMENTS_CHECK_STATE_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_CHECK_STATE_
   end subroutine base_check_state

   ! For backward compatibility only
   subroutine base_set_domain(self _ARG_LOCATION_)
      class (type_model_info),intent(inout) :: self
      _DECLARE_LOCATION_ARG_
   end subroutine base_set_domain
   subroutine base_do_benthos(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_model_info),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
   end subroutine base_do_benthos
   subroutine base_do_benthos_ppdd(self,_ARGUMENTS_DO_BOTTOM_PPDD_)
      class (type_model_info),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_
   end subroutine base_do_benthos_ppdd
   subroutine base_get_surface_exchange(self,_ARGUMENTS_DO_SURFACE_)
      class (type_model_info), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_
   end subroutine base_get_surface_exchange

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes model information.
!
! !INTERFACE:
   subroutine add_named_child_model(parent,modelname,configunit,instancename)
!
! !DESCRIPTION:
!  This function initializes the members of a model information derived type,
!  by setting them to a reasonable default value.
!
! !INPUT/OUTPUT PARAMETER:
      _CLASS_ (type_model_info),target,intent(inout) :: parent
      character(len=*),                intent(in)    :: modelname
      integer,                         intent(in)    :: configunit
      character(len=*),optional,       intent(in)    :: instancename
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),pointer :: child
!
!-----------------------------------------------------------------------
!BOC
      if (present(instancename)) then
         child => factory%create(modelname,instancename,parent,configunit)
      else
         child => factory%create(modelname,modelname,parent,configunit)
      end if

   end subroutine add_named_child_model
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes model information.
!
! !INTERFACE:
   subroutine add_child_model_object(parent,model,configunit,instancename)
!
! !DESCRIPTION:
!  This function initializes the members of a model information derived type,
!  by setting them to a reasonable default value.
!
! !INPUT/OUTPUT PARAMETER:
      _CLASS_ (type_model_info),target,intent(inout) :: parent,model
      integer,                         intent(in)    :: configunit
      character(len=*),                intent(in)    :: instancename
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call set_model_name(model,instancename)
      call add_child_model(parent,model)
      call model%initialize(configunit)

   end subroutine add_child_model_object
!EOC

   subroutine set_variable_property_real(self,variable,name,value)
      class (type_model_info),intent(inout) :: self
      class (type_id),        intent(inout) :: variable
      character(len=*),       intent(in)    :: name
      real(rk),               intent(in)    :: value
      if (.not.associated(variable%metadata)) call fatal_error('set_variable_property_real','variable has not been registered')
      call variable%metadata%properties%set_real(name,value)
   end subroutine

   subroutine set_variable_property_integer(self,variable,name,value)
      class (type_model_info),intent(inout) :: self
      class (type_id),        intent(inout) :: variable
      character(len=*),       intent(in)    :: name
      integer,                intent(in)    :: value
      if (.not.associated(variable%metadata)) call fatal_error('set_variable_property_integer','variable has not been registered')
      call variable%metadata%properties%set_integer(name,value)
   end subroutine

   subroutine set_variable_property_logical(self,variable,name,value)
      class (type_model_info),intent(inout) :: self
      class (type_id),        intent(inout) :: variable
      character(len=*),       intent(in)    :: name
      logical,                intent(in)    :: value
      if (.not.associated(variable%metadata)) call fatal_error('set_variable_property_logical','variable has not been registered')
      call variable%metadata%properties%set_logical(name,value)
   end subroutine

   subroutine add_conserved_quantity_component(self,conserved_quantity,state_variable,scale_factor)
      class (type_model_info),           intent(inout) :: self
      class (type_conserved_quantity_id),intent(inout) :: conserved_quantity
      class (type_state_variable_id),    intent(in)    :: state_variable
      real(rk),optional,                 intent(in)    :: scale_factor
      
      type (type_conserved_quantity_component),pointer :: component

      if (.not.associated(state_variable%metadata)) call fatal_error('add_conserved_quantity_component','variable has not been registered')
      if (.not.associated(conserved_quantity%metadata%components%first)) then
         allocate(conserved_quantity%metadata%components%first)
         component => conserved_quantity%metadata%components%first
      else
         component => conserved_quantity%metadata%components%first
         do while (associated(component%next))
            component => component%next
         end do
         allocate(component%next)
         component => component%next
      end if

      if (present(scale_factor)) component%scale_factor = scale_factor
      select type (internal_variable => state_variable%metadata)
         class is (type_bulk_variable)
            call append_data_pointer(internal_variable%alldata,component%state)
      end select
   end subroutine
#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes model information.
!
! !INTERFACE:
   subroutine initialize_model_info(model,name,parent)
!
! !DESCRIPTION:
!  This function initializes the members of a model information derived type,
!  by setting them to a reasonable default value.
!
! !INPUT/OUTPUT PARAMETER:
      _CLASS_ (type_model_info),target,intent(inout)          :: model
      _CLASS_ (type_model_info),target,intent(inout),optional :: parent
      character(len=*),                intent(in)             :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
      call set_model_name(model,name)
      if (present(parent)) call add_child_model(parent,model)

   end subroutine initialize_model_info
!EOC

   subroutine set_model_name(model,name)
      _CLASS_ (type_model_info),target,intent(inout) :: model
      character(len=*),                intent(in)    :: name

      if (model%name/='') call fatal_error('set_model_name','model name has already been set')

      model%name             = name
      model%name_prefix      = trim(name)//'_'
      model%long_name_prefix = trim(name)//' '
   end subroutine set_model_name

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes model information.
!
! !INTERFACE:
   subroutine add_child_model(parent,model)
!
! !DESCRIPTION:
!  This function initializes the members of a model information derived type,
!  by setting them to a reasonable default value.
!
! !INPUT/OUTPUT PARAMETER:
      _CLASS_ (type_model_info),target,intent(inout) :: parent,model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),pointer :: last_child
!
!-----------------------------------------------------------------------
!BOC
      if (associated(model%parent)) call fatal_error('add_child_model','model has already been connected to a parent')

      model%parent => parent
      if (.not.associated(parent%first_child)) then
         parent%first_child => model
      else
         last_child => parent%first_child
         do while (associated(last_child%next_sibling))
            last_child => last_child%next_sibling
         end do
         last_child%next_sibling => model
      end if

   end subroutine add_child_model
!EOC

subroutine new_bulk_link(model,target,name,merge)
   _CLASS_ (type_model_info),       intent(inout) :: model
   type (type_bulk_variable),target,intent(in)    :: target
   character(len=*),                intent(in)    :: name
   logical,optional,                intent(in)    :: merge

   type (type_bulk_variable_link),pointer :: link

   ! First check if a link with this name exists. If so, merge new target with old target.
   link => model%first_link
   do while (associated(link))
      if (link%name==name) then
         if (.not.present(merge)) call fatal_error('new_bulk_link','Link '//trim(name)//' already exists.')
         if (merge) call merge_variables(link%target,target)
         return
      end if
      link => link%next
   end do

   ! Append a new link to the list.
   if (.not.associated(model%first_link)) then
      allocate(model%first_link)
      link => model%first_link
   else
      link => model%first_link
      do while (associated(link%next))
         link => link%next
      end do
      allocate(link%next)
      link => link%next
   end if
   
   ! Set link attributes.
   link%name = name
   link%target => target
end subroutine new_bulk_link

subroutine new_horizontal_link(first,target,name,merge)
   type (type_horizontal_variable_link),pointer :: first
   type (type_horizontal_variable),target,intent(in) :: target
   character(len=*),                      intent(in) :: name
   logical,optional,                      intent(in) :: merge

   type (type_horizontal_variable_link),pointer :: link

   ! First check if a link with this name exists. If so, merge new target with old target.
   link => first
   do while (associated(link))
      if (link%name==name) then
         if (.not.present(merge)) call fatal_error('new_horizontal_link','Link '//trim(name)//' already exists.')
         if (merge) call merge_variables(link%target,target)
         return
      end if
      link => link%next
   end do

   ! Append a new link to the list.
   if (.not.associated(first)) then
      allocate(first)
      link => first
   else
      link => first
      do while (associated(link%next))
         link => link%next
      end do
      allocate(link%next)
      link => link%next
   end if
   
   ! Set link attributes.
   link%name = name
   link%target => target
end subroutine new_horizontal_link

subroutine new_scalar_link(first,target,name,merge)
   type (type_scalar_variable_link),pointer :: first
   type (type_scalar_variable),target,intent(in) :: target
   character(len=*),                  intent(in) :: name
   logical,optional,                  intent(in) :: merge

   type (type_scalar_variable_link),pointer :: link

   ! First check if a link with this name exists. If so, merge new target with old target.
   link => first
   do while (associated(link))
      if (link%name==name) then
         if (.not.present(merge)) call fatal_error('new_scalar_link','Link '//trim(name)//' already exists.')
         if (merge) call merge_variables(link%target,target)
         return
      end if
      link => link%next
   end do

   ! Append a new link to the list.
   if (.not.associated(first)) then
      allocate(first)
      link => first
   else
      link => first
      do while (associated(link%next))
         link => link%next
      end do
      allocate(link%next)
      link => link%next
   end if
   
   ! Set link attributes.
   link%name = name
   link%target => target
end subroutine new_scalar_link

subroutine new_coupling(model,master,slave)
   _CLASS_ (type_model_info),intent(inout) :: model
   character(len=*),         intent(in)    :: master,slave

   type (type_named_coupling),pointer :: link

   ! Append a new coupling link to the list.
   if (.not.associated(model%first_coupling)) then
      allocate(model%first_coupling)
      link => model%first_coupling
   else
      link => model%first_coupling
      do while (associated(link%next))
         link => link%next
      end do
      allocate(link%next)
      link => link%next
   end if
   
   ! Set coupling attributes.
   link%master = master
   link%slave = slave
end subroutine new_coupling

subroutine append_index(array,index)
   type (type_integer_pointer),dimension(:),_ALLOCATABLE_ :: array
   integer,target :: index
   type (type_integer_pointer),allocatable :: oldarray(:)

   ! Create a new list of integer pointers, or extend it if already allocated.
   if (.not._ALLOCATED_(array)) then
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
   type (type_bulk_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: array
   type (type_bulk_data_pointer),target :: data
   type (type_bulk_data_pointer_pointer),allocatable :: oldarray(:)
   integer :: i

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not._ALLOCATED_(array)) then
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
   type (type_horizontal_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: array
   type (type_horizontal_data_pointer),target :: data
   type (type_horizontal_data_pointer_pointer),allocatable :: oldarray(:)
   integer :: i

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not._ALLOCATED_(array)) then
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
   type (type_scalar_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: array
   type (type_scalar_data_pointer),target :: data
   type (type_scalar_data_pointer_pointer),allocatable :: oldarray(:)
   integer :: i

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not._ALLOCATED_(array)) then
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

subroutine append_string(array,string,exists)
   character(len=attribute_length),dimension(:),_ALLOCATABLE_ :: array
   character(len=*),intent(in) :: string
   logical,intent(out),optional :: exists
   integer :: i
   character(len=attribute_length),allocatable :: oldarray(:)

   if (.not._ALLOCATED_(array)) then
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
   subroutine freeze_model_info(model)
!
! !DESCRIPTION:
!  This function finalizes model initialization. It will resolve all remaining
!  internal dependencies (coupling commands) and generate final authorative lists
!  of state variables, diagnostic variables, conserved quantities and readable
!  variables ("dependencies").
!
! !INPUT/OUTPUT PARAMETER:
      _CLASS_ (type_model_info),intent(inout) :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      if (associated(model%parent)) call fatal_error('freeze_model_info', &
         'freeze_model_info can only operate on the root model.')

      call couple_standard_variables(model)
      call process_coupling_tasks(model)
      call classify_variables(model)

   end subroutine freeze_model_info
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new state variable
!
! !INTERFACE:
   recursive subroutine register_bulk_state_variable(model, id, name, units, long_name, &
                                                     initial_value, vertical_movement, specific_light_extinction, &
                                                     minimum, maximum, missing_value, &
                                                     no_precipitation_dilution, no_river_dilution, &
                                                     standard_variable, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical state variable in the global model database.
!  It returns an identifier that may be used later to retrieve the value of the state variable.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),     intent(inout)                 :: model
      type (type_state_variable_id), intent(inout),target          :: id
      type (type_bulk_variable),     intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                   intent(in)          :: name, long_name, units
      real(rk),                           intent(in),optional :: initial_value,vertical_movement,specific_light_extinction
      real(rk),                           intent(in),optional :: minimum, maximum,missing_value
      logical,                            intent(in),optional :: no_precipitation_dilution,no_river_dilution
      type (type_bulk_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable),pointer :: curinfo
      character(len=256)                :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_bulk_state_variable', &
         'State variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the variable id.
      if (associated(model%parent)) then
         call register_bulk_state_variable(model%parent,id,trim(model%name_prefix)//trim(name), &
            units,trim(model%long_name_prefix)//' '//trim(long_name),       &
            initial_value             = initial_value,             &
            vertical_movement         = vertical_movement,         &
            specific_light_extinction = specific_light_extinction, &
            minimum                   = minimum,                   &
            maximum                   = maximum,                   &
            missing_value             = missing_value,             &
            no_precipitation_dilution = no_precipitation_dilution, &
            no_river_dilution         = no_river_dilution,         &
            standard_variable         = standard_variable,         &
            target                    = curinfo)
      else
         ! Store customized information on state variable.
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         if (present(initial_value))             curinfo%initial_value             = initial_value
         if (present(minimum))                   curinfo%minimum                   = minimum
         if (present(maximum))                   curinfo%maximum                   = maximum
         if (present(missing_value))             curinfo%missing_value             = missing_value
         if (present(vertical_movement))         curinfo%vertical_movement         = vertical_movement
         if (present(specific_light_extinction)) curinfo%specific_light_extinction = specific_light_extinction
         if (present(no_precipitation_dilution)) curinfo%no_precipitation_dilution = no_precipitation_dilution
         if (present(no_river_dilution))         curinfo%no_river_dilution         = no_river_dilution
         if (present(standard_variable))         curinfo%standard_variable         = standard_variable
         call append_index(curinfo%state_indices,id%state_index)
         call append_data_pointer(curinfo%alldata,id%data)
         if (id%name/='') call fatal_error('fabm_types::register_bulk_state_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
#ifdef _FABM_F2003_
         id%metadata => curinfo
#endif

         ! Ensure that initial value falls within prescribed valid range.
         if (curinfo%initial_value<curinfo%minimum .or. curinfo%initial_value>curinfo%maximum) then
            write (text,*) 'Initial value',curinfo%initial_value,'for variable "'//trim(name)//'" lies&
                  &outside allowed range',curinfo%minimum,'to',curinfo%maximum
            call fatal_error('fabm_types::register_bulk_state_variable',text)
         end if
      end if

      call new_link(model,curinfo,name,.not.present(target))

   end subroutine register_bulk_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new state variable
!
! !INTERFACE:
   recursive subroutine register_bottom_state_variable(model, id, name, units, long_name, &
                                                       initial_value, minimum, maximum, missing_value, &
                                                       standard_variable, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical state variable in the global model database.
!  It returns an identifier that may be used later to retrieve the value of the state variable.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),           intent(inout)                 :: model
      type (type_bottom_state_variable_id),intent(inout),target          :: id
      type (type_horizontal_variable),     intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                         intent(in)          :: name, long_name, units
      real(rk),                                 intent(in),optional :: initial_value
      real(rk),                                 intent(in),optional :: minimum, maximum,missing_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable),pointer :: curinfo
      character(len=256)                      :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_bottom_state_variable', &
         'State variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the variable id.
      if (associated(model%parent)) then
         call register_bottom_state_variable(model%parent,id,trim(model%name_prefix)//trim(name),      &
                                             units,trim(model%long_name_prefix)//' '//trim(long_name), &
                                             initial_value             = initial_value,                &
                                             minimum                   = minimum,                      &
                                             maximum                   = maximum,                      &
                                             missing_value             = missing_value,                &
                                             standard_variable         = standard_variable,            &
                                             target                    = curinfo)
      else
         ! Store customized information on state variable.
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         curinfo%domain = domain_bottom
         if (present(initial_value))     curinfo%initial_value     = initial_value
         if (present(minimum))           curinfo%minimum           = minimum
         if (present(maximum))           curinfo%maximum           = maximum
         if (present(missing_value))     curinfo%missing_value     = missing_value
         if (present(standard_variable)) curinfo%standard_variable = standard_variable
         call append_index(curinfo%state_indices,id%bottom_state_index)
         call append_data_pointer(curinfo%alldata,id%horizontal_data)
         if (id%name/='') call fatal_error('fabm_types::register_bottom_state_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
#ifdef _FABM_F2003_
         id%metadata => curinfo
#endif

         ! Ensure that initial value falls within prescribed valid range.
         if (curinfo%initial_value<curinfo%minimum .or. curinfo%initial_value>curinfo%maximum) then
            write (text,*) 'Initial value',curinfo%initial_value,'for variable "'//trim(name)//'" lies&
                  &outside allowed range',curinfo%minimum,'to',curinfo%maximum
            call fatal_error('fabm_types::register_bottom_state_variable',text)
         end if
      end if

      call new_link(model%first_horizontal_link,curinfo,name,.not.present(target))

   end subroutine register_bottom_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new state variable
!
! !INTERFACE:
   recursive subroutine register_surface_state_variable(model, id, name, units, long_name, &
                                                       initial_value, minimum, maximum, missing_value, &
                                                       standard_variable, target)
!
! !DESCRIPTION:
!  This subroutine registers a new surface-bound biogeochemical state variable in the global model database.
!  The identifier "id" may be used later to retrieve the value of the state variable.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),           intent(inout)                 :: model
      type (type_surface_state_variable_id),intent(inout),target         :: id
      type (type_horizontal_variable),     intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                         intent(in)          :: name, long_name, units
      real(rk),                                 intent(in),optional :: initial_value
      real(rk),                                 intent(in),optional :: minimum, maximum,missing_value
      type (type_horizontal_standard_variable), intent(in),optional :: standard_variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable),pointer :: curinfo
      character(len=256)                      :: text
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_surface_state_variable', &
         'State variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the variable id.
      if (associated(model%parent)) then
         call register_surface_state_variable(model%parent,id,trim(model%name_prefix)//trim(name),      &
                                             units,trim(model%long_name_prefix)//' '//trim(long_name), &
                                             initial_value             = initial_value,                &
                                             minimum                   = minimum,                      &
                                             maximum                   = maximum,                      &
                                             missing_value             = missing_value,                &
                                             standard_variable         = standard_variable,            &
                                             target                    = curinfo)
      else
         ! Store customized information on state variable.
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         curinfo%domain = domain_surface
         if (present(initial_value))     curinfo%initial_value     = initial_value
         if (present(minimum))           curinfo%minimum           = minimum
         if (present(maximum))           curinfo%maximum           = maximum
         if (present(missing_value))     curinfo%missing_value     = missing_value
         if (present(standard_variable)) curinfo%standard_variable = standard_variable
         call append_index(curinfo%state_indices,id%surface_state_index)
         call append_data_pointer(curinfo%alldata,id%horizontal_data)
         if (id%name/='') call fatal_error('fabm_types::register_surface_state_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
#ifdef _FABM_F2003_
         id%metadata => curinfo
#endif

         ! Ensure that initial value falls within prescribed valid range.
         if (curinfo%initial_value<curinfo%minimum .or. curinfo%initial_value>curinfo%maximum) then
            write (text,*) 'Initial value',curinfo%initial_value,'for variable "'//trim(name)//'" lies&
                  &outside allowed range',curinfo%minimum,'to',curinfo%maximum
            call fatal_error('fabm_types::register_bottom_state_variable',text)
         end if
      end if

      call new_link(model%first_horizontal_link,curinfo,name,.not.present(target))

   end subroutine register_surface_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   recursive subroutine register_bulk_diagnostic_variable(model, id, name, units, long_name, &
                                                          time_treatment, missing_value, standard_variable, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical diagnostic variable in the global model database.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),         intent(inout)                 :: model
      type (type_diagnostic_variable_id),intent(inout),target          :: id
      type (type_bulk_variable),         intent(inout),target,optional :: target
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
      type (type_bulk_variable),pointer :: curinfo
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_bulk_diagnostic_variable',&
                                             'Diagnostic variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the diagnostic variable id.
      if (associated(model%parent)) then
         call register_bulk_diagnostic_variable(model%parent,id,trim(model%name_prefix)//trim(name), &
                                                units,                                               &
                                                trim(model%long_name_prefix)//' '//trim(long_name),  &
                                                time_treatment=time_treatment,                       &
                                                missing_value=missing_value,                         &
                                                standard_variable=standard_variable,                 &
                                                target = curinfo)
      else
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         if (present(time_treatment))    curinfo%time_treatment    = time_treatment
         if (present(missing_value))     curinfo%missing_value     = missing_value
         if (present(standard_variable)) curinfo%standard_variable = standard_variable
         call append_index(curinfo%write_indices,id%diag_index)
         if (id%name/='') call fatal_error('fabm_types::register_bulk_diagnostic_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
#ifdef _FABM_F2003_
         id%metadata => curinfo
#endif
      end if

      call new_link(model,curinfo,name,.not.present(target))

   end subroutine register_bulk_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   recursive subroutine register_horizontal_diagnostic_variable(model, id, name, units, long_name, &
                                                                time_treatment, missing_value, standard_variable, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical diagnostic variable in the global model database.
!
! !INPUT/OUTPUT PARAMETER:
      _CLASS_ (type_model_info),                    intent(inout)                 :: model
      type (type_horizontal_diagnostic_variable_id),intent(inout),target          :: id
      type (type_horizontal_variable),              intent(inout),target,optional :: target
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
      type (type_horizontal_variable),pointer :: curinfo
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_horizontal_diagnostic_variable',&
                                             'Diagnostic variables may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the diagnostic variable id.
      if (associated(model%parent)) then
         call register_horizontal_diagnostic_variable(model%parent,id,trim(model%name_prefix)//trim(name), &
                                                      units,                                               &
                                                      trim(model%long_name_prefix)//' '//trim(long_name),  &
                                                      time_treatment=time_treatment,                       &
                                                      missing_value=missing_value,                         &
                                                      standard_variable=standard_variable,                 &
                                                      target = curinfo)
      else
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         if (present(time_treatment))    curinfo%time_treatment    = time_treatment
         if (present(missing_value))     curinfo%missing_value     = missing_value
         if (present(standard_variable)) curinfo%standard_variable = standard_variable
         call append_index(curinfo%write_indices,id%horizontal_diag_index)
         if (id%name/='') call fatal_error('fabm_types::register_horizontal_diagnostic_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
#ifdef _FABM_F2003_
         id%metadata => curinfo
#endif
      end if

      call new_link(model%first_horizontal_link,curinfo,name,.not.present(target))

   end subroutine register_horizontal_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a conserved quantity taken from a standard set
!
! !INTERFACE:
   recursive subroutine register_standard_conserved_quantity(model, id, standard_variable, name, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemically conserved quantity in the global
!  model database.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),        intent(inout)                 :: model
      type (type_conserved_quantity_id),intent(inout),target          :: id
      type (type_bulk_variable),        intent(inout),target,optional :: target
      character(len=*),                 intent(in), optional          :: name
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
      type (type_bulk_variable),pointer :: curinfo
      character(len=attribute_length) :: name_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_conserved_quantity',&
         'Conserved quantities may only be registered during initialization.')

      if (present(name)) then
         name_eff = name
      else
         name_eff = standard_variable%name
      end if

      ! Either use the provided variable object, or create a new one.
      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the conserved quantity id.
      if (associated(model%parent)) then
         call register_standard_conserved_quantity(model%parent,id,standard_variable, &
                                                   name=trim(model%name_prefix)//name_eff,target=curinfo)
      else
         curinfo%standard_variable = standard_variable
         curinfo%name      = standard_variable%name
         curinfo%units     = standard_variable%units
         curinfo%long_name = standard_variable%name
         call append_index(curinfo%cons_indices,id%cons_index)
         if (id%name/='') call fatal_error('fabm_types::register_conserved_quantity', &
            'Identifier supplied for '//trim(name_eff)//' is already used by '//trim(id%name)//'.')
         id%name = name_eff
#ifdef _FABM_F2003_
         id%metadata => curinfo
#endif
      end if

      call new_link(model,curinfo,name_eff,.not.present(target))

   end subroutine register_standard_conserved_quantity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new conserved quantity
!
! !INTERFACE:
   recursive subroutine register_custom_conserved_quantity(model, id, name, units, long_name, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemically conserved quantity in the global
!  model database.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),        intent(inout)                 :: model
      type (type_conserved_quantity_id),intent(inout),target          :: id
      type (type_bulk_variable),        intent(inout),target,optional :: target
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
      type (type_bulk_variable),pointer :: curinfo
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_conserved_quantity',&
         'Conserved quantities may only be registered during initialization.')

      ! Either use the provided variable object, or create a new one.
      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the conserved quantity id.
      if (associated(model%parent)) then
         call register_conserved_quantity(model%parent,id,trim(model%name_prefix)//name,      &
                                          units,trim(model%long_name_prefix)//' '//long_name,target=curinfo)
      else
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         call append_index(curinfo%cons_indices,id%cons_index)
         if (id%name/='') call fatal_error('fabm_types::register_conserved_quantity', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
#ifdef _FABM_F2003_
         id%metadata => curinfo
#endif
      end if

      call new_link(model,curinfo,curinfo%name,.not.present(target))

   end subroutine register_custom_conserved_quantity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external state variable
!
! !INTERFACE:
   recursive subroutine register_bulk_state_dependency(model,id,name)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),    intent(inout)        :: model
      type (type_state_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),intent(in) :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      call register_bulk_state_variable(model, id, name, '', name)
      call new_coupling(model,name,name)

   end subroutine register_bulk_state_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external bottom state variable
!
! !INTERFACE:
   recursive subroutine register_bottom_state_dependency(model,id,name)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),           intent(inout)        :: model
      type (type_bottom_state_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),intent(in) :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      call register_bottom_state_variable(model, id, name, '', name)
      call new_coupling(model,name,name)

   end subroutine register_bottom_state_dependency
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on an external surface-bound state variable
!
! !INTERFACE:
   recursive subroutine register_surface_state_dependency(model,id,name)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),            intent(inout)        :: model
      type (type_surface_state_variable_id),intent(inout),target :: id
!
! !INPUT PARAMETERS:
      character(len=*),intent(in) :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      call register_surface_state_variable(model, id, name, '', name)
      call new_coupling(model,name,name)

   end subroutine register_surface_state_dependency
!EOC

   subroutine register_bulk_dependency_sn(model,id,standard_variable)
      _CLASS_ (type_model_info),         intent(inout)        :: model
      type (type_dependency_id),         intent(inout),target :: id
      type (type_bulk_standard_variable),intent(in)           :: standard_variable
      
      call register_bulk_dependency(model,id,standard_variable%name,standard_variable)
   end subroutine register_bulk_dependency_sn

   subroutine register_horizontal_dependency_sn(model,id,standard_variable)
      _CLASS_ (type_model_info),               intent(inout)        :: model
      type (type_horizontal_dependency_id),    intent(inout),target :: id
      type (type_horizontal_standard_variable),intent(in)           :: standard_variable
      
      call register_horizontal_dependency(model,id,standard_variable%name,standard_variable)
   end subroutine register_horizontal_dependency_sn

   subroutine register_global_dependency_sn(model,id,standard_variable)
      _CLASS_ (type_model_info),           intent(inout)        :: model
      type (type_global_dependency_id),    intent(inout),target :: id
      type (type_global_standard_variable),intent(in)           :: standard_variable
      
      call register_global_dependency(model,id,standard_variable%name,standard_variable)
   end subroutine register_global_dependency_sn

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a variable defined on
! the full model domain.
!
! !INTERFACE:
   recursive subroutine register_bulk_dependency(model,id,name,standard_variable,target)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),intent(inout)                 :: model
      type (type_dependency_id),intent(inout),target          :: id
      type (type_bulk_variable),intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                  intent(in)           :: name
      type (type_bulk_standard_variable),intent(in),optional  :: standard_variable
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable),pointer :: curinfo
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_bulk_dependency',&
         'Dependencies may only be registered during initialization.')

      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      if (associated(model%parent)) then
         call register_bulk_dependency(model%parent,id,trim(model%name_prefix)//trim(name),standard_variable,curinfo)
      else
         curinfo%name = name
         if (present(standard_variable)) curinfo%standard_variable = standard_variable
         call append_data_pointer(curinfo%alldata,id%data)
         if (id%name/='') call fatal_error('fabm_types::register_bulk_dependency', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
      end if

      if (associated(model%parent).and.(.not.present(target)).and.(.not.present(standard_variable))) &
         call new_coupling(model,name,name)

      call new_link(model,curinfo,name,.not.present(target))

   end subroutine register_bulk_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a variable defined on
! a horizontal slice of the model domain.
!
! !INTERFACE:
   recursive subroutine register_horizontal_dependency(model,id,name,standard_variable,target)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),           intent(inout)                 :: model
      type (type_horizontal_dependency_id),intent(inout),target          :: id
      type (type_horizontal_variable),     intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                        intent(in)           :: name
      type (type_horizontal_standard_variable),intent(in),optional  :: standard_variable
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_horizontal_variable),pointer :: curinfo
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_horizontal_dependency',&
         'Dependencies may only be registered during initialization.')

      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      if (associated(model%parent)) then
         call register_horizontal_dependency(model%parent,id,trim(model%name_prefix)//trim(name),standard_variable,curinfo)
      else
         curinfo%name = name
         if (present(standard_variable)) curinfo%standard_variable = standard_variable
         call append_data_pointer(curinfo%alldata,id%horizontal_data)
         if (id%name/='') call fatal_error('fabm_types::register_horizontal_dependency', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
      end if

      if (associated(model%parent).and.(.not.present(target)).and.(.not.present(standard_variable))) &
         call new_coupling(model,name,name)

      call new_link(model%first_horizontal_link,curinfo,name,.not.present(target))

   end subroutine register_horizontal_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on a global (space-
! independent) variable.
!
! !INTERFACE:
   recursive subroutine register_global_dependency(model,id,name,standard_variable,target)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),       intent(inout)                 :: model
      type (type_global_dependency_id),intent(inout),target          :: id
      type (type_scalar_variable),     intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                    intent(in)          :: name
      type (type_global_standard_variable),intent(in),optional :: standard_variable
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_scalar_variable),pointer :: curinfo
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (model%frozen) call fatal_error('fabm_types::register_global_dependency',&
         'Dependencies may only be registered during initialization.')

      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if

      if (associated(model%parent)) then
         call register_global_dependency(model%parent,id,trim(model%name_prefix)//trim(name),standard_variable,curinfo)
      else
         curinfo%name = name
         if (present(standard_variable)) curinfo%standard_variable = standard_variable
         call append_data_pointer(curinfo%alldata,id%global_data)
         if (id%name/='') call fatal_error('fabm_types::register_global_dependency', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
      end if

      if (associated(model%parent).and.(.not.present(target)).and.(.not.present(standard_variable))) &
         call new_coupling(model,name,name)

      call new_link(model%first_scalar_link,curinfo,name,.not.present(target))

   end subroutine register_global_dependency
!EOC

#ifdef _FABM_F2003_

subroutine register_bulk_expression_dependency(self,id,expression)
   class (type_model_info),       intent(inout)        :: self
   type (type_dependency_id),     intent(inout),target :: id
   class (type_bulk_expression),  intent(in)           :: expression

   class (type_bulk_expression),allocatable :: copy

   allocate(copy,source=expression)
   copy%out => id%data
   call self%register_dependency(id,trim(copy%output_name))
   copy%output_name = id%name

   call register_expression(self,copy)
   deallocate(copy)
end subroutine

subroutine register_horizontal_expression_dependency(self,id,expression)
   class (type_model_info),             intent(inout)        :: self
   type (type_horizontal_dependency_id),intent(inout),target :: id
   class (type_horizontal_expression),  intent(in)           :: expression

   class (type_horizontal_expression),allocatable :: copy

   allocate(copy,source=expression)
   copy%out => id%horizontal_data
   call self%register_dependency(id,trim(copy%output_name))
   copy%output_name = id%name

   call register_expression(self,copy)
   deallocate(copy)
end subroutine

recursive subroutine register_expression(self,expression)
   class (type_model_info), intent(inout)           :: self
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

recursive subroutine get_real_parameter(self,value,name,units,long_name,scale_factor,default,path)
! !INPUT PARAMETERS:
   class (type_model_info), intent(inout), target :: self
   real(rk),        intent(inout)       :: value
   character(len=*),intent(in)          :: name
   character(len=*),intent(in),optional :: units,long_name,path
   real(rk),        intent(in),optional :: scale_factor,default
!
!EOP
!
! !LOCAL VARIABLES:
   type (type_real_property),pointer :: current_parameter
!
!-----------------------------------------------------------------------
!BOC
   if (present(default)) value = default
   
   if (associated(self%parent)) then
      if (present(path)) then
         call self%parent%get_parameter(value,name,units,long_name,1.0_rk,default,path=trim(self%name)//'/'//trim(path))
      else
         call self%parent%get_parameter(value,name,units,long_name,1.0_rk,default,path=trim(self%name))
      end if
   end if

   ! Store parameter settings
   allocate(current_parameter)
   current_parameter%value = value

   call add_parameter(self,current_parameter,name,units,long_name)

   if (present(scale_factor)) value = value*scale_factor

end subroutine get_real_parameter
!EOC

recursive subroutine get_integer_parameter(self,value,name,units,long_name,default,path)
! !INPUT PARAMETERS:
   class (type_model_info), intent(inout), target :: self
   integer,         intent(inout)       :: value
   character(len=*),intent(in)          :: name
   character(len=*),intent(in),optional :: units,long_name,path
   integer,         intent(in),optional :: default
!
!EOP
!
! !LOCAL VARIABLES:
   type (type_integer_property),pointer :: current_parameter
!
!-----------------------------------------------------------------------
!BOC
   if (present(default)) value = default
   
   if (associated(self%parent)) then
      if (present(path)) then
         call self%parent%get_parameter(value,name,units,long_name,default,path=trim(self%name)//'/'//trim(path))
      else
         call self%parent%get_parameter(value,name,units,long_name,default,path=trim(self%name))
      end if
   end if

   ! Store parameter settings
   allocate(current_parameter)
   current_parameter%value = value

   call add_parameter(self,current_parameter,name,units,long_name)

end subroutine get_integer_parameter
!EOC

recursive subroutine get_logical_parameter(self,value,name,units,long_name,default,path)
! !INPUT PARAMETERS:
   class (type_model_info), intent(inout), target :: self
   logical,         intent(inout)       :: value
   character(len=*),intent(in)          :: name
   character(len=*),intent(in),optional :: units,long_name,path
   logical,         intent(in),optional :: default
!
!EOP
!
! !LOCAL VARIABLES:
   type (type_logical_property),pointer :: current_parameter
!
!-----------------------------------------------------------------------
!BOC
   if (present(default)) value = default
   
   if (associated(self%parent)) then
      if (present(path)) then
         call self%parent%get_parameter(value,name,units,long_name,default,path=trim(self%name)//'/'//trim(path))
      else
         call self%parent%get_parameter(value,name,units,long_name,default,path=trim(self%name))
      end if
   end if

   ! Store parameter settings
   allocate(current_parameter)
   current_parameter%value = value

   call add_parameter(self,current_parameter,name,units,long_name)

end subroutine get_logical_parameter
!EOC

subroutine add_parameter(model,parameter,name,units,long_name)
! !INPUT PARAMETERS:
   class (type_model_info), intent(inout), target :: model
   class (type_property),   intent(inout), target :: parameter
   character(len=*),        intent(in)            :: name
   character(len=*),        intent(in),optional   :: units,long_name
!
!EOP
!
! !LOCAL VARIABLES:
   class (type_property),pointer :: current_parameter
!
!-----------------------------------------------------------------------
!BOC
   ! Store metadata
   parameter%name = name
   if (present(units))     parameter%units     = units
   if (present(long_name)) parameter%long_name = long_name

   ! Append new parameter to model list.
   if (.not.associated(model%first_parameter)) then
      model%first_parameter => parameter
   else
      current_parameter => model%first_parameter
      do while (associated(current_parameter%next))
         current_parameter => current_parameter%next
      end do
      current_parameter%next => parameter
   end if

end subroutine add_parameter
!EOC
#endif

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
      _CLASS_ (type_model_info),intent(inout),target :: model
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_bulk_variable_link),      pointer :: bulk_link,bulk_link2
      type (type_horizontal_variable_link),pointer :: horizontal_link,horizontal_link2
      type (type_scalar_variable_link),    pointer :: scalar_link,scalar_link2
      character(len=attribute_length),_ALLOCATABLE_ :: processed(:)
      logical                                       :: exists
!
!-----------------------------------------------------------------------
!BOC
      bulk_link => model%first_link
      do while (associated(bulk_link))
         if (.not.is_null_standard_variable(bulk_link%target%standard_variable)) then
            call append_string(processed,bulk_link%target%standard_variable%name,exists=exists)
            if (exists) exit
            bulk_link2 => bulk_link%next
            do while (associated(bulk_link2))
               if (compare_standard_variables(bulk_link%target%standard_variable,bulk_link2%target%standard_variable) .and. .not. &
                   is_null_standard_variable(bulk_link2%target%standard_variable)) &
                  call couple_variables(model,bulk_link%target,bulk_link2%target)
               bulk_link2 => bulk_link2%next
            end do
         end if
         bulk_link => bulk_link%next
      end do

      if (_ALLOCATED_(processed)) deallocate(processed)

      horizontal_link => model%first_horizontal_link
      do while (associated(horizontal_link))
         if (.not.is_null_standard_variable(horizontal_link%target%standard_variable)) then
            call append_string(processed,horizontal_link%target%standard_variable%name,exists=exists)
            if (exists) exit
            horizontal_link2 => horizontal_link%next
            do while (associated(horizontal_link2))
               if (compare_standard_variables(horizontal_link%target%standard_variable,horizontal_link2%target%standard_variable).and. .not. &
                   is_null_standard_variable(horizontal_link2%target%standard_variable)) &
                  call couple_variables(model,horizontal_link%target,horizontal_link2%target)
               horizontal_link2 => horizontal_link2%next
            end do
         end if
         horizontal_link => horizontal_link%next
      end do

      if (_ALLOCATED_(processed)) deallocate(processed)

      scalar_link => model%first_scalar_link
      do while (associated(scalar_link))
         if (.not.is_null_standard_variable(scalar_link%target%standard_variable)) then
            call append_string(processed,scalar_link%target%standard_variable%name,exists=exists)
            if (exists) exit
            scalar_link2 => scalar_link%next
            do while (associated(scalar_link2))
               if (compare_standard_variables(scalar_link%target%standard_variable,scalar_link2%target%standard_variable).and. .not. &
                   is_null_standard_variable(scalar_link2%target%standard_variable)) &
                  call couple_variables(model,scalar_link%target,scalar_link2%target)
               scalar_link2 => scalar_link2%next
            end do
         end if
         scalar_link => scalar_link%next
      end do

   end subroutine couple_standard_variables
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process all model-specific coupling tasks.
!
! !INTERFACE:
   recursive subroutine process_coupling_tasks(model)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      _CLASS_ (type_model_info),intent(inout),target   :: model
!
!EOP
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),           pointer :: root,child
      type (type_named_coupling),          pointer :: coupling
      type (type_bulk_variable_link),      pointer :: bulk_link
      type (type_bulk_variable),           pointer :: bulk_master
      type (type_horizontal_variable_link),pointer :: horizontal_link
      type (type_horizontal_variable)     ,pointer :: horizontal_master
      type (type_scalar_variable_link),    pointer :: scalar_link
      type (type_scalar_variable),         pointer :: scalar_master
!
!-----------------------------------------------------------------------
!BOC
      ! Find root model, which will handle the individual coupling tasks.
      root => model
      do while (associated(root%parent))
         root => root%parent
      end do
      
      ! Enumerate named couplings, locate slave and target variables, and couple them.
      coupling => model%first_coupling
      do while (associated(coupling))

         ! Try to find slave variable amongst bulk variables.
         bulk_link => model%first_link
         do while (associated(bulk_link))
            if (bulk_link%name==coupling%slave) then
               bulk_master => find_bulk_variable(model,coupling%master,bulk_link)
               if (associated(bulk_master)) call couple_variables(root,bulk_master,bulk_link%target)
            end if
            bulk_link => bulk_link%next
         end do

         ! Try to find slave variable amongst horizontal variables.
         horizontal_link => model%first_horizontal_link
         do while (associated(horizontal_link))
            if (horizontal_link%name==coupling%slave) then
               horizontal_master => find_horizontal_variable(model,coupling%master,horizontal_link)
               if (associated(horizontal_master)) call couple_variables(root,horizontal_master,horizontal_link%target)
            end if
            horizontal_link => horizontal_link%next
         end do

         ! Try to find slave variable amongst scalar variables.
         scalar_link => model%first_scalar_link
         do while (associated(scalar_link))
            if (scalar_link%name==coupling%slave) then
               scalar_master => find_scalar_variable(model,coupling%master,scalar_link)
               if (associated(scalar_master)) call couple_variables(root,scalar_master,scalar_link%target)
            end if
            scalar_link => scalar_link%next
         end do

         coupling => coupling%next
      end do

      ! Process coupling tasks registered with child models.
      child => model%first_child
      do while (associated(child))
         call process_coupling_tasks(child)
         child => child%next_sibling
      end do
      
   end subroutine process_coupling_tasks
!EOC

recursive subroutine couple_bulk_variables(model,master,slave)
   _CLASS_ (type_model_info),intent(inout),target :: model
   type (type_bulk_variable),intent(inout),target :: master
   type (type_bulk_variable),intent(in),   target :: slave

   type (type_bulk_variable),     pointer :: pslave
   type (type_bulk_variable_link),pointer :: link
   _CLASS_ (type_model_info),     pointer :: child
      
   pslave => slave

   if (associated(pslave,master)) return
   
   ! Process all links and if they used to refer to the specified slave,
   ! redirect them to the specified master.
   link => model%first_link
   do while (associated(link))
      if (associated(pslave,link%target)) then
         link%target => master
         link%coupled = .true.
      end if
      link => link%next
   end do
   
   ! Allow child models to do the same.
   child => model%first_child
   do while (associated(child))
      call couple_variables(child,master,pslave)
      child => child%next_sibling
   end do
   
   if (.not.associated(model%parent)) then
      call merge_variables(master,pslave)
      deallocate(pslave)
   end if
end subroutine couple_bulk_variables

recursive subroutine couple_horizontal_variables(model,master,slave)
   _CLASS_ (type_model_info),      intent(inout),target :: model
   type (type_horizontal_variable),intent(inout),target :: master
   type (type_horizontal_variable),intent(in ),  target :: slave

   type (type_horizontal_variable),     pointer :: pslave
   type (type_horizontal_variable_link),pointer :: link
   _CLASS_ (type_model_info),           pointer :: child

   pslave => slave

   if (associated(pslave,master)) return
   
   ! Process all links and if they used to refer to the specified slave,
   ! redirect them to the specified master.
   link => model%first_horizontal_link
   do while (associated(link))
      if (associated(link%target,pslave)) then
         link%target => master
         link%coupled = .true.
      end if
      link => link%next
   end do
   
   ! Allow child models to do the same.
   child => model%first_child
   do while (associated(child))
      call couple_variables(child,master,pslave)
      child => child%next_sibling
   end do
   
   if (.not.associated(model%parent)) then
      call merge_variables(master,pslave)
      deallocate(pslave)
   end if
end subroutine couple_horizontal_variables

recursive subroutine couple_scalar_variables(model,master,slave)
   _CLASS_ (type_model_info),  intent(inout),target :: model
   type (type_scalar_variable),intent(inout),target :: master
   type (type_scalar_variable),intent(in),   target :: slave

   type (type_scalar_variable),     pointer :: pslave
   type (type_scalar_variable_link),pointer :: link
   _CLASS_ (type_model_info),       pointer :: child

   pslave => slave

   if (associated(pslave,master)) return

   ! Process all links and if they used to refer to the specified slave,
   ! redirect them to the specified master.
   link => model%first_scalar_link
   do while (associated(link))
      if (associated(link%target,pslave)) then
         link%target => master
         link%coupled = .true.
      end if
      link => link%next
   end do

   ! Allow child models to do the same.
   child => model%first_child
   do while (associated(child))
      call couple_variables(child,master,pslave)
      child => child%next_sibling
   end do

   if (.not.associated(model%parent)) then
      call merge_variables(master,pslave)
      deallocate(pslave)
   end if
end subroutine couple_scalar_variables

subroutine merge_bulk_variables(master,slave)
   type (type_bulk_variable),intent(inout) :: master
   type (type_bulk_variable),intent(in)    :: slave
   integer :: i
   type (type_conserved_quantity_component), pointer :: component

   call log_message(trim(slave%name)//' --> '//trim(master%name))
   if (_ALLOCATED_(slave%write_indices)) then
      call fatal_error('merge_bulk_variables','Attempt to couple write-only variable ' &
         //trim(slave%name)//' to '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%state_indices).and..not._ALLOCATED_(master%state_indices)) then
      call fatal_error('merge_bulk_variables','Attempt to couple state variable ' &
         //trim(slave%name)//' to non-state variable '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%cons_indices).and..not._ALLOCATED_(master%cons_indices)) then
      call fatal_error('merge_bulk_variables','Attempt to couple conserved quantity ' &
         //trim(slave%name)//' with non-conserved quantity '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(master%cons_indices).and..not._ALLOCATED_(slave%cons_indices)) then
      call fatal_error('merge_bulk_variables','Attempt to couple non-conserved quantity ' &
         //trim(slave%name)//' with conserved quantity '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%alldata)) then
      do i=1,size(slave%alldata)
         call append_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%state_indices)) then
      do i=1,size(slave%state_indices)
         call append_index(master%state_indices,slave%state_indices(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%cons_indices)) then
      do i=1,size(slave%cons_indices)
         call append_index(master%cons_indices,slave%cons_indices(i)%p)
      end do
   end if
#ifdef _FABM_F2003_
   call master%properties%update(slave%properties,overwrite=.false.)
   if (associated(master%components%first)) then
      component => master%components%first
      do while (associated(component%next))
         component => component%next
      end do
      component%next => slave%components%first
   else
      master%components%first => slave%components%first
   end if
#endif
end subroutine merge_bulk_variables

subroutine merge_horizontal_variables(master,slave)
   type (type_horizontal_variable),intent(inout) :: master
   type (type_horizontal_variable),intent(in)    :: slave
   integer :: i
   
   call log_message(trim(slave%name)//' --> '//trim(master%name))
   if (slave%domain/=master%domain) then
      call fatal_error('merge_horizontal_variables','Domains of coupled variabled ' &
         //trim(slave%name)//' to '//trim(master%name)//' do not match.')
   end if
   if (_ALLOCATED_(slave%write_indices)) then
      call fatal_error('merge_horizontal_variables','Attempt to couple write-only variable ' &
         //trim(slave%name)//' to '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%state_indices).and..not._ALLOCATED_(master%state_indices)) then
      call fatal_error('merge_horizontal_variables','Attempt to couple state variable ' &
         //trim(slave%name)//' to non-state variable '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%cons_indices).and..not._ALLOCATED_(master%cons_indices)) then
      call fatal_error('merge_horizontal_variables','Attempt to couple conserved quantity ' &
         //trim(slave%name)//' with non-conserved quantity '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(master%cons_indices).and..not._ALLOCATED_(slave%cons_indices)) then
      call fatal_error('merge_horizontal_variables','Attempt to couple non-conserved quantity ' &
         //trim(slave%name)//' with conserved quantity '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%alldata)) then
      do i=1,size(slave%alldata)
         call append_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%state_indices)) then
      do i=1,size(slave%state_indices)
         call append_index(master%state_indices,slave%state_indices(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%cons_indices)) then
      do i=1,size(slave%cons_indices)
         call append_index(master%cons_indices,slave%cons_indices(i)%p)
      end do
   end if
#ifdef _FABM_F2003_
   call master%properties%update(slave%properties,overwrite=.false.)
#endif
end subroutine merge_horizontal_variables

subroutine merge_scalar_variables(master,slave)
   type (type_scalar_variable),intent(inout) :: master
   type (type_scalar_variable),intent(in)    :: slave
   integer :: i
   
   call log_message(trim(slave%name)//' --> '//trim(master%name))
   if (_ALLOCATED_(slave%write_indices)) then
      call fatal_error('merge_scalar_variables','Attempt to couple write-only variable ' &
         //trim(slave%name)//' to '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%state_indices).and..not._ALLOCATED_(master%state_indices)) then
      call fatal_error('merge_scalar_variables','Attempt to couple state variable ' &
         //trim(slave%name)//' to non-state variable '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%cons_indices).and..not._ALLOCATED_(master%cons_indices)) then
      call fatal_error('merge_scalar_variables','Attempt to couple conserved quantity ' &
         //trim(slave%name)//' with non-conserved quantity '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(master%cons_indices).and..not._ALLOCATED_(slave%cons_indices)) then
      call fatal_error('merge_scalar_variables','Attempt to couple non-conserved quantity ' &
         //trim(slave%name)//' with conserved quantity '//trim(master%name)//'.')
   end if
   if (_ALLOCATED_(slave%alldata)) then
      do i=1,size(slave%alldata)
         call append_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%state_indices)) then
      do i=1,size(slave%state_indices)
         call append_index(master%state_indices,slave%state_indices(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%cons_indices)) then
      do i=1,size(slave%cons_indices)
         call append_index(master%cons_indices,slave%cons_indices(i)%p)
      end do
   end if
#ifdef _FABM_F2003_
   call master%properties%update(slave%properties,overwrite=.false.)
#endif
end subroutine merge_scalar_variables

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on another variable.
!
! !INTERFACE:
   recursive function find_bulk_variable(model,name,exclude) result(variable)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      _CLASS_ (type_model_info),intent(in),target :: model
      character(len=*),         intent(in)        :: name
      type (type_bulk_variable_link),intent(in),target :: exclude
      type (type_bulk_variable),pointer           :: variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),pointer :: curinfo
      type (type_bulk_variable_link),pointer :: link
!
!-----------------------------------------------------------------------
!BOC
      curinfo => model
      do while (associated(curinfo))
         ! Check all model links
         link => curinfo%first_link
         do while (associated(link))
            if (link%name==name.and..not.associated(link,exclude)) then
               variable => link%target
               return
            end if
            link => link%next
         end do

         ! Variable not found - move to model parent.
         curinfo => curinfo%parent
      end do

      ! Variable not found in model tree.
      nullify(variable)
   end function find_bulk_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on another variable.
!
! !INTERFACE:
   recursive function find_horizontal_variable(model,name,exclude) result(variable)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      _CLASS_ (type_model_info),           intent(in),target :: model
      character(len=*),                    intent(in)        :: name
      type (type_horizontal_variable_link),intent(in),target :: exclude
      type (type_horizontal_variable),pointer                :: variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),pointer :: curinfo
      type (type_horizontal_variable_link),pointer :: link
!
!-----------------------------------------------------------------------
!BOC
      curinfo => model
      do while (associated(curinfo))
         ! Check all model links
         link => curinfo%first_horizontal_link
         do while (associated(link))
            if (link%name==name.and..not.associated(link,exclude)) then
               variable => link%target
               return
            end if
            link => link%next
         end do

         ! Variable not found - move to model parent.
         curinfo => curinfo%parent
      end do

      ! Variable not found in model tree.
      nullify(variable)
   end function find_horizontal_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on another variable.
!
! !INTERFACE:
   recursive function find_scalar_variable(model,name,exclude) result(variable)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      _CLASS_ (type_model_info),       intent(in),target :: model
      character(len=*),                intent(in)        :: name
      type (type_scalar_variable_link),intent(in),target :: exclude
      type (type_scalar_variable),pointer                :: variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),pointer :: curinfo
      type (type_scalar_variable_link),pointer :: link
!
!-----------------------------------------------------------------------
!BOC
      curinfo => model
      do while (associated(curinfo))
         ! Check all model links
         link => curinfo%first_scalar_link
         do while (associated(link))
            if (link%name==name.and..not.associated(link,exclude)) then
               variable => link%target
               return
            end if
            link => link%next
         end do

         ! Variable not found - move to model parent.
         curinfo => curinfo%parent
      end do

      ! Variable not found in model tree.
      nullify(variable)
   end function find_scalar_variable
!EOC

function create_external_bulk_id(model,variable) result(id)
   _CLASS_ (type_model_info),intent(in),   target :: model
   type (type_bulk_variable),intent(inout),target :: variable
   type (type_bulk_variable_id) :: id
   id%variable => variable
   if (_ALLOCATED_(variable%alldata)) then
      allocate(id%alldata(size(variable%alldata)))
      id%alldata = variable%alldata
      id%p => id%alldata(1)%p
   else
      allocate(id%alldata(0))
   end if
   if (_ALLOCATED_(variable%state_indices)) id%state_index = variable%state_indices(1)%p
   if (_ALLOCATED_(variable%write_indices)) id%write_index = variable%write_indices(1)%p
end function create_external_bulk_id

function create_external_horizontal_id(model,variable) result(id)
   _CLASS_ (type_model_info),      intent(in),   target :: model
   type (type_horizontal_variable),intent(inout),target :: variable
   type (type_horizontal_variable_id) :: id
   id%variable => variable
   if (_ALLOCATED_(variable%alldata)) then
      allocate(id%alldata(size(variable%alldata)))
      id%alldata = variable%alldata
      id%p => id%alldata(1)%p
   else
      allocate(id%alldata(0))
   end if
   if (_ALLOCATED_(variable%state_indices)) id%state_index = variable%state_indices(1)%p
   if (_ALLOCATED_(variable%write_indices)) id%write_index = variable%write_indices(1)%p
end function create_external_horizontal_id

function create_external_scalar_id(model,variable) result(id)
   _CLASS_ (type_model_info),  intent(in),   target :: model
   type (type_scalar_variable),intent(inout),target :: variable
   type (type_scalar_variable_id) :: id
   id%variable => variable
   if (_ALLOCATED_(variable%alldata)) then
      allocate(id%alldata(size(variable%alldata)))
      id%alldata = variable%alldata
      id%p => id%alldata(1)%p
   else
      allocate(id%alldata(0))
   end if
   if (_ALLOCATED_(variable%state_indices)) id%state_index = variable%state_indices(1)%p
   if (_ALLOCATED_(variable%write_indices)) id%write_index = variable%write_indices(1)%p
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
   _CLASS_ (type_model_info),         intent(in) :: model
   type (type_bulk_standard_variable),intent(in) :: standard_variable
   type (type_bulk_variable_id)                  :: id

   type (type_bulk_variable_link), pointer :: link
   integer                            :: i

   allocate(id%alldata(0))
   link => model%first_link
   do while (associated(link))
      if (compare_standard_variables(link%target%standard_variable,standard_variable).and. &
          .not.is_null_standard_variable(link%target%standard_variable)) then
         if (_ALLOCATED_(link%target%alldata)) then
            do i=1,size(link%target%alldata)
               call append_data_pointer(id%alldata,link%target%alldata(i)%p)
            end do
         end if
      end if
      link => link%next
   end do
   if (size(id%alldata)>0) id%p => id%alldata(1)%p
end function create_external_bulk_id_for_standard_name

function create_external_horizontal_id_for_standard_name(model,standard_variable) result(id)
   _CLASS_ (type_model_info),               intent(in) :: model
   type (type_horizontal_standard_variable),intent(in) :: standard_variable
   type (type_horizontal_variable_id)                  :: id

   type (type_horizontal_variable_link), pointer :: link
   integer                                       :: i

   allocate(id%alldata(0))
   link => model%first_horizontal_link
   do while (associated(link))
      if (compare_standard_variables(link%target%standard_variable,standard_variable).and. &
          .not.is_null_standard_variable(link%target%standard_variable)) then
         if (_ALLOCATED_(link%target%alldata)) then
            do i=1,size(link%target%alldata)
               call append_data_pointer(id%alldata,link%target%alldata(i)%p)
            end do
         end if
      end if
      link => link%next
   end do
   if (size(id%alldata)>0) id%p => id%alldata(1)%p
end function create_external_horizontal_id_for_standard_name

function create_external_scalar_id_for_standard_name(model,standard_variable) result(id)
   _CLASS_ (type_model_info),           intent(in) :: model
   type (type_global_standard_variable),intent(in) :: standard_variable
   type (type_scalar_variable_id)                  :: id

   type (type_scalar_variable_link), pointer :: link
   integer                                   :: i

   allocate(id%alldata(0))
   link => model%first_scalar_link
   do while (associated(link))
      if (compare_standard_variables(link%target%standard_variable,standard_variable).and. &
          .not.is_null_standard_variable(link%target%standard_variable)) then
         if (_ALLOCATED_(link%target%alldata)) then
            do i=1,size(link%target%alldata)
               call append_data_pointer(id%alldata,link%target%alldata(i)%p)
            end do
         end if
      end if
      link => link%next
   end do
   if (size(id%alldata)>0) id%p => id%alldata(1)%p
end function create_external_scalar_id_for_standard_name

recursive subroutine classify_variables(model)
   _CLASS_ (type_model_info),intent(inout),target :: model

   type (type_bulk_variable_link),                 pointer :: link
   type (type_horizontal_variable_link),           pointer :: horizontal_link
   type (type_scalar_variable_link),               pointer :: scalar_link

   type (type_state_variable_info),                pointer :: statevar
   type (type_horizontal_state_variable_info),     pointer :: hz_statevar
   type (type_diagnostic_variable_info),           pointer :: diagvar
   type (type_horizontal_diagnostic_variable_info),pointer :: hz_diagvar
   type (type_conserved_quantity_info),            pointer :: consvar
   integer                                                 :: nstate,nstate_bot,nstate_surf,ndiag,ndiag_hz,ncons

   integer :: i

   _CLASS_ (type_model_info),pointer :: child
   
   if (model%frozen) call fatal_error('fabm_types::classify_variables', &
      'classify_variables may only be called once (during initialization).')

   ! First process child models
   child => model%first_child
   do while (associated(child))
      call classify_variables(child)
      child => child%next_sibling
   end do

   ! Make sure that no one will register any new variables from this moment onward.
   model%frozen = .true.

   allocate(model%dependencies(0))
   allocate(model%dependencies_hz(0))
   allocate(model%dependencies_scalar(0))

   ! Count number of bulk variables in various categories.
   nstate = 0
   ndiag  = 0
   ncons  = 0
   link => model%first_link
   do while (associated(link))
      if (.not.link%coupled) then
         if (_ALLOCATED_(link%target%write_indices)) then
            ndiag = ndiag+1
            do i=1,size(link%target%write_indices)
               link%target%write_indices(i)%p = ndiag
            end do
         end if
         if (_ALLOCATED_(link%target%state_indices)) then
            nstate = nstate+1
            do i=1,size(link%target%state_indices)
               link%target%state_indices(i)%p = nstate
            end do
         end if
         if (_ALLOCATED_(link%target%cons_indices)) then
            ncons = ncons+1
            do i=1,size(link%target%cons_indices)
               link%target%cons_indices(i)%p = ncons
            end do
         end if
      end if
      if (_ALLOCATED_(link%target%alldata)) then
         call append_string(model%dependencies,link%name)
         if (link%target%standard_variable%name/='') &
            call append_string(model%dependencies,link%target%standard_variable%name)
      end if
      link => link%next
   end do

   ! Count number of horizontal variables in various categories.
   nstate_bot  = 0
   nstate_surf = 0
   ndiag_hz    = 0
   horizontal_link => model%first_horizontal_link
   do while (associated(horizontal_link))
      if (.not.horizontal_link%coupled) then
         if (_ALLOCATED_(horizontal_link%target%write_indices)) then
            ndiag_hz = ndiag_hz+1
            do i=1,size(horizontal_link%target%write_indices)
               horizontal_link%target%write_indices(i)%p = ndiag_hz
            end do
         end if
         if (_ALLOCATED_(horizontal_link%target%state_indices)) then
            select case (horizontal_link%target%domain)
               case (domain_bottom)
                  nstate_bot = nstate_bot+1
                  do i=1,size(horizontal_link%target%state_indices)
                     horizontal_link%target%state_indices(i)%p = nstate_bot
                  end do
               case (domain_surface)
                  nstate_surf = nstate_surf+1
                  do i=1,size(horizontal_link%target%state_indices)
                     horizontal_link%target%state_indices(i)%p = nstate_surf
                  end do
            end select
         end if
      end if
      if (_ALLOCATED_(horizontal_link%target%alldata)) then
         call append_string(model%dependencies_hz,horizontal_link%name)
         if (horizontal_link%target%standard_variable%name/='') &
            call append_string(model%dependencies_hz,horizontal_link%target%standard_variable%name)
      end if
      horizontal_link => horizontal_link%next
   end do

   ! Count number of scalar variables in various categories.
   scalar_link => model%first_scalar_link
   do while (associated(scalar_link))
      if (_ALLOCATED_(scalar_link%target%alldata)) then
         call append_string(model%dependencies_scalar,scalar_link%name)
         if (scalar_link%target%standard_variable%name/='') &
            call append_string(model%dependencies_scalar,scalar_link%target%standard_variable%name)
      end if
      scalar_link => scalar_link%next
   end do

   ! Allocate arrays with variable information that will be accessed by the host model.
   allocate(model%state_variables                (nstate))
   allocate(model%bottom_state_variables         (nstate_bot))
   allocate(model%surface_state_variables        (nstate_surf))
   allocate(model%diagnostic_variables           (ndiag))
   allocate(model%horizontal_diagnostic_variables(ndiag_hz))
   allocate(model%conserved_quantities           (ncons))

   ! Set pointers for backward compatibility (pre 2013-06-15)
   model%state_variables_ben => model%bottom_state_variables
   model%diagnostic_variables_hz => model%horizontal_diagnostic_variables

   ! Classify bulk variables
   link => model%first_link
   do while (associated(link))
      if (.not.link%coupled) then
         ! The model owns this variable (no external master variable has been assigned)
         ! Transfer variable information to the array that will be accessed by the host model.
         if (_ALLOCATED_(link%target%write_indices)) then
            diagvar => model%diagnostic_variables(link%target%write_indices(1)%p)
            diagvar%globalid          = create_external_variable_id(model,link%target)
            diagvar%name              = link%target%name
            diagvar%units             = link%target%units
            diagvar%long_name         = link%target%long_name
            diagvar%standard_variable = link%target%standard_variable
            diagvar%minimum           = link%target%minimum
            diagvar%maximum           = link%target%maximum
            diagvar%missing_value     = link%target%missing_value
            diagvar%time_treatment    = link%target%time_treatment
#ifdef _FABM_F2003_
            call diagvar%properties%update(link%target%properties)
#endif
         end if
         
         if (_ALLOCATED_(link%target%state_indices)) then
            statevar => model%state_variables(link%target%state_indices(1)%p)
            statevar%globalid                  = create_external_variable_id(model,link%target)
            statevar%name                      = link%target%name
            statevar%units                     = link%target%units
            statevar%long_name                 = link%target%long_name
            statevar%standard_variable         = link%target%standard_variable
            statevar%minimum                   = link%target%minimum
            statevar%maximum                   = link%target%maximum
            statevar%missing_value             = link%target%missing_value
            statevar%initial_value             = link%target%initial_value
            statevar%vertical_movement         = link%target%vertical_movement
            statevar%specific_light_extinction = link%target%specific_light_extinction
            statevar%no_precipitation_dilution = link%target%no_precipitation_dilution
            statevar%no_river_dilution         = link%target%no_river_dilution
#ifdef _FABM_F2003_
            call statevar%properties%update(link%target%properties)
#endif
         end if
         
         if (_ALLOCATED_(link%target%cons_indices)) then
            consvar => model%conserved_quantities(link%target%cons_indices(1)%p)
            consvar%globalid          = create_external_variable_id(model,link%target)
            consvar%name              = link%target%name
            consvar%units             = link%target%units
            consvar%long_name         = link%target%long_name
            consvar%standard_variable = link%target%standard_variable
#ifdef _FABM_F2003_
            consvar%components  = link%target%components
            call consvar%properties%update(link%target%properties)
#endif
         end if
      end if
      link => link%next
   end do

   ! Classify horizontal variables
   horizontal_link => model%first_horizontal_link
   do while (associated(horizontal_link))
      if (.not.horizontal_link%coupled) then
         ! The model owns this variable (no external master variable has been assigned)
         ! Transfer variable information to the array that will be accessed by the host model.
         if (_ALLOCATED_(horizontal_link%target%write_indices)) then
            hz_diagvar => model%horizontal_diagnostic_variables(horizontal_link%target%write_indices(1)%p)
            hz_diagvar%globalid          = create_external_variable_id(model,horizontal_link%target)
            hz_diagvar%name              = horizontal_link%target%name
            hz_diagvar%units             = horizontal_link%target%units
            hz_diagvar%long_name         = horizontal_link%target%long_name
            hz_diagvar%standard_variable = horizontal_link%target%standard_variable
            hz_diagvar%minimum           = horizontal_link%target%minimum
            hz_diagvar%maximum           = horizontal_link%target%maximum
            hz_diagvar%missing_value     = horizontal_link%target%missing_value
            hz_diagvar%time_treatment    = horizontal_link%target%time_treatment
#ifdef _FABM_F2003_
            call hz_diagvar%properties%update(horizontal_link%target%properties)
#endif
         end if
         if (_ALLOCATED_(horizontal_link%target%state_indices)) then
            select case (horizontal_link%target%domain)
               case (domain_bottom)
                  hz_statevar => model%bottom_state_variables(horizontal_link%target%state_indices(1)%p)
               case (domain_surface)
                  hz_statevar => model%surface_state_variables(horizontal_link%target%state_indices(1)%p)
            end select
            hz_statevar%globalid          = create_external_variable_id(model,horizontal_link%target)
            hz_statevar%name              = horizontal_link%target%name
            hz_statevar%units             = horizontal_link%target%units
            hz_statevar%long_name         = horizontal_link%target%long_name
            hz_statevar%standard_variable = horizontal_link%target%standard_variable
            hz_statevar%minimum           = horizontal_link%target%minimum
            hz_statevar%maximum           = horizontal_link%target%maximum
            hz_statevar%missing_value     = horizontal_link%target%missing_value
            hz_statevar%initial_value     = horizontal_link%target%initial_value
#ifdef _FABM_F2003_
            call hz_statevar%properties%update(horizontal_link%target%properties)
#endif
         end if
      end if
      horizontal_link => horizontal_link%next
   end do

end subroutine classify_variables


!-----------------------------------------------------------------------

   end module fabm_types

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_types --- Derived types used by 0D biogeochemical modules
!
! !INTERFACE:
   module fabm_types

   use fabm_driver, only: fatal_error, log_message

   implicit none
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
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!
   ! Base data type for biogeochemical models.
   public type_model_info

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
   public type_variable_link
   public type_horizontal_variable_link
   public type_scalar_variable_link
   public type_environment
   public initialize_model_info
   public freeze_model_info
   public create_external_variable_id
   
! !PUBLIC DATA MEMBERS:
!
   integer, parameter, public :: attribute_length = 256
   
   integer, parameter, public :: rk = _FABM_REAL_KIND_

   ! Below a list of names of standard physical-biogeochemical variables that have a well-defined interpretation and unit.
   ! These variables can be used from biogeochemical models by calling register_dependency with
   ! the name of the desired variable during model initialization.
   !
   ! Names are based on the Standard Name Table from the NetCDF Climate and Forecast (CF) Metadata Convention.
   ! See http://cf-pcmdi.llnl.gov/documents/cf-standard-names/.
   ! In deriving names from the CF convention, the following exceptions are made to account for the fact
   ! that FABM handles both marine and limnic systems and has the water column as default domain:
   ! - "sea_water_" prefix is suppressed
   ! - "_in_sea_water" suffix is suppressed
   ! - instead of the "_at_sea_floor" suffix a "bottom_" prefix is used, analogous to the "surface_" prefix used in CF.

   ! Variables defined throughout the water column.
   character(len=attribute_length),parameter,public :: &
     varname_temp    = 'temperature',                                      & ! In-situ temperature (degree_Celsius)
     varname_salt    = 'practical_salinity',                               & ! Salinity on Practical Salinity Scale (1e-3)
     varname_swr     = 'downwelling_shortwave_flux',                       & ! Shortwave [200-4000 nm] radiation (W m-2)
     varname_par     = 'downwelling_photosynthetic_radiative_flux',        & ! Photosynthetically Active [400-700 nm] Radiation (W m-2)
     varname_pres    = 'pressure',                                         & ! Pressure (dbar = 10 kPa)
     varname_dens    = 'density'                                             ! In-situ density (kg m-3)

   ! Variables defined on a horizontal surface (e.g., water surface or bottom).
   character(len=attribute_length),parameter,public :: &
     varname_lon     = 'longitude',                                        & ! Longitude (degree_East)
     varname_lat     = 'latitude',                                         & ! Latitude (degree_North)
     varname_wind_sf = 'wind_speed',                                       & ! Wind speed, defined at 10 m above water surface (m s-1)
     varname_cloud   = 'cloud_area_fraction',                              & ! Cloud cover (1), i.e., a fraction between 0 and 1
     varname_swr_sf  = 'downwelling_shortwave_flux_in_air',                & ! Shortwave [200-4000 nm] radiation, defined at water surface (W m-2)
     varname_par_sf  = 'downwelling_photosynthetic_radiative_flux_in_air', & ! Photosynthetically Active [400-700 nm] Radiation, defined at water surface (W m-2)
     varname_zbot    = 'bottom_depth',                                     & ! Basin floor depth below geoid (approx. mean sea level) (m)
     varname_taub    = 'bottom_stress'                                       ! Bottom stress (Pa)

   ! Added for aed modules
   character(len=attribute_length),parameter,public :: &
     varname_layer_ht= 'cell_thickness',                                                       & ! Layer thickness (m)
     varname_extc    = 'attenuation_coefficient_of_downwelling_photosynthetic_radiative_flux', & ! Attenuation coefficient for Photosynthetically Active [400-700 nm] Radiation (m-1)
     varname_tss     = 'mass_concentration_of_suspended_matter',                               & ! Total suspended matter or suspended solids (g m-3)
     varname_sed_zone= 'env_sed_zone'                                                            ! sedimentation zone

   ! Non-spatial (scalar) variables.
   character(len=attribute_length),parameter,public :: &
     varname_yearday = 'number_of_days_since_start_of_the_year'              ! Decimal day of the year (day), equal to 0.0 at 00:00 1 Jan UTC
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
   ! Variable identifiers used by biogeochemical models.
   ! ====================================================================================================
   
   type type_state_variable_id
      character(len=attribute_length) :: name        = ''
      integer                         :: state_index = -1
      type (type_bulk_data_pointer)   :: data
   end type
      
   type type_bottom_state_variable_id
      character(len=attribute_length)     :: name               = ''
      integer                             :: bottom_state_index = -1
      type (type_horizontal_data_pointer) :: horizontal_data
   end type

   type type_diagnostic_variable_id
      character(len=attribute_length) :: name       = ''
      integer                         :: diag_index = -1
   end type

   type type_horizontal_diagnostic_variable_id
      character(len=attribute_length) :: name              = ''
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

   type type_conserved_quantity_id
      character(len=attribute_length) :: name       = ''
      integer                         :: cons_index = -1
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

   type type_bulk_variable
      ! Metadata
      character(len=attribute_length) :: name                      = ''
      character(len=attribute_length) :: long_name                 = ''
      character(len=attribute_length) :: units                     = ''
      real(rk)                        :: minimum                   = -1.e20_rk
      real(rk)                        :: maximum                   =  1.e20_rk
      real(rk)                        :: missing_value             = -2.e20_rk
      real(rk)                        :: initial_value             = 0.0_rk
      real(rk)                        :: vertical_movement         = 0.0_rk
      real(rk)                        :: specific_light_extinction = 0.0_rk
      logical                         :: no_precipitation_dilution = .false.
      logical                         :: no_river_dilution         = .false.
      integer                         :: time_treatment            = time_treatment_last
      
      ! Arrays with all associated data and index pointers.
      type (type_bulk_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata       _NULL_
      type (type_integer_pointer),          dimension(:),_ALLOCATABLE_ :: state_indices _NULL_
      type (type_integer_pointer),          dimension(:),_ALLOCATABLE_ :: write_indices _NULL_
      type (type_integer_pointer),          dimension(:),_ALLOCATABLE_ :: cons_indices  _NULL_
   end type type_bulk_variable

   type type_horizontal_variable
      ! Metadata
      character(len=attribute_length) :: name           = ''
      character(len=attribute_length) :: long_name      = ''
      character(len=attribute_length) :: units          = ''
      real(rk)                        :: minimum        = -1.e20_rk
      real(rk)                        :: maximum        =  1.e20_rk
      real(rk)                        :: missing_value  = -2.e20_rk
      real(rk)                        :: initial_value  = 0.0_rk
      integer                         :: time_treatment = time_treatment_last
      
      ! Arrays with all associated data and index pointers.
      type (type_horizontal_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata       _NULL_
      type (type_integer_pointer),                dimension(:),_ALLOCATABLE_ :: state_indices _NULL_
      type (type_integer_pointer),                dimension(:),_ALLOCATABLE_ :: write_indices _NULL_
      type (type_integer_pointer),                dimension(:),_ALLOCATABLE_ :: cons_indices  _NULL_
   end type type_horizontal_variable

   type type_scalar_variable
      ! Metadata
      character(len=attribute_length) :: name           = ''
      character(len=attribute_length) :: long_name      = ''
      character(len=attribute_length) :: units          = ''
      real(rk)                        :: minimum        = -1.e20_rk
      real(rk)                        :: maximum        =  1.e20_rk
      real(rk)                        :: missing_value  = -2.e20_rk
      integer                         :: time_treatment = time_treatment_last

      ! Arrays with all associated data and index pointers.
      type (type_scalar_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: alldata       _NULL_
      type (type_integer_pointer),            dimension(:),_ALLOCATABLE_ :: state_indices _NULL_
      type (type_integer_pointer),            dimension(:),_ALLOCATABLE_ :: write_indices _NULL_
      type (type_integer_pointer),            dimension(:),_ALLOCATABLE_ :: cons_indices  _NULL_
   end type type_scalar_variable

   type type_variable_link
      character(len=attribute_length)    :: name    = ''
      type (type_bulk_variable), pointer :: target  => null()
      logical                            :: coupled = .false.
      type (type_variable_link), pointer :: next    => null()
   end type type_variable_link

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
   ! Types to hold variable metadata, used by the external host.
   ! ====================================================================================================

!  Derived type describing a state variable
   type type_state_variable_info
      character(len=attribute_length)   :: name                      = ''
      character(len=attribute_length)   :: long_name                 = ''
      character(len=attribute_length)   :: units                     = ''
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

   type type_horizontal_state_variable_info
      character(len=attribute_length)         :: name          = ''
      character(len=attribute_length)         :: long_name     = ''
      character(len=attribute_length)         :: units         = ''
      real(rk)                                :: initial_value = 0.0_rk
      real(rk)                                :: minimum       = -1.e20_rk
      real(rk)                                :: maximum       =  1.e20_rk
      real(rk)                                :: missing_value = -2.e20_rk
      integer                                 :: externalid    = 0      ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      type (type_horizontal_variable_id)      :: globalid
   end type type_horizontal_state_variable_info

!  Derived type describing a diagnostic variable
   type type_diagnostic_variable_info
      character(len=attribute_length)   :: name           = ''
      character(len=attribute_length)   :: long_name      = ''
      character(len=attribute_length)   :: units          = ''
      real(rk)                          :: minimum        = -1.e20_rk
      real(rk)                          :: maximum        =  1.e20_rk
      real(rk)                          :: missing_value  = -2.e20_rk
      integer                           :: externalid     = 0                    ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      integer                           :: time_treatment = time_treatment_last ! Time treatment: 0=last value, 1=time-integrated, 2=time step-averaged, 3=time step-integrated
      type (type_bulk_variable_id)      :: globalid
   end type type_diagnostic_variable_info

   type type_horizontal_diagnostic_variable_info
      character(len=attribute_length)         :: name           = ''
      character(len=attribute_length)         :: long_name      = ''
      character(len=attribute_length)         :: units          = ''
      real(rk)                                :: minimum        = -1.e20_rk
      real(rk)                                :: maximum        =  1.e20_rk
      real(rk)                                :: missing_value  = -2.e20_rk
      integer                                 :: externalid     = 0                   ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      integer                                 :: time_treatment = time_treatment_last ! Time treatment: 0=last value, 1=time-integrated, 2=time step-averaged, 3=time step-integrated
      type (type_horizontal_variable_id)      :: globalid
   end type type_horizontal_diagnostic_variable_info

!  Derived type describing a conserved quantity
   type type_conserved_quantity_info
      character(len=attribute_length)   :: name       = ''
      character(len=attribute_length)   :: long_name  = ''
      character(len=attribute_length)   :: units      = ''
      integer                           :: externalid = 0       ! Identifier to be used by host (e.g., to hold NetCDF identifier)
      type (type_bulk_variable_id)      :: globalid
   end type type_conserved_quantity_info
   
   ! ====================================================================================================
   ! Base model type, used by biogeochemical models to inherit from, and by external host to
   ! get variable lists and metadata.
   ! ====================================================================================================

   type type_model_info
      ! Flag determining whether the contents of the type are "frozen", i.e., they will not change anymore.
      logical :: frozen = .false.

      ! Arrays with metadata on model variables.
      type (type_state_variable_info),                _ALLOCATABLE_,dimension(:) :: state_variables         _NULL_
      type (type_horizontal_state_variable_info),     _ALLOCATABLE_,dimension(:) :: state_variables_ben     _NULL_
      type (type_diagnostic_variable_info),           _ALLOCATABLE_,dimension(:) :: diagnostic_variables    _NULL_
      type (type_horizontal_diagnostic_variable_info),_ALLOCATABLE_,dimension(:) :: diagnostic_variables_hz _NULL_
      type (type_conserved_quantity_info),            _ALLOCATABLE_,dimension(:) :: conserved_quantities    _NULL_

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

      type (type_variable_link),           pointer :: first_link            => null()
      type (type_horizontal_variable_link),pointer :: first_horizontal_link => null()
      type (type_scalar_variable_link),    pointer :: first_scalar_link     => null()
      type (type_named_coupling),          pointer :: first_coupling        => null()

#ifdef _FABM_F2003_
      contains

      ! Procedures that may be used to register model variables and dependencies during initialization.
      procedure :: register_bulk_state_variable             => register_bulk_state_variable
      procedure :: register_bottom_state_variable           => register_bottom_state_variable
      procedure :: register_bulk_diagnostic_variable        => register_bulk_diagnostic_variable
      procedure :: register_horizontal_diagnostic_variable  => register_horizontal_diagnostic_variable
      procedure :: register_bulk_dependency                 => register_bulk_dependency
      procedure :: register_horizontal_dependency           => register_horizontal_dependency
      procedure :: register_global_dependency               => register_global_dependency
      procedure :: register_conserved_quantity              => register_conserved_quantity
      procedure :: register_bulk_state_dependency           => register_bulk_state_dependency
      procedure :: register_bottom_state_dependency         => register_bottom_state_dependency

      generic :: register_state_variable      => register_bulk_state_variable,register_bottom_state_variable
      generic :: register_diagnostic_variable => register_bulk_diagnostic_variable,register_horizontal_diagnostic_variable
      generic :: register_dependency          => register_bulk_dependency,register_horizontal_dependency,register_global_dependency
      generic :: register_state_dependency    => register_bulk_state_dependency,register_bottom_state_dependency

      ! Procedures that may be overridden by biogeochemical models to provide custom data or functionality.
      procedure :: initialize               => initialize_model_info
      procedure :: set_domain               => base_set_domain
      procedure :: do                       => base_do
      procedure :: do_ppdd                  => base_do_ppdd
      procedure :: do_benthos               => base_do_benthos
      procedure :: do_benthos_ppdd          => base_do_benthos_ppdd
      procedure :: get_light_extinction     => base_get_light_extinction
      procedure :: get_drag                 => base_get_drag
      procedure :: get_albedo               => base_get_albedo
      procedure :: get_conserved_quantities => base_get_conserved_quantities
      procedure :: get_surface_exchange     => base_get_surface_exchange
      procedure :: get_vertical_movement    => base_get_vertical_movement
      procedure :: check_state              => base_check_state
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
   
   interface create_external_variable_id
      module procedure create_external_bulk_id
      module procedure create_external_horizontal_id
      module procedure create_external_scalar_id
   end interface
   
   interface register_state_variable
      module procedure register_bulk_state_variable
      module procedure register_bottom_state_variable
   end interface
   
   interface register_state_dependency
      module procedure register_bulk_state_dependency
      module procedure register_bottom_state_dependency
   end interface

   interface register_dependency
      module procedure register_bulk_dependency
      module procedure register_horizontal_dependency
      module procedure register_global_dependency
   end interface
   
   interface register_diagnostic_variable
      module procedure register_bulk_diagnostic_variable
      module procedure register_horizontal_diagnostic_variable
   end interface

!-----------------------------------------------------------------------

   contains

#ifdef _FABM_F2003_
   subroutine base_set_domain(self _ARG_LOCATION_)
      class (type_model_info),intent(inout) :: self
      _DECLARE_LOCATION_ARG_
   end subroutine base_set_domain
   subroutine base_do(self,_FABM_ARGS_DO_RHS_)
      class (type_model_info),intent(in) ::  self
      _DECLARE_FABM_ARGS_DO_RHS_
   end subroutine base_do
   subroutine base_do_ppdd(self,_FABM_ARGS_DO_PPDD_)
      class (type_model_info),intent(in) :: self
      _DECLARE_FABM_ARGS_DO_PPDD_
   end subroutine base_do_ppdd
   subroutine base_get_light_extinction(self,_FABM_ARGS_GET_EXTINCTION_)
      class (type_model_info), intent(in) :: self
      _DECLARE_FABM_ARGS_GET_EXTINCTION_
   end subroutine base_get_light_extinction
   subroutine base_get_drag(self,_FABM_ARGS_GET_DRAG_)
      class (type_model_info), intent(in) :: self
      _DECLARE_FABM_ARGS_GET_DRAG_
   end subroutine base_get_drag
   subroutine base_get_albedo(self,_FABM_ARGS_GET_ALBEDO_)
      class (type_model_info), intent(in) :: self
      _DECLARE_FABM_ARGS_GET_ALBEDO_
   end subroutine base_get_albedo
   subroutine base_get_conserved_quantities(self,_FABM_ARGS_GET_CONSERVED_QUANTITIES_)
      class (type_model_info), intent(in) :: self
      _DECLARE_FABM_ARGS_GET_CONSERVED_QUANTITIES_
   end subroutine base_get_conserved_quantities
   subroutine base_do_benthos(self,_FABM_ARGS_DO_BENTHOS_RHS_)
      class (type_model_info),intent(in) :: self
      _DECLARE_FABM_ARGS_DO_BENTHOS_RHS_
   end subroutine base_do_benthos
   subroutine base_do_benthos_ppdd(self,_FABM_ARGS_DO_BENTHOS_PPDD_)
      class (type_model_info),intent(in) :: self
      _DECLARE_FABM_ARGS_DO_BENTHOS_PPDD_
   end subroutine base_do_benthos_ppdd
   subroutine base_get_surface_exchange(self,_FABM_ARGS_GET_SURFACE_EXCHANGE_)
      class (type_model_info), intent(in) :: self
      _DECLARE_FABM_ARGS_GET_SURFACE_EXCHANGE_
   end subroutine base_get_surface_exchange
   subroutine base_get_vertical_movement(self,_FABM_ARGS_GET_VERTICAL_MOVEMENT_)
      class (type_model_info), intent(in) :: self
      _DECLARE_FABM_ARGS_GET_VERTICAL_MOVEMENT_
   end subroutine base_get_vertical_movement
   subroutine base_check_state(self,_FABM_ARGS_CHECK_STATE_)
      class (type_model_info), intent(in) :: self
      _DECLARE_FABM_ARGS_CHECK_STATE_
   end subroutine base_check_state
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
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),pointer :: last_child
!
!-----------------------------------------------------------------------
!BOC
      model%name             = name
      model%name_prefix      = trim(name)//'_'
      model%long_name_prefix = trim(name)//' '

      if (present(parent)) then
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
      end if

   end subroutine initialize_model_info
!EOC

subroutine new_link(model,target,name,merge)
   _CLASS_ (type_model_info),       intent(inout) :: model
   type (type_bulk_variable),target,intent(in)    :: target
   character(len=*),                intent(in)    :: name
   logical,optional,                intent(in)    :: merge

   type (type_variable_link),pointer :: link

   ! First check if a link with this name exists. If so, merge new target with old target.
   link => model%first_link
   do while (associated(link))
      if (link%name==name) then
         if (.not.present(merge)) call fatal_error('new_link','Link '//trim(name)//' already exists.')
         if (merge) call merge_bulk_variables(link%target,target)
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
end subroutine new_link

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
         if (.not.present(merge)) call fatal_error('new_link','Link '//trim(name)//' already exists.')
         if (merge) call merge_horizontal_variables(link%target,target)
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
         if (.not.present(merge)) call fatal_error('new_link','Link '//trim(name)//' already exists.')
         if (merge) call merge_scalar_variables(link%target,target)
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

subroutine add_index(array,index)
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
end subroutine add_index

subroutine add_data_pointer(array,data)
   type (type_bulk_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: array
   type (type_bulk_data_pointer),target :: data
   type (type_bulk_data_pointer_pointer),allocatable :: oldarray(:)

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not._ALLOCATED_(array)) then
      allocate(array(1))
   else
      allocate(oldarray(size(array)))
      oldarray = array
      deallocate(array)
      allocate(array(size(oldarray)+1))
      array(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   ! Add pointer to provided data to the list.
   array(size(array))%p => data
end subroutine

subroutine add_horizontal_data_pointer(array,data)
   type (type_horizontal_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: array
   type (type_horizontal_data_pointer),target :: data
   type (type_horizontal_data_pointer_pointer),allocatable :: oldarray(:)

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not._ALLOCATED_(array)) then
      allocate(array(1))
   else
      allocate(oldarray(size(array)))
      oldarray = array
      deallocate(array)
      allocate(array(size(oldarray)+1))
      array(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   ! Add pointer to provided data to the list.
   array(size(array))%p => data
end subroutine

subroutine add_scalar_data_pointer(array,data)
   type (type_scalar_data_pointer_pointer),dimension(:),_ALLOCATABLE_ :: array
   type (type_scalar_data_pointer),target :: data
   type (type_scalar_data_pointer_pointer),allocatable :: oldarray(:)

   ! Create a new list of data pointers, or extend it if already allocated.
   if (.not._ALLOCATED_(array)) then
      allocate(array(1))
   else
      allocate(oldarray(size(array)))
      oldarray = array
      deallocate(array)
      allocate(array(size(oldarray)+1))
      array(1:size(oldarray)) = oldarray
      deallocate(oldarray)
   end if

   ! Add pointer to provided data to the list.
   array(size(array))%p => data
end subroutine

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

      call couple_variables(model)
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
                                    no_precipitation_dilution,no_river_dilution,target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical state variable in the global model database.
!  It returns an identifier that may be used later to retrieve the value of the state variable.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),intent(inout)                 :: model
      _TYPE_STATE_VARIABLE_ID_, intent(inout),target          :: id
      type (type_bulk_variable),intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, long_name, units
      real(rk),              intent(in),optional :: initial_value,vertical_movement,specific_light_extinction
      real(rk),              intent(in),optional :: minimum, maximum,missing_value
      logical,               intent(in),optional :: no_precipitation_dilution,no_river_dilution
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
            target                    = curinfo)
      else
         ! Store customized information on state variable.
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         if (present(initial_value))             curinfo%initial_value = initial_value
         if (present(minimum))                   curinfo%minimum = minimum
         if (present(maximum))                   curinfo%maximum = maximum
         if (present(missing_value))             curinfo%missing_value = missing_value
         if (present(vertical_movement))         curinfo%vertical_movement = vertical_movement
         if (present(specific_light_extinction)) curinfo%specific_light_extinction = specific_light_extinction
         if (present(no_precipitation_dilution)) curinfo%no_precipitation_dilution = no_precipitation_dilution
         if (present(no_river_dilution        )) curinfo%no_river_dilution         = no_river_dilution
         call add_index(curinfo%state_indices,id%state_index)
         call add_data_pointer(curinfo%alldata,id%data)
         if (id%name/='') call fatal_error('fabm_types::register_bulk_state_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name

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
                                    initial_value, minimum, maximum, missing_value, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical state variable in the global model database.
!  It returns an identifier that may be used later to retrieve the value of the state variable.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),      intent(inout)                 :: model
      _TYPE_BOTTOM_STATE_VARIABLE_ID_,intent(inout),target          :: id
      type (type_horizontal_variable),intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, long_name, units
      real(rk),              intent(in),optional :: initial_value
      real(rk),              intent(in),optional :: minimum, maximum,missing_value
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
         call register_bottom_state_variable(model%parent,id,trim(model%name_prefix)//trim(name), &
            units,trim(model%long_name_prefix)//' '//trim(long_name),       &
            initial_value             = initial_value,             &
            minimum                   = minimum,                   &
            maximum                   = maximum,                   &
            missing_value             = missing_value,             &
            target                    = curinfo)
      else
         ! Store customized information on state variable.
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         if (present(initial_value))             curinfo%initial_value = initial_value
         if (present(minimum))                   curinfo%minimum = minimum
         if (present(maximum))                   curinfo%maximum = maximum
         if (present(missing_value))             curinfo%missing_value = missing_value
         call add_index(curinfo%state_indices,id%bottom_state_index)
         call add_horizontal_data_pointer(curinfo%alldata,id%horizontal_data)
         if (id%name/='') call fatal_error('fabm_types::register_bottom_state_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name

         ! Ensure that initial value falls within prescribed valid range.
         if (curinfo%initial_value<curinfo%minimum .or. curinfo%initial_value>curinfo%maximum) then
            write (text,*) 'Initial value',curinfo%initial_value,'for variable "'//trim(name)//'" lies&
                  &outside allowed range',curinfo%minimum,'to',curinfo%maximum
            call fatal_error('fabm_types::register_bottom_state_variable',text)
         end if
      end if

      call new_horizontal_link(model%first_horizontal_link,curinfo,name,.not.present(target))

   end subroutine register_bottom_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   recursive subroutine register_bulk_diagnostic_variable(model, id, name, units, long_name, &
                                                     time_treatment, missing_value, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical diagnostic variable in the global model database.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),    intent(inout)                 :: model
      _TYPE_DIAGNOSTIC_VARIABLE_ID_,intent(inout),target          :: id
      type (type_bulk_variable),    intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, long_name, units
      integer, optional,     intent(in)          :: time_treatment
      real(rk),optional,     intent(in)          :: missing_value
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
                                           units,                                                      &
                                           trim(model%long_name_prefix)//' '//trim(long_name),       &
                                           time_treatment=time_treatment,                              &
                                           missing_value=missing_value,                                &
                                           target = curinfo)
      else
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         if (present(time_treatment)) curinfo%time_treatment = time_treatment
         if (present(missing_value))  curinfo%missing_value = missing_value
         call add_index(curinfo%write_indices,id%diag_index)
         if (id%name/='') call fatal_error('fabm_types::register_bulk_diagnostic_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
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
                                                   time_treatment, missing_value, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical diagnostic variable in the global model database.
!
! !INPUT/OUTPUT PARAMETER:
      _CLASS_ (type_model_info),           intent(inout)                 :: model
      _TYPE_HORIZONTAL_DIAGNOSTIC_VARIABLE_ID_,intent(inout),target      :: id
      type (type_horizontal_variable),     intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, long_name, units
      integer, optional,     intent(in)          :: time_treatment
      real(rk),optional,     intent(in)          :: missing_value
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
                                           units,                                                      &
                                           trim(model%long_name_prefix)//' '//trim(long_name),       &
                                           time_treatment=time_treatment,                              &
                                           missing_value=missing_value,                                &
                                           target = curinfo)
      else
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         if (present(time_treatment)) curinfo%time_treatment = time_treatment
         if (present(missing_value))  curinfo%missing_value = missing_value
         call add_index(curinfo%write_indices,id%horizontal_diag_index)
         if (id%name/='') call fatal_error('fabm_types::register_horizontal_diagnostic_variable', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
      end if

      call new_horizontal_link(model%first_horizontal_link,curinfo,name,.not.present(target))

   end subroutine register_horizontal_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new conserved quantity
!
! !INTERFACE:
   recursive subroutine register_conserved_quantity(model, id, name, units, long_name, target)
!
! !DESCRIPTION:
!  This function registers a new biogeochemically conserved quantity in the global
!  model database.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),   intent(inout)                 :: model
      _TYPE_CONSERVED_QUANTITY_ID_,intent(inout),target          :: id
      type (type_bulk_variable),   intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name
      character(len=*),      intent(in)          :: long_name
      character(len=*),      intent(in)          :: units
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
         call register_conserved_quantity(model%parent,id,trim(model%name_prefix)//name, &
                 units,trim(model%long_name_prefix)//' '//long_name,target=curinfo)
      else
         curinfo%name      = name
         curinfo%units     = units
         curinfo%long_name = long_name
         call add_index(curinfo%cons_indices,id%cons_index)
         if (id%name/='') call fatal_error('fabm_types::register_conserved_quantity', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
      end if

      call new_link(model,curinfo,name,.not.present(target))

   end subroutine register_conserved_quantity
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on another biogeochemical state variable
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
      _CLASS_ (type_model_info),              intent(inout) :: model
      _TYPE_STATE_VARIABLE_ID_,target,        intent(inout) :: id
!
! !INPUT PARAMETERS:
      character(len=*),                       intent(in)    :: name
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
! !IROUTINE: Registers a dependency on another biogeochemical state variable
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
      _CLASS_ (type_model_info),              intent(inout) :: model
      _TYPE_BOTTOM_STATE_VARIABLE_ID_,target, intent(inout) :: id
!
! !INPUT PARAMETERS:
      character(len=*),                       intent(in)    :: name
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
! !IROUTINE: Registers a read-only dependency on another variable.
!
! !INTERFACE:
   recursive subroutine register_bulk_dependency(model,id,name,target)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),intent(inout)                 :: model
      _TYPE_DEPENDENCY_ID_,     intent(inout),target          :: id
      type (type_bulk_variable),intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                       intent(in)    :: name
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
      if (model%frozen) call fatal_error('fabm_types::register_bulk_dependency',&
         'Dependencies may only be registered during initialization.')

      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if
      
      if (associated(model%parent)) then
         call register_bulk_dependency(model%parent,id,trim(model%name_prefix)//trim(name),curinfo)
      else
         curinfo%name = name
         call add_data_pointer(curinfo%alldata,id%data)
         if (id%name/='') call fatal_error('fabm_types::register_bulk_dependency', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
      end if

      if (associated(model%parent).and..not.present(target)) call new_coupling(model,name,name)

      call new_link(model,curinfo,name,.not.present(target))

   end subroutine register_bulk_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on another variable.
!
! !INTERFACE:
   recursive subroutine register_horizontal_dependency(model,id,name,target)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),      intent(inout)                 :: model
      _TYPE_HORIZONTAL_DEPENDENCY_ID_,intent(inout),target          :: id
      type (type_horizontal_variable),intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                       intent(in)    :: name
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
      if (model%frozen) call fatal_error('fabm_types::register_horizontal_dependency',&
         'Dependencies may only be registered during initialization.')

      if (present(target)) then
         curinfo => target
      else
         allocate(curinfo)
      end if
      
      if (associated(model%parent)) then
         call register_horizontal_dependency(model%parent,id,trim(model%name_prefix)//trim(name),curinfo)
      else
         curinfo%name = name
         call add_horizontal_data_pointer(curinfo%alldata,id%horizontal_data)
         if (id%name/='') call fatal_error('fabm_types::register_horizontal_dependency', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
      end if

      if (associated(model%parent).and..not.present(target)) call new_coupling(model,name,name)

      call new_horizontal_link(model%first_horizontal_link,curinfo,name,.not.present(target))

   end subroutine register_horizontal_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on another variable.
!
! !INTERFACE:
   recursive subroutine register_global_dependency(model,id,name,target)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !INPUT/OUTPUT PARAMETERS:
      _CLASS_ (type_model_info),  intent(inout)                 :: model
      _TYPE_GLOBAL_DEPENDENCY_ID_,intent(inout),target          :: id
      type (type_scalar_variable),intent(inout),target,optional :: target
!
! !INPUT PARAMETERS:
      character(len=*),                       intent(in)    :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
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
         call register_global_dependency(model%parent,id,trim(model%name_prefix)//trim(name),curinfo)
      else
         curinfo%name = name
         call add_scalar_data_pointer(curinfo%alldata,id%global_data)
         if (id%name/='') call fatal_error('fabm_types::register_global_dependency', &
            'Identifier supplied for '//trim(name)//' is already used by '//trim(id%name)//'.')
         id%name = name
      end if

      if (associated(model%parent).and..not.present(target)) call new_coupling(model,name,name)

      call new_scalar_link(model%first_scalar_link,curinfo,name,.not.present(target))

   end subroutine register_global_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on another variable.
!
! !INTERFACE:
   recursive subroutine couple_variables(model)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      _CLASS_ (type_model_info),intent(inout),target   :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),           pointer :: root,child

      type (type_named_coupling),          pointer :: coupling

      type (type_variable_link),           pointer :: link
      type (type_bulk_variable),           pointer :: bulk_master,bulk_slave

      type (type_horizontal_variable_link),pointer :: horizontal_link
      type (type_horizontal_variable)     ,pointer :: horizontal_master,horizontal_slave

      type (type_scalar_variable_link),    pointer :: scalar_link
      type (type_scalar_variable),         pointer :: scalar_master,scalar_slave
!
!-----------------------------------------------------------------------
!BOC
      ! Find root model, with which any unresolved target variables will be registered.
      root => model
      do while (associated(root%parent))
         root => root%parent
      end do
      
      ! Enumerate named couplings, locate slave and target variables, and couple them.
      coupling => model%first_coupling
      do while (associated(coupling))

         ! Try to find slave variable amongst bulk variables.
         link => model%first_link
         do while (associated(link))
            if (link%name==coupling%slave) then
               bulk_slave => link%target
               bulk_master => find_bulk_variable(model,coupling%master,link)
               if (.not.associated(bulk_master)) then
                  if ((_ALLOCATED_(bulk_slave%alldata).and..not._ALLOCATED_(bulk_slave%state_indices))) then
                     ! Read-only dependency cannot be resolved in the FABM model tree. Create a placeholder variable at the root level.
                     allocate(bulk_master)
                     bulk_master%name = coupling%master
                     call new_link(root,bulk_master,coupling%master)
                  else
                     call fatal_error('couple_variables','Unable to find target variable ' &
                        //trim(coupling%master)//' for coupled variable '//trim(bulk_slave%name)//'.')
                  end if
               end if
               call couple_bulk_variables(root,bulk_master,bulk_slave)
            end if
            link => link%next
         end do

         ! Try to find slave variable amongst horizontal variables.
         horizontal_link => model%first_horizontal_link
         do while (associated(horizontal_link))
            if (horizontal_link%name==coupling%slave) then
               horizontal_slave => horizontal_link%target
               horizontal_master => find_horizontal_variable(model,coupling%master,horizontal_link)
               if (.not.associated(horizontal_master)) then
                  if (_ALLOCATED_(horizontal_slave%alldata).and..not._ALLOCATED_(horizontal_slave%state_indices)) then
                     ! Read-only dependency cannot be resolved in the FABM model tree. Create a placeholder variable at the root level.
                     allocate(horizontal_master)
                     horizontal_master%name = coupling%master
                     call new_horizontal_link(root%first_horizontal_link,horizontal_master,coupling%master)
                  else
                     call fatal_error('couple_variables','Unable to find target variable ' &
                        //trim(coupling%master)//' for coupled variable '//trim(horizontal_slave%name)//'.')
                  end if
               end if
               call couple_horizontal_variables(root,horizontal_master,horizontal_slave)
            end if
            horizontal_link => horizontal_link%next
         end do

         ! Try to find slave variable amongst scalar variables.
         scalar_link => model%first_scalar_link
         do while (associated(scalar_link))
            if (scalar_link%name==coupling%slave) then
               scalar_master => find_scalar_variable(model,coupling%master,scalar_link)
               scalar_slave => scalar_link%target
               if (.not.associated(scalar_master)) then
                  if ((_ALLOCATED_(scalar_slave%alldata).and..not._ALLOCATED_(scalar_slave%state_indices))) then
                     ! Read-only dependency that cannot be resolved in the FABM model tree. Create a placeholder variable at the root level.
                     allocate(scalar_master)
                     scalar_master%name = coupling%master
                     call new_scalar_link(root%first_scalar_link,scalar_master,coupling%master)
                  else
                     call fatal_error('couple_variables','Unable to find target variable ' &
                        //trim(coupling%master)//' for coupled variable '//trim(scalar_slave%name)//'.')
                  end if
               end if
               call couple_scalar_variables(root,scalar_master,scalar_slave)
            end if
            scalar_link => scalar_link%next
         end do

         coupling => coupling%next
      end do

      ! Process coupling commands registered with child models.
      child => model%first_child
      do while (associated(child))
         call couple_variables(child)
         child => child%next_sibling
      end do
      
   end subroutine couple_variables
!EOC

recursive subroutine couple_bulk_variables(model,master,slave)
   _CLASS_ (type_model_info),intent(inout),target :: model
   type (type_bulk_variable),pointer              :: master,slave

   type (type_variable_link),           pointer :: link
   _CLASS_ (type_model_info),           pointer :: child

   if (.not.associated(master)) call fatal_error('couple_bulk_variables', &
      'Attempt to couple variable '//trim(slave%name)//' to unknown master variable.')
      
   if (associated(master,slave)) return

   ! Process all links and if they used to refer to the specified slave,
   ! redirect them to the specified master.
   link => model%first_link
   do while (associated(link))
      if (associated(link%target,slave)) then
         link%target => master
         link%coupled = .true.
      end if
      link => link%next
   end do
   
   ! Allow child models to do the same.
   child => model%first_child
   do while (associated(child))
      call couple_bulk_variables(child,master,slave)
      child => child%next_sibling
   end do
   
   if (.not.associated(model%parent)) then
      call merge_bulk_variables(master,slave)
      deallocate(slave)
   end if
end subroutine couple_bulk_variables

recursive subroutine couple_horizontal_variables(model,master,slave)
   _CLASS_ (type_model_info),intent(inout),target :: model
   type (type_horizontal_variable),       pointer :: master,slave

   type (type_horizontal_variable_link),pointer :: link
   _CLASS_ (type_model_info),           pointer :: child

   if (.not.associated(master)) call fatal_error('couple_horizontal_variables', &
      'Attempt to couple variable '//trim(slave%name)//' to unknown master variable.')
   if (associated(master,slave).or..not.associated(master)) return

   ! Process all links and if they used to refer to the specified slave,
   ! redirect them to the specified master.
   link => model%first_horizontal_link
   do while (associated(link))
      if (associated(link%target,slave)) then
         link%target => master
         link%coupled = .true.
      end if
      link => link%next
   end do
   
   ! Allow child models to do the same.
   child => model%first_child
   do while (associated(child))
      call couple_horizontal_variables(child,master,slave)
      child => child%next_sibling
   end do
   
   if (.not.associated(model%parent)) then
      call merge_horizontal_variables(master,slave)
      deallocate(slave)
   end if
end subroutine couple_horizontal_variables

recursive subroutine couple_scalar_variables(model,master,slave)
   _CLASS_ (type_model_info),intent(inout),target :: model
   type (type_scalar_variable),           pointer :: master,slave

   type (type_scalar_variable_link),    pointer :: link
   _CLASS_ (type_model_info),           pointer :: child

   if (.not.associated(master)) call fatal_error('couple_scalar_variables', &
      'Attempt to couple variable '//trim(slave%name)//' to unknown master variable.')
   if (associated(master,slave).or..not.associated(master)) return

   ! Process all links and if they used to refer to the specified slave,
   ! redirect them to the specified master.
   link => model%first_scalar_link
   do while (associated(link))
      if (associated(link%target,slave)) then
         link%target => master
         link%coupled = .true.
      end if
      link => link%next
   end do
   
   ! Allow child models to do the same.
   child => model%first_child
   do while (associated(child))
      call couple_scalar_variables(child,master,slave)
      child => child%next_sibling
   end do
   
   if (.not.associated(model%parent)) then
      call merge_scalar_variables(master,slave)
      deallocate(slave)
   end if
end subroutine couple_scalar_variables

subroutine merge_bulk_variables(master,slave)
   type (type_bulk_variable),intent(inout) :: master
   type (type_bulk_variable),intent(in)    :: slave
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
         call add_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%state_indices)) then
      do i=1,size(slave%state_indices)
         call add_index(master%state_indices,slave%state_indices(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%cons_indices)) then
      do i=1,size(slave%cons_indices)
         call add_index(master%cons_indices,slave%cons_indices(i)%p)
      end do
   end if
end subroutine merge_bulk_variables

subroutine merge_horizontal_variables(master,slave)
   type (type_horizontal_variable),intent(inout) :: master
   type (type_horizontal_variable),intent(in)    :: slave
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
         call add_horizontal_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%state_indices)) then
      do i=1,size(slave%state_indices)
         call add_index(master%state_indices,slave%state_indices(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%cons_indices)) then
      do i=1,size(slave%cons_indices)
         call add_index(master%cons_indices,slave%cons_indices(i)%p)
      end do
   end if
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
         call add_scalar_data_pointer(master%alldata,slave%alldata(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%state_indices)) then
      do i=1,size(slave%state_indices)
         call add_index(master%state_indices,slave%state_indices(i)%p)
      end do
   end if
   if (_ALLOCATED_(slave%cons_indices)) then
      do i=1,size(slave%cons_indices)
         call add_index(master%cons_indices,slave%cons_indices(i)%p)
      end do
   end if
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
      type (type_variable_link),intent(in),target :: exclude
      type (type_bulk_variable),pointer           :: variable
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      _CLASS_ (type_model_info),pointer :: curinfo
      type (type_variable_link),pointer :: link
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

function create_external_bulk_id(variable) result(id)
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

function create_external_horizontal_id(variable) result(id)
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

function create_external_scalar_id(variable) result(id)
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

recursive subroutine classify_variables(model)
   _CLASS_ (type_model_info),intent(inout),target :: model

   type (type_variable_link),                      pointer :: link
   type (type_horizontal_variable_link),           pointer :: horizontal_link
   type (type_scalar_variable_link),               pointer :: scalar_link

   type (type_state_variable_info),                pointer :: statevar
   type (type_horizontal_state_variable_info),     pointer :: hz_statevar
   type (type_diagnostic_variable_info),           pointer :: diagvar
   type (type_horizontal_diagnostic_variable_info),pointer :: hz_diagvar
   type (type_conserved_quantity_info),            pointer :: consvar
   integer                                                 :: nstate,nstate_hz,ndiag,ndiag_hz,ncons,nread,nread_hz,nread_scalar
   
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
   
   ! Count number of bulk variables in various categories.
   nstate = 0
   ndiag  = 0
   ncons  = 0
   nread  = 0
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
         if (_ALLOCATED_(link%target%alldata)) nread = nread+1
      end if
      link => link%next
   end do

   ! Count number of horizontal variables in various categories.
   nstate_hz = 0
   ndiag_hz  = 0
   nread_hz  = 0
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
            nstate_hz = nstate_hz+1
            do i=1,size(horizontal_link%target%state_indices)
               horizontal_link%target%state_indices(i)%p = nstate_hz
            end do
         end if
         if (_ALLOCATED_(horizontal_link%target%alldata)) nread_hz = nread_hz+1
      end if
      horizontal_link => horizontal_link%next
   end do

   ! Count number of scalar variables in various categories.
   nread_scalar = 0
   scalar_link => model%first_scalar_link
   do while (associated(scalar_link))
      if (.not.scalar_link%coupled) then
         if (_ALLOCATED_(scalar_link%target%alldata)) nread_scalar = nread_scalar+1
      end if
      scalar_link => scalar_link%next
   end do

   ! Allocate arrays with variable information that will be accessed by the host model.
   allocate(model%state_variables        (nstate))
   allocate(model%state_variables_ben    (nstate_hz))
   allocate(model%diagnostic_variables   (ndiag))
   allocate(model%diagnostic_variables_hz(ndiag_hz))
   allocate(model%conserved_quantities   (ncons))
   allocate(model%dependencies           (nread))
   allocate(model%dependencies_hz        (nread_hz))
   allocate(model%dependencies_scalar    (nread_scalar))

   ! Classify bulk variables
   nread = 0
   link => model%first_link
   do while (associated(link))
      if (.not.link%coupled) then
         ! The model owns this variable (no external master variable has been assigned)
         ! Transfer variable information to the array that will be accessed by the host model.
         if (_ALLOCATED_(link%target%write_indices)) then
            diagvar => model%diagnostic_variables(link%target%write_indices(1)%p)
            diagvar%globalid       = create_external_variable_id(link%target)
            diagvar%name           = link%target%name
            diagvar%units          = link%target%units
            diagvar%long_name      = link%target%long_name
            diagvar%minimum        = link%target%minimum
            diagvar%maximum        = link%target%maximum
            diagvar%missing_value  = link%target%missing_value
            diagvar%time_treatment = link%target%time_treatment
         end if
         
         if (_ALLOCATED_(link%target%state_indices)) then
            statevar => model%state_variables(link%target%state_indices(1)%p)
            statevar%globalid                  = create_external_variable_id(link%target)
            statevar%name                      = link%target%name
            statevar%units                     = link%target%units
            statevar%long_name                 = link%target%long_name
            statevar%minimum                   = link%target%minimum
            statevar%maximum                   = link%target%maximum
            statevar%missing_value             = link%target%missing_value
            statevar%initial_value             = link%target%initial_value
            statevar%vertical_movement         = link%target%vertical_movement
            statevar%specific_light_extinction = link%target%specific_light_extinction
            statevar%no_precipitation_dilution = link%target%no_precipitation_dilution
            statevar%no_river_dilution         = link%target%no_river_dilution
         end if
         
         if (_ALLOCATED_(link%target%cons_indices)) then
            consvar => model%conserved_quantities(link%target%cons_indices(1)%p)
            consvar%globalid    = create_external_variable_id(link%target)
            consvar%name        = link%target%name
            consvar%units       = link%target%units
            consvar%long_name   = link%target%long_name
         end if

         if (_ALLOCATED_(link%target%alldata)) then
            nread = nread+1
            model%dependencies(nread) = link%target%name
         end if
      end if
      link => link%next
   end do

   ! Classify horizontal variables
   nread_hz = 0
   do while (associated(horizontal_link))
      if (.not.horizontal_link%coupled) then
         ! The model owns this variable (no external master variable has been assigned)
         ! Transfer variable information to the array that will be accessed by the host model.
         if (_ALLOCATED_(horizontal_link%target%write_indices)) then
            hz_diagvar => model%diagnostic_variables_hz(horizontal_link%target%write_indices(1)%p)
            hz_diagvar%globalid       = create_external_variable_id(horizontal_link%target)
            hz_diagvar%name           = horizontal_link%target%name
            hz_diagvar%units          = horizontal_link%target%units
            hz_diagvar%long_name      = horizontal_link%target%long_name
            hz_diagvar%minimum        = horizontal_link%target%minimum
            hz_diagvar%maximum        = horizontal_link%target%maximum
            hz_diagvar%missing_value  = horizontal_link%target%missing_value
            hz_diagvar%time_treatment = horizontal_link%target%time_treatment
         end if
         if (_ALLOCATED_(horizontal_link%target%state_indices)) then
            hz_statevar => model%state_variables_ben(horizontal_link%target%state_indices(1)%p)
            hz_statevar%globalid      = create_external_variable_id(horizontal_link%target)
            hz_statevar%name          = horizontal_link%target%name
            hz_statevar%units         = horizontal_link%target%units
            hz_statevar%long_name     = horizontal_link%target%long_name
            hz_statevar%minimum       = horizontal_link%target%minimum
            hz_statevar%maximum       = horizontal_link%target%maximum
            hz_statevar%missing_value = horizontal_link%target%missing_value
            hz_statevar%initial_value = horizontal_link%target%initial_value
         end if
         if (_ALLOCATED_(horizontal_link%target%alldata)) then
            nread_hz = nread_hz+1
            model%dependencies_hz(nread_hz) = horizontal_link%target%name
         end if
      end if
      horizontal_link => horizontal_link%next
   end do

   ! Classify scalar variables
   nread_scalar = 0
   scalar_link => model%first_scalar_link
   do while (associated(scalar_link))
      if (.not.scalar_link%coupled) then
         if (_ALLOCATED_(scalar_link%target%alldata)) then
            nread_scalar = nread_scalar+1
            model%dependencies_scalar(nread_scalar) = scalar_link%target%name
         end if
      end if
      scalar_link => scalar_link%next
   end do

end subroutine classify_variables


!-----------------------------------------------------------------------

   end module fabm_types

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

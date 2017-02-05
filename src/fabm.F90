#include "fabm_driver.h"
#include "fabm_private.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: FABM --- Framework for Aquatic Biogeochemical Models
!
! !INTERFACE:
   module fabm
!
! !DESCRIPTION:
! This is the core module of FABM, serving as the "glue layer" between a
! physical host model (e.g., a general circulation model), and one or more
! specific biogeochemical models. A physical host model will call the interfaces
! of this module to access biogeochemistry.
!
! For more information, see the documentation at http://fabm.net/wiki.
!
! To add new biogeochemical models, edit fabm_library.F90.
!
! !USES:
   use fabm_standard_variables,only: type_bulk_standard_variable, type_horizontal_standard_variable, &
                                     type_global_standard_variable, initialize_standard_variables, type_standard_variable_node
   use fabm_types
   use fabm_library
   use fabm_expressions
   use fabm_driver
   use fabm_properties
   use fabm_builtin_models
   use fabm_coupling
   use fabm_job
!

   implicit none
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public fabm_initialize_library, type_model, fabm_create_model_from_file
   public fabm_initialize, fabm_finalize, fabm_set_domain, fabm_check_ready, fabm_update_time
   public fabm_initialize_state, fabm_initialize_surface_state, fabm_initialize_bottom_state

   public fabm_process_job_all

   ! Process rates and diagnostics for pelagic, surface, bottom.
   public fabm_do, fabm_do_surface, fabm_do_bottom

   ! Vertical movement, light attenuation, feedbacks to drag and albedo
   public fabm_get_vertical_movement, fabm_get_light_extinction, fabm_get_albedo, fabm_get_drag, fabm_get_light

   ! Bookkeeping
   public fabm_check_state, fabm_check_surface_state, fabm_check_bottom_state
   public fabm_get_conserved_quantities, fabm_get_horizontal_conserved_quantities

   ! Management of model variables: retrieve identifiers, get and set data.
   public fabm_get_bulk_variable_id,fabm_get_horizontal_variable_id,fabm_get_scalar_variable_id
   public fabm_get_variable_name, fabm_is_variable_used, fabm_variable_needs_values
   public fabm_link_interior_state_data, fabm_link_bottom_state_data, fabm_link_surface_state_data
   public fabm_link_interior_data, fabm_link_horizontal_data, fabm_link_scalar_data
   public fabm_get_interior_diagnostic_data, fabm_get_horizontal_diagnostic_data

   ! For backward compatibility (pre 11 Dec 2015)
   public fabm_link_bulk_state_data, fabm_get_bulk_diagnostic_data, fabm_link_bulk_data

#ifdef _HAS_MASK_
   ! Set spatial mask
   public fabm_set_mask
#endif

   ! For backward compatibility only (use fabm_do_surface and fabm_do_bottom instead)
   public fabm_get_surface_exchange, fabm_do_benthos

   ! Object with all supported standard variables as its members.
   ! Imported from fabm_types, and made available so hosts only need to "use fabm"
   public standard_variables

   ! Variable identifier types by external physical drivers.
   public type_bulk_variable_id
   public type_horizontal_variable_id
   public type_scalar_variable_id
   public type_external_variable,type_horizontal_state_variable_info

   integer, parameter :: state_none             = 0
   integer, parameter :: state_initialize_done  = 1
   integer, parameter :: state_set_domain_done  = 2
   integer, parameter :: state_check_ready_done = 3

   integer, parameter, public :: data_source_none = 0
   integer, parameter, public :: data_source_host = 1
   integer, parameter, public :: data_source_fabm = 2
   integer, parameter, public :: data_source_user = 3
   integer, parameter, public :: data_source_default = data_source_host

!
! !PUBLIC TYPES:
!

   ! ====================================================================================================
   ! Variable identifiers used by host models.
   ! ====================================================================================================

   type type_external_variable_id
      type (type_internal_variable),pointer :: variable => null()
   end type

   type,extends(type_external_variable_id) :: type_bulk_variable_id
   end type

   type,extends(type_external_variable_id) :: type_horizontal_variable_id
   end type

   type,extends(type_external_variable_id) :: type_scalar_variable_id
   end type

   ! ====================================================================================================
   ! Derived types for variable metadata used by host models.
   ! ====================================================================================================

   type,abstract :: type_external_variable
      character(len=attribute_length) :: name          = ''
      character(len=attribute_length) :: long_name     = ''
      character(len=attribute_length) :: local_long_name = ''
      character(len=attribute_length) :: units         = ''
      character(len=attribute_length) :: path          = ''
      real(rk)                        :: minimum       = -1.e20_rk
      real(rk)                        :: maximum       =  1.e20_rk
      real(rk)                        :: missing_value = -2.e20_rk
      integer                         :: output        = output_instantaneous ! See output_* parameters above
      type (type_property_dictionary) :: properties
      integer                         :: externalid    = 0                    ! Identifier to be used freely by host
      type (type_internal_variable),pointer :: target => null()
   end type

!  Derived type describing a state variable
   type,extends(type_external_variable) :: type_state_variable_info
      type (type_bulk_standard_variable) :: standard_variable
      real(rk)                           :: initial_value             = 0.0_rk
      real(rk)                           :: vertical_movement         = 0.0_rk  ! Vertical movement (m/s, <0: sinking, >0: floating)
      logical                            :: no_precipitation_dilution = .false.
      logical                            :: no_river_dilution         = .false.
      integer                            :: sms_index          = -1
      integer                            :: surface_flux_index = -1
      integer                            :: bottom_flux_index  = -1
      integer                            :: movement_index     = -1
   end type type_state_variable_info

   type,extends(type_external_variable) :: type_horizontal_state_variable_info
      type (type_horizontal_standard_variable) :: standard_variable
      real(rk)                                 :: initial_value = 0.0_rk
      integer                                  :: sms_index = -1
   end type type_horizontal_state_variable_info

!  Derived type describing a diagnostic variable
   type,extends(type_external_variable) :: type_diagnostic_variable_info
      type (type_bulk_standard_variable) :: standard_variable
      integer                            :: time_treatment = time_treatment_last ! See time_treatment_* parameters above
      logical                            :: save           = .false.
      integer                            :: source
   end type type_diagnostic_variable_info

   type,extends(type_external_variable) :: type_horizontal_diagnostic_variable_info
      type (type_horizontal_standard_variable) :: standard_variable
      integer                                  :: time_treatment = time_treatment_last ! See time_treatment_* parameters above
      logical                                  :: save           = .false.
      integer                                  :: source
   end type type_horizontal_diagnostic_variable_info

!  Derived type describing a conserved quantity
   type,extends(type_external_variable) :: type_conserved_quantity_info
      type (type_bulk_standard_variable)    :: standard_variable
      integer                               :: index              = -1
      integer                               :: horizontal_index   = -1
      type (type_internal_variable),pointer :: target_hz => null()
   end type type_conserved_quantity_info

   type type_interior_data_pointer
      real(rk),pointer _DIMENSION_GLOBAL_ :: p => null()
   end type

   type type_horizontal_data_pointer
      real(rk),pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: p => null()
   end type

   type type_scalar_data_pointer
      real(rk),pointer :: p => null()
   end type

   ! Derived type for a single generic biogeochemical model
   type type_model
      type (type_base_model) :: root
      type (type_model_list) :: models

      integer :: state = state_none
      type (type_variable_register) :: variable_register

      class (type_model),pointer :: info => null()  ! For backward compatibility (hosts pre 11/2013); always points to root.

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

      integer                                :: extinction_index = -1
      type (type_internal_variable), pointer :: extinction_target => null()

      type (type_link_list) :: links_postcoupling

      ! Declare the arrays for diagnostic variable values.
      real(rk),allocatable _DIMENSION_GLOBAL_PLUS_1_            :: diag
      real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: diag_hz
      real(rk),allocatable                                      :: diag_missing_value(:)
      real(rk),allocatable                                      :: diag_hz_missing_value(:)

      real(rk),allocatable _DIMENSION_GLOBAL_                   :: zero
      real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_        :: zero_hz
      integer                                                   :: domain_size(_FABM_DIMENSION_COUNT_)
      integer                                                   :: horizontal_domain_size(_HORIZONTAL_DIMENSION_COUNT_)

      type (type_job_manager) :: job_manager

      type (type_job) :: do_interior_job
      type (type_job) :: do_bottom_job
      type (type_job) :: do_surface_job
      type (type_job) :: get_vertical_movement_job
      type (type_job) :: get_conserved_quantities_job
      type (type_job) :: get_horizontal_conserved_quantities_job
      type (type_job) :: get_light_extinction_job
      type (type_job) :: get_drag_job
      type (type_job) :: get_albedo_job
      type (type_job) :: get_diagnostics_job
      type (type_job) :: check_state_job
      type (type_job) :: check_bottom_state_job
      type (type_job) :: check_surface_state_job
      type (type_job) :: initialize_state_job
      type (type_job) :: initialize_bottom_state_job
      type (type_job) :: initialize_surface_state_job

      ! Registry with pointers to global fields of readable variables.
      ! These pointers are accessed to fill the read cache just before individual model instances are called.
      type (type_interior_data_pointer),  allocatable :: data(:)
      type (type_horizontal_data_pointer),allocatable :: data_hz(:)
      type (type_scalar_data_pointer),    allocatable :: data_scalar(:)
      integer, allocatable :: interior_data_sources(:)
      integer, allocatable :: horizontal_data_sources(:)
      integer, allocatable :: scalar_data_sources(:)

#ifdef _HAS_MASK_
#  ifndef _FABM_HORIZONTAL_MASK_
      _FABM_MASK_TYPE_,pointer _DIMENSION_GLOBAL_ :: mask => null()
#  endif
      _FABM_MASK_TYPE_,pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: mask_hz => null()
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==0
      integer :: bottom_index = -1
#  elif _FABM_BOTTOM_INDEX_==-1
      integer,pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: bottom_indices => null()
#  endif
      integer :: surface_index = -1
#endif

      integer :: nscratch = -1
      integer :: nscratch_hz = -1
   contains
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      procedure :: set_bottom_index => fabm_set_bottom_index
      procedure :: set_surface_index => fabm_set_surface_index
#endif

      procedure :: link_interior_data_by_variable => fabm_link_interior_data_by_variable
      procedure :: link_interior_data_by_id   => fabm_link_interior_data_by_id
      procedure :: link_interior_data_by_sn   => fabm_link_interior_data_by_sn
      procedure :: link_interior_data_by_name => fabm_link_interior_data_by_name
      generic :: link_interior_data => link_interior_data_by_variable,link_interior_data_by_id,link_interior_data_by_sn,link_interior_data_by_name

      procedure :: link_horizontal_data_by_variable => fabm_link_horizontal_data_by_variable
      procedure :: link_horizontal_data_by_id       => fabm_link_horizontal_data_by_id
      procedure :: link_horizontal_data_by_sn       => fabm_link_horizontal_data_by_sn
      procedure :: link_horizontal_data_by_name     => fabm_link_horizontal_data_by_name
      generic :: link_horizontal_data => link_horizontal_data_by_variable,link_horizontal_data_by_id,link_horizontal_data_by_sn,link_horizontal_data_by_name

      procedure :: link_scalar_by_id   => fabm_link_scalar_by_id
      procedure :: link_scalar_by_sn   => fabm_link_scalar_by_sn
      procedure :: link_scalar_by_name => fabm_link_scalar_by_name
      generic :: link_scalar => link_scalar_by_id,link_scalar_by_sn,link_scalar_by_name

      procedure :: link_interior_state_data => fabm_link_interior_state_data
      procedure :: link_bottom_state_data   => fabm_link_bottom_state_data
      procedure :: link_surface_state_data  => fabm_link_surface_state_data
      procedure :: link_all_interior_state_data => fabm_link_all_interior_state_data
      procedure :: link_all_bottom_state_data   => fabm_link_all_bottom_state_data
      procedure :: link_all_surface_state_data  => fabm_link_all_surface_state_data

      procedure :: require_interior_data => fabm_require_interior_data
      generic :: require_data => require_interior_data

      procedure :: get_interior_data => fabm_get_interior_data
      procedure :: get_horizontal_data => fabm_get_horizontal_data
      procedure :: get_scalar_data => fabm_get_scalar_data
      generic :: get_data => get_interior_data,get_horizontal_data,get_scalar_data

      procedure :: get_bulk_variable_id_by_name => fabm_get_bulk_variable_id_by_name
      procedure :: get_bulk_variable_id_sn => fabm_get_bulk_variable_id_sn
      generic :: get_bulk_variable_id => get_bulk_variable_id_by_name, get_bulk_variable_id_sn

      procedure :: get_horizontal_variable_id_by_name => fabm_get_horizontal_variable_id_by_name
      procedure :: get_horizontal_variable_id_sn => fabm_get_horizontal_variable_id_sn
      generic :: get_horizontal_variable_id => get_horizontal_variable_id_by_name, get_horizontal_variable_id_sn

      procedure :: get_scalar_variable_id_by_name => fabm_get_scalar_variable_id_by_name
      procedure :: get_scalar_variable_id_sn => fabm_get_scalar_variable_id_sn
      generic :: get_scalar_variable_id => get_scalar_variable_id_by_name, get_scalar_variable_id_sn

      procedure :: interior_variable_needs_values => fabm_interior_variable_needs_values
      procedure :: interior_variable_needs_values_sn => fabm_interior_variable_needs_values_sn
      procedure :: horizontal_variable_needs_values => fabm_horizontal_variable_needs_values
      procedure :: horizontal_variable_needs_values_sn => fabm_horizontal_variable_needs_values_sn
      procedure :: scalar_variable_needs_values => fabm_scalar_variable_needs_values
      procedure :: scalar_variable_needs_values_sn => fabm_scalar_variable_needs_values_sn
      generic :: variable_needs_values => interior_variable_needs_values, interior_variable_needs_values_sn, &
                                          horizontal_variable_needs_values, horizontal_variable_needs_values_sn, &
                                          scalar_variable_needs_values, scalar_variable_needs_values_sn

      ! -----------------------------------------------------------------------------
      ! For backward compatibility (pre 11 Dec 2015)
      procedure :: link_bulk_data_by_variable => fabm_link_interior_data_by_variable
      procedure :: link_bulk_data_by_id   => fabm_link_interior_data_by_id
      procedure :: link_bulk_data_by_sn   => fabm_link_interior_data_by_sn
      procedure :: link_bulk_data_by_name => fabm_link_interior_data_by_name
      generic :: link_bulk_data => link_interior_data_by_variable,link_interior_data_by_id,link_interior_data_by_sn,link_interior_data_by_name

      procedure :: bulk_variable_needs_values => fabm_interior_variable_needs_values
      procedure :: bulk_variable_needs_values_sn => fabm_interior_variable_needs_values_sn
      ! -----------------------------------------------------------------------------

   end type type_model

   type type_integer_list_node
      integer :: value
      type (type_integer_list_node),pointer :: next => null()
   end type

   type,extends(type_base_model) :: type_host_container
      type (type_integer_list_node), pointer :: first => null()
   end type

!
! !PUBLIC INTERFACES:
!
   ! Subroutine calculating local temporal derivatives either as a right-hand side vector,
   ! or production/destruction matrices.
   interface fabm_do
      module procedure fabm_do_rhs
      module procedure fabm_do_ppdd
   end interface

   ! Subroutine calculating local temporal derivatives of bottom layer (benthos & pelagic)
   ! either as a right-hand side vector, or production/destruction matrices.
   interface fabm_do_bottom
      module procedure fabm_do_bottom_rhs
      module procedure fabm_do_bottom_ppdd
   end interface

   interface fabm_link_data
      module procedure fabm_link_interior_data_by_id
      module procedure fabm_link_horizontal_data_by_id
      module procedure fabm_link_scalar_by_id
      module procedure fabm_link_interior_data_by_sn
      module procedure fabm_link_horizontal_data_by_sn
      module procedure fabm_link_scalar_by_sn
   end interface

   ! Subroutine for providing FABM with variable data on the full spatial domain.
   interface fabm_link_interior_data
      module procedure fabm_link_interior_data_by_variable
      module procedure fabm_link_interior_data_by_id
      module procedure fabm_link_interior_data_by_sn
      module procedure fabm_link_interior_data_by_name
   end interface

   ! Subroutine for providing FABM with variable data on horizontal slices of the domain.
   interface fabm_link_horizontal_data
      module procedure fabm_link_horizontal_data_by_variable
      module procedure fabm_link_horizontal_data_by_id
      module procedure fabm_link_horizontal_data_by_sn
      module procedure fabm_link_horizontal_data_by_name
   end interface

   ! Subroutine for providing FABM with scalar variable data.
   interface fabm_link_scalar_data
      module procedure fabm_link_scalar_by_id
      module procedure fabm_link_scalar_by_sn
      module procedure fabm_link_scalar_by_name
   end interface

   interface fabm_get_variable_id
      module procedure fabm_get_bulk_variable_id_sn
      module procedure fabm_get_horizontal_variable_id_sn
      module procedure fabm_get_scalar_variable_id_sn
   end interface

   interface fabm_get_bulk_variable_id
      module procedure fabm_get_bulk_variable_id_by_name
      module procedure fabm_get_bulk_variable_id_sn
   end interface

   interface fabm_get_horizontal_variable_id
      module procedure fabm_get_horizontal_variable_id_by_name
      module procedure fabm_get_horizontal_variable_id_sn
   end interface

   interface fabm_get_scalar_variable_id
      module procedure fabm_get_scalar_variable_id_by_name
      module procedure fabm_get_scalar_variable_id_sn
   end interface

   interface fabm_variable_needs_values
      module procedure fabm_interior_variable_needs_values
      module procedure fabm_horizontal_variable_needs_values
   end interface

   ! For backward compatibility only:
   interface fabm_do_benthos
      module procedure fabm_do_bottom_rhs
      module procedure fabm_do_bottom_ppdd
   end interface
   interface fabm_get_surface_exchange
      module procedure fabm_do_surface
   end interface

   ! Fr backward compatibility (pre 11 Dec 2015)
   ! Subroutine for providing FABM with variable data on the full spatial domain.
   interface fabm_link_bulk_data
      module procedure fabm_link_interior_data_by_variable
      module procedure fabm_link_interior_data_by_id
      module procedure fabm_link_interior_data_by_sn
      module procedure fabm_link_interior_data_by_name
   end interface
   interface fabm_link_bulk_state_data
      module procedure fabm_link_interior_state_data
   end interface fabm_link_bulk_state_data
   interface fabm_get_bulk_diagnostic_data
      module procedure fabm_get_interior_diagnostic_data
   end interface fabm_get_bulk_diagnostic_data
!
!EOP
!-----------------------------------------------------------------------

   contains

   subroutine fabm_initialize_library()
      ! Do nothing if already initialized.
      if (associated(factory)) return

      ! If needed, create default object for communication (e.g., logging, error reporting) with host.
      if (.not.associated(driver)) allocate(type_base_driver::driver)

      ! Create all standard variable objects.
      call initialize_standard_variables()

      ! Create the model factory.
      factory => fabm_model_factory
      call factory%initialize()
   end subroutine fabm_initialize_library

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Create a new model tree from a configuration file.
!
! !INTERFACE:
   function fabm_create_model_from_file(file_unit,file,do_not_initialize) result(model)
!
! !INPUT PARAMETERS:
   character(len=*),optional,intent(in) :: file
   integer,                  intent(in) :: file_unit
   logical,optional,         intent(in) :: do_not_initialize
   type (type_model),pointer            :: model

! !LOCAL VARIABLES:
   logical                   :: isopen,initialize
   character(len=256)        :: file_eff
   integer                   :: i,j,modelcount,ownindex
   character(len=64)         :: models(256),instancename
   class (type_base_model),pointer :: childmodel
   logical,parameter         :: alwayspostfixindex=.false.
   namelist /fabm_nml/ models
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_initialize_library()

   nullify(model)

   ! Determine whether the provided unit has been opened already.
   inquire(file_unit,opened=isopen)

   if (.not.isopen) then
      ! Unit has not been openend - we need to open the configuration file ourselves.
      if (present(file)) then
         ! A file path has been provided - use that.
         file_eff = file
      else
         ! No file path has been provided - use default.
         file_eff = 'fabm.nml'
      end if

      ! Open configuration file.
      open(file_unit,file=file_eff,action='read',status='old',err=98)
   end if

   ! Read main FABM namelist.
   models = ''
   read(file_unit,nml=fabm_nml,err=99,end=100)

   ! Create model tree
   allocate(model)
   do i=1,size(models)
      if (models(i)/='') then
         ! Determine if this model name is used multiple times.
         modelcount = 0
         do j=1,size(models)
            if (models(i)==models(j)) then
               modelcount = modelcount + 1
               if (i==j) ownindex = modelcount
            end if
         end do

         ! If another model uses this name too, append a number to the model name.
         if (alwayspostfixindex .or. modelcount>1) then
            write (unit=instancename,fmt='(a,i2.2)') trim(models(i)),ownindex
         else
            instancename = models(i)
         end if

         ! Ask the factory to create the model.
         call factory%create(trim(models(i)),childmodel)
         if (.not.associated(childmodel)) call fatal_error('fabm_create_model_from_file', &
            '"'//trim(models(i))//'" is not a valid model name.')
         childmodel%user_created = .true.

         call log_message('Initializing '//trim(instancename)//'...')
         call model%root%add_child(childmodel,instancename,configunit=file_unit)
         call log_message( '   initialization succeeded.')
      end if
   end do

   ! Initialize model tree [this freezes the model configuration - no new child models or variables can be added]
   initialize = .not.present(do_not_initialize)
   if (.not.initialize) initialize = .not.do_not_initialize
   if (initialize) call fabm_initialize(model)

   ! If we have opened the configuration file ourselves, close it.
   if (.not.isopen) close(file_unit)

   return

98 call fatal_error('fabm_create_model_from_file','Unable to open FABM configuration file '//trim(file_eff)//'.')
   return

99 call fatal_error('fabm_create_model_from_file','Unable to read namelist "fabm_nml".')
   return

100 call fatal_error('fabm_create_model_from_file','Unable to find namelist "fabm_nml".')
   return

   end function fabm_create_model_from_file
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the biogeochemical model tree.
!
! !INTERFACE:
   subroutine fabm_initialize(self)
!
! !INPUT PARAMETERS:
      class (type_model),target,intent(inout) :: self
!
! !LOCAL VARIABLES:
      type (type_aggregate_variable_access),    pointer :: aggregate_variable_access
      class (type_property),                    pointer :: property => null()
      integer                                           :: islash
      type (type_link),                         pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
      if (self%state>=state_initialize_done) &
         call fatal_error('fabm_initialize','fabm_initialize has already been called on this model object.')

      self%info => self ! For backward compatibility (pre 11/2013 hosts only)

      ! Make sure a variable for light extinction is created at the root level when calling freeze_model_info.
      ! This variable is used from fabm_get_light_extinction.
      aggregate_variable_access => get_aggregate_variable_access(self%root, &
         standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      aggregate_variable_access%interior = ior(aggregate_variable_access%interior,access_read)

      ! Create placeholder variables for zero fields.
      ! Values for these fields will only be provided if actually used by one of the biogeochemical models.
      call self%root%add_interior_variable('zero',act_as_state_variable=.true.,source=source_none)
      call self%root%add_horizontal_variable('zero_hz',act_as_state_variable=.true.,source=source_none)

      ! Filter out expressions that FABM can handle itself.
      ! The remainder, if any, must be handled by the host model.
      call filter_expressions(self)

      ! This will resolve all FABM dependencies and generate final authorative lists of variables of different types.
      call freeze_model_info(self%root)

      ! Raise error for unused coupling commands.
      property => self%root%couplings%first
      do while (associated(property))
         if (.not.self%root%couplings%retrieved%contains(trim(property%name))) then
            islash = index(property%name,'/',.true.)
            call fatal_error('fabm_initialize','model '//property%name(1:islash-1)//' does not contain variable "'//trim(property%name(islash+1:))//'" mentioned in coupling section.')
         end if
         property => property%next
      end do

      ! Build final authorative arrays with variable metadata .
      call classify_variables(self)

      ! Create built-in jobs, which can then be chained by the host/user by calling job%set_next.
      ! (the reason for chaining is to allow later jobs to use results of earlier ones, thus reducing the number of calls needed)
      self%job_manager%default_dependency_handler => self%get_diagnostics_job
      call self%job_manager%create(self%do_interior_job,'do_interior',final_operation=source_do)
      call self%job_manager%create(self%do_surface_job,'do_surface',final_operation=source_do_surface)
      call self%job_manager%create(self%do_bottom_job,'do_bottom',final_operation=source_do_bottom)
      call self%job_manager%create(self%get_vertical_movement_job,'get_vertical_movement',final_operation=source_do,ignore_dependencies=.true.)
      call self%job_manager%create(self%get_conserved_quantities_job,'get_conserved_quantities',final_operation=source_do)
      call self%job_manager%create(self%get_horizontal_conserved_quantities_job,'get_horizontal_conserved_quantities',final_operation=source_do_horizontal)
      call self%job_manager%create(self%get_light_extinction_job,'get_light_extinction',final_operation=source_do,ignore_dependencies=.true.)
      call self%job_manager%create(self%get_diagnostics_job,'get_diagnostics_job',outsource_tasks=.true.)
      call self%job_manager%create(self%initialize_state_job,'initialize_state',final_operation=source_do,ignore_dependencies=.true.)
      call self%job_manager%create(self%initialize_bottom_state_job,'initialize_bottom_state',final_operation=source_do_bottom,ignore_dependencies=.true.)
      call self%job_manager%create(self%initialize_surface_state_job,'initialize_surface_state',final_operation=source_do_surface,ignore_dependencies=.true.)
      call self%job_manager%create(self%check_state_job,'check_state',final_operation=source_do,ignore_dependencies=.true.)
      call self%job_manager%create(self%check_bottom_state_job,'check_bottom_state',final_operation=source_do_bottom,ignore_dependencies=.true.)
      call self%job_manager%create(self%check_surface_state_job,'check_surface_state',final_operation=source_do_surface,ignore_dependencies=.true.)
      call self%job_manager%create(self%get_albedo_job,'get_albedo',final_operation=source_do_surface,ignore_dependencies=.true.)
      call self%job_manager%create(self%get_drag_job,'get_drag',final_operation=source_do_surface,ignore_dependencies=.true.)

      call require_flux_computation(self%do_bottom_job,self%links_postcoupling,domain_bottom)
      call require_flux_computation(self%do_surface_job,self%links_postcoupling,domain_surface)
      call require_flux_computation(self%do_interior_job,self%links_postcoupling,domain_interior)

      ! For vertical movement rates, call all models that access interior state variables,
      ! and explicitly express interest in the movement diagnostics so they will be properly prefilled.
      call require_call_all_with_state(self%get_vertical_movement_job,self%root%links,domain_interior,source_get_vertical_movement)
      link => self%links_postcoupling%first
      do while (associated(link))
         if (associated(link%target%movement_diagnostic)) call self%get_vertical_movement_job%request_variable(link%target%movement_diagnostic%target)
         link => link%next
      end do

      call require_call_all_with_state(self%initialize_state_job,self%root%links,domain_interior,source_initialize_state)
      call require_call_all_with_state(self%initialize_bottom_state_job,self%root%links,domain_bottom,source_initialize_bottom_state)
      call require_call_all_with_state(self%initialize_surface_state_job,self%root%links,domain_surface,source_initialize_surface_state)
      call require_call_all_with_state(self%check_state_job,self%root%links,domain_interior,source_check_state)
      call require_call_all_with_state(self%check_bottom_state_job,self%root%links,domain_bottom,source_check_bottom_state)
      call require_call_all_with_state(self%check_bottom_state_job,self%root%links,domain_interior,source_check_bottom_state)
      call require_call_all_with_state(self%check_surface_state_job,self%root%links,domain_surface,source_check_surface_state)
      call require_call_all_with_state(self%check_surface_state_job,self%root%links,domain_interior,source_check_surface_state)
      call require_call_all(self%get_albedo_job,self%root,source_get_albedo)
      call require_call_all(self%get_drag_job,self%root,source_get_drag)

      self%state = state_initialize_done

   end subroutine fabm_initialize
!EOC

   subroutine filter_expressions(self)
      class (type_model),intent(inout)    :: self
      class (type_expression),    pointer :: current,previous,next
      class (type_depth_integral),pointer :: integral
      logical                             :: filter

      nullify(previous)
      current => self%root%first_expression
      do while (associated(current))
         filter = .false.
         select type (current)
            class is (type_vertical_integral)
               allocate(integral)
               integral%minimum_depth = current%minimum_depth
               integral%maximum_depth = current%maximum_depth
               integral%average       = current%average
               call self%root%add_child(integral,trim(current%output_name)//'_calculator',configunit=-1)
               call integral%request_coupling(integral%id_input,current%input_name)
               call self%root%request_coupling(current%output_name,integral%id_output%link%target%name)
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
   end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   recursive subroutine test_presence(self,model)
!
! !INPUT PARAMETERS:
      class (type_model),     intent(in) :: self
      class (type_base_model),intent(in) :: model
!
! !LOCAL VARIABLES:
      type (type_model_list_node),pointer :: child
      character(len=4) :: strcount
!
!EOP
!-----------------------------------------------------------------------
!BOC
      child => model%children%first
      do while (associated(child))
         if (self%models%count(child%model)/=1) then
            write (strcount,'(i0)') self%models%count(child%model)
            call fatal_error('fabm_initialize::test_presence', &
               'BUG: Model "'//trim(child%model%get_path())//'" is not called exactly one time, but '//trim(strcount)//' times .')
         end if
         call test_presence(self,child%model)
         child => child%next
      end do
   end subroutine test_presence
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Deallocate all data associated with the model object.
!
! !INTERFACE:
   subroutine fabm_finalize(self)
!
! !INPUT PARAMETERS:
   class (type_model),target,intent(inout) :: self
!
!EOP
!-----------------------------------------------------------------------
!BOC
   nullify(self%info)
   self%state = state_none

   ! Deallocate the list of models (this does not deallocate the models themselves!)
   call self%models%finalize()

   ! TODO: this should deallocate the memory of all biogeochemical models

   end subroutine fabm_finalize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with the extents of the spatial domain.
!
! !INTERFACE:
   subroutine fabm_set_domain(self _ARGUMENTS_LOCATION_,seconds_per_time_unit)
!
! !INPUT PARAMETERS:
   class (type_model),target,intent(inout) :: self
   _DECLARE_ARGUMENTS_LOCATION_
   real(rk),optional,        intent(in)    :: seconds_per_time_unit
!
! !LOCAL VARIABLES:
  integer                            :: ivar
  class (type_expression),   pointer :: expression
  integer                            :: n
  type (type_variable_node), pointer :: variable_node
!EOP
!-----------------------------------------------------------------------
!BOC
   if (self%state<state_initialize_done) call fatal_error('fabm_set_domain','fabm_initialize has not yet been called on this model object.')
   if (self%state>=state_set_domain_done) call fatal_error('fabm_set_domain','fabm_set_domain has already been called on this model object.')
   self%state = state_set_domain_done

#if _FABM_DIMENSION_COUNT_>0
   self%domain_size = (/ _LOCATION_ /)
#endif
#if _HORIZONTAL_DIMENSION_COUNT_>0
   self%horizontal_domain_size = (/ _HORIZONTAL_LOCATION_ /)
#endif

   ! Default chaining (temporary; should be done explicitly by host if true)
   call self%do_surface_job%set_previous(self%do_bottom_job)
   call self%do_interior_job%set_previous(self%do_surface_job)

   ! Create job that ensures all diagnostics required by the user are computed.
   call self%get_diagnostics_job%set_previous(self%do_interior_job)
   do ivar=1,size(self%diagnostic_variables)
      if (self%diagnostic_variables(ivar)%save) then
         select case (self%diagnostic_variables(ivar)%target%source)
         case (source_check_state)
            call self%check_state_job%request_variable(self%diagnostic_variables(ivar)%target,copy_to_store=.true.)
         case (source_get_vertical_movement)
            call self%get_vertical_movement_job%request_variable(self%diagnostic_variables(ivar)%target,copy_to_store=.true.)
         case default
            call self%get_diagnostics_job%request_variable(self%diagnostic_variables(ivar)%target,copy_to_store=.true.)
         end select
      end if
   end do
   do ivar=1,size(self%horizontal_diagnostic_variables)
      if (self%horizontal_diagnostic_variables(ivar)%save) &
         call self%get_diagnostics_job%request_variable(self%horizontal_diagnostic_variables(ivar)%target,copy_to_store=.true.)
   end do

   ! Merge write indices when operations can be done in place
   ! This must be done after all variables are requested from the different jobs, so we know which variables
   ! will be retrieved (such variables cannot be merged)
   call merge_indices(self%root)

   ! Note: aggregate variables below are copied to read cache, just in case they are not diagnostics
   ! (they may be equal to a state variable, or to zero)
   do ivar=1,size(self%conserved_quantities)
      call self%get_conserved_quantities_job%request_variable(self%conserved_quantities(ivar)%target,copy_to_cache=.true.)
      call self%get_horizontal_conserved_quantities_job%request_variable(self%conserved_quantities(ivar)%target_hz,copy_to_cache=.true.)
   end do
   call self%get_light_extinction_job%request_variable(self%extinction_target,copy_to_cache=.true.)

   ! Assign write cache indices to all interior diagnostic variables.
   ! Must be done after calls to merge_aggregating_diagnostics.
   self%nscratch = 0
   call assign_write_indices(self%links_postcoupling,domain_interior,self%nscratch)

   ! Assign write cache indices to all horizontal diagnostic variables.
   ! Must be done after calls to merge_aggregating_diagnostics.
   self%nscratch_hz = 0
   call assign_write_indices(self%links_postcoupling,domain_surface,   self%nscratch_hz)
   call assign_write_indices(self%links_postcoupling,domain_horizontal,self%nscratch_hz)
   call assign_write_indices(self%links_postcoupling,domain_bottom,    self%nscratch_hz)

   if (any(self%state_variables(:)%sms_index<=0)) call fatal_error('fabm_set_domain','BUG: sms_index invalid for one or more interior state variables.')
   if (any(self%state_variables(:)%surface_flux_index<=0)) call fatal_error('fabm_set_domain','BUG: surface_flux_index invalid for one or more interior state variables.')
   if (any(self%state_variables(:)%bottom_flux_index<=0)) call fatal_error('fabm_set_domain','BUG: bottom_flux_index invalid for one or more interior state variables.')
   if (any(self%state_variables(:)%movement_index<=0)) call fatal_error('fabm_set_domain','BUG: movement_index invalid for one or more interior state variables.')
   if (any(self%surface_state_variables(:)%sms_index<=0)) call fatal_error('fabm_set_domain','BUG: sms_index invalid for one or more surface state variables.')
   if (any(self%bottom_state_variables(:)%sms_index<=0)) call fatal_error('fabm_set_domain','BUG: sms_index invalid for one or more bottom state variables.')

   ! Initialize all jobs - must be done after the call to filter_readable_variable_registry because that
   ! finalizes the set of variables for which input data is available.
   call self%job_manager%initialize(self%variable_register)

   !call self%job_manager%print()

   ! Allocate arrays with pointers to data.
   allocate(self%data       (self%variable_register%interior_read%count))
   allocate(self%data_hz    (self%variable_register%horizontal_read%count))
   allocate(self%data_scalar(self%variable_register%scalar_read%count))

   ! Allocate and initialize arrays that store the source (host, fabm, user) of all data.
   allocate(self%interior_data_sources  (size(self%data       )))
   allocate(self%horizontal_data_sources(size(self%data_hz    )))
   allocate(self%scalar_data_sources    (size(self%data_scalar)))
   self%interior_data_sources   = data_source_none
   self%horizontal_data_sources = data_source_none
   self%scalar_data_sources     = data_source_none

   ! Provide fields with zeros.
   ! This must be done after read indices are assigned (self%job_manager%initialize), otherwise link_*_data will not store the value!
   allocate(self%zero _INDEX_LOCATION_)
   self%zero = 0.0_rk
   call self%link_interior_data('zero',self%zero)

   allocate(self%zero_hz _INDEX_HORIZONTAL_LOCATION_)
   self%zero_hz = 0.0_rk
   call self%link_horizontal_data('zero_hz',self%zero_hz)

   allocate(self%diag(self%domain_size(1),0:self%variable_register%interior_store%count))
   self%diag = 0.0_rk

   allocate(self%diag_hz(0:self%variable_register%horizontal_store%count))
   self%diag_hz = 0.0_rk

   ! Create arrays with missing values for all stored diagnostics (used as fill value for masked cells)
   allocate(self%diag_missing_value(self%variable_register%interior_store%count))
   allocate(self%diag_hz_missing_value(self%variable_register%horizontal_store%count))

   ! Initialize diagnostic variables to missing value.
   n = 0
   variable_node => self%variable_register%interior_store%first
   do while (associated(variable_node))
      n = n + 1
      self%diag_missing_value(n) = variable_node%target%missing_value
      self%diag(_PREARG_LOCATION_DIMENSIONS_ n) = self%diag_missing_value(n)
      call fabm_link_interior_data(self,variable_node%target,self%diag(_PREARG_LOCATION_DIMENSIONS_ n),source=data_source_fabm)
      variable_node => variable_node%next
   end do

   n = 0
   variable_node => self%variable_register%horizontal_store%first
   do while (associated(variable_node))
      n = n + 1
      self%diag_hz_missing_value(n) = variable_node%target%missing_value
      self%diag_hz(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ n) = self%diag_hz_missing_value(n)
      call fabm_link_horizontal_data(self,variable_node%target,self%diag_hz(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ n),source=data_source_fabm)
      variable_node => variable_node%next
   end do

   call self%variable_register%print()

   if (present(seconds_per_time_unit)) then
      expression => self%root%first_expression
      do while (associated(expression))
         select type (expression)
            class is (type_bulk_temporal_mean)
               expression%period = expression%period/seconds_per_time_unit
               allocate(expression%history(_PREARG_LOCATION_ expression%n+3))
               expression%history = 0.0_rk
               call fabm_link_interior_data(self,expression%output_name, &
                                        expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+3))
            class is (type_horizontal_temporal_mean)
               expression%period = expression%period/seconds_per_time_unit
               allocate(expression%history(_PREARG_HORIZONTAL_LOCATION_ expression%n+3))
               expression%history = 0.0_rk
               call fabm_link_horizontal_data(self,expression%output_name, &
                                              expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+3))
         end select
         expression => expression%next
      end do
   end if

   end subroutine fabm_set_domain
!EOC

   recursive subroutine merge_indices(model)
      class (type_base_model), intent(inout) :: model

      type (type_model_list_node), pointer :: child

      select type (model)
      class is (type_weighted_sum)
         call model%reindex()
      class is (type_horizontal_weighted_sum)
         call model%reindex()
      end select

      ! Process children
      child => model%children%first
      do while (associated(child))
         call merge_indices(child%model)
         child => child%next
      end do
   end subroutine merge_indices

#ifdef _HAS_MASK_
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with the spatial mask.
!
! !INTERFACE:
#  ifdef _FABM_HORIZONTAL_MASK_
   subroutine fabm_set_mask(self, mask_hz)
#  else
   subroutine fabm_set_mask(self, mask, mask_hz)
#  endif
!
! !INPUT PARAMETERS:
   class (type_model),target,intent(inout)                            :: self
#  ifndef _FABM_HORIZONTAL_MASK_
   _FABM_MASK_TYPE_, target, intent(in) _DIMENSION_GLOBAL_            :: mask
#  endif
   _FABM_MASK_TYPE_, target, intent(in) _DIMENSION_GLOBAL_HORIZONTAL_ :: mask_hz
!
! !LOCAL VARIABLES:
   integer :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (self%state<state_set_domain_done) &
      call fatal_error('fabm_set_mask','fabm_set_domain has not yet been called on this model object.')

#  ifndef _FABM_HORIZONTAL_MASK_
#    if !defined(NDEBUG)&&_FABM_DIMENSION_COUNT_>0
   do i=1,size(self%domain_size)
      if (size(mask,i)/=self%domain_size(i)) &
         call fatal_error('fabm_set_mask','dimensions of FABM domain and provided mask do not match.')
   end do
#    endif
   self%mask => mask
#  endif

#  if !defined(NDEBUG)&&_HORIZONTAL_DIMENSION_COUNT_>0
   do i=1,size(self%horizontal_domain_size)
      if (size(mask_hz,i)/=self%horizontal_domain_size(i)) &
         call fatal_error('fabm_set_mask','dimensions of FABM domain and provided horizontal mask do not match.')
   end do
#  endif
   self%mask_hz => mask_hz

   end subroutine fabm_set_mask
!EOC
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==-0
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with the vertical indices of bottommost cells.
!
! !INTERFACE:
   subroutine fabm_set_bottom_index(self, index)
!
! !INPUT PARAMETERS:
   class (type_model),target,intent(inout) :: self
   integer,                  intent(in)    :: index
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (self%state<state_set_domain_done) &
      call fatal_error('fabm_set_bottom_index','fabm_set_domain has not yet been called on this model object.')
   if (index<1) &
      call fatal_error('set_bottom_index','provided index must equal or exceed 1.')
   if (index>self%domain_size(_FABM_DEPTH_DIMENSION_INDEX_)) &
      call fatal_error('set_bottom_index','provided index exceeds size of the depth dimension.')

   self%bottom_index = index

   end subroutine fabm_set_bottom_index
!EOC
#  elif _FABM_BOTTOM_INDEX_==-1
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with the vertical indices of bottommost cells.
!
! !INTERFACE:
   subroutine fabm_set_bottom_index(self, indices)
!
! !INPUT PARAMETERS:
   class (type_model),target,intent(inout)                            :: self
   integer,           target,intent(in) _DIMENSION_GLOBAL_HORIZONTAL_ :: indices
!
! !LOCAL VARIABLES:
   integer :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (self%state<state_set_domain_done) &
      call fatal_error('fabm_set_bottom_index','fabm_set_domain has not yet been called on this model object.')
#    if !defined(NDEBUG)&&_HORIZONTAL_DIMENSION_COUNT_>0
   do i=1,size(self%horizontal_domain_size)
      if (size(indices,i)/=self%horizontal_domain_size(i)) &
         call fatal_error('set_bottom_index','dimensions of FABM domain and provided index field do not match.')
   end do
#    endif

   self%bottom_indices => indices

   end subroutine fabm_set_bottom_index
!EOC
#  endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with the vertical indices of bottommost cells.
!
! !INTERFACE:
   subroutine fabm_set_surface_index(self, index)
!
! !INPUT PARAMETERS:
   class (type_model),target,intent(inout) :: self
   integer,                  intent(in)    :: index
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (self%state<state_set_domain_done) &
      call fatal_error('set_surface_index','fabm_set_domain has not yet been called on this model object.')
   if (index<1) &
      call fatal_error('set_surface_index','provided index must equal or exceed 1.')
   if (index>self%domain_size(_FABM_DEPTH_DIMENSION_INDEX_)) &
      call fatal_error('set_surface_index','provided index exceeds size of the depth dimension.')

   self%surface_index = index

   end subroutine fabm_set_surface_index
!EOC
#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check whether FABM has been provided with all required data.
!
! !INTERFACE:
   subroutine fabm_check_ready(self)
!
! !INPUT PARAMETERS:
   class (type_model),intent(inout),target :: self
!
! !LOCAL VARIABLES:
  logical                           :: ready
  type (type_variable_set)          :: unfulfilled_dependencies
  type (type_variable_node),pointer :: variable_node
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (self%state<state_set_domain_done) &
      call fatal_error('fabm_check_ready','fabm_set_domain has not yet been called on this model object.')

   ready = .true.

#ifdef _HAS_MASK_
#  ifndef _FABM_HORIZONTAL_MASK_
   if (.not.associated(self%mask)) then
      call log_message('spatial mask has not been set. Make sure to call fabm_set_mask.')
      ready = .false.
   end if
#  endif
   if (.not.associated(self%mask_hz)) then
      call log_message('horizontal spatial mask has not been set. Make sure to call fabm_set_mask.')
      ready = .false.
   end if
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==0
   if (self%bottom_index==-1) then
      call log_message('bottom index has not been set. Make sure to call set_bottom_index.')
      ready = .false.
   end if
#  elif _FABM_BOTTOM_INDEX_==-1
   if (.not.associated(self%bottom_indices)) then
      call log_message('bottom indices has not been set. Make sure to call set_bottom_index.')
      ready = .false.
   end if
#  endif
   if (self%surface_index==-1) then
      call log_message('surface index has not been set. Make sure to call set_surface_index.')
      ready = .false.
   end if
#endif

   ! Host has sent all fields - disable all optional fields that were not provided in individual BGC models.
   call filter_readable_variable_registry(self)

   call self%job_manager%process_indices(unfulfilled_dependencies)

   call self%job_manager%print()

   variable_node => unfulfilled_dependencies%first
   do while (associated(variable_node))
      call report_unfulfilled_dependency(variable_node%target)
      variable_node => variable_node%next
   end do

   if (associated(unfulfilled_dependencies%first).or..not.ready) call fatal_error('fabm_check_ready','FABM is lacking required data.')

   self%state = state_check_ready_done

   contains

      subroutine report_unfulfilled_dependency(variable)
         type (type_internal_variable),target :: variable

         type type_model_reference
            class (type_base_model),     pointer :: p    => null()
            type (type_model_reference), pointer :: next => null()
         end type

         type (type_model_reference),pointer :: first,current,next
         type (type_link),           pointer :: link
         character(len=attribute_length)     :: path

         call log_message('UNFULFILLED DEPENDENCY: '//trim(variable%name))
         select case (variable%domain)
         case (domain_interior)
            call log_message('  This is an interior field.')
         case (domain_horizontal,domain_surface,domain_bottom)
            call log_message('  This is a horizontal-only field.')
         case (domain_scalar)
            call log_message('  This is a scalar field (single value valid across the entire domain).')
         end select
         if (variable%units/='') call log_message('  It has units '//trim(variable%units))
         call log_message('  It is needed by the following model instances:')

         first => null()
         link => self%root%links%first
         do while (associated(link))
            if (     associated(link%target,variable)     &                   ! This link points to the target variable,
                .and.associated(link%original%read_index) &                   ! the model that owns the link requests read access for it,
                .and.link%original%presence/=presence_external_optional) then ! and this access is required, not optional
               current => first
               do while (associated(current))
                  if (associated(current%p,link%original%owner)) exit
                  current => current%next
               end do
               if (.not.associated(current)) then
                  ! This model has not been reported before. Do so now and remember that we have done so.
                  allocate(current)
                  current%p => link%original%owner
                  current%next => first
                  first => current
                  path = current%p%get_path()
                  call log_message('    '//trim(path(2:)))
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

   end subroutine fabm_check_ready
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable identifier for the given variable
! name. Returns id\_not\_used if the variable name is unknown.
! The variable identifier can be used later in calls to
! fabm\_link\_data/fabm\_link\_data\_hz.
!
! !INTERFACE:
   function fabm_get_bulk_variable_id_by_name(self,name) result(id)
!
! !INPUT PARAMETERS:
   class (type_model),intent(in) :: self
   character(len=*),  intent(in) :: name
!
! !RETURN VALUE:
   type (type_bulk_variable_id) :: id
!
! !LOCAL VARIABLES:
   type (type_link), pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_interior) then
         if (link%name==name.or.get_safe_name(link%name)==name) then
            id%variable => link%target
            return
         end if
      end if
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_interior.and.link%target%standard_variables%contains(name)) then
         id%variable => link%target
         return
      end if
      link => link%next
   end do

   end function fabm_get_bulk_variable_id_by_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the variable identifier for specified standard
! variable, defined on the full model domain.
!
! !INTERFACE:
   function fabm_get_bulk_variable_id_sn(self,standard_variable) result(id)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(in) :: self
   type (type_bulk_standard_variable),intent(in) :: standard_variable
!
! !RETURN VALUE:
   type (type_bulk_variable_id) :: id
!
! !LOCAL VARIABLES:
   type (type_link), pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%standard_variables%contains(standard_variable)) then
         id%variable => link%target
         return
      end if
      link => link%next
   end do

   end function fabm_get_bulk_variable_id_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable identifier for the given variable
! name.
!
! !INTERFACE:
   function fabm_get_horizontal_variable_id_by_name(self,name) result(id)
!
! !INPUT PARAMETERS:
   class (type_model),intent(in) :: self
   character(len=*),  intent(in) :: name
!
! !RETURN VALUE:
   type (type_horizontal_variable_id) :: id
!
! !LOCAL VARIABLES:
   type (type_link), pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_horizontal.or.link%target%domain==domain_surface.or.link%target%domain==domain_bottom) then
         if (link%name==name.or.get_safe_name(link%name)==name) then
            id%variable => link%target
            return
         end if
      end if
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%links%first
   do while (associated(link))
      if ((link%target%domain==domain_horizontal.or.link%target%domain==domain_surface.or.link%target%domain==domain_bottom).and.link%target%standard_variables%contains(name)) then
         id%variable => link%target
         return
      end if
      link => link%next
   end do

   end function fabm_get_horizontal_variable_id_by_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the variable identifier for specified standard
! variable, defined on a horizontal slice of the model domain.
!
! !INTERFACE:
   function fabm_get_horizontal_variable_id_sn(self,standard_variable) result(id)
!
! !INPUT PARAMETERS:
   class (type_model),                      intent(in) :: self
   type (type_horizontal_standard_variable),intent(in) :: standard_variable
!
! !RETURN VALUE:
   type (type_horizontal_variable_id) :: id
!
! !LOCAL VARIABLES:
   type (type_link), pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%standard_variables%contains(standard_variable)) then
         id%variable => link%target
         return
      end if
      link => link%next
   end do

   end function fabm_get_horizontal_variable_id_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable identifier for the given variable
! name. Returns id\_not\_used if the variable name is unknown.
! The variable identifier can be used later in calls to
! fabm\_link\_data/fabm\_link\_data\_hz.
!
! !INTERFACE:
   function fabm_get_scalar_variable_id_by_name(self,name) result(id)
!
! !INPUT PARAMETERS:
   class (type_model),intent(in)  :: self
   character(len=*),  intent(in)  :: name
!
! !RETURN VALUE:
   type (type_scalar_variable_id) :: id
!
! !LOCAL VARIABLES:
   type (type_link), pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_scalar) then
         if (link%name==name.or.get_safe_name(link%name)==name) then
            id%variable => link%target
            return
         end if
      end if
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_scalar.and.link%target%standard_variables%contains(name)) then
         id%variable => link%target
         return
      end if
      link => link%next
   end do

   end function fabm_get_scalar_variable_id_by_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the variable identifier for specified space-
! independent standard variable.
!
! !INTERFACE:
   function fabm_get_scalar_variable_id_sn(self,standard_variable) result(id)
!
! !INPUT PARAMETERS:
   class (type_model),                  intent(in) :: self
   type (type_global_standard_variable),intent(in) :: standard_variable
!
! !RETURN VALUE:
   type (type_scalar_variable_id) :: id
!
! !LOCAL VARIABLES:
   type (type_link), pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%standard_variables%contains(standard_variable)) then
         id%variable => link%target
         return
      end if
      link => link%next
   end do

   end function fabm_get_scalar_variable_id_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable name for the given variable
! identifier.
!
! !INTERFACE:
   function fabm_get_variable_name(model,id) result(name)
!
! !INPUT PARAMETERS:
   class (type_model),              intent(in) :: model
   class(type_external_variable_id),intent(in) :: id
!
! !RETURN VALUE:
   character(len=attribute_length)             :: name
!
!EOP
!-----------------------------------------------------------------------
!BOC
   name = ''
   if (associated(id%variable)) name = get_safe_name(id%variable%name)

   end function fabm_get_variable_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Determine whether the specified interior variable is used [required]
! by biogeochemical models running in FABM. This does NOT imply its values need
! to be provided by the host; the values may be provided by a FABM module.
!
! !INTERFACE:
   function fabm_is_variable_used(id) result(used)
!
! !INPUT PARAMETERS:
   class(type_external_variable_id), intent(in) :: id
!
! !RETURN VALUE:
   logical                                      :: used
!
!EOP
!-----------------------------------------------------------------------
!BOC
   used = associated(id%variable)
   if (used) used = id%variable%read_indices%value/=-1

   end function fabm_is_variable_used
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Determine whether values for an interior variable are required.
! That is, one or more FABM modules read the variable's value, but no existing module
! provides the value. The variable is specified by its identifier.
!
! !INTERFACE:
   function fabm_interior_variable_needs_values(self,id) result(required)
!
! !INPUT PARAMETERS:
   class (type_model),         intent(in) :: self
   type(type_bulk_variable_id),intent(in) :: id
!
! !RETURN VALUE:
   logical :: required
!
!EOP
!-----------------------------------------------------------------------
!BOC
   required = associated(id%variable)
   if (required) required = id%variable%read_indices%value/=-1
   if (required) required = .not.associated(self%data(id%variable%read_indices%value)%p)

   end function fabm_interior_variable_needs_values
!EOC

!BOP
!
! !IROUTINE: Determine whether values for an interior variable are required.
! That is, one or more FABM modules read the variable's value, but no existing module
! provides the value. The variable is specified by its identity, a "standard
! variable" object.
!
! !INTERFACE:
   function fabm_interior_variable_needs_values_sn(self,standard_variable) result(required)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(in) :: self
   type (type_bulk_standard_variable),intent(in) :: standard_variable
!
! !RETURN VALUE:
   logical :: required
!
! !LOCAL VARIABLES:
   type(type_bulk_variable_id) :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   id = self%get_bulk_variable_id(standard_variable)
   required = self%variable_needs_values(id)

   end function fabm_interior_variable_needs_values_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Determine whether values for a horizontal variable are required.
! That is, one or more FABM modules read the variable's value, but no existing module
! provides the value. The variable is specified by its identifier
!
! !INTERFACE:
   function fabm_horizontal_variable_needs_values(self,id) result(required)
!
! !INPUT PARAMETERS:
   class (type_model),               intent(in) :: self
   type(type_horizontal_variable_id),intent(in) :: id
!
! !RETURN VALUE:
   logical :: required
!
!EOP
!-----------------------------------------------------------------------
!BOC
   required = associated(id%variable)
   if (required) required = id%variable%read_indices%value/=-1
   if (required) required = .not.associated(self%data_hz(id%variable%read_indices%value)%p)

   end function fabm_horizontal_variable_needs_values
!EOC

!BOP
!
! !IROUTINE: Determine whether values for a horizontal variable are required.
! That is, one or more FABM modules read the variable's value, but no existing module
! provides the value. The variable is specified by its identity, a "standard variable"
! object.
!
! !INTERFACE:
   function fabm_horizontal_variable_needs_values_sn(self,standard_variable) result(required)
!
! !INPUT PARAMETERS:
   class (type_model),                      intent(in) :: self
   type (type_horizontal_standard_variable),intent(in) :: standard_variable
!
! !RETURN VALUE:
   logical :: required
!
! !LOCAL VARIABLES:
   type(type_horizontal_variable_id) :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   id = self%get_horizontal_variable_id(standard_variable)
   required = self%variable_needs_values(id)

   end function fabm_horizontal_variable_needs_values_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Determine whether the value for a scalar variable is required.
! That is, one or more FABM modules read the variable's value, but no existing module
! provides the value. The variable is specified by its identifier.
!
! !INTERFACE:
   function fabm_scalar_variable_needs_values(self,id) result(required)
!
! !INPUT PARAMETERS:
   class (type_model),           intent(in) :: self
   type(type_scalar_variable_id),intent(in) :: id
!
! !RETURN VALUE:
   logical :: required
!
!EOP
!-----------------------------------------------------------------------
!BOC
   required = associated(id%variable)
   if (required) required = id%variable%read_indices%value/=-1
   if (required) required = .not.associated(self%data_scalar(id%variable%read_indices%value)%p)

   end function fabm_scalar_variable_needs_values
!EOC

!BOP
!
! !IROUTINE: Determine whether the value for a scalar variable is required.
! That is, one or more FABM modules read the variable's value, but no existing module
! provides the value. The variable is specified by its identity, a "standard variable"
! object.
!
! !INTERFACE:
   function fabm_scalar_variable_needs_values_sn(self,standard_variable) result(required)
!
! !INPUT PARAMETERS:
   class (type_model),                  intent(in) :: self
   type (type_global_standard_variable),intent(in) :: standard_variable
!
! !RETURN VALUE:
   logical :: required
!
! !LOCAL VARIABLES:
   type(type_scalar_variable_id) :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   id = self%get_scalar_variable_id(standard_variable)
   required = self%variable_needs_values(id)

   end function fabm_scalar_variable_needs_values_sn
!EOC

   function get_host_container_model(self) result(host)
      class (type_model), intent(inout) :: self
      class (type_host_container), pointer :: host

      type (type_model_list_node),pointer :: node
      class (type_base_model), pointer :: base_host

      node => self%root%children%find('_host_')
      if (associated(node)) then
         base_host => node%model
         select type (base_host)
         class is (type_host_container)
            host => base_host
         end select
      else
         allocate(host)
         call self%root%add_child(host,'_host_',configunit=-1)
      end if
   end function get_host_container_model

   subroutine fabm_require_interior_data(self,standard_variable,domain)
      class (type_model),                intent(inout) :: self
      type(type_bulk_standard_variable), intent(in)    :: standard_variable
      integer,optional,                  intent(in)    :: domain

      class (type_host_container),  pointer :: host
      type (type_integer_list_node),pointer :: node
      type (type_link),             pointer :: link
      integer                               :: domain_

      if (self%state>=state_initialize_done) &
         call fatal_error('fabm_require_interior_data','model%require_data cannot be called after model initialization.')

      domain_ = domain_interior
      if (present(domain)) domain_ = domain

      host => get_host_container_model(self)

      allocate(node)
      node%next => host%first
      host%first => node
      select case (domain_)
      case (domain_interior)
         call host%add_interior_variable(standard_variable%name,standard_variable%units,standard_variable%name,read_index=node%value,link=link)
      case (domain_horizontal,domain_surface,domain_bottom)
         call host%add_horizontal_variable(standard_variable%name,standard_variable%units,standard_variable%name,read_index=node%value,domain=domain_,link=link)
      case default
         call fatal_error('fabm_require_interior_data','model%require_data called with unknown domain.')
      end select
      call host%request_coupling(link,standard_variable,domain=domain_)
   end subroutine fabm_require_interior_data

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by an internal variable object.
!
! !INTERFACE:
   subroutine fabm_link_interior_data_by_variable(self,variable,dat,source)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: self
   type(type_internal_variable),      intent(in)    :: variable
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
   integer,optional,                  intent(in)    :: source
!
! !LOCAL VARIABLES:
   integer :: i
   integer :: source_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&_FABM_DIMENSION_COUNT_>0
   do i=1,size(self%domain_size)
      if (size(dat,i)/=self%domain_size(i)) then
         call fatal_error('fabm_link_interior_data_by_variable','dimensions of FABM domain and provided array do not match for variable '//trim(variable%name)//'.')
      end if
   end do
#endif

   i = variable%read_indices%value
   if (i/=-1) then
      source_ = data_source_default
      if (present(source)) source_ = source
      if (source_>=self%interior_data_sources(i)) then
         self%data(i)%p => dat
         self%interior_data_sources(i) = source_
      end if
   end if

   end subroutine fabm_link_interior_data_by_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by an external identifier.
!
! !INTERFACE:
   subroutine fabm_link_interior_data_by_id(self,id,dat,source)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: self
   type(type_bulk_variable_id),       intent(in)    :: id
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
   integer,optional,                  intent(in)    :: source
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (associated(id%variable)) call fabm_link_interior_data_by_variable(self,id%variable,dat,source)

   end subroutine fabm_link_interior_data_by_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the interior of the spatial domain. The variable
! is identified by its identity, a "standard variable" object.
!
! !INTERFACE:
   subroutine fabm_link_interior_data_by_sn(model,standard_variable,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: model
   type(type_bulk_standard_variable), intent(in)    :: standard_variable
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
!
! !LOCAL VARIABLES:
   type (type_bulk_variable_id)                             :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_bulk_variable_id_sn(model,standard_variable)

   ! Only link the data if needed (if the variable identifier is valid).
   if (fabm_is_variable_used(id)) call fabm_link_interior_data(model,id,dat)

   end subroutine fabm_link_interior_data_by_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the interior of the spatial domain.
! The variable is identified by its name.
!
! !INTERFACE:
   subroutine fabm_link_interior_data_by_name(model,name,dat)
!
! !INPUT PARAMETERS:
   class (type_model),target,         intent(inout) :: model
   character(len=*),                  intent(in)    :: name
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
!
! !LOCAL VARIABLES:
   type (type_bulk_variable_id) :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_bulk_variable_id(model,name)

   ! Only link the data if needed (if the variable identifier is valid).
   if (fabm_is_variable_used(id)) call fabm_link_interior_data(model,id,dat)

   end subroutine fabm_link_interior_data_by_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by an internal variable object.
!
! !INTERFACE:
   subroutine fabm_link_horizontal_data_by_variable(self,variable,dat,source)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: self
   type (type_internal_variable),                intent(in)    :: variable
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
   integer,optional,                             intent(in)    :: source
!
! !LOCAL VARIABLES:
   integer :: i
   integer :: source_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&_HORIZONTAL_DIMENSION_COUNT_>0
   do i=1,size(self%horizontal_domain_size)
      if (size(dat,i)/=self%horizontal_domain_size(i)) then
         call fatal_error('fabm_link_horizontal_data','dimensions of FABM domain and provided array do not match for variable '//trim(variable%name)//'.')
      end if
   end do
#endif

   i = variable%read_indices%value
   if (i/=-1) then
      source_ = data_source_default
      if (present(source)) source_ = source
      if (source_>=self%horizontal_data_sources(i)) then
         self%data_hz(i)%p => dat
         self%horizontal_data_sources(i) = source_
      end if
   end if

   end subroutine fabm_link_horizontal_data_by_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by an external identifier.
!
! !INTERFACE:
   subroutine fabm_link_horizontal_data_by_id(self,id,dat,source)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: self
   type (type_horizontal_variable_id),           intent(in)    :: id
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
   integer,optional,                             intent(in)    :: source
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (associated(id%variable)) call fabm_link_horizontal_data_by_variable(self,id%variable,dat,source)

   end subroutine fabm_link_horizontal_data_by_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by its identity, a "standard variable" object.
!
! !INTERFACE:
   subroutine fabm_link_horizontal_data_by_sn(model,standard_variable,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: model
   type(type_horizontal_standard_variable),      intent(in)    :: standard_variable
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
!
! !LOCAL VARIABLES:
   type (type_horizontal_variable_id)                             :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_horizontal_variable_id_sn(model,standard_variable)

   ! Only link the data if needed (if the variable identifier is valid).
   if (fabm_is_variable_used(id)) call fabm_link_horizontal_data(model,id,dat)

   end subroutine fabm_link_horizontal_data_by_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by its name.
!
! !INTERFACE:
   subroutine fabm_link_horizontal_data_by_name(model,name,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: model
   character(len=*),                             intent(in)    :: name
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
!
! !LOCAL VARIABLES:
   type(type_horizontal_variable_id) :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_horizontal_variable_id(model,name)

   ! Only link the data if needed (if the variable identifier is valid).
   if (fabm_is_variable_used(id)) call fabm_link_horizontal_data(model,id,dat)

   end subroutine fabm_link_horizontal_data_by_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the scalar that will hold
! data for the specified variable. The variable is identified by an
! external identifier.
!
! !INTERFACE:
   subroutine fabm_link_scalar_by_id(self,id,dat,source)
!
! !INPUT PARAMETERS:
   class (type_model),            intent(inout) :: self
   type (type_scalar_variable_id),intent(in)    :: id
   real(rk),target,               intent(in)    :: dat
   integer,optional,              intent(in)    :: source
!
! !LOCAL VARIABLES:
   integer :: i
   integer :: source_
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (.not.associated(id%variable)) return
   i = id%variable%read_indices%value
   if (i/=-1) then
      source_ = data_source_default
      if (present(source)) source_ = source
      if (source_>=self%scalar_data_sources(i)) then
         self%data_scalar(i)%p => dat
         self%scalar_data_sources(i) = source_
      end if
   end if

   end subroutine fabm_link_scalar_by_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by its identity, a "standard variable" object.
!
! !INTERFACE:
   subroutine fabm_link_scalar_by_sn(model,standard_variable,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                 intent(inout) :: model
   type(type_global_standard_variable),intent(in)    :: standard_variable
   real(rk),target,                    intent(in)    :: dat
!
! !LOCAL VARIABLES:
   type (type_scalar_variable_id) :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_scalar_variable_id_sn(model,standard_variable)

   ! Only link the data if needed (if the variable identifier is valid).
   if (fabm_is_variable_used(id)) call fabm_link_scalar_data(model,id,dat)

   end subroutine fabm_link_scalar_by_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the scalar that will hold
! data for the specified variable. The variable is identified by its name.
!
! !INTERFACE:
   subroutine fabm_link_scalar_by_name(model,name,dat)
!
! !INPUT PARAMETERS:
   class (type_model),intent(inout) :: model
   character(len=*),  intent(in)    :: name
   real(rk),target,   intent(in)    :: dat
!
! !LOCAL VARIABLES:
   type(type_scalar_variable_id) :: id
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_scalar_variable_id(model,name)

   ! Only link the data if needed (if the variable identifier is valid).
   if (fabm_is_variable_used(id)) call fabm_link_scalar_data(model,id,dat)

   end subroutine fabm_link_scalar_by_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! a single pelagic state variable.
!
! !INTERFACE:
   subroutine fabm_link_interior_state_data(self,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: self
   integer,                           intent(in)    :: id
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_link_interior_data(self,self%state_variables(id)%target,dat,source=data_source_fabm)

   end subroutine fabm_link_interior_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! a single benthic state variable.
!
! !INTERFACE:
   subroutine fabm_link_bottom_state_data(self,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: self
   integer,                                      intent(in)    :: id
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_link_horizontal_data(self,self%bottom_state_variables(id)%target,dat,source=data_source_fabm)

   end subroutine fabm_link_bottom_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! a single surface-bound state variable.
!
! !INTERFACE:
   subroutine fabm_link_surface_state_data(self,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: self
   integer,                                      intent(in)    :: id
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_link_horizontal_data(self,self%surface_state_variables(id)%target,dat,source=data_source_fabm)

   end subroutine fabm_link_surface_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with data for all pelagic state variables.
!
! !INTERFACE:
   subroutine fabm_link_all_interior_state_data(self,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                       intent(inout) :: self
   real(rk) _DIMENSION_GLOBAL_PLUS_1_,target,intent(in)    :: dat

   integer :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   if (size(dat,_FABM_DIMENSION_COUNT_+1)/=size(self%state_variables)) &
      call fatal_error('fabm_link_all_interior_state_data','length of last dimension of provided array must match number of interior state variables.')
#endif
   do i=1,size(self%state_variables)
      call fabm_link_interior_state_data(self,i,dat(_PREARG_LOCATION_DIMENSIONS_ i))
   end do

   end subroutine fabm_link_all_interior_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with data for all bottom state variables.
!
! !INTERFACE:
   subroutine fabm_link_all_bottom_state_data(self,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                                  intent(inout) :: self
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_,target,intent(in)    :: dat

   integer :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   if (size(dat,_HORIZONTAL_DIMENSION_COUNT_+1)/=size(self%bottom_state_variables)) &
      call fatal_error('fabm_link_all_bottom_state_data','length of last dimension of provided array must match number of bottom state variables.')
#endif
   do i=1,size(self%bottom_state_variables)
      call fabm_link_bottom_state_data(self,i,dat(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ i))
   end do

   end subroutine fabm_link_all_bottom_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with data for all surface state variables.
!
! !INTERFACE:
   subroutine fabm_link_all_surface_state_data(self,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                                  intent(inout) :: self
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_,target,intent(in)    :: dat

   integer :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   if (size(dat,_HORIZONTAL_DIMENSION_COUNT_+1)/=size(self%surface_state_variables)) &
      call fatal_error('fabm_link_all_surface_state_data','length of last dimension of provided array must match number of surface state variables.')
#endif
   do i=1,size(self%surface_state_variables)
      call fabm_link_surface_state_data(self,i,dat(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ i))
   end do

   end subroutine fabm_link_all_surface_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Returns (a pointer to) the array with data for
! a single diagnostic variable, defined on the full spatial domain.
!
! !INTERFACE:
   function fabm_get_interior_diagnostic_data(self,id) result(dat)
!
! !INPUT PARAMETERS:
   class (type_model),target,         intent(in) :: self
   integer,                           intent(in) :: id
   real(rk) _DIMENSION_GLOBAL_,pointer           :: dat
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   if (.not.allocated(self%diag)) call fatal_error('fabm_get_interior_diagnostic_data','Diagnostics have not been allocated yet.')
#endif

   ! Retrieve a pointer to the array holding the requested data.
   dat => self%diag(_PREARG_LOCATION_DIMENSIONS_ self%diagnostic_variables(id)%target%store_index)

   end function fabm_get_interior_diagnostic_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Returns (a pointer to) the array with data for
! a single diagnostic variable, defined on a horitontal slice of the
! spatial domain.
!
! !INTERFACE:
   function fabm_get_horizontal_diagnostic_data(self,id) result(dat)
!
! !INPUT PARAMETERS:
   class (type_model),target,                  intent(in) :: self
   integer,                                    intent(in) :: id
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,pointer         :: dat
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   if (.not.allocated(self%diag_hz)) call fatal_error('fabm_get_horizontal_diagnostic_data','Diagnostics have not been allocated yet.')
#endif

   ! Retrieve a pointer to the array holding the requested data.
   dat => self%diag_hz(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ self%horizontal_diagnostic_variables(id)%target%store_index)

   end function fabm_get_horizontal_diagnostic_data
!EOC

function fabm_get_interior_data(self,id) result(dat)
   class (type_model),target,         intent(in) :: self
   type(type_bulk_variable_id),       intent(in) :: id
   real(rk) _DIMENSION_GLOBAL_,pointer           :: dat

   nullify(dat)
   if (id%variable%read_indices%value/=-1) dat => self%data(id%variable%read_indices%value)%p
end function fabm_get_interior_data

function fabm_get_horizontal_data(self,id) result(dat)
   class (type_model),target,          intent(in) :: self
   type(type_horizontal_variable_id),  intent(in) :: id
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,pointer :: dat

   nullify(dat)
   if (id%variable%read_indices%value/=-1) dat => self%data_hz(id%variable%read_indices%value)%p
end function fabm_get_horizontal_data

function fabm_get_scalar_data(self,id) result(dat)
   class (type_model),target,          intent(in) :: self
   type(type_scalar_variable_id),      intent(in) :: id
   real(rk),pointer                               :: dat

   nullify(dat)
   if (id%variable%read_indices%value/=-1) dat => self%data_scalar(id%variable%read_indices%value)%p
end function fabm_get_scalar_data

subroutine begin_interior_task(self,task,cache _ARGUMENTS_INTERIOR_IN_)
   type (type_model),intent(inout) :: self
   type (type_task), intent(in)    :: task
   type (type_cache),intent(out)   :: cache
   _DECLARE_ARGUMENTS_INTERIOR_IN_
   _DECLARE_INTERIOR_INDICES_

   integer :: i

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
   _N_ = _STOP_-_START_+1
#  ifdef _HAS_MASK_
   allocate(cache%mask(_N_))
   _DO_CONCURRENT_(_I_,1,_N_)
#    ifdef _FABM_HORIZONTAL_MASK_
      cache%mask _INDEX_SLICE_ = _IS_UNMASKED_(self%mask_hz _INDEX_GLOBAL_HORIZONTAL_(_START_+_I_-1))
#    else
      cache%mask _INDEX_SLICE_ = _IS_UNMASKED_(self%mask _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1))
#    endif
   end do
   _N_ = count(cache%mask)
#  endif
#endif
   
#if defined(_INTERIOR_IS_VECTORIZED_)
   allocate(cache%read(_N_,size(self%data)))
#else
   allocate(cache%read(size(self%data)))
#endif
   do i=1,size(task%load)
      if (task%load(i)) then
         _PACK_GLOBAL_(self%data(i)%p,cache%read,i,cache%mask)
      end if
   end do

#ifdef _HORIZONTAL_IS_VECTORIZED_
   allocate(cache%read_hz(_N_,size(self%data_hz)))
#else
   allocate(cache%read_hz(size(self%data_hz)))
#endif
   do i=1,size(task%load_hz)
      if (task%load_hz(i)) then
         _HORIZONTAL_PACK_GLOBAL_(self%data_hz(i)%p,cache%read_hz,i,cache%mask)
      end if
   end do

   ! Copy global scalars to read cache
   allocate(cache%read_scalar(size(self%data_scalar)))
   do i=1,size(task%load_scalar)
      if (task%load_scalar(i)) cache%read_scalar(i) = self%data_scalar(i)%p
   end do

#if defined(_INTERIOR_IS_VECTORIZED_)
   allocate(cache%write(_N_,self%nscratch))
#else
   allocate(cache%write(self%nscratch))
#endif
   do i=1,size(task%prefill_type)
      if (task%prefill_type(i)==prefill_constant) then
         _CONCURRENT_LOOP_BEGIN_
            cache%write _INDEX_SLICE_PLUS_1_(i) = task%prefill_values(i)
         _LOOP_END_
      elseif (task%prefill_type(i)==prefill_previous_value) then
         _PACK_GLOBAL_PLUS_1_(self%diag,task%prefill_index(i),cache%write,i,cache%mask)
      end if
   end do
end subroutine begin_interior_task

subroutine begin_horizontal_task(self,task,cache _ARGUMENTS_HORIZONTAL_IN_)
   type (type_model),intent(inout) :: self
   type (type_task), intent(in)    :: task
   type (type_cache),intent(out)   :: cache
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   _DECLARE_HORIZONTAL_INDICES_

   integer :: i

#ifdef _HORIZONTAL_IS_VECTORIZED_
   _N_ = _STOP_-_START_+1
#  ifdef _HAS_MASK_
   allocate(cache%mask(_N_))
   _DO_CONCURRENT_(_J_,1,_N_)
      cache%mask _INDEX_HORIZONTAL_SLICE_ = _IS_UNMASKED_(self%mask_hz _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1))
   end do
   _N_ = count(cache%mask)
#  endif
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
   allocate(cache%read_hz(_N_,size(self%data_hz)))
#else
   allocate(cache%read_hz(size(self%data_hz)))
#endif
   do i=1,size(task%load_hz)
      if (task%load_hz(i)) then
         _HORIZONTAL_PACK_GLOBAL_(self%data_hz(i)%p,cache%read_hz,i,cache%mask)
      end if
   end do

   ! Copy global scalars to read cache
   allocate(cache%read_scalar(size(self%data_scalar)))
   do i=1,size(task%load_scalar)
      if (task%load_scalar(i)) cache%read_scalar(i) = self%data_scalar(i)%p
   end do

#ifdef _HORIZONTAL_IS_VECTORIZED_
   allocate(cache%write_hz(_N_,self%nscratch_hz))
#else
   allocate(cache%write_hz(self%nscratch_hz))
#endif
   do i=1,size(task%prefill_type_hz)
      if (task%prefill_type_hz(i)==prefill_constant) then
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = task%prefill_values_hz(i)
         _HORIZONTAL_LOOP_END_
      elseif (task%prefill_type_hz(i)==prefill_previous_value) then
         _HORIZONTAL_PACK_GLOBAL_PLUS_1_(self%diag_hz,task%prefill_index_hz(i),cache%write_hz,i,cache%mask)
      end if
   end do

   ! Also load boundary values for interior fields if performing surface or bottom-specific operations.
   if (task%operation==source_do_surface) then
      call load_surface_data(self,task,cache _ARGUMENTS_HORIZONTAL_IN_)
   elseif (task%operation==source_do_bottom) then
      call load_bottom_data (self,task,cache _ARGUMENTS_HORIZONTAL_IN_)
   end if
end subroutine begin_horizontal_task

subroutine load_surface_data(self,task,cache _ARGUMENTS_HORIZONTAL_IN_)
   type (type_model),intent(inout) :: self
   type (type_task), intent(in)    :: task
   type (type_cache),intent(inout) :: cache
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   _DECLARE_HORIZONTAL_INDICES_

   integer :: i

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
   integer :: _VERTICAL_ITERATOR_
   _VERTICAL_ITERATOR_ = self%surface_index
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
   allocate(cache%read(_N_,size(self%data)))
#elif defined(_INTERIOR_IS_VECTORIZED_)
   allocate(cache%read(1,size(self%data)))
#else
   allocate(cache%read(size(self%data)))
#endif
   do i=1,size(task%load)
      if (task%load(i)) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
#  ifdef _HAS_MASK_
         cache%read(:,i) = pack(self%data(i)%p _INDEX_GLOBAL_INTERIOR_(_START_:_STOP_),cache%mask)
#  else
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            cache%read _INDEX_SLICE_PLUS_1_(i) = self%data(i)%p _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1)
         _HORIZONTAL_LOOP_END_
#  endif
#elif defined(_INTERIOR_IS_VECTORIZED_)
         cache%read(1,i) = self%data(i)%p _INDEX_LOCATION_
#else
         cache%read(i) = self%data(i)%p _INDEX_LOCATION_
#endif
      end if
   end do
end subroutine load_surface_data
   
subroutine load_bottom_data(self,task,cache _ARGUMENTS_HORIZONTAL_IN_)
   type (type_model),intent(inout) :: self
   type (type_task), intent(in)    :: task
   type (type_cache),intent(inout) :: cache
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   _DECLARE_HORIZONTAL_INDICES_

   integer :: i
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
   integer :: _VERTICAL_ITERATOR_
#endif
#if _FABM_BOTTOM_INDEX_==-1
   integer :: j
#endif

#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_BOTTOM_INDEX_==0
   _VERTICAL_ITERATOR_ = self%bottom_index
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
   allocate(cache%read(_N_,size(self%data)))
#elif defined(_INTERIOR_IS_VECTORIZED_)
   allocate(cache%read(1,size(self%data)))
#else
   allocate(cache%read(size(self%data)))
#endif
   do i=1,size(task%load)
      if (task%load(i)) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
#  ifdef _HAS_MASK_
#    if _FABM_BOTTOM_INDEX_==-1
         j = 0
         do _J_=1,_STOP_-_START_+1
            if (cache%mask(_J_)) then
               _VERTICAL_ITERATOR_ = self%bottom_indices _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1)
               j = j + 1
               cache%read(j,i) = self%data(i)%p _INDEX_GLOBAL_INTERIOR_(_START_+_J_-1)
            end if
         end do
#    else
         cache%read(:,i) = pack(self%data(i)%p _INDEX_GLOBAL_INTERIOR_(_START_:_STOP_),cache%mask)
#    endif
#  else
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
#    if _FABM_BOTTOM_INDEX_==-1
            _VERTICAL_ITERATOR_ = self%bottom_indices _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1)
#    endif
            cache%read _INDEX_SLICE_PLUS_1_(i) = self%data(i)%p _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1)
         _HORIZONTAL_LOOP_END_
#  endif
#else
#  if _FABM_BOTTOM_INDEX_==-1
         _VERTICAL_ITERATOR_ = self%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#  endif
#  if defined(_INTERIOR_IS_VECTORIZED_)
         cache%read(1,i) = self%data(i)%p _INDEX_LOCATION_
#  else
         cache%read(i) = self%data(i)%p _INDEX_LOCATION_
#  endif
#endif
      end if
   end do
end subroutine load_bottom_data

subroutine begin_vertical_task(self,task,cache _ARGUMENTS_VERTICAL_IN_)
   type (type_model),intent(inout) :: self
   type (type_task), intent(in)    :: task
   type (type_cache),intent(out)   :: cache
   _DECLARE_ARGUMENTS_VERTICAL_IN_
   _DECLARE_VERTICAL_INDICES_

   integer :: i

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
   _N_ = _VERTICAL_STOP_-_VERTICAL_START_+1
#  ifdef _HAS_MASK_
   allocate(cache%mask(_N_))
   _DO_CONCURRENT_(_I_,1,_N_)
#    ifdef _FABM_HORIZONTAL_MASK_
      cache%mask _INDEX_SLICE_ = _IS_UNMASKED_(self%mask_hz _INDEX_HORIZONTAL_LOCATION_)
#    else
      cache%mask _INDEX_SLICE_ = _IS_UNMASKED_(self%mask _INDEX_GLOBAL_VERTICAL_(_VERTICAL_START_+_I_-1))
#    endif
   end do
   _N_ = count(cache%mask)
#  endif
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
   allocate(cache%read(_N_,size(self%data)))
#elif defined(_INTERIOR_IS_VECTORIZED_)
   allocate(cache%read(1,size(self%data)))
#else
   allocate(cache%read(size(self%data)))
#endif
   do i=1,size(task%load)
      if (task%load(i)) then
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  ifdef _HAS_MASK_
         cache%read(:,i) = pack(self%data(i)%p _INDEX_GLOBAL_VERTICAL_(_VERTICAL_START_:_VERTICAL_STOP_),cache%mask)
#  else
         _CONCURRENT_VERTICAL_LOOP_BEGIN_
            cache%read _INDEX_SLICE_PLUS_1_(i) = self%data(i)%p _INDEX_GLOBAL_VERTICAL_(_VERTICAL_START_+_I_-1)
         _VERTICAL_LOOP_END_
#  endif
#elif defined(_INTERIOR_IS_VECTORIZED_)
         cache%read(1,i) = self%data(i)%p _INDEX_LOCATION_
#else
         cache%read(i) = self%data(i)%p _INDEX_LOCATION_
#endif
      end if
   end do

#ifdef _HORIZONTAL_IS_VECTORIZED_
   allocate(cache%read_hz(1,size(self%data_hz)))
#else
   allocate(cache%read_hz(size(self%data_hz)))
#endif
   do i=1,size(task%load_hz)
      if (task%load_hz(i)) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
         cache%read_hz(1,i) = self%data_hz(i)%p _INDEX_HORIZONTAL_LOCATION_
#else
         cache%read_hz(i) = self%data_hz(i)%p _INDEX_HORIZONTAL_LOCATION_
#endif
      end if
   end do

   ! Copy global scalars to read cache
   allocate(cache%read_scalar(size(self%data_scalar)))
   do i=1,size(task%load_scalar)
      if (task%load_scalar(i)) cache%read_scalar(i) = self%data_scalar(i)%p
   end do

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
   allocate(cache%write(_N_,self%nscratch))
#elif defined(_INTERIOR_IS_VECTORIZED_)
   allocate(cache%write(1,self%nscratch))
#else
   allocate(cache%write(self%nscratch))
#endif
   do i=1,size(task%prefill_type)
      if (task%prefill_type(i)==prefill_constant) then
#if defined(_INTERIOR_IS_VECTORIZED_)
         cache%write(:,i) = task%prefill_values(i)
#else
         cache%write(i) = task%prefill_values(i)
#endif
      end if
   end do

#ifdef _HORIZONTAL_IS_VECTORIZED_
   allocate(cache%write_hz(1,self%nscratch_hz))
#else
   allocate(cache%write_hz(self%nscratch_hz))
#endif
end subroutine begin_vertical_task

subroutine end_interior_task(self,task,cache _ARGUMENTS_INTERIOR_IN_)
   type (type_model),intent(inout) :: self
   type (type_task), intent(in)    :: task
   type (type_cache),intent(inout) :: cache
   _DECLARE_ARGUMENTS_INTERIOR_IN_
   _DECLARE_INTERIOR_INDICES_

   integer :: i

   ! Copy newly written diagnostics that need to be saved to global store.
   do i=1,size(task%save_sources)
      if (task%save_sources(i)/=-1) then
         _UNPACK_TO_GLOBAL_PLUS_1_(cache%write,task%save_sources(i),self%diag,i,cache%mask,self%diag_missing_value(i))
      end if
   end do

#ifdef _HAS_MASK_
   deallocate(cache%mask)
#endif
   deallocate(cache%read)
   deallocate(cache%read_hz)
   deallocate(cache%read_scalar)
   deallocate(cache%write)
end subroutine end_interior_task

subroutine end_horizontal_task(self,task,cache _ARGUMENTS_HORIZONTAL_IN_)
   type (type_model),intent(inout) :: self
   type (type_task), intent(in)    :: task
   type (type_cache),intent(inout) :: cache
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   _DECLARE_HORIZONTAL_INDICES_

   integer :: i

   ! Copy newly written horizontal diagnostics that need to be saved to global store.
   do i=1,size(task%save_sources_hz)
      if (task%save_sources_hz(i)/=-1) then
         _HORIZONTAL_UNPACK_TO_GLOBAL_PLUS_1_(cache%write_hz,task%save_sources_hz(i),self%diag_hz,i,cache%mask,self%diag_hz_missing_value(i))
      end if
   end do

#ifdef _HAS_MASK_
   deallocate(cache%mask)
#endif
   if (allocated(cache%read)) deallocate(cache%read)
   deallocate(cache%read_hz)
   deallocate(cache%read_scalar)
   deallocate(cache%write_hz)
end subroutine end_horizontal_task

subroutine end_vertical_task(self,task,cache _ARGUMENTS_VERTICAL_IN_)
   type (type_model),intent(inout) :: self
   type (type_task), intent(in)    :: task
   type (type_cache),intent(inout) :: cache
   _DECLARE_ARGUMENTS_VERTICAL_IN_
   _DECLARE_VERTICAL_INDICES_

   integer :: i

   ! Copy diagnostics that need to be saved to global store.
   do i=1,size(task%save_sources)
      if (task%save_sources(i)/=-1) then
         _VERTICAL_UNPACK_TO_GLOBAL_PLUS_1_(cache%write,task%save_sources(i),self%diag,i,cache%mask,self%diag_missing_value(i))
      end if
   end do
   do i=1,size(task%save_sources_hz)
      if (task%save_sources_hz(i)/=-1) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
         self%diag_hz(_PREARG_HORIZONTAL_LOCATION_ i) = cache%write_hz(1,task%save_sources_hz(i))
#else
         self%diag_hz(_PREARG_HORIZONTAL_LOCATION_ i) = cache%write_hz(task%save_sources_hz(i))
#endif
      end if
   end do

#ifdef _HAS_MASK_
   deallocate(cache%mask)
#endif
   deallocate(cache%read)
   deallocate(cache%read_hz)
   deallocate(cache%read_scalar)
   deallocate(cache%write)
   deallocate(cache%write_hz)
end subroutine end_vertical_task

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the model state (pelagic).
!
! !INTERFACE:
   subroutine fabm_initialize_state(self _ARGUMENTS_INTERIOR_IN_)
!
! !INPUT PARAMETERS:
   class (type_model),intent(inout) :: self
   _DECLARE_ARGUMENTS_INTERIOR_IN_
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   integer                   :: ivar,read_index
   type (type_call), pointer :: call_node
   logical                   :: set_interior
   _DECLARE_INTERIOR_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_interior_location(self _ARGUMENTS_INTERIOR_IN_,'fabm_initialize_state')
#endif

   call begin_interior_task(self,self%initialize_state_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   do ivar=1,size(self%state_variables)
      read_index = self%state_variables(ivar)%target%read_indices%value
      _CONCURRENT_LOOP_BEGIN_
         cache%read _INDEX_SLICE_PLUS_1_(read_index) = self%state_variables(ivar)%initial_value
      _LOOP_END_
   end do

   ! Allow biogeochemical models to initialize their interior state.
   set_interior = .false.
   call_node => self%initialize_state_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%initialize_state(_ARGUMENTS_INTERIOR_,set_interior)
      call_node => call_node%next
   end do

   ! Copy from cache back to global data store [NB variable values have been set in the *read* cache].
   do ivar=1,size(self%state_variables)
      read_index = self%state_variables(ivar)%target%read_indices%value
      if (self%interior_data_sources(read_index)==data_source_fabm) then
         _UNPACK_TO_GLOBAL_(cache%read,read_index,self%data(read_index)%p,cache%mask,self%state_variables(ivar)%missing_value)
      end if
   end do

   call end_interior_task(self,self%initialize_state_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_initialize_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the bottom model state.
!
! !INTERFACE:
   subroutine fabm_initialize_bottom_state(self _ARGUMENTS_HORIZONTAL_IN_)
!
! !INPUT PARAMETERS:
   class (type_model), intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   integer                   :: ivar,read_index
   type (type_call), pointer :: call_node
   logical                   :: set_horizontal
   _DECLARE_HORIZONTAL_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_initialize_bottom_state')
#endif

   call begin_horizontal_task(self,self%initialize_bottom_state_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   ! Initialize bottom variables
   do ivar=1,size(self%bottom_state_variables)
      read_index = self%bottom_state_variables(ivar)%target%read_indices%value
      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
         cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index) = self%bottom_state_variables(ivar)%initial_value
      _HORIZONTAL_LOOP_END_
   end do

   ! Allow biogeochemical models to initialize their bottom state.
   set_horizontal = .false.
   call_node => self%initialize_bottom_state_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%initialize_bottom_state(_ARGUMENTS_HORIZONTAL_,set_horizontal)
      call_node => call_node%next
   end do

   ! Copy from cache back to global data store [NB variable values have been set in the *read* cache].
   do ivar=1,size(self%bottom_state_variables)
      read_index = self%bottom_state_variables(ivar)%target%read_indices%value
      if (self%horizontal_data_sources(read_index)==data_source_fabm) then
         _HORIZONTAL_UNPACK_TO_GLOBAL_(cache%read_hz,read_index,self%data_hz(read_index)%p,cache%mask,self%bottom_state_variables(ivar)%missing_value)
      end if
   end do

   call end_horizontal_task(self,self%initialize_bottom_state_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_initialize_bottom_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the surface model state.
!
! !INTERFACE:
   subroutine fabm_initialize_surface_state(self _ARGUMENTS_HORIZONTAL_IN_)
!
! !INPUT PARAMETERS:
   class (type_model), intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   integer                   :: ivar,read_index
   type (type_call), pointer :: call_node
   logical                   :: set_horizontal
   _DECLARE_HORIZONTAL_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_initialize_surface_state')
#endif

   call begin_horizontal_task(self,self%initialize_surface_state_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   ! Initialize surface variables
   do ivar=1,size(self%surface_state_variables)
      read_index = self%surface_state_variables(ivar)%target%read_indices%value
      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
         cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index) = self%surface_state_variables(ivar)%initial_value
      _HORIZONTAL_LOOP_END_
   end do

   ! Allow biogeochemical models to initialize their surface state.
   set_horizontal = .false.
   call_node => self%initialize_surface_state_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%initialize_surface_state(_ARGUMENTS_HORIZONTAL_,set_horizontal)
      call_node => call_node%next
   end do

   ! Copy from cache back to global data store [NB variable values have been set in the *read* cache].
   do ivar=1,size(self%surface_state_variables)
      read_index = self%surface_state_variables(ivar)%target%read_indices%value
      if (self%horizontal_data_sources(read_index)==data_source_fabm) then
         _HORIZONTAL_UNPACK_TO_GLOBAL_(cache%read_hz,read_index,self%data_hz(read_index)%p,cache%mask,self%surface_state_variables(ivar)%missing_value)
      end if
   end do

   call end_horizontal_task(self,self%initialize_surface_state_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_initialize_surface_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the biogeochemical
! model tree.
!
! !INTERFACE:
   subroutine fabm_do_rhs(self _ARGUMENTS_INTERIOR_IN_,dy)
!
! !INPUT PARAMETERS:
   class (type_model),                    intent(inout) :: self
   _DECLARE_ARGUMENTS_INTERIOR_IN_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_EXT_SLICE_PLUS_1_, intent(inout) :: dy
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: i,j,k
   _DECLARE_INTERIOR_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_interior_location(self _ARGUMENTS_INTERIOR_IN_,'fabm_do_rhs')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
   call check_extents_2d(dy,_STOP_-_START_+1,size(self%state_variables),'fabm_do_rhs','dy','stop-start+1, # interior state variables')
#  else
   call check_extents_1d(dy,size(self%state_variables),'fabm_do_rhs','dy','# interior state variables')
#  endif
#endif

   call begin_interior_task(self,self%do_interior_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   call_node => self%do_interior_job%final_task%first_call
   do while (associated(call_node))
      if (call_node%source==source_do) then
         call call_node%model%do(_ARGUMENTS_INTERIOR_)
      elseif (call_node%source==source_get_light_extinction) then
         call call_node%model%get_light_extinction(_ARGUMENTS_INTERIOR_)
      end if

      ! Copy outputs of interest to read cache so consecutive models can use it.
      _DO_CONCURRENT_(i,1,size(call_node%copy_commands_int))
         j = call_node%copy_commands_int(i)%read_index
         k = call_node%copy_commands_int(i)%write_index
         _CONCURRENT_LOOP_BEGIN_
            cache%read _INDEX_SLICE_PLUS_1_(j) = cache%write _INDEX_SLICE_PLUS_1_(k)
         _LOOP_END_
      end do

      ! Move to next model
      call_node => call_node%next
   end do

   ! Compose total sources-sinks for each state variable, combining model-specific contributions.
   do i=1,size(self%state_variables)
      k = self%state_variables(i)%sms_index
      _UNPACK_AND_ADD_TO_PLUS_1_(cache%write,k,dy,i,cache%mask,0.0_rk)
   end do

   call end_interior_task(self,self%do_interior_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_do_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the biogeochemical
! model tree in the form of production and destruction matrices.
!
! !INTERFACE:
   subroutine fabm_do_ppdd(self _ARGUMENTS_INTERIOR_IN_,pp,dd)
!
! !INPUT PARAMETERS:
   class (type_model),                   intent(inout) :: self
  _DECLARE_ARGUMENTS_INTERIOR_IN_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_EXT_SLICE_PLUS_2_,intent(inout) :: pp,dd
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: i,j,k
   _DECLARE_INTERIOR_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_interior_location(self _ARGUMENTS_INTERIOR_IN_,'fabm_do_ppdd')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
   call check_extents_3d(pp,_STOP_-_START_+1,size(self%state_variables),size(self%state_variables),'fabm_do_ppdd','pp','stop-start+1, # interior state variables, # interior state variables')
   call check_extents_3d(dd,_STOP_-_START_+1,size(self%state_variables),size(self%state_variables),'fabm_do_ppdd','dd','stop-start+1, # interior state variables, # interior state variables')
#  else
   call check_extents_2d(pp,size(self%state_variables),size(self%state_variables),'fabm_do_ppdd','pp','# interior state variables, # interior state variables')
   call check_extents_2d(dd,size(self%state_variables),size(self%state_variables),'fabm_do_ppdd','dd','# interior state variables, # interior state variables')
#  endif
#endif

   call begin_interior_task(self,self%do_interior_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   call_node => self%do_interior_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%do_ppdd(_ARGUMENTS_INTERIOR_,pp,dd)

      ! Copy outputs of interest to read cache so consecutive models can use it.
      _DO_CONCURRENT_(i,1,size(call_node%copy_commands_int))
         j = call_node%copy_commands_int(i)%read_index
         k = call_node%copy_commands_int(i)%write_index
         _CONCURRENT_LOOP_BEGIN_
            cache%read _INDEX_SLICE_PLUS_1_(j) = cache%write _INDEX_SLICE_PLUS_1_(k)
         _LOOP_END_
      end do

      call_node => call_node%next
   end do

   call end_interior_task(self,self%do_interior_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks whether the current state is valid, and repairs [clips]
! invalid state variables if requested and possible.
!
! !INTERFACE:
   subroutine fabm_check_state(self _ARGUMENTS_INTERIOR_IN_,repair,valid)
!
! !INPUT PARAMETERS:
   class (type_model),     intent(inout) :: self
   _DECLARE_ARGUMENTS_INTERIOR_IN_
   logical,                intent(in)    :: repair
   logical,                intent(out)   :: valid
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   integer                   :: ivar, read_index
   type (type_call), pointer :: call_node
   real(rk)                  :: value,minimum,maximum
   character(len=256)        :: err
   logical                   :: set_interior
   _DECLARE_INTERIOR_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_interior_location(self _ARGUMENTS_INTERIOR_IN_,'fabm_check_state')
#endif

   valid = .true.
   set_interior = .false.

   call begin_interior_task(self,self%check_state_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   ! Allow individual models to check their state for their custom constraints, and to perform custom repairs.
   call_node => self%check_state_job%final_task%first_call
   do while (associated(call_node) .and. valid)
      call call_node%model%check_state(_ARGUMENTS_INTERIOR_,repair,valid,set_interior)
      if (.not. (valid .or. repair)) return
      call_node => call_node%next
   end do

   ! Finally check whether all state variable values lie within their prescribed [constant] bounds.
   ! This is always done, independently of any model-specific checks that may have been called above.

   ! Quick bounds check for the common case where all values are valid.
   do ivar=1,size(self%state_variables)
      read_index = self%state_variables(ivar)%target%read_indices%value
      minimum = self%state_variables(ivar)%minimum
      maximum = self%state_variables(ivar)%maximum
      _LOOP_BEGIN_
         value = cache%read _INDEX_SLICE_PLUS_1_(read_index)
         if (value<minimum.or.value>maximum) valid = .false.
      _LOOP_END_
   end do

   if (.not.valid) then

   ! Check boundaries for pelagic state variables specified by the models.
   ! If repair is permitted, this clips invalid values to the closest boundary.
   do ivar=1,size(self%state_variables)
      ! Shortcuts to variable information - this demonstrably helps the compiler (ifort).
      read_index = self%state_variables(ivar)%target%read_indices%value
      minimum = self%state_variables(ivar)%minimum
      maximum = self%state_variables(ivar)%maximum

      if (repair) then
         _CONCURRENT_LOOP_BEGIN_
            value = cache%read _INDEX_SLICE_PLUS_1_(read_index)
            cache%read _INDEX_SLICE_PLUS_1_(read_index) = max(minimum,min(maximum,value))
         _LOOP_END_
      else
         _LOOP_BEGIN_
            value = cache%read _INDEX_SLICE_PLUS_1_(read_index)
            if (value<minimum) then
               ! State variable value lies below prescribed minimum.
               write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ',trim(self%state_variables(ivar)%name), &
                                                          & ' below minimum value ',minimum
               call log_message(err)
               return
            elseif (value>maximum) then
               ! State variable value exceeds prescribed maximum.
               write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ',trim(self%state_variables(ivar)%name), &
                                                          & ' above maximum value ',maximum
               call log_message(err)
               return
            end if
         _LOOP_END_
      end if
   end do

   end if

   if (set_interior.or..not.valid) then
      do ivar=1,size(self%state_variables)
         read_index = self%state_variables(ivar)%target%read_indices%value
         if (self%interior_data_sources(read_index)==data_source_fabm) then
            _UNPACK_TO_GLOBAL_(cache%read,read_index,self%data(read_index)%p,cache%mask,self%state_variables(ivar)%missing_value)
         end if
      end do
   end if

   call end_interior_task(self,self%check_state_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_check_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks whether the current bottom state is valid, and repairs [clips]
! invalid state variables if requested and possible.
!
! !INTERFACE:
   subroutine fabm_check_bottom_state(self _ARGUMENTS_HORIZONTAL_IN_,repair,valid)
!
! !INPUT PARAMETERS:
   class (type_model),     intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   logical,                intent(in)    :: repair
   logical,                intent(out)   :: valid
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_check_bottom_state')
#endif

   call internal_check_horizontal_state(self,self%check_bottom_state_job _ARGUMENTS_HORIZONTAL_IN_,2,self%bottom_state_variables,repair,valid)

   end subroutine fabm_check_bottom_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks whether the current bottom state is valid, and repairs [clips]
! invalid state variables if requested and possible.
!
! !INTERFACE:
   subroutine fabm_check_surface_state(self _ARGUMENTS_HORIZONTAL_IN_,repair,valid)
!
! !INPUT PARAMETERS:
   class (type_model),     intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   logical,                intent(in)    :: repair
   logical,                intent(out)   :: valid
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_check_surface_state')
#endif

   call internal_check_horizontal_state(self,self%check_surface_state_job _ARGUMENTS_HORIZONTAL_IN_,1,self%info%surface_state_variables,repair,valid)

   end subroutine fabm_check_surface_state
!EOC

subroutine internal_check_horizontal_state(self,job _ARGUMENTS_HORIZONTAL_IN_,flag,state_variables,repair,valid)
   class (type_model),                        intent(inout) :: self
   type (type_job),                           intent(in)    :: job
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   integer,                                   intent(in)    :: flag
   type (type_horizontal_state_variable_info),intent(inout) :: state_variables(:)
   logical,                                   intent(in)    :: repair
   logical,                                   intent(out)   :: valid

   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: ivar,read_index
   real(rk)                  :: value,minimum,maximum
   character(len=256)        :: err
   _DECLARE_HORIZONTAL_INDICES_
   logical                   :: set_horizontal,set_interior
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
   integer :: _VERTICAL_ITERATOR_
#endif
#if _FABM_BOTTOM_INDEX_==-1&&defined(_HORIZONTAL_IS_VECTORIZED_)&&defined(_HAS_MASK_)
   integer :: j
#endif

   call begin_horizontal_task(self,job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   valid = .true.
   set_horizontal = .false.
   set_interior = .false.

   ! Allow individual models to check their state for their custom constraints, and to perform custom repairs.
   call_node => job%final_task%first_call
   do while (associated(call_node) .and. valid)
      if (flag==1) then
         call call_node%model%check_surface_state(_ARGUMENTS_HORIZONTAL_,repair,valid,set_horizontal,set_interior)
      else
         call call_node%model%check_bottom_state(_ARGUMENTS_HORIZONTAL_,repair,valid,set_horizontal,set_interior)
      end if
      if (.not. (valid .or. repair)) return
      call_node => call_node%next
   end do

   ! Check boundaries for horizontal state variables, as prescribed by the owning models.
   ! If repair is permitted, this clips invalid values to the closest boundary.
   do ivar=1,size(state_variables)
      ! Shortcuts to variable information - this demonstrably helps the compiler (ifort).
      read_index = state_variables(ivar)%target%read_indices%value
      minimum = state_variables(ivar)%minimum
      maximum = state_variables(ivar)%maximum

      _HORIZONTAL_LOOP_BEGIN_
         value = cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index)
         if (value<minimum) then
            ! State variable value lies below prescribed minimum.
            valid = .false.
            if (.not.repair) then
               write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ', &
                                                          & trim(state_variables(ivar)%name), &
                                                          & ' below minimum value ',minimum
               call log_message(err)
               return
            end if
            cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index) = minimum
         elseif (value>maximum) then
            ! State variable value exceeds prescribed maximum.
            valid = .false.
            if (.not.repair) then
               write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ', &
                                                          & trim(state_variables(ivar)%name), &
                                                          & ' above maximum value ',maximum
               call log_message(err)
               return
            end if
            cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(read_index) = maximum
         end if
      _HORIZONTAL_LOOP_END_
   end do

   if (set_horizontal.or..not.valid) then
      do ivar=1,size(state_variables)
         read_index = state_variables(ivar)%target%read_indices%value
         if (self%horizontal_data_sources(read_index)==data_source_fabm) then
            _HORIZONTAL_UNPACK_TO_GLOBAL_(cache%read_hz,read_index,self%data_hz(read_index)%p,cache%mask,state_variables(ivar)%missing_value)
         end if
      end do
   end if

   if (set_interior) then
      ! One or more models have provided new values for an interior state variable [at the interface]

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      if (flag==1) then
         _VERTICAL_ITERATOR_ = self%surface_index
      else
#  if _FABM_BOTTOM_INDEX_==0
        _VERTICAL_ITERATOR_ = self%bottom_index
#  elif !defined(_HORIZONTAL_IS_VECTORIZED_)
         _VERTICAL_ITERATOR_ = self%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#  endif
      end if
#endif

      do ivar=1,size(self%state_variables)
         read_index = self%state_variables(ivar)%target%read_indices%value
         if (self%interior_data_sources(read_index)==data_source_fabm) then
#if _FABM_BOTTOM_INDEX_==-1&&defined(_HORIZONTAL_IS_VECTORIZED_)
      if (flag==1) then
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
#  ifdef _HAS_MASK_
         self%data(read_index)%p _INDEX_GLOBAL_INTERIOR_(_START_:_STOP_) = unpack(cache%read(:,read_index),cache%mask,self%state_variables(ivar)%missing_value)
#  else
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            self%data(read_index)%p _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1) = cache%read _INDEX_SLICE_PLUS_1_(read_index)
         _HORIZONTAL_LOOP_END_
#  endif
#elif defined(_INTERIOR_IS_VECTORIZED_)
         self%data(read_index)%p _INDEX_LOCATION_ = cache%read(1,read_index)
#else
         self%data(read_index)%p _INDEX_LOCATION_ = cache%read(read_index)
#endif

#if _FABM_BOTTOM_INDEX_==-1&&defined(_HORIZONTAL_IS_VECTORIZED_)
      else
         ! Special case for bottom if vertical index of bottom point is variable.
#  ifdef _HAS_MASK_
         j = 0
         do _J_=1,_STOP_-_START_+1
            if (cache%mask(_J_)) then
               _VERTICAL_ITERATOR_ = self%bottom_indices _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1)
               j = j + 1
               self%data(read_index)%p _INDEX_GLOBAL_INTERIOR_(_START_+_J_-1) = cache%read(j,read_index)
            end if
         end do
#  else
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            _VERTICAL_ITERATOR_ = self%bottom_indices _INDEX_GLOBAL_HORIZONTAL_(_START_+_J_-1)
            self%data(read_index)%p _INDEX_GLOBAL_INTERIOR_(_START_+_I_-1) = cache%read _INDEX_SLICE_PLUS_1_(read_index)
         _HORIZONTAL_LOOP_END_
#  endif
      end if
#endif
         end if
      end do
   end if

   call end_horizontal_task(self,job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

end subroutine internal_check_horizontal_state

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-water exchange fluxes for all biogeochemical state variables.
! Positive values indicate fluxes into the ocean, negative values indicate fluxes
! out of the ocean. Units are tracer unit * m/s.
!
! !INTERFACE:
   subroutine fabm_do_surface(self _ARGUMENTS_HORIZONTAL_IN_,flux_pel,flux_sf)
!
! !INPUT PARAMETERS:
      class (type_model),                          intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
!
! !INPUT/OUTPUT PARAMETERS:
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(out)          :: flux_pel
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(out),optional :: flux_sf
!
! !LOCAL PARAMETERS:
      type (type_cache)         :: cache
      type (type_call), pointer :: call_node
      integer                   :: i,j,k
      _DECLARE_HORIZONTAL_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_do_surface')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
   call check_extents_2d(flux_pel,_STOP_-_START_+1,size(self%state_variables),'fabm_do_surface','flux_pel','stop-start+1, # interior state variables')
   if (present(flux_sf)) call check_extents_2d(flux_sf,_STOP_-_START_+1,size(self%surface_state_variables),'fabm_do_surface','flux_sf','stop-start+1, # surface state variables')
#  else
   call check_extents_1d(flux_pel,size(self%state_variables),'fabm_do_surface','flux_pel','# interior state variables')
   if (present(flux_sf)) call check_extents_1d(flux_sf,size(self%surface_state_variables),'fabm_do_surface','flux_sf','# surface state variables')
#  endif
#endif
      call fabm_process_job(self,self%do_surface_job _ARGUMENTS_HORIZONTAL_IN_)

      call begin_horizontal_task(self,self%do_surface_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

      call_node => self%do_surface_job%final_task%first_call
      do while (associated(call_node))
         if (call_node%source==source_do_horizontal) then
            call call_node%model%do_horizontal(_ARGUMENTS_HORIZONTAL_)
         else
            call call_node%model%do_surface(_ARGUMENTS_HORIZONTAL_)
         end if

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_hz))
            j = call_node%copy_commands_hz(i)%read_index
            k = call_node%copy_commands_hz(i)%write_index
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
               cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
            _HORIZONTAL_LOOP_END_
         end do

         call_node => call_node%next
      end do

      ! Compose surface fluxes for each interior state variable, combining model-specific contributions.
      flux_pel = 0.0_rk
      do i=1,size(self%state_variables)
         k = self%state_variables(i)%surface_flux_index
         _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(cache%write_hz,k,flux_pel,i,cache%mask,0.0_rk)
      end do

      ! Compose total sources-sinks for each surface-bound state variable, combining model-specific contributions.
      if (present(flux_sf)) then
         flux_sf = 0.0_rk
         do i=1,size(self%surface_state_variables)
            k = self%state_variables(i)%sms_index
            _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(cache%write_hz,k,flux_sf,i,cache%mask,0.0_rk)
         end do
      end if

      call end_horizontal_task(self,self%do_surface_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_do_surface
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process interaction between benthos and bottom layer of the
! pelagic. This calculates the fluxes into all bottom pelagic and benthic variables,
! in variable quantity per surface area per time. This typically implies variable units * m/s
! [bottom fluxes] for the pelagic, and variable units/s [temporal derivatives] for the benthos.
! Positive values denote state variable increases, negative values state variable decreases.
!
! !INTERFACE:
   subroutine fabm_do_bottom_rhs(self _ARGUMENTS_HORIZONTAL_IN_,flux_pel,flux_ben)
!
! !INPUT PARAMETERS:
   class (type_model),                          intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(inout) :: flux_pel,flux_ben
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: i,j,k
   _DECLARE_HORIZONTAL_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_do_bottom_rhs')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
   call check_extents_2d(flux_pel,_STOP_-_START_+1,size(self%state_variables),'fabm_do_bottom_rhs','flux_pel','stop-start+1, # interior state variables')
   call check_extents_2d(flux_ben,_STOP_-_START_+1,size(self%bottom_state_variables),'fabm_do_bottom_rhs','flux_ben','stop-start+1, # bottom state variables')
#  else
   call check_extents_1d(flux_pel,size(self%state_variables),'fabm_do_bottom_rhs','flux_pel','# interior state variables')
   call check_extents_1d(flux_ben,size(self%bottom_state_variables),'fabm_do_bottom_rhs','flux_ben','# bottom state variables')
#  endif
#endif

   call fabm_process_job(self,self%do_bottom_job _ARGUMENTS_HORIZONTAL_IN_)

   call begin_horizontal_task(self,self%do_bottom_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   call_node => self%do_bottom_job%final_task%first_call
   do while (associated(call_node))
      if (call_node%source==source_do_horizontal) then
         call call_node%model%do_horizontal(_ARGUMENTS_HORIZONTAL_)
      else
         call call_node%model%do_bottom(_ARGUMENTS_HORIZONTAL_)
      end if

      ! Copy outputs of interest to read cache so consecutive models can use it.
      _DO_CONCURRENT_(i,1,size(call_node%copy_commands_hz))
         j = call_node%copy_commands_hz(i)%read_index
         k = call_node%copy_commands_hz(i)%write_index
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
         _HORIZONTAL_LOOP_END_
      end do

      call_node => call_node%next
   end do

   ! Compose bottom fluxes for each interior state variable, combining model-specific contributions.
   do i=1,size(self%state_variables)
      k = self%state_variables(i)%bottom_flux_index
      _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(cache%write_hz,k,flux_pel,i,cache%mask,0.0_rk)
   end do

   ! Compose total sources-sinks for each bottom-bound state variable, combining model-specific contributions.
   do i=1,size(self%bottom_state_variables)
      k = self%bottom_state_variables(i)%sms_index
      _HORIZONTAL_UNPACK_AND_ADD_TO_PLUS_1_(cache%write_hz,k,flux_ben,i,cache%mask,0.0_rk)
   end do

   call end_horizontal_task(self,self%do_bottom_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_do_bottom_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process interaction between benthos and bottom layer of the
! pelagic. This calculates the fluxes between all bottom pelagic and benthic variables,
! in variable quantity per surface area per time. This typically imples variable units * m/s
! for the pelagic, and variable units/s for the benthos.
!
! !INTERFACE:
   subroutine fabm_do_bottom_ppdd(self _ARGUMENTS_HORIZONTAL_IN_,pp,dd,benthos_offset)
!
! !INPUT PARAMETERS:
   class (type_model),                          intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   integer,                                     intent(in)    :: benthos_offset
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_2_,intent(inout) :: pp,dd
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: i,j,k
   _DECLARE_HORIZONTAL_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_do_bottom_ppdd')
#endif

   call fabm_process_job(self,self%do_bottom_job _ARGUMENTS_HORIZONTAL_IN_)

   call begin_horizontal_task(self,self%do_bottom_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   call_node => self%do_bottom_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%do_bottom_ppdd(_ARGUMENTS_HORIZONTAL_,pp,dd,benthos_offset)

      ! Copy outputs of interest to read cache so consecutive models can use it.
      _DO_CONCURRENT_(i,1,size(call_node%copy_commands_hz))
         j = call_node%copy_commands_hz(i)%read_index
         k = call_node%copy_commands_hz(i)%write_index
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
         _HORIZONTAL_LOOP_END_
      end do

      call_node => call_node%next
   end do

   call end_horizontal_task(self,self%do_bottom_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_do_bottom_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical movement rates (m/s) for the bio state variables.
! Note that negative values indicate movement towards the bottom, e.g., sinking,
! and positive values indicate movemment towards the surface, e.g., floating.
!
! !INTERFACE:
   subroutine fabm_get_vertical_movement(self _ARGUMENTS_INTERIOR_IN_,velocity)
!
! !INPUT PARAMETERS:
   class (type_model),                   intent(inout) :: self
   _DECLARE_ARGUMENTS_INTERIOR_IN_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_EXT_SLICE_PLUS_1_,intent(out)   :: velocity
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: i,k
   _DECLARE_INTERIOR_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_interior_location(self _ARGUMENTS_INTERIOR_IN_,'fabm_get_vertical_movement')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
   call check_extents_2d(velocity,_STOP_-_START_+1,size(self%state_variables),'fabm_get_vertical_movement','velocity','stop-start+1, # interior state variables')
#  else
   call check_extents_1d(velocity,size(self%state_variables),'fabm_get_vertical_movement','velocity','# interior state variables')
#  endif
#endif

   call begin_interior_task(self,self%get_vertical_movement_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   ! Now allow models to overwrite with spatially-varying sinking rates - if any.
   call_node => self%get_vertical_movement_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%get_vertical_movement(_ARGUMENTS_INTERIOR_)
      call_node => call_node%next
   end do

   ! Compose total sources-sinks for each state variable, combining model-specific contributions.
   do i=1,size(self%state_variables)
      k = self%state_variables(i)%movement_index
      _UNPACK_TO_PLUS_1_(cache%write,k,velocity,i,cache%mask,0.0_rk)
   end do

   call end_interior_task(self,self%get_vertical_movement_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_get_vertical_movement
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine fabm_get_light_extinction(self _ARGUMENTS_INTERIOR_IN_,extinction)
!
! !INPUT PARAMETERS:
   class (type_model),            intent(inout) :: self
   _DECLARE_ARGUMENTS_INTERIOR_IN_
   real(rk) _DIMENSION_EXT_SLICE_,intent(out)   :: extinction
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: i,j,k
   _DECLARE_INTERIOR_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_interior_location(self _ARGUMENTS_INTERIOR_IN_,'fabm_get_light_extinction')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
   call check_extents_1d(extinction,_STOP_-_START_+1,'fabm_get_light_extinction','extinction','stop-start+1')
#  endif
#endif

   call begin_interior_task(self,self%get_light_extinction_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   ! Call all models that calculate extinction components to make sure the extinction diagnostic is up to date.
   call_node => self%get_light_extinction_job%final_task%first_call
   do while (associated(call_node))
      if (call_node%source==source_do) then
         call call_node%model%do(_ARGUMENTS_INTERIOR_)
      elseif (call_node%source==source_get_light_extinction) then
         call call_node%model%get_light_extinction(_ARGUMENTS_INTERIOR_)
      end if

      ! Copy outputs of interest to read cache so consecutive models can use it.
      _DO_CONCURRENT_(i,1,size(call_node%copy_commands_int))
         j = call_node%copy_commands_int(i)%read_index
         k = call_node%copy_commands_int(i)%write_index
         _CONCURRENT_LOOP_BEGIN_
            cache%read _INDEX_SLICE_PLUS_1_(j) = cache%write _INDEX_SLICE_PLUS_1_(k)
         _LOOP_END_
      end do

      ! Move to next model
      call_node => call_node%next
   end do

   ! NB the extinction value is taken from the *read* cache.
   _UNPACK_(cache%read,self%extinction_index,extinction,cache%mask,0.0_rk)

   call end_interior_task(self,self%get_light_extinction_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate the light field
!
! !INTERFACE:
   subroutine fabm_get_light(self _ARGUMENTS_VERTICAL_IN_)
!
! !INPUT PARAMETERS:
   class (type_model), intent(inout) :: self
   _DECLARE_ARGUMENTS_VERTICAL_IN_
!
!EOP
!-----------------------------------------------------------------------
!BOC

   end subroutine fabm_get_light
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get biogeochemistry-induced change in wind drag
! This feedback is represented by a scale factor (0-1), with
! 1 leaving drag unchanged, and 0 suppressing it completely.
!
! !INTERFACE:
   subroutine fabm_get_drag(self _ARGUMENTS_HORIZONTAL_IN_,drag)
!
! !INPUT PARAMETERS:
   class (type_model),                   intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   real(rk) _DIMENSION_HORIZONTAL_SLICE_,intent(out)   :: drag
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   _DECLARE_HORIZONTAL_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_get_drag')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
   call check_extents_1d(drag,_STOP_-_START_+1,'fabm_get_drag','drag','stop-start+1')
#  endif
#endif

   call begin_horizontal_task(self,self%get_drag_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   drag = 1.0_rk
   call_node => self%get_drag_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%get_drag(_ARGUMENTS_HORIZONTAL_,drag)
      call_node => call_node%next
   end do

   call end_horizontal_task(self,self%get_drag_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_get_drag
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get biogeochemistry-induced change in surface albedo
! This feedback is represented by an albedo value (0-1) relating to biogeochemistry.
!
! !INTERFACE:
   subroutine fabm_get_albedo(self _ARGUMENTS_HORIZONTAL_IN_,albedo)
!
! !INPUT PARAMETERS:
   class (type_model),                   intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   real(rk) _DIMENSION_HORIZONTAL_SLICE_,intent(out)   :: albedo
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   _DECLARE_HORIZONTAL_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_get_albedo')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
   call check_extents_1d(albedo,_STOP_-_START_+1,'fabm_get_albedo','albedo','stop-start+1')
#  endif
#endif

   call begin_horizontal_task(self,self%get_albedo_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   albedo = 0.0_rk
   call_node => self%get_albedo_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%get_albedo(_ARGUMENTS_HORIZONTAL_,albedo)
      call_node => call_node%next
   end do

   call end_horizontal_task(self,self%get_albedo_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_get_albedo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities in pelagic domain
!
! !INTERFACE:
   subroutine fabm_get_conserved_quantities(self _ARGUMENTS_INTERIOR_IN_,sums)
!
! !INPUT PARAMETERS:
   class (type_model),                   intent(inout) :: self
   _DECLARE_ARGUMENTS_INTERIOR_IN_
   real(rk) _DIMENSION_EXT_SLICE_PLUS_1_,intent(out)   :: sums
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: i,j,k
   _DECLARE_INTERIOR_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_interior_location(self _ARGUMENTS_INTERIOR_IN_,'fabm_get_conserved_quantities')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
   call check_extents_2d(sums,_STOP_-_START_+1,size(self%conserved_quantities),'fabm_get_conserved_quantities','sums','stop-start+1, # conserved quantities')
#  else
   call check_extents_1d(sums,size(self%conserved_quantities),'fabm_get_conserved_quantities','sums','# conserved quantities')
#  endif
#endif

   call begin_interior_task(self,self%get_conserved_quantities_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   call_node => self%get_conserved_quantities_job%final_task%first_call
   do while (associated(call_node))
      call call_node%model%do(_ARGUMENTS_INTERIOR_)

      ! Copy outputs of interest to read cache so consecutive models can use it.
      _DO_CONCURRENT_(i,1,size(call_node%copy_commands_int))
         j = call_node%copy_commands_int(i)%read_index
         k = call_node%copy_commands_int(i)%write_index
         _CONCURRENT_LOOP_BEGIN_
            cache%read _INDEX_SLICE_PLUS_1_(j) = cache%write _INDEX_SLICE_PLUS_1_(k)
         _LOOP_END_
      end do

      ! Move to next model
      call_node => call_node%next
   end do

   do i=1,size(self%conserved_quantities)
      _UNPACK_TO_PLUS_1_(cache%read,self%conserved_quantities(i)%index,sums,i,cache%mask,0.0_rk)
   end do

   call end_interior_task(self,self%get_conserved_quantities_job%final_task,cache _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities at top and bottom boundaries
!
! !INTERFACE:
   subroutine fabm_get_horizontal_conserved_quantities(self _ARGUMENTS_HORIZONTAL_IN_,sums)
!
! !INPUT PARAMETERS:
   class (type_model),                          intent(inout) :: self
   _DECLARE_ARGUMENTS_HORIZONTAL_IN_
   real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(out)   :: sums
!
! !LOCAL PARAMETERS:
   type (type_cache)         :: cache
   type (type_call), pointer :: call_node
   integer                   :: i,j,k
   _DECLARE_HORIZONTAL_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,'fabm_get_horizontal_conserved_quantities')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
   call check_extents_2d(sums,_STOP_-_START_+1,size(self%conserved_quantities),'fabm_get_horizontal_conserved_quantities','sums','stop-start+1, # conserved quantities')
#  else
   call check_extents_1d(sums,size(self%conserved_quantities),'fabm_get_horizontal_conserved_quantities','sums','# conserved quantities')
#  endif
#endif

   call begin_horizontal_task(self,self%get_horizontal_conserved_quantities_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   call_node => self%get_horizontal_conserved_quantities_job%final_task%first_call
   do while (associated(call_node))
      if (call_node%source==source_do_horizontal) then
         call call_node%model%do_horizontal(_ARGUMENTS_HORIZONTAL_)
      else
         call call_node%model%do_bottom(_ARGUMENTS_HORIZONTAL_)
      end if

      ! Copy outputs of interest to read cache so consecutive models can use it.
      _DO_CONCURRENT_(i,1,size(call_node%copy_commands_hz))
         j = call_node%copy_commands_hz(i)%read_index
         k = call_node%copy_commands_hz(i)%write_index
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
            cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
         _HORIZONTAL_LOOP_END_
      end do

      ! Move to next model
      call_node => call_node%next
   end do

   do i=1,size(self%conserved_quantities)
      _HORIZONTAL_UNPACK_TO_PLUS_1_(cache%read_hz,self%conserved_quantities(i)%horizontal_index,sums,i,cache%mask,0.0_rk)
   end do

   call end_horizontal_task(self,self%get_horizontal_conserved_quantities_job%final_task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_get_horizontal_conserved_quantities
!EOC

   subroutine fabm_process_job_all(self,job _ARGUMENTS_HORIZONTAL_LOCATION_RANGE_)
      class (type_model), intent(inout), target :: self
      type (type_job),    intent(in)            :: job
      _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_RANGE_

      type (type_task),pointer :: task

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      ! Jobs must be applied across the entire depth range (if any),
      ! so we set vertical start and stop indices here.
      integer :: _VERTICAL_START_,_VERTICAL_STOP_
      _VERTICAL_START_ = 1
      _VERTICAL_STOP_ = self%domain_size(_FABM_DEPTH_DIMENSION_INDEX_)
#endif

      task => job%first_task
      do while (associated(task))
         select case (task%operation)
         case (source_do)
            call fabm_process_interior_all(self,task _ARGUMENTS_LOCATION_RANGE_)
         case (source_do_surface,source_do_bottom,source_do_horizontal)
            call fabm_process_horizontal_all(self,task _ARGUMENTS_HORIZONTAL_LOCATION_RANGE_)
         case (source_do_column)
            call fabm_process_vertical_all(self,task _ARGUMENTS_LOCATION_RANGE_)
         end select
         task => task%next
      end do
   end subroutine fabm_process_job_all

   subroutine fabm_process_interior_all(self,job _ARGUMENTS_LOCATION_RANGE_)
      class (type_model),intent(inout), target :: self
      type (type_task),  intent(in)            :: job
      _DECLARE_ARGUMENTS_LOCATION_RANGE_

      _DECLARE_LOCATION_

      _BEGIN_OUTER_INTERIOR_LOOP_
         call fabm_process_interior_slice(self,job _ARGUMENTS_INTERIOR_IN_)
      _END_OUTER_INTERIOR_LOOP_
   end subroutine fabm_process_interior_all

   subroutine fabm_process_horizontal_all(self,job _ARGUMENTS_HORIZONTAL_LOCATION_RANGE_)
      class (type_model),intent(inout), target :: self
      type (type_task),  intent(in)            :: job
      _DECLARE_ARGUMENTS_HORIZONTAL_LOCATION_RANGE_

      _DECLARE_LOCATION_

      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_process_horizontal_slice(self,job _ARGUMENTS_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
   end subroutine fabm_process_horizontal_all

   subroutine fabm_process_vertical_all(self,job _ARGUMENTS_LOCATION_RANGE_)
      class (type_model),intent(inout), target :: self
      type (type_task),  intent(in)            :: job
      _DECLARE_ARGUMENTS_LOCATION_RANGE_

      _DECLARE_LOCATION_

      _BEGIN_OUTER_VERTICAL_LOOP_
         call fabm_process_vertical_slice(self,job _ARGUMENTS_VERTICAL_IN_)
      _END_OUTER_VERTICAL_LOOP_
   end subroutine fabm_process_vertical_all

   subroutine fabm_process_job(self,job _ARGUMENTS_HORIZONTAL_IN_)
      class (type_model),intent(inout) :: self
      type (type_job),   intent(in)    :: job
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_

      type (type_task),pointer :: task

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      ! Jobs must be applied across the entire depth range (if any),
      ! so we set vertical start and stop indices here.
      integer :: _VERTICAL_START_,_VERTICAL_STOP_
#  if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_
      integer :: _ITERATOR_
#  endif
#  if _FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_
      integer :: _VERTICAL_ITERATOR_
#  endif
      _VERTICAL_START_ = 1
      _VERTICAL_STOP_ = self%domain_size(_FABM_DEPTH_DIMENSION_INDEX_)
#endif

      task => job%first_task
      do while (associated(task))
         select case (task%operation)
         case (source_do)
#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_
            do _VERTICAL_ITERATOR_=_VERTICAL_START_,_VERTICAL_STOP_
#endif
               call fabm_process_interior_slice(self,task _ARGUMENTS_INTERIOR_IN_)
#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_
            end do
#endif
         case (source_do_surface,source_do_bottom,source_do_horizontal)
            call fabm_process_horizontal_slice(self,task _ARGUMENTS_HORIZONTAL_IN_)
         case (source_do_column)
#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_
            do _ITERATOR_=_START_,_STOP_
#endif
               call fabm_process_vertical_slice(self,task _ARGUMENTS_VERTICAL_IN_)
#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_
            end do
#endif
         end select
         task => task%next
      end do
   end subroutine fabm_process_job

   subroutine fabm_process_interior_slice(self,task _ARGUMENTS_INTERIOR_IN_)
      class (type_model),intent(inout), target :: self
      type (type_task),  intent(in)            :: task
      _DECLARE_ARGUMENTS_INTERIOR_IN_

      type (type_cache)         :: cache
      type (type_call), pointer :: call_node
      integer                   :: i,j,k
      _DECLARE_INTERIOR_INDICES_

      call begin_interior_task(self,task,cache _ARGUMENTS_INTERIOR_IN_)

      call_node => task%first_call
      do while (associated(call_node))
         call call_node%model%do(_ARGUMENTS_INTERIOR_)

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_int))
            j = call_node%copy_commands_int(i)%read_index
            k = call_node%copy_commands_int(i)%write_index
            _CONCURRENT_LOOP_BEGIN_
               cache%read _INDEX_SLICE_PLUS_1_(j) = cache%write _INDEX_SLICE_PLUS_1_(k)
            _LOOP_END_
         end do

         ! Move to next model
         call_node => call_node%next
      end do

      call end_interior_task(self,task,cache _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_process_interior_slice

   subroutine fabm_process_horizontal_slice(self,task _ARGUMENTS_HORIZONTAL_IN_)
      class (type_model),intent(inout) :: self
      type (type_task),  intent(in)    :: task
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_

      type (type_cache)         :: cache
      type (type_call), pointer :: call_node
      integer                   :: i,j,k
      _DECLARE_HORIZONTAL_INDICES_

      call begin_horizontal_task(self,task,cache _ARGUMENTS_HORIZONTAL_IN_)

      call_node => task%first_call
      do while (associated(call_node))
         select case (call_node%source)
         case (source_do_surface);    call call_node%model%do_surface   (_ARGUMENTS_HORIZONTAL_)
         case (source_do_bottom);     call call_node%model%do_bottom    (_ARGUMENTS_HORIZONTAL_)
         case (source_do_horizontal); call call_node%model%do_horizontal(_ARGUMENTS_HORIZONTAL_)
         end select

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_hz))
            j = call_node%copy_commands_hz(i)%read_index
            k = call_node%copy_commands_hz(i)%write_index
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_
               cache%read_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = cache%write_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
            _HORIZONTAL_LOOP_END_
         end do

         call_node => call_node%next
      end do

      call end_horizontal_task(self,task,cache _ARGUMENTS_HORIZONTAL_IN_)

   end subroutine fabm_process_horizontal_slice

   subroutine fabm_process_vertical_slice(self,task _ARGUMENTS_VERTICAL_IN_)
      class (type_model),intent(inout) :: self
      type (type_task),  intent(in)    :: task
      _DECLARE_ARGUMENTS_VERTICAL_IN_

      type (type_cache)         :: cache
      type (type_call), pointer :: call_node
      integer                   :: i,j,k
      _DECLARE_VERTICAL_INDICES_

      call begin_vertical_task(self,task,cache _ARGUMENTS_VERTICAL_IN_)

      call_node => task%first_call
      do while (associated(call_node))
         call call_node%model%get_light(_ARGUMENTS_VERTICAL_)

         ! Copy outputs of interest to read cache so consecutive models can use it.
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_int))
            j = call_node%copy_commands_int(i)%read_index
            k = call_node%copy_commands_int(i)%write_index
            _CONCURRENT_VERTICAL_LOOP_BEGIN_
               cache%read _INDEX_SLICE_PLUS_1_(j) = cache%write _INDEX_SLICE_PLUS_1_(k)
            _VERTICAL_LOOP_END_
         end do
         _DO_CONCURRENT_(i,1,size(call_node%copy_commands_hz))
            j = call_node%copy_commands_hz(i)%read_index
            k = call_node%copy_commands_hz(i)%write_index
#ifdef _HORIZONTAL_IS_VECTORIZED_
            cache%read_hz(1,j) = cache%write_hz(1,k)
#else
            cache%read_hz(j) = cache%write_hz(k)
#endif
         end do

         call_node => call_node%next
      end do

      call end_vertical_task(self,task,cache _ARGUMENTS_VERTICAL_IN_)

   end subroutine fabm_process_vertical_slice


subroutine fabm_update_time(self,t)
   class (type_model), intent(inout) :: self
   real(rk),           intent(in)    :: t

   class (type_expression),pointer :: expression

   expression => self%root%first_expression
   do while (associated(expression))
      select type (expression)
         class is (type_bulk_temporal_mean)
            call update_bulk_temporal_mean(expression)
         class is (type_horizontal_temporal_mean)
            call update_horizontal_temporal_mean(expression)
      end select
      expression => expression%next
   end do

contains

   subroutine update_bulk_temporal_mean(expression)
      class (type_bulk_temporal_mean), intent(inout) :: expression
      integer  :: i
      real(rk) :: weight_right,frac_outside

      if (expression%ioldest==-1) then
         ! Start of simulation
         expression%next_save_time = t + expression%period/expression%n
         expression%ioldest = 1
      end if
      do while (t>=expression%next_save_time)
         ! Weight for linear interpolation between last stored point and current point, to get at values for desired time.
         weight_right = (expression%next_save_time-expression%last_time)/(t-expression%last_time)

         ! For temporal means:
         ! - remove contribution of oldest point from historical mean (@ n + 2)
         ! - linearly interpolate to desired time (@ ioldest), by computing a weighted mean of the current value (data(expression%in)%p) and the previous value (@ n + 1)
         ! - add contribution of new point to historical mean
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) &
            - expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest)/expression%n
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest) = (1.0_rk-weight_right)*expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+1) &
            + weight_right*self%data(expression%in)%p
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) &
            + expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest)/expression%n

         ! Compute next time for which we want to store output
         expression%next_save_time = expression%next_save_time + expression%period/expression%n

         ! If we just completed the first entire history, compute the running mean and record that it is now valid.
         if (expression%ioldest == expression%n .and. .not. expression%valid) then
            expression%valid = .true.
            expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) = 0.0_rk
            do i=1,expression%n
               expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) + expression%history(_PREARG_LOCATION_DIMENSIONS_ i)
            end do
            expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2)/expression%n
         end if

         ! Increment index for oldest time point
         expression%ioldest = expression%ioldest + 1
         if (expression%ioldest>expression%n) expression%ioldest = 1
      end do

      ! Store current value to enable linear interpolation to next output time in subsequent call.
      expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+1) = self%data(expression%in)%p

      if (expression%valid) then
         ! We have a full history. To compute the temporal mean:
         ! - store values at current time step (@ n + 1)
         ! - set mean (@ n + 3) to historical mean (@ n + 2) but account for change since most recent point in history.

         ! Compute extent of time period outside history
         frac_outside = (t-(expression%next_save_time-expression%period/expression%n))/expression%period

         ! Set corrected running mean (move window by removing part of the start, and appending to the end)
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+3) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) &
            + frac_outside*(-expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest) + expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+1))
      else
         ! We do not have a full history yet; set temporal mean to msising value
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+3) = expression%missing_value
      end if

      expression%last_time = t
   end subroutine

   subroutine update_horizontal_temporal_mean(expression)
      class (type_horizontal_temporal_mean), intent(inout) :: expression
      integer  :: i
      real(rk) :: weight_right,frac_outside

      if (expression%ioldest==-1) then
         ! Start of simulation; set entire history equal to current value.
         do i=1,expression%n+3
            expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ i) = self%data_hz(expression%in)%p
         end do
         expression%next_save_time = t + expression%period/expression%n
         expression%ioldest = 1
      end if
      do while (t>=expression%next_save_time)
         ! Weight for linear interpolation between last stored point and current point, to get at values for desired time.
         weight_right = (expression%next_save_time-expression%last_time)/(t-expression%last_time)

         ! For temporal means:
         ! - remove contribution of oldest point from historical mean
         ! - linearly interpolate to desired time
         ! - add contribution of new point to historical mean
         expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+2) = expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+2) &
            - expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%ioldest)/expression%n
         expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%ioldest) = (1.0_rk-weight_right)*expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+1) &
            + weight_right*self%data_hz(expression%in)%p
         expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+2) = expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+2) &
            + expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%ioldest)/expression%n

         ! Compute next time for which we want to store output
         expression%next_save_time = expression%next_save_time + expression%period/expression%n

         ! Increment index for oldest time point
         expression%ioldest = expression%ioldest + 1
         if (expression%ioldest>expression%n) expression%ioldest = 1
      end do

      ! Compute extent of time period outside history
      frac_outside = (t-(expression%next_save_time-expression%period/expression%n))/expression%period

      ! For temporal means:
      ! - store values at current time step
      ! - for current mean, use historical mean but account for change since most recent point in history.
      expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+1) = self%data_hz(expression%in)%p
      expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+3) = expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+2) &
         + frac_outside*(-expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%ioldest) + expression%history(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ expression%n+1))

      expression%last_time = t
   end subroutine

end subroutine fabm_update_time

subroutine assign_write_indices(link_list,domain,count)
   type (type_link_list),intent(inout) :: link_list
   integer,              intent(in)    :: domain
   integer,              intent(inout) :: count

   type (type_link),pointer :: link

   link => link_list%first
   do while (associated(link))
      if (link%target%domain==domain) then
         if (.not.link%target%write_indices%is_empty().and.link%target%write_indices%value==-1.and..not.associated(link%target%write_owner)) then
            count = count + 1
            call link%target%write_indices%set_value(count)
         end if
      end if
      link => link%next
   end do
end subroutine assign_write_indices

subroutine filter_readable_variable_registry(self)
   class (type_model),intent(inout),target :: self

   integer :: i
   type (type_variable_node),pointer :: variable_node

   i = 0
   variable_node => self%variable_register%interior_read%first
   do while (associated(variable_node))
      i = i + 1
      variable_node%target%in_read_registry = associated(self%data(i)%p)
      variable_node => variable_node%next
   end do

   i = 0
   variable_node => self%variable_register%horizontal_read%first
   do while (associated(variable_node))
      i = i + 1
      variable_node%target%in_read_registry = associated(self%data_hz(i)%p)
      variable_node => variable_node%next
   end do

   i = 0
   variable_node => self%variable_register%scalar_read%first
   do while (associated(variable_node))
      i = i + 1
      variable_node%target%in_read_registry = associated(self%data_scalar(i)%p)
      variable_node => variable_node%next
   end do
end subroutine filter_readable_variable_registry

subroutine classify_variables(self)
   class (type_model),intent(inout),target :: self

   type (type_link),                               pointer :: link,newlink
   type (type_state_variable_info),                pointer :: statevar
   type (type_horizontal_state_variable_info),     pointer :: hz_statevar
   type (type_diagnostic_variable_info),           pointer :: diagvar
   type (type_horizontal_diagnostic_variable_info),pointer :: hz_diagvar
   type (type_conserved_quantity_info),            pointer :: consvar
   type (type_internal_variable),                  pointer :: object
   integer                                                 :: nstate,nstate_bot,nstate_surf,ndiag,ndiag_hz,ncons

   type (type_aggregate_variable_list)         :: aggregate_variable_list
   type (type_aggregate_variable),     pointer :: aggregate_variable
   type (type_set)                             :: dependencies,dependencies_hz,dependencies_scalar
   type (type_model_list_node),        pointer :: model_node
   type (type_standard_variable_node), pointer :: standard_variables_node

   ! Determine the order in which individual biogeochemical models should be called.
   call build_call_list(self%root,self%models)

   ! Consistency check (of find_dependencies): does every model (except the root) appear exactly once in the call list?
   call test_presence(self,self%root)

   ! Create a list of all non-coupled variables, ordered according to the order in which the models will be called.
   ! This promotes sequential memory access and hopefully increases cache hits.
   model_node => self%models%first
   do while (associated(model_node))
      link => model_node%model%links%first
      do while (associated(link))
         if (associated(link%target,link%original).and.index(link%name,'/')==0) &
            newlink => self%links_postcoupling%append(link%target,link%target%name)
         link => link%next
      end do
      model_node => model_node%next
   end do

   ! Count number of conserved quantities and allocate associated array.
   aggregate_variable_list = collect_aggregate_variables(self%root)
   ncons = 0
   aggregate_variable => aggregate_variable_list%first
   do while (associated(aggregate_variable))
      if (aggregate_variable%standard_variable%conserved) ncons = ncons + 1
      aggregate_variable => aggregate_variable%next
   end do
   allocate(self%conserved_quantities(ncons))

   ! Fill list of conserved quantities.
   ! This must be done before building the final authoratitive list of diagnostic variables,
   ! as the calls to append_data_pointer affect the global identifier of diagnostic variables
   ! by adding another pointer that must be set.
   ncons = 0
   aggregate_variable => aggregate_variable_list%first
   do while (associated(aggregate_variable))
      if (aggregate_variable%standard_variable%conserved) then
         ncons = ncons + 1
         consvar => self%conserved_quantities(ncons)
         consvar%standard_variable = aggregate_variable%standard_variable
         consvar%name = trim(consvar%standard_variable%name)
         consvar%units = trim(consvar%standard_variable%units)
         consvar%long_name = trim(consvar%standard_variable%name)
         consvar%path = trim(consvar%standard_variable%name)
         consvar%target => self%root%find_object(trim(aggregate_variable%standard_variable%name))
         if (.not.associated(consvar%target)) call driver%fatal_error('classify_variables', &
            'BUG: conserved quantity '//trim(aggregate_variable%standard_variable%name)//' was not created')
         call consvar%target%read_indices%append(consvar%index)
         consvar%target_hz => self%root%find_object(trim(aggregate_variable%standard_variable%name)//'_at_interfaces')
         if (.not.associated(consvar%target_hz)) call driver%fatal_error('classify_variables', &
            'BUG: conserved quantity '//trim(aggregate_variable%standard_variable%name)//'_at_interfaces was not created')
         call consvar%target_hz%read_indices%append(consvar%horizontal_index)
      end if
      aggregate_variable => aggregate_variable%next
   end do

   ! Get link to extinction variable.
   self%extinction_target => self%root%find_object(trim(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux%name))
   if (.not.associated(self%extinction_target)) call driver%fatal_error('classify_variables', &
      'BUG: variable attenuation_coefficient_of_photosynthetic_radiative_flux was not created')
   call self%extinction_target%read_indices%append(self%extinction_index)

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
      select case (object%domain)
         case (domain_interior)
            if (.not.object%write_indices%is_empty()) then
               ! Interior diagnostic variable.
               ndiag = ndiag+1
            elseif (.not.object%state_indices%is_empty().and..not.object%fake_state_variable) then
               ! Interior state variable.
               select case (object%presence)
                  case (presence_internal)
                     nstate = nstate+1
                     call object%state_indices%set_value(nstate)
                  case (presence_external_required)
                     call fatal_error('classify_variables', &
                        'Variable '//trim(link%name)//' must be coupled to an existing state variable.')
                  case default
                     continue
               end select
            end if
         case (domain_horizontal,domain_surface,domain_bottom)
            if (.not.object%write_indices%is_empty()) then
               ! Horizontal diagnostic variable.
               ndiag_hz = ndiag_hz+1
            elseif (.not.object%state_indices%is_empty().and..not.object%fake_state_variable) then
               ! Horizontal state variable.
               select case (object%presence)
                  case (presence_internal)
                     select case (object%domain)
                        case (domain_bottom)
                           ! Bottom state variable
                           nstate_bot = nstate_bot+1
                           call object%state_indices%set_value(nstate_bot)
                        case (domain_surface)
                           ! Surface state variable
                           nstate_surf = nstate_surf+1
                           call object%state_indices%set_value(nstate_surf)
                     end select
                  case (presence_external_required)
                     call fatal_error('classify_variables', &
                        'Variable '//trim(link%name)//' must be coupled to an existing state variable.')
               end select
            end if
      end select
      link => link%next
   end do

   ! Allocate arrays with variable information that will be accessed by the host model.
   allocate(self%state_variables                (nstate))
   allocate(self%bottom_state_variables         (nstate_bot))
   allocate(self%surface_state_variables        (nstate_surf))
   allocate(self%diagnostic_variables           (ndiag))
   allocate(self%horizontal_diagnostic_variables(ndiag_hz))

   ! Set pointers for backward compatibility (pre 2013-06-15)
   ! Note: this must be done AFTER allocation of the target arrays, above!
   self%state_variables_ben => self%bottom_state_variables
   self%diagnostic_variables_hz => self%horizontal_diagnostic_variables

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
            if (.not.object%write_indices%is_empty()) then
               ! Interior diagnostic variable
               ndiag = ndiag + 1
               diagvar => self%diagnostic_variables(ndiag)
               call copy_variable_metadata(object,diagvar)
               if (associated(object%standard_variables%first)) then
                  select type (standard_variable=>object%standard_variables%first%p)
                     type is (type_bulk_standard_variable)
                        diagvar%standard_variable = standard_variable
                  end select
               end if
               diagvar%time_treatment    = output2time_treatment(diagvar%output)
               diagvar%save = diagvar%output/=output_none
               diagvar%source = object%source
            elseif (object%presence==presence_internal.and..not.object%state_indices%is_empty().and..not.object%fake_state_variable) then
               ! Interior state variable
               nstate = nstate + 1
               statevar => self%state_variables(nstate)
               call copy_variable_metadata(object,statevar)
               if (associated(object%standard_variables%first)) then
                  select type (standard_variable=>object%standard_variables%first%p)
                     type is (type_bulk_standard_variable)
                        statevar%standard_variable = standard_variable
                  end select
               end if
               statevar%initial_value             = object%initial_value
               statevar%vertical_movement         = object%vertical_movement
               statevar%no_precipitation_dilution = object%no_precipitation_dilution
               statevar%no_river_dilution         = object%no_river_dilution
               if (.not.object%sms_sum%target%write_indices%is_empty()) call object%sms_sum%target%write_indices%append(statevar%sms_index)
               if (.not.object%surface_flux_sum%target%write_indices%is_empty()) call object%surface_flux_sum%target%write_indices%append(statevar%surface_flux_index)
               if (.not.object%bottom_flux_sum%target%write_indices%is_empty()) call object%bottom_flux_sum%target%write_indices%append(statevar%bottom_flux_index)
               if (.not.object%movement_diagnostic%target%write_indices%is_empty()) call object%movement_diagnostic%target%write_indices%append(statevar%movement_index)
            end if
         case (domain_horizontal,domain_surface,domain_bottom)
            if (.not.object%write_indices%is_empty()) then
               ! Horizontal diagnostic variable
               ndiag_hz = ndiag_hz + 1
               hz_diagvar => self%horizontal_diagnostic_variables(ndiag_hz)
               call copy_variable_metadata(object,hz_diagvar)
               if (associated(object%standard_variables%first)) then
                  select type (standard_variable=>object%standard_variables%first%p)
                     type is (type_horizontal_standard_variable)
                        hz_diagvar%standard_variable = standard_variable
                  end select
               end if
               hz_diagvar%time_treatment    = output2time_treatment(hz_diagvar%output)
               hz_diagvar%save = hz_diagvar%output/=output_none
               hz_diagvar%source = object%source
            elseif (object%presence==presence_internal.and..not.object%state_indices%is_empty().and..not.object%fake_state_variable) then
               ! Horizontal state variable
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
               if (associated(object%standard_variables%first)) then
                  select type (standard_variable=>object%standard_variables%first%p)
                     type is (type_horizontal_standard_variable)
                        hz_statevar%standard_variable = standard_variable
                  end select
               end if
               hz_statevar%initial_value     = object%initial_value
               if (.not.object%sms_sum%target%write_indices%is_empty()) call object%sms_sum%target%write_indices%append(hz_statevar%sms_index)
            end if
      end select
      link => link%next
   end do

   ! Create lists of variables that may be provided by the host model.
   ! These lists include external dependencies, as well as the model's own state variables,
   ! which may be overridden by the host.
   link => self%root%links%first
   do while (associated(link))
      object =>link%target
      if (.not.object%read_indices%is_empty().and. &
          .not.(object%presence==presence_external_optional.and..not.object%state_indices%is_empty())) then
         select case (object%domain)
            case (domain_interior);                                call dependencies%add(link%name)
            case (domain_horizontal,domain_surface,domain_bottom); call dependencies_hz%add(link%name)
            case (domain_scalar);                                  call dependencies_scalar%add(link%name)
         end select

         standard_variables_node => object%standard_variables%first
         do while (associated(standard_variables_node))
            if (standard_variables_node%p%name/='') then
               select case (object%domain)
                  case (domain_interior);                                call dependencies%add(standard_variables_node%p%name)
                  case (domain_horizontal,domain_surface,domain_bottom); call dependencies_hz%add(standard_variables_node%p%name)
                  case (domain_scalar);                                  call dependencies_scalar%add(standard_variables_node%p%name)
               end select
            end if
            standard_variables_node => standard_variables_node%next
         end do
      end if
      link => link%next
   end do
   call dependencies%to_array(self%dependencies)
   call dependencies_hz%to_array(self%dependencies_hz)
   call dependencies_scalar%to_array(self%dependencies_scalar)
contains
subroutine copy_variable_metadata(internal_variable,external_variable)
   class (type_external_variable),intent(inout) :: external_variable
   type (type_internal_variable),intent(in),target :: internal_variable

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
      external_variable%long_name = trim(owner%long_name)//' '//trim(external_variable%long_name)
      owner => owner%parent
   end do

   call external_variable%properties%update(internal_variable%properties)
end subroutine
end subroutine classify_variables

   recursive subroutine build_call_list(self,list)
      class (type_base_model),intent(in),target   :: self
      type (type_model_list), intent(inout)       :: list

      type (type_model_list_node),pointer :: node

      call find_dependencies(self,list)
      node => self%children%first
      do while (associated(node))
         call build_call_list(node%model,list)
         node => node%next
      end do
   end subroutine build_call_list

   subroutine require_flux_computation(self,link_list,domain)
      type (type_job),      intent(inout) :: self
      type (type_link_list),intent(in)    :: link_list
      integer,              intent(in)    :: domain

      type (type_link), pointer :: link

      link => link_list%first
      do while (associated(link))
         select case (domain)
         case (domain_interior)
            if (link%target%domain==domain_interior.and.associated(link%target%sms_sum)) &
               call self%request_variable(link%target%sms_sum%target)
         case (domain_bottom)
            if (link%target%domain==domain_bottom.and.associated(link%target%sms_sum)) then
               call self%request_variable(link%target%sms_sum%target)
            elseif (link%target%domain==domain_interior.and.associated(link%target%bottom_flux_sum)) then
               call self%request_variable(link%target%bottom_flux_sum%target)
            end if
         case (domain_surface)
            if (link%target%domain==domain_surface.and.associated(link%target%sms_sum)) then
               call self%request_variable(link%target%sms_sum%target)
            elseif (link%target%domain==domain_interior.and.associated(link%target%surface_flux_sum)) then
               call self%request_variable(link%target%surface_flux_sum%target)
            end if
         end select
         link => link%next
      end do
   end subroutine require_flux_computation

   recursive subroutine require_call_all(self,model,source)
      type (type_job),        intent(inout) :: self
      class (type_base_model),intent(in)    :: model
      integer,                intent(in)    :: source

      type (type_model_list_node),pointer :: node

      node => model%children%first
      do while (associated(node))
         call require_call_all(self,node%model,source)
         call self%request_call(node%model,source)
         node => node%next
      end do
   end subroutine require_call_all

   subroutine require_call_all_with_state(self,link_list,domain,source)
      type (type_job),      intent(inout) :: self
      type (type_link_list),intent(in)    :: link_list
      integer,              intent(in)    :: domain
      integer,              intent(in)    :: source

      type (type_link), pointer :: link

      link => link_list%first
      do while (associated(link))
         if (link%target%domain==domain.and..not.link%original%state_indices%is_empty().and..not.link%target%fake_state_variable) &
            call self%request_call(link%original%owner,source)
         link => link%next
      end do
   end subroutine require_call_all_with_state

   subroutine check_interior_location(self _ARGUMENTS_INTERIOR_IN_,routine)
      class (type_model),intent(in) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      character(len=*), intent(in) :: routine

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_loop(_START_,_STOP_,self%domain_size(_FABM_VECTORIZED_DIMENSION_INDEX_),routine)
#endif
#if _FABM_DIMENSION_COUNT_>0&&_FABM_VECTORIZED_DIMENSION_INDEX_!=1
      call check_index(i__,self%domain_size(1),routine,'i')
#endif
#if _FABM_DIMENSION_COUNT_>1&&_FABM_VECTORIZED_DIMENSION_INDEX_!=2
      call check_index(j__,self%domain_size(2),routine,'j')
#endif
#if _FABM_DIMENSION_COUNT_>2&&_FABM_VECTORIZED_DIMENSION_INDEX_!=3
      call check_index(k__,self%domain_size(3),routine,'k')
#endif
   end subroutine check_interior_location

   subroutine check_horizontal_location(self _ARGUMENTS_HORIZONTAL_IN_,routine)
      class (type_model),intent(in) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      character(len=*), intent(in) :: routine

#ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_loop(_START_,_STOP_,self%domain_size(_FABM_VECTORIZED_DIMENSION_INDEX_),routine)
#endif
#if _FABM_DIMENSION_COUNT_>0&&_FABM_VECTORIZED_DIMENSION_INDEX_!=1&&_FABM_DEPTH_DIMENSION_INDEX_!=1
      call check_index(i__,self%domain_size(1),routine,'i')
#endif
#if _FABM_DIMENSION_COUNT_>1&&_FABM_VECTORIZED_DIMENSION_INDEX_!=2&&_FABM_DEPTH_DIMENSION_INDEX_!=2
      call check_index(j__,self%domain_size(2),routine,'j')
#endif
#if _FABM_DIMENSION_COUNT_>2&&_FABM_VECTORIZED_DIMENSION_INDEX_!=3&&_FABM_DEPTH_DIMENSION_INDEX_!=3
      call check_index(k__,self%domain_size(3),routine,'k')
#endif
   end subroutine check_horizontal_location

   subroutine check_vertical_location(self _ARGUMENTS_VERTICAL_IN_,routine)
      class (type_model),intent(in) :: self
      _DECLARE_ARGUMENTS_VERTICAL_IN_
      character(len=*), intent(in) :: routine

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      call check_loop(_VERTICAL_START_,_VERTICAL_STOP_,self%domain_size(_FABM_DEPTH_DIMENSION_INDEX_),routine)
#endif
#if _FABM_DIMENSION_COUNT_>0&&_FABM_DEPTH_DIMENSION_INDEX_!=1
      call check_index(i__,self%domain_size(1),routine,'i')
#endif
#if _FABM_DIMENSION_COUNT_>1&&_FABM_DEPTH_DIMENSION_INDEX_!=2
      call check_index(j__,self%domain_size(2),routine,'j')
#endif
#if _FABM_DIMENSION_COUNT_>2&&_FABM_DEPTH_DIMENSION_INDEX_!=3
      call check_index(k__,self%domain_size(3),routine,'k')
#endif
   end subroutine check_vertical_location

   subroutine check_extents_1d(array,required_size1,routine,array_name,shape_description)
      real(rk),        intent(in) :: array(:)
      integer,         intent(in) :: required_size1
      character(len=*),intent(in) :: routine,array_name,shape_description
      character(len=8) :: actual,required
      if (size(array,1)/=required_size1) then
         write (actual,  '(i0)') size(array,1)
         write (required,'(i0)') required_size1
         call fatal_error(routine,'shape of argument '//trim(array_name)//' is ('//trim(actual)//') but should be ('//trim(required)//') = '//trim(shape_description))
      end if
   end subroutine check_extents_1d

   subroutine check_extents_2d(array,required_size1,required_size2,routine,array_name,shape_description)
      real(rk),        intent(in) :: array(:,:)
      integer,         intent(in) :: required_size1,required_size2
      character(len=*),intent(in) :: routine,array_name,shape_description
      character(len=17) :: actual,required
      if (size(array,1)/=required_size1.or.size(array,2)/=required_size2) then
         write (actual,  '(i0,a,i0)') size(array,1),',',size(array,2)
         write (required,'(i0,a,i0)') required_size1,',',required_size2
         call fatal_error(routine,'shape of argument '//trim(array_name)//' is ('//trim(actual)//') but should be ('//trim(required)//') = '//trim(shape_description))
      end if
   end subroutine check_extents_2d

   subroutine check_extents_3d(array,required_size1,required_size2,required_size3,routine,array_name,shape_description)
      real(rk),        intent(in) :: array(:,:,:)
      integer,         intent(in) :: required_size1,required_size2,required_size3
      character(len=*),intent(in) :: routine,array_name,shape_description
      character(len=26) :: actual,required
      if (size(array,1)/=required_size1.or.size(array,2)/=required_size2.or.size(array,3)/=required_size3) then
         write (actual,  '(i0,a,i0,a,i0)') size(array,1),',',size(array,2),',',size(array,3)
         write (required,'(i0,a,i0,a,i0)') required_size1,',',required_size2,',',required_size3
         call fatal_error(routine,'shape of argument '//trim(array_name)//' is ('//trim(actual)//') but should be ('//trim(required)//') = '//trim(shape_description))
      end if
   end subroutine check_extents_3d

   subroutine check_loop(istart,istop,imax,routine)
      integer,         intent(in) :: istart,istop,imax
      character(len=*),intent(in) :: routine

      character(len=8) :: str1,str2

      if (istart<1) then
         write (str1,'(i0)') istart
         call fatal_error(routine,'Loop start index '//trim(str1)//' is non-positive.')
      end if
      if (istop>imax) then
         write (str1,'(i0)') istop
         write (str2,'(i0)') imax
         call fatal_error(routine,'Loop stop index '//trim(str1)//' exceeds size of vectorized dimension ('//trim(str2)//').')
      end if
      if (istart>istop) then
         write (str1,'(i0)') istart
         write (str2,'(i0)') istop
         call fatal_error(routine,'Loop start index '//trim(str1)//' exceeds stop index '//trim(str2)//'.')
      end if
   end subroutine check_loop

   subroutine check_index(i,i_max,routine,name)
      integer,         intent(in) :: i,i_max
      character(len=*),intent(in) :: routine,name
      character(len=8) :: str1,str2
      if (i<1) then
         write (str1,'(i0)') i
         call fatal_error(routine,'Index '//name//' = '//trim(str1)//' is non-positive.')
      end if
      if (i>i_max) then
         write (str1,'(i0)') i
         write (str2,'(i0)') i_max
         call fatal_error(routine,'Index '//name//' = '//trim(str1)//' exceeds size of associated dimension ('//trim(str2)//').')
      end if
   end subroutine check_index

end module fabm

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

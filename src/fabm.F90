#include "fabm_driver.h"
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
                                     type_global_standard_variable, initialize_standard_variables
   use fabm_types
   use fabm_library
   use fabm_expressions
   use fabm_driver
   use fabm_properties
   use fabm_builtin_models
   use fabm_coupling
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
   public fabm_link_bulk_state_data, fabm_link_bottom_state_data, fabm_link_surface_state_data
   public fabm_link_bulk_data, fabm_link_horizontal_data, fabm_link_scalar_data
   public fabm_get_bulk_diagnostic_data, fabm_get_horizontal_diagnostic_data

#ifdef _FABM_MASK_
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
!
! !PUBLIC TYPES:
!
   ! ====================================================================================================
   ! Variable identifiers used by host models.
   ! ====================================================================================================

   type type_bulk_variable_id
      type (type_internal_variable),pointer :: variable => null()
      integer                            :: state_index = -1
      integer                            :: read_index = -1
   end type

   type type_horizontal_variable_id
      type (type_internal_variable),pointer :: variable => null()
      integer                                  :: state_index = -1
      integer                                  :: read_index = -1
   end type

   type type_scalar_variable_id
      type (type_internal_variable),pointer :: variable => null()
      integer                              :: read_index = -1
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
      type (type_bulk_variable_id)       :: globalid
      integer,allocatable, dimension(:)  :: sms_indices, surface_flux_indices, bottom_flux_indices
   end type type_state_variable_info

   type,extends(type_external_variable) :: type_horizontal_state_variable_info
      type (type_horizontal_standard_variable) :: standard_variable
      real(rk)                                 :: initial_value = 0.0_rk
      type (type_horizontal_variable_id)       :: globalid
      integer,allocatable, dimension(:)        :: sms_indices
   end type type_horizontal_state_variable_info

!  Derived type describing a diagnostic variable
   type,extends(type_external_variable) :: type_diagnostic_variable_info
      type (type_bulk_standard_variable) :: standard_variable
      integer                            :: time_treatment = time_treatment_last ! See time_treatment_* parameters above
      logical                            :: save           = .false.
      integer                            :: save_index     = 0
      integer                            :: source
   end type type_diagnostic_variable_info

   type,extends(type_external_variable) :: type_horizontal_diagnostic_variable_info
      type (type_horizontal_standard_variable) :: standard_variable
      integer                                  :: time_treatment = time_treatment_last ! See time_treatment_* parameters above
      logical                                  :: save           = .false.
      integer                                  :: save_index     = 0
      integer                                  :: source
   end type type_horizontal_diagnostic_variable_info

!  Derived type describing a conserved quantity
   type,extends(type_external_variable) :: type_conserved_quantity_info
      type (type_bulk_standard_variable)            :: standard_variable
      type (type_aggregate_variable),pointer        :: aggregate_variable
      integer                                       :: index              = -1
      integer                                       :: horizontal_index   = -1
      class (type_horizontal_weighted_sum),pointer  :: horizontal_sum     => null()
   end type type_conserved_quantity_info

   ! Derived type for a single generic biogeochemical model
   type type_model
      type (type_base_model)  :: root
      type (type_model_list)  :: models
      type (type_environment) :: environment

      logical :: initialized = .false.

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

      integer                       :: extinction_index = -1
      type (type_model_list)        :: extinction_call_list, conserved_quantity_call_list

      type (type_link_list) :: links_postcoupling

      ! Declare the arrays for diagnostic variable values.
      real(rk),allocatable _DIMENSION_GLOBAL_PLUS_1_            :: diag
      real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: diag_hz
      real(rk),allocatable _DIMENSION_GLOBAL_                   :: zero
      real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_        :: zero_hz
      integer                                                   :: domain_size(_FABM_DIMENSION_COUNT_)
      integer                                                   :: horizontal_domain_size(_FABM_DIMENSION_COUNT_HZ_)
      integer,allocatable,dimension(:)                          :: zero_bulk_indices,zero_surface_indices,zero_bottom_indices
   contains
      procedure :: link_bulk_data_by_variable => fabm_link_bulk_data_by_variable
      procedure :: link_bulk_data_by_id   => fabm_link_bulk_data_by_id
      procedure :: link_bulk_data_by_sn   => fabm_link_bulk_data_by_sn
      procedure :: link_bulk_data_by_name => fabm_link_bulk_data_by_name
      generic :: link_bulk_data => link_bulk_data_by_variable,link_bulk_data_by_id,link_bulk_data_by_sn,link_bulk_data_by_name

      procedure :: link_horizontal_data_by_variable => fabm_link_horizontal_data_by_variable
      procedure :: link_horizontal_data_by_id       => fabm_link_horizontal_data_by_id
      procedure :: link_horizontal_data_by_sn       => fabm_link_horizontal_data_by_sn
      procedure :: link_horizontal_data_by_name     => fabm_link_horizontal_data_by_name
      generic :: link_horizontal_data => link_horizontal_data_by_variable,link_horizontal_data_by_id,link_horizontal_data_by_sn,link_horizontal_data_by_name

      procedure :: link_scalar_by_id   => fabm_link_scalar_by_id
      procedure :: link_scalar_by_sn   => fabm_link_scalar_by_sn
      procedure :: link_scalar_by_name => fabm_link_scalar_by_name
      generic :: link_scalar => link_scalar_by_id,link_scalar_by_sn,link_scalar_by_name

      procedure :: get_bulk_variable_id_by_name => fabm_get_bulk_variable_id_by_name
      procedure :: get_bulk_variable_id_sn => fabm_get_bulk_variable_id_sn
      generic :: get_bulk_variable_id => get_bulk_variable_id_by_name, get_bulk_variable_id_sn

      procedure :: get_horizontal_variable_id_by_name => fabm_get_horizontal_variable_id_by_name
      procedure :: get_horizontal_variable_id_sn => fabm_get_horizontal_variable_id_sn
      generic :: get_horizontal_variable_id => get_horizontal_variable_id_by_name, get_horizontal_variable_id_sn

      procedure :: get_scalar_variable_id_by_name => fabm_get_scalar_variable_id_by_name
      procedure :: get_scalar_variable_id_sn => fabm_get_scalar_variable_id_sn
      generic :: get_scalar_variable_id => get_scalar_variable_id_by_name, get_scalar_variable_id_sn
   end type type_model

   type,extends(type_base_model) :: type_custom_extinction_calculator
      type (type_diagnostic_variable_id) :: id_output
      type (type_model_list)  :: models
   contains
      procedure :: initialize => custom_extinction_calculator_initialize
      procedure :: do => custom_extinction_calculator_do
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
      module procedure fabm_link_bulk_data_by_id
      module procedure fabm_link_horizontal_data_by_id
      module procedure fabm_link_scalar_by_id
      module procedure fabm_link_bulk_data_by_sn
      module procedure fabm_link_horizontal_data_by_sn
      module procedure fabm_link_scalar_by_sn
   end interface

   ! Subroutine for providing FABM with variable data on the full spatial domain.
   interface fabm_link_bulk_data
      module procedure fabm_link_bulk_data_by_variable
      module procedure fabm_link_bulk_data_by_id
      module procedure fabm_link_bulk_data_by_sn
      module procedure fabm_link_bulk_data_by_name
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

   interface fabm_get_variable_name
      module procedure fabm_get_bulk_variable_name
      module procedure fabm_get_horizontal_variable_name
      module procedure fabm_get_scalar_variable_name
   end interface

   interface fabm_is_variable_used
      module procedure fabm_is_bulk_variable_used
      module procedure fabm_is_horizontal_variable_used
      module procedure fabm_is_scalar_variable_used
   end interface

   interface fabm_variable_needs_values
      module procedure fabm_bulk_variable_needs_values
   end interface

   ! For backward compatibility only:
   interface fabm_do_benthos
      module procedure fabm_do_bottom_rhs
      module procedure fabm_do_bottom_ppdd
   end interface
   interface fabm_get_surface_exchange
      module procedure fabm_do_surface
   end interface
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------

   contains

   subroutine fabm_initialize_library()
      logical, save :: initialized = .false.

      ! Do nothing if already initialized.
      if (initialized) return

      ! If needed, create default object for communication (e.g., logging, error reporting) with host.
      if (.not.associated(driver)) allocate(type_base_driver::driver)

      ! Create all standard variable objects.
      call initialize_standard_variables()

      ! Create the model factory.
      call fabm_create_model_factory()

      initialized = .true.
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
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman

   logical                   :: isopen,initialize
   character(len=256)        :: file_eff
   integer                   :: i,j,modelcount,ownindex
   character(len=64)         :: models(256),instancename
   class (type_base_model),pointer :: childmodel
   logical,parameter         :: alwayspostfixindex=.false.
   namelist /fabm_nml/ models
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

         call log_message('Initializing biogeochemical model "'//trim(instancename)//'"...')
         call model%root%add_child(childmodel,instancename,configunit=file_unit)
         call log_message('model "'//trim(instancename)//'" initialized successfully.')
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
      type (type_aggregate_variable),pointer :: aggregate_variable
      class (type_custom_extinction_calculator),pointer :: extinction_calculator
!EOP
!-----------------------------------------------------------------------
!BOC
      self%info => self ! For backward compatibility (pre 11/2013 hosts only)

      ! Make sure a variable for light extinction is created at the root level when calling freeze_model_info.
      ! This variable is used from fabm_get_light_extinction.
      aggregate_variable => get_aggregate_variable(self%root, &
         standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      aggregate_variable%bulk_required = .true.

      ! Create placeholder variables for zero fields.
      ! Values for these fields will only be provided if actually used by one of the biogeochemical models.
      call self%root%add_bulk_variable('zero')
      call self%root%add_horizontal_variable('zero_hz')

      allocate(extinction_calculator)
      call self%root%add_child(extinction_calculator,'custom_extinction_calculator',configunit=-1)

      ! Filter out expressions that FABM can handle itself.
      ! The remainder, if any, must be handled by the host model.
      call filter_expressions(self)

      ! This will resolve all FABM dependencies and generate final authorative lists of variables of different types.
      call freeze_model_info(self%root)

      ! Build final authorative arrays with variable metadata .
      call classify_variables(self)

      call extinction_calculator%models%extend(self%models)
      call add_to_model_call_list(self,'attenuation_coefficient_of_photosynthetic_radiative_flux',self%extinction_call_list)

      self%initialized = .true.

   end subroutine fabm_initialize
!EOC

   subroutine filter_expressions(self)
      class (type_model),intent(inout)           :: self
      class (type_expression),           pointer :: current,previous,next
      class (type_simple_depth_integral),pointer :: integral
      logical                                    :: filter

      nullify(previous)
      current => self%root%first_expression
      do while (associated(current))
         filter = .false.
         select type (current)
            class is (type_vertical_integral)
#ifndef _FABM_DEPTH_DIMENSION_INDEX_
               ! For models without depth dimension, FABM can calculate depth averages and integrals itself.
               allocate(integral)
               integral%minimum_depth = current%minimum_depth
               integral%maximum_depth = current%maximum_depth
               integral%average       = current%average
               call self%root%add_child(integral,trim(current%output_name)//'_calculator',configunit=-1)
               call integral%request_coupling(integral%id_input,current%input_name)
               call self%root%request_coupling(current%output_name,integral%id_output%link%target%name)
               filter = .true.
#endif
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
      type (type_model_list_node),pointer :: child
      character(len=4) :: strcount
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   nullify(self%info)
   self%initialized = .false.

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
   subroutine fabm_set_domain(self _ARG_LOCATION_,seconds_per_time_unit)
!
! !INPUT PARAMETERS:
   class (type_model),target,intent(inout) :: self
   _DECLARE_LOCATION_ARG_
   real(rk),optional,        intent(in)    :: seconds_per_time_unit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  type (type_model_list_node), pointer :: node
  integer                              :: ivar,n,index
  class (type_expression),pointer      :: expression
  type (type_link),pointer             :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if _FABM_DIMENSION_COUNT_>0
   self%domain_size = (/ _LOCATION_ /)
#endif
#if _FABM_DIMENSION_COUNT_HZ_>0
   self%horizontal_domain_size = (/ _LOCATION_HZ_ /)
#endif

   ! Forward domain to individual biogeochemical models.
   ! Also determine whether one of the models needs conservation checks for its sink and source terms.
   node => self%models%first
   do while (associated(node))
      call node%model%set_domain(_LOCATION_)
      node => node%next
   end do

   ! Allocate memory to hold readable variables for a single domain slice.
   allocate(self%environment%prefetch _INDEX_SLICE_PLUS_1_(size(self%environment%data)))
   allocate(self%environment%prefetch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(size(self%environment%data_hz)))
   allocate(self%environment%prefetch_scalar(size(self%environment%data_scalar)))
   self%environment%prefetch = 0.0_rk
   self%environment%prefetch_hz = 0.0_rk
   self%environment%prefetch_scalar = 0.0_rk

   allocate(self%zero _INDEX_LOCATION_)
   self%zero = 0.0_rk
   call self%link_bulk_data('zero',self%zero)

   allocate(self%zero_hz _INDEX_HORIZONTAL_LOCATION_)
   self%zero_hz = 0.0_rk
   call self%link_horizontal_data('zero_hz',self%zero_hz)

   ! Merge diagnostic variables that contribute to the same sink/source terms,
   ! surface fluxes, bottom fluxes, provided no other variables depend on them.
   link => self%links_postcoupling%first
   do while (associated(link))
      call merge_aggregating_diagnostics(link%target%sms_list)
      call merge_aggregating_diagnostics(link%target%surface_flux_list)
      call merge_aggregating_diagnostics(link%target%bottom_flux_list)
      link => link%next
   end do

   ! Create arrays with diagnostic indices contributing to sink/source terms,
   ! surface fluxes, bottom fluxes. Must be done after calls to merge_aggregating_diagnostics.
   do ivar=1,size(self%state_variables)
      call copy_write_indices(self%state_variables(ivar)%target%sms_list,         self%state_variables(ivar)%sms_indices)
      call copy_write_indices(self%state_variables(ivar)%target%surface_flux_list,self%state_variables(ivar)%surface_flux_indices)
      call copy_write_indices(self%state_variables(ivar)%target%bottom_flux_list, self%state_variables(ivar)%bottom_flux_indices)
   end do
   do ivar=1,size(self%surface_state_variables)
      call copy_write_indices(self%surface_state_variables(ivar)%target%sms_list, self%surface_state_variables(ivar)%sms_indices)
   end do
   do ivar=1,size(self%bottom_state_variables)
      call copy_write_indices(self%bottom_state_variables(ivar)%target%sms_list, self%bottom_state_variables(ivar)%sms_indices)
   end do

   ! Assign write indices in scratch space to all bulk diagnostic variables.
   ! Must be done after calls to merge_aggregating_diagnostics.
   n = 0
   link => self%links_postcoupling%first
   do while (associated(link))
      select case (link%target%domain)
         case (domain_bulk)
            if (.not.link%target%write_indices%is_empty().and.link%target%write_indices%value==-1) then
               n = n + 1
               call link%target%write_indices%set_value(n)
            end if
      end select
      link => link%next
   end do

   ! Allocate scratch memory for bulk variables.
   allocate(self%environment%scratch _INDEX_SLICE_PLUS_1_(n))

   ! Assign write indices in scratch space to all bulk diagnostic variables.
   ! Must be done after calls to merge_aggregating_diagnostics.
   n = 0
   link => self%links_postcoupling%first
   do while (associated(link))
      select case (link%target%domain)
         case (domain_horizontal,domain_surface,domain_bottom)
            if (.not.link%target%write_indices%is_empty().and.link%target%write_indices%value==-1) then
               n = n + 1
               call link%target%write_indices%set_value(n)
            end if
      end select
      link => link%next
   end do

   ! Allocate scratch memory for horizontal variables.
   allocate(self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(n))

   ! Create lists of scratch variables that require zeroing before calling biogeochemical models.
   link => self%links_postcoupling%first
   do while (associated(link))
      select case (link%target%domain)
         case (domain_bulk)
            call copy_write_indices(link%target%sms_list,         self%zero_bulk_indices)
            call copy_write_indices(link%target%bottom_flux_list, self%zero_bottom_indices)
            call copy_write_indices(link%target%surface_flux_list,self%zero_surface_indices)
         case (domain_bottom)
            call copy_write_indices(link%target%sms_list,self%zero_bottom_indices)
         case (domain_surface)
            call copy_write_indices(link%target%sms_list,self%zero_surface_indices)
      end select
      link => link%next
   end do

#ifdef _FABM_MASK_
   allocate(self%environment%prefetch_mask _INDEX_SLICE_)
#endif

#ifdef _FULL_DOMAIN_IN_SLICE_
   ! Fill memory for diagnostic variable with missing values, and save pointers for data retrieval.
   do ivar=1,size(self%diagnostic_variables)
      index = self%diagnostic_variables(ivar)%target%write_indices%value
      if (index>0) then
         self%environment%scratch(_PREARG_LOCATION_DIMENSIONS_ index) = self%diagnostic_variables(ivar)%missing_value
         call fabm_link_bulk_data(self,self%diagnostic_variables(ivar)%target, &
            self%environment%scratch(_PREARG_LOCATION_DIMENSIONS_ index))
      end if
   end do
   do ivar=1,size(self%horizontal_diagnostic_variables)
      index = self%horizontal_diagnostic_variables(ivar)%target%write_indices%value
      if (index>0) then
         self%environment%scratch_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ index) = self%horizontal_diagnostic_variables(ivar)%missing_value
         call fabm_link_horizontal_data(self,self%horizontal_diagnostic_variables(ivar)%target, &
            self%environment%scratch_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ index))
      end if
   end do
#else
   ! Allocate memory for full-domain storage of bulk diagnostics
   n = 0
   do ivar=1,size(self%diagnostic_variables)
      if (self%diagnostic_variables(ivar)%save.or.self%diagnostic_variables(ivar)%target%read_indices%value/=-1) then
         n = n + 1
         self%diagnostic_variables(ivar)%save_index = n
      end if
   end do
   allocate(self%diag(_PREARG_LOCATION_ 0:n))
   self%diag = 0.0_rk

   ! Allocate memory for full-domain storage of horizontal diagnostics
   n = 0
   do ivar=1,size(self%horizontal_diagnostic_variables)
      if (self%horizontal_diagnostic_variables(ivar)%save.or.self%horizontal_diagnostic_variables(ivar)%target%read_indices%value/=-1) then
         n = n + 1
         self%horizontal_diagnostic_variables(ivar)%save_index = n
      end if
   end do
   allocate(self%diag_hz(_PREARG_LOCATION_HZ_ 0:n))
   self%diag_hz = 0.0_rk

   ! Initialize diagnostic variables to missing value.
   do ivar=1,size(self%diagnostic_variables)
      index = self%diagnostic_variables(ivar)%save_index
      if (index>0) then
         self%diag(_PREARG_LOCATION_DIMENSIONS_ index) = self%diagnostic_variables(ivar)%missing_value
         call fabm_link_bulk_data(self,self%diagnostic_variables(ivar)%target, self%diag(_PREARG_LOCATION_DIMENSIONS_ index))
      end if
   end do
   do ivar=1,size(self%horizontal_diagnostic_variables)
      index = self%horizontal_diagnostic_variables(ivar)%save_index
      if (index>0) then
         self%diag_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ index) = self%horizontal_diagnostic_variables(ivar)%missing_value
         call fabm_link_horizontal_data(self,self%horizontal_diagnostic_variables(ivar)%target, self%diag_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ index))
      end if
   end do
#endif

   if (present(seconds_per_time_unit)) then
      expression => self%root%first_expression
      do while (associated(expression))
         select type (expression)
            class is (type_bulk_temporal_mean)
               expression%period = expression%period/seconds_per_time_unit
               allocate(expression%history(_PREARG_LOCATION_ expression%n+3))
               expression%history = 0.0_rk
               call fabm_link_bulk_data(self,expression%output_name, &
                                        expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+3))
            class is (type_horizontal_temporal_mean)
               expression%period = expression%period/seconds_per_time_unit
               allocate(expression%history(_PREARG_LOCATION_HZ_ expression%n+3))
               expression%history = 0.0_rk
               call fabm_link_horizontal_data(self,expression%output_name, &
                                              expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+3))
         end select
         expression => expression%next
      end do
   end if

   end subroutine fabm_set_domain
!EOC

   subroutine merge_aggregating_diagnostics(list)
      type (type_link_list),intent(inout) :: list

      type (type_link),             pointer :: link
      type (type_internal_variable),pointer :: free_target

      nullify(free_target)
      link => list%first
      do while (associated(link))
         if (link%target%read_indices%is_empty()) then
            ! This diagnostic is only used as increment of source-sink terms, surface flux or bulk flux.
            ! It can be merged.
            if (associated(free_target)) then
               ! We already found a previous variable that can be merged - merge current and previous.
               call free_target%write_indices%extend(link%target%write_indices)
               call link%target%write_indices%clear()
            else
               ! This is the first variable that can be merged. Record it as such and move on.
               free_target => link%target
            end if
         end if
         link => link%next
      end do
   end subroutine merge_aggregating_diagnostics

#ifdef _FABM_MASK_
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with the spatial mask.
!
! !INTERFACE:
   subroutine fabm_set_mask(self, mask)
!
! !INPUT PARAMETERS:
   class (type_model),target,intent(inout)                    :: self
   _FABM_MASK_TYPE_, target, intent(in)    _DIMENSION_GLOBAL_ :: mask
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   self%environment%mask => mask

   end subroutine fabm_set_mask
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  type (type_link),pointer :: link
  logical                  :: ready
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ready = .true.

#ifdef _FABM_MASK_
   if (.not.associated(self%environment%mask)) call log_message('spatial mask has not been set.')
#endif

   link => self%links_postcoupling%first
   do while (associated(link))
      if (.not.link%target%read_indices%is_empty().and..not.link%target%presence==presence_external_optional) then
         select case (link%target%domain)
            case (domain_bulk)
               if (.not.associated(self%environment%data(link%target%read_indices%value)%p)) then
                  call log_message('data for dependency "'//trim(link%name)// &
                     & '", defined on the full model domain, have not been provided.')
                  ready = .false.
               end if
            case (domain_horizontal,domain_surface,domain_bottom)
               if (.not.associated(self%environment%data_hz(link%target%read_indices%value)%p)) then
                  call log_message('data for dependency "'//trim(link%name)// &
                     &  '", defined on a horizontal slice of the model domain, have not been provided.')
                  ready = .false.
               end if
            case (domain_scalar)
               if (.not.associated(self%environment%data_scalar(link%target%read_indices%value)%p)) then
                  call log_message('data for dependency "'//trim(link%name)// &
                     &  '", defined as global scalar quantity, have not been provided.')
                  ready = .false.
               end if
         end select
      end if
      link => link%next
   end do

   if (.not.ready) call fatal_error('fabm_check_ready','FABM is lacking required data.')

   ! Host has sent all fields - disbale all optional fields that were not provided in individual BGC models.
   call filter_readable_variable_registry(self)

   end subroutine fabm_check_ready
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the model state (pelagic).
!
! !INTERFACE:
   subroutine fabm_initialize_state(self _ARG_LOCATION_ND_)
!
! !INPUT PARAMETERS:
   class (type_model),     intent(inout) :: self
   _DECLARE_LOCATION_ARG_ND_
!
! !LOCAL PARAMETERS:
   integer                               :: ivar
   type (type_model_list_node), pointer  :: node
   real(rk)                              :: initial_value
   type (type_bulk_data_pointer)         :: p
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   call prefetch(self%environment _ARG_LOCATION_ND_)

   do ivar=1,size(self%state_variables)
      ! Shortcuts to variable information - this demonstrably helps the compiler with vectorization (ifort).
      p = self%environment%data(self%state_variables(ivar)%globalid%read_index)
      initial_value = self%state_variables(ivar)%initial_value
      _CONCURRENT_LOOP_BEGIN_EX_(self%environment)
         p%p _INDEX_LOCATION_ = initial_value
      _CONCURRENT_LOOP_END_
   end do

   ! Allow biogeochemical models to initialize their bulk state.
   node => self%models%first
   do while (associated(node))
      call node%model%initialize_state(_ARGUMENTS_ND_IN_)
      node => node%next
   end do

   end subroutine fabm_initialize_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the bottom model state.
!
! !INTERFACE:
   subroutine fabm_initialize_bottom_state(self _ARG_LOCATION_VARS_HZ_)
!
! !INPUT PARAMETERS:
   class (type_model), intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
!
! !LOCAL PARAMETERS:
   integer                               :: ivar
   type (type_model_list_node), pointer  :: node
   real(rk)                              :: initial_value
   type (type_horizontal_data_pointer)   :: p
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   ! Initialize bottom variables
   do ivar=1,size(self%bottom_state_variables)
      p = self%environment%data_hz(self%bottom_state_variables(ivar)%globalid%read_index)
      initial_value = self%bottom_state_variables(ivar)%initial_value
      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
         p%p _INDEX_HORIZONTAL_LOCATION_ = initial_value
      _CONCURRENT_HORIZONTAL_LOOP_END_
   end do

   ! Allow biogeochemical models to initialize their bottom state.
   node => self%models%first
   do while (associated(node))
      call node%model%initialize_bottom_state(_ARGUMENTS_IN_HZ_)
      node => node%next
   end do

   end subroutine fabm_initialize_bottom_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the surface model state.
!
! !INTERFACE:
   subroutine fabm_initialize_surface_state(self _ARG_LOCATION_VARS_HZ_)
!
! !INPUT PARAMETERS:
   class (type_model),     intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
!
! !LOCAL PARAMETERS:
   integer                               :: ivar
   type (type_model_list_node), pointer  :: node
   real(rk)                              :: initial_value
   type (type_horizontal_data_pointer)   :: p
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   ! Initialize surface variables
   do ivar=1,size(self%surface_state_variables)
      p = self%environment%data_hz(self%surface_state_variables(ivar)%globalid%read_index)
      initial_value = self%surface_state_variables(ivar)%initial_value
      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
         p%p _INDEX_HORIZONTAL_LOCATION_ = initial_value
      _CONCURRENT_HORIZONTAL_LOOP_END_
   end do

   ! Allow biogeochemical models to initialize their surface state.
   node => self%models%first
   do while (associated(node))
      call node%model%initialize_surface_state(_ARGUMENTS_IN_HZ_)
      node => node%next
   end do

   end subroutine fabm_initialize_surface_state
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
   type (type_link),       pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_bulk) then
         if (link%target%name==name.or.get_safe_name(link%target%name)==name) then
            id = create_external_bulk_id(link%target)
            return
         end if
      end if
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_bulk.and.associated(link%target%standard_variable)) then
         if (link%target%standard_variable%name==name) then
            id = create_external_bulk_id(link%target)
            return
         end if
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
!EOP
!-----------------------------------------------------------------------
!BOC
   id = create_external_bulk_id_for_standard_name(self,standard_variable)

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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
   type (type_link),       pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_horizontal.or.link%target%domain==domain_surface.or.link%target%domain==domain_bottom) then
         if (link%target%name==name.or.get_safe_name(link%target%name)==name) then
            id = create_external_horizontal_id(link%target)
            return
         end if
      end if
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%links%first
   do while (associated(link))
      if ((link%target%domain==domain_horizontal.or.link%target%domain==domain_surface.or.link%target%domain==domain_bottom) &
          .and.associated(link%target%standard_variable)) then
         if (link%target%standard_variable%name==name) then
            id = create_external_horizontal_id(link%target)
            return
         end if
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
!EOP
!-----------------------------------------------------------------------
!BOC
   id = create_external_horizontal_id_for_standard_name(self,standard_variable)

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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
   type (type_link),       pointer :: link
!
!EOP
!-----------------------------------------------------------------------
!BOC
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_scalar) then
         if (link%target%name==name.or.get_safe_name(link%target%name)==name) then
            id = create_external_scalar_id(link%target)
            return
         end if
      end if
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%links%first
   do while (associated(link))
      if (link%target%domain==domain_scalar.and.associated(link%target%standard_variable)) then
         if (link%target%standard_variable%name==name) then
            id = create_external_scalar_id(link%target)
            return
         end if
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
!EOP
!-----------------------------------------------------------------------
!BOC
   id = create_external_scalar_id_for_standard_name(self,standard_variable)

   end function fabm_get_scalar_variable_id_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable name for the given variable
! identifier.
!
! !INTERFACE:
   function fabm_get_bulk_variable_name(model,id) result(name)
!
! !INPUT PARAMETERS:
   class (type_model),            intent(in)  :: model
   type(type_bulk_variable_id),   intent(in)  :: id
!
! !RETURN VALUE:
   character(len=attribute_length)            :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   name = ''
   if (associated(id%variable)) name = get_safe_name(id%variable%name)

   end function fabm_get_bulk_variable_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable name for the given variable
! identifier.
!
! !INTERFACE:
   function fabm_get_horizontal_variable_name(model,id) result(name)
!
! !INPUT PARAMETERS:
   class (type_model),               intent(in) :: model
   type(type_horizontal_variable_id),intent(in) :: id
!
! !RETURN VALUE:
   character(len=attribute_length)              :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   name = ''
   if (associated(id%variable)) name = get_safe_name(id%variable%name)

   end function fabm_get_horizontal_variable_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable name for the given variable
! identifier.
!
! !INTERFACE:
   function fabm_get_scalar_variable_name(model,id) result(name)
!
! !INPUT PARAMETERS:
   class (type_model),            intent(in) :: model
   type(type_scalar_variable_id), intent(in) :: id
!
! !RETURN VALUE:
   character(len=attribute_length)           :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   name = ''
   if (associated(id%variable)) name = get_safe_name(id%variable%name)

   end function fabm_get_scalar_variable_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Determine whether a bulk variable is used [required] by biogeochemical models running in FABM.
!
! !INTERFACE:
   function fabm_is_bulk_variable_used(id) result(used)
!
! !INPUT PARAMETERS:
   type(type_bulk_variable_id),   intent(in)  :: id
!
! !RETURN VALUE:
   logical                                    :: used
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   used = id%read_index/=-1

   end function fabm_is_bulk_variable_used
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Determine whether a bulk variable is available to biogeochemical models running in FABM.
!
! !INTERFACE:
   function fabm_bulk_variable_needs_values(self,id) result(required)
!
! !INPUT PARAMETERS:
   class (type_model),         intent(inout) :: self
   type(type_bulk_variable_id),intent(in)    :: id
!
! !RETURN VALUE:
   logical                                :: required
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   required = id%read_index/=-1
   if (required) required = .not.associated(self%environment%data(id%read_index)%p)

   end function fabm_bulk_variable_needs_values
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable name for the given variable
! identifier.
!
! !INTERFACE:
   function fabm_is_horizontal_variable_used(id) result(used)
!
! !INPUT PARAMETERS:
   type(type_horizontal_variable_id),intent(in):: id
!
! !RETURN VALUE:
   logical                                     :: used
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   used = id%read_index/=-1

   end function fabm_is_horizontal_variable_used
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable name for the given variable
! identifier.
!
! !INTERFACE:
   function fabm_is_scalar_variable_used(id) result(used)
!
! !INPUT PARAMETERS:
   type(type_scalar_variable_id), intent(in)  :: id
!
! !RETURN VALUE:
   logical                                    :: used
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   used = id%read_index/=-1

   end function fabm_is_scalar_variable_used
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by an internal variable object.
!
! !INTERFACE:
   subroutine fabm_link_bulk_data_by_variable(self,variable,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: self
   type(type_internal_variable),      intent(in)    :: variable
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&_FABM_DIMENSION_COUNT_>0
   do i=1,size(self%domain_size)
      if (size(dat,i)/=self%domain_size(i)) then
         call fatal_error('fabm_link_bulk_data_by_variable','dimensions of FABM domain and provided array do not match for variable '//trim(variable%name)//'.')
      end if
   end do
#endif

   if (variable%read_indices%value/=-1) self%environment%data(variable%read_indices%value)%p => dat

   end subroutine fabm_link_bulk_data_by_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by an external identifier.
!
! !INTERFACE:
   subroutine fabm_link_bulk_data_by_id(self,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: self
   type(type_bulk_variable_id),       intent(inout) :: id
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                                :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&_FABM_DIMENSION_COUNT_>0
   do i=1,size(self%domain_size)
      if (size(dat,i)/=self%domain_size(i)) then
         call fatal_error('fabm_link_bulk_data','dimensions of FABM domain and provided array do not match for variable '//trim(id%variable%name)//'.')
      end if
   end do
#endif

   if (id%read_index/=-1) self%environment%data(id%read_index)%p => dat

   end subroutine fabm_link_bulk_data_by_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by its standard name.
!
! !INTERFACE:
   subroutine fabm_link_bulk_data_by_sn(model,standard_variable,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: model
   type(type_bulk_standard_variable), intent(in)    :: standard_variable
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
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
   if (fabm_is_variable_used(id)) call fabm_link_bulk_data(model,id,dat)

   end subroutine fabm_link_bulk_data_by_sn
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by its name.
!
! !INTERFACE:
   subroutine fabm_link_bulk_data_by_name(model,name,dat)
!
! !INPUT PARAMETERS:
   class (type_model),target,         intent(inout) :: model
   character(len=*),                  intent(in)    :: name
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   type (type_bulk_variable_id)                             :: id
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_bulk_variable_id(model,name)

   ! Only link the data if needed (if the variable identifier is valid).
   if (fabm_is_variable_used(id)) call fabm_link_bulk_data(model,id,dat)

   end subroutine fabm_link_bulk_data_by_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by an internal variable object.
!
! !INTERFACE:
   subroutine fabm_link_horizontal_data_by_variable(self,variable,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: self
   type (type_internal_variable),                intent(inout) :: variable
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                                :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&_FABM_DIMENSION_COUNT_HZ_>0
   do i=1,size(self%horizontal_domain_size)
      if (size(dat,i)/=self%horizontal_domain_size(i)) then
         call fatal_error('fabm_link_horizontal_data','dimensions of FABM domain and provided array do not match for variable '//trim(variable%name)//'.')
      end if
   end do
#endif

   if (variable%read_indices%value/=-1) self%environment%data_hz(variable%read_indices%value)%p => dat

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
   subroutine fabm_link_horizontal_data_by_id(self,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: self
   type (type_horizontal_variable_id),           intent(inout) :: id
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                                :: i
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&_FABM_DIMENSION_COUNT_HZ_>0
   do i=1,size(self%horizontal_domain_size)
      if (size(dat,i)/=self%horizontal_domain_size(i)) then
         call fatal_error('fabm_link_horizontal_data','dimensions of FABM domain and provided array do not match for variable '//trim(id%variable%name)//'.')
      end if
   end do
#endif

   if (id%read_index/=-1) self%environment%data_hz(id%read_index)%p => dat

   end subroutine fabm_link_horizontal_data_by_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by its standard name.
!
! !INTERFACE:
   subroutine fabm_link_horizontal_data_by_sn(model,standard_variable,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: model
   type(type_horizontal_standard_variable),      intent(in)    :: standard_variable
   real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   type(type_horizontal_variable_id) :: id
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
   subroutine fabm_link_scalar_by_id(self,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),            intent(inout) :: self
   type (type_scalar_variable_id),intent(inout) :: id
   real(rk),target,               intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   if (id%read_index/=-1) self%environment%data_scalar(id%read_index)%p => dat

   end subroutine fabm_link_scalar_by_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by its standard name.
!
! !INTERFACE:
   subroutine fabm_link_scalar_by_sn(model,standard_variable,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                 intent(inout) :: model
   type(type_global_standard_variable),intent(in)    :: standard_variable
   real(rk),target,                    intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   type(type_scalar_variable_id) :: id
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
   subroutine fabm_link_bulk_state_data(self,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: self
   integer,                           intent(in)    :: id
   real(rk) _DIMENSION_GLOBAL_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_link_bulk_data(self,self%state_variables(id)%globalid,dat)

   end subroutine fabm_link_bulk_state_data
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_link_horizontal_data(self,self%bottom_state_variables(id)%globalid,dat)

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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_link_horizontal_data(self,self%surface_state_variables(id)%globalid,dat)

   end subroutine fabm_link_surface_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Returns (a pointer to) the array with data for
! a single diagnostic variable, defined on the full spatial domain.
!
! !INTERFACE:
   function fabm_get_bulk_diagnostic_data(self,id) result(dat)
!
! !INPUT PARAMETERS:
   class (type_model),target,         intent(in) :: self
   integer,                           intent(in) :: id
   real(rk) _DIMENSION_GLOBAL_,pointer           :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef _FULL_DOMAIN_IN_SLICE_
   ! Retrieve a pointer to the array holding the requested data.
   dat => self%environment%scratch(_PREARG_LOCATION_DIMENSIONS_ self%diagnostic_variables(id)%target%write_indices%value)
#else
   ! Retrieve a pointer to the array holding the requested data.
   dat => self%diag(_PREARG_LOCATION_DIMENSIONS_ self%diagnostic_variables(id)%save_index)
#endif

   end function fabm_get_bulk_diagnostic_data
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
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef _FULL_DOMAIN_IN_SLICE_
   ! Retrieve a pointer to the array holding the requested data.
   dat => self%environment%scratch_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ self%horizontal_diagnostic_variables(id)%target%write_indices%value)
#else
   ! Retrieve a pointer to the array holding the requested data.
   dat => self%diag_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ self%horizontal_diagnostic_variables(id)%save_index)
#endif

   end function fabm_get_horizontal_diagnostic_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the biogeochemical
! model tree.
!
! !INTERFACE:
   subroutine fabm_do_rhs(self _ARG_LOCATION_ND_,dy)
!
! !INPUT PARAMETERS:
   class (type_model),           intent(inout) :: self
   _DECLARE_LOCATION_ARG_ND_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_SLICE_PLUS_1_, intent(inout) :: dy
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
   integer :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&defined(_FABM_USE_1D_LOOP_)
   if (size(dy,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_do_rhs','Size of first dimension of dy should match loop length.')
#endif
   call prefetch(self%environment _ARG_LOCATION_ND_)

   ! Initialize scratch variables that will hold source-sink terms to zero,
   ! because they will be incremented.
   do i=1,size(self%zero_bulk_indices)
      j = self%zero_bulk_indices(i)
      _CONCURRENT_LOOP_BEGIN_EX_(self%environment)
         self%environment%scratch _INDEX_SLICE_PLUS_1_(j) = 0.0_rk
      _CONCURRENT_LOOP_END_
   end do

   node => self%models%first
   do while (associated(node))
      call node%model%do(_ARGUMENTS_ND_IN_)

      ! Copy newly written diagnostics to prefetch
      do i=1,size(node%model%reused_diag)
         if (node%model%reused_diag(i)%source==source_do) then
            j = node%model%reused_diag(i)%read_index
            k = node%model%reused_diag(i)%write_index
            _CONCURRENT_LOOP_BEGIN_EX_(self%environment)
               self%environment%prefetch _INDEX_SLICE_PLUS_1_(j) = self%environment%scratch _INDEX_SLICE_PLUS_1_(k)
            _CONCURRENT_LOOP_END_
         end if
      end do

      ! Move to next model
      node => node%next
   end do

   ! Compose total sources-sinks for each state variable, combining model-specific contributions.
   do i=1,size(self%state_variables)
      do j=1,size(self%state_variables(i)%sms_indices)
         k = self%state_variables(i)%sms_indices(j)
         _CONCURRENT_LOOP_BEGIN_EX_(self%environment)
            dy _INDEX_SLICE_PLUS_1_(i) = dy _INDEX_SLICE_PLUS_1_(i) + self%environment%scratch _INDEX_SLICE_PLUS_1_(k)
         _CONCURRENT_LOOP_END_
      end do
   end do

#ifndef _FULL_DOMAIN_IN_SLICE_
   ! Copy newly written diagnostics that need to be saved to global store.
   do i=1,size(self%diagnostic_variables)
      k = self%diagnostic_variables(i)%save_index
      if (self%diagnostic_variables(i)%source==source_do.and.k/=0) then
         j = self%diagnostic_variables(i)%target%write_indices%value
         _CONCURRENT_LOOP_BEGIN_EX_(self%environment)
            self%diag(_PREARG_LOCATION_ k) = self%environment%scratch _INDEX_SLICE_PLUS_1_(j)
         _CONCURRENT_LOOP_END_
      end if
   end do
#endif

   end subroutine fabm_do_rhs
!EOC

subroutine prefetch(environment _ARG_LOCATION_ND_)
   type (type_environment),intent(inout) :: environment
  _DECLARE_LOCATION_ARG_ND_

   integer :: i

#ifdef _FABM_MASK_
   _CONCURRENT_LOOP_BEGIN_EX_NOMASK_(environment)
      environment%prefetch_mask _INDEX_SLICE_ = _FABM_IS_UNMASKED_(environment%mask _INDEX_LOCATION_)
   _CONCURRENT_LOOP_END_
#endif

   do i=1,size(environment%data)
      if (associated(environment%data(i)%p)) then
         _CONCURRENT_LOOP_BEGIN_EX_(environment)
            environment%prefetch _INDEX_SLICE_PLUS_1_(i) = environment%data(i)%p _INDEX_LOCATION_
         _CONCURRENT_LOOP_END_
      end if
   end do
   do i=1,size(environment%data_hz)
      if (associated(environment%data_hz(i)%p)) then
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(environment)
            environment%prefetch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = environment%data_hz(i)%p _INDEX_HORIZONTAL_LOCATION_
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end if
   end do
   do i=1,size(environment%data_scalar)
      if (associated(environment%data_scalar(i)%p)) environment%prefetch_scalar(i) = environment%data_scalar(i)%p
   end do
end subroutine

subroutine prefetch_hz(environment _ARG_LOCATION_VARS_HZ_)
   type (type_environment),intent(inout) :: environment
  _DECLARE_LOCATION_ARG_HZ_

   integer :: i

#ifdef _FABM_MASK_
   _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_NOMASK_(environment)
      environment%prefetch_mask _INDEX_SLICE_ = _FABM_IS_UNMASKED_(environment%mask _INDEX_LOCATION_)
   _CONCURRENT_HORIZONTAL_LOOP_END_
#endif

   do i=1,size(environment%data)
      if (associated(environment%data(i)%p)) then
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(environment)
            environment%prefetch _INDEX_SLICE_PLUS_1_(i) = environment%data(i)%p _INDEX_LOCATION_
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end if
   end do
   do i=1,size(environment%data_hz)
      if (associated(environment%data_hz(i)%p)) then
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(environment)
            environment%prefetch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = environment%data_hz(i)%p _INDEX_HORIZONTAL_LOCATION_
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end if
   end do
   do i=1,size(environment%data_scalar)
      if (associated(environment%data_scalar(i)%p)) environment%prefetch_scalar(i) = environment%data_scalar(i)%p
   end do
end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the biogeochemical
! model tree in the form of production and destruction matrices.
!
! !INTERFACE:
   subroutine fabm_do_ppdd(self _ARG_LOCATION_ND_,pp,dd)
!
! !INPUT PARAMETERS:
   class (type_model),          intent(inout) :: self
  _DECLARE_LOCATION_ARG_ND_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_SLICE_PLUS_2_,intent(inout) :: pp,dd
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
   integer                              :: i,j,k
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&defined(_FABM_USE_1D_LOOP_)
   if (size(pp,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_do_ppdd','Size of first dimension of pp should match loop length.')
   if (size(dd,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_do_ppdd','Size of first dimension of dd should match loop length.')
#endif

   call prefetch(self%environment _ARG_LOCATION_ND_)

   node => self%models%first
   do while (associated(node))
      call node%model%do_ppdd(_ARGUMENTS_ND_IN_,pp,dd)

      ! Copy newly written diagnostics to prefetch
      do i=1,size(node%model%reused_diag)
         if (node%model%reused_diag(i)%source==source_do) then
            j = node%model%reused_diag(i)%read_index
            k = node%model%reused_diag(i)%write_index
            _CONCURRENT_LOOP_BEGIN_EX_(self%environment)
               self%environment%prefetch _INDEX_SLICE_PLUS_1_(j) = self%environment%scratch _INDEX_SLICE_PLUS_1_(k)
            _CONCURRENT_LOOP_END_
         end if
      end do

      node => node%next
   end do

#ifndef _FULL_DOMAIN_IN_SLICE_
   ! Copy newly written diagnostics that need to be saved to global store.
   do i=1,size(self%diagnostic_variables)
      k = self%diagnostic_variables(i)%save_index
      if (self%diagnostic_variables(k)%source==source_do.and.k/=0) then
         j = self%diagnostic_variables(i)%target%write_indices%value
         _CONCURRENT_LOOP_BEGIN_EX_(self%environment)
            self%diag(_PREARG_LOCATION_ k) = self%environment%scratch _INDEX_SLICE_PLUS_1_(j)
         _CONCURRENT_LOOP_END_
      end if
   end do
#endif

   end subroutine fabm_do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks whether the current state is valid, and repairs [clips]
! invalid state variables if requested and possible.
!
! !INTERFACE:
   subroutine fabm_check_state(self _ARG_LOCATION_ND_,repair,valid)
!
! !INPUT PARAMETERS:
   class (type_model),     intent(inout) :: self
   _DECLARE_LOCATION_ARG_ND_
   logical,                intent(in)    :: repair
   logical,                intent(out)   :: valid
!
! !LOCAL PARAMETERS:
   integer                               :: ivar
   type (type_model_list_node), pointer :: node
   real(rk)                              :: value,minimum,maximum
   character(len=256)                    :: err
   type (type_bulk_data_pointer)         :: p
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   valid = .true.

   call prefetch(self%environment _ARG_LOCATION_ND_)

   ! Allow individual models to check their state for their custom constraints, and to perform custom repairs.
   node => self%models%first
   do while (associated(node) .and. valid)
      call node%model%check_state(_ARGUMENTS_ND_IN_,repair,valid)
      if (.not. (valid .or. repair)) return
      node => node%next
   end do

   ! Finally check whether all state variable values lie within their prescribed [constant] bounds.
   ! This is always done, independently of any model-specific checks that may have been called above.

   ! Quick bounds check for the common case where all values are valid.
   do ivar=1,size(self%state_variables)
      minimum = self%state_variables(ivar)%minimum
      maximum = self%state_variables(ivar)%maximum
      _LOOP_BEGIN_EX_(self%environment)
         value = self%environment%prefetch _INDEX_SLICE_PLUS_1_(self%state_variables(ivar)%globalid%read_index)
         if (value<minimum.or.value>maximum) valid = .false.
      _LOOP_END_
   end do
   if (valid) return

   ! Check boundaries for pelagic state variables specified by the models.
   ! If repair is permitted, this clips invalid values to the closest boundary.
   do ivar=1,size(self%state_variables)
      ! Shortcuts to variable information - this demonstrably helps the compiler (ifort).
      p = self%environment%data(self%state_variables(ivar)%globalid%read_index)
      minimum = self%state_variables(ivar)%minimum
      maximum = self%state_variables(ivar)%maximum

      _LOOP_BEGIN_EX_(self%environment)
         value = self%environment%prefetch _INDEX_SLICE_PLUS_1_(self%state_variables(ivar)%globalid%read_index)
         if (value<minimum) then
            ! State variable value lies below prescribed minimum.
            valid = .false.
            if (.not.repair) then
               write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ',trim(self%state_variables(ivar)%name), &
                                                          & ' below minimum value ',minimum
               call log_message(err)
               return
            end if
            p%p _INDEX_LOCATION_ = minimum
         elseif (value>maximum) then
            ! State variable value exceeds prescribed maximum.
            valid = .false.
            if (.not.repair) then
               write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ',trim(self%state_variables(ivar)%name), &
                                                          & ' above maximum value ',maximum
               call log_message(err)
               return
            end if
            p%p _INDEX_LOCATION_ = maximum
         end if
      _LOOP_END_
   end do

   end subroutine fabm_check_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks whether the current bottom state is valid, and repairs [clips]
! invalid state variables if requested and possible.
!
! !INTERFACE:
   subroutine fabm_check_bottom_state(self _ARG_LOCATION_VARS_HZ_,repair,valid)
!
! !INPUT PARAMETERS:
   class (type_model),     intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
   logical,                intent(in)    :: repair
   logical,                intent(out)   :: valid
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   valid = .true.

   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   ! Allow individual models to check their state for their custom constraints, and to perform custom repairs.
   node => self%models%first
   do while (associated(node) .and. valid)
      call node%model%check_bottom_state(_ARGUMENTS_IN_HZ_,repair,valid)
      if (.not. (valid .or. repair)) return
      node => node%next
   end do

   ! Finally check whether all state variable values lie within their prescribed [constant] bounds.
   ! This is always done, independently of any model-specific checks that may have been called above.
   call internal_check_horizontal_state(self _ARG_LOCATION_VARS_HZ_,self%bottom_state_variables,repair,valid)

   end subroutine fabm_check_bottom_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks whether the current bottom state is valid, and repairs [clips]
! invalid state variables if requested and possible.
!
! !INTERFACE:
   subroutine fabm_check_surface_state(self _ARG_LOCATION_VARS_HZ_,repair,valid)
!
! !INPUT PARAMETERS:
   class (type_model),     intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
   logical,                intent(in)    :: repair
   logical,                intent(out)   :: valid
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   valid = .true.

   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   ! Allow individual models to check their state for their custom constraints, and to perform custom repairs.
   node => self%models%first
   do while (associated(node) .and. valid)
      call node%model%check_surface_state(_ARGUMENTS_IN_HZ_,repair,valid)
      if (.not. (valid .or. repair)) return
      node => node%next
   end do

   ! Finally check whether all state variable values lie within their prescribed [constant] bounds.
   ! This is always done, independently of any model-specific checks that may have been called above.
   call internal_check_horizontal_state(self _ARG_LOCATION_VARS_HZ_,self%info%surface_state_variables,repair,valid)

   end subroutine fabm_check_surface_state
!EOC

subroutine internal_check_horizontal_state(self _ARG_LOCATION_VARS_HZ_,state_variables,repair,valid)
   class (type_model),                        intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
   type (type_horizontal_state_variable_info),intent(inout) :: state_variables(:)
   logical,                                   intent(in)    :: repair
   logical,                                   intent(out)   :: valid

   integer                               :: ivar
   real(rk)                              :: value,minimum,maximum
   character(len=256)                    :: err
   type (type_horizontal_data_pointer)   :: p_hz

   ! Check boundaries for horizontal state variables, as prescribed by the owning models.
   ! If repair is permitted, this clips invalid values to the closest boundary.
   do ivar=1,size(state_variables)
      ! Shortcuts to variable information - this demonstrably helps the compiler (ifort).
      p_hz = self%environment%data_hz(state_variables(ivar)%globalid%read_index)
      minimum = state_variables(ivar)%minimum
      maximum = state_variables(ivar)%maximum

      _HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
         value = self%environment%prefetch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(state_variables(ivar)%globalid%read_index)
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
            p_hz%p _INDEX_HORIZONTAL_LOCATION_ = minimum
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
            p_hz%p _INDEX_HORIZONTAL_LOCATION_ = maximum
         end if
      _HORIZONTAL_LOOP_END_
   end do
end subroutine internal_check_horizontal_state

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-water exchange fluxes for all biogeochemical state variables.
! Positive values indicate fluxes into the ocean, negative values indicate fluxes
! out of the ocean. Units are tracer unit * m/s.
!
! !INTERFACE:
   subroutine fabm_do_surface(self _ARG_LOCATION_VARS_HZ_,flux_pel,flux_sf)
!
! !INPUT PARAMETERS:
      class (type_model),                          intent(inout) :: self
      _DECLARE_LOCATION_ARG_HZ_
!
! !INPUT/OUTPUT PARAMETERS:
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(out)          :: flux_pel
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(out),optional :: flux_sf
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
      type (type_model_list_node), pointer :: node
      integer                              :: i,j,k
!-----------------------------------------------------------------------
!BOC
      call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

      ! Initialize scratch variables that will hold surface-attached source-sink terms
      ! or bulk surface fluxes to zero, because they will be incremented.
      do i=1,size(self%zero_surface_indices)
         j = self%zero_surface_indices(i)
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
            self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = 0.0_rk
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end do

#ifndef _FULL_DOMAIN_IN_SLICE_
      ! Copy the value of any horizontal diagnostics kept in global store to scratch space.
      ! This is needed to preserve the value of bottom diagnostics.
      do i=1,size(self%horizontal_diagnostic_variables)
         k = self%horizontal_diagnostic_variables(i)%save_index
         if (k/=0.and.self%horizontal_diagnostic_variables(i)%source==source_unknown) then
            j = self%horizontal_diagnostic_variables(i)%target%write_indices%value
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
               self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = self%diag_hz(_PREARG_LOCATION_HZ_ k)
            _CONCURRENT_HORIZONTAL_LOOP_END_
         end if
      end do
#endif

      node => self%models%first
      do while (associated(node))
         call node%model%do_surface(_ARGUMENTS_IN_HZ_)

         ! Copy newly written diagnostics to prefetch
         do i=1,size(node%model%reused_diag_hz)
            if (node%model%reused_diag_hz(i)%source==source_do_surface.or.node%model%reused_diag_hz(i)%source==source_unknown) then
               j = node%model%reused_diag_hz(i)%read_index
               k = node%model%reused_diag_hz(i)%write_index
               _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
                  self%environment%prefetch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
               _CONCURRENT_HORIZONTAL_LOOP_END_
            end if
         end do

         node => node%next
      end do

      ! Compose surface fluxes for each bulk state variable, combining model-specific contributions.
      flux_pel = 0.0_rk
      do i=1,size(self%state_variables)
         do j=1,size(self%state_variables(i)%surface_flux_indices)
            k = self%state_variables(i)%surface_flux_indices(j)
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
               flux_pel _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = flux_pel _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) + self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
            _CONCURRENT_HORIZONTAL_LOOP_END_
         end do
      end do

      ! Compose total sources-sinks for each surface-bound state variable, combining model-specific contributions.
      if (present(flux_sf)) then
         flux_sf = 0.0_rk
         do i=1,size(self%surface_state_variables)
            do j=1,size(self%surface_state_variables(i)%sms_indices)
               k = self%surface_state_variables(i)%sms_indices(j)
               _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
                  flux_sf _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = flux_sf _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) + self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
               _CONCURRENT_HORIZONTAL_LOOP_END_
            end do
         end do
      end if

#ifndef _FULL_DOMAIN_IN_SLICE_
      ! Copy newly written horizontal diagnostics that need to be saved to global store.
      do i=1,size(self%horizontal_diagnostic_variables)
         k = self%horizontal_diagnostic_variables(i)%save_index
         if (k/=0.and.(self%horizontal_diagnostic_variables(i)%source==source_do_surface.or.self%horizontal_diagnostic_variables(i)%source==source_unknown)) then
            j = self%horizontal_diagnostic_variables(i)%target%write_indices%value
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
               self%diag_hz(_PREARG_LOCATION_HZ_ k) = self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j)
            _CONCURRENT_HORIZONTAL_LOOP_END_
         end if
      end do
#endif

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
   subroutine fabm_do_bottom_rhs(self _ARG_LOCATION_VARS_HZ_,flux_pel,flux_ben)
!
! !INPUT PARAMETERS:
   class (type_model),                          intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(inout) :: flux_pel,flux_ben
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model_list_node), pointer :: node
   integer                              :: i,j,k
!-----------------------------------------------------------------------
!BOC
   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   ! Initialize scratch variables that will hold bottom-attached source-sink terms
   ! or bulk bottom fluxes to zero, because they will be incremented.
   do i=1,size(self%zero_bottom_indices)
      j = self%zero_bottom_indices(i)
      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
         self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = 0.0_rk
      _CONCURRENT_HORIZONTAL_LOOP_END_
   end do

#ifndef _FULL_DOMAIN_IN_SLICE_
   ! Copy the value of any horizontal diagnostics kept in global store to scratch space.
   ! This is needed to preserve the value of surface diagnostics.
   do i=1,size(self%horizontal_diagnostic_variables)
      k = self%horizontal_diagnostic_variables(i)%save_index
      if (k/=0.and.self%horizontal_diagnostic_variables(i)%source==source_unknown) then
         j = self%horizontal_diagnostic_variables(i)%target%write_indices%value
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
            self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = self%diag_hz(_PREARG_LOCATION_HZ_ k)
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end if
   end do
#endif

   node => self%models%first
   do while (associated(node))
      call node%model%do_bottom(_ARGUMENTS_IN_HZ_)

      ! Copy newly written diagnostics to prefetch
      do i=1,size(node%model%reused_diag_hz)
         if (node%model%reused_diag_hz(i)%source==source_do_bottom.or.node%model%reused_diag_hz(i)%source==source_unknown) then
            j = node%model%reused_diag_hz(i)%read_index
            k = node%model%reused_diag_hz(i)%write_index
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
               self%environment%prefetch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
            _CONCURRENT_HORIZONTAL_LOOP_END_
         end if
      end do

      node => node%next
   end do

   ! Compose bottom fluxes for each bulk state variable, combining model-specific contributions.
   do i=1,size(self%state_variables)
      do j=1,size(self%state_variables(i)%bottom_flux_indices)
         k = self%state_variables(i)%bottom_flux_indices(j)
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
            flux_pel _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = flux_pel _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) + self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end do
   end do

   ! Compose total sources-sinks for each bottom-bound state variable, combining model-specific contributions.
   do i=1,size(self%bottom_state_variables)
      do j=1,size(self%bottom_state_variables(i)%sms_indices)
         k = self%bottom_state_variables(i)%sms_indices(j)
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
            flux_ben _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = flux_ben _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) + self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end do
   end do

#ifndef _FULL_DOMAIN_IN_SLICE_
   ! Copy newly written horizontal diagnostics that need to be saved to global store.
   do i=1,size(self%horizontal_diagnostic_variables)
      k = self%horizontal_diagnostic_variables(i)%save_index
      if (k/=0.and.(self%horizontal_diagnostic_variables(i)%source==source_do_bottom.or.self%horizontal_diagnostic_variables(i)%source==source_unknown)) then
         j = self%horizontal_diagnostic_variables(i)%target%write_indices%value
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
            self%diag_hz(_PREARG_LOCATION_HZ_ k) = self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j)
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end if
   end do
#endif

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
   subroutine fabm_do_bottom_ppdd(self _ARG_LOCATION_VARS_HZ_,pp,dd,benthos_offset)
!
! !INPUT PARAMETERS:
   class (type_model),        intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
   integer,                   intent(in)    :: benthos_offset
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_2_,intent(inout) :: pp,dd
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model_list_node), pointer :: node
   integer                              :: i,j,k
!-----------------------------------------------------------------------
!BOC

#ifndef _FULL_DOMAIN_IN_SLICE_
   ! Copy the value of any horizontal diagnostics kept in global store to scratch space.
   ! This is needed to preserve the value of surface diagnostics.
   do i=1,size(self%horizontal_diagnostic_variables)
      k = self%horizontal_diagnostic_variables(i)%save_index
      if (k/=0.and.self%horizontal_diagnostic_variables(i)%source==source_unknown) then
         j = self%horizontal_diagnostic_variables(i)%target%write_indices%value
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
            self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = self%diag_hz(_PREARG_LOCATION_HZ_ k)
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end if
   end do
#endif

   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   node => self%models%first
   do while (associated(node))
      call node%model%do_bottom_ppdd(_ARGUMENTS_IN_HZ_,pp,dd,benthos_offset)

      ! Copy newly written diagnostics to prefetch
      do i=1,size(node%model%reused_diag_hz)
         if (node%model%reused_diag_hz(i)%source==source_do_bottom.or.node%model%reused_diag_hz(i)%source==source_unknown) then
            j = node%model%reused_diag_hz(i)%read_index
            k = node%model%reused_diag_hz(i)%write_index
            _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
               self%environment%prefetch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j) = self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(k)
            _CONCURRENT_HORIZONTAL_LOOP_END_
         end if
      end do

      node => node%next
   end do

#ifndef _FULL_DOMAIN_IN_SLICE_
   ! Copy newly written horizontal diagnostics that need to be saved to global store.
   do i=1,size(self%horizontal_diagnostic_variables)
      k = self%horizontal_diagnostic_variables(i)%save_index
      if (k/=0.and.(self%horizontal_diagnostic_variables(i)%source==source_do_bottom.or.self%horizontal_diagnostic_variables(i)%source==source_unknown)) then
         j = self%horizontal_diagnostic_variables(i)%target%write_indices%value
         _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
            self%diag_hz(_PREARG_LOCATION_HZ_ k) = self%environment%scratch_hz _INDEX_HORIZONTAL_SLICE_PLUS_1_(j)
         _CONCURRENT_HORIZONTAL_LOOP_END_
      end if
   end do
#endif

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
   subroutine fabm_get_vertical_movement(self _ARG_LOCATION_ND_,velocity)
!
! !INPUT PARAMETERS:
   class (type_model),               intent(inout) :: self
   _DECLARE_LOCATION_ARG_ND_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_SLICE_PLUS_1_,intent(out)  :: velocity
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model_list_node), pointer :: node
   integer                              :: i
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&defined(_FABM_USE_1D_LOOP_)
   if (size(velocity,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_get_vertical_movement','Size of first dimension of velocity should match loop length.')
#endif

   call prefetch(self%environment _ARG_LOCATION_ND_)

   ! First set constant sinking rates.
   do i=1,size(self%state_variables)
      ! Use variable-specific constant vertical velocities.
      ! Automatic vectorization by compiler is possible without further modification (ifort 14.0, 2013-12-06)
      _LOOP_BEGIN_EX_(self%environment)
         velocity _INDEX_SLICE_PLUS_1_(i) = self%state_variables(i)%vertical_movement
      _LOOP_END_
   end do

   ! Now allow models to overwrite with spatially-varying sinking rates - if any.
   node => self%models%first
   do while (associated(node))
      call node%model%get_vertical_movement(_ARGUMENTS_ND_IN_,velocity)
      node => node%next
   end do

   end subroutine fabm_get_vertical_movement
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine fabm_get_light_extinction(self _ARG_LOCATION_ND_,extinction)
!
! !INPUT PARAMETERS:
   class (type_model),        intent(inout) :: self
   _DECLARE_LOCATION_ARG_ND_
   real(rk) _DIMENSION_SLICE_,intent(out)   :: extinction
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(NDEBUG)&&defined(_FABM_USE_1D_LOOP_)
   if (size(extinction,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_get_light_extinction','Size of first dimension of extinction should match loop length.')
#endif

   call prefetch(self%environment _ARG_LOCATION_ND_)

   ! Call all models that calculate extinction components to make sure the extinction diagnostic is up to date.
   node => self%extinction_call_list%first
   do while (associated(node))
      call node%model%do(_ARGUMENTS_ND_IN_)
      node => node%next
   end do

   _LOOP_BEGIN_EX_(self%environment)
      extinction _INDEX_SLICE_ = self%environment%data(self%extinction_index)%p _INDEX_LOCATION_
   _LOOP_END_

   end subroutine fabm_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calculate the light field
!
! !INTERFACE:
   subroutine fabm_get_light(self _ARG_LOCATION_VERT_)
!
! !INPUT PARAMETERS:
   class (type_model), intent(inout) :: self
   _DECLARE_LOCATION_ARG_VERT_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
   integer                              :: i,j,k
!EOP
!-----------------------------------------------------------------------
!BOC
   call prefetch(self%environment _ARG_LOCATION_ND_)

   node => self%models%first
   do while (associated(node))
      call node%model%get_light(_ARGUMENTS_VERT_IN_)
      node => node%next
   end do

#ifndef _FULL_DOMAIN_IN_SLICE_
   ! Copy diagnostics that need to be saved to global store.
   do i=1,size(self%diagnostic_variables)
      k = self%diagnostic_variables(i)%save_index
      if (self%diagnostic_variables(i)%source==source_do_column.and.k/=0) then
         j = self%diagnostic_variables(i)%target%write_indices%value
         _CONCURRENT_VERTICAL_LOOP_BEGIN_EX_(self%environment)
            self%diag(_PREARG_LOCATION_ k) = self%environment%scratch _INDEX_SLICE_PLUS_1_(j)
         _CONCURRENT_VERTICAL_LOOP_END_
      end if
   end do
#endif

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
   subroutine fabm_get_drag(self _ARG_LOCATION_VARS_HZ_,drag)
!
! !INPUT PARAMETERS:
   class (type_model),                   intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
   real(rk) _DIMENSION_HORIZONTAL_SLICE_,intent(out)   :: drag
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
!EOP
!-----------------------------------------------------------------------
!BOC
   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   drag = 1.0_rk
   node => self%models%first
   do while (associated(node))
      call node%model%get_drag(_ARGUMENTS_IN_HZ_,drag)
      node => node%next
   end do

   end subroutine fabm_get_drag
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get biogeochemistry-induced change in surface albedo
! This feedback is represented by an albedo value (0-1) relating to biogeochemistry.
!
! !INTERFACE:
   subroutine fabm_get_albedo(self _ARG_LOCATION_VARS_HZ_,albedo)
!
! !INPUT PARAMETERS:
   class (type_model),                   intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
   real(rk) _DIMENSION_HORIZONTAL_SLICE_,intent(out)   :: albedo
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
!EOP
!-----------------------------------------------------------------------
!BOC
   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   albedo = 0.0_rk
   node => self%models%first
   do while (associated(node))
      call node%model%get_albedo(_ARGUMENTS_IN_HZ_,albedo)
      node => node%next
   end do

   end subroutine fabm_get_albedo
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities in pelagic domain
!
! !INTERFACE:
   subroutine fabm_get_conserved_quantities(self _ARG_LOCATION_ND_,sums)
!
! !INPUT PARAMETERS:
   class (type_model),               intent(inout) :: self
   _DECLARE_LOCATION_ARG_ND_
   real(rk) _DIMENSION_SLICE_PLUS_1_,intent(out)   :: sums
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model_list_node), pointer :: node
   integer :: i
!-----------------------------------------------------------------------
!BOC
   call prefetch(self%environment _ARG_LOCATION_ND_)

   node => self%conserved_quantity_call_list%first
   do while (associated(node))
      call node%model%do(_ARGUMENTS_ND_IN_)
      node => node%next
   end do

   do i=1,size(self%conserved_quantities)
      _CONCURRENT_LOOP_BEGIN_EX_(self%environment)
         sums _INDEX_SLICE_PLUS_1_(i) = self%environment%data(self%conserved_quantities(i)%index)%p _INDEX_LOCATION_
      _CONCURRENT_LOOP_END_
   end do

   end subroutine fabm_get_conserved_quantities
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities at top and bottom boundaries
!
! !INTERFACE:
   subroutine fabm_get_horizontal_conserved_quantities(self _ARG_LOCATION_VARS_HZ_,sums)
!
! !INPUT PARAMETERS:
   class (type_model),                          intent(inout) :: self
   _DECLARE_LOCATION_ARG_HZ_
   real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(out)   :: sums
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   integer :: i

!-----------------------------------------------------------------------
!BOC
   call prefetch_hz(self%environment _ARG_LOCATION_VARS_HZ_)

   do i=1,size(self%conserved_quantities)
      if (associated(self%conserved_quantities(i)%horizontal_sum)) &
         call self%conserved_quantities(i)%horizontal_sum%evaluate_horizontal(_ARGUMENTS_IN_HZ_)
      _CONCURRENT_HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
         sums _INDEX_HORIZONTAL_SLICE_PLUS_1_(i) = self%environment%data_hz(self%conserved_quantities(i)%horizontal_index)%p _INDEX_HORIZONTAL_LOCATION_
      _CONCURRENT_HORIZONTAL_LOOP_END_
   end do

   end subroutine fabm_get_horizontal_conserved_quantities
!EOC

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
         ! Start of simulation; set entire history equal to current value.
         do i=1,expression%n+3
            expression%history(_PREARG_LOCATION_DIMENSIONS_ i) = self%environment%data(expression%in)%p
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
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) &
            - expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest)/expression%n
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest) = (1.0_rk-weight_right)*expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+1) &
            + weight_right*self%environment%data(expression%in)%p
         expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) &
            + expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest)/expression%n

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
      expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+1) = self%environment%data(expression%in)%p
      expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+3) = expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+2) &
         + frac_outside*(-expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%ioldest) + expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+1))

      expression%last_time = t
   end subroutine

   subroutine update_horizontal_temporal_mean(expression)
      class (type_horizontal_temporal_mean), intent(inout) :: expression
      integer  :: i
      real(rk) :: weight_right,frac_outside

      if (expression%ioldest==-1) then
         ! Start of simulation; set entire history equal to current value.
         do i=1,expression%n+3
            expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ i) = self%environment%data_hz(expression%in)%p
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
         expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+2) = expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+2) &
            - expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%ioldest)/expression%n
         expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%ioldest) = (1.0_rk-weight_right)*expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+1) &
            + weight_right*self%environment%data_hz(expression%in)%p
         expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+2) = expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+2) &
            + expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%ioldest)/expression%n

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
      expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+1) = self%environment%data_hz(expression%in)%p
      expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+3) = expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+2) &
         + frac_outside*(-expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%ioldest) + expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+1))

      expression%last_time = t
   end subroutine

end subroutine

subroutine copy_write_indices(source,target)
   type (type_link_list),intent(in) :: source
   integer,allocatable              :: target(:)

   integer                  :: n,nold
   type (type_link),pointer :: link
   integer,allocatable      :: old(:)

   n = 0
   link => source%first
   do while (associated(link))
      if (.not.link%target%write_indices%is_empty()) n = n + 1
      link => link%next
   end do

   if (allocated(target)) then
      nold = size(target)
      allocate(old(nold))
      old(:) = target
      deallocate(target)
   else
      nold = 0
   end if

   allocate(target(nold+n))
   if (allocated(old)) target(1:nold) = old

   n = nold
   link => source%first
   do while (associated(link))
      if (.not.link%target%write_indices%is_empty()) then
         n = n + 1
         call link%target%write_indices%append(target(n))
      end if
      link => link%next
   end do
end subroutine copy_write_indices

function create_external_bulk_id(variable) result(id)
   type (type_internal_variable),intent(inout),target :: variable
   type (type_bulk_variable_id) :: id

   if (variable%domain/=domain_bulk) call driver%fatal_error('create_external_bulk_id','BUG: called on non-bulk variable.')
   id%variable => variable
   if (.not.variable%read_indices%is_empty())  id%read_index = variable%read_indices%value
   if (.not.variable%state_indices%is_empty()) id%state_index = variable%state_indices%value
end function create_external_bulk_id

function create_external_horizontal_id(variable) result(id)
   type (type_internal_variable),intent(inout),target :: variable
   type (type_horizontal_variable_id) :: id

   if (variable%domain/=domain_horizontal.and.variable%domain/=domain_surface.and.variable%domain/=domain_bottom) &
      call driver%fatal_error('create_external_horizontal_id','BUG: called on non-horizontal variable.')
   id%variable => variable
   if (.not.variable%read_indices%is_empty())  id%read_index = variable%read_indices%value
   if (.not.variable%state_indices%is_empty()) id%state_index = variable%state_indices%value
end function create_external_horizontal_id

function create_external_scalar_id(variable) result(id)
   type (type_internal_variable),intent(inout),target :: variable
   type (type_scalar_variable_id) :: id
   if (variable%domain/=domain_scalar) call driver%fatal_error('create_external_scalar_id','BUG: called on non-scalar variable.')
   id%variable => variable
   if (.not.variable%read_indices%is_empty()) id%read_index = variable%read_indices%value
end function create_external_scalar_id

function create_external_bulk_id_for_standard_name(self,standard_variable) result(id)
   class (type_model),                intent(in) :: self
   type (type_bulk_standard_variable),intent(in) :: standard_variable
   type (type_bulk_variable_id)                  :: id

   type (type_link), pointer :: link

   link => self%root%links%first
   do while (associated(link))
      if (associated(link%target%standard_variable)) then
         select type (variable=>link%target%standard_variable)
            class is (type_bulk_standard_variable)
               if (variable%compare(standard_variable)) then
                  id = create_external_bulk_id(link%target)
                  return
               end if
         end select
      end if
      link => link%next
   end do
end function create_external_bulk_id_for_standard_name

function create_external_horizontal_id_for_standard_name(self,standard_variable) result(id)
   class (type_model),                      intent(in) :: self
   type (type_horizontal_standard_variable),intent(in) :: standard_variable
   type (type_horizontal_variable_id)                  :: id

   type (type_link), pointer :: link

   link => self%root%links%first
   do while (associated(link))
      if (associated(link%target%standard_variable)) then
         select type (variable=>link%target%standard_variable)
            class is (type_horizontal_standard_variable)
               if (variable%compare(standard_variable)) then
                  id = create_external_horizontal_id(link%target)
                  return
               end if
         end select
      end if
      link => link%next
   end do
end function create_external_horizontal_id_for_standard_name

function create_external_scalar_id_for_standard_name(self,standard_variable) result(id)
   class (type_model),                  intent(in) :: self
   type (type_global_standard_variable),intent(in) :: standard_variable
   type (type_scalar_variable_id)                  :: id

   type (type_link), pointer :: link

   link => self%root%links%first
   do while (associated(link))
      if (associated(link%target%standard_variable)) then
         select type (variable=>link%target%standard_variable)
            class is (type_global_standard_variable)
               if (variable%compare(standard_variable)) then
                  id = create_external_scalar_id(link%target)
                  return
               end if
         end select
      end if
      link => link%next
   end do
end function create_external_scalar_id_for_standard_name

recursive subroutine set_diagnostic_indices(self)
   class (type_base_model),intent(inout) :: self

   type (type_link), pointer :: link
   integer                   :: n,n_hz
   type (type_model_list_node),pointer :: node

   n = 0
   n_hz = 0
   link => self%links%first
   do while (associated(link))
      if (index(link%name,'/')==0.and.associated(link%original%write_index).and..not.link%target%read_indices%is_empty()) then
         ! Variable is a diagnostic written to by current model, and read by at least one model.
         select case (link%target%domain)
            case (domain_bulk)
               n = n + 1
            case (domain_horizontal,domain_surface,domain_bottom)
               n_hz = n_hz + 1
         end select
      end if
      link => link%next
   end do

   allocate(self%reused_diag(n))
   allocate(self%reused_diag_hz(n_hz))

   n = 0
   n_hz = 0
   link => self%links%first
   do while (associated(link))
      if (index(link%name,'/')==0.and.associated(link%original%write_index).and..not.link%target%read_indices%is_empty()) then
         ! Variable is written to by current model, and read by someone.
         select case (link%target%domain)
            case (domain_bulk)
               n = n + 1
               self%reused_diag(n)%source = link%target%source
               call link%target%write_indices%append(self%reused_diag(n)%write_index)
               call link%target%read_indices%append (self%reused_diag(n)%read_index)
            case (domain_horizontal,domain_surface,domain_bottom)
               n_hz = n_hz + 1
               self%reused_diag_hz(n_hz)%source = link%target%source
               call link%target%write_indices%append(self%reused_diag_hz(n_hz)%write_index)
               call link%target%read_indices%append (self%reused_diag_hz(n_hz)%read_index)
         end select
      end if
      link => link%next
   end do

   node => self%children%first
   do while (associated(node))
      call set_diagnostic_indices(node%model)
      node => node%next
   end do
end subroutine set_diagnostic_indices

subroutine create_readable_variable_registry(self)
   class (type_model),intent(inout),target :: self

   integer :: nread,nread_hz,nread_scalar
   type (type_link), pointer :: link

   nread = 0
   nread_hz = 0
   nread_scalar = 0
   link => self%links_postcoupling%first
   do while (associated(link))
      if (.not.link%target%read_indices%is_empty()) then
         select case (link%target%domain)
            case (domain_bulk)
               ! Bulk variable read by one or more models
               nread = nread+1
               call link%target%read_indices%set_value(nread)
            case (domain_horizontal,domain_surface,domain_bottom)
               ! Horizontal variable read by one or more models
               nread_hz = nread_hz+1
               call link%target%read_indices%set_value(nread_hz)
            case (domain_scalar)
               ! Scalar variable read by one or more models
               nread_scalar = nread_scalar+1
               call link%target%read_indices%set_value(nread_scalar)
         end select
      end if   
      link => link%next
   end do

   allocate(self%environment%data       (nread))
   allocate(self%environment%data_hz    (nread_hz))
   allocate(self%environment%data_scalar(nread_scalar))
end subroutine create_readable_variable_registry

function variable_from_data_index(self,index,domain) result(variable)
   class (type_model),intent(inout),target :: self
   integer,           intent(in)           :: index,domain
   type (type_internal_variable),pointer   :: variable

   type (type_link), pointer :: link

   nullify(variable)
   link => self%links_postcoupling%first
   do while (associated(link))
      if (link%target%read_indices%value==index .and. iand(link%target%domain,domain)/=0) then
         variable => link%target
         return
      end if   
      link => link%next
   end do
end function variable_from_data_index

subroutine filter_readable_variable_registry(self)
   class (type_model),intent(inout),target :: self

   integer :: nread,nread_hz,nread_scalar
   type (type_link), pointer :: link

   ! For all readable variables that are designated as "optional" by the biogeochemical model, and have not been
   ! assigned a value, set their index in the biogeochemical model to -1 so the model is aware that it is not used/unavailable.
   nread = 0
   nread_hz = 0
   nread_scalar = 0
   link => self%links_postcoupling%first
   do while (associated(link))
      if (.not.link%target%read_indices%is_empty()) then
         select case (link%target%domain)
            case (domain_bulk)
               nread = nread+1
               if (link%target%presence==presence_external_optional &
                   .and..not.associated(self%environment%data(nread)%p)) &
                  call link%target%read_indices%set_value(-1)
            case (domain_horizontal,domain_surface,domain_bottom)
               nread_hz = nread_hz+1
               if (link%target%presence==presence_external_optional &
                   .and..not.associated(self%environment%data_hz(nread_hz)%p)) &
                  call link%target%read_indices%set_value(-1)
            case (domain_scalar)
               nread_scalar = nread_scalar+1
               if (link%target%presence==presence_external_optional &
                   .and..not.associated(self%environment%data_scalar(nread_scalar)%p)) &
                  call link%target%read_indices%set_value(-1)
         end select
      end if   
      link => link%next
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
   class (type_base_model),                        pointer :: model

   type (type_aggregate_variable),    pointer :: aggregate_variable
   type (type_set) :: dependencies,dependencies_hz,dependencies_scalar
   type (type_model_list_node),   pointer :: model_node

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

   call set_diagnostic_indices(self%root)

   ! Count number of conserved quantities and allocate associated array.
   ncons = 0
   aggregate_variable => self%root%first_aggregate_variable
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
   aggregate_variable => self%root%first_aggregate_variable
   do while (associated(aggregate_variable))
      if (aggregate_variable%standard_variable%conserved) then
         ncons = ncons + 1
         consvar => self%conserved_quantities(ncons)
         consvar%standard_variable = aggregate_variable%standard_variable
         consvar%name = trim(consvar%standard_variable%name)
         consvar%units = trim(consvar%standard_variable%units)
         consvar%long_name = trim(consvar%standard_variable%name)
         consvar%path = trim(consvar%standard_variable%name)
         consvar%aggregate_variable => aggregate_variable

         ! Get read index of the variable that will keep the total of conserved quantity in bulk domain.
         ! If this variable is a diagnostic, add the source model and any dependencies to the call list
         ! for conserved quantity calculations.
         object => self%root%find_object(trim(aggregate_variable%standard_variable%name))
         if (.not.associated(object)) call driver%fatal_error('classify_variables', &
            'BUG: conserved quantity '//trim(aggregate_variable%standard_variable%name)//' was not created')
         call object%read_indices%append(consvar%index)
         if (.not.object%write_indices%is_empty()) call find_dependencies(object%owner,self%conserved_quantity_call_list)

         ! Store pointer to total of conserved quantity at surface + bottom.
         object => self%root%find_object(trim(aggregate_variable%standard_variable%name)//'_at_interfaces')
         if (.not.associated(object)) call driver%fatal_error('classify_variables', &
            'BUG: conserved quantity '//trim(aggregate_variable%standard_variable%name)//'_at_interfaces was not created')
         call object%read_indices%append(consvar%horizontal_index)

         ! Store pointer to model that computes total of conserved quantity at surface + bottom, so we can force recomputation.
         model => self%root%find_model(trim(aggregate_variable%standard_variable%name)//'_at_interfaces_calculator')
         if (associated(model)) then
            select type (model)
               class is (type_horizontal_weighted_sum)
                  consvar%horizontal_sum => model
            end select
         end if
      end if
      aggregate_variable => aggregate_variable%next
   end do

   ! Get link to extinction variable.
   object => self%root%find_object(trim(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux%name))
   if (.not.associated(object)) call driver%fatal_error('classify_variables', &
      'BUG: variable attenuation_coefficient_of_photosynthetic_radiative_flux was not created')
   call object%read_indices%append(self%extinction_index)

   ! From this point on, variables will stay as they are.
   ! Coupling is done, and the framework will not add further read indices.

   call create_readable_variable_registry(self)

   ! Count number of bulk variables in various categories.
   nstate = 0
   ndiag  = 0
   nstate_bot  = 0
   nstate_surf = 0
   ndiag_hz    = 0
   link => self%links_postcoupling%first
   do while (associated(link))
      object => link%target
      select case (object%domain)
         case (domain_bulk)
            if (.not.object%write_indices%is_empty()) then
               ! Bulk diagnostic variable.
               ndiag = ndiag+1
            elseif (.not.object%state_indices%is_empty()) then
               ! Bulk state variable.
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
            elseif (.not.object%state_indices%is_empty()) then
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
         case (domain_bulk)
            if (.not.object%write_indices%is_empty()) then
               ! Bulk diagnostic variable
               ndiag = ndiag + 1
               diagvar => self%diagnostic_variables(ndiag)
               call copy_variable_metadata(object,diagvar)
               if (associated(object%standard_variable)) then
                  select type (standard_variable=>object%standard_variable)
                     type is (type_bulk_standard_variable)
                        diagvar%standard_variable = standard_variable
                  end select
               end if
               diagvar%time_treatment    = output2time_treatment(diagvar%output)
               diagvar%save = diagvar%output/=output_none
               diagvar%source = object%source
            elseif (object%presence==presence_internal.and..not.object%state_indices%is_empty()) then
               ! Bulk state variable
               nstate = nstate + 1
               statevar => self%state_variables(nstate)
               call copy_variable_metadata(object,statevar)
               statevar%globalid                  = create_external_bulk_id(object)
               if (associated(object%standard_variable)) then
                  select type (standard_variable=>object%standard_variable)
                     type is (type_bulk_standard_variable)
                        statevar%standard_variable = standard_variable
                  end select
               end if
               statevar%initial_value             = object%initial_value
               statevar%vertical_movement         = object%vertical_movement
               statevar%no_precipitation_dilution = object%no_precipitation_dilution
               statevar%no_river_dilution         = object%no_river_dilution
            end if
         case (domain_horizontal,domain_surface,domain_bottom)
            if (.not.object%write_indices%is_empty()) then
               ! Horizontal diagnostic variable
               ndiag_hz = ndiag_hz + 1
               hz_diagvar => self%horizontal_diagnostic_variables(ndiag_hz)
               call copy_variable_metadata(object,hz_diagvar)
               if (associated(object%standard_variable)) then
                  select type (standard_variable=>object%standard_variable)
                     type is (type_horizontal_standard_variable)
                        hz_diagvar%standard_variable = standard_variable
                  end select
               end if
               hz_diagvar%time_treatment    = output2time_treatment(hz_diagvar%output)
               hz_diagvar%save = hz_diagvar%output/=output_none
               hz_diagvar%source = object%source
            elseif (object%presence==presence_internal.and..not.object%state_indices%is_empty()) then
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
               hz_statevar%globalid          = create_external_horizontal_id(object)
               if (associated(object%standard_variable)) then
                  select type (standard_variable=>object%standard_variable)
                     type is (type_horizontal_standard_variable)
                        hz_statevar%standard_variable = standard_variable
                  end select
               end if
               hz_statevar%initial_value     = object%initial_value
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
            case (domain_bulk);                                    call dependencies%add(link%name)
            case (domain_horizontal,domain_surface,domain_bottom); call dependencies_hz%add(link%name)
            case (domain_scalar);                                  call dependencies_scalar%add(link%name)
         end select
         if (associated(object%standard_variable)) then
            if (object%standard_variable%name/='') then
               select case (object%domain)
                  case (domain_bulk);                                    call dependencies%add(object%standard_variable%name)
                  case (domain_horizontal,domain_surface,domain_bottom); call dependencies_hz%add(object%standard_variable%name)
                  case (domain_scalar);                                  call dependencies_scalar%add(object%standard_variable%name)
               end select
            end if
         end if
      end if
      link => link%next
   end do
   call dependencies%to_array(self%dependencies)
   call dependencies_hz%to_array(self%dependencies_hz)
   call dependencies_scalar%to_array(self%dependencies_scalar)
end subroutine classify_variables

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

   subroutine custom_extinction_calculator_initialize(self,configunit)
      class (type_custom_extinction_calculator), intent(inout), target :: self
      integer,                                   intent(in)            :: configunit

      call self%register_diagnostic_variable(self%id_output,'total','1/m','total custom light extinction',output=output_none)
      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
         self%id_output)
   end subroutine

   subroutine custom_extinction_calculator_do(self,_ARGUMENTS_DO_)
      class (type_custom_extinction_calculator),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      type (type_model_list_node), pointer :: node
      real(rk) _DIMENSION_SLICE_AUTOMATIC_ :: extinction

      ! Ask all active biogeochemical models for custom light extinction values.
      ! (custom as in: not computed by FABM from concentration+specific_light_extincton,
      ! but by the model's own implementation of the get_light_extinction routine)
      extinction = 0.0_rk
      node => self%models%first
      do while (associated(node))
         call node%model%get_light_extinction(_ARGUMENTS_ND_,extinction)
         node => node%next
      end do

      ! Transfer computed light extinction values to the output diagnostic.
      _CONCURRENT_LOOP_BEGIN_
         _SET_DIAGNOSTIC_(self%id_output,extinction _INDEX_SLICE_)
      _CONCURRENT_LOOP_END_
   end subroutine

   subroutine add_to_model_call_list(self,variable_name,list)
      class (type_model),     intent(inout),target :: self
      character(len=*),       intent(in)           :: variable_name
      type (type_model_list), intent(inout)        :: list

      type (type_internal_variable),pointer :: object

      object => self%root%find_object(variable_name)
      if (.not.object%write_indices%is_empty()) call find_dependencies(object%owner,list)
   end subroutine add_to_model_call_list

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

end module fabm

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

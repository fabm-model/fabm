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
! of this module to access biogeochemistry. Specific biogeochemical models must
! be referenced in this module to be available in FABM. Locations where
! specific biogeochemical models may be referenced are indicated by
! ADD\_NEW\_MODEL\_HERE strings in the code comments.
!
! For more information, see the documentation at /doc/documentation.pdf.
!
! !USES:
   use fabm_standard_variables,only: type_bulk_standard_variable, type_horizontal_standard_variable, &
                                     type_global_standard_variable, initialize_standard_variables
   use fabm_types
   use fabm_library
   use fabm_expressions
   use fabm_driver
   use fabm_properties
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
   public fabm_get_conserved_quantities, fabm_get_horizontal_conserved_quantities, fabm_state_to_conserved_quantities

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
   public create_external_variable_id
!
! !PUBLIC TYPES:
!
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
   ! Derived types for variable metadata used by host models.
   ! ====================================================================================================

   type,abstract :: type_external_variable
      character(len=attribute_length) :: name          = ''
      character(len=attribute_length) :: long_name     = ''
      character(len=attribute_length) :: units         = ''
      character(len=attribute_length) :: path          = ''
      real(rk)                        :: minimum       = -1.e20_rk
      real(rk)                        :: maximum       =  1.e20_rk
      real(rk)                        :: missing_value = -2.e20_rk
      integer                         :: output        = output_instantaneous ! See output_* parameters above
      type (type_property_dictionary) :: properties
      integer                         :: externalid    = 0                    ! Identifier to be used freely by host
   end type

!  Derived type describing a state variable
   type,extends(type_external_variable) :: type_state_variable_info
      type (type_bulk_standard_variable) :: standard_variable
      real(rk)                           :: initial_value             = 0.0_rk
      real(rk)                           :: vertical_movement         = 0.0_rk  ! Vertical movement (m/s, <0: sinking, >0: floating)
      logical                            :: no_precipitation_dilution = .false.
      logical                            :: no_river_dilution         = .false.
      type (type_bulk_variable_id)       :: globalid
   end type type_state_variable_info

   type,extends(type_external_variable) :: type_horizontal_state_variable_info
      type (type_horizontal_standard_variable) :: standard_variable
      real(rk)                                 :: initial_value = 0.0_rk
      type (type_horizontal_variable_id)       :: globalid
   end type type_horizontal_state_variable_info

!  Derived type describing a diagnostic variable
   type,extends(type_external_variable) :: type_diagnostic_variable_info
      type (type_bulk_standard_variable) :: standard_variable
      integer                            :: time_treatment = time_treatment_last ! See time_treatment_* parameters above
      type (type_bulk_variable_id)       :: globalid
   end type type_diagnostic_variable_info

   type,extends(type_external_variable) :: type_horizontal_diagnostic_variable_info
      type (type_horizontal_standard_variable) :: standard_variable
      integer                                  :: time_treatment = time_treatment_last ! See time_treatment_* parameters above
      type (type_horizontal_variable_id)       :: globalid
   end type type_horizontal_diagnostic_variable_info

!  Derived type describing a conserved quantity
   type,extends(type_external_variable) :: type_conserved_quantity_info
      type (type_bulk_standard_variable)            :: standard_variable
      type (type_aggregate_variable),pointer        :: aggregate_variable
      type (type_bulk_data_pointer)                 :: data
      type (type_horizontal_data_pointer)           :: horizontal_data
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

      type (type_bulk_data_pointer) :: extinction_data

      ! Derived types that belong to specific biogeochemical models.
      !type (type_metu_mnemiopsis)           :: metu_mnemiopsis
      ! ADD_NEW_MODEL_HERE - required if the model groups its data in a custom derived type
   contains
      procedure :: link_bulk_data_by_id   => fabm_link_bulk_data_by_id
      procedure :: link_bulk_data_by_sn   => fabm_link_bulk_data_by_sn
      procedure :: link_bulk_data_by_name => fabm_link_bulk_data_by_name
      generic :: link_bulk_data => link_bulk_data_by_id,link_bulk_data_by_sn,link_bulk_data_by_name

      procedure :: link_horizontal_data_by_id   => fabm_link_horizontal_data_by_id
      procedure :: link_horizontal_data_by_sn   => fabm_link_horizontal_data_by_sn
      procedure :: link_horizontal_data_by_name => fabm_link_horizontal_data_by_name
      generic :: link_horizontal_data => link_horizontal_data_by_id,link_horizontal_data_by_sn,link_horizontal_data_by_name

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
      module procedure fabm_link_bulk_data_by_id
      module procedure fabm_link_bulk_data_by_sn
      module procedure fabm_link_bulk_data_by_name
   end interface

   ! Subroutine for providing FABM with variable data on horizontal slices of the domain.
   interface fabm_link_horizontal_data
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

   interface fabm_do_surface
      module procedure fabm_do_surface_nosfstate
      module procedure fabm_do_surface_sfstate
   end interface

   ! For backward compatibility only:
   interface fabm_do_benthos
      module procedure fabm_do_bottom_rhs
      module procedure fabm_do_bottom_ppdd
   end interface
   interface fabm_get_surface_exchange
      module procedure fabm_do_surface_nosfstate
      module procedure fabm_do_surface_sfstate
   end interface

   interface create_external_variable_id
      module procedure create_external_bulk_id
      module procedure create_external_horizontal_id
      module procedure create_external_scalar_id
      module procedure create_external_bulk_id_for_standard_name
      module procedure create_external_horizontal_id_for_standard_name
      module procedure create_external_scalar_id_for_standard_name
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
      type (type_model_list_node),   pointer :: node
      type (type_aggregate_variable),pointer :: aggregate_variable
!EOP
!-----------------------------------------------------------------------
!BOC
      self%info => self ! For backward compatibility (pre 11/2013 hosts only)

      ! Make sure a variable for light extinction is created at the root level when calling freeze_model_info.
      ! This variable is used from fabm_get_light_extinction.
      aggregate_variable => get_aggregate_variable(self%root,standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      aggregate_variable%bulk_required = .true.

      ! Filter out expressions that FABM can handle itself.
      ! The remainder, if any, must be handled by the host model.
      call filter_expressions(self)

      ! This will resolve all FABM dependencies and generate final authorative lists of variables of different types.
      call freeze_model_info(self%root)

      ! Build final authorative arrays with variable metadata .
      call classify_variables(self)

      ! Determine the order in which individual biogeochemical models should be called.
      node => self%root%children%first
      do while (associated(node))
         call find_dependencies(node%model,self%models)
         node => node%next
      end do

      ! Consistency check (of find_dependencies): does every model (except the root) appear exactly once in the call list?
      call test_presence(self,self%root)

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
         select type (current)
            class is (type_vertical_integral)
#ifndef _FABM_DEPTH_DIMENSION_INDEX_
               ! For models without depth dimension, FABM can calculate depth averages and integrals itself.
               allocate(integral)
               integral%minimum_depth = current%minimum_depth
               integral%maximum_depth = current%maximum_depth
               integral%average       = current%average
               call self%root%add_child(integral,current%output_name,configunit=-1)
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
!EOP
!-----------------------------------------------------------------------
!BOC
      child => model%children%first
      do while (associated(child))
         if (self%models%count(child%model)/=1) call fatal_error('fabm_initialize::test_presence', &
            'BUG: Model "'//trim(model%name)//'" is not called exactly one time.')
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
  integer                              :: ivar
  class (type_expression),pointer      :: expression
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Forward domain to individual biogeochemical models.
   ! Also determine whether one of the models needs conservation checks for its sink and source terms.
   node => self%models%first
   do while (associated(node))
      call node%model%set_domain(_LOCATION_)
      node => node%next
   end do

   ! Allocate arrays for diagnostic variables defined on the full domain and on horizontonal slices.
   allocate(self%environment%diag   (_PREARG_LOCATION_ size(self%diagnostic_variables)))
   allocate(self%environment%diag_hz(_PREARG_LOCATION_HZ_ size(self%horizontal_diagnostic_variables)))
   self%environment%nstate = size(self%state_variables)

   if (_VARIABLE_REGISTERED_(self%root%id_zero)) then
      allocate(self%environment%zero _INDEX_LOCATION_)
      self%environment%zero = 0.0_rk
      call fabm_link_bulk_data(self,self%root%id_zero%link%name,self%environment%zero)
   end if
   if (_VARIABLE_REGISTERED_(self%root%id_zero_hz)) then
      allocate(self%environment%zero_hz _INDEX_HORIZONTAL_LOCATION_)
      self%environment%zero_hz = 0.0_rk
      call fabm_link_horizontal_data(self,self%root%id_zero_hz%link%name,self%environment%zero_hz)
   end if

   ! Initialize diagnostic variables to missing value.
   do ivar=1,size(self%diagnostic_variables)
      self%environment%diag(_PREARG_LOCATION_DIMENSIONS_ ivar) = self%diagnostic_variables(ivar)%missing_value
   end do
   do ivar=1,size(self%horizontal_diagnostic_variables)
      self%environment%diag_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ ivar) = self%horizontal_diagnostic_variables(ivar)%missing_value
   end do

   ! If diagnostic variables also appear as dependency, send the corresponding array slice for generic read-only access.
   do ivar=1,size(self%diagnostic_variables)
      call fabm_link_bulk_data(self,self%diagnostic_variables(ivar)%globalid, &
                               self%environment%diag(_PREARG_LOCATION_DIMENSIONS_ ivar))
   end do
   do ivar=1,size(self%horizontal_diagnostic_variables)
      call fabm_link_horizontal_data(self,self%horizontal_diagnostic_variables(ivar)%globalid, &
                                     self%environment%diag_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ ivar))
   end do

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
   class (type_model),target,intent(inout)                            :: self
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
   class (type_model),intent(in),target :: self
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  type (type_link),        pointer :: link
  logical                          :: ready
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ready = .true.

#ifdef _FABM_MASK_
   if (.not.associated(self%environment%mask)) call log_message('spatial mask has not been set.')
#endif

   link => self%root%first_link
   do while (associated(link))
      if (link%owner) then
         select type (object=>link%target)
            class is (type_bulk_variable)
               if (allocated(object%alldata).and..not.object%presence==presence_external_optional) then
                  if (.not.associated(object%alldata(1)%p%p)) then
                     call log_message('data for dependency "'//trim(link%name)// &
                        & '", defined on the full model domain, have not been provided.')
                     ready = .false.
                  end if
               end if
            class is (type_horizontal_variable)
               if (allocated(object%alldata).and..not.object%presence==presence_external_optional) then
                  if (.not.associated(object%alldata(1)%p%p)) then
                     call log_message('data for dependency "'//trim(link%name)// &
                        &  '", defined on a horizontal slice of the model domain, have not been provided.')
                     ready = .false.
                  end if
               end if
            class is (type_scalar_variable)
               if (allocated(object%alldata).and..not.object%presence==presence_external_optional) then
                  if (.not.associated(object%alldata(1)%p%p)) then
                     call log_message('data for dependency "'//trim(link%name)// &
                        &  '", defined as global scalar quantity, have not been provided.')
                     ready = .false.
                  end if
               end if
         end select
      end if
      link => link%next
   end do

   if (.not.ready) call fatal_error('fabm_check_ready','FABM is lacking required data.')

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
   do ivar=1,size(self%state_variables)
      ! Shortcuts to variable information - this demonstrably helps the compiler with vectorization (ifort).
      p = self%state_variables(ivar)%globalid%p
      initial_value = self%state_variables(ivar)%initial_value
      _LOOP_BEGIN_EX_(self%environment)
         _SET_EX_(p,initial_value)
      _LOOP_END_
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
   ! Initialize bottom variables
   do ivar=1,size(self%bottom_state_variables)
      p = self%bottom_state_variables(ivar)%globalid%p
      initial_value = self%bottom_state_variables(ivar)%initial_value
      _HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
         _SET_HORIZONTAL_EX_(p,initial_value)
      _HORIZONTAL_LOOP_END_
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
   ! Initialize surface variables
   do ivar=1,size(self%surface_state_variables)
      p = self%surface_state_variables(ivar)%globalid%p
      initial_value = self%surface_state_variables(ivar)%initial_value
      _HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
         _SET_HORIZONTAL_EX_(p,initial_value)
      _HORIZONTAL_LOOP_END_
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
   link => self%root%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_bulk_variable)
            if (object%name==name.or.get_safe_name(object%name)==name) then
               id = create_external_variable_id(object)
               return
            end if
      end select
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_bulk_variable)
            if (object%standard_variable%name==name) then
               id = create_external_variable_id(self%root,object%standard_variable)
               return
            end if
      end select
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
   id = create_external_variable_id(self%root,standard_variable)

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
   link => self%root%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_horizontal_variable)
            if (object%name==name.or.get_safe_name(object%name)==name) then
               id = create_external_variable_id(object)
               return
            end if
      end select
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_horizontal_variable)
            if (object%standard_variable%name==name) then
               id = create_external_variable_id(self%root,object%standard_variable)
               return
            end if
      end select
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
   id = create_external_variable_id(self%root,standard_variable)

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
   link => self%root%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_scalar_variable)
            if (object%name==name.or.get_safe_name(object%name)==name) then
               id = create_external_variable_id(object)
               return
            end if
      end select
      link => link%next
   end do

   ! Name not found among variable names. Now try standard names that are in use.
   link => self%root%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_scalar_variable)
            if (object%standard_variable%name==name) then
               id = create_external_variable_id(self%root,object%standard_variable)
               return
            end if
      end select
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
   id = create_external_variable_id(self%root,standard_variable)

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
   used = associated(id%p)

   end function fabm_is_bulk_variable_used
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Determine whether a bulk variable is available to biogeochemical models running in FABM.
!
! !INTERFACE:
   function fabm_bulk_variable_needs_values(id) result(required)
!
! !INPUT PARAMETERS:
   type(type_bulk_variable_id),intent(in) :: id
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
   required = associated(id%p)
   if (required) required = .not.associated(id%p%p)

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
   used = associated(id%p)

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
   used = associated(id%p)

   end function fabm_is_scalar_variable_used
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by its integer id.
!
! !INTERFACE:
   subroutine fabm_link_bulk_data_by_id(model,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                intent(inout) :: model
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
#if defined(DEBUG)&&_FABM_DIMENSION_COUNT_>0
   do i=1,size(shape(dat))
      if (size(dat,i)/=size(model%environment%diag,i)) then
         call fatal_error('fabm_link_bulk_data','dimensions of FABM domain and provided array do not match.')
      end if
   end do
#endif

   do i=1,size(id%alldata)
      id%alldata(i)%p%p => dat
   end do

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
! The variable is identified by its integer id.
!
! !INTERFACE:
   subroutine fabm_link_horizontal_data_by_id(model,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),                           intent(inout) :: model
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
#ifdef DEBUG
#ifndef _FABM_HORIZONTAL_IS_SCALAR_
   do i=1,size(shape(dat))
      if (size(dat,i)/=size(model%environment%diag,i)) then
         call fatal_error('fabm_link_horizontal_data','dimensions of FABM domain and provided array do not match.')
      end if
   end do
#endif
#endif

   do i=1,size(id%alldata)
      id%alldata(i)%p%p => dat
   end do

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
! data for the specified variable. The variable is identified by its integer id.
!
! !INTERFACE:
   subroutine fabm_link_scalar_by_id(model,id,dat)
!
! !INPUT PARAMETERS:
   class (type_model),            intent(inout) :: model
   type (type_scalar_variable_id),intent(inout) :: id
   real(rk),target,               intent(in)    :: dat
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
   do i=1,size(id%alldata)
      id%alldata(i)%p%p => dat
   end do

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
   ! Retrieve a pointer to the array holding the requested data.
   dat => self%environment%diag(_PREARG_LOCATION_DIMENSIONS_ self%diagnostic_variables(id)%globalid%write_index)

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
   ! Retrieve a pointer to the array holding the requested data.
   dat => self%environment%diag_hz(_PREARG_LOCATION_DIMENSIONS_HZ_ self%horizontal_diagnostic_variables(id)%globalid%write_index)

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
!EOP
!-----------------------------------------------------------------------
!BOC
#if defined(DEBUG)&&defined(_FABM_USE_1D_LOOP_)
   if (size(dy,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_do_rhs','Size of first dimension of dy should match loop length.')
#endif
   node => self%models%first
   do while (associated(node))
      if (node%model%check_conservation) then
         call node%model%do_check_conservation(_ARGUMENTS_ND_IN_,dy)
      else
         call node%model%do(_ARGUMENTS_ND_IN_,dy)
      end if
      node => node%next
   end do

   end subroutine fabm_do_rhs
!EOC


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
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
#if defined(DEBUG)&&defined(_FABM_USE_1D_LOOP_)
   if (size(pp,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_do_ppdd','Size of first dimension of pp should match loop length.')
   if (size(dd,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_do_ppdd','Size of first dimension of dd should match loop length.')
#endif

   node => self%models%first
   do while (associated(node))
      call node%model%do_ppdd(_ARGUMENTS_ND_IN_,pp,dd)
      node => node%next
   end do

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

   ! Allow individual models to check their state for their custom constraints, and to perform custom repairs.
   node => self%models%first
   do while (associated(node) .and. valid)
      call node%model%check_state(_ARGUMENTS_ND_IN_,repair,valid)
      if (.not. (valid .or. repair)) return
      node => node%next
   end do

   ! Finally check whether all state variable values lie within their prescribed [constant] bounds.
   ! This is always done, independently of any model-specific checks that may have been called above.

   ! Check boundaries for pelagic state variables specified by the models.
   ! If repair is permitted, this clips invalid values to the closest boundary.
   do ivar=1,size(self%state_variables)
      ! Shortcuts to variable information - this demonstrably helps the compiler (ifort).
      p = self%state_variables(ivar)%globalid%p
      minimum = self%state_variables(ivar)%minimum
      maximum = self%state_variables(ivar)%maximum

      _LOOP_BEGIN_EX_(self%environment)
         _GET_STATE_EX_(self%environment,p,value)
         if (value<minimum) then
            ! State variable value lies below prescribed minimum.
            valid = .false.
            if (.not.repair) then
               write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ',trim(self%state_variables(ivar)%name), &
                                                          & ' below minimum value ',minimum
               call log_message(err)
               return
            end if
            _SET_STATE_EX_(self%environment,p,minimum)
         elseif (value>maximum) then
            ! State variable value exceeds prescribed maximum.
            valid = .false.
            if (.not.repair) then
               write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',value,' of variable ',trim(self%state_variables(ivar)%name), &
                                                          & ' above maximum value ',maximum
               call log_message(err)
               return
            end if
            _SET_STATE_EX_(self%environment,p,maximum)
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
      p_hz = state_variables(ivar)%globalid%p
      minimum = state_variables(ivar)%minimum
      maximum = state_variables(ivar)%maximum

      _HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
       _GET_STATE_BEN_EX_(self%environment,p_hz,value)
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
            _SET_STATE_BEN_EX_(root%environment,p_hz,minimum)
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
            _SET_STATE_BEN_EX_(root%environment,p_hz,maximum)
         end if
      _HORIZONTAL_LOOP_END_
   end do
end subroutine

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-water exchange fluxes for all biogeochemical state variables.
! Positive values indicate fluxes into the ocean, negative values indicate fluxes
! out of the ocean. Units are tracer unit * m/s.
!
! !INTERFACE:
   subroutine fabm_do_surface_nosfstate(self _ARG_LOCATION_VARS_HZ_,flux_pel)
!
! !INPUT PARAMETERS:
      class (type_model),                          intent(inout) :: self
      _DECLARE_LOCATION_ARG_HZ_
!
! !INPUT/OUTPUT PARAMETERS:
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(out)   :: flux_pel
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
      type (type_model_list_node), pointer :: node
      real(rk),allocatable _DIMENSION_HORIZONTAL_SLICE_PLUS_1_ALLOCATABLE_ :: flux_sf
!-----------------------------------------------------------------------
!BOC
      allocate(flux_sf(_HORIZONTAL_SLICE_SHAPE_ 0))
      flux_pel = 0.0_rk
      node => self%models%first
      do while (associated(node))
         call node%model%do_surface(_ARGUMENTS_IN_HZ_,flux_pel,flux_sf)
         node => node%next
      end do

   end subroutine
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-water exchange fluxes for all biogeochemical state variables.
! Positive values indicate fluxes into the ocean, negative values indicate fluxes
! out of the ocean. Units are tracer unit * m/s.
!
! !INTERFACE:
   subroutine fabm_do_surface_sfstate(self _ARG_LOCATION_VARS_HZ_,flux_pel,flux_sf)
!
! !INPUT PARAMETERS:
      class (type_model),                          intent(inout) :: self
      _DECLARE_LOCATION_ARG_HZ_
!
! !INPUT/OUTPUT PARAMETERS:
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_,intent(out)   :: flux_pel,flux_sf
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
      type (type_model_list_node), pointer :: node
!-----------------------------------------------------------------------
!BOC
      flux_pel = 0.0_rk
      flux_sf = 0.0_rk
      node => self%models%first
      do while (associated(node))
         call node%model%do_surface(_ARGUMENTS_IN_HZ_,flux_pel,flux_sf)
         node => node%next
      end do

   end subroutine
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
!-----------------------------------------------------------------------
!BOC
   node => self%models%first
   do while (associated(node))
      call node%model%do_bottom(_ARGUMENTS_IN_HZ_,flux_pel,flux_ben)
      node => node%next
   end do

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
!-----------------------------------------------------------------------
!BOC
   node => self%models%first
   do while (associated(node))
      call node%model%do_bottom_ppdd(_ARGUMENTS_IN_HZ_,pp,dd,benthos_offset)
      node => node%next
   end do

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
#if defined(DEBUG)&&defined(_FABM_USE_1D_LOOP_)
   if (size(velocity,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_get_vertical_movement','Size of first dimension of velocity should match loop length.')
#endif

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
#if defined(DEBUG)&&defined(_FABM_USE_1D_LOOP_)
   if (size(extinction,1)/=fabm_loop_stop-fabm_loop_start+1) &
      call fatal_error('fabm_get_light_extinction','Size of first dimension of extinction should match loop length.')
#endif

   _LOOP_BEGIN_EX_(self%environment)
      _GET_EX_(self%extinction_data,extinction _INDEX_SLICE_)
   _LOOP_END_
   node => self%models%first
   do while (associated(node))
      call node%model%get_light_extinction(_ARGUMENTS_ND_IN_,extinction)
      node => node%next
   end do

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
!EOP
!-----------------------------------------------------------------------
!BOC
   node => self%models%first
   do while (associated(node))
      call node%model%get_light(_ARGUMENTS_VERT_IN_)
      node => node%next
   end do

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
   integer :: i
!-----------------------------------------------------------------------
!BOC
   do i=1,size(self%conserved_quantities)
      if (associated(self%conserved_quantities(i)%aggregate_variable%sum)) &
         call self%conserved_quantities(i)%aggregate_variable%sum%evaluate(_ARGUMENTS_ND_IN_)
      _LOOP_BEGIN_EX_(self%environment)
         _GET_EX_(self%conserved_quantities(i)%data,sums _INDEX_SLICE_PLUS_1_(i))
      _LOOP_END_
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
   do i=1,size(self%conserved_quantities)
      if (associated(self%conserved_quantities(i)%aggregate_variable%horizontal_sum)) &
         call self%conserved_quantities(i)%aggregate_variable%horizontal_sum%evaluate_horizontal(_ARGUMENTS_IN_HZ_)
      _HORIZONTAL_LOOP_BEGIN_EX_(self%environment)
         _GET_HORIZONTAL_EX_(self%conserved_quantities(i)%horizontal_data,sums _INDEX_HORIZONTAL_SLICE_PLUS_1_(i))
      _HORIZONTAL_LOOP_END_
   end do

   end subroutine fabm_get_horizontal_conserved_quantities
!EOC

subroutine fabm_state_to_conserved_quantities(self _ARG_LOCATION_ND_,y,sums)
   class (type_model),               intent(inout)  :: self
   _DECLARE_LOCATION_ARG_ND_
   real(rk) _DIMENSION_SLICE_PLUS_1_,intent(in)     :: y
   real(rk) _DIMENSION_SLICE_PLUS_1_,intent(out)    :: sums
   call self%root%state_to_conserved_quantities(_ARGUMENTS_ND_IN_,y,sums)
end subroutine

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
            expression%history(_PREARG_LOCATION_DIMENSIONS_ i) = expression%in%p
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
            + weight_right*expression%in%p
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
      expression%history(_PREARG_LOCATION_DIMENSIONS_ expression%n+1) = expression%in%p
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
            expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ i) = expression%in%p
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
            + weight_right*expression%in%p
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
      expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+1) = expression%in%p
      expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+3) = expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+2) &
         + frac_outside*(-expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%ioldest) + expression%history(_PREARG_LOCATION_DIMENSIONS_HZ_ expression%n+1))

      expression%last_time = t
   end subroutine

end subroutine

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
   if (.not.variable%state_indices%is_empty()) id%state_index = variable%state_indices%pointers(1)%p
   if (.not.variable%write_indices%is_empty()) id%write_index = variable%write_indices%pointers(1)%p
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
   if (.not.variable%state_indices%is_empty()) id%state_index = variable%state_indices%pointers(1)%p
   if (.not.variable%write_indices%is_empty()) id%write_index = variable%write_indices%pointers(1)%p
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
   if (.not.variable%state_indices%is_empty()) id%state_index = variable%state_indices%pointers(1)%p
   if (.not.variable%write_indices%is_empty()) id%write_index = variable%write_indices%pointers(1)%p
end function create_external_scalar_id

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
            if (object%standard_variable%compare(standard_variable).and..not.object%standard_variable%is_null()) then
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
            if (object%standard_variable%compare(standard_variable).and..not.object%standard_variable%is_null()) then
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
            if (object%standard_variable%compare(standard_variable).and..not.object%standard_variable%is_null()) then
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

subroutine classify_variables(self)
   class (type_model),intent(inout),target :: self

   type (type_link),pointer :: link

   type (type_state_variable_info),                pointer :: statevar
   type (type_horizontal_state_variable_info),     pointer :: hz_statevar
   type (type_diagnostic_variable_info),           pointer :: diagvar
   type (type_horizontal_diagnostic_variable_info),pointer :: hz_diagvar
   type (type_conserved_quantity_info),            pointer :: consvar
   class (type_internal_object),                   pointer :: object
   integer                                                 :: nstate,nstate_bot,nstate_surf,ndiag,ndiag_hz,ncons

   type (type_aggregate_variable),    pointer :: aggregate_variable
   type (type_set) :: dependencies,dependencies_hz,dependencies_scalar

   ! Count number of conserved quantities and allocate associated array.
   ncons = 0
   aggregate_variable => self%root%first_aggregate_variable
   do while (associated(aggregate_variable))
      ncons = ncons + 1
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
      ncons = ncons + 1
      consvar => self%conserved_quantities(ncons)
      consvar%standard_variable = aggregate_variable%standard_variable
      consvar%name = trim(consvar%standard_variable%name)
      consvar%units = trim(consvar%standard_variable%units)
      consvar%long_name = trim(consvar%standard_variable%name)
      consvar%path = trim(consvar%standard_variable%name)
      consvar%aggregate_variable => aggregate_variable
      object => self%root%find_object(trim(aggregate_variable%standard_variable%name))
      select type (object)
         class is (type_bulk_variable)
            call append_data_pointer(object%alldata,consvar%data)
      end select
      object => self%root%find_object(trim(aggregate_variable%standard_variable%name)//'_at_interfaces')
      select type (object)
         class is (type_horizontal_variable)
            call append_data_pointer(object%alldata,consvar%horizontal_data)
      end select
      aggregate_variable => aggregate_variable%next
   end do

   ! Get link to extinction variable.
   object => self%root%find_object(trim(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux%name))
   select type (object)
      class is (type_bulk_variable)
         call append_data_pointer(object%alldata,self%extinction_data)
   end select

   ! Count number of bulk variables in various categories.
   nstate = 0
   ndiag  = 0
   nstate_bot  = 0
   nstate_surf = 0
   ndiag_hz    = 0
   link => self%root%first_link
   do while (associated(link))
      if (link%owner) then
         select type (object => link%target)
            class is (type_bulk_variable)
               ! This is a variable owned by the model.
               if (.not.object%write_indices%is_empty()) then
                  ! This is a diagnostic variable.
                  ndiag = ndiag+1
                  call object%write_indices%set_value(ndiag)
               end if
               if (.not.object%state_indices%is_empty()) then
                  ! This is a state variable.
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
            class is (type_horizontal_variable)
               if (.not.object%write_indices%is_empty()) then
                  ! This is a diagnostic variable.
                  ndiag_hz = ndiag_hz+1
                  call object%write_indices%set_value(ndiag_hz)
               end if
               if (.not.object%state_indices%is_empty()) then
                  ! This is a state variable.
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
      end if
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
   link => self%root%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_bulk_variable)
            if (link%owner) then
               ! The model owns this variable (no external master variable has been assigned)
               ! Transfer variable information to the array that will be accessed by the host model.
               if (.not.object%write_indices%is_empty()) then
                  ndiag = ndiag + 1
                  diagvar => self%diagnostic_variables(ndiag)
                  call copy_variable_metadata(object,diagvar)
                  diagvar%globalid          = create_external_variable_id(object)
                  diagvar%standard_variable = object%standard_variable
                  diagvar%time_treatment    = output2time_treatment(diagvar%output)
                  call diagvar%properties%update(object%properties)
               end if

               if (object%presence==presence_internal.and..not.object%state_indices%is_empty()) then
                  nstate = nstate + 1
                  statevar => self%state_variables(nstate)
                  call copy_variable_metadata(object,statevar)
                  statevar%globalid                  = create_external_variable_id(object)
                  statevar%standard_variable         = object%standard_variable
                  statevar%initial_value             = object%initial_value
                  statevar%vertical_movement         = object%vertical_movement
                  statevar%no_precipitation_dilution = object%no_precipitation_dilution
                  statevar%no_river_dilution         = object%no_river_dilution
               end if
            end if
         class is (type_horizontal_variable)
            if (link%owner) then
               ! The model owns this variable (no external master variable has been assigned)
               ! Transfer variable information to the array that will be accessed by the host model.
               if (.not.object%write_indices%is_empty()) then
                  ndiag_hz = ndiag_hz + 1
                  hz_diagvar => self%horizontal_diagnostic_variables(ndiag_hz)
                  call copy_variable_metadata(object,hz_diagvar)
                  hz_diagvar%globalid          = create_external_variable_id(object)
                  hz_diagvar%standard_variable = object%standard_variable
                  hz_diagvar%time_treatment    = output2time_treatment(hz_diagvar%output)
               end if
               if (object%presence==presence_internal.and..not.object%state_indices%is_empty()) then
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
                  hz_statevar%globalid          = create_external_variable_id(object)
                  hz_statevar%standard_variable = object%standard_variable
                  hz_statevar%initial_value     = object%initial_value
               end if
            end if
      end select
      link => link%next
   end do

   ! Create lists of variables that may be provided by the host model.
   ! These lists include external dependencies, as well as the model's own state variables,
   ! which may be overridden by the host.
   link => self%root%first_link
   do while (associated(link))
      select type (object=>link%target)
         class is (type_bulk_variable)
            if (allocated(object%alldata).and. &
                .not.(object%presence==presence_external_optional.and..not.object%state_indices%is_empty())) then
               call dependencies%add(link%name)
               if (object%standard_variable%name/='') call dependencies%add(object%standard_variable%name)
            end if
         class is (type_horizontal_variable)
            if (allocated(object%alldata).and. &
                .not.(object%presence==presence_external_optional.and..not.object%state_indices%is_empty())) then
               call dependencies_hz%add(link%name)
               if (object%standard_variable%name/='') call dependencies_hz%add(object%standard_variable%name)
            end if
         class is (type_scalar_variable)
            if (allocated(object%alldata)) then
               call dependencies_scalar%add(link%name)
               if (object%standard_variable%name/='') call dependencies_scalar%add(object%standard_variable%name)
            end if
      end select
      link => link%next
   end do
   call dependencies%to_array(self%dependencies)
   call dependencies_hz%to_array(self%dependencies_hz)
   call dependencies_scalar%to_array(self%dependencies_scalar)
end subroutine classify_variables

subroutine copy_variable_metadata(internal_variable,external_variable)
   class (type_external_variable),intent(inout) :: external_variable
   class (type_internal_variable),intent(in)    :: internal_variable

   external_variable%name          = get_safe_name(internal_variable%name)
   external_variable%long_name     = internal_variable%long_name
   external_variable%units         = internal_variable%units
   external_variable%path          = internal_variable%name
   external_variable%minimum       = internal_variable%minimum
   external_variable%maximum       = internal_variable%maximum
   external_variable%missing_value = internal_variable%missing_value
   external_variable%output        = internal_variable%output
   call external_variable%properties%update(internal_variable%properties)
end subroutine

end module fabm

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

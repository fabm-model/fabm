!$Id$
#include "rmbm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RMBM --- Repository of Marine Biogeochemcial Models
!
! !INTERFACE:
   module rmbm
!
! !DESCRIPTION:
! This module encapsulates various specific biogeochemical models.
!
! ------------------------------------------------------------------------------------------
! How to register a new biogeochemical model:
! ------------------------------------------------------------------------------------------
!
! A new biogeochemical model only needs to be registered in this file
! in order to be usable by GOTM, MOM4 and any other physical drivers for RMBM.
!
! For the impatient: throughout this file, locations where additional bio models may be
! referenced are indicated by the comment "ADD_NEW_MODEL_HERE". You can search for this
! string to get started quickly.
!
! In detail:
!
! 1) Add a use statement that references your new model (cf. "use rmbm_npzd" below)
!
! 2) Define a unique integer identifier for your new model (cf. "npzd_id" below)
!
! 3) If your model uses model-specific data (e.g., parameter values), it is strongly recommended
!    that these are grouped in a (model-specific) derived type. This is a *requirement*
!    if you want to be able to run multiple instances of your model side-by-side!
!    (because in that case, any module-level variables would be shared between instances)
!
!    If your model groups data in a derived type, add an instance of this type as member to the
!    derived type "type_model" defined below.
!
! 4) Add the model as option to the "select" statements in the following subroutines:
!
!    "get_model_name" - to provide a short name (string) for your model.
!    "rmbm_init"      - to initialize your model: register prognostic and diagnostic variables,
!                       read in namelists, etc.
!    "rmbm_do_rhs"    - to provide the local temporal derivatives of your model.
!
! The following steps are optional:
!
! 5) If the sinking rate of any of the model variables varies in time and/or space, a subroutine
!    that provides the sinking rates (m/s) must be added as option to the "select" statement in
!    "rmbm_get_vertical_movement".
!
!    If this function is not provided, sinking rates will be set to the constant sinking
!    rates specified by the model upon its calls to register_state_variable.
!
! 6) If any of the model variables attenuate light, and this effect cannot be captured by a 
!    specific extinction coefficient for these model variables, a function that calculates the light
!    extinction coefficient (/m) from the current model state must be added as option to the "select"
!    statement in rmbm_get_light_extinction. This allows for complete customization of the bio
!    extinction.
!
!    If this function is not provided, the bio extinction will be calculated from the specific
!    extinction coefficients for each state variable as specified in the model information,
!    specifically, the member "specific_light_extinction" of the state variable information.
!    (note: specific extinction coefficients default to 0: no attenuation due to biogeochemical components) 
!
! 7) If (part of) the model state variable are composed of conserved quantities (energy and/or
!    chemical elements), a function that provides the sum of these quantities given the model
!    state can be added as option to the "select" statement in get_conserved_quantities_rmbm.
!    In that case, the totals of these conserved quantities can be diagnosed at run-time, and
!    included in the model output.
!
!    Note that your model should register any conserved quantities during initialization.
!
! ------------------------------------------------------------------------------------------
! How to use this library of biogeochemical models in a physical host model
! ------------------------------------------------------------------------------------------
!
! There are no explicit requirements on the physical model, except for the use of the
! interfaces documented below. However, in order to use RMBM efficiently, data for 
! environmental and biotic variables should be stored in a particular way by the host.
!
! Specifically:
!
! At any given point in time, all variable values (at all points in space) of any single
! environmental or biotic variable must be accessible through a single F90 pointer. For
! spatial models, all data for a single variable must thus be stored in a single array (e.g., a
! one-dimensional array for column models, or a three-dimensionsal array for global
! circulation models). Ultimately, the variable will be represented by a F90 array slice, which
! means that the variable for storing may use additional dimensions (e.g., for time),
! which are then set to some index chosen by the host, and ignored by RMBM.
!
! Note that this does *not* require that the values for all variables combined are
! stored in a single contiguous block of memory. In practice, it does mean that the
! values for a *single* variable are best stored in a contiguous block of memory.
!
! Variable values will be accessed through F90 pointers. Therefore, the FORTRAN variables 
! at the side of the host that contain environmental or biotic data must be able to serve
! as the target of a F90 pointer. That means that it should either have the "target" or
! "pointer" attribute, or be a member of a derived type with one of these attributes.
!
! The arrays on the side of the host could thus be defined with any of the following:
!
! for a 3D model:
!
! real,dimension(:,:,:),allocatable,target  :: temp,salt
! real,dimension(nx,ny,nz),         target  :: temp,salt
! real,dimension(:,:,:),            pointer :: temp,salt
!
! for a 1D model:
!
! real,dimension(:),allocatable,target  :: temp,salt
! real,dimension(nz),           target  :: temp,salt
! real,dimension(:),            pointer :: temp,salt
!
! They can also be stored as members of a derived type:
!
! type type_data
!    real,dimension(:,:,:),allocatable  :: temp,salt
! end type
! type (type_data),pointer :: dat
!
! Note that the type instance "dat" needs to have the "target" or "pointer" attribute,
! but the member variables "temp" and "salt" do not.
!
! Also, the data arrays may contain additional dimensions that should be set to a fixed
! index when used by RMBM, for instance:
!
! real,dimension(nx,ny,nz,3),target  :: temp,salt
!
! Later the last dimension can be ignored by providing RMBM with a slices such as
! temp(:,:,:,1) and salt(:,:,:,1)
!
! This can be needed for models that use an additional dimension for the time step index,
! e.g., mom4.
!
! If one or more environmental or biotic variables are not stored in the required manner,
! it will be upon the coupling layer between RMBM and its physical host to create arrays
! of the required structure for these variables, and make sure that the contained values
! are current upon each call to RMBM.

! A typical example of this would be a variable that stores bottom values on a 3D cartesian
! model grid: there the bottom index will vary in the horizontal, making it impossible to
! get at bottom values with a simple array slice (with fixed depth index). Before calls to
! RMBM, the coupling layer should therefore read bottom values from the physical host, and
! put them in the "bottom value" array that has been provided to RMBM.
!
! It may be clear that this is feasible for variables that have one
! or more dimensions fewer than the main grid, such as bottom/surface layers, but it will
! be (very!) computationally expensive if it is to be done for many variables that exist
! on the full grid.
!
! How to use the RMBM in a physical host:
!
! 1) Add a use statement that references this module, e.g.:
!
!    use rmbm
!
! 2) Define a pointer to an instance of the derived type "type_model":
!
!    type (model_type), pointer :: model
!
!    This instance will hold all information on the selected biogeochemical model(s),
!    including descriptive strings, state variable information, and parameter values.
!
! 3) Create the root model instance by calling rmbm_create_model. Typically, the root model
!    will be a model container that is created by calling create_model without arguments:
!
!    model => rmbm_create_model()
!
!    After the root model is created, it is empty. Child models can be added to it by
!    repeatedly calling rmbm_create_model, with the existing root model as parent argument:
!
!    childmodel => rmbm_create_model(1,parent=model)
!    childmodel => rmbm_create_model(2,parent=model)
!
!    This adds child models with identifiers 1 (currently NPZD) and 2 (currently Mnemiopsis)
!    to the root.
!
!    Note 1: it is not required to first construct an empty root model - if you know that you
!    will use one model from RMBM only, you can directly construct the desired model by calling
!    rmbm_create_model with a model identifier and no parent argument. This might render a slight
!    performance increase because RMBM does not need to enumerate the child models at each call.
!
!    Note 2: it is also possible to create a tree structure of models, nesting deeper than the
!    one level demonstarted above. By creating empty child models (omitting the model identifier
!    in rmbm_create_model), and then using these child models as parent in other calls to
!    rmbm_create_model, arbitrary model trees can be created. This is useful, because
!    namespaces at
!
! 4) Initialize the model tree by calling "init_rmbm" on the root model:
!
!    init_rmbm(model,nmlunit)
!
!    with "nmlunit" being a unit specifier (integer) of a file that has already been openend.
!    The RMBM models will read namelists from this file during initialization.
!    After initialization the caller is responsible for closing the file.
!
! 5) Provide the model with pointers to the arrays that will hold the current values of environmental
!    variables, and biotic state variables. These arrays normally are spatially explicit, i.e., they
!    contain all spatial dimensions.
!
!    For instance, a 3D model might have an array with temperature and salinity values allocated as
!    
!    real,dimension(nx,ny,nz),target :: temp,salt
!
!    These should be provided to RMBM as follows:
!
!    model%environment%temp => temp(:,:,:)
!    model%environment%salt => salt(:,:,:)
!
!    For an overview of all environmental variables that should be provided in this manner,
!    see the definition of the environment type ("type_environment") in rmbm_types.F90.
!
!    Similarly, if the n model state variables are stored by the host in variable
!
!    real,dimension(nx,ny,nz,n),target :: state
!
!    then arrays for all biological state variables must be provided as follows:
!
!    do ivar=1,n
!       model%state(ivar) => state(:,:,:,ivar)
!    end do
!
! 6) Access the model by the following subroutines:
!    - rmbm_do:                       to get local temporal derivatives
!    - rmbm_get_vertical_movement:    to get current vertcial movement rates for the state variables
!    - rmbm_get_light_extinction:     to get the combined light extinction coefficient due to
!                                     biogeochemical components
!    - rmbm_get_conserved_quantities: to get the sums of the conserved quantities described
!                                     by the model
!    - rmbm_get_surface_exchange:     to get updated fluxes over the air-sea interface
!    - rmbm_check_state:              to check the validity of current state variable values,
!                                     and repair (clip) these if desired
!
!    Additional information on the model (e.g. long names, units for its variables) is present in the
!    "info" member of the model, after the model has been initialized. For an overview of available metadata,
!    see the definition of the info type ("type_model_info") in rmbm_types.F90, as well as the definitions
!    of its contained types ("type_state_variable_info", "type_diagnostic_variable_info", etc.)
!
! !USES:
   use rmbm_types
   use rmbm_driver
!   
!  Reference specific biogeochemical models:
   use rmbm_npzd
   use rmbm_mnemiopsis
   use rmbm_co2sys
   ! ADD_NEW_MODEL_HERE
   
   implicit none
!   
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_model, rmbm_create_model, rmbm_init, rmbm_set_domain, rmbm_do, &
          rmbm_link_benthos_state_data,rmbm_link_state_data, &
          rmbm_link_data,rmbm_link_data_hz,rmbm_get_variable_id, &
          rmbm_get_diagnostic_data, rmbm_get_diagnostic_data_hz, &
          rmbm_check_state, rmbm_get_vertical_movement, rmbm_get_light_extinction, &
          rmbm_get_conserved_quantities, rmbm_get_surface_exchange, rmbm_do_benthos

#ifndef RMBM_MANAGE_DIAGNOSTICS
   public rmbm_link_diagnostic_data_hz,rmbm_link_diagnostic_data
#endif
!
! !PRIVATE DATA MEMBERS:

!  Identifiers for specific biogeochemical models.
!  Note: identifiers <=100 are reserved for models ported from the General Ocean Turbulence Model
   integer, parameter :: model_container_id = -1
   integer, parameter :: npzd_id            =  1
   integer, parameter :: carbonate_id       =  101
   integer, parameter :: mnemiopsis_id      =  102
   ! ADD_NEW_MODEL_HERE

! !PUBLIC TYPES:
!
!  Single generic biogeochemical model
   type type_model
   
      ! Model identifier, name (= variable name prefix), and metadata.
      integer                :: id
      type (type_model_info) :: info
      
      ! Pointers for linking to parent and child models
      ! (models can be linked to form a tree structure)
      type (type_model),pointer :: parent
      type (type_model),pointer :: firstchild
      type (type_model),pointer :: nextsibling
      
      ! Pointer to next non-container model.
      ! To speed up iterations over all models in the tree, non-container models are
      ! connected in a single linked list.
      ! For the root model, nextmodel points to the first non-container model in the list.
      ! For nested (non-root) container models, the pointer will remain disassociated(!)
      type (type_model),pointer :: nextmodel

      ! Derived types that belong to specific biogeochemical models.
      type (type_npzd)       :: npzd
      type (type_mnemiopsis) :: mnemiopsis
      type (type_co2_sys)    :: carbonate
      ! ADD_NEW_MODEL_HERE
      
      ! Pointers to the current state and environment.
      ! These are allocated upon initialization by RMBM,
      ! and should be set by the physical host model.
      type (type_environment),pointer :: environment
   end type type_model
   
   ! Arrays for integer identifiers and names of all available biogeochemical models.
   integer,          pointer,dimension(:) :: modelids   => null()
   character(len=64),pointer,dimension(:) :: modelnames => null()
!
! !PUBLIC INTERFACES:
!
!  Access to temporal derivatives
   interface rmbm_do
      module procedure rmbm_do_rhs
      module procedure rmbm_do_ppdd
      module procedure rmbm_do_rhs_1d
   end interface
   
   interface rmbm_link_data
      module procedure rmbm_link_data
      module procedure rmbm_link_data_char
   end interface

   interface rmbm_link_data_hz
      module procedure rmbm_link_data_hz
      module procedure rmbm_link_data_hz_char
   end interface
   
   interface rmbm_create_model
      module procedure rmbm_create_model_by_id
      module procedure rmbm_create_model_by_name
   end interface
!
! !REVISION HISTORY:!
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
   
   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:Register all biogeochemical models integer identifier and
!  short model name (letters, numbers and underscores only)
!
! !INTERFACE:
   subroutine register_models()
!
! !USES:
   implicit none
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
!-----------------------------------------------------------------------
!BOC
   call register_model(model_container_id,'')
   call register_model(npzd_id           ,'npzd')
   call register_model(mnemiopsis_id     ,'mnemiopsis')
   call register_model(carbonate_id      ,'co2sys')
   ! ADD_NEW_MODEL_HERE - required
   
   end subroutine register_models
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register integer identifier and short name (letters,
! numbers and underscores only) of a new biogeochemical model.
!
! !INTERFACE:
   subroutine register_model(id,name)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in)         :: id
   character(len=*),intent(in) :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
   integer :: oldcount,i
   integer,          pointer,dimension(:) :: modelids_new
   character(len=64),pointer,dimension(:) :: modelnames_new
   character(len=256) :: text
!
!-----------------------------------------------------------------------
!BOC
   ! Determine original model count.
   oldcount = 0
   if (associated(modelids)) oldcount = ubound(modelids,1)
   
   ! Create new arrays for identifiers and names, one longer to accomodate the new model.
   allocate(modelids_new  (oldcount+1))
   allocate(modelnames_new(oldcount+1))
   
   if (associated(modelids)) then
      ! First check whether the provided model identifier or name are not in use yet.
      do i=1,ubound(modelids,1)
         if (modelids(i)==id) then
            write (text,fmt='(a,i4,a)') 'model identifier ',id,' has already been registered.'
            call fatal_error('rmbm::register_model_name',text)
         end if
         if (modelnames(i)==name) &
            call fatal_error('rmbm::register_model_name','model name "'//trim(name)//'" has already been registered.')
      end do
      
      ! Copy the old identifiers and names to the new extended array.
      modelids_new  (1:oldcount) = modelids(:)
      modelnames_new(1:oldcount) = modelnames(:)
      
      ! Deallocate the old arrays.
      deallocate(modelids)
      deallocate(modelnames)
   end if
   
   ! Add the new identifier and name.
   modelids_new  (oldcount+1) = id
   modelnames_new(oldcount+1) = name
   
   ! Assign the new arrays to the module-level pointers.
   modelids   => modelids_new
   modelnames => modelnames_new
   
   end subroutine register_model
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the short model name (letters, numbers and underscores
! only) for the biogeochemcial model with the specified integer identifier.
!
! !INTERFACE:
   function get_model_name(id) result(name)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer, intent(in) :: id
   character(len=64)   :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   character(len=256) :: text
   integer :: i
!-----------------------------------------------------------------------
!BOC
   if (.not.associated(modelids)) call register_models()
   do i=1,ubound(modelids,1)
      if (modelids(i)==id) then
         name = modelnames(i)
         return
      end if
   end do
   write (text,fmt='(i4,a)') id,' is not a valid model identifier.'
   call fatal_error('rmbm::get_model_name',text)
   
   end function get_model_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the integer model identifier from the specified short
! model name (letters, numbers and underscores only)
!
! !INTERFACE:
   function get_model_id(name) result(id)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   character(len=*),intent(in) :: name
   integer                     :: id
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   integer :: i
!-----------------------------------------------------------------------
!BOC
   if (.not.associated(modelnames)) call register_models()
   do i=1,ubound(modelnames,1)
      if (modelnames(i)==name) then
         id = modelids(i)
         return
      end if
   end do
   call fatal_error('rmbm::get_model_id',trim(name)//' is not a valid model name.')
   
   end function get_model_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Create a new model. This will be a specific model if the
! model identifier is provided, or a container of child models if the
! identifier is omitted.
!
! !INTERFACE:
   function rmbm_create_model_by_id(modelid,parent) result(model)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   integer,         optional,         intent(in)    :: modelid
   type (type_model),pointer,optional               :: parent
   type (type_model),pointer                        :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
   type (type_model),pointer                        :: curmodel
!-----------------------------------------------------------------------
!BOC
   allocate(model)

   ! Set the model identifier.
   if (present(modelid)) then
      model%id = modelid
   else
      model%id = model_container_id
   end if

   ! Initialize model info
   call init_model_info(model%info)

   ! Set the model name
   model%info%name = get_model_name(model%id)

   ! Make sure the pointers to parent and sibling models are explicitly dissociated.
   nullify(model%parent)
   nullify(model%firstchild)
   nullify(model%nextsibling)
   nullify(model%nextmodel)
   
   ! Make sure the pointers to the current environment are explicitly dissociated.
   nullify(model%environment)

   ! Connect to parent container if provided.
   if (present(parent)) then
      ! Make sure the provided parent model is a container.
      if (parent%id.ne.model_container_id) &
         call fatal_error('rmbm_create_model','A child model can only be added to a container, not to an existing model.')

      ! Link parent container to child.
      if (associated(parent%firstchild)) then
         ! The target container already contains one or more children.
         ! Find the last, then append the new model to the list.
         curmodel => parent%firstchild
         do while (associated(curmodel%nextsibling))
            curmodel => curmodel%nextsibling
         end do
         curmodel%nextsibling => model
         curmodel%info%nextsibling => model%info
      else
         ! The target container does not contain any children yet.
         ! Add the current model as first child.
         parent%firstchild => model
         parent%info%firstchild => model%info
      end if
      
      ! Link the child model to its parent container.
      model%parent => parent
      model%info%parent => parent%info

      if (model%id.ne.model_container_id) then
         ! This model is not a container.
         ! Add the model to the flattened list of non-container models,
         ! which is used at runtime to iterate over all models without needing recursion.
      
         ! First non-container model will be pointed to by the root of the tree - find it.
         curmodel => parent
         do while (associated(curmodel%parent))
            curmodel => curmodel%parent
         end do
         
         ! Find the last model in the flattened list of non-container models.
         ! (first model is pointed to by the root of the tree)
         do while (associated(curmodel%nextmodel))
            curmodel => curmodel%nextmodel
         end do
         
         ! Add current model to the flattened list of non-container models.
         curmodel%nextmodel => model
      end if

   end if
   
   end function rmbm_create_model_by_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Create a new model by name. This will be a specific model if the
! model name is provided, or a container of child models if the
! identifier is omitted.
!
! !INTERFACE:
   function rmbm_create_model_by_name(modelname,parent) result(model)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   character(len=*),    intent(in)    :: modelname
   type (type_model),pointer,optional :: parent
   type (type_model),pointer          :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   if (present(parent)) then
      model => rmbm_create_model_by_id(get_model_id(modelname),parent=parent)
   else
      model => rmbm_create_model_by_id(get_model_id(modelname))
   end if
   
   end function rmbm_create_model_by_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the biogeochemical model tree.
!
! !INTERFACE:
   subroutine rmbm_init(root,nmlunit)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),target,               intent(inout) :: root
   integer,                                intent(in)    :: nmlunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  logical                    :: isopen
  integer                    :: ivar
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Check whether we are operating on the root of a model tree.
   if (associated(root%parent)) call fatal_error('rmbm_init','rmbm_init can only be called on the root of a model tree.')

   ! Check whether the unit provided by the host actually refers to an open file.
   inquire(nmlunit,opened=isopen)
   if (.not.isopen) call fatal_error('rmbm_init','input configuration file has not been opened yet!')
   
   call init_model(root,nmlunit)

   ! Allocate arrays for storage of (references to) data.
   allocate(root%environment)

   ! Set all pointers to external data to dissociated.
   allocate(root%environment%var_hz(ubound(root%info%dependencies_hz,1)))
   allocate(root%environment%var   (ubound(root%info%dependencies,   1)))
   do ivar=1,ubound(root%environment%var_hz,1)
      nullify(root%environment%var_hz(ivar)%data)
   end do
   do ivar=1,ubound(root%environment%var,1)
      nullify(root%environment%var(ivar)%data)
   end do

#ifdef RMBM_SINGLE_STATE_VARIABLE_ARRAY
   nullify(root%environment%state)
   nullify(root%environment%state_ben)
#else
   allocate(root%environment%state    (ubound(root%info%state_variables,    1)))
   allocate(root%environment%state_ben(ubound(root%info%state_variables_ben,1)))
   do ivar=1,ubound(root%environment%state,1)
      nullify(root%environment%state(ivar)%data)
   end do
   do ivar=1,ubound(root%environment%state_ben,1)
      nullify(root%environment%state_ben(ivar)%data)
   end do
#endif

   ! Transfer pointer to environment to all child models.
   call set_model_data_members(root,root%environment)

   end subroutine rmbm_init
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the provided biogeochemical model
!
! !INTERFACE:
   recursive subroutine init_model(model,nmlunit)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),target,               intent(inout) :: model
   integer,                                intent(in)    :: nmlunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer                    :: count,ownindex
  character(len= 64)         :: modelname
  type (type_model), pointer :: curchild,curchild2
  logical                    :: isopen
  logical,parameter          :: alwayspostfixindex=.false.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Retrieve model name based on its integer identifier.
   modelname = get_model_name(model%id)
   call log_message('Initializing biogeochemical model '//trim(modelname)//'...')
   
   ! Allow the selected model to initialize
   select case (model%id)
      case (npzd_id)
         call npzd_init(model%npzd,model%info,nmlunit)
      case (mnemiopsis_id)
         call mnemiopsis_init(model%mnemiopsis,model%info,nmlunit)
      case (carbonate_id)
         call co2sys_init(model%carbonate,model%info,nmlunit)
      ! ADD_NEW_MODEL_HERE - required
      
      case (model_container_id)
         ! This is a container of models. Loop over all children and allow each to initialize.
         curchild => model%firstchild
         do while (associated(curchild))
            ! Find the number of the child model. This will be used in the prefix of model variable names.
            count = 0
            curchild2 => model%firstchild
            do while (associated(curchild2))
               if (curchild2%info%name.eq.curchild%info%name) then
                  count = count + 1
                  if (associated(curchild,curchild2)) ownindex = count
               end if
               curchild2 => curchild2%nextsibling
            end do

            ! Create prefixes for names of variables of the child model.
            if (alwayspostfixindex .or. count>1) then
               write (unit=curchild%info%nameprefix,     fmt='(a,i2.2,a)') trim(curchild%info%name),ownindex,'_'
               write (unit=curchild%info%longnameprefix, fmt='(a,i2.2,a)') trim(curchild%info%name),ownindex,' '
            else
               curchild%info%nameprefix     = trim(curchild%info%name)//'_'
               curchild%info%longnameprefix = trim(curchild%info%name)//' '
            end if

            ! Initialize child model.
            call init_model(curchild,nmlunit)
            
            ! Move to next child model.
            curchild => curchild%nextsibling
         end do
         
      case default
         call fatal_error('rmbm_init','no valid biogeochemical model specified!')
         
   end select
   
   call log_message('model '//trim(modelname)//' initialized successfully.')

   ! Check whether the unit provided by the host has not been closed by the biogeochemical model.
   inquire(nmlunit,opened=isopen)
   if (.not.isopen) call fatal_error('rmbm_init','input configuration file was closed by model '//trim(modelname)//'.')

   end subroutine init_model
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tell RMBM about the extents of the spatial domain.
! This allows it to create spatially-explicit arrays for internal use.
!
! !INTERFACE:
   subroutine rmbm_set_domain(root,LOCATION)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),target,               intent(inout) :: root
   LOCATION_TYPE,                          intent(in)    :: LOCATION
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer                    :: i
!EOP
!-----------------------------------------------------------------------
!BOC
#ifdef RMBM_MANAGE_DIAGNOSTICS
   ! RMBM will manage and store current values of diagnostic variables.

   ! Allocate arrays for diagnostic variables defined on the full domain and on horizontonal slices.
   allocate(root%environment%diag   (ubound(root%info%diagnostic_variables,1) ARG_LOCATION))
   allocate(root%environment%diag_hz(ubound(root%info%diagnostic_variables_hz,1) ARG_LOCATION_HZ))

   ! Initialize diagnostic variables to zero.
   root%environment%diag    = _ZERO_
   root%environment%diag_hz = _ZERO_
   
   ! If diagnostic variables also appear as dependency, send the corresponding array slice for generic read-only access.
   do i=1,ubound(root%info%diagnostic_variables,1)
      if (root%info%diagnostic_variables(i)%dependencyid.ne.id_not_used) &
         call rmbm_link_data(root,root%info%diagnostic_variables(i)%dependencyid,root%environment%diag(i ARG_LOCATION_DIMENSIONS))
   end do
   do i=1,ubound(root%info%diagnostic_variables_hz,1)
      if (root%info%diagnostic_variables_hz(i)%dependencyid.ne.id_not_used) &
         call rmbm_link_data_hz(root,root%info%diagnostic_variables_hz(i)%dependencyid,root%environment%diag_hz(i ARG_LOCATION_DIMENSIONS_HZ))
   end do
#endif
end subroutine rmbm_set_domain

recursive subroutine set_model_data_members(model,environment)
   type (type_model),             intent(inout) :: model
   type (type_environment),target,intent(in)    :: environment
   
   type (type_model), pointer :: curchild
   
   model%environment => environment
   curchild => model%firstchild
   do while (associated(curchild))
      call set_model_data_members(curchild,environment)
      curchild => curchild%nextsibling
   end do
end subroutine set_model_data_members

function rmbm_get_variable_id(model,name,shape) result(id)
   type (type_model),             intent(in)  :: model
   character(len=*),              intent(in)  :: name
   integer,                       intent(in)  :: shape
   integer                                    :: id
   character(len=64),pointer                  :: source(:)
   
   select case (shape)
      case (shape_hz)
         source => model%info%dependencies_hz
      case (shape_full)
         source => model%info%dependencies
   end select

   if (associated(source)) then
      do id=1,ubound(source,1)
         if (source(id)==name) return
      end do
   end if
   
   id = id_not_used
end function rmbm_get_variable_id

subroutine rmbm_link_data(model,id,dat)
   type (type_model),                       intent(inout) :: model
   integer,                                 intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS,target,intent(in)    :: dat
   
   model%environment%var(id)%data => dat
end subroutine rmbm_link_data

subroutine rmbm_link_data_char(model,name,dat)
   type (type_model),                       intent(inout) :: model
   character(len=*),                        intent(in)    :: name
   REALTYPE ATTR_LOCATION_DIMENSIONS,target,intent(in)    :: dat
   
   integer                                                :: id
   
   id = rmbm_get_variable_id(model,name,shape_full)
   if (id.ne.id_not_used) call rmbm_link_data(model,id,dat)
end subroutine rmbm_link_data_char

subroutine rmbm_link_data_hz(model,id,dat)
   type (type_model),                          intent(inout) :: model
   integer,                                    intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS_HZ,target,intent(in)    :: dat
   
   model%environment%var_hz(id)%data => dat
end subroutine rmbm_link_data_hz

subroutine rmbm_link_data_hz_char(model,name,dat)
   type (type_model),                          intent(inout) :: model
   character(len=*),                           intent(in)    :: name
   REALTYPE ATTR_LOCATION_DIMENSIONS_HZ,target,intent(in)    :: dat
   
   integer                                                   :: id
   
   id = rmbm_get_variable_id(model,name,shape_hz)
   if (id.ne.id_not_used) call rmbm_link_data_hz(model,id,dat)
end subroutine rmbm_link_data_hz_char

#ifdef RMBM_SINGLE_STATE_VARIABLE_ARRAY

subroutine rmbm_link_state_data(model,dat)
   type (type_model),                                intent(inout) :: model
   REALTYPE ATTR_LOCATION_DIMENSIONS_PLUS_ONE,target,intent(in)    :: dat
   
   integer                                                         :: id
   
   if (ubound(dat,1).ne.ubound(model%info%state_variables,1)) &
      call fatal_error('rmbm::rmbm_link_state_data','The length of the first dimension of the state variable&
      & array does not match the number of state variables.')
   
   model%environment%state => dat
   do id=1,ubound(model%info%state_variables,1)
      if (model%info%state_variables(id)%dependencyid.ne.id_not_used) &
         call rmbm_link_data(model,model%info%state_variables(id)%dependencyid,dat(id ARG_LOCATION_DIMENSIONS))
   end do
end subroutine rmbm_link_state_data

subroutine rmbm_link_benthos_state_data(model,dat)
   type (type_model),                                   intent(inout) :: model
   REALTYPE ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE,target,intent(in)    :: dat

   integer                                                            :: id
   
   if (ubound(dat,1).ne.ubound(model%info%state_variables_ben,1)) &
      call fatal_error('rmbm::rmbm_link_benthos_state_data','The length of the first dimension of the benthic&
      & state variable array does not match the number of state variables.')

   model%environment%state_ben => dat
   do id=1,ubound(model%info%state_variables_ben,1)
      if (model%info%state_variables_ben(id)%dependencyid.ne.id_not_used) &
         call rmbm_link_data_hz(model,model%info%state_variables_ben(id)%dependencyid,dat(id ARG_LOCATION_DIMENSIONS_HZ))
   end do
end subroutine rmbm_link_benthos_state_data

#else

subroutine rmbm_link_state_data(model,id,dat)
   type (type_model),                       intent(inout) :: model
   integer,                                 intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS,target,intent(in)    :: dat
   
   model%environment%state(id)%data => dat
   if (model%info%state_variables(id)%dependencyid.ne.id_not_used) &
      call rmbm_link_data(model,model%info%state_variables(id)%dependencyid,dat)
end subroutine rmbm_link_state_data

subroutine rmbm_link_benthos_state_data(model,id,dat)
   type (type_model),                          intent(inout) :: model
   integer,                                    intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS_HZ,target,intent(in)    :: dat
   
   model%environment%state_ben(id)%data => dat
   if (model%info%state_variables_ben(id)%dependencyid.ne.id_not_used) &
      call rmbm_link_data_hz(model,model%info%state_variables_ben(id)%dependencyid,dat)
end subroutine rmbm_link_benthos_state_data

#endif

function rmbm_get_diagnostic_data(model,id) result(dat)
   type (type_model),                       intent(inout) :: model
   integer,                                 intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS,pointer              :: dat
   
#ifdef RMBM_MANAGE_DIAGNOSTICS
   dat => model%environment%diag(id ARG_LOCATION_DIMENSIONS)
#else
   dat => model%environment%var(model%info%diagnostic_variables(id)%dependencyid)%data
#endif
end function rmbm_get_diagnostic_data

function rmbm_get_diagnostic_data_hz(model,id) result(dat)
   type (type_model),                          intent(inout) :: model
   integer,                                    intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS_HZ,pointer              :: dat
   
#ifdef RMBM_MANAGE_DIAGNOSTICS
   dat => model%environment%diag_hz(id ARG_LOCATION_DIMENSIONS_HZ)
#else
   dat => model%environment%var_hz(model%info%diagnostic_variables_hz(id)%dependencyid)%data
#endif
end function rmbm_get_diagnostic_data_hz

#ifndef RMBM_MANAGE_DIAGNOSTICS

subroutine rmbm_link_diagnostic_data(model,id,dat)
   type (type_model),                       intent(inout) :: model
   integer,                                 intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS,target,intent(in)    :: dat
   
   call rmbm_link_data(model,model%info%diagnostic_variables(id)%dependencyid,dat)
end subroutine rmbm_link_diagnostic_data

subroutine rmbm_link_diagnostic_data_hz(model,id,dat)
   type (type_model),                          intent(inout) :: model
   integer,                                    intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS_HZ,target,intent(in)    :: dat
   
   call rmbm_link_data_hz(model,model%info%diagnostic_variables_hz(id)%dependencyid,dat)
end subroutine rmbm_link_diagnostic_data_hz

#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the biogeochemical
! model tree.
!
! !INTERFACE:
   subroutine rmbm_do_rhs(root,LOCATION,dy)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),      intent(inout) :: root
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, dimension(:), intent(inout) :: dy
!
! !LOCAL PARAMETERS:
   type (type_model), pointer            :: curmodel
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   if (associated(root%parent)) &
      call fatal_error('rmbm_do_rhs','rmbm_do_rhs may only be called on the root of the model tree, or on non-container models.')

   curmodel => root%nextmodel
   do while (associated(curmodel))
      select case (curmodel%id)
         case (npzd_id)
            call npzd_do(curmodel%npzd,dy RMBM_ARGS_IN)
         case (mnemiopsis_id)
            call mnemiopsis_do(curmodel%mnemiopsis,dy RMBM_ARGS_IN)
         case (carbonate_id)
            call co2sys_do(curmodel%carbonate,dy RMBM_ARGS_IN)
         ! ADD_NEW_MODEL_HERE - required, unless the model provides production/destruction
         ! matrices instead of a temporal derivative vector. In that case, add the model to
         ! do_ppdd_to_rhs.
         case default
           call do_ppdd_to_rhs(curmodel,LOCATION,dy)
      end select
      curmodel => curmodel%nextmodel
   end do

   end subroutine rmbm_do_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the biogeochemical
! model tree, by caluclating and summing production/destruction matrices.
!
! !INTERFACE:
   subroutine do_ppdd_to_rhs(root,LOCATION,dy)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),      intent(inout) :: root
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, dimension(:), intent(inout) :: dy
!
! !LOCAL PARAMETERS:
   REALTYPE,allocatable                  :: pp(:,:),dd(:,:)
   integer                               :: i,j
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   ! The model does not provide the temporal derivatives itself.
   ! In that case, it provides production/destruction matrices instead.
   ! Retrieve those, and calculate the temporal derivatives based on their contents.
   !
   ! NB pp/dd are allocated on the spot here, which is really expensive.
   ! However, this could only be improved by making RMBM aware of the full model domain,
   ! and making biogeochemical models tell RMB whther they provide pp/dd. Only then can RMBM
   ! allocate pp/dd for the full domain. (allocating the arrays once for a local
   ! point in space seems attractive, but that would rule out parallelization)
   allocate(pp(1:ubound(dy,1),1:ubound(dy,1)))
   allocate(dd(1:ubound(dy,1),1:ubound(dy,1)))
   pp = _ZERO_
   dd = _ZERO_
   select case (root%id)
      case (npzd_id)
         call npzd_do_ppdd(root%npzd,pp,dd RMBM_ARGS_IN)
      ! ADD_NEW_MODEL_HERE - only needed if the model is not added to rmbm_do_rhs.
      case default
         call fatal_error('rmbm::do_ppdd_to_rhs','model '//trim(root%info%name)//' does not provide a&
         & subroutine for calculating the temporal derivative vector.')
   end select
   do i=1,ubound(dy,1)
      do j=1,ubound(dy,1)
         dy(i) = dy(i) + pp(i,j)-dd(i,j)
      end do
   end do
   deallocate(pp)
   deallocate(dd)
   
   end subroutine do_ppdd_to_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives in 1D
!
! !INTERFACE:
   subroutine rmbm_do_rhs_1d(root ARG_LOCATION_1DLOOP,istart,istop,dy)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: root
   integer,                intent(in)    :: istart,istop
   DEFINE_LOCATION_1DLOOP
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: dy(:,:)
!
! !LOCAL PARAMETERS:
   type (type_model), pointer            :: curmodel
   LOCATION_TYPE                         :: VARIABLE_1DLOOP
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   curmodel = root%nextmodel
   do while (associated(curmodel))
      select case (curmodel%id)
         ! ADD_NEW_MODEL_HERE - optional
         case default
            do VARIABLE_1DLOOP=istart,istop
               call rmbm_do_rhs(curmodel,LOCATION,dy(VARIABLE_1DLOOP,:))
            end do
      end select
      curmodel = curmodel%nextmodel
   end do

   end subroutine rmbm_do_rhs_1d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the biogeochemical
! model tree in the form of production and destruction matrices.
!
! !INTERFACE:
   recursive subroutine rmbm_do_ppdd(root,LOCATION,pp,dd)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),      intent(inout) :: root
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(:,:),intent(inout) :: pp,dd
!
! !LOCAL PARAMETERS:
   REALTYPE,allocatable                  :: dy(:)
   integer                               :: i
   type (type_model), pointer            :: curmodel
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   curmodel => root%nextmodel
   do while (associated(curmodel))
      select case (curmodel%id)
         case (npzd_id)
            call npzd_do_ppdd(curmodel%npzd,pp,dd RMBM_ARGS_IN)
         ! ADD_NEW_MODEL_HERE - optional
         case default
            ! The model does not provide production/destruction matrices itself.
            ! In that case, it provides temporal derivatives instead.
            ! Retrieve those, and create degenerate production/destruction from their contents.
            allocate(dy(ubound(pp,1)))
            dy = _ZERO_
            call rmbm_do_rhs(curmodel,LOCATION,dy)
            do i=1,ubound(pp,1)
               pp(i,i) = pp(i,i) + dy(i)
            end do
            deallocate(dy)
      end select
      curmodel => curmodel%nextmodel
   end do

   end subroutine rmbm_do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks whether the current state is valid, and repairs [clips]
! invalid state variables if requested and possible.
!
! !INTERFACE:
   function rmbm_check_state(root,LOCATION,repair) result(valid)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),      intent(inout) :: root
   LOCATION_TYPE,          intent(in)    :: LOCATION
   logical,                intent(in)    :: repair

   logical                               :: valid
!
! !LOCAL PARAMETERS:
   integer                               :: i
   type (type_model), pointer            :: curmodel
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   valid = .true.
   
   ! Allow models to perform custom repairs.
   curmodel => root%nextmodel
   do while (associated(curmodel) .and. valid)
      select case (curmodel%id)
         ! ADD_NEW_MODEL_HERE - optional
         ! If a subroutine or function is provided:
         ! Set "valid" to .false. if the state variable values are invalid,
         ! repair them if they are invalid and "repair" is .true. (but still keep "valid" set to .false.)
      end select
      curmodel => curmodel%nextmodel
   end do
   
   ! If the present values are invalid and repair is not permitted, we are done.
   if (.not. (valid .or. repair)) return

   ! Check absolute variable boundaries specified by the models.
   ! If repair is permitted, this clips invalid values to the closest boundary.
   do i=1,ubound(root%info%state_variables,1)
      if (root%_GET_STATE_(root%info%state_variables(i)%globalid)<root%info%state_variables(i)%minimum) then
         valid = .false.
         if (repair) root%_GET_STATE_(root%info%state_variables(i)%globalid) =&
             root%info%state_variables(i)%minimum
      elseif (root%_GET_STATE_(root%info%state_variables(i)%globalid)> &
               root%info%state_variables(i)%maximum) then
         valid = .false.
         if (repair) root%_GET_STATE_(root%info%state_variables(i)%globalid) =&
             root%info%state_variables(i)%maximum
      end if
   end do

   end function rmbm_check_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-water exchange fluxes for all biogeochemical state variables.
! Positive values indicate fluxes into the ocean, negative values indicate fluxes
! out of the ocean. Units are tracer unit * m/s.
!
! !INTERFACE:
   subroutine rmbm_get_surface_exchange(root,LOCATION,flux)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),      intent(inout) :: root
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: flux(:)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer            :: curmodel
!-----------------------------------------------------------------------
!BOC
   curmodel => root%nextmodel
   do while (associated(curmodel))
      select case (curmodel%id)
         case (carbonate_id)
            call co2sys_get_surface_exchange(curmodel%carbonate,flux RMBM_ARGS_IN)
         ! ADD_NEW_MODEL_HERE - optional
      end select
      curmodel => curmodel%nextmodel
   end do

   end subroutine rmbm_get_surface_exchange
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process interaction between benthos and bottom layer of the
! pelagic. This calculates the fluxes into all bottom pelagic and benthic variables,
! in variable units * m/s (similar to surface fluxes). Positive values denote state
! variable increases, negative values state variable decreases, for both pelagic and
! benthic variables. Horizontal-only diagnostic variables can also be set in this routine.
!
! !INTERFACE:
   subroutine rmbm_do_benthos(root,LOCATION,flux_ben,flux_pel)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: root
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(:),  intent(inout) :: flux_ben,flux_pel
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer            :: curmodel
!-----------------------------------------------------------------------
!BOC
   curmodel => root%nextmodel
   do while (associated(curmodel))
      select case (curmodel%id)
         ! ADD_NEW_MODEL_HERE - optional
      end select
      curmodel => curmodel%nextmodel
   end do

   end subroutine rmbm_do_benthos
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical movement rates (m/s) for the bio state variables.
! Note that negative values indidcate movement towards the bottom, e.g., sinking,
! and positive values indicate movemment towards the surface, e.g., floating.
!
! !INTERFACE:
   subroutine rmbm_get_vertical_movement(root,LOCATION,vertical_movement)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: root
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: vertical_movement(:)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer            :: curmodel
   integer                               :: i
!-----------------------------------------------------------------------
!BOC
   curmodel => root%nextmodel
   do while (associated(curmodel))
      select case (curmodel%id)
         ! ADD_NEW_MODEL_HERE - optional
         case default
            ! Default: use the constant sinking rates specified in state variable properties.
            do i=1,ubound(curmodel%info%state_variables,1)
               vertical_movement(curmodel%info%state_variables(i)%globalid) =&
                   curmodel%info%state_variables(i)%vertical_movement
            end do
      end select
      curmodel => curmodel%nextmodel
   end do

   end subroutine rmbm_get_vertical_movement
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   function rmbm_get_light_extinction(root,LOCATION) result(extinction)
!
! !INPUT PARAMETERS:
   type (type_model), intent(in) :: root
   LOCATION_TYPE,     intent(in) :: LOCATION
   REALTYPE                      :: extinction
   
   integer                       :: i
   type (type_model), pointer    :: curmodel
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   curmodel => root%nextmodel
   do while (associated(curmodel))
      select case (curmodel%id)
         case (npzd_id)
            extinction = npzd_get_light_extinction(curmodel%npzd,curmodel%environment,LOCATION)
         ! ADD_NEW_MODEL_HERE - optional
         case default
            ! Default: use constant specific light extinction values specified in the state variable properties
            extinction = _ZERO_
            do i=1,ubound(curmodel%info%state_variables,1)
               if (curmodel%info%state_variables(i)%specific_light_extinction.ne._ZERO_) &
                  extinction = extinction + root%_GET_STATE_(curmodel%info%state_variables(i)%globalid) &
                                    & * curmodel%info%state_variables(i)%specific_light_extinction
            end do
      end select
      curmodel => curmodel%nextmodel
   end do

   end function rmbm_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities
!
! !INTERFACE:
   subroutine rmbm_get_conserved_quantities(root,LOCATION,sums)
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
   type (type_model), intent(inout)   :: root
   LOCATION_TYPE,     intent(in)      :: LOCATION
   REALTYPE,          intent(inout)   :: sums(:)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer    :: curmodel

!-----------------------------------------------------------------------
!BOC
   curmodel => root%nextmodel
   do while (associated(curmodel))
      select case (curmodel%id)
         case (npzd_id)
            call npzd_get_conserved_quantities(curmodel%npzd,sums RMBM_ARGS_IN)
         ! ADD_NEW_MODEL_HERE - optional
         case default
            ! Default: the model does not describe any conserved quantities.
            if (ubound(curmodel%info%conserved_quantities,1).gt.0) &
               call fatal_error('rmbm_get_conserved_quantities','the model specifies that it describes one or more conserved &
                    &quantities, but a function that provides sums of these quantities has not been specified.')
      end select
      curmodel => curmodel%nextmodel
   end do

   end subroutine rmbm_get_conserved_quantities
!EOC

end module rmbm

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

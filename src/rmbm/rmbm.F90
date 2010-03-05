!$Id: bio_0d_gen.F90,v 1.8 2009-05-11 13:57:31 jorn Exp $
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
! A new biogeochemical models only need to be registered in this file
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
! 2) Define an integer identifier for your new model (cf. "npzd_id" below)
!
! 3) If your model uses model-specific data (e.g. parameter values) that are grouped in a
!    (model-specific) derived type, add an instance of this type as member to the derived
!    type "type_model" defined below.
!
!    Note: if you want to be able to run multiple instances of your model
!    side-by-side, your model cannot use any module-level variables.
!    Grouping the model parameters in a derived type, then, is a *requirement*.
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
! 6) If any of the model variables attenuate light, a function that calculates the light
!    extinction coefficient (/m) from the current model state must be added as option to the "select"
!    statement in rmbm_get_bio_extinction. This allows for complete customization of the bio
!    extinction.
!
!    If this function is not provided, the bio extinction will be calculated from the specific
!    extinction coefficients for each state variable as specified in the model information,
!    specifically, the member "specific_light_extinction" of the state variable information.
!    (note: specific extinction coefficients default to 0: no attenuation due to biogeochemical components) 
!
! 7) If (part of) the model state variable are composed of conserved quantities (energy and/or
!    chemical elements), a function that provides the sum of these quantities given the model
!    state must be added as option to the "select" statement in get_conserved_quantities_rmbm.
!    Note that your model should reigster these conserved quantities during initialization.
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
! circulation models). Ultimately, the variable will be represented by an array slice, which
! means that the variable for storing is allowed to use additional dimensions (e.g., for time),
! which will be ignored (by being set to some index chosen by the host) by RMBM.
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
! Later thwe last dimension can be ignored by providing RMBM with a slices such as
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
! get at bottom values with a simple aray slice (with fixed depth index). Before calls to
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
!    This adds child models with identifiers 1 (currently NPZD) and 2 (currently jellyfish)
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
!    - rmbm_do:                   to get local temporal derivatives
!    - rmbm_get_vertical_movement:    to get current vertcial movement rates for the state variables
!    - rmbm_get_bio_extinction:       to get the combined light extinction coefficient due to
!                                     biogeochemical components
!    - rmbm_get_conserved_quantities: to get the sums of the conserved quantities described
!                                     by the model
!    - rmbm_update_air_sea_exchange:  to get updated fluxes over the air-sea interface
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
   use rmbm_jellyfish
   use rmbm_co2sys
   ! ADD_NEW_MODEL_HERE
!   
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_model, rmbm_init, rmbm_create_model, rmbm_do, &
          rmbm_link_variable_data,rmbm_get_variable_id, &
          rmbm_check_state, rmbm_get_vertical_movement, rmbm_get_bio_extinction, &
          rmbm_get_conserved_quantities, rmbm_update_air_sea_exchange
!
! !PRIVATE DATA MEMBERS:

!  Identifiers for specific biogeochemical models:
   integer, parameter :: model_container_id = -1
   integer, parameter :: npzd_id            =  1
   integer, parameter :: jellyfish_id       =  2
   integer, parameter :: carbonate_id       =  3
   ! ADD_NEW_MODEL_HERE

! !PUBLIC TYPES:
!
!  Single generic biogeochemical model
   type type_model
   
      ! Model identifier, name (= variable name prefix), and metadata.
      integer                :: id
      character(len=64)      :: name
      type (type_model_info) :: info
      
      ! Pointers for linking to parent and child models
      ! (models can be linked to form a tree structure)
      type (type_model),pointer :: parent      => NULL()
      type (type_model),pointer :: firstchild  => NULL()
      type (type_model),pointer :: nextsibling => NULL()

      ! Pointers to the current state and environment.
      ! These are allocated upon initialization by RMBM,
      ! and should be set by the physical host model.
      type (type_environment),pointer :: environment
      type (type_state),      pointer :: state(:)

      ! Derived types that belong to specific biogeochemical models.
      type (type_npzd)      :: npzd
      type (type_jellyfish) :: jellyfish
      type (type_co2_sys)   :: carbonate
      ! ADD_NEW_MODEL_HERE
      
   end type type_model
!
! !PUBLIC INTERFACES:
!
!  Access to temporal derivatives
   interface rmbm_do
      module procedure rmbm_do_rhs
      module procedure rmbm_do_ppdd
      module procedure rmbm_do_rhs_1d
   end interface rmbm_do
   
   interface rmbm_link_variable_data
      module procedure rmbm_supply_variable_data_2d
      module procedure rmbm_supply_variable_data_3d
   end interface rmbm_link_variable_data
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
! !IROUTINE:Get short model name (letters, numbers and underscores only)
!
! !INTERFACE:
   function get_model_name(bio_model) result(name)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in) :: bio_model
   character(len=64)   :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   character(len=256) :: text
!-----------------------------------------------------------------------
!BOC
   select case (bio_model)
      case (model_container_id)
         name = ''
      case (npzd_id)
         name = 'npzd'
      case (jellyfish_id)
         name = 'jellyfish'
      case (carbonate_id)
         name = 'co2sys'
      ! ADD_NEW_MODEL_HERE
      case default
         write (text,fmt='(i,a)') bio_model,' is not a valid bio model identifier.'
         call fatal_error('rmbm::get_model_name',text)
   end select
   
   end function get_model_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Create a new model. This will be a specific model if the
! model identifier is provided, or a contained of child models if the
! identifier is omitted.
!
! !INTERFACE:
   function rmbm_create_model(modelid,name,parent) result(model)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,         optional,         intent(in)    :: modelid
   character(len=*),optional,         intent(in)    :: name
   type (type_model),pointer,optional,intent(inout) :: parent
   type (type_model),pointer                        :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
   type (type_model),pointer                        :: previoussibling
!-----------------------------------------------------------------------
!BOC
   allocate(model)

   ! Set the model identifier.
   if (present(modelid)) then
      model%id = modelid
   else
      model%id = model_container_id
   end if

   ! Set the model name
   if (present(name)) then
      model%name = name
   else
      model%name = get_model_name(model%id)
   end if

   ! Initialize model info
   call init_model_info(model%info)

   ! Connect to parent model if specified.
   if (present(parent)) then
      if (parent%id.ne.model_container_id) &
         call fatal_error('rmbm_create_model','A child model can only be added to a container, not to an existing model.')

      if (associated(parent%firstchild)) then
         previoussibling => parent%firstchild
         do while (associated(previoussibling%nextsibling))
            previoussibling => previoussibling%nextsibling
         end do
         previoussibling%nextsibling => model
         previoussibling%info%nextsibling => model%info
      else
         parent%firstchild => model
         parent%info%firstchild => model%info
      end if
      model%parent => parent
      model%info%parent => parent%info
   end if
   
   end function rmbm_create_model
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the selected 0d biogeochemical model
!
! !INTERFACE:
   recursive subroutine rmbm_init(model,nmlunit)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),target,               intent(inout) :: model
   integer,                                intent(in)    :: nmlunit
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer                    :: i,ichild
  character(len= 64)         :: modelname
  type (type_model), pointer :: curchild,curchild2
  character(len=256)         :: childname
!EOP
!-----------------------------------------------------------------------
!BOC
   modelname = get_model_name(model%id)
   call log_message('Initializing biogeochemical model '//trim(modelname))

   ! Allow the selected model to initialize
   select case (model%id)
      case (npzd_id)
         call init_bio_npzd_0d(model%npzd,model%info,nmlunit)
      case (jellyfish_id)
         call init_bio_jellyfish_0d(model%jellyfish,model%info,nmlunit)
      case (carbonate_id)
         call init_bio_co2_sys_0d(model%carbonate,model%info,nmlunit)
      ! ADD_NEW_MODEL_HERE
      
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            ! Find the number of the child model
            ichild = 1
            curchild2 => model%firstchild
            do while (.not. associated(curchild,curchild2))
               if (curchild2%name.eq.curchild%name) ichild = ichild + 1
               curchild2 => curchild2%nextsibling
            end do

            ! Create prefixes for names of variables of the child model.
            write (unit=curchild%info%nameprefix,     fmt='(a,i2.2,a)') trim(curchild%name),ichild,'_'
            write (unit=curchild%info%longnameprefix, fmt='(a,i2.2,a)') trim(curchild%name),ichild,' '

            ! Initialize child model.
            call rmbm_init(curchild,nmlunit)
            
            ! Move to next child model.
            curchild => curchild%nextsibling
         end do
         
      case default
         call fatal_error('rmbm_init','no valid biogeochemical model specified!')
         
   end select
   
   if (.not. associated(model%parent)) then
      allocate(model%environment)
      allocate(model%environment%var2d(ubound(model%info%dependencies2d,1)))
      allocate(model%environment%var3d(ubound(model%info%dependencies3d,1)))
      allocate(model%state(model%info%state_variable_count))
      call set_state_data(model,model%state,model%environment)
   end if
   
   call log_message('model '//trim(modelname)//' initialized successfully.')

   end subroutine rmbm_init
!EOC

recursive subroutine set_state_data(model,state,environment)
   type (type_model),             intent(inout) :: model
   type (type_state),      target,intent(in)    :: state(:)
   type (type_environment),target,intent(in)    :: environment
   
   type (type_model), pointer :: curchild
   
   model%state       => state
   model%environment => environment
   curchild => model%firstchild
   do while (associated(curchild))
      call set_state_data(curchild,state,environment)
      curchild => curchild%nextsibling
   end do
end subroutine set_state_data

function rmbm_get_variable_id(model,name,shape) result(id)
   type (type_model),             intent(in)  :: model
   character(len=*),              intent(in)  :: name
   integer,                       intent(in)  :: shape
   integer                                    :: id
   character(len=64),pointer                  :: source(:)
   
   select case (shape)
      case (shape2d)
         source => model%info%dependencies2d
      case (shape3d)
         source => model%info%dependencies3d
   end select

   if (associated(source)) then
      do id=1,ubound(source,1)
         if (source(id)==name) return
      end do
   end if
   
   id = -1
end function rmbm_get_variable_id

recursive subroutine rmbm_supply_variable_data_3d(model,id,dat)
   type (type_model),                       intent(inout) :: model
   integer,                                 intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS,target,intent(in)    :: dat
   
   model%environment%var3d(id)%data => dat
end subroutine rmbm_supply_variable_data_3d

recursive subroutine rmbm_supply_variable_data_2d(model,id,dat)
   type (type_model),                          intent(inout) :: model
   integer,                                    intent(in)    :: id
   REALTYPE ATTR_LOCATION_DIMENSIONS_HZ,target,intent(in)    :: dat
   
   model%environment%var2d(id)%data => dat
end subroutine rmbm_supply_variable_data_2d


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the selected 0d
! biogeochemical model
!
! !INTERFACE:
   recursive subroutine rmbm_do_rhs(model,LOCATION,dy,diag)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE, dimension(:), intent(inout) :: dy,diag
!
! !LOCAL PARAMETERS:
   REALTYPE,allocatable                  :: pp(:,:),dd(:,:)
   integer                               :: i,j
   type (type_model), pointer            :: curchild
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (npzd_id)
         call do_bio_npzd_0d(model%npzd,model%state,model%environment,LOCATION,dy,diag)
      case (jellyfish_id)
         call do_bio_jellyfish_0d(model%jellyfish,model%state,model%environment,LOCATION,dy,diag)
      case (carbonate_id)
         call do_bio_co2_sys_0d(model%carbonate,model%state,model%environment,LOCATION,dy,diag)
      ! ADD_NEW_MODEL_HERE
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            call rmbm_do_rhs(curchild,LOCATION,dy,diag)
            curchild => curchild%nextsibling
         end do
      case default
         allocate(pp(1:ubound(dy,1),1:ubound(dy,1)))
         allocate(dd(1:ubound(dy,1),1:ubound(dy,1)))
         pp = _ZERO_
         dd = _ZERO_
         call rmbm_do_ppdd(model,LOCATION,pp,dd,diag)
         do i=1,ubound(dy,1)
            do j=1,ubound(dy,1)
               dy(i) = dy(i) + pp(i,j)-dd(i,j)
            end do
         end do
         deallocate(pp)
         deallocate(dd)
   end select

   end subroutine rmbm_do_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives in 1D
!
! !INTERFACE:
   recursive subroutine rmbm_do_rhs_1d(model ARG_LOCATION_1DLOOP,istart,istop,dy,diag)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   integer,                intent(in)    :: istart,istop
   DEFINE_LOCATION_1DLOOP
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: dy(:,:)
   REALTYPE,               intent(inout) :: diag(:,:)
!
! !LOCAL PARAMETERS:
   type (type_model), pointer            :: curchild
   LOCATION_TYPE                         :: VARIABLE_1DLOOP
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            call rmbm_do_rhs_1d(curchild ARG_LOCATION_1DLOOP,istart,istop,dy,diag)
            curchild => curchild%nextsibling
         end do
      case default
         do VARIABLE_1DLOOP=istart,istop
            call rmbm_do_rhs(curchild,LOCATION,dy(VARIABLE_1DLOOP,:),diag(VARIABLE_1DLOOP,:))
         end do
   end select

   end subroutine rmbm_do_rhs_1d
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the selected 0d
! biogeochemical model
!
! !INTERFACE:
   recursive subroutine rmbm_do_ppdd(model,LOCATION,pp,dd,diag)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,dimension(:,:),intent(inout) :: pp,dd
   REALTYPE,dimension(:),  intent(inout) :: diag
!
! !LOCAL PARAMETERS:
   REALTYPE,allocatable                  :: dy(:)
   integer                               :: i
   type (type_model), pointer            :: curchild
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (npzd_id)
         call do_bio_npzd_0d_ppdd(model%npzd,model%state,model%environment,LOCATION,pp,dd,diag)
      ! ADD_NEW_MODEL_HERE - optional
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            call rmbm_do_ppdd(curchild,LOCATION,pp,dd,diag)
            curchild => curchild%nextsibling
         end do
      case default
         allocate(dy(ubound(pp,1)))
         dy = _ZERO_
         call rmbm_do_rhs(model,LOCATION,dy,diag)
         do i=1,ubound(pp,1)
            pp(i,i) = pp(i,i) + dy(i)
         end do
         deallocate(dy)
   end select

   end subroutine rmbm_do_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Checks whether the current state is valid, and repairs [clips]
! invalid state variables if requested and possible.
!
! !INTERFACE:
   recursive function rmbm_check_state(model,LOCATION,repair) result(valid)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATION_TYPE,          intent(in)    :: LOCATION
   logical,                intent(in)    :: repair

   logical                               :: valid
!
! !LOCAL PARAMETERS:
   integer                               :: i
   type (type_model), pointer            :: curchild
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   valid = .true.
   select case (model%id)
      ! ADD_NEW_MODEL_HERE - optional
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            if (.not. rmbm_check_state(curchild,LOCATION,repair)) valid = .false.
            curchild => curchild%nextsibling
         end do
      case default
         do i=1,model%info%state_variable_count
            if (model%state(model%info%variables(i)%globalid)%data(LOCATION)<model%info%variables(i)%minimum) then
               valid = .false.
               if (repair) model%state(model%info%variables(i)%globalid)%data(LOCATION) = model%info%variables(i)%minimum
            elseif (model%state(model%info%variables(i)%globalid)%data(LOCATION)>model%info%variables(i)%maximum) then
               valid = .false.
               if (repair) model%state(model%info%variables(i)%globalid)%data(LOCATION) = model%info%variables(i)%maximum
            end if
         end do
   end select

   end function rmbm_check_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-sea exchange fluxes for all bio state variables.
! Positive values indicate fluxes into the ocean, negative values indicate fluxes
! out of the ocean. Units are m/s * tracer concentration.
!
! !INTERFACE:
   recursive subroutine rmbm_update_air_sea_exchange(model,LOCATION,flux)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: flux(:)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer            :: curchild
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (carbonate_id)
         call update_air_sea_co2_sys_0d(model%carbonate,model%state,model%environment,LOCATION,flux)
      ! ADD_NEW_MODEL_HERE - optional
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            call rmbm_update_air_sea_exchange(curchild,LOCATION,flux)
            curchild => curchild%nextsibling
         end do
      case default
   end select

   end subroutine rmbm_update_air_sea_exchange
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical movement rates (m/s) for the bio state variables.
! Note that negative values indidcate movement towards the bottom, e.g., sinking,
! and positive values indicate movemment towards the surface, e.g., floating.
!
! !INTERFACE:
   recursive subroutine rmbm_get_vertical_movement(model,LOCATION,vertical_movement)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATION_TYPE,          intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: vertical_movement(:)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer            :: curchild
   integer                               :: i
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      ! ADD_NEW_MODEL_HERE - optional
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            call rmbm_get_vertical_movement(curchild,LOCATION,vertical_movement)
            curchild => curchild%nextsibling
         end do
      case default
         ! Default: use the constant sinking rates specified in state variable properties.
         do i=1,model%info%state_variable_count
            vertical_movement(model%info%variables(i)%globalid) = model%info%variables(i)%vertical_movement
         end do
   end select

   end subroutine rmbm_get_vertical_movement
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   recursive function rmbm_get_bio_extinction(model,LOCATION) result(extinction)
!
! !INPUT PARAMETERS:
   type (type_model), intent(in) :: model
   LOCATION_TYPE,     intent(in) :: LOCATION
   REALTYPE                      :: extinction
   
   integer                       :: i
   type (type_model), pointer    :: curchild
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (npzd_id)
         extinction = get_bio_extinction_npzd_0d(model%npzd,model%state,model%environment,LOCATION)
      ! ADD_NEW_MODEL_HERE - optional
      case (model_container_id)
         extinction = _ZERO_
         curchild => model%firstchild
         do while (associated(curchild))
            extinction = extinction + rmbm_get_bio_extinction(curchild,LOCATION)
            curchild => curchild%nextsibling
         end do
      case default
         ! Default: use constant specific light extinction values specified in the state variable properties
         extinction = _ZERO_
         do i=1,model%info%state_variable_count
            if (model%info%variables(i)%specific_light_extinction.ne._ZERO_) &
               extinction = extinction + model%state(model%info%variables(i)%globalid)%data(LOCATION) &
                                 & * model%info%variables(i)%specific_light_extinction
         end do
   end select

   end function rmbm_get_bio_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities
!
! !INTERFACE:
   recursive subroutine rmbm_get_conserved_quantities(model,LOCATION,sums)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model), intent(in)      :: model
   LOCATION_TYPE,     intent(in)      :: LOCATION
   REALTYPE,          intent(inout)   :: sums(:)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer    :: curchild

!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (npzd_id)
         call get_conserved_quantities_npzd_0d(model%npzd,model%state,model%environment,LOCATION,sums)
      ! ADD_NEW_MODEL_HERE - optional
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            call rmbm_get_conserved_quantities(curchild,LOCATION,sums)
            curchild => curchild%nextsibling
         end do
      case default
         ! Default: the model does not describe any conserved quantities.
         if (model%info%conserved_quantity_count.gt.0) &
            call fatal_error('rmbm_get_conserved_quantities','the model specifies that it describes one or more conserved &
                 &quantities, but a function that provides sums of these quantities has not been specified.')
   end select

   end subroutine rmbm_get_conserved_quantities
!EOC

end module rmbm

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

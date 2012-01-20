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
! ADD_NEW_MODEL_HERE strings in the code comments.
!
! For more information, see the documentation at /doc/documentation.pdf.
!
! !USES:
   use fabm_types
   use fabm_driver
   use fabm_library
!
!  Reference modules of specific biogeochemical models
   use fabm_gotm_npzd
   use fabm_gotm_fasham
   use fabm_metu_mnemiopsis
   use fabm_pml_ersem
   use fabm_pml_carbonate
   use fabm_examples_benthic_predator
   use fabm_examples_npzd_det
   use fabm_examples_npzd_nut
   use fabm_examples_npzd_phy
   use fabm_examples_npzd_zoo
   use fabm_iow_ergom
   use fabm_bb_passive
   ! ADD_NEW_MODEL_HERE - required if the model is contained in a Fortran 90 module

   implicit none
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_model, fabm_create_model, fabm_create_model_from_file, fabm_init, fabm_set_domain, fabm_check_ready, &
          fabm_get_variable_id,fabm_link_benthos_state_data, fabm_link_state_data, &
          fabm_link_data,fabm_link_data_hz,fabm_link_scalar, &
          fabm_get_diagnostic_data, fabm_get_diagnostic_data_hz, &
          fabm_do, fabm_check_state, fabm_get_vertical_movement, fabm_get_light_extinction, &
          fabm_get_conserved_quantities, fabm_get_surface_exchange, fabm_do_benthos

#ifndef _FABM_MANAGE_DIAGNOSTICS_
   public fabm_link_diagnostic_data_hz,fabm_link_diagnostic_data
#endif

! !PUBLIC TYPES:
!
   ! Derived type for a single generic biogeochemical model
   type type_model
      ! Model identifier and metadata.
      integer                :: id
      _CLASS_ (type_model_info),pointer :: info

      ! Pointers for linking to parent and child models
      ! (models can be linked to form a tree structure)
      type (type_model),pointer :: parent
      type (type_model),pointer :: firstchild
      type (type_model),pointer :: nextsibling

      ! Pointer to next non-container model.
      ! To speed up iterations over all models in the tree, non-container models are
      ! connected in a singly linked list.
      ! For the root model, nextmodel points to the first non-container model in the list.
      ! For nested (non-root) container models, the pointer remains disassociated(!)
      type (type_model),pointer :: nextmodel

      ! Derived types that belong to specific biogeochemical models.
      type (type_gotm_npzd)                 :: gotm_npzd
      type (type_gotm_fasham)               :: gotm_fasham
      type (type_metu_mnemiopsis)           :: metu_mnemiopsis
      type (type_pml_ersem)                 :: pml_ersem
      type (type_pml_carbonate)             :: pml_carbonate
      type (type_examples_benthic_predator) :: examples_benthic_predator
      type (type_examples_npzd_det)         :: examples_npzd_det
      type (type_examples_npzd_nut)         :: examples_npzd_nut
      type (type_examples_npzd_phy)         :: examples_npzd_phy
      type (type_examples_npzd_zoo)         :: examples_npzd_zoo
      type (type_iow_ergom)                 :: iow_ergom
      type (type_bb_passive)                :: bb_passive
      ! ADD_NEW_MODEL_HERE - required if the model groups its data in a custom derived type

      ! Pointer to the current spatially explicit environment.
      ! It is assigned to during initialization by FABM.
      type (type_environment),pointer :: environment

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
   interface fabm_do_benthos
      module procedure fabm_do_benthos_rhs
      module procedure fabm_do_benthos_ppdd
   end interface

   ! Subroutine for providing FABM with variable data on the full spatial domain.
   interface fabm_link_data
      module procedure fabm_link_data
      module procedure fabm_link_data_char
   end interface

   ! Subroutine for providing FABM with variable data on horizontal slices of the domain.
   interface fabm_link_data_hz
      module procedure fabm_link_data_hz
      module procedure fabm_link_data_hz_char
   end interface

   ! Subroutine for providing FABM with scalar variable data.
   interface fabm_link_scalar
      module procedure fabm_link_scalar
      module procedure fabm_link_scalar_char
   end interface

   ! Function for creating new models based on integer id [deprecated] or name.
   interface fabm_create_model
      module procedure fabm_create_model_by_id
      module procedure fabm_create_model_by_name
   end interface
!
! !PRIVATE DATA MEMBERS:

!  Identifiers for specific biogeochemical models.
   integer, parameter :: model_container_id           = -1
#ifdef _FABM_F2003_
   integer, parameter :: model_f2003_id               =  0
#endif
   integer, parameter :: gotm_npzd_id                 =  1
   integer, parameter :: gotm_fasham_id               =  4
   integer, parameter :: pml_ersem_id                 =  99
   integer, parameter :: pml_carbonate_id             =  101
   integer, parameter :: metu_mnemiopsis_id           =  102
   integer, parameter :: examples_benthic_predator_id =  103
   integer, parameter :: examples_npzd_det_id         =  104
   integer, parameter :: examples_npzd_nut_id         =  105
   integer, parameter :: examples_npzd_phy_id         =  106
   integer, parameter :: examples_npzd_zoo_id         =  107
   integer, parameter :: iow_ergom_id                 =  108
   integer, parameter :: bb_passive_id                =  109
   ! ADD_NEW_MODEL_HERE - required. Identifier values are arbitrary, but they must be unique.
   ! Note: values <=100 are reserved for models ported from the General Ocean Turbulence Model.

   ! Arrays for integer identifiers and names of all available biogeochemical models.
   ! Allocated and set on demand when model names are first used.
   integer,          allocatable,dimension(:) :: modelids
   character(len=64),allocatable,dimension(:) :: modelnames
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
! !IROUTINE:Register the integer identifiers and short names (letters,
! numbers and underscores only) of all biogeochemical models.
!
! !INTERFACE:
   subroutine register_models()
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
!-----------------------------------------------------------------------
!BOC
   ! Start with empty arrays with model identifiers and names.
   allocate(modelids  (0))
   allocate(modelnames(0))

   ! Register specific biogeochemical models.
   call register_model(model_container_id,          '')
   call register_model(gotm_npzd_id,                'gotm_npzd')
   call register_model(gotm_fasham_id,              'gotm_fasham')
   call register_model(metu_mnemiopsis_id,          'metu_mnemiopsis')
   call register_model(pml_ersem_id,                'pml_ersem')
   call register_model(pml_carbonate_id,            'pml_carbonate')
   call register_model(examples_benthic_predator_id,'examples_benthic_predator')
   call register_model(examples_npzd_det_id,        'examples_npzd_det')
   call register_model(examples_npzd_nut_id,        'examples_npzd_nut')
   call register_model(examples_npzd_phy_id,        'examples_npzd_phy')
   call register_model(examples_npzd_zoo_id,        'examples_npzd_zoo')
   call register_model(iow_ergom_id,                'iow_ergom')
   call register_model(bb_passive_id,               'bb_passive')
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
! !INPUT PARAMETERS:
   integer, intent(in)         :: id
   character(len=*),intent(in) :: name
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
   integer :: oldcount,i
   integer,          allocatable,dimension(:) :: modelids_old
   character(len=64),allocatable,dimension(:) :: modelnames_old
   character(len=256) :: text
!
!-----------------------------------------------------------------------
!BOC
   ! Determine original model count.
   oldcount = ubound(modelids,1)

   ! First check whether the provided model identifier or name are not in use yet.
   do i=1,ubound(modelids,1)
      if (modelids(i)==id) then
         write (text,fmt='(a,i4,a)') 'model identifier ',id,' has already been registered.'
         call fatal_error('fabm::register_model_name',text)
      end if
      if (modelnames(i)==name) &
         call fatal_error('fabm::register_model_name','model name "'//trim(name)//'" has already been registered.')
   end do

   ! Copy current identifiers and names to temporary storage.
   allocate(modelids_old  (oldcount))
   allocate(modelnames_old(oldcount))
   modelids_old  (:) = modelids
   modelnames_old(:) = modelnames

   ! Create extended arrays
   deallocate(modelids)
   deallocate(modelnames)
   allocate(modelids  (oldcount+1))
   allocate(modelnames(oldcount+1))

   ! Copy the old identifiers and names to the new extended array.
   modelids  (1:oldcount) = modelids_old  (:)
   modelnames(1:oldcount) = modelnames_old(:)

   ! Deallocate the temporary arrays with old values.
   deallocate(modelids_old)
   deallocate(modelnames_old)

   ! Add the new identifier and name.
   modelids  (oldcount+1) = id
   modelnames(oldcount+1) = name

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
   ! If models have not been registered yet, do this first.
   if (.not.allocated(modelids)) call register_models()

   ! Enumerate all models and compare their integer id with the supplied one.
   ! Return the name if a match is found.
   do i=1,ubound(modelids,1)
      if (modelids(i)==id) then
         name = modelnames(i)
         return
      end if
   end do

   ! Model identifier was not found - throw an error.
   write (text,fmt='(i4,a)') id,' is not a valid model identifier registered in fabm::register_models.'
   call fatal_error('fabm::get_model_name',text)

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
   ! Default id: model not found.
   id = id_not_used

   ! If models have not been registered yet, do this first.
   if (.not.allocated(modelnames)) call register_models()

   ! Enumerate all models and compare their name with the supplied one.
   ! Return the integer id if a match is found.
   do i=1,ubound(modelnames,1)
      if (modelnames(i)==name) then
         id = modelids(i)
         return
      end if
   end do

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
   function fabm_create_model_by_id(modelid,parent,info) result(model)
!
! !INPUT PARAMETERS:
   integer,                 optional,intent(in)    :: modelid
   type (type_model),target,optional,intent(inout) :: parent
   _CLASS_ (type_model_info),pointer,optional         :: info
   type (type_model),pointer                       :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
   type (type_model),pointer                       :: curmodel
   integer                                         :: modelid_eff
   _CLASS_ (type_model_info),pointer               :: info_eff
!-----------------------------------------------------------------------
!BOC
   ! Determine effective model identifier.
   if (present(modelid)) then
      modelid_eff = modelid
   else
      modelid_eff = model_container_id
   end if

   ! Determine effective model info (null if not present).
   if (present(info)) then
      info_eff => info
   else
      nullify(info_eff)
   end if

   ! Allocate storage space for the model.
   allocate(model)

   if (associated(info_eff)) then
      ! Use user-supplied model info.
      model%info => info_eff
   else
      ! Allocate and initialize model info.
      allocate(model%info)
      call init_model_info(model%info)
      model%info%name = get_model_name(modelid_eff)
   end if

   ! Make sure the pointers to parent and sibling models are dissociated.
   nullify(model%parent)
   nullify(model%firstchild)
   nullify(model%nextsibling)
   nullify(model%nextmodel)

   ! Make sure the pointer to the current environment is dissociated.
   nullify(model%environment)

   ! Set the model identifier.
   model%id = modelid_eff

   ! Connect to parent container if provided.
   if (present(parent)) then
      ! Make sure the provided parent model is a container.
      if (parent%id.ne.model_container_id) &
         call fatal_error('fabm_create_model','A child model can only be added to a container, not to an existing model.')

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
   else
      ! No parent provided - ensure that the created model is a container.
      if (model%id.ne.model_container_id) &
         call fatal_error('fabm_create_model_by_id','Non-container models must be created as children of an existing container.')
   end if

   if (model%id.ne.model_container_id) then
      ! This model is not a container.
      ! Add the model to the flattened list of non-container models,
      ! which is used at runtime to iterate over all models without needing recursion.

      ! First non-container model will be pointed to by the root of the tree - find it.
      ! Note that the code above ensures that the model has a parent.
      curmodel => model%parent
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

   end function fabm_create_model_by_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Create a new model by name. This will be a specific model if the
! model name is provided, or a container of child models if the
! identifier is omitted.
!
! !INTERFACE:
   function fabm_create_model_by_name(modelname,parent) result(model)
!
! !INPUT PARAMETERS:
   character(len=*),                 intent(in)    :: modelname
   type (type_model),target,optional,intent(inout) :: parent
   type (type_model),pointer                       :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
   integer :: modelid
   _CLASS_ (type_model_info),pointer :: modelinfo
!-----------------------------------------------------------------------
!BOC
   nullify(modelinfo)

   ! Resolve model name.
   modelid = get_model_id(modelname)

#ifdef _FABM_F2003_
   ! Try to get model from Fortran 2003 model library
   if (modelid.eq.id_not_used) then
      modelinfo => fabm_library_create_model(modelname)
      if (associated(modelinfo)) then
         modelinfo%name = modelname
         modelid = model_f2003_id
      end if
   end if
#endif

   ! If we have not found the model now, throw an error.
   if (modelid.eq.id_not_used) &
      call fatal_error('fabm_create_model_by_name','"'//trim(modelname)//'" is not a valid model name registered in fabm::register_models.')

   ! Obtain the integer model identifier, and redirect to the function operating on that.
   if (present(parent)) then
      model => fabm_create_model_by_id(modelid,parent=parent,info=modelinfo)
   else
      model => fabm_create_model_by_id(modelid,info=modelinfo)
   end if

   end function fabm_create_model_by_name
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Create a new model tree from a configuration file.
!
! !INTERFACE:
   function fabm_create_model_from_file(file_unit,file) result(model)
!
! !INPUT PARAMETERS:
   character(len=*),optional,intent(in) :: file
   integer,                  intent(in) :: file_unit
   type (type_model),pointer            :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman

   logical                   :: isopen
   character(len=256)        :: file_eff
   integer                   :: i
   character(len=64)         :: models(256)
   type (type_model),pointer :: childmodel
   namelist /fabm_nml/ models
!EOP
!-----------------------------------------------------------------------
!BOC
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
   model => fabm_create_model()
   do i=1,ubound(models,1)
      if (models(i).ne.'') childmodel => fabm_create_model(trim(models(i)),parent=model)
   end do

   ! Initialize model tree
   call fabm_init(model,file_unit)

   if (.not.isopen) then
      ! We have opened the configuration file ourselves - close it.
      close(file_unit)
   end if

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
   subroutine fabm_init(root,nmlunit)
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
   if (associated(root%parent)) call fatal_error('fabm_init','fabm_init can only be called on the root of a model tree.')

   ! Check whether the unit provided by the host actually refers to an open file.
   inquire(nmlunit,opened=isopen)
   if (.not.isopen) call fatal_error('fabm_init','input configuration file has not been opened yet.')

   ! Initialize the model (this automatically initializes all contained models)
   call init_model(root,nmlunit)

   ! Allocate arrays for storage of (references to) data.
   allocate(root%environment)

   ! Set all pointers to external data to dissociated.
   allocate(root%environment%var       (ubound(root%info%dependencies,       1)))
   allocate(root%environment%var_hz    (ubound(root%info%dependencies_hz,    1)))
   allocate(root%environment%var_scalar(ubound(root%info%dependencies_scalar,1)))
   do ivar=1,ubound(root%environment%var,1)
      nullify(root%environment%var(ivar)%data)
   end do
   do ivar=1,ubound(root%environment%var_hz,1)
      nullify(root%environment%var_hz(ivar)%data)
   end do
   do ivar=1,ubound(root%environment%var_scalar,1)
      nullify(root%environment%var_scalar(ivar)%data)
   end do

   ! Transfer pointer to environment to all child models.
   call set_model_data_members(root,root%environment)

   end subroutine fabm_init
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with the extents of the spatial domain.
!
! !INTERFACE:
   subroutine fabm_set_domain(root _ARG_LOCATION_)
!
! !INPUT PARAMETERS:
   type (type_model),target,               intent(inout) :: root
   _DECLARE_LOCATION_ARG_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  type (type_model), pointer :: model
#ifdef _FABM_MANAGE_DIAGNOSTICS_
  integer                    :: ivar
#endif
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Check whether we are operating on the root of a model tree.
   if (associated(root%parent)) &
      call fatal_error('fabm_set_domain','fabm_set_domain can only be called on the root of a model tree.')

   ! Enumerate all non-container models in the tree.
   model => root%nextmodel
   do while (associated(model))
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%set_domain(_LOCATION_)
#endif
         case (pml_ersem_id)
            call pml_ersem_set_domain(model%pml_ersem,_DOMAIN_SIZE_1D_)
         ! ADD_NEW_MODEL_HERE - optional, only needed if the model needs to be informed about the spatial domain.
         !
         ! Typical model call:
         ! call MODELNAME_set_domain(model%MODELNAME,_DOMAIN_SIZE_1D_)
      end select
      model => model%nextmodel
   end do

#ifdef _FABM_MANAGE_DIAGNOSTICS_
   ! FABM will manage and store current values of diagnostic variables.
   ! Allocate memory for this, and link this memory to FABM's variable data pointers.

   ! Allocate arrays for diagnostic variables defined on the full domain and on horizontonal slices.
   allocate(root%environment%diag   (ubound(root%info%diagnostic_variables,1) _ARG_LOCATION_))
   allocate(root%environment%diag_hz(ubound(root%info%diagnostic_variables_hz,1) _ARG_LOCATION_HZ_))

   ! Initialize diagnostic variables to zero.
   root%environment%diag    = _ZERO_
   root%environment%diag_hz = _ZERO_

   ! If diagnostic variables also appear as dependency, send the corresponding array slice for generic read-only access.
   do ivar=1,ubound(root%info%diagnostic_variables,1)
      call fabm_link_data(root,root%info%diagnostic_variables(ivar)%globalid%dependencyid,root%environment%diag(ivar _ARG_LOCATION_DIMENSIONS_))
   end do
   do ivar=1,ubound(root%info%diagnostic_variables_hz,1)
      call fabm_link_data_hz(root,root%info%diagnostic_variables_hz(ivar)%globalid%dependencyid,root%environment%diag_hz(ivar _ARG_LOCATION_DIMENSIONS_HZ_))
   end do
#endif

   end subroutine fabm_set_domain
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check whether FABM has been provided with all required data.
!
! !INTERFACE:
   subroutine fabm_check_ready(root)
!
! !INPUT PARAMETERS:
   type (type_model),intent(in) :: root
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer                    :: ivar
  logical                    :: ready
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Check whether we are operating on the root of a model tree.
   if (associated(root%parent)) &
      call fatal_error('fabm_check_ready','fabm_check_ready can only be called on the root of a model tree.')

   ready = .true.
   do ivar=1,ubound(root%environment%var,1)
      if (.not.associated(root%environment%var(ivar)%data)) then
         call log_message('data for dependency "'//trim(root%info%dependencies(ivar))// &
            & '", defined on the full model domain, have not been provided.')
         ready = .false.
      end if
   end do
   do ivar=1,ubound(root%environment%var_hz,1)
      if (.not.associated(root%environment%var_hz(ivar)%data)) then
         call log_message('data for dependency "'//trim(root%info%dependencies_hz(ivar))// &
            &  '", defined on a horizontal slice of the model domain, have not been provided.')
         ready = .false.
      end if
   end do
   do ivar=1,ubound(root%environment%var_scalar,1)
      if (.not.associated(root%environment%var_scalar(ivar)%data)) then
         call log_message('data for dependency "'//trim(root%info%dependencies_scalar(ivar))// &
            &  '", defined as global scalar quantity, have not been provided.')
         ready = .false.
      end if
   end do
   if (.not.ready) call fatal_error('fabm_check_ready','FABM is lacking required data.')

   end subroutine fabm_check_ready
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the provided biogeochemical model
!
! !INTERFACE:
   recursive subroutine init_model(model,nmlunit)
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
  type (type_model), pointer :: curchild,curchild2
  logical                    :: isopen
  logical,parameter          :: alwayspostfixindex=.false.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Log start of initialization of this model, unless it is a container only.
   if (model%id.ne.model_container_id) &
      call log_message('Initializing biogeochemical model "'//trim(model%info%name)//'"...')

   ! Allow the selected model to initialize
   select case (model%id)
#ifdef _FABM_F2003_
      case (model_f2003_id)
         call model%info%initialize(nmlunit)
#endif
      case (gotm_npzd_id)
         call gotm_npzd_init(model%gotm_npzd,model%info,nmlunit)
      case (gotm_fasham_id)
         call gotm_fasham_init(model%gotm_fasham,model%info,nmlunit)
      case (pml_ersem_id)
         call pml_ersem_init(model%pml_ersem,model%info,nmlunit)
      case (metu_mnemiopsis_id)
         call metu_mnemiopsis_init(model%metu_mnemiopsis,model%info,nmlunit)
      case (pml_carbonate_id)
         call pml_carbonate_init(model%pml_carbonate,model%info,nmlunit)
      case (examples_benthic_predator_id)
         call examples_benthic_predator_init(model%examples_benthic_predator,model%info,nmlunit)
      case (examples_npzd_det_id)
         call examples_npzd_det_init(model%examples_npzd_det,model%info,nmlunit)
      case (examples_npzd_nut_id)
         call examples_npzd_nut_init(model%examples_npzd_nut,model%info,nmlunit)
      case (examples_npzd_phy_id)
         call examples_npzd_phy_init(model%examples_npzd_phy,model%info,nmlunit)
      case (examples_npzd_zoo_id)
         call examples_npzd_zoo_init(model%examples_npzd_zoo,model%info,nmlunit)
      case (iow_ergom_id)
         call iow_ergom_init(model%iow_ergom,model%info,nmlunit)
      case (bb_passive_id)
         call bb_passive_init(model%bb_passive,model%info,nmlunit)
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
         call fatal_error('init_model','model "'//trim(model%info%name)//'" has not been registered in init_model.')

   end select

   ! Log successful initialization of this model, unless it is a container only.
   if (model%id.ne.model_container_id) &
      call log_message('model "'//trim(model%info%name)//'" initialized successfully.')

   ! Debug check: make sure the unit provided by the host has not been closed by the biogeochemical model.
   inquire(nmlunit,opened=isopen)
   if (.not.isopen) call fatal_error('init_model','input configuration file was closed by model "'//trim(model%info%name)//'".')

   end subroutine init_model
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide all descendant models with a pointer to the environment.
!
! !INTERFACE:
   recursive subroutine set_model_data_members(model,environment)
!
! !INPUT PARAMETERS:
   type (type_model),             intent(inout) :: model
   type (type_environment),target,intent(in)    :: environment
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   type (type_model), pointer :: curchild
!EOP
!-----------------------------------------------------------------------
!BOC
   model%environment => environment
   call freeze_model_info(model%info)
   curchild => model%firstchild
   do while (associated(curchild))
      call set_model_data_members(curchild,environment)
      curchild => curchild%nextsibling
   end do

   end subroutine set_model_data_members
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Obtain the integer variable identifier for the given variable
! name. Returns id_not_used if the variable name is unknown.
! The variable identifier can be used later in calls to fabm_link_data/fabm_link_data_hz.
!
! !INTERFACE:
   function fabm_get_variable_id(model,name,shape) result(id)
!
! !INPUT PARAMETERS:
   type (type_model),             intent(in)  :: model
   character(len=*),              intent(in)  :: name
   integer,                       intent(in)  :: shape
!
! !RETURN VALUE:
   integer                                    :: id
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   character(len=64),pointer                  :: source(:)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain a reference to the appropriate array with variable names.
   select case (shape)
      case (shape_hz)
         source => model%info%dependencies_hz
      case (shape_full)
         source => model%info%dependencies
   end select

   ! Try to locate the variable in the array with variable names.
   ! Return the identifier when found.
   if (associated(source)) then
      do id=1,ubound(source,1)
         if (source(id)==name) return
      end do
   end if

   ! Return default identifier: variable not found.
   id = id_not_used

   end function fabm_get_variable_id
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on the full spatial domain. The variable
! is identified by its integer id.
!
! !INTERFACE:
   subroutine fabm_link_data(model,id,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                       intent(inout) :: model
   integer,                                 intent(in)    :: id
   REALTYPE _ATTR_LOCATION_DIMENSIONS_,target,intent(in)  :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store a pointer to the provided array.
   model%environment%var(id)%data => dat

   end subroutine fabm_link_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by its name.
!
! !INTERFACE:
   subroutine fabm_link_data_char(model,name,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                         intent(inout) :: model
   character(len=*),                          intent(in)    :: name
   REALTYPE _ATTR_LOCATION_DIMENSIONS_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                                :: id
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_variable_id(model,name,shape_full)

   ! Only link the data if needed (if the variable identifier is valid).
   if (id.ne.id_not_used) call fabm_link_data(model,id,dat)

   end subroutine fabm_link_data_char
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by its integer id.
!
! !INTERFACE:
   subroutine fabm_link_data_hz(model,id,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                            intent(inout) :: model
   integer,                                      intent(in)    :: id
   REALTYPE _ATTR_LOCATION_DIMENSIONS_HZ_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store a pointer to the provided array.
   model%environment%var_hz(id)%data => dat

   end subroutine fabm_link_data_hz
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! the specified variable, defined on a horizontal slice of the spatial domain.
! The variable is identified by its name.
!
! !INTERFACE:
   subroutine fabm_link_data_hz_char(model,name,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                          intent(inout) :: model
   character(len=*),                           intent(in)    :: name
   REALTYPE _ATTR_LOCATION_DIMENSIONS_HZ_,target,intent(in)  :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                                   :: id
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_variable_id(model,name,shape_hz)

   ! Only link the data if needed (if the variable identifier is valid).
   if (id.ne.id_not_used) call fabm_link_data_hz(model,id,dat)

   end subroutine fabm_link_data_hz_char
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the scalar that will hold
! data for the specified variable. The variable is identified by its integer id.
!
! !INTERFACE:
   subroutine fabm_link_scalar(model,id,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                            intent(inout) :: model
   integer,                                      intent(in)    :: id
   REALTYPE,target,                              intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store a pointer to the provided array.
   model%environment%var_scalar(id)%data => dat

   end subroutine fabm_link_scalar
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the scalar that will hold
! data for the specified variable. The variable is identified by its name.
!
! !INTERFACE:
   subroutine fabm_link_scalar_char(model,name,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                          intent(inout) :: model
   character(len=*),                           intent(in)    :: name
   REALTYPE,target,                            intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
   integer                                                   :: id
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Obtain integer identifier of the variable.
   id = fabm_get_variable_id(model,name,shape_scalar)

   ! Only link the data if needed (if the variable identifier is valid).
   if (id.ne.id_not_used) call fabm_link_scalar(model,id,dat)

   end subroutine fabm_link_scalar_char
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! a single pelagic state variable.
!
! !INTERFACE:
   subroutine fabm_link_state_data(model,id,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                         intent(inout) :: model
   integer,                                   intent(in)    :: id
   REALTYPE _ATTR_LOCATION_DIMENSIONS_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_link_data(model,model%info%state_variables(id)%globalid%dependencyid,dat)

   end subroutine fabm_link_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! a single benthic state variable.
!
! !INTERFACE:
   subroutine fabm_link_benthos_state_data(model,id,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                            intent(inout) :: model
   integer,                                      intent(in)    :: id
   REALTYPE _ATTR_LOCATION_DIMENSIONS_HZ_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   call fabm_link_data_hz(model,model%info%state_variables_ben(id)%globalid%dependencyid,dat)

   end subroutine fabm_link_benthos_state_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Returns (a pointer to) the array with data for
! a single diagnostic variable, defined on the full spatial domain.
!
! !INTERFACE:
   function fabm_get_diagnostic_data(model,id) result(dat)
!
! !INPUT PARAMETERS:
   type (type_model),                       intent(in) :: model
   integer,                                 intent(in) :: id
   REALTYPE _ATTR_LOCATION_DIMENSIONS_,pointer         :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Retrieve a pointer to the array holding the requested data.
   dat => model%environment%var(model%info%diagnostic_variables(id)%globalid%dependencyid)%data

   end function fabm_get_diagnostic_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Returns (a pointer to) the array with data for
! a single diagnostic variable, defined on a horitontal slice of the
! spatial domain.
!
! !INTERFACE:
   function fabm_get_diagnostic_data_hz(model,id) result(dat)
!
! !INPUT PARAMETERS:
   type (type_model),                          intent(in) :: model
   integer,                                    intent(in) :: id
   REALTYPE _ATTR_LOCATION_DIMENSIONS_HZ_,pointer         :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Retrieve a pointer to the array holding the requested data.
   dat => model%environment%var_hz(model%info%diagnostic_variables_hz(id)%globalid%dependencyid)%data

   end function fabm_get_diagnostic_data_hz
!EOC

#ifndef _FABM_MANAGE_DIAGNOSTICS_

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! a single diagnostic state variable, defined on the full spatial domain.
!
! !INTERFACE:
   subroutine fabm_link_diagnostic_data(model,id,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                         intent(inout) :: model
   integer,                                   intent(in)    :: id
   REALTYPE _ATTR_LOCATION_DIMENSIONS_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Diagnostic data is managed by the host, which means that FABM treats it like a generic
   ! variable (e.g., an external dependency). Redirect to the generic function.
   call fabm_link_data(model,model%info%diagnostic_variables(id)%globalid%dependencyid,dat)

   end subroutine fabm_link_diagnostic_data
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Provide FABM with (a pointer to) the array with data for
! a single diagnostic state variable, defined on a horizontal slice of the
! spatial domain.
!
! !INTERFACE:
   subroutine fabm_link_diagnostic_data_hz(model,id,dat)
!
! !INPUT PARAMETERS:
   type (type_model),                            intent(inout) :: model
   integer,                                      intent(in)    :: id
   REALTYPE _ATTR_LOCATION_DIMENSIONS_HZ_,target,intent(in)    :: dat
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Diagnostic data is managed by the host, which means that FABM treats it like a generic
   ! variable (e.g., an external dependency). Redirect to the generic function.
   call fabm_link_data_hz(model,model%info%diagnostic_variables_hz(id)%globalid%dependencyid,dat)

   end subroutine fabm_link_diagnostic_data_hz
!EOC

#endif

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the biogeochemical
! model tree.
!
! !INTERFACE:
   subroutine fabm_do_rhs(root _ARG_LOCATION_ND_,dy)
!
! !INPUT PARAMETERS:
   type (type_model),      intent(inout) :: root
   _DECLARE_LOCATION_ARG_ND_
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE _ATTR_DIMENSIONS_1_, intent(inout) :: dy
!
! !LOCAL PARAMETERS:
   type (type_model), pointer            :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_DO_RHS_ _FABM_ARGS_ND_IN_,dy

   ! Ensure that this subroutine is called on the root of the model tree only.
   if (associated(root%parent)) &
      call fatal_error('fabm_do_rhs','fabm_do_rhs may only be called on the root of the model tree, or on non-container models.')

   ! Enumerate all non-container models in the tree.
   model => root%nextmodel
   do while (associated(model))
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%do(_INPUT_ARGS_DO_RHS_)
#endif
         case (gotm_npzd_id)
            call gotm_npzd_do(model%gotm_npzd,_INPUT_ARGS_DO_RHS_)
         case (gotm_fasham_id)
            call gotm_fasham_do(model%gotm_fasham,_INPUT_ARGS_DO_RHS_)
         case (pml_ersem_id)
            call pml_ersem_do(model%pml_ersem,_INPUT_ARGS_DO_RHS_)
         case (metu_mnemiopsis_id)
            call metu_mnemiopsis_do(model%metu_mnemiopsis,_INPUT_ARGS_DO_RHS_)
         case (pml_carbonate_id)
            call pml_carbonate_do(model%pml_carbonate,_INPUT_ARGS_DO_RHS_)
         case (examples_benthic_predator_id)
         case (examples_npzd_det_id)
            call examples_npzd_det_do(model%examples_npzd_det,_INPUT_ARGS_DO_RHS_)
         case (examples_npzd_nut_id)
            call examples_npzd_nut_do(model%examples_npzd_nut,_INPUT_ARGS_DO_RHS_)
         case (examples_npzd_phy_id)
            call examples_npzd_phy_do(model%examples_npzd_phy,_INPUT_ARGS_DO_RHS_)
         case (examples_npzd_zoo_id)
            call examples_npzd_zoo_do(model%examples_npzd_zoo,_INPUT_ARGS_DO_RHS_)
         case(iow_ergom_id)
            call iow_ergom_do(model%iow_ergom,_INPUT_ARGS_DO_RHS_)
         ! ADD_NEW_MODEL_HERE - required if the model features one or more active pelagic state variables,
         ! unless the model provides production/destruction matrices instead of a temporal derivative vector.
         ! In that case, add the model to fabm_do_ppdd.
         !
         ! Typical model call:
         ! call MODELNAME_do(model%MODELNAME,_INPUT_ARGS_DO_RHS_)
      end select
      model => model%nextmodel
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
   recursive subroutine fabm_do_ppdd(root _ARG_LOCATION_ND_,pp,dd)
!
! !INPUT PARAMETERS:
   type (type_model),           intent(inout) :: root
  _DECLARE_LOCATION_ARG_ND_
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE _ATTR_DIMENSIONS_2_,intent(inout) :: pp,dd
!
! !LOCAL PARAMETERS:
   type (type_model), pointer                 :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_DO_PPDD_ _FABM_ARGS_ND_IN_,pp,dd

   ! Enumerate all non-container models in the tree.
   model => root%nextmodel
   do while (associated(model))
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%do_ppdd(_INPUT_ARGS_DO_PPDD_)
#endif
         case (gotm_npzd_id)
            call gotm_npzd_do_ppdd(model%gotm_npzd,_INPUT_ARGS_DO_PPDD_)
         case (gotm_fasham_id)
            call gotm_fasham_do_ppdd(model%gotm_fasham,_INPUT_ARGS_DO_PPDD_)
         ! ADD_NEW_MODEL_HERE - optional, only if the model provides a subroutine for calculating local
         ! production/destruction matrices. This is required for certain temporal integration schemes,
         ! e.g., Patankar, Modified Patankar.
         !
         ! Typical model call:
         ! call MODELNAME_do_ppdd(model%MODELNAME,_INPUT_ARGS_DO_PPDD_)

         case default
           call fatal_error('fabm_do_ppdd','model "'//trim(model%info%name)//'" does not provide a subroutine &
              &that calculates local production/destruction matrices.')
      end select
      model => model%nextmodel
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
   subroutine fabm_check_state(root _ARG_LOCATION_ND_,repair,valid)
!
! !INPUT PARAMETERS:
   type (type_model),      intent(inout) :: root
   _DECLARE_LOCATION_ARG_ND_
   logical,                intent(in)    :: repair
   logical,                intent(out)   :: valid
!
! !LOCAL PARAMETERS:
   integer                               :: ivar
   type (type_model), pointer            :: model
   REALTYPE                              :: val
   character(len=256)                    :: err
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_CHECK_STATE_ _FABM_ARGS_ND_IN_,repair,valid

   valid = .true.

   ! Enumerate all non-container models in the model tree, and allow them to perform custom repairs.
   model => root%nextmodel
   do while (associated(model) .and. valid)
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%check_state(_INPUT_ARGS_CHECK_STATE_)
#endif
         ! ADD_NEW_MODEL_HERE - optional, only if the validity of state variable values cannot
         ! be checked simply by ascertaining whether its values lie within prescribed [constant] bounds.
         !
         ! Typical model call:
         ! call MODELNAME_check_state(model%MODELNAME,_INPUT_ARGS_CHECK_STATE_)
         !
         ! If the biogeochemical model provides a subroutine for checking the model state:
         ! Set "valid" to .false. if the state variable values are invalid,
         ! repair them if they are invalid and "repair" is .true. (but still keep "valid" set to .false.)
      end select

      ! If the present values are invalid and repair is not permitted, we are done.
      if (.not. (valid .or. repair)) return

      model => model%nextmodel
   end do

   ! Finally check whether all state variable values lie within their prescribed [constant] bounds.
   ! This is always done, independently of any model-specific checks that may have been called above.

   ! Enter spatial loops (if any)
   _FABM_LOOP_BEGIN_

   ! Check boundaries for pelagic state variables specified by the models.
   ! If repair is permitted, this clips invalid values to the closest boundary.
   do ivar=1,ubound(root%info%state_variables,1)
      _GET_STATE_EX_(root%environment,root%info%state_variables(ivar)%globalid,val)
      if (val<root%info%state_variables(ivar)%minimum) then
         ! State variable value lies below prescribed minimum.
         valid = .false.
         if (.not.repair) then
            write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',val,' of variable ',trim(root%info%state_variables(ivar)%name), &
                                                       & ' below minimum value ',root%info%state_variables(ivar)%minimum
            call log_message(err)
            return
         end if
         _SET_STATE_EX_(root%environment,root%info%state_variables(ivar)%globalid,root%info%state_variables(ivar)%minimum)
      elseif (val>root%info%state_variables(ivar)%maximum) then
         ! State variable value exceeds prescribed maximum.
         valid = .false.
         if (.not.repair) then
            write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',val,' of variable ',trim(root%info%state_variables(ivar)%name), &
                                                       & ' above maximum value ',root%info%state_variables(ivar)%maximum
            call log_message(err)
            return
         end if
         _SET_STATE_EX_(root%environment,root%info%state_variables(ivar)%globalid,root%info%state_variables(ivar)%maximum)
      end if
   end do

   ! Check boundaries for benthic state variables specified by the models.
   ! If repair is permitted, this clips invalid values to the closest boundary.
   do ivar=1,ubound(root%info%state_variables_ben,1)
      _GET_STATE_BEN_EX_(root%environment,root%info%state_variables_ben(ivar)%globalid,val)
      if (val<root%info%state_variables_ben(ivar)%minimum) then
         ! State variable value lies below prescribed minimum.
         valid = .false.
         if (.not.repair) then
            write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',val,' of variable ',trim(root%info%state_variables_ben(ivar)%name), &
                                                       & ' below minimum value ',root%info%state_variables_ben(ivar)%minimum
            call log_message(err)
            return
         end if
         _SET_STATE_BEN_EX_(root%environment,root%info%state_variables_ben(ivar)%globalid,root%info%state_variables_ben(ivar)%minimum)
      elseif (val>root%info%state_variables_ben(ivar)%maximum) then
         ! State variable value exceeds prescribed maximum.
         valid = .false.
         if (.not.repair) then
            write (unit=err,fmt='(a,e12.4,a,a,a,e12.4)') 'Value ',val,' of variable ',trim(root%info%state_variables_ben(ivar)%name), &
                                                       & ' above maximum value ',root%info%state_variables_ben(ivar)%maximum
            call log_message(err)
            return
         end if
         _SET_STATE_BEN_EX_(root%environment,root%info%state_variables_ben(ivar)%globalid,root%info%state_variables_ben(ivar)%maximum)
      end if
   end do

   ! Leave spatial loops (if any)
   _FABM_LOOP_END_

   end subroutine fabm_check_state
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-water exchange fluxes for all biogeochemical state variables.
! Positive values indicate fluxes into the ocean, negative values indicate fluxes
! out of the ocean. Units are tracer unit * m/s.
!
! !INTERFACE:
   subroutine fabm_get_surface_exchange(root _ARG_LOCATION_VARS_HZ_,flux)
!
! !INPUT PARAMETERS:
   type (type_model), intent(inout) :: root
   _DECLARE_LOCATION_ARG_HZ_
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE _ATTR_DIMENSIONS_1_HZ_,intent(out)   :: flux
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer       :: model
!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_GET_SURFACE_EXCHANGE_ _FABM_ARGS_IN_HZ_,flux
   flux = _ZERO_
   model => root%nextmodel
   do while (associated(model))
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%get_surface_exchange(_INPUT_ARGS_GET_SURFACE_EXCHANGE_)
#endif
         case (pml_carbonate_id)
            call pml_carbonate_get_surface_exchange(model%pml_carbonate,_INPUT_ARGS_GET_SURFACE_EXCHANGE_)
         case (iow_ergom_id)
            call iow_ergom_get_surface_exchange(model%iow_ergom,_INPUT_ARGS_GET_SURFACE_EXCHANGE_)
         ! ADD_NEW_MODEL_HERE - optional, only if the model specifies fluxes of one or
         ! more of its state variables across the air-water interface.
         !
         ! Typical model call:
         ! call MODELNAME_get_surface_exchange(model%MODELNAME,_INPUT_ARGS_GET_SURFACE_EXCHANGE_)
      end select
      model => model%nextmodel
   end do

   end subroutine fabm_get_surface_exchange
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Process interaction between benthos and bottom layer of the
! pelagic. This calculates the fluxes into all bottom pelagic and benthic variables,
! in variable quantity per surface area per time. This typically imples variable units * m/s
! [bottom fluxes] for the pelagic, and variable units/s [temporal derivatives] for the benthos.
! Positive values denote state variable increases, negative values state variable decreases.
!
! !INTERFACE:
   subroutine fabm_do_benthos_rhs(root _ARG_LOCATION_VARS_HZ_,flux_pel,flux_ben)
!
! !INPUT PARAMETERS:
   type (type_model),         intent(inout)    :: root
   _DECLARE_LOCATION_ARG_HZ_
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE _ATTR_DIMENSIONS_1_HZ_,intent(inout) :: flux_pel,flux_ben
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer               :: model
!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_DO_BENTHOS_RHS_ _FABM_ARGS_IN_HZ_,flux_pel,flux_ben

   model => root%nextmodel
   do while (associated(model))
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%do_benthos(_INPUT_ARGS_DO_BENTHOS_RHS_)
#endif
         case (examples_benthic_predator_id)
            call examples_benthic_predator_do_benthos(model%examples_benthic_predator,_INPUT_ARGS_DO_BENTHOS_RHS_)
         case(iow_ergom_id)
            call iow_ergom_do_benthos(model%iow_ergom,_INPUT_ARGS_DO_BENTHOS_RHS_)
         ! ADD_NEW_MODEL_HERE - optional, only if the model has benthic state variables,
         ! or specifies bottom fluxes for its pelagic state variables.
         !
         ! Typical model call:
         ! call MODELNAME_do_benthos(model%MODELNAME,_INPUT_ARGS_DO_BENTHOS_RHS_)
      end select
      model => model%nextmodel
   end do

   end subroutine fabm_do_benthos_rhs
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
   subroutine fabm_do_benthos_ppdd(root _ARG_LOCATION_VARS_HZ_,pp,dd,benthos_offset)
!
! !INPUT PARAMETERS:
   type (type_model),         intent(inout) :: root
   _DECLARE_LOCATION_ARG_HZ_
   integer,                   intent(in)    :: benthos_offset
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE _ATTR_DIMENSIONS_2_HZ_,intent(inout) :: pp,dd
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer               :: model
!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_DO_BENTHOS_PPDD_ _FABM_ARGS_IN_HZ_,pp,dd,benthos_offset

   model => root%nextmodel
   do while (associated(model))
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%do_benthos_ppdd(_INPUT_ARGS_DO_BENTHOS_PPDD_)
#endif
         ! ADD_NEW_MODEL_HERE - optional, only if the model has benthic state variables,
         ! or specifies bottom fluxes for its pelagic state variables.
         !
         ! Typical model call:
         ! call MODELNAME_do_benthos(model%MODELNAME,_INPUT_ARGS_DO_BENTHOS_PPDD_)
      end select
      model => model%nextmodel
   end do

   end subroutine fabm_do_benthos_ppdd
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the vertical movement rates (m/s) for the bio state variables.
! Note that negative values indicate movement towards the bottom, e.g., sinking,
! and positive values indicate movemment towards the surface, e.g., floating.
!
! !INTERFACE:
   subroutine fabm_get_vertical_movement(root _ARG_LOCATION_ND_,velocity)
!
! !INPUT PARAMETERS:
   type (type_model),           intent(inout) :: root
   _DECLARE_LOCATION_ARG_ND_
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE _ATTR_DIMENSIONS_1_,intent(out)  :: velocity
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer               :: model
   integer                                  :: i
   _TYPE_STATE_VARIABLE_ID_                 :: varid
!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_GET_VERTICAL_MOVEMENT_ _FABM_ARGS_ND_IN_,velocity

   model => root%nextmodel
   do while (associated(model))
      ! First set constant sinking rates.

      do i=1,ubound(model%info%state_variables,1)
         varid = model%info%state_variables(i)%globalid

         ! Enter spatial loops (if any)
         _FABM_LOOP_BEGIN_

         ! Use variable-specific constant vertical velocities.
         _SET_VERTICAL_MOVEMENT_(varid,model%info%state_variables(i)%vertical_movement)

         ! Leave spatial loops (if any)
         _FABM_LOOP_END_
      end do

      ! Now allow models to overwrite with spatially-varying sinking rates - if any.

      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%get_vertical_movement(_INPUT_ARGS_GET_VERTICAL_MOVEMENT_)
#endif
         ! ADD_NEW_MODEL_HERE - optional, only if the model specifies time- and/or space
         ! varying vertical velocities for one or more state variables.
         !
         ! Typical model call:
         ! call MODELNAME_get_vertical_movement(model%MODELNAME,_INPUT_ARGS_GET_VERTICAL_MOVEMENT_)
      end select
      model => model%nextmodel
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
   subroutine fabm_get_light_extinction(root _ARG_LOCATION_ND_,extinction)
!
! !INPUT PARAMETERS:
   type (type_model),           intent(inout) :: root
   _DECLARE_LOCATION_ARG_ND_
   REALTYPE _ATTR_DIMENSIONS_0_,intent(out)   :: extinction
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL PARAMETERS:
   REALTYPE                                   :: curext,val
   integer                                    :: i
   type (type_model), pointer                 :: model
!EOP
!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_GET_LIGHT_EXTINCTION_ _FABM_ARGS_ND_IN_,extinction

   extinction = _ZERO_
   model => root%nextmodel
   do while (associated(model))
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%get_light_extinction(_INPUT_ARGS_GET_LIGHT_EXTINCTION_)
#endif
         case (gotm_npzd_id)
            call gotm_npzd_get_light_extinction(model%gotm_npzd,_INPUT_ARGS_GET_LIGHT_EXTINCTION_)
         case (gotm_fasham_id)
            call gotm_fasham_get_light_extinction(model%gotm_fasham,_INPUT_ARGS_GET_LIGHT_EXTINCTION_)
         case (iow_ergom_id)
            call iow_ergom_get_light_extinction(model%iow_ergom,_INPUT_ARGS_GET_LIGHT_EXTINCTION_)
         ! ADD_NEW_MODEL_HERE - optional, only if light attenuation in the model cannot be captured by
         ! state variable specific extinction coefficients.
         !
         ! Typical model call:
         ! call MODELNAME_get_light_extinction(model%MODELNAME,_INPUT_ARGS_GET_LIGHT_EXTINCTION_)

         case default
            ! Default: use constant specific light extinction values specified in the state variable properties

            ! Enter spatial loops (if any)
            _FABM_LOOP_BEGIN_

            ! Use variable-specific light extinction coefficients.
            do i=1,ubound(model%info%state_variables,1)
               curext = model%info%state_variables(i)%specific_light_extinction
               if (curext.ne._ZERO_) then
                  _GET_STATE_EX_(root%environment,model%info%state_variables(i)%globalid,val)
                  _SET_EXTINCTION_(val*curext)
               end if
            end do

            ! Enter spatial loops (if any)
            _FABM_LOOP_END_
      end select
      model => model%nextmodel
   end do

   end subroutine fabm_get_light_extinction
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities
!
! !INTERFACE:
   subroutine fabm_get_conserved_quantities(root _ARG_LOCATION_ND_,sums)
!
! !INPUT PARAMETERS:
   type (type_model), intent(inout)         :: root
   _DECLARE_LOCATION_ARG_ND_
   REALTYPE _ATTR_DIMENSIONS_1_,intent(out) :: sums
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer               :: model

!-----------------------------------------------------------------------
!BOC
#define _INPUT_ARGS_GET_CONSERVED_QUANTITIES_ _FABM_ARGS_ND_IN_,sums

   sums = _ZERO_
   model => root%nextmodel
   do while (associated(model))
      select case (model%id)
#ifdef _FABM_F2003_
         case (model_f2003_id)
            call model%info%get_conserved_quantities(_INPUT_ARGS_GET_CONSERVED_QUANTITIES_)
#endif
         case (gotm_npzd_id)
            call gotm_npzd_get_conserved_quantities(model%gotm_npzd,_INPUT_ARGS_GET_CONSERVED_QUANTITIES_)
         case (gotm_fasham_id)
            call gotm_fasham_get_conserved_quantities(model%gotm_fasham,_INPUT_ARGS_GET_CONSERVED_QUANTITIES_)
         ! ADD_NEW_MODEL_HERE - optional, required only if the model exports one or more
         ! conserved quantities.
         !
         ! Typical model call:
         ! call MODELNAME_get_conserved_quantities(model%MODELNAME,_INPUT_ARGS_GET_CONSERVED_QUANTITIES_)

         case default
            ! Default: the model does not describe any conserved quantities.
            if (ubound(model%info%conserved_quantities,1).gt.0) &
               call fatal_error('fabm_get_conserved_quantities','model '//trim(model%info%name)//' registered &
                    &one or more conserved quantities, but a function that provides sums of these &
                    &quantities has not been provided.')
      end select
      model => model%nextmodel
   end do

   end subroutine fabm_get_conserved_quantities
!EOC

end module fabm

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

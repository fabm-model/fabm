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
! For the impatient: throughout this file, location where additional bio models may be
! referenced are indicated by the comment "ADD_NEW_MODEL_HERE". You can search for this
! string to get started quickly.
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
!    Grouping the model parameters in a derived type then is a *requirement*.
!
! 4) Add the model as option to the "select" statements in the following subroutines:
!    "get_model_name" - to provide a short name (string) for your model.
!    "init_rmbm"      - to initialize your model: register prognostic and diagnostic variables,
!                       read in namelists, etc.
!    "do_rmbm_rhs"    - to provide the local temporal derivatives of your model.
!
! The following steps are optional:
!
! 5) If the sinking rate of any of the model variables varies in time and/or space, a subroutine
!    that provides the sinking rates (m/s) must be added as option to the "select" statement in
!    "get_vertical_movement_rmbm".

!    Note that in that case, your model must also set member "dynamic_vertical_movement" in the
!    model information to a value of 2 during initialization (otherwise your function will
!    not be called).
!
!    If a function is not provided (dynamic_vertical_movement=0), sinking rates are assumed to be constant
!    in time and space; they will be taken from the vertical_movement member of the respective
!    type_state_variable_info derived type (see bio_types.F90).
!
! 6) If any of the model variables attenuate light, a function that calculates the light
!    extinction coefficient (/m) from the current model state must be added as option to the "select"
!    statement in get_bio_extinction_rmbm. This allows for complete customization of the bio
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
! 1) Add a use statement that references this module, e.g.:
!
!    "use rmbm"
!
! 2) Define a pointer to an instance of the derived type "type_model":
!
!    "type (model_type), pointer :: model"
!
!    This instance will hold all information on the selected biogeochemical model,
!    including descriptive strings, state variable information, and parameter values.
!
! 3) Create the root model instance by calling create_model. Typically, the root model
!    will be a model container that is created by calling create_model without arguments:
!
!    "model => rmbm_create_model()
!
!    After the root model is created, submodels can be added to it by repeatedly calling
!    create_model, with the existing root model as parent argument:
!
!    "childmodel => create_model(1,parent=model)"
!    "childmodel => create_model(2,parent=model)"
!
!    This adds submodels with identifiers 1 and 2 to the root.
!
! 4) Initialize the model tree by calling "init_rmbm" on the root model:
!
!    "init_rmbm(model,nmlunit)"
!
!    with "nmlunit" being a unit specifier (integer) of a file that has already been openend.
!    The RMBM models will read namelists from this file during initialization.
!    After initialization the caller is responsible for closing the file.
!
! 5) Access the model by the following subroutines:
!    - do_rmbm_rhs:                   to get local temporal derivatives
!    - get_vertical_movement_rmbm:    to get current sinking rates for the state variables
!    - get_bio_extinction_rmbm:       to get the combined light extinction coefficient due to
!                                     biogeochemical components
!    - get_conserved_quantities_rmbm: to get the sums of the conserved quantities described
!                                     by the model
!
!    Additional information on the model (e.g. descriptive string for its variables) is present in the
!    "info" member of the derived type, after the model (collection) has been initialized.
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
   
      integer                :: id
      character(len=64)      :: name
      type (type_model_info) :: info
      
      ! Pointers for linking to parent and child models
      type (type_model),pointer :: parent      => NULL()
      type (type_model),pointer :: firstchild  => NULL()
      type (type_model),pointer :: nextsibling => NULL()

      ! Pointers to the current state and environment.
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
         call fatal_error('bio_0d_gen::get_model_name',text)
   end select
   
   end function get_model_name
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Create a new model. This will be a specific model if the
! model identifier is provided, or a contained of child models if the
! idnetifier is omitted.
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

   ! Connect to parent model if specified.
   if (present(parent)) then
      if (parent%id.ne.model_container_id) &
         call fatal_error('rmbm::add_child_model','A child model can only be added to a container, not to an existing model.')

      if (associated(parent%firstchild)) then
         previoussibling => parent%firstchild
         do while (associated(previoussibling%nextsibling))
            previoussibling => previoussibling%nextsibling
         end do
         previoussibling%nextsibling => model
      else
         parent%firstchild => model
      end if
      model%parent => parent
   end if
   
   end function rmbm_create_model
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the selected 0d biogeochemical model
!
! !INTERFACE:
   recursive subroutine rmbm_init(model,nmlunit,nameprefix,longnameprefix,master)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),                      intent(inout) :: model
   integer,                                intent(in)    :: nmlunit
   character(len=*),      optional,        intent(in)    :: nameprefix,longnameprefix
   type (type_model_info),optional,target, intent(in)    :: master
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer                    :: i,ichild
  character(len= 64)         :: modelname
  type (type_model), pointer :: curchild,curchild2
  character(len=256)         :: childnameprefix,childlongnameprefix,childname
!EOP
!-----------------------------------------------------------------------
!BOC
   modelname = get_model_name(model%id)
   call log_message('Initializing biogeochemical model '//trim(modelname))

   call init_model_info(model%info)
   if (present(master        )) model%info%master => master
   if (present(nameprefix    )) model%info%nameprefix = nameprefix
   if (present(longnameprefix)) model%info%longnameprefix = longnameprefix

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
            end do
            write (unit=childname, fmt='(a,i2.2)') trim(curchild%name),ichild
         
            ! Create prefixes for names of variables of the child model.
            if (present(nameprefix)) then
              childnameprefix = nameprefix//trim(childname)//'_'
            else
              childnameprefix = trim(childname)//'_'
            end if
            if (present(longnameprefix)) then
              childlongnameprefix = longnameprefix//trim(childname)//' '
            else
              childlongnameprefix = trim(childname)//' '
            end if

            ! Initialize child model.
            call rmbm_init(curchild,nmlunit,childnameprefix,childlongnameprefix,model%info)
            
            ! Move to next child model.
            curchild => curchild%nextsibling
         end do
         
      case default
         call fatal_error('bio_0d_gen::init_bio_single','no valid biogeochemical model specified!')
         
   end select
   
   if (.not. associated(model%parent)) then
      allocate(model%environment)
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
   LOCATIONTYPE,           intent(in)    :: LOCATION
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
   recursive subroutine rmbm_do_rhs_1d(model,LOCATION_1DLOOP,dy,diag)
!
! !USES:
   LOOP1D_USE
   
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATIONTYPE,           intent(in)    :: LOCATION_1DLOOP
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: dy(:,:)
   REALTYPE,               intent(inout) :: diag(:,:)
!
! !LOCAL PARAMETERS:
   type (type_model), pointer            :: curchild
   LOCATIONTYPE                          :: VARIABLE_1DLOOP
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
            call rmbm_do_rhs_1d(curchild,LOCATION_1DLOOP,dy,diag)
            curchild => curchild%nextsibling
         end do
      case default
         do VARIABLE_1DLOOP=1,LENGTH_1DLOOP
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
   LOCATIONTYPE,           intent(in)    :: LOCATION
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
   LOCATIONTYPE,           intent(in)    :: LOCATION
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
! !IROUTINE: Get the air-sea exchange fluxes
!
! !INTERFACE:
   recursive subroutine rmbm_update_air_sea_exchange(model,LOCATION,numc,flux)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATIONTYPE,           intent(in)    :: LOCATION
   integer,                intent(in)    :: numc
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: flux(1:numc)
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
         call update_air_sea_co2_sys_0d(model%carbonate,model%state,model%environment,LOCATION,numc,flux)
      ! ADD_NEW_MODEL_HERE - optional
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            call rmbm_update_air_sea_exchange(curchild,LOCATION,numc,flux)
            curchild => curchild%nextsibling
         end do
      case default
   end select

   end subroutine rmbm_update_air_sea_exchange
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the sinking rates for the state variables of the
! selected 0d biogeochemical model
!
! !INTERFACE:
   recursive subroutine rmbm_get_vertical_movement(model,LOCATION,vertical_movement)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATIONTYPE,           intent(in)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: vertical_movement(:)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!
   type (type_model), pointer            :: curchild
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
         if (model%info%dynamic_vertical_movement.ne.0) &
            call fatal_error('get_vertical_movement_single','the 0d model specifies that vertical movement is time- and/or &
                 &space-dependent, but a function that provides these sinking rates has not been specified.')
         vertical_movement = model%info%variables%vertical_movement
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
   LOCATIONTYPE,      intent(in) :: LOCATION
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
   recursive subroutine rmbm_get_conserved_quantities(model,LOCATION,num,sums)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model), intent(in)      :: model
   LOCATIONTYPE,      intent(in)      :: LOCATION
   integer,           intent(in)      :: num
   REALTYPE,          intent(inout)   :: sums(1:num)
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
         call get_conserved_quantities_npzd_0d(model%npzd,model%state,model%environment,LOCATION,num,sums)
      ! ADD_NEW_MODEL_HERE - optional
      case (model_container_id)
         curchild => model%firstchild
         do while (associated(curchild))
            call rmbm_get_conserved_quantities(curchild,LOCATION,num,sums)
            curchild => curchild%nextsibling
         end do
      case default
         ! Default: the model does not describe any conserved quantities.
         if (model%info%conserved_quantity_count.gt.0) &
            call fatal_error('get_conserved_quantities_single','the model specifies that it describes one or more conserved &
                 &quantities, but a function that provides sums of these quantities has not been specified.')
   end select

   end subroutine rmbm_get_conserved_quantities
!EOC

end module rmbm

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

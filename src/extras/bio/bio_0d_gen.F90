!$Id: bio_0d_gen.F90,v 1.8 2009-05-11 13:57:31 jorn Exp $
#include"cppdefs.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_0d_gen --- Generic 0D biogeochemical model
!
! !INTERFACE:
   module bio_0d_gen
!
! !DESCRIPTION:
! This module encapsulates all specific 0d biogeochemical models. Each
! subroutine and function is called with an identifier [bio_model] specifying
! the desired specific biogeochemical model; based on this identifier
! the correct subroutine or function of the specific model will be called.
!
! A new 0D biogeochemical models only need to be registered in this file
! in order to be usable by GOTM, GETM and the independent 0D driver.
!
! How to register a new biogeochemical model:
!
! 1) Add a use statement that references your new model (cf. "use bio_npzd_0d" below)
!
! 2) Define an integer identifier > 1000 for your new model (cf. "npzd_0d_id" below)
!
! 3) If your model uses model-specific data (e.g. parameter values) that are grouped in a
!    (model-specific) derived type, add an instance of this type as member to "type_model"
!    defined below. Note: if you want to be able to run multiple instances of your model
!    side-by-side, grouping the model parameters in a derived type is a *requirement*.
!
! 4) Add the model as option to the "select" statements in the following subroutines:
!    "get_model_name", "init_bio_single", "do_bio_single".
!
! The following steps are optional:
!
! 5) If the sinking rate of any of the model variables varies in time and/or space, a subroutine
!    that provides the sinking rates (m/s) must be added as option to the "select" statement in
!    "get_vertical_movement_single". Note that in that case, your model must also set member
!    "dynamic_vertical_movement" in the model information to a value of 2 (otherwise your function will
!    not be called).
!    If a function is not provided (dynamic_vertical_movement=0), sinking rates are assumed to be constant
!    in time and space; they will be taken from the vertical_movement member of the respective
!    type_state_variable_info derived type (see bio_types.F90).
!
! 6) If any of the model variables attentuate light, a function that calculate the light
!    extinction coefficient (/m) from the current model state may be added as option to the "select"
!    statement in get_bio_extinction_single. This allows for complete customization of the bio
!    extinction.
!    If this function is not provided, the bio extinction will be calculated from the specific
!    extinction coefficients for each state variable as specified in the model information,
!    specifically member "specific_light_extinction" of the state variable information.
!    (note: specific extinction coefficients default to 0: no attenuation due to biogeochemical components) 
!
! 7) If (part of) the model state variable are composed of conserved quantities (energy and/or
!    chemical elements), a function that provides the sum of these quantities given the model
!    state must be added as option to the "select" statement in get_conserved_quantities_single.
!
! How to use this library of biogeochemical models:
!
! 1) Add a use statement that references this module.
!
! 2) Define an instance of the derived type "type_model" if you want to use one single model, or
!    of type "type_model_collection" if you want to use multiple models side-by-side (N.B. these
!    models will not be aware of each other).
!    In both cases this instance will hold all information on the selected biogeochemical model,
!    including descriptive strings, state variable information, and parameter values.
!
! 3) Initialize the derived type instance by calling "init_bio_0d_generic"
!
! 4) Access the model by the following subroutines:
!    - do_bio_0d_generic: to get local temporal derivatives
!    - get_vertical_movement_bio_0d_generic: to get current sinking rates for the state variables
!    - get_bio_extinction_bio_0d_generic: to get the combined light extinction coefficient due to
!      biogeochemical components
!    - get_conserved_quantities_bio_0d_generic: to get the sums of the conserved quantities described
!      by the model
!    Additional information on the model (e.g. descriptive string for its variables) is present in the
!    "info" member of the derived type, after the model (collection) has been initialized.
!
! !USES:
   use bio_types
   use bio_driver
!   
!  Reference specific biogeochemical models:
   use bio_npzd_0d
   use bio_jellyfish_0d
   use bio_co2_sys_0d
!   
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_model, type_model_collection, &
          get_model_name, init_bio_0d_generic,do_bio_0d_generic, &
          get_vertical_movement_bio_0d_generic, get_bio_extinction_bio_0d_generic, &
          get_conserved_quantities_bio_0d_generic, update_airsea_exchange_bio_0d_generic
!
! !PRIVATE DATA MEMBERS:

!  Identifiers for specific biogeochemical models:
   integer, parameter :: no_model_id = -1
   integer, parameter :: npzd_0d_id = 1001
   integer, parameter :: jellyfish_0d_id = 1002
   integer, parameter :: carbonate_id = 1003

! !PUBLIC TYPES:
!
!  Single generic biogeochemical model
   type type_model
      integer :: id
      type (type_model_info) :: info
      
      ! Derived types that below to specific biogeochemical models.
      type (type_npzd)      :: npzd
      type (type_jellyfish) :: jellyfish
      type (type_co2_sys)   :: carbonate
   end type type_model
!
!  Collection of generic biogeochemical models.
   type type_model_collection
      integer :: count
      type (type_model_info) :: info
      type (type_model),allocatable :: models(:)
   end type type_model_collection
!
! !PUBLIC INTERFACES:
!
!  Model initialization
   interface init_bio_0d_generic
      module procedure init_bio_single
      module procedure init_bio_collection
   end interface init_bio_0d_generic
!
!  Access to temporal derivatives
   interface do_bio_0d_generic
      module procedure do_bio_single
      module procedure do_bio_single_rhs
      module procedure do_bio_collection
      module procedure do_bio_collection_rhs
   end interface do_bio_0d_generic
!
!  Access to sinking rates
   interface get_vertical_movement_bio_0d_generic
      module procedure get_vertical_movement_single
      module procedure get_vertical_movement_collection
   end interface get_vertical_movement_bio_0d_generic
!
!  Access to light extinction coefficient
   interface get_bio_extinction_bio_0d_generic
      module procedure get_bio_extinction_single
      module procedure get_bio_extinction_collection
   end interface get_bio_extinction_bio_0d_generic
!
!  Access to sums of conserved quantities
   interface get_conserved_quantities_bio_0d_generic
      module procedure get_conserved_quantities_single
      module procedure get_conserved_quantities_collection
   end interface get_conserved_quantities_bio_0d_generic
!
!  Access to air-sea exchange
   interface update_airsea_exchange_bio_0d_generic
      module procedure update_airsea_exchange_single
      module procedure update_airsea_exchange_collection
   end interface update_airsea_exchange_bio_0d_generic
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
!-----------------------------------------------------------------------
!BOC
   select case (bio_model)
      case (npzd_0d_id)
         name = 'npzd'
      case (jellyfish_0d_id)
         name = 'jellyfish'
      case (carbonate_id)
         name = 'co2sys'
      case default
         stop 'bio_0d_gen::get_model_name: no valid biogeochemical model specified!'
   end select
   
   end function get_model_name
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the selected 0d biogeochemical model
!
! !INTERFACE:
   function init_bio_single(bio_model,nmlunit,nmlfilename,nameprefix,longnameprefix,master) result(model)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,                                intent(in) :: bio_model,nmlunit
   character(len=*),                       intent(in) :: nmlfilename
   character(len=*),      optional,        intent(in) :: nameprefix,longnameprefix
   type (type_model_info),optional,target, intent(in) :: master
   type (type_model)                                  :: model
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer           :: i
  character(len=64) :: modelname
!EOP
!-----------------------------------------------------------------------
!BOC
   modelname = get_model_name(bio_model)
   LEVEL2 'Initializing biogeochemical model '//trim(modelname)

   call init_model_info(model%info)
   if (present(master        )) model%info%master => master
   if (present(nameprefix    )) model%info%nameprefix = nameprefix
   if (present(longnameprefix)) model%info%longnameprefix = longnameprefix

   ! Preset model identifier to none, in case initialization
   model%id = no_model_id

   ! Allow the selected model to initialize
   open(nmlunit,file=nmlfilename,action='read',status='old',err=98)
   select case (bio_model)
      case (no_model_id)
      case (npzd_0d_id)
         call init_bio_npzd_0d(model%npzd,model%info,nmlunit)
      case (jellyfish_0d_id)
         call init_bio_jellyfish_0d(model%jellyfish,model%info,nmlunit)
      case (carbonate_id)
         call init_bio_co2_sys_0d(model%carbonate,model%info,nmlunit)
      case default
         stop 'bio_0d_gen::init_bio_single: no valid biogeochemical model specified!'
   end select
   close(nmlunit)
   LEVEL3 'model '//trim(modelname)//' initialized successfully from '//trim(nmlfilename)
   
   ! Store the identifier for the selected model.
   model%id = bio_model
   
   return
   
98 LEVEL2 'I could not open '//trim(nmlfilename)
   LEVEL2 'If thats not what you want you have to supply '//trim(nmlfilename)
   LEVEL2 'See the bio example on www.gotm.net for a working '//trim(nmlfilename)
   return

   end function init_bio_single
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the selected 0d
! biogeochemical model
!
! !INTERFACE:
   subroutine do_bio_single(model,LOCATIONVARIABLE,numc,num_diag,pp,dd,diag)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATIONTYPE,           intent(in)    :: LOCATIONVARIABLE
   integer,                intent(in)    :: numc,num_diag
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: pp(1:numc,1:numc)
   REALTYPE,               intent(inout) :: dd(1:numc,1:numc)
   REALTYPE,               intent(inout) :: diag(1:num_diag)
!
! !LOCAL PARAMETERS:
   REALTYPE,allocatable                  :: dy(:)
   integer                               :: i
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (no_model_id)
      case (npzd_0d_id)
         call do_bio_npzd_0d_ppdd(model%npzd,LOCATIONVARIABLE,numc,num_diag,pp,dd,diag)
      case default
         allocate(dy(numc))
         dy = _ZERO_
         call do_bio_single_rhs(model,LOCATIONVARIABLE,numc,num_diag,dy,diag)
         do i=1,numc
            pp(i,i) = pp(i,i) + dy(i)
         end do
         deallocate(dy)
   end select

   end subroutine do_bio_single
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the selected 0d
! biogeochemical model
!
! !INTERFACE:
   subroutine do_bio_single_rhs(model,LOCATIONVARIABLE,numc,num_diag,dy,diag)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATIONTYPE,           intent(in)    :: LOCATIONVARIABLE
   integer,                intent(in)    :: numc,num_diag
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: dy(1:numc)
   REALTYPE,               intent(inout) :: diag(1:num_diag)
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
   select case (model%id)
      case (no_model_id)
      case (npzd_0d_id)
         call do_bio_npzd_0d(model%npzd,LOCATIONVARIABLE,numc,num_diag,dy,diag)
      case (jellyfish_0d_id)
         call do_bio_jellyfish_0d(model%jellyfish,LOCATIONVARIABLE,numc,num_diag,dy,diag)
      case (carbonate_id)
         call do_bio_co2_sys_0d(model%carbonate,LOCATIONVARIABLE,numc,num_diag,dy,diag)
      case default
         allocate(pp(1:numc,1:numc))
         allocate(dd(1:numc,1:numc))
         pp = _ZERO_
         dd = _ZERO_
         call do_bio_single(model,LOCATIONVARIABLE,numc,num_diag,pp,dd,diag)
         do i=1,numc
            do j=1,numc
               dy(i) = dy(i) + pp(i,j)-dd(i,j)
            end do
         end do
         deallocate(pp)
         deallocate(dd)
   end select

   end subroutine do_bio_single_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-sea exchange fluxes
!
! !INTERFACE:
   subroutine update_airsea_exchange_single(model,LOCATIONVARIABLE,numc,flux)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATIONTYPE,           intent(in)    :: LOCATIONVARIABLE
   integer,                intent(in)    :: numc
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: flux(1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (carbonate_id)
         call update_airsea_co2_sys_0d(model%carbonate,LOCATIONVARIABLE,numc,flux)
   end select

   end subroutine update_airsea_exchange_single
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the sinking rates for the state variables of the
! selected 0d biogeochemical model
!
! !INTERFACE:
   subroutine get_vertical_movement_single(model,LOCATIONVARIABLE,numc,vertical_movement)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model),      intent(in)    :: model
   LOCATIONTYPE,           intent(in)    :: LOCATIONVARIABLE
   integer,                intent(in)    :: numc
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,               intent(inout) :: vertical_movement(1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case default
         ! Default: use the constant sinking rates specified in state variable properties.
         if (model%info%dynamic_vertical_movement.ne.0) &
            stop 'get_vertical_movement_single: the 0d model specifies that vertical movement is time- and/or &
                 &space-dependent, but a function that provides these sinking rates has not been specified.'
         vertical_movement = model%info%variables%vertical_movement
   end select

   end subroutine get_vertical_movement_single
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   function get_bio_extinction_single(model,LOCATIONVARIABLE) result(extinction)
!
! !INPUT PARAMETERS:
   type (type_model), intent(in) :: model
   LOCATIONTYPE,      intent(in) :: LOCATIONVARIABLE
   REALTYPE                      :: extinction
   
   integer                       :: i
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (npzd_0d_id)
         extinction = get_bio_extinction_npzd_0d(model%npzd,LOCATIONVARIABLE)
      case default
         ! Default: use constant specific light extinction values specified in the state variable properties
         extinction = _ZERO_
         do i=1,model%info%state_variable_count
            if (model%info%variables(i)%specific_light_extinction.ne._ZERO_) &
               extinction = extinction + getbiovar(model%info%variables(i)%globalid,LOCATIONVARIABLE) &
                                 & * model%info%variables(i)%specific_light_extinction
         end do
   end select

   end function get_bio_extinction_single
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities
!
! !INTERFACE:
   subroutine get_conserved_quantities_single(model,LOCATIONVARIABLE,num,sums)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model), intent(in)      :: model
   LOCATIONTYPE,      intent(in)      :: LOCATIONVARIABLE
   integer,           intent(in)      :: num
   REALTYPE,          intent(inout)   :: sums(1:num)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (model%id)
      case (npzd_0d_id)
         call get_conserved_quantities_npzd_0d(model%npzd,LOCATIONVARIABLE,num,sums)
      case default
         ! Default: the model does not describe any conserved quantities.
         if (model%info%conserved_quantity_count.gt.0) &
            stop 'get_conserved_quantities_single: the model specifies that it describes one or more conserved &
                 &quantities, but a function that provides sums of these quantities has not been specified.'
   end select

   end subroutine get_conserved_quantities_single
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise a collection of 0d biogeochemical models
!
! !INTERFACE:
   function init_bio_collection(count,bio_model,nmlunit,nmlfilename,nameprefixes,longnameprefixes) result(collection)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer,                     intent(in) :: count,bio_model(count),nmlunit
   character(len=*),            intent(in) :: nmlfilename(count)
   character(len=*),optional,   intent(in) :: nameprefixes(count),longnameprefixes(count)
   type (type_model_collection),target     :: collection
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer :: i,nstate,ndiag,ncons
  character(len=64) :: nameprefix,longnameprefix,modelname
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Allocate the memory for the set of biogeochemical models in the collection.
   collection%count = count
   allocate(collection%models(1:count))
   
   ! Initialize empty model info for the collection.
   call init_model_info(collection%info)

   ! Allow each biogeochemical model to initialize.
   do i = 1,count
      ! Generate prefix strings if not specified
      modelname = get_model_name(bio_model(i))
      if (present(nameprefixes)) then
         nameprefix = nameprefixes(i)
      else
         write (unit=nameprefix, fmt='(a,i2.2,a)') trim(modelname),i,'_'
      end if
      if (present(longnameprefixes)) then
         longnameprefix = longnameprefixes(i)
      else
         write (unit=longnameprefix, fmt='(a,a,i2.2)') trim(modelname),' ',i
      end if
      
      ! Initialize the current model
      collection%models(i) = init_bio_0d_generic(bio_model(i),nmlunit,nmlfilename(i),nameprefix,longnameprefix,master=collection%info)
   end do

   ! Create table for information on the entire collection of models.
   collection%info%dynamic_vertical_movement = maxval(collection%models%info%dynamic_vertical_movement)
   
   end function init_bio_collection
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the collection of biogeochemical models.
!
! !INTERFACE:
   subroutine do_bio_collection(collection,LOCATIONVARIABLE,numc,num_diag,pp,dd,diag)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model_collection), intent(in)    :: collection
   LOCATIONTYPE,                 intent(in)    :: LOCATIONVARIABLE
   integer,                      intent(in)    :: numc,num_diag
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,                     intent(inout) :: pp(1:numc,1:numc)
   REALTYPE,                     intent(inout) :: dd(1:numc,1:numc)
   REALTYPE,                     intent(inout) :: diag(1:num_diag)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   if (numc.ne.collection%info%state_variable_count) then
      FATAL numc,collection%info%state_variable_count
      stop 'do_bio_collection'
   end if
   
   do i = 1,collection%count
      call do_bio_single(collection%models(i),LOCATIONVARIABLE,numc,num_diag,pp,dd,diag)
   end do

   end subroutine do_bio_collection
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the local temporal derivatives for the collection of biogeochemical models.
!
! !INTERFACE:
   subroutine do_bio_collection_rhs(collection,LOCATIONVARIABLE,numc,num_diag,dy,diag)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model_collection), intent(in)    :: collection
   LOCATIONTYPE,                 intent(in)    :: LOCATIONVARIABLE
   integer,                      intent(in)    :: numc,num_diag
!
! !INPUT/OUTPUT PARAMETERS:
   REALTYPE,                     intent(inout) :: dy(1:numc)
   REALTYPE,                     intent(inout) :: diag(1:num_diag)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   if (numc.ne.collection%info%state_variable_count) then
      FATAL numc,collection%info%state_variable_count
      stop 'do_bio_collection'
   end if
   
   do i = 1,collection%count
      call do_bio_single_rhs(collection%models(i),LOCATIONVARIABLE,numc,num_diag,dy,diag)
   end do

   end subroutine do_bio_collection_rhs
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the sinking rates for the state variables of the
! collection of biogeochemical models
!
! !INTERFACE:
   subroutine get_vertical_movement_collection(collection,LOCATIONVARIABLE,vertical_movement)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model_collection), intent(in)  :: collection
   LOCATIONTYPE,                 intent(in)  :: LOCATIONVARIABLE
   REALTYPE,                     intent(out) :: vertical_movement(1:collection%info%state_variable_count)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   do i = 1,collection%count
      call get_vertical_movement_single(collection%models(i),LOCATIONVARIABLE,collection%info%state_variable_count,vertical_movement)
   end do

   end subroutine get_vertical_movement_collection
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   function get_bio_extinction_collection(collection,LOCATIONVARIABLE) result(extinction)
!
! !INPUT PARAMETERS:
   type (type_model_collection), intent(in) :: collection
   LOCATIONTYPE,                 intent(in) :: LOCATIONVARIABLE
   REALTYPE                                 :: extinction
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!EOP
!-----------------------------------------------------------------------
!BOC
   extinction = _ZERO_
   do i = 1,collection%count
      extinction = extinction + get_bio_extinction_bio_0d_generic(collection%models(i),LOCATIONVARIABLE)
   end do

   end function get_bio_extinction_collection
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the total of all conserved quantities for a collection
! of models.
!
! !INTERFACE:
   subroutine get_conserved_quantities_collection(collection,LOCATIONVARIABLE,sums)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model_collection), intent(in)  :: collection
   LOCATIONTYPE,                 intent(in)  :: LOCATIONVARIABLE
   REALTYPE,                     intent(out) :: sums(1:collection%info%conserved_quantity_count)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   do i = 1,collection%count
      call get_conserved_quantities_single(collection%models(i),LOCATIONVARIABLE, &
              collection%info%conserved_quantity_count,sums)
   end do

   end subroutine get_conserved_quantities_collection
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the air-sea exchange fluxes for a collection of models.
!
! !INTERFACE:
   subroutine update_airsea_exchange_collection(collection,LOCATIONVARIABLE,numc,flux)
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   type (type_model_collection), intent(in)    :: collection
   LOCATIONTYPE,                 intent(in)    :: LOCATIONVARIABLE
   integer,                      intent(in)    :: numc
   REALTYPE,                     intent(inout) :: flux(1:numc)
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
! !LOCAL VARIABLES:
  integer :: i
!EOP
!-----------------------------------------------------------------------
!BOC
   do i = 1,collection%count
      call update_airsea_exchange_single(collection%models(i),LOCATIONVARIABLE,numc,flux)
   end do

   end subroutine update_airsea_exchange_collection
!EOC

   end module bio_0d_gen

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

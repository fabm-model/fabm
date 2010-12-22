!$Id$
#include "rmbm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: rmbm_types --- Derived types used by 0D biogeochemical modules
!
! !INTERFACE:
   module rmbm_types
   
   use rmbm_driver, only: fatal_error
   
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
   public type_model_info
   public type_state_variable_info,type_diagnostic_variable_info,type_conserved_quantity_info
   public init_model_info
   public register_state_variable, register_diagnostic_variable, register_conserved_quantity, &
          register_state_dependency, register_dependency
   public type_environment,type_state_hz,type_state,type_state_variable_id
!
! !PUBLIC DERIVED TYPES:
!

   type type_state_variable_id
      integer :: id
      REALTYPE :: scale
   end type type_state_variable_id

   ! Properties of a single state variable
   type type_state_variable_info
      character(len=64) :: name, longname, units

      REALTYPE :: initial_value              ! Initial state variable value
      REALTYPE :: minimum,maximum            ! Valid range
      REALTYPE :: vertical_movement          ! Vertical movement (m/s) due to e.g. sinking, floating or activity. Note: positive for upward movement!
      REALTYPE :: specific_light_extinction  ! Specific light extinction (/m/state variable unit)
#if 0
      REALTYPE :: mussels_inhale             ! Whether this variable can be consumed by mussels
#endif
      logical :: no_precipitation_dilution,no_river_dilution
      
      integer  :: dependencyid      ! This is a globally unique identifier for the variable that can be used to retrieve values.
      type (type_state_variable_id) :: globalid
      integer  :: id
   end type type_state_variable_info

   ! Properties of a diagnostic variable
   type type_diagnostic_variable_info
      character(len=64) :: name, longname, units
      integer           :: id,dependencyid
      
      ! Time treatment:
      ! 0: last value
      ! 1: time-integrated
      ! 2: time step-averaged
      ! 3: time step-integrated
      integer           :: time_treatment
   end type type_diagnostic_variable_info
   
   integer, parameter,public  :: time_treatment_last=0,time_treatment_integrated=1, &
                                 time_treatment_averaged=2,time_treatment_step_integrated=3

   ! Properties of a conserved quantity
   type type_conserved_quantity_info
      character(len=64) :: name, longname, units
      integer           :: globalid,id
   end type type_conserved_quantity_info
   
   ! Global 0D model properties
   type type_model_info
      type (type_state_variable_info),     pointer,dimension(:) :: state_variables_ben,state_variables
      type (type_diagnostic_variable_info),pointer,dimension(:) :: diagnostic_variables_hz,diagnostic_variables
      type (type_conserved_quantity_info), pointer,dimension(:) :: conserved_quantities
      
      type (type_model_info),pointer :: parent
      type (type_model_info),pointer :: firstchild
      type (type_model_info),pointer :: nextsibling
      
      character(len=64) :: name,nameprefix,longnameprefix
      
      character(len=64),pointer :: dependencies(:)
      character(len=64),pointer :: dependencies_hz(:)
   end type type_model_info
   
   ! Parameters
   integer, parameter, public         :: shape_hz=2,shape_full=3
   
   integer, parameter, public         :: id_not_used=-1
   
   character(len=64),parameter,public :: &
     varname_temp    = 'env_temp',    & ! Temperature (degrees Celsius)
     varname_salt    = 'env_salt',    & ! Salinity (psu)
     varname_par     = 'env_par',     & ! Photosynthetically Active Radiation (W/m^2)
     varname_pres    = 'env_pres',    & ! Pressure (dbar = 10 kPa)
     varname_dens    = 'env_dens',    & ! Density (kg/m^3)
     varname_wind_sf = 'env_wind_sf', & ! Wind speed at 10 m above surface (m/s)
     varname_par_sf  = 'env_par_sf'     ! Photosynthetically Active Radiation at surface (W/m^2)
                                  
   type type_state
      REALTYPE,pointer _ATTR_LOCATION_DIMENSIONS_ :: data
   end type type_state

   type type_state_hz
      REALTYPE,pointer _ATTR_LOCATION_DIMENSIONS_HZ_ :: data
   end type type_state_hz

   type type_environment

      ! Pointer(s) to arrays that will hold state variable values.
#ifdef RMBM_SINGLE_STATE_VARIABLE_ARRAY
      ! Data for all state variables are stored in a single array.
      REALTYPE,pointer _ATTR_LOCATION_DIMENSIONS_PLUS_ONE_    :: state     ! pointer to data of pelagic state variables
      REALTYPE,pointer _ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE_ :: state_ben ! pointer to data of benthic state variables
#else
      ! Data for the state variables may be stored in different arrays.
      type (type_state   ), dimension(:), _ALLOCATABLE_ :: state     _NULL_ ! array of pointers to data of pelagic state variables
      type (type_state_hz), dimension(:), _ALLOCATABLE_ :: state_ben _NULL_ ! array of pointers to data of benthic state variables
#endif

      ! Pointer(s) to arrays that will hold values of "generic" variables, that is,
      ! all internal and external dependencies and, if _RMBM_MANAGE_DIAGNOSTICS_ is not set, diagnsotic variables as well.
      type (type_state   ), dimension(:), _ALLOCATABLE_ :: var    _NULL_ ! array of pointers to data of all pelagic variables (state and diagnostic)
      type (type_state_hz), dimension(:), _ALLOCATABLE_ :: var_hz _NULL_ ! array of pointers to data of all horizontal variables (state and diagnostic, surface and bottom)
      
#ifdef _RMBM_MANAGE_DIAGNOSTICS_
      ! RMBM will manage the current value of diagnostic variables itself.
      ! Declare the arrays for this purpose.
      REALTYPE,_ALLOCATABLE_ _ATTR_LOCATION_DIMENSIONS_PLUS_ONE_    :: diag    _NULL_
      REALTYPE,_ALLOCATABLE_ _ATTR_LOCATION_DIMENSIONS_HZ_PLUS_ONE_ :: diag_hz _NULL_
#endif

   end type type_environment

!-----------------------------------------------------------------------

   contains
   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes model information.
!
! !INTERFACE:
   subroutine init_model_info(modelinfo)
!
! !DESCRIPTION:
!  This function initializes the members of a model information derived type,
!  by setting them to a reasonable default value.
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETER:
      type (type_model_info),intent(inout) :: modelinfo
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      allocate(modelinfo%state_variables_ben(0))
      allocate(modelinfo%state_variables(0))
      allocate(modelinfo%diagnostic_variables_hz(0))
      allocate(modelinfo%diagnostic_variables(0))
      allocate(modelinfo%conserved_quantities(0))
      
      nullify(modelinfo%parent)
      nullify(modelinfo%firstchild)
      nullify(modelinfo%nextsibling)

      allocate(modelinfo%dependencies_hz(0))
      allocate(modelinfo%dependencies(0))
      
      modelinfo%nameprefix     = ''
      modelinfo%longnameprefix = ''
            
   end subroutine init_model_info
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes a state variable.
!
! !INTERFACE:
   subroutine init_state_variable_info(varinfo)
!
! !DESCRIPTION:
!  This function initializes the members of a state variable derived type,
!  by setting them to a reasonable default value.
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETER:
      type (type_state_variable_info), intent(inout) :: varinfo
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      varinfo%name = ''
      varinfo%units = ''
      varinfo%longname = ''
      varinfo%initial_value = _ZERO_
      varinfo%minimum = -1.e20
      varinfo%maximum = 1.e20
      varinfo%vertical_movement = _ZERO_
      varinfo%specific_light_extinction = _ZERO_
#if 0
      varinfo%mussels_inhale = .false.
#endif
      varinfo%no_precipitation_dilution = .false.
      varinfo%no_river_dilution         = .false.
      varinfo%dependencyid = id_not_used
      varinfo%globalid%id = id_not_used
   end subroutine init_state_variable_info
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes a diagnostic variable.
!
! !INTERFACE:
   subroutine init_diagnostic_variable_info(varinfo)
!
! !DESCRIPTION:
!  This function initializes the members of a diagnostic variable derived type,
!  by setting them to a reasonable default value.
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETER:
      type (type_diagnostic_variable_info), intent(inout) :: varinfo
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      varinfo%name = ''
      varinfo%units = ''
      varinfo%longname = ''
      varinfo%time_treatment = time_treatment_last
      varinfo%id = id_not_used
      varinfo%dependencyid = id_not_used
   end subroutine init_diagnostic_variable_info
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initializes a conserved quantity.
!
! !INTERFACE:
   subroutine init_conserved_quantity_info(conservedinfo)
!
! !DESCRIPTION:
!  This function initializes the members of a conserved quantity derived type,
!  by setting them to a reasonable default value.
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETER:
      type (type_conserved_quantity_info), intent(inout) :: conservedinfo
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      conservedinfo%name = ''
      conservedinfo%units = ''
      conservedinfo%longname = ''
      conservedinfo%id = id_not_used
   end subroutine init_conserved_quantity_info
!EOC
   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new state variable
!
! !INTERFACE:
   recursive function register_state_variable(modelinfo, name, units, longname, &
                                    initial_value, vertical_movement, specific_light_extinction, &
                                    mussels_inhale, minimum, maximum, &
                                    no_precipitation_dilution,no_river_dilution,benthic) &
                                    result(id)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical state variable in the global model database.
!  It returns an identifier that may be used later to retrieve the value of the state variable.
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETER:
      type (type_model_info),intent(inout)       :: modelinfo
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, longname, units
      REALTYPE,              intent(in),optional :: initial_value,vertical_movement,specific_light_extinction
      REALTYPE,              intent(in),optional :: minimum, maximum
      logical,               intent(in),optional :: mussels_inhale,no_precipitation_dilution,no_river_dilution
      logical,               intent(in),optional :: benthic
!
! !OUTPUT PARAMETER:
      type (type_state_variable_id)              :: id
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_state_variable_info),pointer :: variables_old(:),variables_new(:),curinfo
      character(len=256)                      :: text
      logical                                 :: benthic_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Determine whether this is a benthic variable (.false. by default)
      ! If so, select the corresponding array of state variables instead of the normal one.
      benthic_eff = .false.
      if (present(benthic)) benthic_eff = benthic
      if (benthic_eff) then
         variables_old => modelinfo%state_variables_ben
      else
         variables_old => modelinfo%state_variables
      end if

      ! Extend the state variable array and copy over old values.
      allocate(variables_new(ubound(variables_old,1)+1))
      variables_new(1:ubound(variables_old,1)) = variables_old(:)
      deallocate(variables_old)
      
      ! Assign new state variable array.
      if (benthic_eff) then
         modelinfo%state_variables_ben => variables_new
      else
         modelinfo%state_variables => variables_new
      end if
      
      curinfo => variables_new(ubound(variables_new,1))
      
      ! Initialize state variable info.
      call init_state_variable_info(curinfo)
      
      ! Store customized information on state variable.
      curinfo%name     = name
      curinfo%units    = units
      curinfo%longname = longname
      if (present(initial_value))             curinfo%initial_value = initial_value
      if (present(minimum))                   curinfo%minimum = minimum
      if (present(maximum))                   curinfo%maximum = maximum
      if (present(vertical_movement))         curinfo%vertical_movement = vertical_movement
      if (present(specific_light_extinction)) curinfo%specific_light_extinction = specific_light_extinction
#if 0
      if (present(mussels_inhale))            curinfo%mussels_inhale = mussels_inhale
#endif
      if (present(no_precipitation_dilution)) curinfo%no_precipitation_dilution = no_precipitation_dilution
      if (present(no_river_dilution        )) curinfo%no_river_dilution         = no_river_dilution

      ! Check for positive definiteness of initial value, if positive definiteness is specified.
      ! NB: although it would make sense to allow a value of exactly zero for positive definite
      ! state variables, this is not accepted by some integration schemes (notably: Patankar-related
      ! ones, causing NaNs). Therefore, initial values of zero are forbidden here.
      if (curinfo%initial_value<curinfo%minimum .or. curinfo%initial_value>curinfo%maximum) then
         write (text,*) 'Initial value',curinfo%initial_value,'for variable "'//trim(name)//'" lies&
               &outside allowed range',curinfo%minimum,'to',curinfo%maximum
         call fatal_error('rmbm_types::register_state_variable',text)
      end if
                  
      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the variable id.
      if (associated(modelinfo%parent)) then
         id = register_state_variable(modelinfo%parent,trim(modelinfo%nameprefix)//trim(curinfo%name), &
                trim(curinfo%units),trim(modelinfo%longnameprefix)//' '//trim(curinfo%longname),       &
                initial_value             = curinfo%initial_value,             &
                vertical_movement         = curinfo%vertical_movement,         &
                specific_light_extinction = curinfo%specific_light_extinction, &
                minimum                   = curinfo%minimum,                   &
                maximum                   = curinfo%maximum,                   &
                no_precipitation_dilution = curinfo%no_precipitation_dilution, &
                no_river_dilution         = curinfo%no_river_dilution,         &
                benthic                   = benthic_eff)
      else
         id%id = ubound(variables_new,1)
      end if
      
      ! Save the state variable's global id
      ! (index into state variable array of the root of the model tree).
      curinfo%globalid = id
      id%scale = _ONE_
      
   end function register_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   recursive function register_diagnostic_variable(modelinfo, name, units, longname, shape, time_treatment) result(id)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical diagnostic variable in the global model database.
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETER:
      type (type_model_info),intent(inout)       :: modelinfo
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, longname, units
      integer, optional,     intent(in)          :: time_treatment,shape
!
! !OUTPUT PARAMETER:
      integer                                    :: id
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_diagnostic_variable_info),pointer :: variables_old(:),variables_new(:),curinfo
      integer                                      :: shape_eff
!
!-----------------------------------------------------------------------
!BOC
      ! Determine whether this is a benthic state variable (.false. by default)
      shape_eff = shape_full
      if (present(shape)) shape_eff = shape
      select case (shape_eff)
         case (shape_hz)
            variables_old => modelinfo%diagnostic_variables_hz
         case (shape_full)
            variables_old => modelinfo%diagnostic_variables
         case default
            call fatal_error('rmbm_types::register_diagnostic_variable','invalid shape argument provided.')
      end select

      ! Extend the state variable array and copy over old values.
      allocate(variables_new(ubound(variables_old,1)+1))
      variables_new(1:ubound(variables_old,1)) = variables_old(:)
      deallocate(variables_old)
      
      ! Assign new state variable array.
      select case (shape_eff)
         case (shape_hz)
            modelinfo%diagnostic_variables_hz => variables_new
         case (shape_full)
            modelinfo%diagnostic_variables => variables_new
      end select
      
      curinfo => variables_new(ubound(variables_new,1))

      ! Initialize diagnostic variable info.
      call init_diagnostic_variable_info(curinfo)
      
      ! Store customized information on diagnostic variable.
      curinfo%name     = name
      curinfo%units    = units
      curinfo%longname = longname
      if (present(time_treatment)) curinfo%time_treatment = time_treatment
                  
      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the diagnostic variable id.
      if (associated(modelinfo%parent)) then
         id = register_diagnostic_variable(modelinfo%parent,trim(modelinfo%nameprefix)//trim(curinfo%name), &
                                           trim(curinfo%units),                                             &
                                           trim(modelinfo%longnameprefix)//' '//trim(curinfo%longname),     &
                                           time_treatment=curinfo%time_treatment,                           &
                                           shape = shape_eff)
      else
#ifdef _RMBM_MANAGE_DIAGNOSTICS_
         ! RMBM manages diagnostic variables - the identifier will be the number of the
         ! diagnostic variable in the global list maintained by the root of the model tree.
         id = ubound(variables_new,1)
#else
         ! The host manages diagnostic variables - the identifier will be the number of the
         ! diagnostic variable in the global list of dependencies, maintained by the root of the model tree.
         id = register_dependency(modelinfo,curinfo%name,shape)
#endif
      end if

#ifndef _RMBM_MANAGE_DIAGNOSTICS_
      curinfo%dependencyid = id
#endif
   end function register_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new conserved quantity
!
! !INTERFACE:
   recursive function register_conserved_quantity(modelinfo, name, units, longname) result(id)
!
! !DESCRIPTION:
!  This function registers a new biogeochemically conserved quantity in the global
!  model database.
!
! !USES:
   implicit none
!
! !INPUT/OUTPUT PARAMETER:
      type (type_model_info),intent(inout)       :: modelinfo
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, longname, units
!
! !OUTPUT PARAMETER:
      integer                                    :: id
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_conserved_quantity_info),pointer :: quantities_new(:)
!
!-----------------------------------------------------------------------
!BOC
      ! Extend array with conserved quantities, if needed
      allocate(quantities_new(ubound(modelinfo%conserved_quantities,1)+1))
      quantities_new(1:ubound(modelinfo%conserved_quantities,1)) = modelinfo%conserved_quantities(:)
      deallocate(modelinfo%conserved_quantities)
      modelinfo%conserved_quantities => quantities_new
      
      ! By default, the conserved quantity id is its index within the model.
      id = ubound(quantities_new,1)

      ! Initialize conserved quantity info.
      call init_conserved_quantity_info(quantities_new(id))
      
      ! Store customized information on conserved quantity.
      quantities_new(id)%name     = name
      quantities_new(id)%units    = units
      quantities_new(id)%longname = longname
                  
      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the conserved quantity id.
      if (associated(modelinfo%parent)) &
         id = register_conserved_quantity(modelinfo%parent,trim(modelinfo%nameprefix)//name, &
                 units,trim(modelinfo%longnameprefix)//' '//longname)
                
      ! Save the conserved quantity's global id.
      quantities_new(ubound(quantities_new,1))%globalid = id
      
   end function register_conserved_quantity
!EOC
   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on another biogeochemical state variable
!
! !INTERFACE:
   function register_state_dependency(modelinfo,name,benthic,mustexist) result(id)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
      type (type_model_info),target,          intent(in) :: modelinfo
      character(len=*),                       intent(in) :: name
      logical,optional,                       intent(in) :: benthic,mustexist
!
! !OUTPUT PARAMETER:
      type (type_state_variable_id)                      :: id
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      integer                                 :: i
      logical                                 :: benthic_eff,mustexist_eff
      type (type_model_info),pointer          :: curinfo
      type (type_state_variable_info),pointer :: variables(:)
!
!-----------------------------------------------------------------------
!BOC
      id%id = id_not_used
      
      ! If this model does not have a parent, there is no context to search variables in.
      if (.not. associated(modelinfo%parent)) return
      
      ! Determine whether this must be a benthos state varible (default: no).
      benthic_eff = .false.
      if (present(benthic)) benthic_eff = benthic

      ! Determine whether to throw an error if the variable is not found (default: yes).
      mustexist_eff = .true.
      if (present(mustexist)) mustexist_eff = mustexist

      ! Search every ancestor, starting at the parent,
      ! and moving one level up every time until the variable is found (or we reach the root).
      curinfo => modelinfo%parent
      do while (associated(curinfo))
         if (benthic_eff) then
            variables => curinfo%state_variables_ben
         else
            variables => curinfo%state_variables
         end if
         do i = 1,ubound(variables,1)
            if (variables(i)%name==name) then
               id = variables(i)%globalid
               return
            end if
         end do
         curinfo => curinfo%parent
      end do
      
      ! If we reaached this point, the variable was not found.
      ! Throw an error if the variable must exist.
      if (mustexist_eff) call fatal_error('rmbm_types::register_state_dependency', &
         'state variable dependency '//trim(name)//' of model '//trim(modelinfo%name)//' was not found.')
      
   end function register_state_dependency
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a read-only dependency on another variable.
!
! !INTERFACE:
   recursive function register_dependency(modelinfo,name,shape) result(id)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !USES:
   implicit none
!
! !INPUT PARAMETERS:
      type (type_model_info),target,          intent(inout) :: modelinfo
      character(len=*),                       intent(in)    :: name
      integer,optional,                       intent(in)    :: shape
!
! !OUTPUT PARAMETER:
      integer                           :: id
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      integer                                      :: n,i,realshape
      type (type_model_info),pointer               :: proot
      character(len=64),pointer                    :: dependencies_new(:)
      character(len=64),pointer                    :: source(:)
      type (type_state_variable_info),     pointer :: statevarinfo(:)
      type (type_diagnostic_variable_info),pointer :: diagvarinfo(:)
!
!-----------------------------------------------------------------------
!BOC
      ! Get effective shape argument - defaults to the full domain (instead of a horizontal slice).
      realshape = shape_full
      if (present(shape)) realshape = shape

      ! Find the root of the model tree
      proot => modelinfo
      do while (associated(proot%parent))
         proot => proot%parent
      end do
      
      ! Get pointer to array with dependencies (use shape argument to determine which)
      select case (realshape)
         case (shape_hz)
            source => proot%dependencies_hz
            statevarinfo => proot%state_variables_ben
            diagvarinfo => proot%diagnostic_variables_hz
         case (shape_full)
            source => proot%dependencies
            statevarinfo => proot%state_variables
            diagvarinfo => proot%diagnostic_variables
         case default
            call fatal_error('rmbm_types::register_dependency','Invalid shape argument given')
      end select

      n = ubound(source,1)

      ! Search existing dependencies and return the corresponding id if found.
      do id=1,n
         if (source(id).eq.name) return
      end do
      
      ! Dependency was not registered yet - create extended array to hold new dependency.
      allocate(dependencies_new(n+1))
      dependencies_new(1:n) = source(:)
      deallocate(source)
      dependencies_new(n+1) = name
      
      ! Assign new array with dependencies (variable to assign to depends on shape argument)
      select case (realshape)
         case (shape_hz)
            proot%dependencies_hz => dependencies_new
         case (shape_full)
            proot%dependencies => dependencies_new
      end select

      ! Variable identifier equal the previous number of dependencies plus 1.
      id = n+1

      ! Determine whether this new dependency matches a registered state variable.
      ! (in that case RMBM can provide the value internally)
      do i=1,ubound(statevarinfo,1)
         if (statevarinfo(i)%name==name) then
            statevarinfo(i)%dependencyid = id
         end if
      end do

#ifdef _RMBM_MANAGE_DIAGNOSTICS_
      ! Determine whether this new dependency matches a registered diagnostic variable.
      ! (in that case RMBM can provide the value internally)
      do i=1,ubound(diagvarinfo,1)
         if (diagvarinfo(i)%name==name) then
            diagvarinfo(i)%dependencyid = id
         end if
      end do
#endif
      
   end function register_dependency
!EOC

!-----------------------------------------------------------------------
   
   end module rmbm_types

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

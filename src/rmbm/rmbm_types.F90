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
   public type_environment,type_state_2d,type_state
!
! !PUBLIC DERIVED TYPES:
!

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
      
      integer  :: globalid                   ! This is a globally unique identifier for the variable that can be used to retrieve values.
      integer  :: id
   end type type_state_variable_info

   ! Properties of a diagnostic variable
   type type_diagnostic_variable_info
      character(len=64) :: name, longname, units
      integer           :: id,globalid
      
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
      type (type_state_variable_info),     pointer,dimension(:) :: state_variables_2d,state_variables_3d
      type (type_diagnostic_variable_info),pointer,dimension(:) :: diagnostic_variables_2d,diagnostic_variables_3d
      type (type_conserved_quantity_info), pointer,dimension(:) :: conserved_quantities
      
      type (type_model_info),pointer :: parent
      type (type_model_info),pointer :: firstchild
      type (type_model_info),pointer :: nextsibling
      
      character(len=64) :: nameprefix,longnameprefix
      
      character(len=64),pointer :: dependencies3d(:)
      character(len=64),pointer :: dependencies2d(:)
   end type type_model_info
   
   ! Parameters
   integer, parameter, public         :: shape2d=2,shape3d=3
   
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
      REALTYPE,pointer ATTR_LOCATION_DIMENSIONS :: data
   end type type_state

   type type_state_2d
      REALTYPE,pointer ATTR_LOCATION_DIMENSIONS_HZ :: data
   end type type_state_2d

   type type_environment
      type (type_state   ), dimension(:), _ALLOCATABLE :: state3d _NULL ! array of pointers to data of pelagic state variables
      type (type_state   ), dimension(:), _ALLOCATABLE :: var3d   _NULL ! array of pointers to data of all pelagic variables (state and diagnostic)
      type (type_state_2d), dimension(:), _ALLOCATABLE :: state2d _NULL ! array of pointers to data of benthic state variables
      type (type_state_2d), dimension(:), _ALLOCATABLE :: var2d   _NULL ! array of pointers to data of all horizontal variables (state and diagnostic, surface and bottom)
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
      allocate(modelinfo%state_variables_2d(0))
      allocate(modelinfo%state_variables_3d(0))
      allocate(modelinfo%diagnostic_variables_2d(0))
      allocate(modelinfo%diagnostic_variables_3d(0))
      allocate(modelinfo%conserved_quantities(0))
      
      nullify(modelinfo%parent)
      nullify(modelinfo%firstchild)
      nullify(modelinfo%nextsibling)

      allocate(modelinfo%dependencies2d(0))
      allocate(modelinfo%dependencies3d(0))
      
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
      varinfo%id = id_not_used
      varinfo%time_treatment = time_treatment_last
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
      conservedinfo%id = -1
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
      integer                                    :: id
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
      ! Determine whether this is a benthic state variable (.false. by default)
      ! If so, select the array of 2d state variables instead of the normal array with 3d variables.
      benthic_eff = .false.
      if (present(benthic)) benthic_eff = benthic
      if (benthic_eff) then
         variables_old => modelinfo%state_variables_2d
      else
         variables_old => modelinfo%state_variables_3d
      end if

      ! Extend the state variable array and copy over old values.
      allocate(variables_new(ubound(variables_old,1)+1))
      variables_new(1:ubound(variables_old,1)) = variables_old(:)
      deallocate(variables_old)
      
      ! Assign new state variable array.
      if (benthic_eff) then
         modelinfo%state_variables_2d => variables_new
      else
         modelinfo%state_variables_3d => variables_new
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
         id = ubound(variables_new,1)
      end if
      
      ! Save the state variable's global id
      ! (index into state variable array of the root of the model tree).
      curinfo%globalid = id
      
   end function register_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   recursive function register_diagnostic_variable(modelinfo, name, units, longname, benthic, time_treatment) result(id)
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
      integer, optional,     intent(in)          :: time_treatment
      logical, optional,     intent(in)          :: benthic
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
      logical                                      :: benthic_eff
      integer                                      :: shape
!
!-----------------------------------------------------------------------
!BOC
      ! Determine whether this is a benthic state variable (.false. by default)
      benthic_eff = .false.
      if (present(benthic)) benthic_eff = benthic
      if (benthic_eff) then
         variables_old => modelinfo%diagnostic_variables_2d
      else
         variables_old => modelinfo%diagnostic_variables_3d
      end if

      ! Extend the state variable array and copy over old values.
      allocate(variables_new(ubound(variables_old,1)+1))
      variables_new(1:ubound(variables_old,1)) = variables_old(:)
      deallocate(variables_old)
      
      ! Assign new state variable array.
      if (benthic_eff) then
         modelinfo%diagnostic_variables_2d => variables_new
      else
         modelinfo%diagnostic_variables_3d => variables_new
      end if
      
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
                                           benthic = benthic_eff)
      else
         shape = shape3d
         if (benthic_eff) shape = shape2d
         id = register_dependency(modelinfo,curinfo%name,shape)
      end if
                
      ! Save the diagnostic variable's global id.
      ! (index into dependencies array of the root of the tree)
      curinfo%globalid = id
      
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
   recursive function register_state_dependency(modelinfo,name,benthic) result(id)
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
      logical,optional,                       intent(in) :: benthic
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
      integer                                 :: i
      logical                                 :: benthic_eff
      type (type_model_info),pointer          :: curinfo
      type (type_state_variable_info),pointer :: variables(:)
!
!-----------------------------------------------------------------------
!BOC
      id = -1
      
      ! If this model does not have a parent, there is no context to search variables in.
      if (.not. associated(modelinfo%parent)) return
      
      benthic_eff = .false.
      if (present(benthic)) benthic_eff = benthic

      ! Search amongst ancestors.
      curinfo => modelinfo%parent
      do while (associated(curinfo))
         if (benthic_eff) then
            variables => curinfo%state_variables_2d
         else
            variables => curinfo%state_variables_3d
         end if
         do i = 1,ubound(variables,1)
            if (variables(i)%name==name) then
               id = variables(i)%globalid
               return
            end if
         end do
         curinfo => curinfo%parent
      end do
      
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
      integer                                 :: n,realshape
      type (type_model_info),pointer          :: proot
      character(len=64),pointer               :: dependencies_new(:)
      character(len=64),pointer               :: source(:)
!
!-----------------------------------------------------------------------
!BOC
      ! Get effective shape argument - defaults to 3D.
      realshape = shape3d
      if (present(shape)) realshape = shape

      ! Find the root of the model tree
      proot => modelinfo
      do while (associated(proot%parent))
         proot => proot%parent
      end do
      
      ! Get pointer to array with dependencies (use shape argument to determine which)
      select case (realshape)
         case (shape2d)
            source => proot%dependencies2d
         case (shape3d)
            source => proot%dependencies3d
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
         case (shape2d)
            proot%dependencies2d => dependencies_new
         case (shape3d)
            proot%dependencies3d => dependencies_new
      end select
      
      ! Variable identifier equal the previous number of dependencies plus 1.
      id = n+1
      
   end function register_dependency
!EOC

!-----------------------------------------------------------------------
   
   end module rmbm_types

!-----------------------------------------------------------------------
! Copyright under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

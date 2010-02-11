!$Id: bio_types.F90,v 1.6 2009-05-11 13:41:41 jorn Exp $
#include"cppdefs.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_types --- Derived types used by 0D biogeochemical modules
!
! !INTERFACE:
   module bio_types
   
   implicit none
!
! !DESCRIPTION:
! This module contains the derived types that are used for communication between
! the 0D biogeochemical models and a hosting physical environment (e.g. GOTM).
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
          register_state_variable_dependency
!
! !PUBLIC DERIVED TYPES:
!

   ! Properties of a single state variable
   type type_state_variable_info
      character(len=64) :: name, longname, unit

      REALTYPE :: initial_value              ! Initial state variable value
      REALTYPE :: vertical_movement          ! Vertical movement (m/s) due to e.g. sinking, floating or activity. Note: positive for upward movement!
      REALTYPE :: specific_light_extinction  ! Specific light extinction (/m/state variable unit)
#if 0
      REALTYPE :: mussels_inhale             ! Whether this variable can be consumed by mussels
#endif
      logical  :: positive_definite          ! Whether this variable is positive definite (negative and zero values are invalid)
      
      integer  :: globalid                   ! This is a globally unique identifier for the variable that can be used in getbiovar
   end type type_state_variable_info

   ! Properties of a diagnostic variable
   type type_diagnostic_variable_info
      character(len=64) :: name, longname, unit
      integer           :: id,globalid
      
      ! Time treatment:
      ! 0: last value
      ! 1: time-integrated
      ! 2: time step-averaged
      ! 3: time step-integrated
      integer           :: time_treatment
   end type type_diagnostic_variable_info

   ! Properties of a conserved quantity
   type type_conserved_quantity_info
      character(len=64) :: name, longname, unit
      integer           :: globalid,id
   end type type_conserved_quantity_info

   ! Global 0D model properties
   type type_model_info
      ! Number of state variables
      integer :: state_variable_count, diagnostic_variable_count, conserved_quantity_count

      ! Vertical movement type
      ! 0: vertical movement of all variables is constant in time and space
      ! 2: vertical movement of one or more variables depends on time and space
      !    (optionally including the current state and local enviroment)
      ! NB: other values (e.g., 1) are reserved for semi-variable vertical movement:
      ! e.g., in depth only, but not in time. Currently this is not implemented and
      ! therefore requires dynamic_vertical_movement=2.
      integer  :: dynamic_vertical_movement
      
      type (type_state_variable_info),     pointer :: variables(:)
      type (type_diagnostic_variable_info),pointer :: diagnostic_variables(:)
      type (type_conserved_quantity_info), pointer :: conserved_quantities(:)
      
      type (type_model_info),pointer :: master
      
      character(len=64) :: nameprefix,longnameprefix
   end type type_model_info

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
   IMPLICIT NONE
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
      modelinfo%state_variable_count      = 0
      modelinfo%diagnostic_variable_count = 0
      modelinfo%conserved_quantity_count  = 0
      
      modelinfo%dynamic_vertical_movement = 0

      modelinfo%variables => null()
      modelinfo%diagnostic_variables => null()
      modelinfo%conserved_quantities => null()
      
      modelinfo%master => null()
      
      modelinfo%nameprefix = ''
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
   IMPLICIT NONE
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
      varinfo%unit = ''
      varinfo%longname = ''
      varinfo%initial_value = _ZERO_
      varinfo%vertical_movement = _ZERO_
      varinfo%specific_light_extinction = _ZERO_
      varinfo%positive_definite = .false.
#if 0
      varinfo%mussels_inhale = .false.
#endif
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
   IMPLICIT NONE
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
      varinfo%unit = ''
      varinfo%longname = ''
      varinfo%id = -1
      varinfo%time_treatment = 0
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
   IMPLICIT NONE
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
      conservedinfo%unit = ''
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
   recursive function register_state_variable(modelinfo, name, unit, longname, &
                                    initial_value, vertical_movement, specific_light_extinction, &
                                    positive_definite, mussels_inhale) result(id)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical state variable in the global model database.
!  It returns an identifier that may be used later to retrieve the value of the state variable.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETER:
      type (type_model_info),intent(inout)       :: modelinfo
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, longname, unit
      REALTYPE,              intent(in),optional :: initial_value,vertical_movement,specific_light_extinction
      logical,               intent(in),optional :: positive_definite,mussels_inhale
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
      type (type_state_variable_info),pointer :: variables_new(:)
!
!-----------------------------------------------------------------------
!BOC
      ! Check for positive definiteness of initial value, if positive definiteness is specified.
      ! NB: although it would make sense to allow a value of exactly zero for positive definite
      ! state variables, this is not accepted by some integration schemes (notably: Patankar-related
      ! ones, causing NaNs). Therefore, initial values of zero are forbidden here.
      if (present(positive_definite)) then
         if (positive_definite .and. initial_value.le._ZERO_) then
            FATAL 'Variable "'//trim(name)//'" is specified to be positive definite, and must therefore be&
                 & initialized to a positive, non-zero value. However, the initial value is currently set to', &
                 & initial_value
            stop 'bio_types::register_state_variable'
         end if
      end if

      ! Extend the state variable array.
      allocate(variables_new(modelinfo%state_variable_count+1))
      if (associated(modelinfo%variables)) then
         variables_new(1:modelinfo%state_variable_count) = modelinfo%variables(:)
         deallocate(modelinfo%variables)
      end if
      modelinfo%variables => variables_new
      modelinfo%state_variable_count = modelinfo%state_variable_count+1
      
      ! By default, the variable id is the index of the state variable.
      id = modelinfo%state_variable_count

      ! Initialize state variable info.
      call init_state_variable_info(modelinfo%variables(id))
      
      ! Store customized information on state variable.
      modelinfo%variables(id)%name = name
      modelinfo%variables(id)%unit = unit
      modelinfo%variables(id)%longname = longname
      if (present(initial_value))             modelinfo%variables(id)%initial_value = initial_value
      if (present(vertical_movement))         modelinfo%variables(id)%vertical_movement = vertical_movement
      if (present(specific_light_extinction)) modelinfo%variables(id)%specific_light_extinction = specific_light_extinction
      if (present(positive_definite))         modelinfo%variables(id)%positive_definite = positive_definite
#if 0
      if (present(mussels_inhale))            modelinfo%variables(id)%mussels_inhale = mussels_inhale
#endif
                  
      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the variable id.
      if (associated(modelinfo%master)) &
         id = register_state_variable(modelinfo%master,trim(modelinfo%nameprefix)//name, &
                unit,trim(modelinfo%longnameprefix)//' '//longname, &
                initial_value             = modelinfo%variables(id)%initial_value, &
                vertical_movement         = modelinfo%variables(id)%vertical_movement, &
                specific_light_extinction = modelinfo%variables(id)%specific_light_extinction, &
                positive_definite         = modelinfo%variables(id)%positive_definite)
                
      ! Save the variable's global id.
      modelinfo%variables(modelinfo%state_variable_count)%globalid = id
      
   end function register_state_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new diagnostic variable
!
! !INTERFACE:
   recursive function register_diagnostic_variable(modelinfo, name, unit, longname, time_treatment) result(id)
!
! !DESCRIPTION:
!  This function registers a new biogeochemical diagnostic variable in the global model database.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETER:
      type (type_model_info),intent(inout)       :: modelinfo
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, longname, unit
      integer, optional,     intent(in)          :: time_treatment
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
      type (type_diagnostic_variable_info),pointer :: variables_new(:)
!
!-----------------------------------------------------------------------
!BOC
      ! Extend array with diagnostic variables, if needed
      allocate(variables_new(modelinfo%diagnostic_variable_count+1))
      if (associated(modelinfo%diagnostic_variables)) then
         variables_new(1:modelinfo%diagnostic_variable_count) = modelinfo%diagnostic_variables(:)
         deallocate(modelinfo%diagnostic_variables)
      end if
      modelinfo%diagnostic_variables => variables_new
      modelinfo%diagnostic_variable_count = modelinfo%diagnostic_variable_count+1
      
      ! By default, the variable id is the index of the diagnostic variable.
      id = modelinfo%diagnostic_variable_count

      ! Initialize diagnostic variable info.
      call init_diagnostic_variable_info(modelinfo%diagnostic_variables(id))
      
      ! Store customized information on diagnostic variable.
      modelinfo%diagnostic_variables(id)%name = name
      modelinfo%diagnostic_variables(id)%unit = unit
      modelinfo%diagnostic_variables(id)%longname = longname
      if (present(time_treatment)) modelinfo%diagnostic_variables(id)%time_treatment = time_treatment
                  
      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the diagnostic variable id.
      if (associated(modelinfo%master)) &
         id = register_diagnostic_variable(modelinfo%master,trim(modelinfo%nameprefix)//name, &
                 unit,trim(modelinfo%longnameprefix)//' '//longname, &
                 modelinfo%diagnostic_variables(id)%time_treatment)
                
      ! Save the diagnostic variable's global id.
      modelinfo%diagnostic_variables(modelinfo%diagnostic_variable_count)%globalid = id
      
   end function register_diagnostic_variable
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a new conserved quantity
!
! !INTERFACE:
   recursive function register_conserved_quantity(modelinfo, name, unit, longname) result(id)
!
! !DESCRIPTION:
!  This function registers a new biogeochemically conserved quantity in the global
!  model database.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETER:
      type (type_model_info),intent(inout)       :: modelinfo
!
! !INPUT PARAMETERS:
      character(len=*),      intent(in)          :: name, longname, unit
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
      allocate(quantities_new(modelinfo%conserved_quantity_count+1))
      if (associated(modelinfo%conserved_quantities)) then
         quantities_new(1:modelinfo%conserved_quantity_count) = modelinfo%conserved_quantities(:)
         deallocate(modelinfo%conserved_quantities)
      end if
      modelinfo%conserved_quantities => quantities_new
      modelinfo%conserved_quantity_count = modelinfo%conserved_quantity_count+1
      
      ! By default, the conserved quantity id is its index within the model.
      id = modelinfo%conserved_quantity_count

      ! Initialize conserved quantity info.
      call init_conserved_quantity_info(modelinfo%conserved_quantities(id))
      
      ! Store customized information on conserved quantity.
      modelinfo%conserved_quantities(id)%name = name
      modelinfo%conserved_quantities(id)%unit = unit
      modelinfo%conserved_quantities(id)%longname = longname
                  
      ! If this model runs as part of a larger collection,
      ! the collection (the "master") determines the conserved quantity id.
      if (associated(modelinfo%master)) &
         id = register_conserved_quantity(modelinfo%master,trim(modelinfo%nameprefix)//name, &
                 unit,trim(modelinfo%longnameprefix)//' '//longname)
                
      ! Save the conserved quantity's global id.
      modelinfo%conserved_quantities(modelinfo%conserved_quantity_count)%globalid = id
      
   end function register_conserved_quantity
!EOC
   
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on another biogeochemical state variable
!
! !INTERFACE:
   recursive function register_state_variable_dependency(modelinfo,name) result(id)
!
! !DESCRIPTION:
!  This function searches for a biogeochemical state variable by the user-supplied name
!  in the global model database. It returns the identifier of the variable (or -1 if
!  the variable is not found), which may be used to retrieve the variable value at a later stage.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
      type (type_model_info),intent(in) :: modelinfo
      character(len=*),      intent(in) :: name
!
! !OUTPUT PARAMETER:
      integer                           :: id
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
!-----------------------------------------------------------------------
!BOC
      if (associated(modelinfo%master)) then
         ! There is a master model attached - ask that model for the variable id instead.
         id = register_state_variable_dependency(modelinfo%master,name)
      else
         ! No master model attached - look up the variable in our ptivate list of state variables.
         do id = 1,modelinfo%state_variable_count
            if (modelinfo%variables(id)%name==name) exit
         end do
         if (id>modelinfo%state_variable_count) then
            FATAL 'Could not locate variable "'//trim(name)//'".'
            stop 'bio_types::register_state_variable_dependency'
            id = -1
         end if
      end if
      
   end function register_state_variable_dependency
!EOC

!-----------------------------------------------------------------------
   
   end module bio_types

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

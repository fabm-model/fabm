! ====================================================================================================
! This module declares objects ("identities") for standard physical-biogeochemical variables
! that have a well-defined interpretation and unit.
!
! To use these identities, "use" the fabm_types module, which provides a single "standard_variables"
! variable (derived type) that has all standard identities as its members. This can then be used as
! standard_variables%temperature, standard_variables%wind_speed, etc.
! For a list of all supported variables, please see:
! http://sourceforge.net/apps/mediawiki/fabm/index.php?title=List_of_standard_variables
!
! Biogeochemical models can use these "identity" objects in two ways. First, they can access the value
! of the corresponding variable by registering them as dependency. To do so, call register_dependency
! with the "identity" object (e.g., standard_variables%temperature) as argument.
! Additionally, biogeochemical models can assign a standard identities to their own variables during
! registration. FABM will couple all variables that have been assigned the same identity. Thus, using
! standard identities enables implicit variable coupling.
! ====================================================================================================
! The names of standard variables are based on the Standard Name Table from the NetCDF Climate and
! Forecast (CF) Metadata Convention. See http://cf-pcmdi.llnl.gov/documents/cf-standard-names/.
! In deriving names from the CF convention, the following exceptions are made to account for the fact
! that FABM handles both marine and limnic systems and has the water column as default domain:
! - "sea_water_" prefix is suppressed.
! - "_in_sea_water" suffix is suppressed.
! - instead of the "_at_sea_floor" suffix a "bottom_" prefix is used, analogous to the "surface_"
!   prefix used in CF.
! - the "sea_floor_" prefix is replaced by a "bottom_" prefix.
! ====================================================================================================

module fabm_standard_variables

   private

   public type_standard_variable, type_bulk_standard_variable, type_horizontal_standard_variable, type_global_standard_variable
   public standard_variables, initialize_standard_variables

   ! ====================================================================================================
   ! Data types that contain all metadata needed to describe standard variables.
   ! ====================================================================================================

   type type_standard_variable
      character(len=256) :: name  = ''    ! Name
      character(len=64)  :: units = ''    ! Units
      character(len=512) :: cf_names = '' ! Comma-separated list of standard names defined in the NetCDF CF convention
                                          ! (http://cf-pcmdi.llnl.gov/documents/cf-standard-names/)
      logical            :: aggregate_variable = .false. ! Whether biogeochemical models can contribute (add to) this variable.
                                                         ! If .true., this variable is always available with a default value of 0.
   contains
      procedure :: is_null => standard_variable_is_null
   end type

   type,extends(type_standard_variable) :: type_bulk_standard_variable
   contains
      procedure :: compare => bulk_standard_variable_compare
   end type

   type,extends(type_standard_variable) :: type_horizontal_standard_variable
   contains
      procedure :: compare => horizontal_standard_variable_compare
   end type

   type,extends(type_standard_variable) :: type_global_standard_variable
   contains
      procedure :: compare => global_standard_variable_compare
   end type

   type type_standard_variable_collection
#include "standard_variables.h"
   end type

   ! Single instance of the collection that contains all standard variables.
   type (type_standard_variable_collection),save :: standard_variables

contains

   subroutine initialize_standard_variables()
#include "standard_variable_assignments.h"
   end subroutine

   logical function standard_variable_is_null(variable)
      class (type_standard_variable),intent(in) :: variable
      standard_variable_is_null = (variable%name==''.and. variable%units=='')
   end function

   logical function standard_variable_compare(variable1,variable2)
      class (type_standard_variable),intent(in) :: variable1,variable2
      standard_variable_compare = (variable1%name ==''.or.variable2%name ==''.or.variable1%name ==variable2%name ) &
                            .and. (variable1%units==''.or.variable2%units==''.or.variable1%units==variable2%units)
   end function

   logical function bulk_standard_variable_compare(variable1,variable2)
      class (type_bulk_standard_variable),intent(in) :: variable1,variable2
      bulk_standard_variable_compare = standard_variable_compare(variable1,variable2)
   end function

   logical function horizontal_standard_variable_compare(variable1,variable2)
      class (type_horizontal_standard_variable),intent(in) :: variable1,variable2
      horizontal_standard_variable_compare = standard_variable_compare(variable1,variable2)
   end function

   logical function global_standard_variable_compare(variable1,variable2)
      class (type_global_standard_variable),intent(in) :: variable1,variable2
      global_standard_variable_compare = standard_variable_compare(variable1,variable2)
   end function

end module fabm_standard_variables
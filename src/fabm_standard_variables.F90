! ====================================================================================================
! This module declares objects ("identities") for standard physical-biogeochemical variables
! that have a well-defined interpretation and unit.
!
! To use these identities, "use" the fabm_types module, which declares a single "standard_variables"
! variable that has all standard identities as its members.
! (for a list, see the definion of type_standard_variable_collection below)
!
! Biogeochemical models can use these "identity" objects in two ways. First, they can access the value
! of the corresponding variable by registering them as dependency. To do so, call register_dependency
! with the "identity" object as argument.
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
   public type_standard_variable_collection, standard_variables
   
   ! ====================================================================================================
   ! Data types that contain all metadata needed to describe standard variables.
   ! ====================================================================================================

   type type_standard_variable
      character(len=256) :: name  = ''    ! Name
      character(len=64)  :: units = ''    ! Units
      character(len=512) :: cf_names = '' ! Comma-separated list of standard names defined in the NetCDF CF convention
                                          ! (http://cf-pcmdi.llnl.gov/documents/cf-standard-names/)
   end type
   
   type,extends(type_standard_variable) :: type_bulk_standard_variable
   end type
   
   type,extends(type_standard_variable) :: type_horizontal_standard_variable
   end type
   
   type,extends(type_standard_variable) :: type_global_standard_variable
   end type

#include "standard_variables.h"

   end module fabm_standard_variables
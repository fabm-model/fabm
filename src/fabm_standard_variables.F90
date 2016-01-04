! ====================================================================================================
! This module declares objects ("identities") for standard physical-biogeochemical variables
! that have a well-defined interpretation and unit.
!
! To use these identities, "use" the fabm_types module, which provides a single "standard_variables"
! variable (derived type) that has all standard identities as its members. This can then be used as
! standard_variables%temperature, standard_variables%wind_speed, etc.
! For a list of all supported variables, please see:
! http://sourceforge.net/p/fabm/wiki/List_of_standard_variables
!
! Biogeochemical models can use these "identity" objects in two ways. First, they can access the value
! of the corresponding variable by registering them as dependency. To do so, call register_dependency
! with the "identity" object (e.g., standard_variables%temperature) as argument.
! Additionally, biogeochemical models can assign a standard identities to their own variables during
! registration. FABM will couple all variables that have been assigned the same identity. Thus, using
! standard identities enables implicit variable coupling.
! ====================================================================================================
! The names of standard variables are based on the Standard Name Table from the NetCDF Climate and
! Forecast (CF) Metadata Convention. See http://cfconventions.org/.
! In deriving names from the CF convention, the following exceptions are made to account for the fact
! that FABM handles both marine and limnic systems and has the water column as default domain:
! - "sea_water_" prefix is suppressed.
! - "_in_sea_water" suffix is suppressed.
! - instead of the "_at_sea_floor" suffix a "bottom_" prefix is used, analogous to the "surface_"
!   prefix used in CF.
! - the "sea_floor_" prefix is replaced by a "bottom_" prefix.
! ====================================================================================================

module fabm_standard_variables

   implicit none

   private

   public type_base_standard_variable, type_bulk_standard_variable, type_horizontal_standard_variable, type_global_standard_variable
   public type_standard_variable_node, type_standard_variable_set
   public standard_variables, initialize_standard_variables

   ! ====================================================================================================
   ! Data types that contain all metadata needed to describe standard variables.
   ! ====================================================================================================

   type type_base_standard_variable
      character(len=256) :: name  = ''    ! Name
      character(len=64)  :: units = ''    ! Units
      character(len=512) :: cf_names = '' ! Comma-separated list of standard names defined in the NetCDF CF convention
                                          ! (http://cf-pcmdi.llnl.gov/documents/cf-standard-names/)
      logical            :: aggregate_variable = .false. ! Whether biogeochemical models can contribute (add to) this variable.
                                                         ! If .true., this variable is always available with a default value of 0.
      logical            :: conserved = .false.          ! Whether this variable shoudl be included in lists of conserved quantities.
   contains
      procedure :: is_null => standard_variable_is_null
      procedure :: compare => standard_variable_compare
   end type

   type,extends(type_base_standard_variable) :: type_bulk_standard_variable
   end type

   type,extends(type_base_standard_variable) :: type_horizontal_standard_variable
   end type

   type,extends(type_base_standard_variable) :: type_global_standard_variable
   end type

   type type_standard_variable_node
      class (type_base_standard_variable),pointer :: p    => null()
      type (type_standard_variable_node), pointer :: next => null()
   end type

   type type_standard_variable_set
      type (type_standard_variable_node), pointer :: first => null()
   contains
      procedure :: contains_variable => standard_variable_set_contains_variable
      procedure :: contains_name     => standard_variable_set_contains_name
      generic   :: contains => contains_variable,contains_name
      procedure :: add      => standard_variable_set_add
      procedure :: update   => standard_variable_set_update
      procedure :: finalize => standard_variable_set_finalize
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
      class (type_base_standard_variable),intent(in) :: variable
      standard_variable_is_null = (variable%name==''.and. variable%units=='')
   end function

   logical function standard_variable_compare(variable1,variable2)
      class (type_base_standard_variable),intent(in) :: variable1,variable2
      standard_variable_compare = .false.

      ! First test whether the types match.
      select type (variable1)
         class is (type_bulk_standard_variable)
            select type (variable2)
               class is (type_bulk_standard_variable)
                  standard_variable_compare = .true.
            end select
         class is (type_horizontal_standard_variable)
            select type (variable2)
               class is (type_horizontal_standard_variable)
                  standard_variable_compare = .true.
            end select
         class is (type_global_standard_variable)
            select type (variable2)
               class is (type_global_standard_variable)
                  standard_variable_compare = .true.
            end select
      end select

      ! If types do not match, the standard variables are not equal - we're done.
      if (.not.standard_variable_compare) return

      ! Compare the metadata of the standard variables.
      standard_variable_compare = (variable1%name ==''.or.variable2%name ==''.or.variable1%name ==variable2%name ) &
                            .and. (variable1%units==''.or.variable2%units==''.or.variable1%units==variable2%units)
   end function standard_variable_compare

   logical function standard_variable_set_contains_variable(self,standard_variable)
      class (type_standard_variable_set), intent(in) :: self
      class (type_base_standard_variable),intent(in) :: standard_variable

      type (type_standard_variable_node), pointer :: node

      standard_variable_set_contains_variable = .true.
      node => self%first
      do while (associated(node))
         if (standard_variable%compare(node%p)) return
         node => node%next
      end do
      standard_variable_set_contains_variable = .false.
   end function standard_variable_set_contains_variable

   logical function standard_variable_set_contains_name(self,name)
      class (type_standard_variable_set),intent(in) :: self
      character(len=*),                  intent(in) :: name

      type (type_standard_variable_node), pointer :: node

      standard_variable_set_contains_name = .true.
      node => self%first
      do while (associated(node))
         if (node%p%name==name) return
         node => node%next
      end do
      standard_variable_set_contains_name = .false.
   end function standard_variable_set_contains_name

   subroutine standard_variable_set_add(self,standard_variable)
      class (type_standard_variable_set), intent(inout) :: self
      class (type_base_standard_variable),intent(in)    :: standard_variable

      type (type_standard_variable_node), pointer :: node

      if (self%contains(standard_variable)) return

      if (.not.associated(self%first)) then
         allocate(self%first)
         node => self%first
      else
         node => self%first
         do while (associated(node%next))
            node => node%next
         end do
         allocate(node%next)
         node => node%next
      end if
      allocate(node%p,source=standard_variable)
   end subroutine standard_variable_set_add

   subroutine standard_variable_set_update(self,other)
      class (type_standard_variable_set),intent(inout) :: self
      class (type_standard_variable_set),intent(in)    :: other

      type (type_standard_variable_node), pointer :: node

      node => other%first
      do while (associated(node))
         call self%add(node%p)
         node => node%next
      end do
   end subroutine standard_variable_set_update

   subroutine standard_variable_set_finalize(self)
      class (type_standard_variable_set),intent(inout) :: self

      type (type_standard_variable_node), pointer :: node,next_node

      node => self%first
      do while (associated(node))
         next_node => node%next
         deallocate(node%p)
         deallocate(node)
         node => next_node
      end do
      self%first => null()
   end subroutine standard_variable_set_finalize

end module fabm_standard_variables
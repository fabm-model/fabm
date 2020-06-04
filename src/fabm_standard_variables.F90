! ====================================================================================================
! This module declares objects ("identities") for standard physical-biogeochemical variables
! that have a well-defined interpretation and unit.
!
! To use these identities, "use" the fabm_types module, which provides a single "standard_variables"
! variable (derived type) that has all standard identities as its members. This can then be used as
! standard_variables%temperature, standard_variables%wind_speed, etc.
! For a list of all supported variables, please see:
! https://github.com/fabm-model/fabm/wiki/List-of-standard-variables
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

   public type_base_standard_variable
   public type_universal_standard_variable, type_domain_specific_standard_variable, type_interior_standard_variable, type_horizontal_standard_variable, type_surface_standard_variable, type_bottom_standard_variable, type_global_standard_variable
   public type_standard_variable_node, type_standard_variable_set
   public standard_variables, initialize_standard_variables

   ! ====================================================================================================
   ! Data types that contain all metadata needed to describe standard variables.
   ! ====================================================================================================

   type type_base_standard_variable
      character(len=256) :: name  = ''    ! Name
      character(len=64)  :: units = ''    ! Units
      character(len=512) :: cf_names = '' ! Comma-separated list of standard names defined in the NetCDF CF convention
                                          ! (http://cfconventions.org/standard-names.html)
      logical            :: aggregate_variable = .false. ! Whether biogeochemical models can contribute (add to) this variable.
                                                         ! If .true., this variable is always available with a default value of 0.
      logical, private   :: resolved = .false.
   contains
      procedure :: resolve         => base_standard_variable_resolve
      procedure :: assert_resolved => base_standard_variable_assert_resolved
   end type

   type, extends(type_base_standard_variable) :: type_domain_specific_standard_variable
      type (type_universal_standard_variable), pointer :: universal => null()
   contains
      procedure :: typed_resolve => domain_specific_standard_variable_typed_resolve
   end type

   type, extends(type_domain_specific_standard_variable) :: type_interior_standard_variable
   end type

   type, extends(type_domain_specific_standard_variable) :: type_horizontal_standard_variable
   end type

   type, extends(type_horizontal_standard_variable) :: type_surface_standard_variable
   end type

   type, extends(type_horizontal_standard_variable) :: type_bottom_standard_variable
   end type

   type, extends(type_domain_specific_standard_variable) :: type_global_standard_variable
   end type

   type, extends(type_base_standard_variable) :: type_universal_standard_variable
      logical                                                    :: conserved = .false.          ! Whether this variable should be included in lists of conserved quantities.
      type (type_interior_standard_variable),   pointer, private :: pin_interior   => null()
      type (type_horizontal_standard_variable), pointer, private :: pat_interfaces => null()
      type (type_surface_standard_variable),    pointer, private :: pat_surface    => null()
      type (type_bottom_standard_variable),     pointer, private :: pat_bottom     => null()
   contains
      procedure :: typed_resolve => universal_standard_variable_typed_resolve
      procedure :: in_interior   => universal_standard_variable_in_interior
      procedure :: at_interfaces => universal_standard_variable_at_interfaces
      procedure :: at_surface    => universal_standard_variable_at_surface
      procedure :: at_bottom     => universal_standard_variable_at_bottom
   end type

   type type_standard_variable_node
      class (type_base_standard_variable), pointer :: p    => null()
      type (type_standard_variable_node),  pointer :: next => null()
   end type

   type type_standard_variable_set
      type (type_standard_variable_node), pointer :: first => null()
   contains
      procedure :: contains_variable => standard_variable_set_contains_variable
      procedure :: contains_name     => standard_variable_set_contains_name
      generic   :: contains => contains_variable, contains_name
      procedure :: add      => standard_variable_set_add
      procedure :: update   => standard_variable_set_update
      procedure :: finalize => standard_variable_set_finalize
   end type

   type type_standard_variable_collection
      type (type_standard_variable_node), pointer :: first => null()
#include "standard_variables.h"
   end type

   ! Single instance of the collection that contains all standard variables.
   type (type_standard_variable_collection), save :: standard_variables

contains

   function base_standard_variable_resolve(self) result(p)
      class (type_base_standard_variable), intent(in), target :: self
      class (type_base_standard_variable), pointer            :: p

      type (type_standard_variable_node), pointer :: node

      if (self%resolved) then
         p => self
         return
      end if

      node => standard_variables%first
      do while (associated(node))
         if (compare(self, node%p)) then
            p => node%p
            return
         end if
         node => node%next
      end do

      allocate(p, source=self)
      call add(p)

   contains

      logical function compare(variable1, variable2)
         class (type_base_standard_variable), intent(in) :: variable1, variable2

         ! Compare the type of the standard variables.
         compare = same_type_as(variable1, variable2) .or. extends_type_of(variable1, variable2) .or. extends_type_of(variable2, variable1)

         ! Compare the metadata of the standard variables.
         if (compare) compare = (variable1%name  == '' .or. variable2%name  == '' .or. variable1%name  == variable2%name ) &
                          .and. (variable1%units == '' .or. variable2%units == '' .or. variable1%units == variable2%units)
      end function

   end function base_standard_variable_resolve

   subroutine add_child(standard_variable, name, units, universal)
      class (type_domain_specific_standard_variable), target :: standard_variable
      character(len=*), intent(in)                           :: name, units
      type (type_universal_standard_variable),        target :: universal
      standard_variable%name = name
      standard_variable%units = units
      standard_variable%aggregate_variable = universal%aggregate_variable
      standard_variable%universal => universal
      call add(standard_variable)
   end subroutine

   recursive subroutine add(standard_variable)
      class (type_base_standard_variable), target, intent(inout) :: standard_variable

      type (type_standard_variable_node), pointer :: node

      select type (standard_variable)
      class is (type_universal_standard_variable)
         allocate(standard_variable%pin_interior, standard_variable%pat_surface, standard_variable%pat_bottom, standard_variable%pat_interfaces)
         call add_child(standard_variable%pin_interior,   trim(standard_variable%name),                     standard_variable%units, standard_variable)
         call add_child(standard_variable%pat_surface,    trim(standard_variable%name) // '_at_surface',    trim(standard_variable%units) // '*m', standard_variable)
         call add_child(standard_variable%pat_bottom,     trim(standard_variable%name) // '_at_bottom',     trim(standard_variable%units) // '*m', standard_variable)
         call add_child(standard_variable%pat_interfaces, trim(standard_variable%name) // '_at_interfaces', trim(standard_variable%units) // '*m', standard_variable)
      end select

      allocate(node)
      node%p => standard_variable
      node%next => standard_variables%first
      standard_variables%first => node
      standard_variable%resolved = .true.
   end subroutine

   function universal_standard_variable_typed_resolve(self) result(p)
      class (type_universal_standard_variable), target  :: self
      class (type_universal_standard_variable), pointer :: p
      class (type_base_standard_variable), pointer :: presolved
      presolved => self%resolve()
      select type (presolved)
      class is (type_universal_standard_variable)
         p => presolved
#ifndef NDEBUG
      class default
         write (*,*) 'universal_standard_variable_typed_resolve: BUG wrong type returned'
         stop 1
#endif
      end select
   end function

   function universal_standard_variable_in_interior(self) result(p)
      class (type_universal_standard_variable), target  :: self
      class (type_interior_standard_variable), pointer :: p
      class (type_universal_standard_variable), pointer  :: presolved
      presolved => self%typed_resolve()
      p => presolved%pin_interior
   end function

   function universal_standard_variable_at_surface(self) result(p)
      class (type_universal_standard_variable), target  :: self
      class (type_surface_standard_variable), pointer :: p
      class (type_universal_standard_variable), pointer  :: presolved
      presolved => self%typed_resolve()
      p => presolved%pat_surface
   end function

   function universal_standard_variable_at_bottom(self) result(p)
      class (type_universal_standard_variable), target  :: self
      class (type_bottom_standard_variable), pointer :: p
      class (type_universal_standard_variable), pointer  :: presolved
      presolved => self%typed_resolve()
      p => presolved%pat_bottom
   end function

   function universal_standard_variable_at_interfaces(self) result(p)
      class (type_universal_standard_variable), target  :: self
      class (type_horizontal_standard_variable), pointer :: p
      class (type_universal_standard_variable), pointer  :: presolved
      presolved => self%typed_resolve()
      p => presolved%pat_interfaces
   end function

   function domain_specific_standard_variable_typed_resolve(self) result(p)
      class (type_domain_specific_standard_variable), target  :: self
      class (type_domain_specific_standard_variable), pointer :: p
      class (type_base_standard_variable), pointer :: pbase
      pbase => self%resolve()
      select type (pbase)
      class is (type_domain_specific_standard_variable)
         p => pbase
#ifndef NDEBUG
      class default
         write (*,*) 'domain_specific_standard_variable_typed_resolve: BUG wrong type returned'
         stop 1
#endif
      end select
   end function

   subroutine base_standard_variable_assert_resolved(self)
      class (type_base_standard_variable), intent(in) :: self

      if (self%resolved) return
      write (*,*) 'FATAL ERROR: standard_variable_collection_assert_contains: "' // trim(self%name) // '" not in standard variable collection."'
      stop 1
   end subroutine

   subroutine initialize_standard_variables()
#include "standard_variable_assignments.h"
   end subroutine

   logical function standard_variable_set_contains_variable(self, standard_variable)
      class (type_standard_variable_set),  intent(in) :: self
      class (type_base_standard_variable), target     :: standard_variable

      type (type_standard_variable_node), pointer :: node

#ifndef NDEBUG
      call standard_variable%assert_resolved()
#endif
      standard_variable_set_contains_variable = .true.
      node => self%first
      do while (associated(node))
#ifndef NDEBUG
         call node%p%assert_resolved()
#endif
         if (associated(node%p, standard_variable)) return
         node => node%next
      end do
      standard_variable_set_contains_variable = .false.
   end function standard_variable_set_contains_variable

   logical function standard_variable_set_contains_name(self, name)
      class (type_standard_variable_set), intent(in) :: self
      character(len=*),                   intent(in) :: name

      type (type_standard_variable_node), pointer :: node

      standard_variable_set_contains_name = .true.
      node => self%first
      do while (associated(node))
         if (node%p%name == name) return
         node => node%next
      end do
      standard_variable_set_contains_name = .false.
   end function standard_variable_set_contains_name

   subroutine standard_variable_set_add(self, standard_variable)
      class (type_standard_variable_set),  intent(inout) :: self
      class (type_base_standard_variable), target        :: standard_variable

      type (type_standard_variable_node), pointer :: node

#ifndef NDEBUG
      call standard_variable%assert_resolved()
#endif

      if (self%contains(standard_variable)) return

      if (.not. associated(self%first)) then
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
      node%p => standard_variable
   end subroutine standard_variable_set_add

   subroutine standard_variable_set_update(self, other)
      class (type_standard_variable_set), intent(inout) :: self
      class (type_standard_variable_set), intent(in)    :: other

      type (type_standard_variable_node), pointer :: node

      node => other%first
      do while (associated(node))
         call self%add(node%p)
         node => node%next
      end do
   end subroutine standard_variable_set_update

   subroutine standard_variable_set_finalize(self)
      class (type_standard_variable_set), intent(inout) :: self

      type (type_standard_variable_node), pointer :: node, next_node

      node => self%first
      do while (associated(node))
         next_node => node%next
         deallocate(node)
         node => next_node
      end do
      self%first => null()
   end subroutine standard_variable_set_finalize

end module fabm_standard_variables

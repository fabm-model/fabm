#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_global_types

   use fabm_types, only: rke, type_base_model, type_variable_id

   implicit none

   private

   public type_domain, type_catalog, type_store
   public type_global_model
   public type_global_interior_dependency_id, type_global_horizontal_dependency_id
   public type_global_interior_variable_id, type_global_horizontal_variable_id, type_global_scalar_variable_id

   ! --------------------------------------------------------------------------
   ! Derived type for model domain
   ! (spatial extent, masks, indices of surface and bottom layers)
   ! --------------------------------------------------------------------------

   type type_domain
      ! Information about the model domain
      integer :: shape(_FABM_DIMENSION_COUNT_)
      integer :: start(_FABM_DIMENSION_COUNT_)
      integer :: stop(_FABM_DIMENSION_COUNT_)
      integer :: horizontal_shape(_HORIZONTAL_DIMENSION_COUNT_)

#ifdef _HAS_MASK_
#  ifndef _FABM_HORIZONTAL_MASK_
      _FABM_MASK_TYPE_, pointer _ATTRIBUTES_GLOBAL_ :: mask => null()
#  endif
      _FABM_MASK_TYPE_, pointer _ATTRIBUTES_GLOBAL_HORIZONTAL_ :: mask_hz => null()
#endif

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==-1
      integer, pointer _ATTRIBUTES_GLOBAL_HORIZONTAL_ :: bottom_indices => null()
#  endif
#endif
   end type

   ! --------------------------------------------------------------------------
   ! Derived types for catalog with pointers to all available fields
   ! --------------------------------------------------------------------------

   type type_interior_data_pointer
      real(rke), pointer _ATTRIBUTES_GLOBAL_ :: p => null()
   end type

   type type_horizontal_data_pointer
      real(rke), pointer _ATTRIBUTES_GLOBAL_HORIZONTAL_ :: p => null()
   end type

   type type_scalar_data_pointer
      real(rke), pointer :: p => null()
   end type

   type type_catalog
      type (type_interior_data_pointer),   allocatable :: interior(:)
      type (type_horizontal_data_pointer), allocatable :: horizontal(:)
      type (type_scalar_data_pointer),     allocatable :: scalar(:)
      integer, allocatable :: interior_sources(:)
      integer, allocatable :: horizontal_sources(:)
      integer, allocatable :: scalar_sources(:)
   end type

   ! --------------------------------------------------------------------------
   ! Derived type for variable store
   ! (spatially explicit model outputs needed by other BGC modules or host)
   ! --------------------------------------------------------------------------

   type type_store
      real(rke), allocatable _DIMENSION_GLOBAL_PLUS_1_            :: interior
      real(rke), allocatable _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: horizontal
      real(rke), allocatable, dimension(:)                        :: scalar
      real(rke), allocatable                                      :: interior_fill_value(:)
      real(rke), allocatable                                      :: horizontal_fill_value(:)
      real(rke), allocatable                                      :: scalar_fill_value(:)
      real(rke), allocatable                                      :: interior_missing_value(:)
      real(rke), allocatable                                      :: horizontal_missing_value(:)
   end type


   type, extends(type_variable_id) :: type_global_interior_dependency_id
      integer :: icatalog = -1
   end type

   type, extends(type_variable_id) :: type_global_horizontal_dependency_id
      integer :: icatalog = -1
   end type

   type, extends(type_variable_id) :: type_global_interior_variable_id
      real(rke), pointer _ATTRIBUTES_GLOBAL_CONTIGUOUS_ :: p => null()
   end type

   type, extends(type_variable_id) :: type_global_horizontal_variable_id
      real(rke), pointer _ATTRIBUTES_GLOBAL_HORIZONTAL_CONTIGUOUS_ :: p => null()
   end type

   type, extends(type_variable_id) :: type_global_scalar_variable_id
      real(rke), pointer :: p => null()
   end type

   type, extends(type_base_model) :: type_global_model
   contains
      procedure :: set_data
      procedure :: update
   end type

contains

   subroutine set_data(self, store, seconds_per_time_unit)
      class (type_global_model), intent(inout) :: self
      type (type_store), target                :: store
      real(rke), optional, intent(in)          :: seconds_per_time_unit
   end subroutine

   subroutine update(self, catalog _POSTARG_LOCATION_RANGE_, time)
      class (type_global_model), intent(in) :: self
      type (type_catalog),       intent(in) :: catalog
      _DECLARE_ARGUMENTS_LOCATION_RANGE_
      real(rke), optional,       intent(in) :: time
   end subroutine

end module

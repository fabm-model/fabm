#include "fabm_driver.h"

module fabm_builtin_models
   use fabm_types
   use fabm_builtin_scale
   use fabm_builtin_sum
   use fabm_builtin_constant
   use fabm_builtin_time_filter
   use fabm_builtin_depth_integral
   use fabm_builtin_source
   use fabm_builtin_relaxation
   use fabm_builtin_depth_mapping
   use fabm_builtin_tracer

   implicit none

   private

   public type_weighted_sum, type_horizontal_weighted_sum
   public type_depth_integral, type_bounded_depth_integral
   public copy_fluxes, copy_horizontal_fluxes
   public type_surface_source

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: builtin_factory

   type, extends(type_base_model) :: type_horizontal_layer
      type (type_dependency_id)                     :: id_source
      type (type_horizontal_diagnostic_variable_id) :: id_result
   contains
      procedure :: after_coupling  => horizontal_layer_after_coupling
   end type

   type,extends(type_horizontal_layer) :: type_bottom_field
   contains
      procedure :: initialize => bottom_field_initialize
      procedure :: do_bottom  => bottom_field_do_bottom
   end type

   type, extends(type_horizontal_layer) :: type_surface_field
   contains
      procedure :: initialize => surface_field_initialize
      procedure :: do_surface => surface_field_do_surface
   end type

   type, extends(type_base_model) :: type_column_projection
      type (type_horizontal_dependency_id) :: id_source
      type (type_diagnostic_variable_id)   :: id_result
   contains
      procedure :: initialize => column_projection_initialize
      procedure :: do         => column_projection_do
   end type

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('tracer');                   allocate(type_tracer::model)
         case ('bulk_constant');            allocate(type_interior_constant::model)
         case ('interior_constant');        allocate(type_interior_constant::model)
         case ('horizontal_constant');      allocate(type_horizontal_constant::model)
         case ('bottom_constant');          allocate(type_bottom_constant::model)
         case ('surface_constant');         allocate(type_surface_constant::model)
         case ('global_constant');          allocate(type_global_constant::model)
         case ('surface_flux');             allocate(type_constant_surface_flux::model)
         case ('constant_surface_flux');    allocate(type_constant_surface_flux::model)
         case ('external_surface_flux');    allocate(type_external_surface_flux::model)
         case ('external_bottom_flux');     allocate(type_external_bottom_flux::model)
         case ('interior_source');          allocate(type_interior_source::model)
         case ('bottom_source');            allocate(type_bottom_source::model)
         case ('surface_source');           allocate(type_surface_source::model)
         case ('interior_relaxation');      allocate(type_interior_relaxation::model)
         case ('column_projection');        allocate(type_column_projection::model)
         case ('weighted_sum');             allocate(type_weighted_sum::model)
         case ('horizontal_weighted_sum');  allocate(type_horizontal_weighted_sum::model)
         case ('bottom_layer');             allocate(type_bottom_field::model)
         case ('surface_layer');            allocate(type_surface_field::model)
         case ('vertical_integral');        allocate(type_depth_integral::model)
         case ('bounded_vertical_integral');allocate(type_bounded_depth_integral::model)
         case ('interior_temporal_mean');   allocate(type_interior_temporal_mean::model)
         case ('surface_temporal_mean');    allocate(type_surface_temporal_mean::model)
         case ('bottom_temporal_mean');     allocate(type_bottom_temporal_mean::model)
         case ('surface_temporal_maximum'); allocate(type_surface_temporal_maximum::model)
         case ('depth_integrated_particle_override'); allocate(type_depth_integrated_particle_override::model)
         case ('vertical_depth_range');     allocate(type_vertical_depth_range::model)
         ! Add new examples models here
      end select

   end subroutine

   subroutine horizontal_layer_after_coupling(self)
      class (type_horizontal_layer), intent(inout) :: self

      if (associated(self%id_result%link%target, self%id_result%link%original)) self%id_result%link%target%units = self%id_source%link%target%units
   end subroutine horizontal_layer_after_coupling

   subroutine bottom_field_initialize(self,configunit)
      class (type_bottom_field), intent(inout), target :: self
      integer,                   intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_bottom/))
      call self%register_diagnostic_variable(self%id_result, 'result', '', 'bottom values', source=source_do_bottom)
      call self%register_dependency(self%id_source, 'source', '', 'interior values')
   end subroutine bottom_field_initialize

   subroutine bottom_field_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_bottom_field), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: value

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_source, value)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_result, value)
      _HORIZONTAL_LOOP_END_
   end subroutine bottom_field_do_bottom

   subroutine surface_field_initialize(self,configunit)
      class (type_surface_field), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_surface/))
      call self%register_diagnostic_variable(self%id_result, 'result', '', 'surface values', source=source_do_surface)
      call self%register_dependency(self%id_source, 'source', '', 'interior values')
   end subroutine surface_field_initialize

   subroutine surface_field_do_surface(self,_ARGUMENTS_DO_SURFACE_)
      class (type_surface_field), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: value

      _HORIZONTAL_LOOP_BEGIN_
         _GET_(self%id_source, value)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_result, value)
      _HORIZONTAL_LOOP_END_
   end subroutine surface_field_do_surface

   subroutine column_projection_initialize(self,configunit)
      class (type_column_projection),intent(inout),target :: self
      integer,                       intent(in)           :: configunit

      call self%register_implemented_routines((/source_do/))
      call self%register_dependency(self%id_source,'source', '', 'horizontal source')
      call self%register_diagnostic_variable(self%id_result,'result', '', 'interior result')
   end subroutine column_projection_initialize

   subroutine column_projection_do(self,_ARGUMENTS_DO_)
      class (type_column_projection), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: value

      _LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_source,value)
         _SET_DIAGNOSTIC_(self%id_result,value)
      _LOOP_END_
   end subroutine column_projection_do

end module fabm_builtin_models

#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_v0_compatibility

   use fabm, type_bulk_variable_id => type_fabm_interior_variable_id, type_horizontal_variable_id => type_fabm_horizontal_variable_id, type_scalar_variable_id => type_fabm_scalar_variable_id, type_external_variable => type_fabm_variable, type_horizontal_state_variable_info => type_fabm_horizontal_state_variable
   use fabm_config, only: fabm_configure_model
   use fabm_properties, only: type_property_dictionary
   use fabm_types, only: rke, attribute_length, source_get_light_extinction, source_get_albedo, source_get_drag, type_base_model, type_model_list_node, type_bulk_standard_variable, type_interior_standard_variable, standard_variables
   use fabm_debug
   use fabm_job, only: type_job, type_call
   use fabm_driver, only: driver
   use fabm_work

   implicit none

   public fabm_create_model_from_yaml_file

   public type_bulk_variable_id, type_external_variable, type_model, standard_variables

   public fabm_initialize, fabm_finalize, fabm_set_domain, fabm_check_ready, fabm_update_time
   public fabm_initialize_state, fabm_initialize_surface_state, fabm_initialize_bottom_state

   ! Process rates and diagnostics for pelagic, surface, bottom.
   public fabm_do, fabm_do_surface, fabm_do_bottom

   ! Vertical movement, light attenuation, feedbacks to drag and albedo
   public fabm_get_vertical_movement, fabm_get_light_extinction, fabm_get_albedo, fabm_get_drag, fabm_get_light

   ! Bookkeeping
   public fabm_check_state, fabm_check_surface_state, fabm_check_bottom_state
   public fabm_get_conserved_quantities, fabm_get_horizontal_conserved_quantities

   ! Management of model variables: retrieve identifiers, get and set data.
   public fabm_get_variable_name, fabm_is_variable_used
   public fabm_get_interior_diagnostic_data, fabm_get_horizontal_diagnostic_data

#ifdef _HAS_MASK_
   ! Set spatial mask
   public fabm_set_mask
#endif

   type, extends(type_fabm_model) :: type_model
      ! Short names for interior state/diagnostic variables
      type (type_fabm_interior_state_variable),      pointer, dimension(:) :: state_variables      => null()
      type (type_fabm_interior_diagnostic_variable), pointer, dimension(:) :: diagnostic_variables => null()

      type (type_job) :: get_light_extinction_job
      type (type_job) :: get_albedo_job
      type (type_job) :: get_drag_job
      type (type_bulk_variable_id)       :: extinction_id
      type (type_horizontal_variable_id) :: albedo_id
      type (type_horizontal_variable_id) :: surface_drag_id
   contains
      procedure :: initialize => fabm_initialize
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      procedure :: set_surface_index
#  if _FABM_BOTTOM_INDEX_==0
      procedure :: set_bottom_index
#  endif
#endif
      procedure :: get_bulk_variable_id_by_name => fabm_get_bulk_variable_id_by_name
      procedure :: get_bulk_variable_id_sn => fabm_get_bulk_variable_id_sn
      generic :: get_bulk_variable_id => get_bulk_variable_id_by_name, get_bulk_variable_id_sn
   end type

   ! Subroutine calculating local temporal derivatives either as a right-hand side vector,
   ! or production/destruction matrices.
   interface fabm_do
      module procedure fabm_do_rhs
      module procedure fabm_do_ppdd
   end interface

   ! Subroutine calculating local temporal derivatives of bottom layer (benthos & pelagic)
   ! either as a right-hand side vector, or production/destruction matrices.
   interface fabm_do_bottom
      module procedure fabm_do_bottom_rhs
      module procedure fabm_do_bottom_ppdd
   end interface

   interface fabm_update_time
      module procedure fabm_update_time1
      module procedure fabm_update_time2
   end interface

contains

   subroutine fabm_create_model_from_yaml_file(model, path, do_not_initialize, parameters, unit)
      type (type_model),                         intent(out) :: model
      character(len=*),                optional, intent(in)  :: path
      logical,                         optional, intent(in)  :: do_not_initialize
      type (type_property_dictionary), optional, intent(in)  :: parameters
      integer,                         optional, intent(in)  :: unit

      logical :: initialize

      ! Make sure the library is initialized.
      call fabm_initialize_library()

      call fabm_configure_model(model%root, model%schedules, model%log, path, parameters, unit)

      ! Initialize model tree
      initialize = .true.
      if (present(do_not_initialize)) initialize = .not. do_not_initialize
      if (initialize) call model%initialize()
   end subroutine

   subroutine fabm_initialize(self)
      class (type_model), target, intent(inout) :: self

      call self%type_fabm_model%initialize()

      call self%job_manager%create(self%get_light_extinction_job, 'get_light_extinction', source=source_get_light_extinction, previous=self%prepare_inputs_job)
      call self%job_manager%create(self%get_albedo_job,'get_albedo', source=source_get_albedo, previous=self%finalize_outputs_job)
      call self%job_manager%create(self%get_drag_job,'get_drag', source=source_get_drag, previous=self%finalize_outputs_job)
      self%extinction_id = self%get_interior_variable_id(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      self%albedo_id = self%get_horizontal_variable_id(standard_variables%surface_albedo)
      self%surface_drag_id = self%get_horizontal_variable_id(standard_variables%surface_drag_coefficient_in_air)
      _ASSERT_(associated(self%extinction_id%variable), 'fabm_initialize', 'BUG: variable attenuation_coefficient_of_photosynthetic_radiative_flux not found')
      _ASSERT_(associated(self%albedo_id%variable), 'fabm_initialize', 'BUG: variable surface_albedo not found')
      _ASSERT_(associated(self%surface_drag_id%variable), 'fabm_initialize', 'BUG: variable surface_drag_coefficient_in_air not found')
      call self%get_light_extinction_job%request_variable(self%extinction_id%variable)
      call self%get_albedo_job%request_variable(self%albedo_id%variable)
      call self%get_drag_job%request_variable(self%surface_drag_id%variable)

      ! Short names for interior variables
      self%state_variables => self%interior_state_variables
      self%diagnostic_variables => self%interior_diagnostic_variables
   end subroutine

   subroutine fabm_finalize(self)
      class (type_model), target, intent(inout) :: self
      call self%finalize()
   end subroutine

   subroutine fabm_set_domain(self _POSTARG_LOCATION_, seconds_per_time_unit)
      class (type_model), target, intent(inout) :: self
      _DECLARE_ARGUMENTS_LOCATION_
      real(rke), optional,        intent(in)    :: seconds_per_time_unit
      call self%set_domain(_PREARG_LOCATION_ seconds_per_time_unit)
   end subroutine

#ifdef _FABM_DEPTH_DIMENSION_INDEX_

   ! --------------------------------------------------------------------------
   ! set_surface_index: set vertical index of surface layer
   ! --------------------------------------------------------------------------
   subroutine set_surface_index(self, index)
      class (type_model), intent(inout) :: self
      integer,            intent(in)    :: index

      if (self%status < status_set_domain_done) &
         call driver%fatal_error('set_surface_index', 'set_domain has not yet been called on this model object.')
      if (index < 1) &
         call driver%fatal_error('set_surface_index', 'provided index must equal or exceed 1.')
      if (index > self%domain%shape(_FABM_DEPTH_DIMENSION_INDEX_)) &
         call driver%fatal_error('set_surface_index', 'provided index exceeds size of the depth dimension.')

#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
      self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_) = index
#    else
      self%domain%start(_FABM_DEPTH_DIMENSION_INDEX_) = index
#    endif
   end subroutine set_surface_index

#  if _FABM_BOTTOM_INDEX_==0

   ! --------------------------------------------------------------------------
   ! set_bottom_index: set vertical index of bottom layer
   ! --------------------------------------------------------------------------
   subroutine set_bottom_index(self, index)
      class (type_model), intent(inout) :: self
      integer,            intent(in)    :: index

      if (self%status < status_set_domain_done) &
         call driver%fatal_error('set_bottom_index', 'set_domain has not yet been called on this model object.')
      if (index < 1) &
         call driver%fatal_error('set_bottom_index', 'provided index must equal or exceed 1.')
      if (index > self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)) &
         call driver%fatal_error('set_bottom_index', 'provided index exceeds size of the depth dimension.')

#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
      self%domain%start(_FABM_DEPTH_DIMENSION_INDEX_) = index
#    else
      self%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_) = index
#    endif
   end subroutine set_bottom_index
#  endif

#endif

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
   subroutine fabm_set_mask(self, mask_hz)
      class (type_model), target, intent(inout)                             :: self
      _FABM_MASK_TYPE_,   target, intent(in) _ATTRIBUTES_GLOBAL_HORIZONTAL_ :: mask_hz
      call self%set_mask(mask_hz)
   end subroutine
#  else
   subroutine fabm_set_mask(self, mask, mask_hz)
      class (type_model), target, intent(inout)                             :: self
      _FABM_MASK_TYPE_,   target, intent(in) _ATTRIBUTES_GLOBAL_            :: mask
      _FABM_MASK_TYPE_,   target, intent(in) _ATTRIBUTES_GLOBAL_HORIZONTAL_ :: mask_hz
      call self%set_mask(mask, mask_hz)
   end subroutine
#  endif
#endif

   subroutine fabm_check_ready(self)
      class (type_model), intent(inout), target :: self
      call self%start()
   end subroutine

   subroutine fabm_initialize_state(self _POSTARG_INTERIOR_IN_)
      class (type_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      call self%initialize_interior_state(_ARG_INTERIOR_IN_)
   end subroutine

   subroutine fabm_initialize_bottom_state(self _POSTARG_HORIZONTAL_IN_)
      class (type_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      call self%initialize_bottom_state(_ARG_HORIZONTAL_IN_)
   end subroutine

   subroutine fabm_initialize_surface_state(self _POSTARG_HORIZONTAL_IN_)
      class (type_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      call self%initialize_surface_state(_ARG_HORIZONTAL_IN_)
   end subroutine

   subroutine fabm_check_state(self _POSTARG_INTERIOR_IN_, repair, valid)
      class (type_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      logical,            intent(in)    :: repair
      logical,            intent(out)   :: valid
      call self%check_interior_state(_PREARG_INTERIOR_IN_ repair, valid)
   end subroutine

   subroutine fabm_check_bottom_state(self _POSTARG_HORIZONTAL_IN_, repair, valid)
      class (type_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      logical,            intent(in)    :: repair
      logical,            intent(out)   :: valid
      call self%check_bottom_state(_PREARG_HORIZONTAL_IN_ repair, valid)
   end subroutine

   subroutine fabm_check_surface_state(self _POSTARG_HORIZONTAL_IN_, repair, valid)
      class (type_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      logical,            intent(in)    :: repair
      logical,            intent(out)   :: valid
      call self%check_surface_state(_PREARG_HORIZONTAL_IN_ repair, valid)
   end subroutine

   subroutine fabm_get_vertical_movement(self _POSTARG_INTERIOR_IN_, velocity)
      class (type_model),                     intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_1_, intent(out)   :: velocity
      call self%get_vertical_movement(_PREARG_INTERIOR_IN_ velocity)
   end subroutine

   subroutine fabm_get_conserved_quantities(self _POSTARG_INTERIOR_IN_, sums)
      class (type_model),                     intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_1_, intent(out)   :: sums
      call self%get_interior_conserved_quantities(_PREARG_INTERIOR_IN_ sums)
   end subroutine

   subroutine fabm_get_horizontal_conserved_quantities(self _POSTARG_HORIZONTAL_IN_, sums)
      class (type_model),                                intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(out)   :: sums
      call self%get_horizontal_conserved_quantities(_PREARG_HORIZONTAL_IN_ sums)
   end subroutine

   subroutine fabm_do_rhs(self _POSTARG_INTERIOR_IN_, dy)
      class (type_model),                     intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_1_, intent(inout) :: dy
      call self%get_interior_sources(_PREARG_INTERIOR_IN_ dy)
   end subroutine

   subroutine fabm_do_ppdd(self _POSTARG_INTERIOR_IN_, pp, dd)
      class (type_model),                     intent(inout) :: self
     _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_2_, intent(inout) :: pp, dd
      call self%get_interior_sources(_PREARG_INTERIOR_IN_ pp, dd)
   end subroutine

   subroutine fabm_do_bottom_rhs(self _POSTARG_HORIZONTAL_IN_, flux_pel, flux_ben)
      class (type_model),                                intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(inout) :: flux_pel, flux_ben
      call self%get_bottom_sources(_PREARG_HORIZONTAL_IN_ flux_pel, flux_ben)
   end subroutine

   subroutine fabm_do_bottom_ppdd(self _POSTARG_HORIZONTAL_IN_, pp, dd, benthos_offset)
      class (type_model),                                intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      integer,                                           intent(in)    :: benthos_offset
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_2_, intent(inout) :: pp, dd
      call self%get_bottom_sources(_PREARG_HORIZONTAL_IN_ pp, dd, benthos_offset)
   end subroutine

   subroutine fabm_do_surface(self _POSTARG_HORIZONTAL_IN_, flux_pel, flux_sf)
      class (type_model),                                intent(inout)         :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(out)           :: flux_pel
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(out), optional :: flux_sf
      call self%get_surface_sources(_PREARG_HORIZONTAL_IN_ flux_pel, flux_sf)
   end subroutine

   function fabm_get_interior_diagnostic_data(self, index) result(dat)
      class (type_model), intent(in)         :: self
      integer,            intent(in)         :: index
      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: dat
      dat => self%get_interior_diagnostic_data(index)
   end function

   function fabm_get_horizontal_diagnostic_data(self, index) result(dat)
      class (type_model), intent(in)                    :: self
      integer,            intent(in)                    :: index
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: dat
      dat => self%get_horizontal_diagnostic_data(index)
   end function

   function fabm_is_variable_used(id) result(used)
      class (type_fabm_variable_id), intent(in) :: id
      logical                                   :: used
      used = associated(id%variable)
      if (used) used = .not. id%variable%read_indices%is_empty()
   end function

   function fabm_get_variable_name(self, id) result(name)
      class (type_model),            intent(in) :: self
      class (type_fabm_variable_id), intent(in) :: id
      character(len=attribute_length)           :: name
      name = self%get_variable_name(id)
   end function

   subroutine fabm_get_light(self _POSTARG_VERTICAL_IN_)
      class (type_model), intent(inout) :: self
      _DECLARE_ARGUMENTS_VERTICAL_IN_
   end subroutine

   subroutine fabm_get_light_extinction(self _POSTARG_INTERIOR_IN_, extinction)
      class (type_model),              intent(inout) :: self
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_, intent(out)   :: extinction

      _DECLARE_INTERIOR_INDICES_

#ifndef NDEBUG
      call check_interior_location(self%domain%start, self%domain%stop _POSTARG_INTERIOR_IN_, 'fabm_get_light_extinction')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_extents_1d(extinction, _STOP_ - _START_ + 1, 'fabm_get_light_extinction', 'extinction', 'stop-start+1')
#  endif
#endif

      call process_interior_slice(self%get_light_extinction_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_int _POSTARG_INTERIOR_IN_)

      ! Copy light extinction from write cache to output array provided by host
      _UNPACK_(self%cache_int%write, self%extinction_id%variable%write_indices%value, extinction, self%cache_int, 0.0_rke)
   end subroutine fabm_get_light_extinction

   subroutine fabm_get_drag(self _POSTARG_HORIZONTAL_IN_, drag)
      class (type_model),                         intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_, intent(out)   :: drag

      _DECLARE_HORIZONTAL_INDICES_

#ifndef NDEBUG
       call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'fabm_get_drag')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
       call check_extents_1d(drag, _STOP_ - _START_ + 1, 'fabm_get_drag', 'drag', 'stop-start+1')
#  endif
#endif

      call process_horizontal_slice(self%get_drag_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      _HORIZONTAL_UNPACK_(self%cache_hz%write_hz, self%surface_drag_id%variable%write_indices%value, drag, self%cache_hz, 0.0_rke)
      drag = drag + 1.0_rke
   end subroutine fabm_get_drag

   subroutine fabm_get_albedo(self _POSTARG_HORIZONTAL_IN_, albedo)
      class (type_model),                         intent(inout) :: self
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_, intent(out)   :: albedo

      _DECLARE_HORIZONTAL_INDICES_

#ifndef NDEBUG
      call check_horizontal_location(self%domain%start, self%domain%stop _POSTARG_HORIZONTAL_IN_, 'fabm_get_albedo')
#  ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_extents_1d(albedo, _STOP_ - _START_ + 1, 'fabm_get_albedo', 'albedo', 'stop-start+1')
#  endif
#endif

      call process_horizontal_slice(self%get_albedo_job%first_task, self%domain, self%catalog, self%cache_fill_values, self%store, self%cache_hz _POSTARG_HORIZONTAL_IN_)

      _HORIZONTAL_UNPACK_(self%cache_hz%write_hz, self%albedo_id%variable%write_indices%value, albedo, self%cache_hz, 0.0_rke)
   end subroutine fabm_get_albedo

   subroutine fabm_update_time1(self, t)
      class (type_model),  intent(inout) :: self
      real(rke), optional, intent(in)    :: t
      call self%prepare_inputs(t)
   end subroutine

   subroutine fabm_update_time2(self, t, year, month, day, seconds)
      class (type_model), intent(inout) :: self
      real(rke),          intent(in)    :: t
      integer,            intent(in)    :: year, month, day
      real(rke),          intent(in)    :: seconds
      call self%prepare_inputs(t, year, month, day, seconds)
   end subroutine

   function fabm_get_bulk_variable_id_by_name(self, name) result(id)
      class (type_model),      intent(in)   :: self
      character(len=*),        intent(in)   :: name
      type (type_bulk_variable_id)          :: id
      id = self%get_interior_variable_id(name)
   end function

   function fabm_get_bulk_variable_id_sn(self, standard_variable) result(id)
      class (type_model),                     intent(in) :: self
      type (type_interior_standard_variable), intent(in) :: standard_variable
      type (type_bulk_variable_id)                       :: id
      id = self%get_interior_variable_id(standard_variable)
   end function

end module

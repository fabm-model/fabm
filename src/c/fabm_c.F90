#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_c

   use iso_c_binding, only: c_double, c_int, c_char, C_NULL_CHAR, c_f_pointer, c_loc, c_ptr, c_null_ptr

   !DIR$ ATTRIBUTES DLLEXPORT :: STATE_VARIABLE,DIAGNOSTIC_VARIABLE,CONSERVED_QUANTITY

   use fabm, only: type_fabm_model, type_fabm_variable, fabm_get_version, status_start_done, fabm_create_model
   use fabm_types, only: rk => rke, attribute_length, type_model_list_node, type_base_model, &
                         factory, type_link, type_link_list, type_internal_variable, type_variable_list, type_variable_node, &
                         domain_interior, domain_horizontal, domain_scalar
   use fabm_driver, only: type_base_driver, driver
   use fabm_properties, only: type_property, type_property_dictionary
   use fabm_python_helper
   use fabm_c_helper

   implicit none

   public

   integer, parameter :: INTERIOR_STATE_VARIABLE        = 1
   integer, parameter :: SURFACE_STATE_VARIABLE         = 2
   integer, parameter :: BOTTOM_STATE_VARIABLE          = 3
   integer, parameter :: INTERIOR_DIAGNOSTIC_VARIABLE   = 4
   integer, parameter :: HORIZONTAL_DIAGNOSTIC_VARIABLE = 5
   integer, parameter :: CONSERVED_QUANTITY             = 6
   integer, parameter :: INTERIOR_DEPENDENCY            = 7
   integer, parameter :: HORIZONTAL_DEPENDENCY          = 8
   integer, parameter :: SCALAR_DEPENDENCY              = 9

   logical, save :: error_occurred = .false.
   character(len=:), allocatable, save :: error_message

   type, extends(type_base_driver) :: type_python_driver
   contains
      procedure :: fatal_error => python_driver_fatal_error
      procedure :: log_message => python_driver_log_message
   end type

   type type_model_wrapper
      class (type_fabm_model), pointer :: p => null()
      type (type_variable_list)        :: environment
      type (type_link_list)            :: coupling_link_list
      type (type_property_dictionary)  :: forced_parameters, forced_couplings
   end type

contains

   subroutine get_version(length, version_string) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_version
      integer(c_int), value, intent(in) :: length
      character(kind=c_char)            :: version_string(length)

      character(len=length-1) :: string

      call fabm_get_version(string)
      call copy_to_c_string(string, version_string)
   end subroutine get_version


   subroutine get_driver_settings(ndim, idepthdim) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_driver_settings
      integer(c_int), intent(out) :: ndim, idepthdim

      ndim = _FABM_DIMENSION_COUNT_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      idepthdim = _FABM_DEPTH_DIMENSION_INDEX_
#else
      idepthdim = -1
#endif
   end subroutine get_driver_settings

   function create_model(path _POSTARG_LOCATION_) result(ptr) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: create_model
      character(kind=c_char), target, intent(in) :: path(*)
#if _FABM_DIMENSION_COUNT_ > 0
      integer (kind=c_int), value, intent(in) :: _LOCATION_
#endif
      type(c_ptr)                                :: ptr

      type (type_model_wrapper),       pointer :: model
      character(len=attribute_length), pointer :: ppath
      class (type_property),           pointer :: property

      ! Initialize driver object used by FABM for logging/error reporting.
      if (.not. associated(driver)) allocate(type_python_driver::driver)

      ! Build FABM model tree (configuration will be read from file specified as argument).
      allocate(model)
      call c_f_pointer(c_loc(path), ppath)
      model%p => fabm_create_model(path=ppath(:index(ppath, C_NULL_CHAR) - 1), parameters=model%forced_parameters)

      ! Get a list of all parameters that had an explicit value specified.
      property => model%p%root%parameters%first
      do while (associated(property))
         if (.not. model%p%root%parameters%missing%contains(property%name)) &
            call model%forced_parameters%set_property(property)
         property => property%next
      end do

      ! Get a list of all active couplings.
      property => model%p%root%couplings%first
      do while (associated(property))
         if (.not. model%p%root%couplings%missing%contains(property%name)) &
            call model%forced_couplings%set_property(property)
         property => property%next
      end do

      ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
      call model%p%set_domain(_PREARG_LOCATION_ 1._rk)

      ! Retrieve arrays to hold values for environmental variables and corresponding metadata.
      call get_environment_metadata(model%p, model%environment)

      call get_couplings(model%p, model%coupling_link_list)

      ptr = c_loc(model)
   end function create_model

   subroutine reinitialize(model)
      type (type_model_wrapper), intent(inout) :: model

      class (type_fabm_model),     pointer :: newmodel
      type (type_model_list_node), pointer :: node
      class (type_base_model),     pointer :: childmodel
      class (type_property),       pointer :: property, next
      _DECLARE_LOCATION_

#  if _FABM_DIMENSION_COUNT_ > 0
      i__ = model%p%domain%shape(1)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 1
      j__ = model%p%domain%shape(2)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 2
      k__ = model%p%domain%shape(3)
#  endif

      ! Create new model object.
      allocate(newmodel)

      ! Transfer forced parameters to root of the model.
      call newmodel%root%parameters%update(model%forced_parameters)
      call newmodel%root%couplings%update(model%forced_couplings)

      ! Re-create original models
      node => model%p%root%children%first
      do while (associated(node))
         if (node%model%user_created) then
            call factory%create(node%model%type_name, childmodel)
            childmodel%user_created = .true.
            call newmodel%root%add_child(childmodel, node%model%name, node%model%long_name, configunit=-1)
         end if
         node => node%next
      end do

      ! Clean up old model
      call finalize(model, keep_forced=.true.)
      model%p => newmodel

      ! Initialize new model
      call model%p%initialize()

      ! Removed unused forced parameters from root model.
      property => model%p%root%parameters%first
      do while (associated(property))
         if (.not. model%p%root%parameters%retrieved%contains(property%name)) then
            next => property%next
            call model%p%root%parameters%delete(property%name)
            property => next
         else
            property => property%next
         end if
      end do

      ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
      call model%p%set_domain(_PREARG_LOCATION_ 1._rk)

      ! Retrieve arrays to hold values for environmental variables and corresponding metadata.
      call get_environment_metadata(model%p, model%environment)

      call get_couplings(model%p, model%coupling_link_list)
   end subroutine reinitialize

   subroutine start(pmodel) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: start
      type (c_ptr), intent(in), value :: pmodel

      type (type_model_wrapper), pointer :: model

      call c_f_pointer(pmodel, model)
      call model%p%start()
   end subroutine start

   integer(c_int) function get_error_state() bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_error_state
      get_error_state = logical2int(error_occurred)
   end function get_error_state

   subroutine reset_error_state() bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: reset_error_state
      error_occurred = .false.
      error_message = ''
   end subroutine reset_error_state

   subroutine get_error(length, message) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_error
      integer(c_int),        intent(in), value             :: length
      character(kind=c_char),intent(out),dimension(length) :: message
      call copy_to_c_string(error_message, message)
   end subroutine get_error

   integer(c_int) function model_count(pmodel) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: model_count
      type (c_ptr), intent(in), value :: pmodel

      type (type_model_wrapper),   pointer :: model
      type (type_model_list_node), pointer :: node

      call c_f_pointer(pmodel, model)
      model_count = 0
      node => model%p%root%children%first
      do while (associated(node))
         model_count = model_count + 1
         node => node%next
      end do
   end function model_count

   subroutine get_counts(pmodel, nstate_interior, nstate_surface, nstate_bottom, ndiagnostic_interior, ndiagnostic_horizontal, &
      ndependencies_interior, ndependencies_horizontal, ndependencies_scalar, nconserved, nparameters, ncouplings) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_counts
      type (c_ptr),   intent(in), value :: pmodel
      integer(c_int), intent(out)       :: nstate_interior, nstate_surface, nstate_bottom
      integer(c_int), intent(out)       :: ndiagnostic_interior, ndiagnostic_horizontal
      integer(c_int), intent(out)       :: ndependencies_interior, ndependencies_horizontal, ndependencies_scalar
      integer(c_int), intent(out)       :: nconserved, nparameters, ncouplings

      type (type_model_wrapper), pointer :: model
      type (type_variable_node), pointer :: node

      call c_f_pointer(pmodel, model)
      nstate_interior = size(model%p%interior_state_variables)
      nstate_surface = size(model%p%surface_state_variables)
      nstate_bottom = size(model%p%bottom_state_variables)
      ndiagnostic_interior = size(model%p%interior_diagnostic_variables)
      ndiagnostic_horizontal = size(model%p%horizontal_diagnostic_variables)

      ndependencies_interior = 0
      ndependencies_horizontal = 0
      ndependencies_scalar = 0
      node => model%environment%first
      do while (associated(node))
         select case (node%target%domain)
         case (domain_interior); ndependencies_interior = ndependencies_interior + 1
         case (domain_scalar);   ndependencies_scalar = ndependencies_scalar + 1
         case default;           ndependencies_horizontal = ndependencies_horizontal + 1
         end select
         node => node%next
      end do

      nconserved = size(model%p%conserved_quantities)
      nparameters = model%p%root%parameters%size()
      ncouplings = model%coupling_link_list%count()
   end subroutine get_counts

   subroutine get_variable_metadata(pmodel, category, index, length, name, units, long_name, path) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_variable_metadata
      type (c_ptr),           intent(in), value              :: pmodel
      integer(c_int),         intent(in), value              :: category, index, length
      character(kind=c_char), intent(out), dimension(length) :: name, units, long_name, path

      type (type_model_wrapper),  pointer :: model
      class (type_fabm_variable), pointer :: variable

      call c_f_pointer(pmodel, model)

      ! Get a pointer to the target variable
      select case (category)
      case (INTERIOR_STATE_VARIABLE)
         variable => model%p%interior_state_variables(index)
      case (SURFACE_STATE_VARIABLE)
         variable => model%p%surface_state_variables(index)
      case (BOTTOM_STATE_VARIABLE)
         variable => model%p%bottom_state_variables(index)
      case (INTERIOR_DIAGNOSTIC_VARIABLE)
         variable => model%p%interior_diagnostic_variables(index)
      case (HORIZONTAL_DIAGNOSTIC_VARIABLE)
         variable => model%p%horizontal_diagnostic_variables(index)
      case (CONSERVED_QUANTITY)
         variable => model%p%conserved_quantities(index)
      end select
      call copy_to_c_string(variable%name,            name)
      call copy_to_c_string(variable%units,           units)
      call copy_to_c_string(variable%local_long_name, long_name)
      call copy_to_c_string(variable%path,            path)
   end subroutine get_variable_metadata

   function get_variable(pmodel, category, index) bind(c) result(pvariable)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_variable
      type (c_ptr),   intent(in), value :: pmodel
      integer(c_int), intent(in), value :: category, index
      type (c_ptr)                      :: pvariable

      type (type_model_wrapper),     pointer :: model
      type (type_internal_variable), pointer :: variable
      type (type_variable_node),     pointer :: node
      integer                                :: i, domain

      call c_f_pointer(pmodel, model)

      ! Get a pointer to the target variable
      select case (category)
      case (INTERIOR_STATE_VARIABLE)
         variable => model%p%interior_state_variables(index)%target
      case (SURFACE_STATE_VARIABLE)
         variable => model%p%surface_state_variables(index)%target
      case (BOTTOM_STATE_VARIABLE)
         variable => model%p%bottom_state_variables(index)%target
      case (INTERIOR_DIAGNOSTIC_VARIABLE)
         variable => model%p%interior_diagnostic_variables(index)%target
      case (HORIZONTAL_DIAGNOSTIC_VARIABLE)
         variable => model%p%horizontal_diagnostic_variables(index)%target
      case (CONSERVED_QUANTITY)
         variable => model%p%conserved_quantities(index)%target
      case (INTERIOR_DEPENDENCY, HORIZONTAL_DEPENDENCY, SCALAR_DEPENDENCY)
         select case (category)
         case (INTERIOR_DEPENDENCY);   domain = domain_interior
         case (HORIZONTAL_DEPENDENCY); domain = domain_horizontal
         case (SCALAR_DEPENDENCY);     domain = domain_scalar
         end select
         i = 0
         variable => null()
         node => model%environment%first
         do while (associated(node))
            if (iand(node%target%domain, domain) /= 0) then
               i = i + 1
               if (i == index) variable => node%target
            end if
            node => node%next
         end do
      end select
      pvariable = c_null_ptr
      if (associated(variable)) pvariable = c_loc(variable)
   end function get_variable

   subroutine get_parameter_metadata(pmodel, index, length, name, units, long_name, typecode, has_default) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_parameter_metadata
      type (c_ptr),           intent(in), value              :: pmodel
      integer(c_int),         intent(in), value              :: index, length
      character(kind=c_char), intent(out), dimension(length) :: name, units, long_name
      integer(c_int),         intent(out)                    :: typecode, has_default

      type (type_model_wrapper), pointer :: model
      integer                            :: i
      class (type_property),     pointer :: property

      call c_f_pointer(pmodel, model)

      i = 1
      property => model%p%root%parameters%first
      do while (associated(property))
         if (index == i) exit
         property => property%next
         i = i + 1
      end do

      call copy_to_c_string(property%name, name)
      call copy_to_c_string(property%units, units)
      call copy_to_c_string(property%long_name, long_name)
      typecode = property%typecode()
      has_default = logical2int(property%has_default)
   end subroutine get_parameter_metadata

   subroutine get_coupling(pmodel, index, slave, master) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_coupling
      type (c_ptr),   intent(in), value :: pmodel
      integer(c_int), intent(in), value :: index
      type (c_ptr),   intent(out)       :: slave, master

      type (type_model_wrapper), pointer :: model
      type (type_link),pointer :: link_slave
      integer                  :: i

      call c_f_pointer(pmodel, model)
      link_slave => model%coupling_link_list%first
      do i=2,index
         link_slave => link_slave%next
      end do
      slave = c_loc(link_slave%original)
      master = c_loc(link_slave%target)
   end subroutine get_coupling

   function variable_get_suitable_masters(pmodel, pvariable) result(plist) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: variable_get_suitable_masters
      type (c_ptr), intent(in), value :: pmodel
      type (c_ptr), intent(in), value :: pvariable
      type (c_ptr)                    :: plist

      type (type_model_wrapper),     pointer :: model
      type (type_internal_variable), pointer :: variable
      type (type_link_list),         pointer :: list

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pvariable, variable)
      list => get_suitable_masters(model%p, variable)
      plist = c_loc(list)
   end function variable_get_suitable_masters

   subroutine get_model_metadata(pmodel, name, length, long_name, user_created) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_model_metadata
      type (c_ptr),           intent(in), value  :: pmodel
      character(kind=c_char), intent(in), target :: name(*)
      integer(c_int),         intent(in), value  :: length
      character(kind=c_char), intent(out)        :: long_name(length)
      integer(c_int),         intent(out)        :: user_created

      type (type_model_wrapper),       pointer :: model
      character(len=attribute_length), pointer :: pname
      class (type_base_model),         pointer :: found_model

      call c_f_pointer(pmodel, model)
      call c_f_pointer(c_loc(name), pname)
      found_model => model%p%root%find_model(pname(:index(pname, C_NULL_CHAR) - 1))
      if (.not.associated(found_model)) call driver%fatal_error('get_model_metadata', &
         'model "'//pname(:index(pname, C_NULL_CHAR) - 1) // '" not found.')
      call copy_to_c_string(found_model%long_name, long_name)
      user_created = logical2int(found_model%user_created)
   end subroutine get_model_metadata

   function c_f_pointer_interior(model, ptr) result(pdata)
      type (type_model_wrapper), intent(in) :: model
      real(c_double), intent(in), target :: ptr(*)
      real(c_double) _ATTRIBUTES_GLOBAL_, pointer :: pdata
#if _FABM_DIMENSION_COUNT_ > 0
      call c_f_pointer(c_loc(ptr), pdata, model%p%domain%shape)
#else
      call c_f_pointer(c_loc(ptr), pdata)
#endif
   end function

   function c_f_pointer_horizontal(model, ptr) result(pdata)
      type (type_model_wrapper), intent(in) :: model
      real(c_double), intent(in), target :: ptr(*)
      real(c_double) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: pdata
#if  _HORIZONTAL_DIMENSION_COUNT_ > 0
      call c_f_pointer(c_loc(ptr), pdata, model%p%domain%horizontal_shape)
#else
      call c_f_pointer(c_loc(ptr), pdata)
#endif
   end function

   subroutine link_interior_data(pmodel, pvariable, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_interior_data
      type (c_ptr),   intent(in), value  :: pmodel
      type (c_ptr),   intent(in), value  :: pvariable
      real(c_double), intent(in), target :: dat(*)

      type (type_model_wrapper),          pointer :: model
      type (type_internal_variable),      pointer :: variable
      real(c_double) _ATTRIBUTES_GLOBAL_, pointer :: interior_data

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pvariable, variable)
      interior_data => c_f_pointer_interior(model, dat)
      call model%p%link_interior_data(variable, interior_data)
   end subroutine link_interior_data

   subroutine link_horizontal_data(pmodel, pvariable, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_horizontal_data
      type (c_ptr),   intent(in), value  :: pmodel
      type (c_ptr),   intent(in), value  :: pvariable
      real(c_double), intent(in), target :: dat(*)

      type (type_model_wrapper),                     pointer :: model
      type (type_internal_variable),                 pointer :: variable
      real(c_double) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: horizontal_data

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pvariable, variable)
      horizontal_data => c_f_pointer_horizontal(model, dat)
      call model%p%link_horizontal_data(variable, horizontal_data)
   end subroutine link_horizontal_data

   subroutine link_scalar(pmodel, pvariable, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_scalar
      type (c_ptr),   intent(in), value  :: pmodel
      type (c_ptr),   intent(in), value  :: pvariable
      real(c_double), intent(in), target :: dat(*)

      type (type_model_wrapper),     pointer :: model
      type (type_internal_variable), pointer :: variable
      real(c_double),                pointer :: scalar_data

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pvariable, variable)
      call c_f_pointer(c_loc(dat), scalar_data)
      call model%p%link_scalar(variable, scalar_data)
   end subroutine link_scalar

   subroutine link_interior_state_data(pmodel, index, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_interior_state_data
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), value  :: index
      real(c_double), intent(in), target :: dat(*)

      type (type_model_wrapper),          pointer :: model
      real(c_double) _ATTRIBUTES_GLOBAL_, pointer :: dat_

      call c_f_pointer(pmodel, model)
      dat_ => c_f_pointer_interior(model, dat)
      dat_ = model%p%interior_state_variables(index)%initial_value
      call model%p%link_interior_state_data(index, dat_)
   end subroutine link_interior_state_data

   subroutine link_surface_state_data(pmodel, index, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_surface_state_data
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), value  :: index
      real(c_double), intent(in), target :: dat(*)

      type (type_model_wrapper),                     pointer :: model
      real(c_double) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: dat_

      call c_f_pointer(pmodel, model)
      dat_ => c_f_pointer_horizontal(model, dat)
      dat_ = model%p%surface_state_variables(index)%initial_value
      call model%p%link_surface_state_data(index, dat_)
   end subroutine link_surface_state_data

   subroutine link_bottom_state_data(pmodel, index, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_bottom_state_data
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), value  :: index
      real(c_double), intent(in), target :: dat(*)

      type (type_model_wrapper),                     pointer :: model
      real(c_double) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: dat_

      call c_f_pointer(pmodel, model)
      dat_ => c_f_pointer_horizontal(model, dat)
      dat_ = model%p%bottom_state_variables(index)%initial_value
      call model%p%link_bottom_state_data(index, dat_)
   end subroutine link_bottom_state_data

   subroutine get_sources(pmodel, t, sources_interior, sources_surface, sources_bottom, do_surface, do_bottom, cell_thickness) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_sources
      type (c_ptr),   value,  intent(in) :: pmodel
      real(c_double), value,  intent(in) :: t
      real(c_double), target, intent(in) :: sources_interior(*), sources_surface(*), sources_bottom(*)
      integer(c_int), value,  intent(in) :: do_surface, do_bottom
      real(c_double), target, intent(in) :: cell_thickness(*)

      logical :: surface, bottom
      type (type_model_wrapper), pointer :: model
      real(c_double) _DIMENSION_GLOBAL_PLUS_1_, pointer :: sources_interior_
      real(c_double) _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_, pointer :: sources_surface_, sources_bottom_
      real(c_double) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_, allocatable :: fluxes
      real(c_double) _ATTRIBUTES_GLOBAL_, pointer :: cell_thickness_
      _DECLARE_LOCATION_
#  if _FABM_DIMENSION_COUNT_ > 0
      integer :: _LOCATION_RANGE_
#  endif

      call c_f_pointer(pmodel, model)
      if (model%p%status < status_start_done) then
         call driver%fatal_error('get_sources', 'start has not been called yet.')
         return
      end if

#  if _FABM_DIMENSION_COUNT_ > 0
      istart__ = model%p%domain%start(1)
      istop__ = model%p%domain%stop(1)
      i__ = model%p%domain%shape(1)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 1
      jstart__ = model%p%domain%start(2)
      jstop__ = model%p%domain%stop(2)
      j__ = model%p%domain%shape(2)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 2
      kstart__ = model%p%domain%start(3)
      kstop__ = model%p%domain%stop(3)
      k__ = model%p%domain%shape(3)
#  endif

      surface = int2logical(do_surface)
      bottom = int2logical(do_bottom)
      if ((surface .or. bottom) .and. size(model%p%interior_state_variables) > 0) cell_thickness_ => c_f_pointer_interior(model, cell_thickness)

      call c_f_pointer(c_loc(sources_interior), sources_interior_, (/_PREARG_LOCATION_ size(model%p%interior_state_variables)/))
      call c_f_pointer(c_loc(sources_surface), sources_surface_, (/_PREARG_HORIZONTAL_LOCATION_ size(model%p%surface_state_variables)/))
      call c_f_pointer(c_loc(sources_bottom), sources_bottom_, (/_PREARG_HORIZONTAL_LOCATION_ size(model%p%bottom_state_variables)/))
#ifdef _HORIZONTAL_IS_VECTORIZED_
      allocate(fluxes(_ITERATOR_, size(model%p%interior_state_variables)))
#else
      allocate(fluxes(size(model%p%interior_state_variables)))
#endif

      if (t < 0) then
         call model%p%prepare_inputs()
      else
         call model%p%prepare_inputs(t)
      end if

      sources_interior_ = 0.0_rk
      sources_surface_ = 0.0_rk
      sources_bottom_ = 0.0_rk

      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%p%get_interior_sources(_PREARG_INTERIOR_IN_ sources_interior_ _INDEX_GLOBAL_INTERIOR_PLUS_1_(_START_:_STOP_,:))
      _END_OUTER_INTERIOR_LOOP_

      if (surface) then
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         _VERTICAL_ITERATOR_ = model%p%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)
#  else
         _VERTICAL_ITERATOR_ = model%p%domain%start(_FABM_DEPTH_DIMENSION_INDEX_)
#  endif
#endif

         _BEGIN_OUTER_HORIZONTAL_LOOP_
            fluxes = 0.0_rk
            call model%p%get_surface_sources(_PREARG_HORIZONTAL_IN_ fluxes, sources_surface_ _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_:_STOP_,:))
            if (size(model%p%interior_state_variables) > 0) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
               _DO_CONCURRENT_(_ITERATOR_,_START_,_STOP_)
                  sources_interior_(_PREARG_LOCATION_ :) = sources_interior_(_PREARG_LOCATION_ :) &
                     + fluxes(_ITERATOR_,:) / cell_thickness_ _INDEX_LOCATION_
               end do
#else
               sources_interior_(_PREARG_LOCATION_ :) = sources_interior_(_PREARG_LOCATION_ :) &
                  + fluxes(:) / cell_thickness_ _INDEX_LOCATION_
#endif
            end if
         _END_OUTER_HORIZONTAL_LOOP_
      end if

      if (bottom) then
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==0
#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         _VERTICAL_ITERATOR_ = model%p%domain%start(_FABM_DEPTH_DIMENSION_INDEX_)
#    else
         _VERTICAL_ITERATOR_ = model%p%domain%stop(_FABM_DEPTH_DIMENSION_INDEX_)
#    endif
#  elif !defined(_HORIZONTAL_IS_VECTORIZED_)
         _VERTICAL_ITERATOR_ = model%p%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#  endif
#endif

         _BEGIN_OUTER_HORIZONTAL_LOOP_
            fluxes = 0.0_rk
            call model%p%get_bottom_sources_rhs(_PREARG_HORIZONTAL_IN_ fluxes, sources_bottom_ _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_:_STOP_,:))
            if (size(model%p%interior_state_variables) > 0) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
               _DO_CONCURRENT_(_ITERATOR_,_START_,_STOP_)
#if _FABM_BOTTOM_INDEX_==-1
                  _VERTICAL_ITERATOR_ = model%p%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#endif
                  sources_interior_(_PREARG_LOCATION_ :) = sources_interior_(_PREARG_LOCATION_ :) &
                     + fluxes(_ITERATOR_,:) / cell_thickness_ _INDEX_LOCATION_
               end do
#else
#if _FABM_BOTTOM_INDEX_==-1
               _VERTICAL_ITERATOR_ = model%p%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#endif
               sources_interior_(_PREARG_LOCATION_ :) = sources_interior_(_PREARG_LOCATION_ :) &
                  + fluxes(:) / cell_thickness_ _INDEX_LOCATION_
#endif
            end if
         _END_OUTER_HORIZONTAL_LOOP_
      end if

      call model%p%finalize_outputs()
   end subroutine get_sources

   function check_state(pmodel, repair_) bind(c) result(valid_)
      !DIR$ ATTRIBUTES DLLEXPORT :: check_state
      type (c_ptr),         intent(in), value :: pmodel
      integer(c_int),value, intent(in)        :: repair_
      integer(c_int)                          :: valid_

      type (type_model_wrapper), pointer :: model
      logical :: repair, all_valid, interior_valid, surface_valid, bottom_valid
      _DECLARE_LOCATION_

#  if _FABM_DIMENSION_COUNT_ > 0
      integer :: _LOCATION_RANGE_
      istart__ = model%p%domain%start(1)
      istop__ = model%p%domain%stop(1)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 1
      jstart__ = model%p%domain%start(2)
      jstop__ = model%p%domain%stop(2)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 2
      kstart__ = model%p%domain%start(3)
      kstop__ = model%p%domain%stop(3)
#  endif

      call c_f_pointer(pmodel, model)
      if (model%p%status < status_start_done) then
         call driver%fatal_error('check_state', 'start has not been called yet.')
         return
      end if

      repair = int2logical(repair_)

      ! Check interior state everywhere
      all_valid = .true.
      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%p%check_interior_state(_PREARG_INTERIOR_IN_ repair, interior_valid)
         all_valid = all_valid .and. interior_valid
      _END_OUTER_INTERIOR_LOOP_

      ! Check surface state everywhere
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%p%check_surface_state(_PREARG_HORIZONTAL_IN_ repair, surface_valid)
         all_valid = all_valid .and. surface_valid
      _END_OUTER_HORIZONTAL_LOOP_

      ! Check bottom state everywhere
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%p%check_bottom_state(_PREARG_HORIZONTAL_IN_ repair, bottom_valid)
         all_valid = all_valid .and. bottom_valid
      _END_OUTER_HORIZONTAL_LOOP_

      valid_ = logical2int(all_valid)
   end function check_state

   function get_interior_diagnostic_data(pmodel, index) result(ptr) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_interior_diagnostic_data
      type (c_ptr),   intent(in), value :: pmodel
      integer(c_int), intent(in), value :: index
      type(c_ptr)                       :: ptr
      real(rk) _ATTRIBUTES_GLOBAL_, pointer :: pvalue

      type (type_model_wrapper), pointer :: model

      call c_f_pointer(pmodel, model)
      ptr = c_null_ptr
      pvalue => model%p%get_interior_diagnostic_data(index)
      if (associated(pvalue)) ptr = c_loc(pvalue)
   end function get_interior_diagnostic_data

   function get_horizontal_diagnostic_data(pmodel, index) result(ptr) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_horizontal_diagnostic_data
      type (c_ptr),   intent(in), value :: pmodel
      integer(c_int), intent(in), value :: index
      type(c_ptr)                       :: ptr
      real(rk) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: pvalue

      type (type_model_wrapper), pointer :: model

      call c_f_pointer(pmodel, model)
      ptr = c_null_ptr
      pvalue => model%p%get_horizontal_diagnostic_data(index)
      if (associated(pvalue)) ptr = c_loc(pvalue)
   end function get_horizontal_diagnostic_data

   subroutine finalize(model, keep_forced)
      type (type_model_wrapper), intent(inout) :: model
      logical, optional,         intent(in)    :: keep_forced

      logical :: keep_forced_

      keep_forced_ = .false.
      if (present(keep_forced)) keep_forced_ = keep_forced

      call model%p%finalize()
      call model%environment%finalize()
      if (.not. keep_forced) then
         call model%forced_parameters%finalize()
         call model%forced_couplings%finalize()
      end if
      deallocate(model%p)
   end subroutine finalize

   subroutine python_driver_fatal_error(self, location, message)
      class (type_python_driver), intent(inout) :: self
      character(len=*),           intent(in)    :: location, message

      if (error_occurred) return
      error_occurred = .true.
      error_message = trim(location) // ': ' // trim(message)
   end subroutine python_driver_fatal_error

   subroutine python_driver_log_message(self, message)
      class (type_python_driver), intent(inout) :: self
      character(len=*),           intent(in)    :: message

      !write (*,*) trim(message)
   end subroutine python_driver_log_message

end module fabm_c

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - Public License - www.gnu.org
!-----------------------------------------------------------------------

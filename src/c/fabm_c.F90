#include "fabm_driver.h"
#include "fabm_private.h"

module fabm_c

   use iso_c_binding, only: c_int, c_char, C_NULL_CHAR, c_f_pointer, c_loc, c_ptr, c_null_ptr, c_funptr, c_f_procpointer

   use fabm, only: type_fabm_model, type_fabm_variable, fabm_get_version, status_start_done, fabm_create_model
   use fabm_types, only: rke, attribute_length, type_model_list_node, type_base_model, &
                         factory, type_link, type_link_list, type_internal_variable, type_variable_list, type_variable_node, &
                         domain_interior, domain_horizontal, domain_scalar, get_free_unit, &
                         type_interior_standard_variable, type_horizontal_standard_variable
   use fabm_driver, only: type_base_driver, driver
   use fabm_c_variable, only: type_standard_variable_wrapper
   use fabm_python_helper
   use fabm_c_helper
   use yaml_settings, only: type_settings, type_key_value_pair, type_scalar_value, type_real_setting, type_integer_setting, &
      type_logical_setting, type_string_setting

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
      class (type_fabm_model),   pointer :: p => null()
      type (type_variable_list)          :: environment
      type (type_link_list)              :: coupling_link_list
   end type

   interface
      subroutine log_callback_interface(msg) bind(c)
         import c_char
         character(kind=c_char), intent(in) :: msg(*)
      end subroutine
   end interface
   procedure (log_callback_interface), pointer, save :: log_callback => null()

contains

   subroutine get_version(length, version_string) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_version
      integer(c_int), value, intent(in) :: length
      character(kind=c_char)            :: version_string(length)

      character(len=length-1) :: string

      ! Initialize driver object used by FABM for logging/error reporting.
      if (.not. associated(driver)) allocate(type_python_driver::driver)

      call fabm_get_version(string)
      call copy_to_c_string(string, version_string)
   end subroutine get_version

   subroutine set_log_callback(cb) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_log_callback
      type(c_funptr), intent(in), value :: cb
      call c_f_procpointer(cb, log_callback)
   end subroutine

   subroutine get_driver_settings(ndim, idepthdim, mask_type, variable_bottom_index) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_driver_settings
      integer(c_int), intent(out) :: ndim, idepthdim, mask_type, variable_bottom_index

      ndim = _FABM_DIMENSION_COUNT_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      ! We convert from 1-based Fortran index to 0-based C index, accounting for reversal of dimension order!
      idepthdim = _FABM_DIMENSION_COUNT_ - _FABM_DEPTH_DIMENSION_INDEX_
#else
      idepthdim = -1
#endif
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      mask_type = 1
#  else
      mask_type = 2
#  endif
#else
      mask_type = 0
#endif
#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_BOTTOM_INDEX_==-1
      variable_bottom_index = logical2int(.true.)
#else
      variable_bottom_index = logical2int(.false.)
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

      ! Initialize driver object used by FABM for logging/error reporting.
      if (.not. associated(driver)) allocate(type_python_driver::driver)

      ! Build FABM model tree (configuration will be read from file specified as argument).
      allocate(model)
      call c_f_pointer(c_loc(path), ppath)
      model%p => fabm_create_model(path=ppath(:index(ppath, C_NULL_CHAR) - 1))

      ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
      call model%p%set_domain(_PREARG_LOCATION_ 1._rke)

      ! Retrieve arrays to hold values for environmental variables and corresponding metadata.
      call get_environment_metadata(model%p, model%environment)

      call get_couplings(model%p, model%coupling_link_list)

      ptr = c_loc(model)
   end function create_model

   subroutine save_settings(pmodel, path, display) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: save_settings
      type (c_ptr), value,            intent(in) :: pmodel
      character(kind=c_char), target, intent(in) :: path(*)
      integer(c_int), value,          intent(in) :: display

      type (type_model_wrapper),       pointer :: model
      character(len=attribute_length), pointer :: ppath

      call c_f_pointer(pmodel, model)
      call c_f_pointer(c_loc(path), ppath)
      call model%p%settings%save(ppath(:index(ppath, C_NULL_CHAR) - 1), unit=get_free_unit(), display=display)
   end subroutine

#  if _FABM_DIMENSION_COUNT_ > 0
   subroutine set_domain_start(pmodel _POSTARG_LOCATION_) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_domain_start
      type (c_ptr),   intent(in), value  :: pmodel
      integer (kind=c_int), value, intent(in) :: _LOCATION_
      type (type_model_wrapper), pointer :: model
      call c_f_pointer(pmodel, model)
      call model%p%set_domain_start(_LOCATION_)
   end subroutine

   subroutine set_domain_stop(pmodel _POSTARG_LOCATION_) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_domain_stop
      type (c_ptr),   intent(in), value  :: pmodel
      integer (kind=c_int), value, intent(in) :: _LOCATION_
      type (type_model_wrapper), pointer :: model
      call c_f_pointer(pmodel, model)
      call model%p%set_domain_stop(_LOCATION_)
   end subroutine
#endif

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
   subroutine set_mask(pmodel, horizontal_mask_) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_mask
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), target :: horizontal_mask_(*)

      type (type_model_wrapper),                     pointer :: model
      integer(c_int) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: horizontal_mask

      call c_f_pointer(pmodel, model)
#if  _HORIZONTAL_DIMENSION_COUNT_ > 0
      call c_f_pointer(c_loc(horizontal_mask_), horizontal_mask, model%p%domain%horizontal_shape)
#else
      call c_f_pointer(c_loc(horizontal_mask_), horizontal_mask)
#endif
      call model%p%set_mask(horizontal_mask)
   end subroutine
#  else
   subroutine set_mask(pmodel, interior_mask_, horizontal_mask_) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_mask
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), target :: interior_mask_(*), horizontal_mask_(*)

      type (type_model_wrapper),                     pointer :: model
      integer(c_int) _ATTRIBUTES_GLOBAL_,            pointer :: interior_mask
      integer(c_int) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: horizontal_mask

      call c_f_pointer(pmodel, model)
#if  _HORIZONTAL_DIMENSION_COUNT_ > 0
      call c_f_pointer(c_loc(horizontal_mask_), horizontal_mask, model%p%domain%horizontal_shape)
#else
      call c_f_pointer(c_loc(horizontal_mask_), horizontal_mask)
#endif
#if  _FABM_DIMENSION_COUNT_ > 0
      call c_f_pointer(c_loc(interior_mask_), interior_mask, model%p%domain%shape)
#else
      call c_f_pointer(c_loc(interior_mask_), interior_mask)
#endif
      call model%p%set_mask(interior_mask, horizontal_mask)
   end subroutine
#  endif
#endif

#if _FABM_BOTTOM_INDEX_==-1
   subroutine set_bottom_index(pmodel, bottom_index_) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_bottom_index
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), target :: bottom_index_(*)

      type (type_model_wrapper),                     pointer :: model
      integer(c_int) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: bottom_index

      call c_f_pointer(pmodel, model)
#if  _HORIZONTAL_DIMENSION_COUNT_ > 0
      call c_f_pointer(c_loc(bottom_index_), bottom_index, model%p%domain%horizontal_shape)
#else
      call c_f_pointer(c_loc(bottom_index_), bottom_index)
#endif
      call model%p%set_bottom_index(bottom_index)
   end subroutine
#endif

   subroutine reinitialize(model)
      type (type_model_wrapper), intent(inout) :: model

      class (type_fabm_model),     pointer :: newmodel
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

      newmodel => fabm_create_model(settings=model%p%settings)

      ! Clean up old model
      call finalize(model)
      model%p => newmodel

      ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
      call model%p%set_domain(_PREARG_LOCATION_ 1._rke)

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
      nparameters = count_parameters(model%p%root)
      ncouplings = model%coupling_link_list%count()

   contains

      function count_parameters(root) result(n)
         class (type_base_model), intent(in) :: root
         integer                             :: n

         type (type_model_list_node), pointer :: instance
         type (type_key_value_pair),  pointer :: pair

         n = 0
         instance => root%children%first
         do while (associated(instance))
            if (instance%model%user_created) then
               pair => instance%model%parameters%first
               do while (associated(pair))
                  select type (scalar_value => pair%value)
                  class is (type_scalar_value)
                     n = n + 1
                  end select
                  pair => pair%next
               end do
            end if
            instance => instance%next
         end do
      end function

   end subroutine get_counts

   function get_parameter_by_index(root, i, name) result(scalar_value)
      class (type_base_model),    intent(inout) :: root
      integer,                    intent(in)    :: i
      character(len=*), optional, intent(out)   :: name
      class (type_scalar_value), pointer   :: scalar_value

      integer                              :: n
      type (type_model_list_node), pointer :: instance
      type (type_key_value_pair),  pointer :: pair

      n = 0
      instance => root%children%first
      do while (associated(instance))
         if (instance%model%user_created) then
            pair => instance%model%parameters%first
            do while (associated(pair))
               select type (value => pair%value)
               class is (type_scalar_value)
                  n = n + 1
                  if (n == i) then
                     scalar_value => value
                     if (present(name)) name = trim(instance%model%name) // '/' // trim(pair%name)
                     return
                  end if
               end select
               pair => pair%next
            end do
         end if
         instance => instance%next
      end do
   end function

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

   subroutine set_variable_save(pmodel, category, index, value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_variable_save
      type (c_ptr),           intent(in), value :: pmodel
      integer(c_int),         intent(in), value :: category, index, value

      type (type_model_wrapper),  pointer :: model

      call c_f_pointer(pmodel, model)

      select case (category)
      case (INTERIOR_DIAGNOSTIC_VARIABLE)
         model%p%interior_diagnostic_variables(index)%save = int2logical(value)
      case (HORIZONTAL_DIAGNOSTIC_VARIABLE)
         model%p%horizontal_diagnostic_variables(index)%save = int2logical(value)
      end select
   end subroutine set_variable_save

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
         variable => model%p%interior_diagnostic_variables(index)%original
      case (HORIZONTAL_DIAGNOSTIC_VARIABLE)
         variable => model%p%horizontal_diagnostic_variables(index)%original
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
      class (type_scalar_value), pointer :: scalar_value
      character(len=length)              :: name_

      call c_f_pointer(pmodel, model)

      scalar_value => get_parameter_by_index(model%p%root, index, name_)
      call copy_to_c_string(name_, name)
      if (allocated(scalar_value%units)) then
         call copy_to_c_string(scalar_value%units, units)
      else
         call copy_to_c_string('', units)
      end if
      call copy_to_c_string(scalar_value%long_name, long_name)
      select type (scalar_value)
      class is (type_real_setting)
         typecode = typecode_real
      class is (type_integer_setting)
         typecode = typecode_integer
      class is (type_logical_setting)
         typecode = typecode_logical
      class is (type_string_setting)
         typecode = typecode_string
      class default
         typecode = typecode_unknown
      end select
      has_default = logical2int(scalar_value%has_default)
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
      real(rke), target,         intent(in) :: ptr(*)
      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: pdata
#if _FABM_DIMENSION_COUNT_ > 0
      call c_f_pointer(c_loc(ptr), pdata, model%p%domain%shape)
#else
      call c_f_pointer(c_loc(ptr), pdata)
#endif
   end function

   function c_f_pointer_horizontal(model, ptr) result(pdata)
      type (type_model_wrapper), intent(in) :: model
      real(rke), target,         intent(in) :: ptr(*)
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: pdata
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
      real(rke),    intent(in), target :: dat(*)

      type (type_model_wrapper),          pointer :: model
      type (type_internal_variable),      pointer :: variable
      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: interior_data

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pvariable, variable)
      interior_data => c_f_pointer_interior(model, dat)
      call model%p%link_interior_data(variable, interior_data)
   end subroutine link_interior_data

   subroutine link_horizontal_data(pmodel, pvariable, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_horizontal_data
      type (c_ptr),   intent(in), value  :: pmodel
      type (c_ptr),   intent(in), value  :: pvariable
      real(rke),    intent(in), target :: dat(*)

      type (type_model_wrapper),                     pointer :: model
      type (type_internal_variable),                 pointer :: variable
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: horizontal_data

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pvariable, variable)
      horizontal_data => c_f_pointer_horizontal(model, dat)
      call model%p%link_horizontal_data(variable, horizontal_data)
   end subroutine link_horizontal_data

   subroutine link_scalar(pmodel, pvariable, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_scalar
      type (c_ptr),   intent(in), value  :: pmodel
      type (c_ptr),   intent(in), value  :: pvariable
      real(rke),    intent(in), target :: dat(*)

      type (type_model_wrapper),     pointer :: model
      type (type_internal_variable), pointer :: variable
      real(rke),                     pointer :: scalar_data

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pvariable, variable)
      call c_f_pointer(c_loc(dat), scalar_data)
      call model%p%link_scalar(variable, scalar_data)
   end subroutine link_scalar

   subroutine link_interior_state_data(pmodel, index, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_interior_state_data
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), value  :: index
      real(rke),      intent(in), target :: dat(*)

      type (type_model_wrapper),          pointer :: model
      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: dat_

      call c_f_pointer(pmodel, model)
      dat_ => c_f_pointer_interior(model, dat)
      dat_ = model%p%interior_state_variables(index)%initial_value
      call model%p%link_interior_state_data(index, dat_)
   end subroutine link_interior_state_data

   subroutine link_surface_state_data(pmodel, index, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_surface_state_data
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), value  :: index
      real(rke),      intent(in), target :: dat(*)

      type (type_model_wrapper),                     pointer :: model
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: dat_

      call c_f_pointer(pmodel, model)
      dat_ => c_f_pointer_horizontal(model, dat)
      dat_ = model%p%surface_state_variables(index)%initial_value
      call model%p%link_surface_state_data(index, dat_)
   end subroutine link_surface_state_data

   subroutine link_bottom_state_data(pmodel, index, dat) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_bottom_state_data
      type (c_ptr),   intent(in), value  :: pmodel
      integer(c_int), intent(in), value  :: index
      real(rke),      intent(in), target :: dat(*)

      type (type_model_wrapper),                     pointer :: model
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: dat_

      call c_f_pointer(pmodel, model)
      dat_ => c_f_pointer_horizontal(model, dat)
      dat_ = model%p%bottom_state_variables(index)%initial_value
      call model%p%link_bottom_state_data(index, dat_)
   end subroutine link_bottom_state_data

   subroutine get_sources(pmodel, t, sources_interior, sources_surface, sources_bottom, &
      do_surface, do_bottom, cell_thickness) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_sources
      type (c_ptr),   value,  intent(in) :: pmodel
      real(rke),      value,  intent(in) :: t
      real(rke),      target, intent(in) :: sources_interior(*), sources_surface(*), sources_bottom(*)
      integer(c_int), value,  intent(in) :: do_surface, do_bottom
      real(rke),      target, intent(in) :: cell_thickness(*)

      logical :: surface, bottom
      type (type_model_wrapper), pointer :: model
      real(rke) _DIMENSION_GLOBAL_PLUS_1_, pointer :: sources_interior_
      real(rke) _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_, pointer :: sources_surface_, sources_bottom_
      real(rke) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_, allocatable :: fluxes
      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: cell_thickness_
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
      i__ = model%p%domain%shape(1)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 1
      j__ = model%p%domain%shape(2)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 2
      k__ = model%p%domain%shape(3)
#  endif

      call c_f_pointer(c_loc(sources_interior), sources_interior_, &
         (/_PREARG_LOCATION_ size(model%p%interior_state_variables)/))
      call c_f_pointer(c_loc(sources_surface), sources_surface_, &
         (/_PREARG_HORIZONTAL_LOCATION_ size(model%p%surface_state_variables)/))
      call c_f_pointer(c_loc(sources_bottom), sources_bottom_, &
         (/_PREARG_HORIZONTAL_LOCATION_ size(model%p%bottom_state_variables)/))

      surface = int2logical(do_surface)
      bottom = int2logical(do_bottom)
      if ((surface .or. bottom) .and. size(model%p%interior_state_variables) > 0) &
         cell_thickness_ => c_f_pointer_interior(model, cell_thickness)

#  if _FABM_DIMENSION_COUNT_ > 0
      istart__ = model%p%domain%start(1)
      istop__ = model%p%domain%stop(1)
      i__ = istop__ - istart__ + 1
#  endif
#  if _FABM_DIMENSION_COUNT_ > 1
      jstart__ = model%p%domain%start(2)
      jstop__ = model%p%domain%stop(2)
      j__ = jstop__ - jstart__ + 1
#  endif
#  if _FABM_DIMENSION_COUNT_ > 2
      kstart__ = model%p%domain%start(3)
      kstop__ = model%p%domain%stop(3)
      k__ = kstop__ - kstart__ + 1
#  endif


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

      sources_interior_ = 0.0_rke
      sources_surface_ = 0.0_rke
      sources_bottom_ = 0.0_rke

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
            fluxes = 0.0_rke   ! get_surface_sources increments fluxes, so zero it first
            call model%p%get_surface_sources(_PREARG_HORIZONTAL_IN_ fluxes, sources_surface_ _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_:_STOP_,:))
            if (size(model%p%interior_state_variables) > 0) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
               _DO_CONCURRENT_(_ITERATOR_,_START_,_STOP_)
#endif

#ifdef _HAS_MASK_
               if (_IS_UNMASKED_(model%p%domain%mask_hz _INDEX_HORIZONTAL_LOCATION_)) then
#endif   

                  sources_interior_(_PREARG_LOCATION_ :) = sources_interior_(_PREARG_LOCATION_ :) &
#ifdef _HORIZONTAL_IS_VECTORIZED_
                  + fluxes(_ITERATOR_ - _START_ + 1,:) &
#else
                  + fluxes(:) &
#endif
                      / cell_thickness_ _INDEX_LOCATION_

#ifdef _HAS_MASK_
               end if   ! if unmasked
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
               end do   ! inner loop
#endif
            end if   ! if interior state variables
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
            fluxes = 0.0_rke   ! get_bottom_sources increments fluxes, so zero it first
            call model%p%get_bottom_sources_rhs(_PREARG_HORIZONTAL_IN_ fluxes, sources_bottom_ _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_:_STOP_,:))
            if (size(model%p%interior_state_variables) > 0) then
#ifdef _HORIZONTAL_IS_VECTORIZED_
               _DO_CONCURRENT_(_ITERATOR_,_START_,_STOP_)
#endif

#ifdef _HAS_MASK_
               if (_IS_UNMASKED_(model%p%domain%mask_hz _INDEX_HORIZONTAL_LOCATION_)) then
#endif

#if _FABM_BOTTOM_INDEX_==-1
                  _VERTICAL_ITERATOR_ = model%p%domain%bottom_indices _INDEX_HORIZONTAL_LOCATION_
#endif

           sources_interior_(_PREARG_LOCATION_ :) = sources_interior_(_PREARG_LOCATION_ :) &
#ifdef _HORIZONTAL_IS_VECTORIZED_
                     + fluxes(_ITERATOR_ - _START_ + 1,:) &
#else
                     + fluxes(:) &
#endif
                     / cell_thickness_ _INDEX_LOCATION_

#ifdef _HAS_MASK_
               end if   ! if unmasked
#endif

#ifdef _HORIZONTAL_IS_VECTORIZED_
               end do   ! inner loop
#endif
            end if   ! if interior state variables
         _END_OUTER_HORIZONTAL_LOOP_
      end if   ! if bottom

      call model%p%finalize_outputs()
   end subroutine get_sources

   subroutine get_vertical_movement(pmodel, velocity) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_vertical_movement
      type (c_ptr),   value,  intent(in) :: pmodel
      real(rke),      target, intent(in) :: velocity(*)

      type (type_model_wrapper), pointer :: model
      real(rke) _DIMENSION_GLOBAL_PLUS_1_, pointer :: velocity_
      _DECLARE_LOCATION_
#  if _FABM_DIMENSION_COUNT_ > 0
      integer :: _LOCATION_RANGE_
#  endif

      call c_f_pointer(pmodel, model)
      if (model%p%status < status_start_done) then
         call driver%fatal_error('get_vertical_movement', 'start has not been called yet.')
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

      call c_f_pointer(c_loc(velocity), velocity_, (/_PREARG_LOCATION_ size(model%p%interior_state_variables)/))

      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%p%get_vertical_movement(_PREARG_INTERIOR_IN_ velocity_ _INDEX_GLOBAL_INTERIOR_PLUS_1_(_START_:_STOP_,:))
      _END_OUTER_INTERIOR_LOOP_
   end subroutine get_vertical_movement

   subroutine get_conserved_quantities(pmodel, sums, cell_thickness) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_conserved_quantities
      type (c_ptr),   value,  intent(in) :: pmodel
      real(rke),      target, intent(in) :: sums(*)
      real(rke),      target, intent(in) :: cell_thickness(*)

      type (type_model_wrapper), pointer :: model
      real(rke) _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_, pointer :: sums_
      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: cell_thickness_
      real(rke) _DIMENSION_SLICE_PLUS_1_, allocatable :: sums_int
      integer :: ivar
      _DECLARE_LOCATION_
#  if _FABM_DIMENSION_COUNT_ > 0
      integer :: _LOCATION_RANGE_
#  endif

      call c_f_pointer(pmodel, model)
      if (model%p%status < status_start_done) then
         call driver%fatal_error('get_conserved_quantities', 'start has not been called yet.')
         return
      end if

#  if _FABM_DIMENSION_COUNT_ > 0
      i__ = model%p%domain%shape(1)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 1
      j__ = model%p%domain%shape(2)
#  endif
#  if _FABM_DIMENSION_COUNT_ > 2
      k__ = model%p%domain%shape(3)
#  endif

      cell_thickness_ => c_f_pointer_interior(model, cell_thickness)
      call c_f_pointer(c_loc(sums), sums_, &
         (/_PREARG_HORIZONTAL_LOCATION_ size(model%p%conserved_quantities)/))

#  if _FABM_DIMENSION_COUNT_ > 0
      istart__ = model%p%domain%start(1)
      istop__ = model%p%domain%stop(1)
      i__ = istop__ - istart__ + 1
#  endif
#  if _FABM_DIMENSION_COUNT_ > 1
      jstart__ = model%p%domain%start(2)
      jstop__ = model%p%domain%stop(2)
      j__ = jstop__ - jstart__ + 1
#  endif
#  if _FABM_DIMENSION_COUNT_ > 2
      kstart__ = model%p%domain%start(3)
      kstop__ = model%p%domain%stop(3)
      k__ = kstop__ - kstart__ + 1
#  endif

#ifdef _INTERIOR_IS_VECTORIZED_
      allocate(sums_int(_ITERATOR_, size(model%p%conserved_quantities)))
#else
      allocate(sums_int(size(model%p%conserved_quantities)))
#endif

      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%p%get_horizontal_conserved_quantities(_PREARG_HORIZONTAL_IN_ sums_ _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_START_:_STOP_,:))
      _END_OUTER_HORIZONTAL_LOOP_

      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%p%get_interior_conserved_quantities(_PREARG_INTERIOR_IN_ sums_int)
         _DO_CONCURRENT_(ivar,1,size(model%p%conserved_quantities))
#ifdef _INTERIOR_IS_VECTORIZED_
            _DO_CONCURRENT_(_ITERATOR_,_START_,_STOP_)
#endif

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
               if (_IS_UNMASKED_(model%p%domain%mask_hz _INDEX_HORIZONTAL_LOCATION_)) then
#  else
               if (_IS_UNMASKED_(model%p%domain%mask _INDEX_LOCATION_)) then
#  endif
#endif

                  sums_ _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_ITERATOR_,ivar) = sums_ _INDEX_GLOBAL_HORIZONTAL_PLUS_1_(_ITERATOR_,ivar) &
                     + cell_thickness_ _INDEX_LOCATION_ * &
#ifdef _INTERIOR_IS_VECTORIZED_
                     sums_int(_ITERATOR_ - _START_ + 1,ivar)
#else
                     sums_int(ivar)
#endif

#ifdef _HAS_MASK_
               end if   ! if not masked
#endif

#ifdef _INTERIOR_IS_VECTORIZED_
            end do   ! inner interior loop
#endif

         end do   ! state variable loop
      _END_OUTER_INTERIOR_LOOP_
   end subroutine get_conserved_quantities

   function check_state(pmodel, repair_) bind(c) result(valid_)
      !DIR$ ATTRIBUTES DLLEXPORT :: check_state
      type (c_ptr),         intent(in), value :: pmodel
      integer(c_int),value, intent(in)        :: repair_
      integer(c_int)                          :: valid_

      type (type_model_wrapper), pointer :: model
      logical :: repair, all_valid, interior_valid, surface_valid, bottom_valid
      _DECLARE_LOCATION_
#if _FABM_DIMENSION_COUNT_ > 0
      integer :: _LOCATION_RANGE_
#endif

      call c_f_pointer(pmodel, model)
      if (model%p%status < status_start_done) then
         call driver%fatal_error('check_state', 'start has not been called yet.')
         return
      end if

#if _FABM_DIMENSION_COUNT_ > 0
      istart__ = model%p%domain%start(1)
      istop__ = model%p%domain%stop(1)
#endif
#if _FABM_DIMENSION_COUNT_ > 1
      jstart__ = model%p%domain%start(2)
      jstop__ = model%p%domain%stop(2)
#endif
#if _FABM_DIMENSION_COUNT_ > 2
      kstart__ = model%p%domain%start(3)
      kstop__ = model%p%domain%stop(3)
#endif

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

      real(rke) _ATTRIBUTES_GLOBAL_, pointer :: pvalue

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

      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: pvalue

      type (type_model_wrapper), pointer :: model

      call c_f_pointer(pmodel, model)
      ptr = c_null_ptr
      pvalue => model%p%get_horizontal_diagnostic_data(index)
      if (associated(pvalue)) ptr = c_loc(pvalue)
   end function get_horizontal_diagnostic_data

   function get_standard_variable_data(pmodel, pstandard_variable, horizontal) result(ptr) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_standard_variable_data
      type (c_ptr),   intent(in), value :: pmodel, pstandard_variable
      integer(c_int), intent(out)       :: horizontal
      type(c_ptr)                       :: ptr

      type (type_model_wrapper),                pointer :: model
      type (type_standard_variable_wrapper),    pointer :: standard_variable
      real(rke) _ATTRIBUTES_GLOBAL_,            pointer :: pvalue
      real(rke) _ATTRIBUTES_GLOBAL_HORIZONTAL_, pointer :: pvalue_hz

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pstandard_variable, standard_variable)
      select type (p => standard_variable%p)
      class is (type_interior_standard_variable)
         pvalue => model%p%get_data(model%p%get_interior_variable_id(p))
         if (associated(pvalue)) ptr = c_loc(pvalue)
         horizontal = 0
      class is (type_horizontal_standard_variable)
         pvalue_hz => model%p%get_data(model%p%get_horizontal_variable_id(p))
         if (associated(pvalue_hz)) ptr = c_loc(pvalue_hz)
         horizontal = 1
      end select
   end function get_standard_variable_data

   subroutine require_data(pmodel, pstandard_variable) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: require_data
      type (c_ptr),   intent(in), value :: pmodel, pstandard_variable

      type (type_model_wrapper),             pointer :: model
      type (type_standard_variable_wrapper), pointer :: standard_variable

      call c_f_pointer(pmodel, model)
      call c_f_pointer(pstandard_variable, standard_variable)
      select type (p => standard_variable%p)
      class is (type_interior_standard_variable)
         call model%p%require_data(p)
      class is (type_horizontal_standard_variable)
         call model%p%require_data(p)
      end select
   end subroutine require_data

   subroutine finalize(model)
      type (type_model_wrapper), intent(inout) :: model

      call model%p%finalize()
      call model%environment%finalize()
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

      character(kind=c_char) :: cmessage(len(message) + 1)

      if (associated(log_callback)) then
         call copy_to_c_string(message, cmessage)
         call log_callback(cmessage)
      else
         write (*,*) trim(message)
      end if
   end subroutine python_driver_log_message

end module fabm_c

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - Public License - www.gnu.org
!-----------------------------------------------------------------------

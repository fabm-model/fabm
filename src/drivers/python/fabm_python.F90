#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: Python interface to the Framework for Aquatic Biogeochemical Models (FABM)
!
! !INTERFACE:
   module fabm_python
!
! !DESCRIPTION:
! TODO
!
! !USES:
   use iso_c_binding, only: c_double, c_int, c_char, C_NULL_CHAR, c_f_pointer, c_loc, c_ptr

   !DIR$ ATTRIBUTES DLLEXPORT :: STATE_VARIABLE,DIAGNOSTIC_VARIABLE,CONSERVED_QUANTITY

   use fabm
   use fabm_config
   use fabm_types, only:rk,attribute_length,type_model_list_node
   use fabm_driver, only: type_base_driver, driver
   use fabm_properties, only: type_property, type_property_dictionary, type_real_property
   use fabm_python_helper

   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
   public

   integer,parameter :: BULK_STATE_VARIABLE            = 1
   integer,parameter :: SURFACE_STATE_VARIABLE         = 2
   integer,parameter :: BOTTOM_STATE_VARIABLE          = 3
   integer,parameter :: BULK_DIAGNOSTIC_VARIABLE       = 4
   integer,parameter :: HORIZONTAL_DIAGNOSTIC_VARIABLE = 5
   integer,parameter :: CONSERVED_QUANTITY             = 6

   type (type_model),private,target,save :: model
   real(8),dimension(:),pointer :: state
   character(len=1024),dimension(:),allocatable :: environment_names,environment_units
   character(len=1024) :: yaml_path

   type (type_property_dictionary),save,private :: forced_parameters

   type,extends(type_base_driver) :: type_python_driver
   contains
      procedure :: fatal_error => python_driver_fatal_error
      procedure :: log_message => python_driver_log_message
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the model
!
! !INTERFACE:
   subroutine initialize(path) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: initialize
      character(kind=c_char),target,intent(in) :: path(*)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
      character(len=attribute_length),pointer :: ppath
!-----------------------------------------------------------------------
!BOC
      call c_f_pointer(c_loc(path), ppath)
      yaml_path = ppath(:index(ppath,C_NULL_CHAR)-1)

      if (model%initialized) call finalize()

      if (.not.associated(driver)) allocate(type_python_driver::driver)

      ! Build FABM model tree (configuration will be read from fabm.yaml).
      call fabm_create_model_from_yaml_file(model,path=yaml_path,parameters=forced_parameters)

      ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
      call fabm_set_domain(model)

      ! Retrieve arrays to hold values for environmental variables and corresponding metadata.
      call get_environment_metadata(model,environment_names,environment_units)

   end subroutine initialize
!EOC

   subroutine check_ready()
      !DIR$ ATTRIBUTES DLLEXPORT :: check_ready
      call fabm_check_ready(model)
   end subroutine

   integer(c_int) function model_count() bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: model_count

      type (type_model_list_node), pointer :: node

      model_count = 0
      node => model%root%children%first
      do while (associated(node))
         model_count = model_count + 1
         node => node%next
      end do
   end function

   subroutine get_variable_counts(nstate_bulk,nstate_surface,nstate_bottom,ndiagnostic_bulk,ndiagnostic_horizontal,nconserved, &
      ndependencies,nparameters) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_variable_counts
      integer(c_int),intent(out) :: nstate_bulk,nstate_surface,nstate_bottom
      integer(c_int),intent(out) :: ndiagnostic_bulk,ndiagnostic_horizontal
      integer(c_int),intent(out) :: nconserved,ndependencies,nparameters
      nstate_bulk = size(model%state_variables)
      nstate_surface = size(model%surface_state_variables)
      nstate_bottom = size(model%bottom_state_variables)
      ndiagnostic_bulk = size(model%diagnostic_variables)
      ndiagnostic_horizontal = size(model%horizontal_diagnostic_variables)
      nconserved = size(model%conserved_quantities)
      ndependencies = size(environment_names)
      nparameters = model%root%parameters%size()
   end subroutine

   subroutine get_variable_metadata(category,index,length,name,units,long_name,path) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_variable_metadata
      integer(c_int),intent(in), value             :: category,index,length
      character(kind=c_char),intent(out),dimension(length) :: name,units,long_name,path

      class (type_external_variable),pointer :: variable

      ! Get a pointer to the target variable
      select case (category)
      case (BULK_STATE_VARIABLE)
         variable => model%state_variables(index)
      case (SURFACE_STATE_VARIABLE)
         variable => model%surface_state_variables(index)
      case (BOTTOM_STATE_VARIABLE)
         variable => model%bottom_state_variables(index)
      case (BULK_DIAGNOSTIC_VARIABLE)
         variable => model%diagnostic_variables(index)
      case (HORIZONTAL_DIAGNOSTIC_VARIABLE)
         variable => model%horizontal_diagnostic_variables(index)
      case (CONSERVED_QUANTITY)
         variable => model%conserved_quantities(index)
      end select
      call copy_to_c_string(variable%name,     name)
      call copy_to_c_string(variable%units,    units)
      call copy_to_c_string(variable%long_name,long_name)
      call copy_to_c_string(variable%path,     path)
   end subroutine

   subroutine get_parameter_metadata(index,length,name,units,long_name,typecode) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_parameter_metadata
      integer(c_int),intent(in), value             :: index,length
      character(kind=c_char),intent(out),dimension(length) :: name,units,long_name
      integer(c_int),intent(out)                   :: typecode

      integer                       :: i
      class (type_property),pointer :: property

      i = 1
      property => model%root%parameters%first
      do while (associated(property))
         if (index==i) exit
         property => property%next
         i = i + 1
      end do

      call copy_to_c_string(property%name,     name)
      call copy_to_c_string(property%units,    units)
      call copy_to_c_string(property%long_name,long_name)
      typecode = property%typecode()
   end subroutine

   subroutine get_dependency_metadata(index,length,name,units) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_dependency_metadata
      integer(c_int),intent(in), value             :: index,length
      character(kind=c_char),intent(out),dimension(length) :: name,units

      call copy_to_c_string(environment_names(index),name)
      call copy_to_c_string(environment_units(index),units)
   end subroutine

   subroutine link_dependency_data(index,value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_dependency_data
      integer(c_int),intent(in),value  :: index
      real(c_double),intent(in),target :: value

      call fabm_link_bulk_data(model,environment_names(index),value)
      call fabm_link_horizontal_data(model,environment_names(index),value)
      call fabm_link_scalar_data(model,environment_names(index),value)
   end subroutine

   subroutine link_bulk_state_data(index,value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_bulk_state_data
      integer(c_int),intent(in),   value  :: index
      real(c_double),intent(inout),target :: value

      value = model%state_variables(index)%initial_value
      call fabm_link_bulk_state_data(model,index,value)
   end subroutine

   subroutine link_surface_state_data(index,value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_surface_state_data
      integer(c_int),intent(in),   value  :: index
      real(c_double),intent(inout),target :: value

      value = model%surface_state_variables(index)%initial_value
      call fabm_link_surface_state_data(model,index,value)
   end subroutine
   
   subroutine link_bottom_state_data(index,value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: link_bottom_state_data
      integer(c_int),intent(in),   value  :: index
      real(c_double),intent(inout),target :: value

      value = model%bottom_state_variables(index)%initial_value
      call fabm_link_bottom_state_data(model,index,value)
   end subroutine

   subroutine get_rates(pelagic_rates_) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_rates
      real(c_double),target,intent(in) :: pelagic_rates_(*)

      real(c_double),pointer :: pelagic_rates(:)

      call c_f_pointer(c_loc(pelagic_rates_),pelagic_rates, &
        (/size(model%state_variables)+size(model%surface_state_variables)+size(model%bottom_state_variables)/))
      pelagic_rates = 0.0_rk
      call fabm_do(model,pelagic_rates)

      ! Compute rate of change in conserved quantities
      !call fabm_state_to_conserved_quantities(model,pelagic_rates,conserved_rates)

      ! Normalize rate of change in conserved quantities to sum of absolute rates of change.
      !call fabm_state_to_conserved_quantities(model,abs(pelagic_rates),abs_conserved_rates)
      !where (abs_conserved_rates>0.0_rk) conserved_rates = conserved_rates/abs_conserved_rates
   end subroutine

   subroutine get_bulk_diagnostic_data(index,ptr) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_bulk_diagnostic_data
      integer(c_int),intent(in),value :: index
      type(c_ptr),   intent(out)      :: ptr
      ptr = c_loc(fabm_get_bulk_diagnostic_data(model,index))
   end subroutine

   subroutine get_horizontal_diagnostic_data(index,ptr) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_horizontal_diagnostic_data
      integer(c_int),intent(in),value :: index
      type(c_ptr),   intent(out)      :: ptr
      ptr = c_loc(fabm_get_horizontal_diagnostic_data(model,index))
   end subroutine

   subroutine finalize() bind(c)
      call fabm_finalize(model)
      if (allocated(environment_names)) deallocate(environment_names)
      if (allocated(environment_units)) deallocate(environment_units)
   end subroutine

   subroutine reset_parameter(name) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: reset_parameter
      character(kind=c_char),target,intent(in) :: name(*)

      character(len=attribute_length),pointer :: pname

      call c_f_pointer(c_loc(name), pname)
      call forced_parameters%delete(pname(:index(pname,C_NULL_CHAR)-1))

      ! Re-initialize the model using updated parameter values
      call initialize(yaml_path)
   end subroutine

   subroutine set_real_parameter(name,value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_real_parameter
      character(kind=c_char),target,    intent(in) :: name(*)
      real(c_double),value,intent(in) :: value

      character(len=attribute_length),pointer :: pname

      call c_f_pointer(c_loc(name), pname)
      call forced_parameters%set_real(pname(:index(pname,C_NULL_CHAR)-1),value)

      ! Re-initialize the model using updated parameter values
      call initialize(yaml_path)
   end subroutine

   function get_real_parameter(name,default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_real_parameter
      character(kind=c_char),target,    intent(in) :: name(*)
      real(c_double),value,intent(in) :: default
      real(c_double)                  :: value
      character(len=attribute_length),pointer :: pname
      call c_f_pointer(c_loc(name), pname)
      value = model%root%parameters%get_real(pname(:index(pname,C_NULL_CHAR)-1),default)
   end function

   subroutine set_integer_parameter(name,value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_integer_parameter
      character(kind=c_char),target,    intent(in) :: name(*)
      integer(c_int),value,intent(in) :: value

      character(len=attribute_length),pointer :: pname

      call c_f_pointer(c_loc(name), pname)
      call forced_parameters%set_integer(pname(:index(pname,C_NULL_CHAR)-1),value)

      ! Re-initialize the model using updated parameter values
      call initialize(yaml_path)
   end subroutine

   function get_integer_parameter(name,default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_integer_parameter
      character(kind=c_char),target,    intent(in) :: name(*)
      integer(c_int),value,intent(in) :: default
      integer(c_int)                  :: value
      character(len=attribute_length),pointer :: pname
      call c_f_pointer(c_loc(name), pname)
      value = model%root%parameters%get_integer(pname(:index(pname,C_NULL_CHAR)-1),default)
   end function

   subroutine set_logical_parameter(name,value) bind(c)
      !DIR$ ATTRIBUTES DLLEXPORT :: set_logical_parameter
      character(kind=c_char),target,    intent(in) :: name(*)
      integer(c_int),value,intent(in) :: value

      character(len=attribute_length),pointer :: pname

      call c_f_pointer(c_loc(name), pname)
      call forced_parameters%set_logical(pname(:index(pname,C_NULL_CHAR)-1),value/=0)

      ! Re-initialize the model using updated parameter values
      call initialize(yaml_path)
   end subroutine

   function get_logical_parameter(name,default) bind(c) result(value)
      !DIR$ ATTRIBUTES DLLEXPORT :: get_logical_parameter
      character(kind=c_char),target,    intent(in) :: name(*)
      integer(c_int),value,intent(in) :: default
      integer(c_int)                  :: value
      character(len=attribute_length),pointer :: pname
      call c_f_pointer(c_loc(name), pname)
      if (model%root%parameters%get_logical(pname(:index(pname,C_NULL_CHAR)-1),default/=0)) then
         value = 1
      else
         value = 0
      end if
   end function

   subroutine set_string_parameter(name,value)
      character(len=*),intent(in) :: name,value
      call forced_parameters%set_string(name,value)

      ! Re-initialize the model using updated parameter values
      call initialize(yaml_path)
   end subroutine

   function get_string_parameter(name,default) result(value)
      character,target,intent(in) :: name(*),default(*)
      character            :: value(1024)
      class (type_property),pointer :: property
      character(len=attribute_length),pointer :: pname,pdefault
      call c_f_pointer(c_loc(name), pname)
      call c_f_pointer(c_loc(default), pdefault)
      value = pdefault(:index(pdefault,C_NULL_CHAR)-1)
      property => model%root%parameters%get_property(pname(:index(pname,C_NULL_CHAR)-1))
      value = property%to_string()
   end function

   subroutine python_driver_fatal_error(self,location,message)
      class (type_python_driver),intent(inout) :: self
      character(len=*),          intent(in)    :: location,message

      write (*,*) trim(location)//': '//trim(message)
      stop 1
   end subroutine

   subroutine python_driver_log_message(self,message)
      class (type_python_driver),intent(inout) :: self
      character(len=*),          intent(in)    :: message

      !write (*,*) trim(message)
   end subroutine

   subroutine copy_to_c_string(string,cstring)
      character(len=*),intent(in)  :: string
      character,       intent(out) :: cstring(:)
      integer i,n
      n = min(len_trim(string),size(cstring)-1)
      do i=1,n
         cstring(i) = string(i:i)
      end do
      cstring(n+1) = C_NULL_CHAR 
   end subroutine

   end module fabm_python

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

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
   use fabm
   use fabm_config
   use fabm_types, only:rk,attribute_length
   use fabm_driver, only: type_base_driver, driver
   use fabm_properties, only: type_property, type_property_dictionary
   use fabm_python_helper

   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
   public

   type (type_model),private,target,save :: model
   real(8),dimension(:),allocatable :: state,environment
   character(len=1024),dimension(:),allocatable :: state_names,state_units
   character(len=1024),dimension(:),allocatable :: environment_names,environment_units
   character(len=1024),dimension(:),allocatable :: conserved_quantity_names,conserved_quantity_units
   character(len=1024),dimension(:),allocatable :: parameter_long_names,parameter_names,parameter_units
   integer,            dimension(:),allocatable :: parameter_types

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
   subroutine initialize()
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman
!
!EOP
!
! !LOCAL VARIABLES:
      integer           :: i,n
      class (type_property),pointer :: property
!
!-----------------------------------------------------------------------
!BOC
      if (model%initialized) call finalize()

      if (.not.associated(driver)) allocate(type_python_driver::driver)

      ! Build FABM model tree (configuration will be read from fabm.yaml).
      call fabm_create_model_from_yaml_file(model,parameters=forced_parameters)

      ! Send information on spatial domain to FABM (this also allocates memory for diagnostics)
      call fabm_set_domain(model)

      ! Create arrays to hold state variable values, use the initial values specified by the model,
      ! link state data to FABM.
      n = size(model%state_variables) + size(model%bottom_state_variables) + size(model%surface_state_variables)
      allocate(state(n))
      allocate(state_names(n))
      allocate(state_units(n))
      n = 0
      do i=1,size(model%state_variables)
         n = n + 1
         state(n) = model%state_variables(i)%initial_value
         state_names(n) = model%state_variables(i)%name
         state_units(n) = model%state_variables(i)%units
         call fabm_link_bulk_state_data(model,i,state(n))
      end do
      do i=1,size(model%bottom_state_variables)
         n = n + 1
         state(n) = model%bottom_state_variables(i)%initial_value
         state_names(n) = model%bottom_state_variables(i)%name
         state_units(n) = model%bottom_state_variables(i)%units
         call fabm_link_bulk_state_data(model,i,state(n))
      end do
      do i=1,size(model%surface_state_variables)
         n = n + 1
         state(n) = model%surface_state_variables(i)%initial_value
         state_names(n) = model%surface_state_variables(i)%name
         state_units(n) = model%surface_state_variables(i)%units
         call fabm_link_bulk_state_data(model,i,state(n))
      end do

      allocate(conserved_quantity_names(size(model%conserved_quantities)))
      allocate(conserved_quantity_units(size(model%conserved_quantities)))
      conserved_quantity_names = model%conserved_quantities(:)%name
      conserved_quantity_units = model%conserved_quantities(:)%units

      ! Create arrays for parameter names, units, data types.
      call model%root%parameters%keys(parameter_names)
      allocate(parameter_types(size(parameter_names)))
      allocate(parameter_units(size(parameter_names)))
      allocate(parameter_long_names(size(parameter_names)))
      do i = 1,size(parameter_names)
         property => model%root%parameters%get_property(parameter_names(i))
         parameter_types(i) = property%typecode()
         parameter_units(i) = property%units
         parameter_long_names(i) = property%long_name
      end do

      ! Retrieve arrays to hold values for environmental variables and corresponding metadata.
      call get_environment(model,environment_names,environment_units,environment)

      call fabm_check_ready(model)

   end subroutine initialize
!EOC

   subroutine get_rates(n,m,pelagic_rates,conserved_rates)
      integer, intent(in)  :: n,m
      real(8),intent(inout) :: pelagic_rates(n),conserved_rates(m)

      real(8) :: abs_conserved_rates(m)

      pelagic_rates = 0.0_rk
      call fabm_do(model,pelagic_rates)

      ! Compute rate of change in conserved quantities
      call fabm_state_to_conserved_quantities(model,pelagic_rates,conserved_rates)

      ! Normalize rate of change in conserved quantities to sum of absolute rates of change.
      call fabm_state_to_conserved_quantities(model,abs(pelagic_rates),abs_conserved_rates)
      where (abs_conserved_rates>0.0_rk) conserved_rates = conserved_rates/abs_conserved_rates
   end subroutine

   subroutine finalize()
      call fabm_finalize(model)
      if (allocated(state)) deallocate(state)
      if (allocated(state_names)) deallocate(state_names)
      if (allocated(state_units))  deallocate(state_units)
      if (allocated(environment)) deallocate(environment)
      if (allocated(environment_names)) deallocate(environment_names)
      if (allocated(environment_units)) deallocate(environment_units)
      if (allocated(parameter_names)) deallocate(parameter_names)
      if (allocated(parameter_long_names)) deallocate(parameter_long_names)
      if (allocated(parameter_types)) deallocate(parameter_types)
      if (allocated(parameter_units)) deallocate(parameter_units)
      if (allocated(conserved_quantity_names)) deallocate(conserved_quantity_names)
      if (allocated(conserved_quantity_units)) deallocate(conserved_quantity_units)
   end subroutine

   subroutine reset_parameter(name)
      character(len=*),intent(in) :: name
      call forced_parameters%delete(name)

      ! Re-initialize the model using updated parameter values
      call initialize()
   end subroutine

   subroutine set_real_parameter(name,value)
      character(len=*),intent(in) :: name
      real(8),         intent(in) :: value
      call forced_parameters%set_real(name,value)

      ! Re-initialize the model using updated parameter values
      call initialize()
   end subroutine

   function get_real_parameter(name,default) result(value)
      character(len=*),intent(in) :: name
      real(8),         intent(in) :: default
      real(8)                     :: value
      class (type_property),pointer :: property
      value = default
      property => model%root%parameters%get_property(name)
      value = property%to_real()
   end function

   subroutine set_integer_parameter(name,value)
      character(len=*),intent(in) :: name
      integer,         intent(in) :: value
      call forced_parameters%set_integer(name,value)

      ! Re-initialize the model using updated parameter values
      call initialize()
   end subroutine

   function get_integer_parameter(name,default) result(value)
      character(len=*),intent(in) :: name
      integer,         intent(in) :: default
      integer                     :: value
      class (type_property),pointer :: property
      value = default
      property => model%root%parameters%get_property(name)
      value = property%to_integer()
   end function

   subroutine set_logical_parameter(name,value)
      character(len=*),intent(in) :: name
      logical,         intent(in) :: value
      call forced_parameters%set_logical(name,value)

      ! Re-initialize the model using updated parameter values
      call initialize()
   end subroutine

   function get_logical_parameter(name,default) result(value)
      character(len=*),intent(in) :: name
      logical,         intent(in) :: default
      logical                     :: value
      class (type_property),pointer :: property
      value = default
      property => model%root%parameters%get_property(name)
      value = property%to_logical()
   end function

   subroutine set_string_parameter(name,value)
      character(len=*),intent(in) :: name,value
      call forced_parameters%set_string(name,value)

      ! Re-initialize the model using updated parameter values
      call initialize()
   end subroutine

   function get_string_parameter(name,default) result(value)
      character(len=*),intent(in) :: name,default
      character(len=1024)         :: value
      class (type_property),pointer :: property
      value = default
      property => model%root%parameters%get_property(name)
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

   end module fabm_python

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

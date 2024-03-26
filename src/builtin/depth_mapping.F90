#include "fabm_driver.h"

module fabm_builtin_depth_mapping
   use fabm_types
   use fabm_builtin_depth_integral, only: type_depth_integral
   use fabm_particle

   implicit none

   private

   public type_vertical_distribution, type_vertical_depth_range
   public set_vertical_distribution, register_depth_explicit_dependency, register_depth_explicit_state_dependency
   public type_depth_integrated_particle
   public type_depth_integrated_particle_override

   type, extends(type_base_model) :: type_vertical_distribution
      type (type_diagnostic_variable_id) :: id_w
   end type

   type, extends(type_vertical_distribution) :: type_vertical_depth_range
      type (type_dependency_id) :: id_z, id_h
      real(rk) :: minimum_depth = 0.0_rk
      real(rk) :: maximum_depth = 1e5_rk
   contains
      procedure :: initialize => vertical_depth_range_initialize
      procedure :: do         => vertical_depth_range_do
   end type

   type, extends(type_base_model), public :: type_weighted_depth_integral
      type (type_dependency_id)                     :: id_source  ! depth-explicit variable to integrate
      type (type_dependency_id)                     :: id_h       ! layer thicknesses (m)
      type (type_dependency_id)                     :: id_w       ! vertical distribution weights
      type (type_horizontal_dependency_id)          :: id_w_int   ! depth-integrated vertical distribution weights (for averaging)
      type (type_horizontal_diagnostic_variable_id) :: id_result  ! result (depth integral or average)

      logical :: average               = .false.
      logical :: act_as_state_variable = .false.
      logical :: proportional_change   = .false.
      integer :: domain = domain_bottom
   contains
      procedure :: initialize => weighted_depth_integral_initialize
      procedure :: do_column  => weighted_depth_integral_do_column
   end type

   type, extends(type_base_model), public :: type_projector
      type (type_horizontal_dependency_id) :: id_source
      type (type_dependency_id)            :: id_w
      type (type_horizontal_dependency_id) :: id_w_int
      type (type_diagnostic_variable_id)   :: id_result

      logical :: act_as_state_variable = .false.
   contains
      procedure :: initialize => projector_initialize
      procedure :: do         => projector_do
   end type

   ! This model type takes a depth-integrated source term and distributes it over a pelagic state variable.
   ! The fraction injected into each layer scales with the vertical distribution weights
   ! (as well as layer thickness)
   type, extends(type_base_model) :: type_absolute_rate_distributor
      type (type_horizontal_dependency_id) :: id_int        ! Depth-integrated variable for which sources-minus-sinks are specified
      type (type_horizontal_dependency_id) :: id_sms_int    ! Depth-integrated sources-minus-sinks
      type (type_state_variable_id)        :: id_target     ! Depth-explicit variable that should absorb the sources-minus-sinks
      type (type_dependency_id)            :: id_w          ! Weights for the vertical distribution of the sinks and sources
      type (type_horizontal_dependency_id) :: id_w_int      ! Depth integral of weights
   contains
      procedure :: initialize => absolute_rate_distributor_initialize
      procedure :: do         => absolute_rate_distributor_do
   end type

   ! This model type takes a depth-integrated source term and distributes it over a pelagic state variable.
   ! The fraction injected into each layer scales with the vertical distribution weights and
   ! the local concentration of the targeted pelagic variable (as well as layer thickness)
   type, extends(type_absolute_rate_distributor) :: type_relative_rate_distributor
   contains
      procedure :: do         => relative_rate_distributor_do
   end type

   ! This model type takes a depth-explicit source term and distributes it over a depth-integrated state variable.
   type, extends(type_base_model) :: type_projected_rate_distributor
      type (type_dependency_id)            :: id_sms
      type (type_dependency_id)            :: id_h
      type (type_bottom_state_variable_id) :: id_target
   contains
      procedure :: initialize => projected_rate_distributor_initialize
      procedure :: do_column  => projected_rate_distributor_do_column
   end type

   ! This model type inherits from the base particle model and adds several procedures
   ! for coupling to depth-integrate and depth-averaged pelagic variables.
   type, extends(type_particle_model) :: type_depth_integrated_particle
      type (type_horizontal_dependency_id) :: id_w_int
   contains
      procedure :: initialize => depth_integrated_particle_initialize
      procedure :: register_mapped_model_dependency
      procedure :: request_mapped_coupling1 => depth_integrated_particle_request_mapped_coupling1
      procedure :: request_mapped_coupling2 => depth_integrated_particle_request_mapped_coupling2
      procedure :: request_mapped_coupling3 => depth_integrated_particle_request_mapped_coupling3
      procedure :: request_mapped_coupling4 => depth_integrated_particle_request_mapped_coupling4
      generic :: request_mapped_coupling => request_mapped_coupling1, request_mapped_coupling2, &
         request_mapped_coupling3, request_mapped_coupling4
      procedure :: request_mapped_coupling_to_model1 => depth_integrated_particle_request_mapped_coupling_to_model1
      procedure :: request_mapped_coupling_to_model2 => depth_integrated_particle_request_mapped_coupling_to_model2
      procedure :: request_mapped_coupling_to_model3 => depth_integrated_particle_request_mapped_coupling_to_model3
      procedure :: request_mapped_coupling_to_model4 => depth_integrated_particle_request_mapped_coupling_to_model4
      generic :: request_mapped_coupling_to_model => request_mapped_coupling_to_model1, request_mapped_coupling_to_model2, &
         request_mapped_coupling_to_model3, request_mapped_coupling_to_model4
   end type

   type, extends(type_depth_integrated_particle) :: type_depth_integrated_particle_override
      integer :: domain = domain_bottom
   contains
      procedure :: initialize                  => depth_integrated_particle_override_initialize
      procedure :: complete_internal_variables => depth_integrated_particle_override_complete_internal_variables
   end type

   ! This model type couples to a depth-explicit particle and creates depth-integrated equivalents
   ! for each of its state variables, which then accept depth-integrated source terms.
   ! These are redistributed over the original depth-explicit variables.
   ! Integration and redistribution take into account a set of weights that can be used
   ! to specify the vertical distribution/habitat of the model creating the particle integrator.
   type, extends(type_particle_model) :: type_particle_integrator
      integer :: domain = domain_bottom
      logical :: proportional_change = .false.
   contains
      procedure :: initialize                  => particle_integrator_initialize
      procedure :: complete_internal_variables => particle_integrator_complete_internal_variables
   end type

   ! This model type couples to a depth-integrated particle and creates depth-explicit equivalents
   ! for each of its state variables, which then accept depth-explicit source terms.
   ! These source terms are depth-integrated and applied to the original depth-integrated variables.
   ! Projection of depth-integrated variables takes into account a set of weights that can be used
   ! to specify the vertical distribution/habitat.
   type, extends(type_particle_model) :: type_projected_particle
   contains
      procedure :: initialize                  => projected_particle_initialize
      procedure :: complete_internal_variables => projected_particle_complete_internal_variables
   end type

contains

   subroutine vertical_depth_range_initialize(self, configunit)
      class (type_vertical_depth_range), intent(inout), target :: self
      integer,                           intent(in)            :: configunit

      call self%register_implemented_routines((/source_do/))
      call self%get_parameter(self%minimum_depth, 'minimum_depth', 'm', 'minimum depth (distance from surface)', &
         default=self%minimum_depth)
      call self%get_parameter(self%maximum_depth, 'maximum_depth', 'm', 'maximum depth (distance from surface)', &
         default=self%maximum_depth)
      call self%register_diagnostic_variable(self%id_w, 'w', '1', 'weights')
      call self%register_dependency(self%id_z, standard_variables%depth)
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)
   end subroutine

   subroutine vertical_depth_range_do(self, _ARGUMENTS_DO_)
      class (type_vertical_depth_range), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real :: z_center, h, z_top, z_bottom, h_overlap

      _LOOP_BEGIN_
         _GET_(self%id_z, z_center)
         _GET_(self%id_h, h)
         z_top = z_center - 0.5_rk * h
         z_bottom = z_center + 0.5_rk * h
         h_overlap = max(0.0_rk, min(z_bottom, self%maximum_depth) - max(z_top, self%minimum_depth))
         _SET_DIAGNOSTIC_(self%id_w, h_overlap / h)
      _LOOP_END_
   end subroutine

   subroutine set_vertical_distribution(self, w_calculator, w_int_link)
      class (type_base_model),                      intent(inout) :: self
      class (type_vertical_distribution), optional, intent(inout) :: w_calculator
      type (type_link), pointer,          optional                :: w_int_link

      class (type_depth_integral), pointer :: w_int_calculator
      type (type_link), pointer            :: w_int_link_

      ! Register a dependency on the vertical distribution weights
      ! We use add_interior_variable instead of register_dependency because we do not
      ! need read access to the weights ourselves (they are only for child models)
      call self%add_interior_variable('w', '1', 'vertical distribution weights')

      if (present(w_calculator)) then
         call self%add_child(w_calculator, 'w_calculator')
         call self%request_coupling('w', 'w_calculator/w')
      end if

      ! Create a child model that calculates the depth integral of the weights
      ! This is used to calculate depth averages form depth integrals, among others.
      allocate(w_int_calculator)
      call self%add_child(w_int_calculator, 'w_int_calculator')
      call w_int_calculator%request_coupling('source', '../w')

      ! Make the depth integrated available as a local variable
      if (present(w_int_link)) then
         ! A variable has already been registered to hold the depth-integrated weights
         w_int_link_ => w_int_link
      else
         ! Create a local alias for the depth-integrated weights
         w_int_link_ => null()
         call self%add_horizontal_variable('w_int', 'm', 'integral of vertical distribution weights', link=w_int_link_)
      end if
      call self%request_coupling(w_int_link_, 'w_int_calculator/result')
   end subroutine

   subroutine register_depth_explicit_dependency(self, id, name, units, long_name, average, link)
      class (type_base_model),               intent(inout) :: self
      class (type_horizontal_dependency_id), intent(inout) :: id
      character(len=*),                      intent(in)    :: name, units, long_name
      logical, optional,                     intent(in)    :: average
      type (type_link), optional, pointer :: link

      class (type_weighted_depth_integral), pointer :: depth_integral
      character(len=4)                              :: postfix

      ! Create a placeholder for the depth-explicit source variable
      if (present(link)) link => null()
      call self%add_interior_variable(name, units, long_name, link=link)

      ! Determine postfix of the name of the dept-integrated/averaged variable
      postfix = '_int'
      if (present(average)) postfix = '_ave'

      ! Create a child model that will depth-integrate or depth-average the source variable
      ! Couple this to the vertical distribution weights and our placeholder for the source variable
      allocate(depth_integral)
      if (present(average)) depth_integral%average = average
      call self%add_child(depth_integral, trim(name) // postfix // '_calculator')
      call depth_integral%request_coupling('source', '../' // trim(name))
      call depth_integral%request_coupling('w', '../w')
      if (depth_integral%average) call depth_integral%request_coupling('w_int', '../w_int')

      ! If the provided identifer has not yet been registered, do so first
      ! Then couple the identifier to the result of the depth integration
      if (.not. associated(id%link)) call self%register_dependency(id, trim(name) // postfix, units, long_name)
      call self%request_coupling(id%link, trim(name) // postfix // '_calculator/result')
   end subroutine

   subroutine register_depth_explicit_state_dependency(self, link_int, name, units, long_name, proportional_change, domain, &
         link, depth_integral_out)
      class (type_base_model),              intent(inout) :: self
      type (type_link), target,             intent(inout) :: link_int
      character(len=*),                     intent(in)    :: name, units, long_name
      logical, optional,                    intent(in)    :: proportional_change
      integer, optional,                    intent(in)    :: domain
      type (type_link),                     optional, pointer :: link
      class (type_weighted_depth_integral), optional, pointer :: depth_integral_out

      class (type_weighted_depth_integral), pointer :: depth_integral

      ! Create a placeholder for the depth-explicit source variable
      if (present(link)) link => null()
      call self%add_interior_variable(name, units, long_name, act_as_state_variable=.true., link=link)

      ! Create a child model that will depth-integrate the source variable
      ! Couple this to the vertical distribution weights and our placeholder for the source variable
      allocate(depth_integral)
      depth_integral%act_as_state_variable = .true.
      if (present(domain)) depth_integral%domain = domain
      if (present(proportional_change)) depth_integral%proportional_change = proportional_change
      call self%add_child(depth_integral, name // '_int_calculator')
      call depth_integral%request_coupling('source', '../' // name)
      call depth_integral%request_coupling('w', '../w')

      ! Couple the provided variable link to the result of the depth integration
      call self%request_coupling(link_int, name // '_int_calculator/result')

      ! If requested, return the model instance that performs the depth integration
      if (present(depth_integral_out)) depth_integral_out => depth_integral
   end subroutine

   subroutine weighted_depth_integral_initialize(self, configunit)
      class (type_weighted_depth_integral), intent(inout), target :: self
      integer,                              intent(in)            :: configunit

      class (type_absolute_rate_distributor), pointer :: rate_distributor

      call self%register_implemented_routines((/source_do_column/))

      call self%get_parameter(self%average, 'average', '', 'average', default=self%average)
      call self%get_parameter(self%act_as_state_variable, 'act_as_state_variable', '', 'act as state variable', &
         default=self%act_as_state_variable)

      call self%register_diagnostic_variable(self%id_result, 'result', 'm', 'integral', source=source_do_column, &
         act_as_state_variable=self%act_as_state_variable, domain=self%domain)

      call self%register_dependency(self%id_w, 'w', '1', 'weights')
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)
      call self%register_dependency(self%id_source, 'source', '', 'source')
      if (self%average) call self%register_dependency(self%id_w_int, 'w_int', 'm', 'depth-integrated weights')

      if (self%act_as_state_variable) then
         ! Other model instances can provide source terms for the depth-integrated variable.
         ! Distribute these over the column and apply the result to the original depth-explicit state variable.
         call self%get_parameter(self%proportional_change, 'proportional_change', '', 'distribute sources as relative rate', &
            default=self%proportional_change)
         if (self%proportional_change) then
            allocate(type_relative_rate_distributor::rate_distributor)
         else
            allocate(rate_distributor)
         end if
         self%id_source%link%target%fake_state_variable = .true.
         call self%add_child(rate_distributor, 'sms_distributor')
         call rate_distributor%request_coupling('int', '../result')
         call rate_distributor%request_coupling('sms_int', '../result_sms_tot')
         call rate_distributor%request_coupling('w', '../w')
         call rate_distributor%request_coupling('w_int', '../w_int')
         call rate_distributor%request_coupling('target', '../source')
      end if
   end subroutine

   subroutine weighted_depth_integral_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_weighted_depth_integral), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: h, w, w_int, local_value, result

      result = 0.0_rk
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_h, h)
         _GET_(self%id_w, w)
         _GET_(self%id_source, local_value)
         result = result + w * h * local_value
      _VERTICAL_LOOP_END_

      if (self%average) then
         _GET_HORIZONTAL_(self%id_w_int, w_int)
         if (w_int > 0.0_rk) result = result / w_int
      end if
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_result, result)
   end subroutine

   subroutine projector_initialize(self, configunit)
      class (type_projector), intent(inout), target :: self
      integer,                intent(in)            :: configunit

      class (type_projected_rate_distributor), pointer :: rate_distributor

      call self%register_implemented_routines((/source_do/))

      call self%get_parameter(self%act_as_state_variable, 'act_as_state_variable', '', 'act as state variable', &
         default=self%act_as_state_variable)

      call self%register_diagnostic_variable(self%id_result, 'result', 'm', 'projection', &
         act_as_state_variable=self%act_as_state_variable)

      call self%register_dependency(self%id_source, 'source', '', 'source')
      call self%register_dependency(self%id_w, 'w', '1', 'weights')
      call self%register_dependency(self%id_w_int, 'w_int', 'm', 'depth-integrated weights')

      if (self%act_as_state_variable) then
         ! Other model instances can provide source terms for the projected variable.
         ! Make sure these are depth-integrated and applied to the original depth-integrated state variable.
         allocate(rate_distributor)
         self%id_source%link%target%fake_state_variable = .true.
         call self%add_child(rate_distributor, 'sms_distributor')
         call rate_distributor%request_coupling(rate_distributor%id_sms, '../result_sms_tot')
         call rate_distributor%request_coupling('target', '../source')
      end if
   end subroutine

   subroutine projector_do(self, _ARGUMENTS_DO_)
      class (type_projector), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: w, w_int, int_value

      _LOOP_BEGIN_
         _GET_(self%id_w, w)
         _GET_HORIZONTAL_(self%id_w_int, w_int)
         _GET_HORIZONTAL_(self%id_source, int_value)
         _SET_DIAGNOSTIC_(self%id_result, w / w_int * int_value)
      _LOOP_END_
   end subroutine

   subroutine absolute_rate_distributor_initialize(self, configunit)
      class (type_absolute_rate_distributor), intent(inout), target :: self
      integer,                                intent(in)            :: configunit

      call self%register_implemented_routines((/source_do/))
      call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
      call self%register_dependency(self%id_w, 'w', '1', 'weights for vertical distribution')
      call self%register_dependency(self%id_w_int, 'w_int', 'm', 'depth-integrated weights for vertical distribution')
      call self%register_dependency(self%id_int, 'int', '', 'depth-integrated target variable')
      call self%register_dependency(self%id_sms_int, 'sms_int', '', 'depth-integrated sources-sinks')
   end subroutine

   subroutine absolute_rate_distributor_do(self, _ARGUMENTS_DO_)
      class (type_absolute_rate_distributor), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sms_int, w_int, w, sms_per_w

      _LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_sms_int, sms_int)
         _GET_HORIZONTAL_(self%id_w_int, w_int)
         _GET_(self%id_w, w)

         ! The fraction of the depth-integrated source going to this cell is (w * h) / sum(w * h)
         ! To obtain the local depth-averaged source we divide by h, and thus get w / sum(w * h) = w / w_int
         sms_per_w = sms_int / w_int
         _ADD_SOURCE_(self%id_target, sms_per_w * w)
      _LOOP_END_
   end subroutine

   subroutine relative_rate_distributor_do(self, _ARGUMENTS_DO_)
      class (type_relative_rate_distributor), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sms_int, w, w_value_int, relative_change, local_value

      _LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_sms_int, sms_int)
         _GET_HORIZONTAL_(self%id_int, w_value_int)
         if (w_value_int > 0.0_rk) then
            relative_change = sms_int / w_value_int
         else
            relative_change = 0.0_rk
         end if

         _GET_(self%id_w, w)
         _GET_(self%id_target, local_value)
         _ADD_SOURCE_(self%id_target, relative_change * local_value * w)
      _LOOP_END_
   end subroutine

   subroutine projected_rate_distributor_initialize(self, configunit)
      class (type_projected_rate_distributor), intent(inout), target :: self
      integer,                                 intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_column/))
      call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
      call self%register_dependency(self%id_h, standard_variables%cell_thickness)
      call self%register_dependency(self%id_sms, 'sms', '', 'sources-sinks')
      self%id_target%bottom_sms%link%target%source = source_do_column
   end subroutine

   subroutine projected_rate_distributor_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_projected_rate_distributor), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: sms, h, sms_int

      sms_int = 0.0_rk
      _VERTICAL_LOOP_BEGIN_
         _GET_(self%id_sms, sms)
         _GET_(self%id_h, h)
         sms_int = sms_int + sms * h
      _VERTICAL_LOOP_END_
      _ADD_BOTTOM_SOURCE_(self%id_target, sms_int)
   end subroutine

   subroutine register_mapped_model_dependency(self, id, name, domain, proportional_change)
      class (type_depth_integrated_particle), intent(inout) :: self
      type (type_model_id), target, optional, intent(inout) :: id
      character(len=*),                       intent(in)    :: name
      integer, optional,                      intent(in)    :: domain
      logical, optional,                      intent(in)    :: proportional_change

      class (type_particle_integrator), pointer :: integrator

      ! Set up a child instance that depth-integrates all state variables of the
      ! depth-explicit (pelagic) model instance identified by the "name" argument
      allocate(integrator)
      if (present(domain)) integrator%domain = domain
      if (present(proportional_change)) integrator%proportional_change = proportional_change
      call self%add_child(integrator, name // '_integrator')

      ! Allow the user to coupled the requesting model instance (the parent of the integrator)
      ! to the named depth-explicit instance. Make sure the integrator uses that.
      call self%register_model_dependency(name)
      call integrator%request_coupling('source', '../' // name)

      ! Link the provided model identifier to the depth-integrated model instance
      call self%register_model_dependency(id, name // '_int')
      call self%request_coupling(name // '_int', './' // name // '_integrator')
   end subroutine

   subroutine depth_integrated_particle_initialize(self, configunit)
      class (type_depth_integrated_particle), intent(inout), target :: self
      integer,                                intent(in)            :: configunit

      call self%register_dependency(self%id_w_int, 'w_int', 'm', 'depth-integrated vertical distribution weights')
      call set_vertical_distribution(self, w_int_link=self%id_w_int%link)
   end subroutine

   subroutine depth_integrated_particle_request_mapped_coupling1(self, id, standard_variable, average)
      class (type_depth_integrated_particle),   intent(inout) :: self
      type (type_bottom_dependency_id), target, intent(inout) :: id
      type (type_interior_standard_variable),   intent(in)    :: standard_variable
      logical, optional,                        intent(in)    :: average

      type (type_link), pointer :: link

      call register_depth_explicit_dependency(self, id, trim(standard_variable%name), trim(standard_variable%units), &
         trim(standard_variable%name), link=link, average=average)
      call self%request_coupling(link, standard_variable)
   end subroutine

   subroutine depth_integrated_particle_request_mapped_coupling2(self, id, standard_variable, average)
      class (type_depth_integrated_particle),    intent(inout) :: self
      type (type_surface_dependency_id), target, intent(inout) :: id
      type (type_interior_standard_variable),    intent(in)    :: standard_variable
      logical, optional,                         intent(in)    :: average

      type (type_link), pointer :: link

      call register_depth_explicit_dependency(self, id, trim(standard_variable%name), trim(standard_variable%units), &
         trim(standard_variable%name), link=link, average=average)
      call self%request_coupling(link, standard_variable)
   end subroutine

   subroutine depth_integrated_particle_request_mapped_coupling3(self, id, standard_variable, average)
      class (type_depth_integrated_particle),   intent(inout) :: self
      type (type_bottom_dependency_id), target, intent(inout) :: id
      type (type_universal_standard_variable),  intent(in)    :: standard_variable
      logical, optional,                        intent(in)    :: average

      call self%request_mapped_coupling(id, standard_variable%in_interior(), average)
   end subroutine

   subroutine depth_integrated_particle_request_mapped_coupling4(self, id, standard_variable, average)
      class (type_depth_integrated_particle),    intent(inout) :: self
      type (type_surface_dependency_id), target, intent(inout) :: id
      type (type_universal_standard_variable),   intent(in)    :: standard_variable
      logical, optional,                         intent(in)    :: average

      call self%request_mapped_coupling(id, standard_variable%in_interior(), average)
   end subroutine

   subroutine depth_integrated_particle_request_mapped_coupling_to_model1(self, id, target_model_name, standard_variable, average)
      class (type_depth_integrated_particle),   intent(inout) :: self
      type (type_bottom_dependency_id), target, intent(inout) :: id
      character(len=*),                         intent(in)    :: target_model_name
      type (type_universal_standard_variable),  intent(in)    :: standard_variable
      logical, optional,                        intent(in)    :: average

      type (type_link), pointer :: link

      call register_depth_explicit_dependency(self, id, target_model_name // '_' // trim(standard_variable%name), &
         trim(standard_variable%units), target_model_name // ' ' // trim(standard_variable%name), link=link, average=average)
      call self%request_coupling_to_model_generic(link, target_model_name=target_model_name, target_standard_variable=standard_variable)
   end subroutine

   subroutine depth_integrated_particle_request_mapped_coupling_to_model2(self, id, target_model_name, standard_variable, average)
      class (type_depth_integrated_particle),    intent(inout) :: self
      type (type_surface_dependency_id), target, intent(inout) :: id
      character(len=*),                          intent(in)    :: target_model_name
      type (type_universal_standard_variable),   intent(in)    :: standard_variable
      logical, optional,                         intent(in)    :: average

      type (type_link), pointer :: link

      call register_depth_explicit_dependency(self, id, target_model_name // '_' // trim(standard_variable%name), &
         trim(standard_variable%units), target_model_name // ' ' // trim(standard_variable%name), link=link, average=average)
      call self%request_coupling_to_model_generic(link, target_model_name=target_model_name, target_standard_variable=standard_variable)
   end subroutine

   subroutine depth_integrated_particle_request_mapped_coupling_to_model3(self, id, target_model_name, standard_variable, proportional_change)
      class (type_depth_integrated_particle),  intent(inout)         :: self
      type (type_bottom_state_variable_id),    intent(inout), target :: id
      character(len=*),                        intent(in)            :: target_model_name
      type (type_universal_standard_variable), intent(in)            :: standard_variable
      logical, optional,                       intent(in)            :: proportional_change

      type (type_link), pointer :: link
      class (type_weighted_depth_integral),   pointer :: depth_integral

      call register_depth_explicit_state_dependency(self, id%link, target_model_name // '_' // trim(standard_variable%name), &
         trim(standard_variable%units), target_model_name // ' ' // trim(standard_variable%name), link=link, &
         proportional_change=proportional_change, depth_integral_out=depth_integral, domain=domain_bottom)
      call self%add_to_aggregate_variable(standard_variable, depth_integral%id_result)
      call self%request_coupling_to_model_generic(link, target_model_name=target_model_name, target_standard_variable=standard_variable)
   end subroutine

   subroutine depth_integrated_particle_request_mapped_coupling_to_model4(self, id, target_model_name, standard_variable, proportional_change)
      class (type_depth_integrated_particle),   intent(inout)         :: self
      type (type_surface_state_variable_id),    intent(inout), target :: id
      character(len=*),                         intent(in)            :: target_model_name
      type (type_universal_standard_variable),  intent(in)            :: standard_variable
      logical, optional,                        intent(in)            :: proportional_change

      type (type_link), pointer :: link
      class (type_weighted_depth_integral),   pointer :: depth_integral

      call register_depth_explicit_state_dependency(self, id%link, target_model_name // '_' // trim(standard_variable%name), &
         trim(standard_variable%units), target_model_name // ' ' // trim(standard_variable%name), link=link, &
         proportional_change=proportional_change, depth_integral_out=depth_integral, domain=domain_surface)
      call self%add_to_aggregate_variable(standard_variable, depth_integral%id_result)
      call self%request_coupling_to_model_generic(link, target_model_name=target_model_name, target_standard_variable=standard_variable)
   end subroutine

   subroutine particle_integrator_initialize(self, configunit)
      class (type_particle_integrator), intent(inout), target :: self
      integer,                          intent(in)            :: configunit

      ! Register a dependency on the depth-explicit model instance for which we will provide depth integrals.
      call self%register_model_dependency('source', 'model instance to depth-integrate')
   end subroutine

   recursive subroutine particle_integrator_complete_internal_variables(self)
      class (type_particle_integrator), intent(inout) :: self

      type (type_base_model),               pointer :: source
      type (type_link),                     pointer :: link, link_int
      class (type_weighted_depth_integral), pointer :: depth_integral
      type (type_contribution),             pointer :: contribution

      ! Retrieve the coupled depth-explcit model instance, ensuring that all its variables have been registered.
      source => self%resolve_model_dependency('source', require_internal_variables=.true.)

      ! Enumerate all state variables (and diagnostics acting as state variable) of the
      ! original depth explicit model instance, and for each create a depth-integrated equivalent.
      link => source%links%first
      do while (associated(link))
         if (index(link%name,'/') == 0 .and. link%target%presence == presence_internal &
             .and. (link%target%source == source_state .or. link%target%fake_state_variable)) then
            select case (link%target%domain)
            case (domain_interior)
               ! Create a child instance that depth-integrates the original depth-explicit state variable
               allocate(depth_integral)
               depth_integral%act_as_state_variable = .true.
               depth_integral%proportional_change = self%proportional_change
               depth_integral%domain = self%domain
               call self%add_child(depth_integral, trim(link%name) // '_int_calculator')

               ! Couple the weights and source variable of the depth integral
               call depth_integral%request_coupling(depth_integral%id_source, link)
               call depth_integral%request_coupling(depth_integral%id_w, '../w')

               ! Set up an alias for the depth-integrated variable (coupled to the result of the depth integrator)
               ! Note: the variable must be flagged as state variable (source=source_state or act_as_state_variable=.true.) in
               ! order to be included in state variable id lists set up when a particle-based coupling to this instance is made.
               link_int => null()
               call self%add_horizontal_variable(trim(link%name), trim(link%target%units) // ' m', 'depth-integrated ' &
                  // trim(link%target%long_name), source=source_state, link=link_int, domain=self%domain)
               call self%request_coupling(link_int, trim(link%name) // '_int_calculator/result')

               ! For the depth-integrated fake state variable, make sure it contributes to the same aggregate variable(s)
               ! that its depth-explicit source variable contributes to. This ensures conservation checks will work for
               ! models that provide source terms for the depth-integrated fake state variable.
               contribution => link%target%contributions%first
               do while (associated(contribution))
                  call depth_integral%add_to_aggregate_variable(contribution%target%universal, depth_integral%id_result, &
                     contribution%scale_factor, contribution%include_background)
                  contribution => contribution%next
               end do
            end select
         end if
         link => link%next
      end do
   end subroutine

   subroutine projected_particle_initialize(self, configunit)
      class (type_projected_particle), target, intent(inout) :: self
      integer,                                 intent(in)    :: configunit

      ! Register a dependency on the depth-integrated model instance for which we will
      ! create depth-explicit projections of each state variable
      call self%register_model_dependency('source', 'model instance to project over water column')
   end subroutine

   subroutine projected_particle_complete_internal_variables(self)
      class (type_projected_particle), intent(inout) :: self

      type (type_base_model),   pointer :: source
      type (type_link),         pointer :: link, link_int
      class (type_projector),   pointer :: projector
      type (type_contribution), pointer :: contribution

      ! Retrieve the coupled depth-explicit model, ensuring that all its variables have been registered.
      source => self%resolve_model_dependency('source', require_internal_variables=.true.)

      ! Enumerate all depth-integrated state variables (and diagnostics acting as state variable) of the
      ! coupled model instance, and for each create a depth-explicit equivalent.
      ! Projection over the water column is done by a child model
      link_int => source%links%first
      do while (associated(link_int))
         if (index(link_int%name,'/') == 0 .and. link_int%target%presence == presence_internal &
             .and. (link_int%target%source == source_state .or. link_int%target%fake_state_variable)) then
            select case (link_int%target%domain)
            case (domain_bottom, domain_surface)
               ! Create a child instance that projects the depth-integrated state variable over the pelagic
               allocate(projector)
               projector%act_as_state_variable = .true.
               call self%add_child(projector, trim(link_int%name) // '_calculator')

               ! Couple the weights and source variable of the depth integral
               call projector%request_coupling(projector%id_source, link_int)
               call projector%request_coupling(projector%id_w, '../w')
               call projector%request_coupling(projector%id_w_int, '../w_int')

               ! Set up a alias for the depth-explicit variable (coupled to the result of the child model that projects over depth)
               ! Note: the variable must be flagged as state variable (source=source_state or act_as_state_variable=.true.) in
               ! order to be included in state variable id lists set up when a particle-based coupling to this instance is made.
               link => null()
               call self%add_interior_variable(trim(link_int%name), trim(link_int%target%units) // ' m-1', 'depth-explicit ' &
                  // trim(link_int%target%long_name), source=source_state, link=link)
               call self%request_coupling(link, trim(link_int%name) // '_calculator/result')

               ! For the depth-explicit fake state variable, make sure it contributes to the same aggregate variable(s)
               ! that its depth-integrated source variable contributes to. This ensures conservation checks will work for
               ! models that provide source terms for the projected state variable.
               contribution => link_int%target%contributions%first
               do while (associated(contribution))
                  call projector%add_to_aggregate_variable(contribution%target%universal, projector%id_result, &
                     contribution%scale_factor, contribution%include_background)
                  contribution => contribution%next
               end do
            end select
         end if
         link_int => link_int%next
      end do
   end subroutine
   
   subroutine depth_integrated_particle_override_initialize(self, configunit)
      class (type_depth_integrated_particle_override), target, intent(inout) :: self
      integer,                                                 intent(in)    :: configunit

      type (type_projected_particle), pointer :: projection

      ! Alllow the base class to initialize
      ! This sets up a dependency on vertial distribution weights w
      call self%type_depth_integrated_particle%initialize(configunit)

      ! Register a dependence on the depth-explicit model instance for which we will
      ! replace interior state variables with depth-integrated ones
      call self%register_model_dependency('target', 'model instance in which to replace depth-explicit state variables with projected depth-integrated ones')

      ! Add a child model that will project our depth-integrated variables over the water column
      allocate(projection)
      call self%add_child(projection, 'projection')
      call projection%request_coupling('source', '..')
   end subroutine

   subroutine depth_integrated_particle_override_complete_internal_variables(self)
      class (type_depth_integrated_particle_override), intent(inout) :: self

      type (type_base_model),   pointer :: target_instance
      type (type_link),         pointer :: link, link_int
      type (type_contribution), pointer :: contribution

      ! Retrieve the coupled depth-explicit model instance, ensuring that all its variables have been registered.
      target_instance => self%resolve_model_dependency('target', require_internal_variables=.true.)

      ! Enumerate all depth-explicit state variables from coupled model instance,
      ! and for each create a depth-integrated equivalent.
      link => target_instance%links%first
      do while (associated(link))
         if (index(link%name,'/') == 0 .and. link%target%presence == presence_internal &
             .and. link%target%source == source_state .and. link%target%domain == domain_interior) then

            ! Create the depth-integrated state variable that will replace the original depth-explicit state variable
            link_int => null()
            call self%add_horizontal_variable(trim(link%name), trim(link%target%units) // ' m', 'depth-integrated ' &
               // trim(link%target%long_name), source=source_state, link=link_int, domain=self%domain)

            ! Ensure the depth-integrated state variable contributes to the same aggregate variable(s) that the original
            ! depth-explicit source variable contributes to. Make sure to register contributions to bottom/surface as
            ! well as to interfaces in general. The latter is used to compute totals of conserved quantities.
            contribution => link%target%contributions%first
            do while (associated(contribution))
               call link_int%target%contributions%add(contribution%target%universal%at_interfaces(), &
                  contribution%scale_factor, contribution%include_background)
               if (self%domain == domain_bottom) then
                  call link_int%target%contributions%add(contribution%target%universal%at_bottom(), &
                     contribution%scale_factor, contribution%include_background)
               else
                  call link_int%target%contributions%add(contribution%target%universal%at_surface(), &
                     contribution%scale_factor, contribution%include_background)
               end if
               contribution => contribution%next
            end do

            ! Couple the original depth-explicit state variable to the projected diagnostic
            call target_instance%request_coupling(link, trim(self%get_path()) // '/projection/' // trim(link%name))
         end if
         link => link%next
      end do
   end subroutine

end module

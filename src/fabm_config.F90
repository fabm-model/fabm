#include "fabm_driver.h"

module fabm_config

   use fabm_types
   use fabm_properties, only: type_property_dictionary, type_property, type_set
   use fabm_driver
   use fabm_schedule

   use yaml_settings

   implicit none

   private

   public fabm_configure_model

   type, extends(type_dictionary_populator) :: type_instances_populator
      class (type_base_model), pointer :: root => null()
      class (type_schedules), pointer :: schedules => null()
      logical :: check_conservation, require_initialization, require_all_parameters
   contains
      procedure :: create => create_instance
   end type

contains

   subroutine fabm_configure_model(root, settings, schedules, log, path, parameters, unit)
      class (type_base_model),                   intent(inout) :: root
      type (type_settings),                      intent(out)   :: settings
      class (type_schedules),                    intent(inout) :: schedules
      logical,                                   intent(out)   :: log
      character(len=*),                optional, intent(in)    :: path
      type (type_property_dictionary), optional, intent(in)    :: parameters
      integer,                         optional, intent(in)    :: unit

      integer                          :: unit_eff
      character(len=256)               :: path_eff

      ! Determine the path to use for YAML file.
      if (present(path)) then
          path_eff = trim(path)
      else
          path_eff = 'fabm.yaml'
      end if

      ! Determine the unit to use for YAML file.
      if (present(unit)) then
          unit_eff = unit
      else
          unit_eff = get_free_unit()
      end if

      ! Parse YAML file.
      call settings%load(trim(path_eff), unit_eff)
      ! TODO :check for errors

      ! If custom parameter values were provided, transfer these to the root model.
      !TODO if (present(parameters)) call root%parameters%update(parameters)

      log = settings%get_logical('log', 'write log files for debugging FABM', default=.false.)
      ! TODO :check for errors
      call create_model_tree_from_dictionary(root, settings, schedules)
   end subroutine fabm_configure_model

   subroutine create_model_tree_from_dictionary(root, settings, schedules)
      class (type_base_model), intent(inout), target :: root
      class (type_settings),   intent(inout) :: settings
      class (type_schedules),  intent(inout), target :: schedules

      type (type_instances_populator) :: instances_populator
      class (type_settings), pointer  :: instances

      instances_populator%root => root
      instances_populator%schedules => schedules
      instances_populator%check_conservation = settings%get_logical('check_conservation', 'add diagnostics for the rate of change of conserved quantities', default=.false.)
      ! TODO :check for errors
      instances_populator%require_initialization = settings%get_logical('require_initialization', 'require initial values for all state variables', default=.false.)
      ! TODO :check for errors
      instances_populator%require_all_parameters = settings%get_logical('require_all_parameters', 'require values for all parameters', default=.false.)
      ! TODO :check for errors

      instances => settings%get_child('instances', populator=instances_populator)
      ! TODO :check for errors
   end subroutine create_model_tree_from_dictionary

   recursive subroutine create_instance(self, pair)
      class (type_instances_populator), intent(inout) :: self
      type (type_key_value_pair),       intent(inout) :: pair

      class (type_settings), pointer :: instance_settings
      logical                             :: use_model
      character(len=:), allocatable       :: modelname
      character(len=:), allocatable       :: long_name
      class (type_settings), pointer :: parameters, schedule
      type (type_dictionary)              :: parametermap
      class (type_dictionary),    pointer :: childmap
      class (type_property),      pointer :: property
      type (type_key_value_pair), pointer :: pair
      type (type_set)                     :: initialized_set,background_set
      type (type_link),           pointer :: link
      type (type_error),          pointer :: config_error
      integer                             :: schedule_pattern, source
      character(len=64)                   :: pattern
      class (type_base_model), pointer :: model

      instance_settings => type_settings_create(pair)

      use_model = instance_settings%get_logical('use', 'use this model', default=.true.)
      ! TODO :check for errors
      if (.not. use_model) then
         call log_message('SKIPPING model instance ' // pair%name // ' because it has use=false set.')
         return
      end if

      ! Retrieve model name (default to instance name if not provided).
      modelname = instance_settings%get_string('model', 'model type', default=pair%name)
      ! TODO :check for errors

      ! Retrieve descriptive name for the model instance (default to instance name if not provided).
      long_name = instance_settings%get_string('long_name', 'descriptive name for model instance', default=pair%name)
      ! TODO :check for errors

      ! Try to create the model based on name.
      call factory%create(modelname, model)
      if (.not. associated(model)) then
         call fatal_error('create_model_from_dictionary', &
            pair%name // ': "' // modelname // '" is not a valid model name.')
         return
      end if
      model%user_created = .true.

      ! Transfer user-specified parameter values to the model.
      model%parameters => instance_settings%get_child('parameters')
      model%couplings => instance_settings%get_child('coupling')
      ! TODO :check for errors

      ! Add the model to its parent.
      call log_message('Initializing ' // trim(instancename) // '...')
      call log_message('   model type: ' // trim(modelname))
      call parent%add_child(model, instancename, long_name)

      ! TODO Check for parameters requested by the model, but not present in the configuration file.

      call log_message('   initialization succeeded.')

      ! Interpret coupling links specified in configuration file.
      ! These override any couplings requested by the models during initialization.
      ! This step must therefore occur after model initialization [from parent%add_child] has completed.
      childmap => node%get_dictionary('coupling', required=.false., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_from_dictionary', config_error%message)
      if (associated(childmap)) then
         pair => childmap%first
         do while (associated(pair))
            select type (value => pair%value)
            class is (type_scalar)
               ! Register couplings at the root level, so they override whatever the models themselves request.
               call parent%couplings%set_string(trim(instancename) // '/' // trim(pair%key), trim(value%string))
            class is (type_node)
               call fatal_error('create_model_from_dictionary', 'The value of ' // trim(value%path) // &
                  ' must be a string, not a nested dictionary.')
            end select
            pair => pair%next
         end do
      end if

      ! Parse scheduling instructions
      schedule => instance_settings%get_child('schedule')
      ! TODO: check errors
      if (associated(childmap)) then
         pair => childmap%first
         do while (associated(pair))
            select type (value => pair%value)
            class is (type_dictionary)
               select case (pair%key)
               case ('interior')
                  source = source_do
               case default
                  call fatal_error('create_model_from_dictionary', 'Scheduler currently only supports "interior" &
                     &(not "' // trim(pair%key) // '").')
               end select
               pattern = trim(value%get_string('pattern', error=config_error))
               if (associated(config_error)) call fatal_error('create_model_from_dictionary', config_error%message)
               select case (pattern)
               case ('monthly')
                  schedule_pattern = schedule_pattern_monthly
               case default
                  call fatal_error('create_model_from_dictionary', 'Scheduler currently only supports "monthly" &
                     &as a pattern (not "' // trim(pattern) // '").')
               end select
               call schedules%add(model, source, schedule_pattern)
            class is (type_node)
               call fatal_error('create_model_from_dictionary','The value of ' // trim(value%path) // &
                  ' must be a dictionary.')
            end select
            pair => pair%next
         end do
      end if

      ! Transfer user-specified initial state to the model.
      childmap => node%get_dictionary('initialization', required=.false., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_from_dictionary', config_error%message)
      if (associated(childmap)) call parse_initialization(model, childmap, initialized_set, get_background=.false.)

      ! Transfer user-specified background value to the model.
      childmap => node%get_dictionary('background', required=.false., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_from_dictionary', config_error%message)
      if (associated(childmap)) call parse_initialization(model, childmap, background_set, get_background=.true.)

      ! Verify whether all state variables have been provided with an initial value.
      !link => model%first_link
      !do while (associated(link))
      !   if (link%owner.and..not.associated(model%couplings%find(link%name))) then
      !      ! This link is our own: not coupled to another variable, not an alias, and no intention to couple was registered.
      !            if (object%source == source_state) then
      !               if (object%presence/=presence_external_optional .and. .not.initialized_set%contains(trim(link%name))) then
      !                  ! State variable is required, but initial value is not explicitly provided.
      !                  if (require_initialization) then
      !                     call fatal_error('parse_initialization','model '//trim(model%name) &
      !                        //': no initial value provided for variable "'//trim(link%name)//'".')
      !                  else
      !                     call log_message('WARNING: no initial value provided for state variable "'//trim(link%name)// &
      !                        '" of model "'//trim(model%name)//'" - using default.')
      !                  end if
      !               elseif (object%presence==presence_external_optional .and. initialized_set%contains(trim(link%name))) then
      !                  ! Optional state variable is not used, but an initial value is explicitly provided.
      !                  call fatal_error('parse_initialization','model '//trim(model%name) &
      !                     //': initial value provided for variable "'//trim(link%name)//'", but this variable is not used.')
      !               end if
      !            end if
      !   end if
      !   link => link%next
      !end do
      call initialized_set%finalize()

      model%check_conservation = node%get_logical('check_conservation', default=self%check_conservation)
      ! TODO: check errors
   end subroutine create_model_from_dictionary

   subroutine parse_initialization(model, node, initialized_set, get_background)
      class (type_base_model), intent(inout) :: model
      class (type_dictionary), intent(in)    :: node
      type (type_set),         intent(out)   :: initialized_set
      logical,                 intent(in)    :: get_background

      type (type_key_value_pair),    pointer :: pair
      type (type_internal_variable), pointer :: object
      logical                                :: success
      real(rk)                               :: realvalue

      ! Transfer user-specified initial state to the model.
      pair => node%first
      do while (associated(pair))
         select type (value => pair%value)
         class is (type_scalar)
            object => model%find_object(trim(pair%key))
            if (.not. associated(object)) then
               call fatal_error('parse_initialization', &
                  trim(value%path) // ': "' // trim(pair%key) // '" is not a member of model "' // trim(model%name) // '".')
               return
            end if
            if (object%source /= source_state) call fatal_error('parse_initialization', &
               trim(value%path) // ': "' // trim(pair%key) // '" is not a state variable of model "' // trim(model%name) // '".')
            realvalue = value%to_real(default=real(0, real_kind), success=success)
            if (.not. success) call fatal_error('parse_initialization', &
               trim(value%path) // ': "' // trim(value%string) // '" is not a real number.')
            if (get_background) then
               call object%background_values%set_value(realvalue)
            else
               object%initial_value = realvalue
            end if
            call initialized_set%add(trim(pair%key))
         class is (type_null)
            call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to null.')
         class is (type_dictionary)
            call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to a dictionary.')
         end select
         pair => pair%next
      end do
   end subroutine

end module fabm_config

module fabm_config_types
   use yaml_types
   implicit none
   public
end module

module fabm_yaml
   use yaml
   implicit none
   public
end module

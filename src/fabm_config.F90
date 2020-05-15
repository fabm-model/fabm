#include "fabm_driver.h"

module fabm_config

   use fabm_types
   use fabm_properties, only: type_property_dictionary, type_property, type_set
   use fabm_driver
   use fabm_schedule

   use yaml_types
   use yaml, yaml_parse=>parse, yaml_error_length=>error_length

   implicit none

   private

   public fabm_configure_model

contains

   subroutine fabm_configure_model(root, schedules, log, path, parameters, unit)
      class (type_base_model),                   intent(inout) :: root
      class (type_schedules),                    intent(inout) :: schedules
      logical,                                   intent(out)   :: log
      character(len=*),                optional, intent(in)    :: path
      type (type_property_dictionary), optional, intent(in)    :: parameters
      integer,                         optional, intent(in)    :: unit

      class (type_node), pointer       :: node
      character(len=yaml_error_length) :: yaml_error
      integer                          :: unit_eff
      character(len=256)               :: path_eff
      type (type_error), pointer       :: config_error

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
      node => yaml_parse(trim(path_eff), unit_eff, yaml_error)
      if (yaml_error /= '') call fatal_error('fabm_configure_model', trim(yaml_error))
      if (.not. associated(node)) call fatal_error('fabm_configure_model', &
         'No configuration information found in ' // trim(path_eff) // '.')
      !call node%dump(output_unit,0)

      ! If custom parameter values were provided, transfer these to the root model.
      if (present(parameters)) call root%parameters%update(parameters)

      ! Create model tree from YAML root node.
      select type (node)
      class is (type_dictionary)
         log = node%get_logical('log', default=.false., error=config_error)
         if (associated(config_error)) call fatal_error('fabm_configure_model', config_error%message)
         call create_model_tree_from_dictionary(root, node, schedules)
      class is (type_node)
         call fatal_error('fabm_configure_model', trim(path_eff) // ' must contain a dictionary &
            &at the root (non-indented) level, not a single value. Are you missing a trailing colon?')
      end select
   end subroutine fabm_configure_model

   subroutine create_model_tree_from_dictionary(root, mapping, schedules)
      class (type_base_model), intent(inout) :: root
      class (type_dictionary), intent(in)    :: mapping
      class (type_schedules),  intent(inout) :: schedules

      class (type_node),          pointer :: node
      type (type_dictionary)              :: empty_dict
      character(len=64)                   :: instancename
      type (type_key_value_pair), pointer :: pair
      logical                             :: check_conservation, require_initialization, require_all_parameters
      type (type_error),          pointer :: config_error

      config_error => null()
      check_conservation = mapping%get_logical('check_conservation', default=.false., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_tree_from_dictionary', config_error%message)
      require_initialization = mapping%get_logical('require_initialization', default=.false., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_tree_from_dictionary', config_error%message)
      require_all_parameters = mapping%get_logical('require_all_parameters', default=.false., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_tree_from_dictionary', config_error%message)

      node => mapping%get('instances')
      if (.not. associated(node)) &
         call fatal_error('create_model_tree_from_dictionary', 'No "instances" dictionary found at root level.')
      pair => null()
      select type (node)
      class is (type_dictionary)
         pair => node%first
      class is (type_null)
      class is (type_node)
         call fatal_error('create_model_tree_from_dictionary', trim(node%path) // &
            ' must be a dictionary with (model name : information) pairs, not a single value.')
      end select

      if (.not. associated(pair)) &
         call log_message('WARNING: no model instances specified. FABM is effectively disabled.')

      ! Iterate over all models (key:value pairs below "instances" node at root level) and
      ! create corresponding objects.
      do while (associated(pair))
         instancename = trim(pair%key)
         select type (dict => pair%value)
         class is (type_dictionary)
            call create_model_from_dictionary(instancename, dict, root, &
                                              require_initialization, require_all_parameters, check_conservation, schedules)
         class is (type_null)
            call create_model_from_dictionary(instancename, empty_dict, root, &
                                              require_initialization, require_all_parameters, check_conservation, schedules)
         class is (type_node)
            call fatal_error('create_model_tree_from_dictionary', 'Configuration information for model "' // &
               trim(instancename) // '" must be a dictionary, not a single value.')
         end select
         pair => pair%next
      end do

      ! Check whether any keys at the root level remain unused.
      pair => mapping%first
      do while (associated(pair))
         if (.not. pair%accessed) call fatal_error('create_model_tree_from_dictionary', 'Unrecognized option "' // &
            trim(pair%key) // '" found at root level.')
         pair => pair%next
      end do
   end subroutine create_model_tree_from_dictionary

   subroutine create_model_from_dictionary(instancename, node, parent, &
                                           require_initialization, require_all_parameters, check_conservation, schedules)
      character(len=*),        intent(in)            :: instancename
      class (type_dictionary), intent(in)            :: node
      class (type_base_model), intent(inout), target :: parent
      logical,                 intent(in)            :: require_initialization, require_all_parameters, check_conservation
      class (type_schedules),  intent(inout)         :: schedules
      class (type_base_model), pointer               :: model

      logical                             :: use_model
      character(len=64)                   :: modelname
      character(len=256)                  :: long_name
      type (type_dictionary)              :: parametermap
      class (type_dictionary),    pointer :: childmap
      class (type_property),      pointer :: property
      type (type_key_value_pair), pointer :: pair
      type (type_set)                     :: initialized_set,background_set
      type (type_link),           pointer :: link
      type (type_error),          pointer :: config_error
      integer                             :: schedule_pattern, source
      character(len=64)                   :: pattern

      config_error => null()

      use_model = node%get_logical('use', default=.true., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_from_dictionary', config_error%message)
      if (.not. use_model) then
         call log_message('SKIPPING model instance ' // trim(instancename) // ' because it has use=false set.')
         return
      end if

      ! Retrieve model name (default to instance name if not provided).
      modelname = trim(node%get_string('model', default=instancename, error=config_error))
      if (associated(config_error)) then
         call fatal_error('create_model_from_dictionary', config_error%message)
         return
      end if

      ! Retrieve descriptive name for the model instance (default to instance name if not provided).
      long_name = trim(node%get_string('long_name', default=instancename, error=config_error))
      if (associated(config_error)) then
         call fatal_error('create_model_from_dictionary', config_error%message)
         return
      end if

      ! Try to create the model based on name.
      call factory%create(trim(modelname), model)
      if (.not. associated(model)) then
         call fatal_error('create_model_from_dictionary', &
            trim(instancename) // ': "' // trim(modelname) // '" is not a valid model name.')
         return
      end if
      model%user_created = .true.

      ! Transfer user-specified parameter values to the model.
      childmap => node%get_dictionary('parameters', required=.false., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_from_dictionary', config_error%message)
      if (associated(childmap)) then
         call childmap%flatten(parametermap, '')
         pair => parametermap%first
         do while (associated(pair))
            select type (value => pair%value)
            class is (type_scalar)
               call model%parameters%set_string(trim(pair%key), trim(value%string))
            class is (type_node)
               call fatal_error('create_model_from_dictionary', 'BUG: "flatten" should &
                  &have ensured that the value of ' // trim(value%path) // ' is scalar, not a nested dictionary.')
            end select
            pair => pair%next
         end do
         call parametermap%finalize()
      end if

      ! Add the model to its parent.
      call log_message('Initializing ' // trim(instancename) // '...')
      call log_message('   model type: ' // trim(modelname))
      call parent%add_child(model, instancename, long_name, configunit=-1)

      ! Check for parameters requested by the model, but not present in the configuration file.
      if (require_all_parameters .and. associated(model%parameters%missing%first)) then
         call fatal_error('create_model_from_dictionary', 'Value for parameter "'// &
            trim(model%parameters%missing%first%string) // '" of model "' // trim(instancename) // '" is not provided.')
         return
      end if

      ! Check for parameters present in configuration file, but not interpreted by the models.
      property => model%parameters%first
      do while (associated(property))
         if (.not. model%parameters%retrieved%contains(property%name)) then
            call fatal_error('create_model_from_dictionary', &
               'Unrecognized parameter "' // trim(property%name) // '" found below ' // trim(childmap%path) // '.')
            return
         end if
         property => property%next
      end do

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
      childmap => node%get_dictionary('schedule', required=.false., error=config_error)
      if (associated(config_error)) call fatal_error('create_model_from_dictionary', config_error%message)
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

      model%check_conservation = node%get_logical('check_conservation', default=check_conservation, error=config_error)
      if (associated(config_error)) call fatal_error('create_model_from_dictionary', config_error%message)

      ! Check whether any keys at the model level remain unused.
      pair => node%first
      do while (associated(pair))
         if (.not. pair%accessed) call fatal_error('create_model_from_dictionary', &
            'Unrecognized option "' // trim(pair%key) // '" found below ' // trim(node%path) // '.')
         pair => pair%next
      end do

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

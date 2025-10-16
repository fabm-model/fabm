#include "fabm_driver.h"

module fabm_config

   use fabm_types
   use fabm_driver
   use fabm_schedule

   use yaml_settings
   use yaml_types, only: yaml_rk => real_kind

   implicit none

   private

   public fabm_configure_model, fabm_load_settings

   type, extends(type_dictionary_populator) :: type_instances_populator
      class (type_base_model), pointer :: root      => null()
      class (type_schedules),  pointer :: schedules => null()
      logical :: check_conservation
      logical :: require_all_parameters
   contains
      procedure :: create => create_instance
   end type

contains

   subroutine fabm_load_settings(settings, path, unit)
      class (type_fabm_settings),               intent(inout) :: settings
      character(len=*),               optional, intent(in)    :: path
      integer,                        optional, intent(in)    :: unit

      integer :: unit_

      ! Determine the unit to use for YAML file.
      if (present(unit)) then
          unit_ = unit
      else
          unit_ = get_free_unit()
      end if

      ! Parse YAML file.
      if (present(path)) then
          call settings%load(path, unit_, error_reporter=yaml_settings_error_reporter)
      else
          call settings%load('fabm.yaml', unit_, error_reporter=yaml_settings_error_reporter)
      end if
   end subroutine fabm_load_settings

   subroutine fabm_configure_model(root, settings, schedules, log, require_initialization)
      class (type_base_model),    target, intent(inout) :: root
      class (type_fabm_settings),         intent(inout) :: settings
      class (type_schedules),     target, intent(inout) :: schedules
      logical,                    target, intent(inout) :: log
      logical,                    target, intent(inout) :: require_initialization

      type (type_instances_populator) :: instances_populator
      class (type_settings), pointer  :: instances

      call settings%get(log, 'log', 'write log files for debugging FABM', default=.false., display=display_advanced)
      call settings%get(require_initialization, 'require_initialization', 'require initial values for all state variables', default=.false., display=display_advanced)

      instances_populator%root => root
      instances_populator%schedules => schedules
      instances_populator%check_conservation = settings%get_logical('check_conservation', 'add diagnostics for the rate of change of conserved quantities', default=.false., display=display_advanced)
      instances_populator%require_all_parameters = settings%get_logical('require_all_parameters', 'require values for all parameters', default=.false., display=display_hidden)
      instances => settings%get_child('instances', populator=instances_populator)
   end subroutine fabm_configure_model

   recursive subroutine create_instance(self, pair)
      class (type_instances_populator), intent(inout) :: self
      type (type_key_value_pair),       intent(inout) :: pair

      class (type_settings),      pointer :: subsettings
      class (type_fabm_settings), pointer :: instance_settings
      logical                             :: use_model, ignored
      character(len=:), allocatable       :: modelname
      character(len=:), allocatable       :: long_name
      type (type_link),           pointer :: link
      class (type_base_model),    pointer :: model
      integer                             :: schedule_pattern
      real(rk)                            :: realvalue

      subsettings => type_settings_create(pair)
      select type (subsettings)
      class is (type_fabm_settings)
         instance_settings => subsettings
      end select

      use_model = instance_settings%get_logical('use', 'use this model', default=.true., display=display_advanced)
      if (.not. use_model) then
         call log_message('SKIPPING model instance ' // pair%name // ' because it has use=false set.')
         ignored = instance_settings%ignore()
         return
      end if

      ! Retrieve model name (default to instance name if not provided).
      modelname = instance_settings%get_string('model', 'model type', default=pair%name)

      ! Retrieve descriptive name for the model instance (default to instance name if not provided).
      long_name = instance_settings%get_string('long_name', 'descriptive name for model instance', default=pair%name, display=display_advanced)

      ! Try to create the model based on name.
      call factory%create(modelname, model)
      if (.not. associated(model)) then
         call fatal_error('create_model_from_dictionary', &
            pair%name // ': "' // modelname // '" is not a valid model name.')
         return
      end if
      model%user_created = .true.

      call instance_settings%get(model%check_conservation, 'check_conservation', 'check_conservation', default=self%check_conservation, display=display_advanced)

      ! Transfer user-specified parameter values to the model.
      call instance_settings%attach_child(model%parameters, 'parameters', display=display_advanced)
      call instance_settings%attach_child(model%initialization, 'initialization', display=display_advanced)
      call instance_settings%attach_child(model%couplings, 'coupling', display=display_advanced)

      ! Add the model to its parent.
      call log_message('Initializing ' // pair%name // '...')
      call log_message('   model type: ' // modelname)
      call self%root%add_child(model, pair%name, long_name)
      call log_message('   initialization succeeded.')

      subsettings => instance_settings%get_child('schedule', display=display_advanced)
      schedule_pattern = subsettings%get_integer('interior', 'interior', options=(/option(0, 'always', 'always'), &
         option(1, 'monthly', 'monthly')/), default=0)
      if (schedule_pattern /= 0) call self%schedules%add(model, source_do, schedule_pattern)

      ! Transfer user-specified background value to the model.
      subsettings => instance_settings%get_child('background', display=display_advanced)
      link => model%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. link%target%source == source_state .and. link%target%presence == presence_internal) then
            realvalue = subsettings%get_real(trim(link%name), trim(link%target%long_name), trim(link%target%units), &
               default=real(link%target%background_values%value, yaml_rk), display=display_advanced)
            call link%target%background_values%set_value(realvalue)
         end if
         link => link%next
      end do

   end subroutine create_instance

   subroutine yaml_settings_error_reporter(message)
      character(len=*), intent(in) :: message
      call driver%fatal_error('yaml_settings', message)
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

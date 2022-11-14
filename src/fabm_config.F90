#include "fabm_driver.h"

module fabm_config

   use fabm_types
   use fabm_driver
   use fabm_schedule

   use yaml_settings

   implicit none

   private

   public fabm_configure_model, fabm_load_settings

   type, extends(type_dictionary_populator) :: type_instances_populator
      class (type_base_model), pointer :: root      => null()
      class (type_schedules),  pointer :: schedules => null()
      logical :: check_conservation
      logical :: require_initialization
      logical :: require_all_parameters
   contains
      procedure :: create => create_instance
   end type

contains

   subroutine fabm_load_settings(settings, path, unit, error_reporter)
      class (type_fabm_settings),               intent(inout) :: settings
      character(len=*),               optional, intent(in)    :: path
      integer,                        optional, intent(in)    :: unit
      procedure(error_reporter_proc), optional                :: error_reporter

      integer :: unit_

      ! Determine the unit to use for YAML file.
      if (present(unit)) then
          unit_ = unit
      else
          unit_ = get_free_unit()
      end if

      ! Parse YAML file.
      if (present(path)) then
          call settings%load(path, unit_, error_reporter=error_reporter)
      else
          call settings%load('fabm.yaml', unit_, error_reporter=error_reporter)
      end if
   end subroutine fabm_load_settings

   subroutine fabm_configure_model(root, settings, schedules, log)
      class (type_base_model),    target, intent(inout) :: root
      class (type_fabm_settings),         intent(inout) :: settings
      class (type_schedules),     target, intent(inout) :: schedules
      logical,                    target, intent(inout) :: log

      type (type_instances_populator) :: instances_populator
      class (type_settings), pointer  :: instances

      call settings%get(log, 'log', 'write log files for debugging FABM', default=.false.)

      instances_populator%root => root
      instances_populator%schedules => schedules
      instances_populator%check_conservation = settings%get_logical('check_conservation', 'add diagnostics for the rate of change of conserved quantities', default=.false.)
      instances_populator%require_initialization = settings%get_logical('require_initialization', 'require initial values for all state variables', default=.false.)
      instances_populator%require_all_parameters = settings%get_logical('require_all_parameters', 'require values for all parameters', default=.false.)
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

      use_model = instance_settings%get_logical('use', 'use this model', default=.true.)
      if (.not. use_model) then
         call log_message('SKIPPING model instance ' // pair%name // ' because it has use=false set.')
         ignored = instance_settings%ignore()
         return
      end if

      ! Retrieve model name (default to instance name if not provided).
      modelname = instance_settings%get_string('model', 'model type', default=pair%name)

      ! Retrieve descriptive name for the model instance (default to instance name if not provided).
      long_name = instance_settings%get_string('long_name', 'descriptive name for model instance', default=pair%name)

      ! Try to create the model based on name.
      call factory%create(modelname, model)
      if (.not. associated(model)) then
         call fatal_error('create_model_from_dictionary', &
            pair%name // ': "' // modelname // '" is not a valid model name.')
         return
      end if
      model%user_created = .true.

      call instance_settings%get(model%check_conservation, 'check_conservation', 'check_conservation', default=self%check_conservation)

      ! Transfer user-specified parameter values to the model.
      call instance_settings%attach_child(model%parameters, 'parameters')
      call instance_settings%attach_child(model%couplings, 'coupling')

      ! Add the model to its parent.
      call log_message('Initializing ' // pair%name // '...')
      call log_message('   model type: ' // modelname)
      call self%root%add_child(model, pair%name, long_name)
      call log_message('   initialization succeeded.')

      subsettings => instance_settings%get_child('schedule')
      schedule_pattern = subsettings%get_integer('interior', 'interior', options=(/option(0, 'always', 'always'), &
         option(1, 'monthly', 'monthly')/), default=0, display=display_advanced)
      if (schedule_pattern /= 0) call self%schedules%add(model, source_do, schedule_pattern)

      ! Transfer user-specified initial state to the model.
      subsettings => instance_settings%get_child('initialization')
      link => model%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. link%target%source == source_state .and. link%target%presence == presence_internal) then
            if (self%require_initialization) then
               call subsettings%get(link%target%initial_value, trim(link%name), trim(link%target%long_name), &
                  trim(link%target%units), minimum=link%target%minimum, maximum=link%target%maximum)
            else
               call subsettings%get(link%target%initial_value, trim(link%name), trim(link%target%long_name), &
                  trim(link%target%units), minimum=link%target%minimum, maximum=link%target%maximum, &
                  default=link%target%initial_value)
            end if
         end if
         link => link%next
      end do

      ! Transfer user-specified background value to the model.
      subsettings => instance_settings%get_child('background')
      link => model%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. link%target%source == source_state &
             .and. allocated(link%target%background_values%pointers)) then
            realvalue = subsettings%get_real(trim(link%name), trim(link%target%long_name), trim(link%target%units), &
               default=link%target%background_values%pointers(1)%p)
            call link%target%background_values%set_value(realvalue)
         end if
         link => link%next
      end do

   end subroutine create_instance

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

#ifdef _FABM_F2003_

#include "fabm_driver.h"

module fabm_config

   use fabm_types
   use fabm_properties,only:type_property_dictionary,type_property,type_set
   use fabm,only:type_model,fabm_initialize

   use fabm_config_types
   use fabm_yaml,yaml_parse=>parse,yaml_error_length=>error_length

   implicit none

   private

   public fabm_create_model_from_yaml_file

contains

   subroutine fabm_create_model_from_yaml_file(model,do_not_initialize,parameters)
      type (type_model),                       intent(out) :: model
      logical,                        optional,intent(in)  :: do_not_initialize
      type (type_property_dictionary),optional,intent(in)  :: parameters

      class (type_node),pointer        :: node
      character(len=yaml_error_length) :: yaml_error
      character(len=*),parameter       :: path = 'fabm.yaml'

      if (.not.associated(driver)) allocate(type_base_driver::driver)

      ! Parse YAML file.
      node => yaml_parse(path,get_free_unit(),yaml_error)
      if (yaml_error/='') call driver%fatal_error('fabm_create_model_from_yaml_file',trim(yaml_error))
      if (.not.associated(node)) call driver%fatal_error('fabm_create_model_from_yaml_file','No configuration information found in '//trim(path)//'.')
      !call node%dump(output_unit,0)

      ! Create model tree from YAML root node.
      select type (node)
         class is (type_dictionary)
            ! Create F2003 model tree.
            call create_model_tree_from_dictionary(model,node,do_not_initialize,parameters)
         class default
            call driver%fatal_error('fabm_create_model_from_yaml_file',trim(path)//' must contain "MODEL:" entries at the root (non-indented) level, not a single value. Are you missing a trailing colon?')
      end select

   end subroutine fabm_create_model_from_yaml_file

   subroutine create_model_tree_from_dictionary(model,mapping,do_not_initialize,parameters)
      type (type_model),                       intent(out) :: model
      class (type_dictionary),                 intent(in)  :: mapping
      logical,                        optional,intent(in)  :: do_not_initialize
      type (type_property_dictionary),optional,intent(in)  :: parameters

      class (type_node),pointer          :: node
      character(len=64)                  :: instancename
      type (type_key_value_pair),pointer :: pair
      class (type_base_model),   pointer :: childmodel
      logical                            :: initialize,check_conservation,success

      ! If custom parameter values were provided, transfer these to the root model.
      if (present(parameters)) call model%root%parameters%update(parameters)

      check_conservation = .false.
      node => mapping%get('check_conservation')
      if (associated(node)) then
         select type (node)
            class is (type_scalar)
               check_conservation = node%to_logical(.false.,success)
               if (.not.success) call driver%fatal_error('create_model_tree_from_dictionary','/check_conservation is set to "'//trim(node%string)//'", which cannot be interpreted as a Boolean value.')
            class default
               call driver%fatal_error('create_model_tree_from_dictionary','/check_conservation must be a key, not a dictionary.')
         end select
      end if

      node => mapping%get('models')
      if (.not.associated(node)) call driver%fatal_error('create_model_tree_from_dictionary','No "models" dictionary found at root level.')
      select type (node)
         class is (type_dictionary)
            pair => node%first
         class default
            nullify(pair)
            call driver%fatal_error('create_model_tree_from_dictionary','"/models" must be a dictionary with (model name : information) pairs, not a single value.')
      end select

      ! Iterate over all models (key:value pairs at root level) and
      ! create corresponding objects.
      do while (associated(pair))
         instancename = trim(pair%key)
         select type (dict=>pair%value)
            class is (type_dictionary)
               childmodel => create_model_from_dictionary(instancename,dict,model%root,require_initialization=.false.)
               childmodel%check_conservation = check_conservation
            class default
               call driver%fatal_error('create_model_tree_from_dictionary','Configuration information for model "'//trim(instancename)//'" must be a dictionary, not a single value.')
         end select
         pair => pair%next
      end do

      ! Check whether any keys at the root level remain unused.
      pair => mapping%first
      do while (associated(pair))
         if (.not.pair%accessed) call driver%fatal_error('create_model_tree_from_dictionary','Unrecognized option "'//trim(pair%key)//'" found at root level.')
         pair => pair%next
      end do

      ! Initialize model tree
      initialize = .true.
      if (present(do_not_initialize)) initialize = .not.do_not_initialize
      if (initialize) call fabm_initialize(model)

   end subroutine create_model_tree_from_dictionary

   function create_model_from_dictionary(instancename,node,parent,require_initialization) result(model)
      character(len=*),       intent(in)           :: instancename
      class (type_dictionary),intent(in)           :: node
      class (type_base_model),intent(inout),target :: parent
      logical,                intent(in)           :: require_initialization
      class (type_base_model),pointer              :: model

      character(len=64)                  :: modelname
      character(len=256)                 :: long_name
      type (type_dictionary)             :: parametermap
      class (type_node),pointer          :: childnode
      class (type_property),pointer      :: property
      type (type_key_value_pair),pointer :: pair
      type (type_set)                    :: initialized_set
      type (type_link),pointer           :: link
      logical                            :: success

      ! Retrieve model name (default to instance name if not provided).
      modelname = trim(node%get_string('model',default=instancename))

      ! Retrieve descriptive name for the model instance (default to instance name if not provided).
      long_name = trim(node%get_string('long_name',default=instancename))

      ! If FABM does not have its model factory yet, create it.
      call fabm_create_model_factory()

      ! Try to create the model based on name.
      call factory%create(trim(modelname),model)
      if (.not.associated(model)) call driver%fatal_error('create_model_from_dictionary',trim(instancename)//': "'//trim(modelname)//'" is not a valid model name.')

      ! Transfer user-specified parameter values to the model.
      childnode => node%get('parameters')
      if (associated(childnode)) then
         select type (childnode)
            class is (type_dictionary)
               call childnode%flatten(parametermap,'')
               pair => parametermap%first
               do while (associated(pair))
                  select type (value=>pair%value)
                     class is (type_scalar)
                        call model%parameters%set_string(trim(pair%key),trim(value%string))
                     class default
                        call driver%fatal_error('create_model_from_dictionary','BUG: "flatten" should have ensured that the value of /models/'//trim(instancename)//'/parameters/'//trim(pair%key)//' is scalar, not a nested dictionary.')
                  end select
                  pair => pair%next   
               end do
            class default
               call driver%fatal_error('create_model_from_dictionary','/models/'//trim(instancename)//'/parameters must be a dictionary, not a string.')
         end select
      end if

      ! Add the model to its parent.
      call model%parameters%reset_accessed()
      call driver%log_message('Initializing biogeochemical model "'//trim(instancename)//'" (type "'//trim(modelname)//'")...')
      call parent%add_child(model,instancename,long_name,configunit=-1)
      call driver%log_message('model "'//trim(instancename)//'" initialized successfully.')

      ! Check for parameters present in configuration file, but not interpreted by the models.
      property => model%parameters%first
      do while (associated(property))
         if (.not.property%accessed) call driver%fatal_error('create_model_from_dictionary','Unrecognized parameter '//trim(property%name)//' found below /models/'//trim(instancename)//'.')
         property => property%next
      end do

      ! Interpret coupling links specified in configuration file.
      ! These override any couplings requested by the models during initialization.
      ! This step must therefore occur after model initialization [from parent%add_child] has completed.
      childnode => node%get('coupling')
      if (associated(childnode)) then
         select type (childnode)
            class is (type_dictionary)
               pair => childnode%first
               do while (associated(pair))
                  select type (value=>pair%value)
                     class is (type_scalar)
                        call model%request_coupling(trim(get_safe_name(pair%key)),trim(get_safe_name(value%string)),required=.true.)
                     class default
                        call driver%fatal_error('create_model_from_dictionary','The value of '//trim(instancename)//'/coupling/'//trim(pair%key)//' must be a string, not a nested dictionary.')
                  end select
                  pair => pair%next   
               end do
            class default
               call driver%fatal_error('create_model_from_dictionary','/models/'//trim(instancename)//'/coupling must be a dictionary, not a single value.')
         end select
      end if

      ! Transfer user-specified initial state to the model.
      childnode => node%get('initialization')
      if (associated(childnode)) then
         select type (childnode)
            class is (type_dictionary)
               call parse_initialization(model,childnode,initialized_set)
            class default
               call driver%fatal_error('create_model_from_dictionary','/models/'//trim(instancename)//'/parameters must be a dictionary, not a string.')
         end select
      end if

      ! Verify whether all state variables have been provided with an initial value.
      link => model%first_link
      do while (associated(link))
         if (link%owner.and.link%coupling%master=='') then
            ! This link is our own: not coupled to another variable, not an alias, and no intention to couple was registered.
            select type (object=>link%target)
               class is (type_internal_variable)
                  if (allocated(object%state_indices)) then
                     if (object%presence/=presence_external_optional .and. .not.initialized_set%contains(trim(link%name))) then
                        ! State variable is required, but initial value is not explicitly provided.
                        if (require_initialization) then
                           call driver%fatal_error('parse_initialization','/models/'//trim(model%name) &
                              //': no initial value provided for variable "'//trim(link%name)//'".')
                        else
                           call driver%log_message('WARNING: no explicit initial value provided for "'//trim(model%name)//'/'//trim(link%name)//'" - using default.')
                        end if
                     elseif (object%presence==presence_external_optional .and. initialized_set%contains(trim(link%name))) then
                        ! Optional state variable is not used, but an initial value is explicitly provided.
                        call driver%fatal_error('parse_initialization','/models/'//trim(model%name) &
                           //': initial value provided for variable "'//trim(link%name)//'", but this variable is not used.')
                     end if
                  end if
            end select
         end if
         link => link%next
      end do

      childnode => node%get('check_conservation')
      if (associated(childnode)) then
         select type (childnode)
            class is (type_scalar)
               model%check_conservation = childnode%to_logical(.false.,success)
               if (.not.success) call driver%fatal_error('create_model_from_dictionary','/models/'//trim(instancename)//'/check_conservation is set to "'//trim(childnode%string)//'", which cannot be interpreted as a Boolean value.')
            class default
               call driver%fatal_error('create_model_from_dictionary','/models/'//trim(instancename)//'/check_conservation must be a key, not a dictionary.')
         end select
      end if

      ! Check whether any keys at the model level remain unused.
      pair => node%first
      do while (associated(pair))
         if (.not.pair%accessed) call driver%fatal_error('create_model_from_dictionary','Unrecognized option "'//trim(pair%key)//'" found at /models/'//trim(instancename)//'.')
         pair => pair%next
      end do

   end function create_model_from_dictionary

   subroutine parse_initialization(model,node,initialized_set)
      class (type_base_model),intent(inout) :: model
      class (type_dictionary),intent(in)    :: node
      type (type_set),        intent(out)   :: initialized_set

      type (type_key_value_pair),  pointer :: pair
      class (type_internal_object),pointer :: object
      logical                              :: is_state_variable

      ! Transfer user-specified initial state to the model.
      pair => node%first
      do while (associated(pair))
         select type (value=>pair%value)
            class is (type_scalar)
               object => model%find_object(trim(pair%key))
               if (.not.associated(object)) call driver%fatal_error('parse_initialization','/models/'//trim(model%name)//'/initialization/'//trim(pair%key)//': "'//trim(pair%key)//'" is not a member of model "'//trim(model%name)//'".')
               is_state_variable = .false.
               select type (object)
                  class is (type_internal_variable)
                     if (allocated(object%state_indices)) then
                        read (value%string,*,err=99,end=99) object%initial_value
                        call initialized_set%add(trim(pair%key))
                        is_state_variable = .true.
                     end if
               end select
               if (.not.is_state_variable) call driver%fatal_error('parse_initialization','/models/'//trim(model%name)//'/initialization/'//trim(pair%key)//': "'//trim(pair%key)//'" is not a state variable of model "'//trim(model%name)//'".')
            class default
               call driver%fatal_error('parse_initialization','/models/'//trim(model%name)//'/initialization/'//trim(pair%key)//' must be set to a scalar value, not to a nested dictionary.')
         end select
         pair => pair%next   
      end do

      return

99    select type (value=>pair%value)
         class is (type_scalar)
            call driver%fatal_error('parse_initialization','/models/'//trim(model%name)//'/initialization/'//trim(pair%key)//': "'//trim(value%string)//'" is not a real number.')
      end select
   end subroutine

end module fabm_config

#endif

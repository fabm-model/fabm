! This module defined a new model type, type_particle_model, which may be used as base
! model type for other biogeochemical models. This model type is seen as a base type for models
! that can be considered "integral physical particles" in the context of Lagrangian models.
! For instance: a phytoplankton species or type, a zooplankton species or type, a specific
! (size) class of organic matter, etc.
!
! The new model type inherits from type_base_model, and adds functionality for accessing
! other active models by declaring a "model dependency", which may be coupled to an actual model
! at run-time just as variables are.
! By registering a model dependency (self%register_model_dependency) and using the returned model
! identifier, model variables can be coupled to named variables of the other model
! (self%request_coupling_to_model) and a list of the identifiers of all state variables in the other model
! can be obtained (model_id%state). In turn, those identifiers can be used to e.g. prescribe a
! specific rate of change to all state variables of the other model.
!
! This functionality enables model-to-model coupling in addition to variable-to-variable
! coupling. Model-to-model coupling is convenient when each model features state variables
! that always feature together in biogeochemical processes (i.e., the process operates on
! all of the model's state variables simultaneously).
!
! For instance, a phytoplankton model could serve as prey for a zooplankton model,
! in which case zooplankton predation operates on all state variables (e.g., all chemical constituents)
! of the phytoplankton model. This type of interaction can be achieved by declaring a "prey" model
! dependency in the zooplankton model, and linking that at run-time to the phytoplankton model.
! The zooplankton model can then directly query the variables of the prey model (e.g., total carbon),
! and subsequently apply a specific loss rate to all of the prey model's state variables.

module fabm_particle

   use fabm_types

   implicit none

   private

   public type_particle_model, type_model_id

   integer, parameter :: not_done = 0, busy = 1, done   = 2

   type type_model_reference
      logical                                 :: resolving = .false.
      character(len=attribute_length)         :: name  = ''
      character(len=attribute_length)         :: long_name  = ''
      class (type_base_model),        pointer :: model => null()
      type (type_model_id),           pointer :: id    => null()
      type (type_model_reference),    pointer :: next  => null()
   end type

   type type_model_id
      type (type_model_reference), pointer              :: reference => null()
      type (type_state_variable_id),        allocatable :: state(:)
      type (type_bottom_state_variable_id), allocatable :: bottom_state(:)
      type (type_surface_state_variable_id),allocatable :: surface_state(:)
   end type

   type, extends(type_coupling_task) :: type_coupling_from_model
      type (type_model_reference), pointer :: model_reference => null()
      integer                              :: access          = access_read
      class (type_particle_model), pointer :: owner           => null()
   contains
      procedure :: resolve => coupling_from_model_resolve
   end type

   type, extends(type_base_model) :: type_particle_model
      type (type_model_reference), pointer :: first_model_reference   => null()
      integer                              :: internal_variable_state = not_done
   contains
      ! Used by inheriting biogeochemical models:
      procedure :: register_model_dependency1
      procedure :: register_model_dependency2
      generic :: register_model_dependency => register_model_dependency1, register_model_dependency2

      procedure :: request_named_coupling_to_model
      procedure :: request_standard_coupling_to_model
      procedure :: request_named_coupling_to_named_model
      procedure :: request_standard_coupling_to_named_model
      procedure :: request_coupling_to_model_generic
      generic :: request_coupling_to_model => request_named_coupling_to_model,request_standard_coupling_to_model, &
                                              request_named_coupling_to_named_model,request_standard_coupling_to_named_model
      procedure :: resolve_model_dependency
      procedure :: request_coupling_nn

      ! Hooks called by FABM:
      procedure :: before_coupling
      procedure :: resolve_model_reference
      procedure :: complete_internal_variables

      procedure :: finalize
   end type

contains

   subroutine register_model_dependency1(self, id, name, long_name)
      class (type_particle_model),  intent(inout)           :: self
      type (type_model_id), target, intent(inout)           :: id
      character(len=*),             intent(in)              :: name
      character(len=*), optional,   intent(in)              :: long_name

      type (type_model_reference), pointer :: reference

      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_model_dependency',&
         'Dependencies may only be registered during initialization.')

      reference => add_model_reference(self, name, long_name, require_empty_id=.true.)
      reference%id => id
      id%reference => reference
   end subroutine register_model_dependency1

   subroutine register_model_dependency2(self, name, long_name)
      class (type_particle_model),  intent(inout) :: self
      character(len=*),             intent(in)    :: name
      character(len=*), optional,   intent(in)    :: long_name

      type (type_model_reference), pointer :: reference

      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_model_dependency',&
         'Dependencies may only be registered during initialization.')

      reference => add_model_reference(self, name, long_name)
   end subroutine register_model_dependency2

   function add_model_reference(self, name, long_name, require_empty_id) result(reference)
      class (type_particle_model), intent(inout) :: self
      character(len=*),            intent(in)    :: name
      character(len=*), optional,  intent(in)    :: long_name
      logical,          optional,  intent(in)    :: require_empty_id

      type (type_model_reference), pointer :: reference

      ! First search for existing model reference.
      reference => find_model_reference(self, name, require_empty_id)
      if (associated(reference)) return

      ! Register this model reference.
      if (.not. associated(self%first_model_reference)) then
         ! No model references yet - register the first.
         allocate(self%first_model_reference)
         reference => self%first_model_reference
      else
         ! One or more model references exist - append to the linked list.
         reference => self%first_model_reference
         do while (associated(reference%next))
            reference => reference%next
         end do
         allocate(reference%next)
         reference => reference%next
      end if

      ! Set attributes of the model reference.
      reference%name = name
      if (present(long_name)) then
         reference%long_name = long_name
      else
         reference%long_name = name
      end if
   end function add_model_reference

   function find_model_reference(self, name, require_empty_id) result(reference)
      class (type_particle_model), intent(inout) :: self
      character(len=*),            intent(in)    :: name
      logical,          optional,  intent(in)    :: require_empty_id

      logical                              :: require_empty_id_
      type (type_model_reference), pointer :: reference

      require_empty_id_ = .false.
      if (present(require_empty_id)) require_empty_id_ = require_empty_id

      reference => self%first_model_reference
      do while (associated(reference))
         if (reference%name == name .and. .not. (associated(reference%id) .and. require_empty_id_)) return
         reference => reference%next
      end do
   end function find_model_reference

   function resolve_model_dependency(self, name, require_internal_variables) result(model)
      class (type_particle_model), intent(inout) :: self
      character(len=*),            intent(in)    :: name
      logical, optional,           intent(in)    :: require_internal_variables
      class (type_base_model), pointer :: model

      type (type_model_reference), pointer :: reference

      model => null()
      reference => find_model_reference(self, name)
      if (associated(reference)) model => self%resolve_model_reference(reference, require_internal_variables)
   end function resolve_model_dependency

   recursive subroutine request_coupling_nn(self, name, target_name)
      class (type_particle_model), intent(inout), target :: self
      character(len=*),            intent(in)            :: name, target_name

      type (type_model_reference), pointer :: reference

      reference => find_model_reference(self, name)
      if (associated(reference)) then
         ! The provided name is a local model reference that needs to be coupled to another model instance
         call self%couplings%set_string(name, target_name)
      else
         ! The provided name does not match that of a local model reference.
         ! It therefore must be a variable. Forward to base type
         call self%type_base_model%request_coupling_nn(name, target_name)
      end if
   end subroutine request_coupling_nn

   subroutine request_coupling_to_model_generic(self, link, target_model, target_model_name, target_name, target_standard_variable)
      class (type_particle_model), target, intent(inout)           :: self
      type (type_link),target                                      :: link
      type (type_model_id),                intent(inout), optional :: target_model
      character(len=*),                    intent(in),    optional :: target_name, target_model_name
      class (type_base_standard_variable), intent(in),    optional :: target_standard_variable

      type (type_coupling_from_model),     pointer :: coupling
      class (type_coupling_task),          pointer :: base_coupling
      class (type_base_standard_variable), pointer :: standard_variable

      if (self%coupling_task_list%includes_custom) call self%fatal_error('request_coupling_to_model_generic', &
         'Too late to add custom couplings (of type type_coupling_from_model')

      ! Create object describing the coupling, and send it to FABM.
      ! This must be a pointer, because FABM will manage its memory and deallocate when appropriate.
      allocate(coupling)
      base_coupling => coupling
      call self%request_coupling(link, base_coupling)
      if (.not. associated(base_coupling)) return
      coupling%owner => self
      if (present(target_name)) coupling%target_name = target_name
      if (present(target_standard_variable)) then
         standard_variable => target_standard_variable%resolve()
         select type (standard_variable)
         class is (type_universal_standard_variable)
            select case (coupling%link%target%domain)
            case (domain_interior);   coupling%target_standard_variable => standard_variable%in_interior()
            case (domain_surface);    coupling%target_standard_variable => standard_variable%at_surface()
            case (domain_bottom);     coupling%target_standard_variable => standard_variable%at_bottom()
            case (domain_horizontal); coupling%target_standard_variable => standard_variable%at_interfaces()
            end select
         class is (type_domain_specific_standard_variable)
            coupling%target_standard_variable => standard_variable
         end select
         if (coupling%link%target%source == source_state .or. coupling%link%target%fake_state_variable) &
            coupling%access = access_state
      end if
      if (present(target_model)) then
         coupling%model_reference => target_model%reference
      else
         if (.not. present(target_model_name)) call self%fatal_error('request_coupling_to_model_generic', &
            'BUG: target_model or target_model_name must be provided.')
         coupling%model_reference => add_model_reference(self, target_model_name)
      end if
   end subroutine request_coupling_to_model_generic

   subroutine request_named_coupling_to_model(self, id, target_model, target_name)
      class (type_particle_model), target, intent(inout) :: self
      class (type_variable_id),            intent(in)    :: id
      type (type_model_id),                intent(inout) :: target_model
      character(len=*),                    intent(in)    :: target_name

      if (.not. associated(id%link)) &
         call self%fatal_error('request_named_coupling_to_model', 'variable must be registered before it can be coupled.')
      call request_coupling_to_model_generic(self, id%link, target_model=target_model, target_name=target_name)
   end subroutine request_named_coupling_to_model

   subroutine request_standard_coupling_to_model(self, id, target_model, target_standard_variable)
      class (type_particle_model), target, intent(inout) :: self
      class (type_variable_id),            intent(in)    :: id
      type (type_model_id),                intent(inout) :: target_model
      class (type_base_standard_variable), intent(in)    :: target_standard_variable

      if (.not. associated(id%link)) &
         call self%fatal_error('request_standard_coupling_to_model', 'variable must be registered before it can be coupled.')
      call request_coupling_to_model_generic(self, id%link, target_model=target_model, target_standard_variable=target_standard_variable)
   end subroutine request_standard_coupling_to_model

   subroutine request_named_coupling_to_named_model(self, id, target_model_name, target_name)
      class (type_particle_model), target, intent(inout) :: self
      class (type_variable_id),            intent(in)    :: id
      character(len=*),                    intent(in)    :: target_model_name
      character(len=*),                    intent(in)    :: target_name

      if (.not. associated(id%link)) &
         call self%fatal_error('request_named_coupling_to_named_model', 'variable must be registered before it can be coupled.')
      call request_coupling_to_model_generic(self, id%link, target_model_name=target_model_name, target_name=target_name)
   end subroutine request_named_coupling_to_named_model

   subroutine request_standard_coupling_to_named_model(self, id, target_model_name, target_standard_variable)
      class (type_particle_model), target, intent(inout) :: self
      class (type_variable_id),            intent(in)    :: id
      character(len=*),                    intent(in)    :: target_model_name
      class (type_base_standard_variable), intent(in)    :: target_standard_variable

      if (.not. associated(id%link)) &
         call self%fatal_error('request_named_coupling_to_named_model', 'variable must be registered before it can be coupled.')
      call request_coupling_to_model_generic(self, id%link, target_model_name=target_model_name, target_standard_variable=target_standard_variable)
   end subroutine request_standard_coupling_to_named_model

   recursive function resolve_model_reference(self, reference, require_internal_variables) result(model)
      class (type_particle_model), intent(inout),target :: self
      type (type_model_reference), intent(inout)        :: reference
      logical, optional,           intent(in)           :: require_internal_variables
      class (type_base_model), pointer :: model

      type (type_model_reference), pointer :: reference2
      character(len=:), allocatable        :: model_target_name
      integer                              :: istart
      class (type_particle_model), pointer :: source_model
      logical                              :: require_internal_variables_

      ! If the model reference has been resolved before, we do not return immediately but continue
      ! below in case this is the first time that require_internal_variables was set.

      if (.not. associated(reference%model)) then
         ! If this model reference is *already* in the process of being resolved,
         ! that can only be because it is pointing to itself, either directly or indirectly.
         if (reference%resolving) then
            call self%fatal_error('resolve_model_reference', 'Circular dependency')
            return
         end if

         ! Flag this references as "in the process of being resolved",
         ! in order to support detection of circular dependencies
         reference%resolving = .true.

         ! First find a coupling for the referenced model.
         model_target_name = self%couplings%get_string(trim(reference%name), trim(reference%long_name))

         ! Try to find referenced model among actual models (as opposed to among model dependencies)
         reference%model => self%find_model(model_target_name, recursive=.true.)

         if (.not. associated(reference%model)) then
            ! Search among named model dependencies instead

            ! Find starting position of local name (excluding any preprended path components)
            istart = index(model_target_name, '/', .true.) + 1

            if (istart == 1) then
               ! No slash in path; search model references within current model
               source_model => self
            else
               ! Try model references in specified other model
               source_model => null()
               model => self%find_model(model_target_name(:istart - 1), recursive=.true.)
               if (associated(model)) then
                  select type (model)
                  class is (type_particle_model)
                     source_model => model
                  end select
               end if
            end if

            ! Search references in source model
            if (associated(source_model)) then
               reference2 => find_model_reference(source_model, model_target_name(istart:))
               if (associated(reference2)) reference%model => source_model%resolve_model_reference(reference2)
            end if
         end if

         if (.not. associated(reference%model)) then
            call self%fatal_error('resolve_model_reference', &
               'Referenced model instance "' // model_target_name // '" not found.')
            return
         end if

         reference%resolving = .false.
      end if

      model => reference%model

      require_internal_variables_ = .false.
      if (present(require_internal_variables)) require_internal_variables_ = require_internal_variables
      if (require_internal_variables_) then
         select type (model)
         class is (type_particle_model)
            call complete_internal_variables_if_needed(model)
         end select
      end if
   end function resolve_model_reference

   recursive subroutine complete_internal_variables_if_needed(self)
      class (type_particle_model), intent(inout) :: self

      select case (self%internal_variable_state)
      case (done)
         return
      case (busy)
         call self%fatal_error('complete_internal_variables_if_needed', 'Circular dependency')
      end select
      self%internal_variable_state = busy
      call self%complete_internal_variables()
      self%internal_variable_state = done
   end subroutine complete_internal_variables_if_needed

   recursive subroutine complete_internal_variables(self)
      class (type_particle_model), intent(inout) :: self
   end subroutine

   recursive subroutine before_coupling(self)
      class (type_particle_model), intent(inout) :: self

      type (type_model_reference), pointer :: reference
      class (type_base_model),     pointer :: model

      ! For model references that include the model state, buid the corresponding lists
      ! of state variable identifiers now, and request their coupling to the referenced model.
      ! This must be done before the actual [variable] coupling starts, and therefore has
      ! to happen here, in before_coupling.
      reference => self%first_model_reference
      do while (associated(reference))
         if (associated(reference%id)) then
            model => resolve_model_reference(self, reference, require_internal_variables=.true.)
            if (.not. associated(model)) return
            call build_state_id_list(self, reference, domain_interior)
            call build_state_id_list(self, reference, domain_surface)
            call build_state_id_list(self, reference, domain_bottom)
         end if
         reference => reference%next
      end do

      call complete_internal_variables_if_needed(self)
   end subroutine before_coupling

   function coupling_from_model_resolve(self) result(link)
      class (type_coupling_from_model), intent(inout) :: self
      type (type_link), pointer :: link

      class (type_base_model), pointer :: model

      model => resolve_model_reference(self%owner, self%model_reference)
      if (.not. associated(model)) then
         ! Model not found. A fatal error will already have been reported.
         ! Just return a harmless result so the host gets the opportunity for error handling
         link => null()
      elseif (associated(self%target_standard_variable)) then
         ! Coupling to a standard [aggregate] variable
         link => get_aggregate_variable_access(model, self%target_standard_variable, self%access)
      else
         ! Coupling to a named variable
         link => model%find_link(trim(self%target_name))
      end if
   end function

   subroutine build_state_id_list(self, reference, domain)
      class (type_particle_model), intent(inout) :: self
      type (type_model_reference), intent(inout) :: reference
      integer,                     intent(in)    :: domain

      type (type_link), pointer :: link
      integer                   :: n
      character(len=10)         :: strindex

      ! Count number of state variables in target model.
      n = 0
      link => reference%model%links%first
      do while (associated(link))
         if (index(link%name, '/') == 0 .and. link%target%domain == domain .and. link%original%presence == presence_internal &
             .and.(link%original%source == source_state .or. link%original%fake_state_variable)) n = n + 1
         link => link%next
      end do

      ! Allocate array to hold state variable identifiers.
      select case (domain)
      case (domain_interior)
         allocate(reference%id%state(n))
      case (domain_bottom)
         allocate(reference%id%bottom_state(n))
      case (domain_surface)
         allocate(reference%id%surface_state(n))
      end select

      ! Connect target state variables to identifiers in model reference.
      n = 0
      link => reference%model%links%first
      do while (associated(link))
         if (index(link%name,'/') == 0 .and. link%target%domain == domain .and. link%original%presence == presence_internal &
             .and. (link%original%source == source_state .or. link%original%fake_state_variable)) then
            n = n + 1
            write (strindex,'(i0)') n
            select case (domain)
            case (domain_interior)
               call self%register_state_dependency(reference%id%state(n), trim(reference%name) // '_state' // trim(strindex), &
                  link%target%units, trim(reference%name) // ' state variable ' // trim(strindex))
               call self%request_coupling(reference%id%state(n), link)
            case (domain_bottom)
               call self%register_state_dependency(reference%id%bottom_state(n), trim(reference%name) // '_bottom_state' // trim(strindex), &
                  link%target%units, trim(reference%name) // ' bottom state variable ' // trim(strindex))
               call self%request_coupling(reference%id%bottom_state(n), link)
            case (domain_surface)
               call self%register_state_dependency(reference%id%surface_state(n), trim(reference%name) // '_surface_state' // trim(strindex), &
                  link%target%units, trim(reference%name) // ' surface state variable ' // trim(strindex))
               call self%request_coupling(reference%id%surface_state(n), link)
            end select
         end if
         link => link%next
      end do

   end subroutine build_state_id_list

   recursive subroutine finalize(self)
      class (type_particle_model),  intent(inout) :: self

      type (type_model_reference), pointer :: reference, next_reference

      reference => self%first_model_reference
      do while (associated(reference))
         next_reference => reference%next
         deallocate(reference)
         reference => next_reference
      end do
      self%first_model_reference => null()

      call self%type_base_model%finalize()
   end subroutine

end module

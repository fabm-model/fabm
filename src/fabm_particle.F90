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
   use fabm_standard_variables
   use fabm_properties

   implicit none

   private

   public type_particle_model, type_model_id

   type type_model_reference
      character(len=attribute_length)         :: name  = ''
      class (type_base_model),        pointer :: model => null()
      type (type_model_id),           pointer :: id    => null()
      type (type_model_reference),    pointer :: next  => null()
   end type

   type type_model_id
      type (type_model_reference),pointer               :: reference => null()
      type (type_state_variable_id),        allocatable :: state(:)
      type (type_bottom_state_variable_id), allocatable :: bottom_state(:)
      type (type_surface_state_variable_id),allocatable :: surface_state(:)
   end type

   type,extends(type_coupling_task) :: type_coupling_from_model
      type (type_model_reference), pointer :: model_reference => null()
      integer                              :: access          = access_read
   end type

   type,extends(type_base_model) :: type_particle_model
      type (type_model_reference),    pointer :: first_model_reference => null()
   contains
      ! Used by inheriting biogeochemical models:
      procedure :: register_model_dependency
      procedure :: request_named_coupling_to_model
      procedure :: request_standard_coupling_to_model
      procedure :: request_named_coupling_to_named_model
      procedure :: request_standard_coupling_to_named_model
      generic :: request_coupling_to_model => request_named_coupling_to_model,request_standard_coupling_to_model, &
                                              request_named_coupling_to_named_model,request_standard_coupling_to_named_model

      ! Hooks called by FABM:
      procedure :: before_coupling
   end type

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Registers a dependency on another model.
!
! !INTERFACE:
   subroutine register_model_dependency(self,id,name)
!
! !DESCRIPTION:
!  This subroutine registers a dependency on a named external model.
!  It initializes the caller-supplied identifier, which may subsequently be
!  used to request variable coupling with the referenced model, or to
!  access all of the model's state variables.
!
! !INPUT/OUTPUT PARAMETERS:
      class (type_particle_model), intent(inout) :: self
      type (type_model_id),target, intent(inout) :: id
!
! !INPUT PARAMETERS:
      character(len=*),            intent(in)    :: name
!
!EOP
!
! !LOCAL VARIABLES:
      type (type_model_reference), pointer :: reference
!
!-----------------------------------------------------------------------
!BOC
      ! Check whether the model information may be written to (only during initialization)
      if (self%frozen) call self%fatal_error('register_model_dependency',&
         'Dependencies may only be registered during initialization.')

      reference => add_model_reference(self,name,require_empty_id=.true.)
      reference%id => id
      id%reference => reference

   end subroutine register_model_dependency
!EOC

   function add_model_reference(self,name,require_empty_id) result(reference)
      class (type_particle_model), intent(inout) :: self
      character(len=*),            intent(in)    :: name
      logical,                     intent(in)    :: require_empty_id

      type (type_model_reference), pointer :: reference

      ! First search for existing model reference.
      reference => self%first_model_reference
      do while (associated(reference))
         if (reference%name==name.and..not.(associated(reference%id).and.require_empty_id)) return
         reference => reference%next
      end do

      ! Register this model reference.
      if (.not.associated(self%first_model_reference)) then
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
   end function add_model_reference

   subroutine request_coupling_to_model_generic(self,slave_variable,master_model,master_model_name,master_name,master_standard_variable)
      class (type_particle_model),       intent(inout)          :: self
      class (type_variable_id),          intent(in)             :: slave_variable
      class (type_model_id),             intent(inout),optional :: master_model
      character(len=*),                  intent(in),   optional :: master_name, master_model_name
      type (type_bulk_standard_variable),intent(in),   optional :: master_standard_variable

      type (type_coupling_from_model), pointer :: coupling
      class (type_coupling_task),      pointer :: base_coupling

      if (.not.associated(slave_variable%link)) &
         call self%fatal_error('request_coupling_to_model_generic','slave variable must be registered before it is coupled.')

      ! Create object describing the coupling, and send it to FABM.
      ! This must be a pointer, because FABM will manage its memory and deallocate when appropriate.
      allocate(coupling)
      coupling%slave => slave_variable%link
      if (present(master_name)) coupling%master_name = master_name
      if (present(master_standard_variable)) then
         allocate(coupling%master_standard_variable,source=master_standard_variable)
         select type (slave_variable)
         class is (type_dependency_id)
            coupling%domain = domain_interior
         class is (type_horizontal_dependency_id)
            coupling%domain = domain_horizontal
         class is (type_state_variable_id)
            coupling%domain = domain_interior
            coupling%access = access_state
         class is (type_bottom_state_variable_id)
            coupling%domain = domain_bottom
            coupling%access = access_state
         class is (type_surface_state_variable_id)
            coupling%domain = domain_surface
            coupling%access = access_state
         class default
            call self%fatal_error('request_coupling_to_model_sn','Provided variable id must be of ones of the following types: &
               &type type_dependency_id, type_horizontal_dependency_id, type_state_variable_id, type_bottom_state_variable_id, type_surface_state_variable_id.')
         end select
      end if
      if (present(master_model)) then
         coupling%model_reference => master_model%reference
      else
         if (.not.present(master_model_name)) call self%fatal_error('request_coupling_to_model_generic','BUG: master_model or master_model_name must be provided.')
         coupling%model_reference => add_model_reference(self,master_model_name,require_empty_id=.false.)
      end if
      base_coupling => coupling
      if (.not.self%coupling_task_list%add(base_coupling,.false.)) deallocate(coupling)
   end subroutine request_coupling_to_model_generic

   subroutine request_named_coupling_to_model(self,slave_variable,master_model,master_variable)
      class (type_particle_model),intent(inout) :: self
      class (type_variable_id),   intent(in)    :: slave_variable
      class (type_model_id),      intent(inout) :: master_model
      character(len=*),           intent(in)    :: master_variable

      call request_coupling_to_model_generic(self,slave_variable,master_model=master_model,master_name=master_variable)
   end subroutine request_named_coupling_to_model

   subroutine request_standard_coupling_to_model(self,slave_variable,master_model,master_variable)
      class (type_particle_model),       intent(inout) :: self
      class (type_variable_id),          intent(in)    :: slave_variable
      class (type_model_id),             intent(inout) :: master_model
      type (type_bulk_standard_variable),intent(in)    :: master_variable

      call request_coupling_to_model_generic(self,slave_variable,master_model=master_model,master_standard_variable=master_variable)
   end subroutine request_standard_coupling_to_model

   subroutine request_named_coupling_to_named_model(self,slave_variable,master_model,master_variable)
      class (type_particle_model),intent(inout) :: self
      class (type_variable_id),   intent(in)    :: slave_variable
      character(len=*),           intent(in)    :: master_model
      character(len=*),           intent(in)    :: master_variable

      call request_coupling_to_model_generic(self,slave_variable,master_model_name=master_model,master_name=master_variable)
   end subroutine request_named_coupling_to_named_model

   subroutine request_standard_coupling_to_named_model(self,slave_variable,master_model,master_variable)
      class (type_particle_model),       intent(inout) :: self
      class (type_variable_id),          intent(in)    :: slave_variable
      character(len=*),                  intent(in)    :: master_model
      type (type_bulk_standard_variable),intent(in)    :: master_variable

      call request_coupling_to_model_generic(self,slave_variable,master_model_name=master_model,master_standard_variable=master_variable)
   end subroutine request_standard_coupling_to_named_model

   subroutine before_coupling(self)
      class (type_particle_model),intent(inout) :: self

      type (type_model_reference),           pointer :: reference,reference2
      class (type_property),                 pointer :: model_master_name
      type (type_aggregate_variable_access), pointer :: aggregate_variable_access
      class (type_coupling_task),            pointer :: coupling, next_coupling
      character(len=attribute_length)                :: master_name

      reference => self%first_model_reference
      do while (associated(reference))
         ! First find a coupling for the referenced model.
         model_master_name => self%couplings%find_in_tree(reference%name)
         if (.not.associated(model_master_name)) call self%fatal_error('before_coupling','Model reference "'//trim(reference%name)//'" was not coupled.')

         ! First try among earlier model references
         reference2 => self%first_model_reference
         do while (associated(reference2))
            select type (model_master_name)
               class is (type_string_property)
                  if (model_master_name%value==reference2%name) reference%model => reference2%model
            end select
            reference2 => reference2%next
         end do

         if (.not.associated(reference%model)) then
            ! Now find the actual model that was coupled.
            select type (model_master_name)
               class is (type_string_property)
                  reference%model => self%find_model(model_master_name%value,recursive=.true.)
                  if (.not.associated(reference%model)) call self%fatal_error('before_coupling','Referenced model "'//trim(model_master_name%value)//'" not found.')
            end select
         end if

         reference => reference%next
      end do

      reference => self%first_model_reference
      do while (associated(reference))
         if (associated(reference%id)) then
            call build_state_id_list(self,reference,domain_interior)
            call build_state_id_list(self,reference,domain_surface)
            call build_state_id_list(self,reference,domain_bottom)
         end if
         reference => reference%next
      end do

      ! Resolve all couplings to variables in specific other models.
      coupling => self%coupling_task_list%first
      do while (associated(coupling))
         next_coupling => coupling%next
         select type (model_coupling=>coupling)
         class is (type_coupling_from_model)
            if (coupling%master_name/='') then
               ! Coupling to a named variable
               master_name = trim(model_coupling%model_reference%model%get_path())//'/'//trim(coupling%master_name)
            else
               ! Coupling to a standard [aggregate] variable
               select type (standard_variable=>coupling%master_standard_variable)
               class is (type_bulk_standard_variable)
                  aggregate_variable_access => get_aggregate_variable_access(model_coupling%model_reference%model,standard_variable)
                  select case (coupling%domain)
                  case (domain_interior)
                     aggregate_variable_access%interior = ior(aggregate_variable_access%interior,model_coupling%access)
                     master_name = trim(model_coupling%model_reference%model%get_path())//'/'//trim(standard_variable%name)
                  case (domain_horizontal)
                     aggregate_variable_access%horizontal = ior(aggregate_variable_access%horizontal,model_coupling%access)
                     master_name = trim(model_coupling%model_reference%model%get_path())//'/'//trim(standard_variable%name)//'_at_interfaces'
                  case (domain_bottom)
                     aggregate_variable_access%bottom = ior(aggregate_variable_access%bottom,model_coupling%access)
                     master_name = trim(model_coupling%model_reference%model%get_path())//'/'//trim(standard_variable%name)//'_at_bottom'
                  case (domain_surface)
                     aggregate_variable_access%surface = ior(aggregate_variable_access%surface,model_coupling%access)
                     master_name = trim(model_coupling%model_reference%model%get_path())//'/'//trim(standard_variable%name)//'_at_surface'
                  end select
               end select
            end if
            call self%request_coupling(coupling%slave,master_name)
         end select
         coupling => next_coupling
      end do
   end subroutine before_coupling

   subroutine build_state_id_list(self,reference,domain)
      class (type_particle_model),intent(inout) :: self
      type (type_model_reference),intent(inout) :: reference
      integer,                    intent(in)    :: domain

      type (type_link),           pointer :: link
      integer                             :: n
      character(len=10)                   :: strindex

      ! Count number of state variables in target model.
      n = 0
      link => reference%model%links%first
      do while (associated(link))
         if (index(link%name,'/')==0.and.link%target%domain==domain.and.link%original%presence==presence_internal &
             .and.(.not.link%original%state_indices%is_empty().or.link%original%fake_state_variable)) n = n + 1
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
         if (index(link%name,'/')==0.and.link%target%domain==domain.and.link%original%presence==presence_internal &
             .and.(.not.link%original%state_indices%is_empty().or.link%original%fake_state_variable)) then
            n = n + 1
            write (strindex,'(i0)') n
            select case (domain)
               case (domain_interior)
                  call self%register_state_dependency(reference%id%state(n),trim(reference%name)//'_state'//trim(strindex), &
                     link%target%units,trim(reference%name)//' state variable '//trim(strindex))
                  call self%request_coupling_to_model(reference%id%state(n),reference%id,link%name)
               case (domain_bottom)
                  call self%register_state_dependency(reference%id%bottom_state(n),trim(reference%name)//'_bottom_state'//trim(strindex), &
                     link%target%units,trim(reference%name)//' bottom state variable '//trim(strindex))
                  call self%request_coupling_to_model(reference%id%bottom_state(n),reference%id,link%name)
               case (domain_surface)
                  call self%register_state_dependency(reference%id%surface_state(n),trim(reference%name)//'_surface_state'//trim(strindex), &
                     link%target%units,trim(reference%name)//' surface state variable '//trim(strindex))
                  call self%request_coupling_to_model(reference%id%surface_state(n),reference%id,link%name)
            end select
         end if
         link => link%next
      end do

   end subroutine build_state_id_list

end module
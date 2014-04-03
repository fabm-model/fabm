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

   type type_model_reference
      character(len=attribute_length)         :: name           = ''
      class (type_base_model),        pointer :: model          => null()
      type (type_coupling_from_model),pointer :: first_coupling => null()
      type (type_model_id),           pointer :: id             => null()
      type (type_model_reference),    pointer :: next           => null()
   end type

   type type_model_id
      type (type_model_reference),pointer :: reference => null()
      type (type_state_variable_id),allocatable :: state(:)
   end type

   type,extends(type_coupling) :: type_coupling_from_model
      type (type_model_reference),pointer :: model_reference => null()
   end type
   
   type,extends(type_base_model) :: type_particle_model
      type (type_model_reference),pointer :: first_model_reference => null()
   contains
      ! Used by inheriting biogeochemical models:
      procedure :: register_model_dependency
      procedure :: request_coupling_to_model

      ! Hooks called by FABM:
      procedure :: before_coupling
      procedure :: after_coupling
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
      type (type_model_id),target,      intent(inout) :: id
!
! !INPUT PARAMETERS:
      character(len=*),                 intent(in)    :: name
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
      reference%id => id

      id%reference => reference

   end subroutine register_model_dependency
!EOC

   subroutine request_coupling_to_model(self,slave_variable,master_model,master_variable)
      class (type_particle_model),intent(inout) :: self
      class (type_variable_id),   intent(in)    :: slave_variable
      class (type_model_id),      intent(inout) :: master_model
      character(len=*),           intent(in)    :: master_variable

      class (type_coupling_from_model), pointer :: coupling
      class (type_coupling), pointer :: plain_coupling

      ! Create object describing the coupling, and send it to FABM.
      ! This must be a pointer, because FABM will manage its memory and deallocate when appropriate.
      allocate(coupling)
      coupling%slave = slave_variable%link%name
      coupling%master = master_variable
      coupling%model_reference => master_model%reference
      plain_coupling => coupling
      call self%couplings%add(plain_coupling)
   end subroutine

   subroutine before_coupling(self)
      class (type_particle_model),intent(inout) :: self

      type (type_model_reference),pointer :: reference
      class (type_coupling),      pointer :: coupling

      reference => self%first_model_reference
      do while (associated(reference))
         ! First find a coupling for the referenced model.
         coupling => self%couplings%find(reference%name,remove=.true.)
         if (.not.associated(coupling)) call self%fatal_error('before_coupling','Model reference "'//trim(reference%name)//'" was not coupled.')

         ! Now find the actual model that was coupled.
         reference%model => self%find_model(coupling%master,recursive=.true.)
         if (.not.associated(reference%model)) call self%fatal_error('before_coupling','Referenced model "'//trim(reference%name)//'" not found.')

         ! Now deallocate coupling (self%couplings%find rendered control to us)
         deallocate(coupling)
         
         reference => reference%next
      end do

      ! Resolve all couplings to variables in specific other models.
      coupling => self%couplings%first
      do while (associated(coupling))
         select type (coupling)
            class is (type_coupling_from_model)
               coupling%master = trim(coupling%model_reference%model%name)//'_'//trim(coupling%master)
         end select
         coupling => coupling%next
      end do
   end subroutine
   
   subroutine after_coupling(self)
      class (type_particle_model),intent(inout) :: self

      type (type_model_reference),pointer :: reference
      type (type_link),           pointer :: link
      integer                             :: n

      reference => self%first_model_reference
      do while (associated(reference))
         ! Count number of state variables in target model.
         n = 0
         link => reference%model%first_link
         do while (associated(link))
            select type (object=>link%target)
               class is (type_bulk_variable)
                  if (object%presence==presence_internal.and.link%owner.and..not.object%state_indices%is_empty()) n = n + 1
            end select
            link => link%next
         end do

         ! Allocate array to hold state variable identifiers.
         allocate(reference%id%state(n))

         ! Connect target state variables to identifiers in model reference.
         n = 0
         link => reference%model%first_link
         do while (associated(link))
            select type (object=>link%target)
               class is (type_bulk_variable)
                  if (object%presence==presence_internal.and.link%owner.and..not.object%state_indices%is_empty()) then
                     n = n + 1
                     call connect_bulk_state_variable_id(self,object,reference%id%state(n))
                  end if
            end select
            link => link%next
         end do

         reference => reference%next
      end do

   end subroutine after_coupling

end module
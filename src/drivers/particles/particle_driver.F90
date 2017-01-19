module fabm_particle_driver

   use fabm
   use fabm_config
   use fabm_types, only: rk, attribute_length, output_none
   use fabm_driver, only: driver

   implicit none

   private

   public type_fabm_particle_state, type_output

   type type_output
      character(len=attribute_length) :: name      = ''
      character(len=attribute_length) :: units     = ''
      character(len=attribute_length) :: long_name = ''
      real(rk), pointer               :: data(:)   => null()
      logical                         :: save      = .false.
      type (type_output), pointer     :: specific_to => null()
      integer                         :: id        = -1
      type (type_output), pointer     :: next      => null()
   end type

   type type_fabm_particle_state
      integer                     :: npar
      integer                     :: nstate
      real(rk), allocatable       :: y(:,:)
      logical,  allocatable       :: active(:)
      logical                     :: active_hz
      type (type_model), private  :: model
      type (type_output), pointer :: first_output => null()
      real(rk), allocatable       :: particle_count(:)
      real(rk), allocatable       :: work(:,:)
   contains
      procedure :: configure
      procedure :: initialize
      procedure :: start
      procedure :: is_variable_used
      procedure :: send_data
      procedure :: send_horizontal_data
      procedure :: send_scalar
      procedure :: get_vertical_movement
      procedure :: get_sources
      procedure :: check_state
   end type

   contains

   subroutine configure(self, dictionary, npar)
!DEC$ ATTRIBUTES DLLEXPORT :: configure
      use yaml_types, only: type_node, type_dictionary, type_key_value_pair, type_scalar, type_error

      class (type_fabm_particle_state), target, intent(inout) :: self
      class (type_dictionary),                  intent(inout) :: dictionary
      integer,          optional,               intent(in)    :: npar

      integer                             :: ivar
      type (type_output),         pointer :: output
      type (type_output),         pointer :: particle_count_output
      type (type_error),          pointer :: yaml_error
      type (type_key_value_pair), pointer :: pair
      class (type_scalar),        pointer :: config_path

      ! Determine particle count (provided as argument, or otherwise read from yaml).
      if (present(npar)) then
         self%npar = npar
      else
         self%npar = dictionary%get_integer('particle_count', error=yaml_error)
         if (associated(yaml_error)) call driver%fatal_error('particle_driver::configure', trim(yaml_error%message))
      end if

      config_path => dictionary%get_scalar('fabm_configuration', required=.false., error=yaml_error)
      if (associated(yaml_error)) call driver%fatal_error('particle_driver::configure', trim(yaml_error%message))
      if (associated(config_path)) then
         ! Explicitly specified FABM configuration file
         call fabm_create_model_from_yaml_file(self%model, trim(config_path%string))
      else
         ! Default FABM configuration file (fabm.yaml)
         call fabm_create_model_from_yaml_file(self%model)
      end if
      self%nstate = size(self%model%state_variables)

      ! Allocate the particle state.
      allocate(self%y(self%npar, self%nstate))
      allocate(self%work(self%npar, self%nstate))

      particle_count_output => null()
      do ivar = 1, size(self%model%state_variables)
         output => add_output(self%model%state_variables(ivar), .true.)
         output%data => self%y(:, ivar)
      end do
      do ivar = 1, size(self%model%diagnostic_variables)
         output => add_output(self%model%diagnostic_variables(ivar), self%model%diagnostic_variables(ivar)%output /= output_none)
         output%id = ivar
      end do

      allocate(self%active(self%npar))
      self%active = .true.
      self%active_hz = .false.

   contains

      function add_output(variable, save) result(output)
         class (type_external_variable), intent(in) :: variable
         logical,                        intent(in) :: save

         type (type_output), pointer :: output

         allocate(output)
         output%name      = variable%name
         output%units     = variable%units
         output%long_name = variable%long_name
         output%save      = save
         if (variable%target%specific_to == -2) then
            ! This is a variable that describes a property specific to water parcels (e.g., age).
            ! To average this over a volume, we need to sum its value over all particles in the volume,
            ! then divide by the number of particles in the volume. The averging infrastucture
            ! does this by computing \Sigma c w / \Sigma w, for property value c and weight w, and the
            ! sums taken over all particles. For the current case, w is therefore a constant (typically 1).
            ! We explicitly create a dummy weight field with this value below.
            if (.not.associated(particle_count_output)) then
               ! Particle count variable does not exist yet - create it.
               allocate(self%particle_count(self%npar))
               self%particle_count(:) = 1
               allocate(particle_count_output)
               particle_count_output%name = 'count'
               particle_count_output%units = '# m-3'
               particle_count_output%long_name = 'count'
               particle_count_output%save = .true.
               particle_count_output%data => self%particle_count
               particle_count_output%next => self%first_output
               self%first_output => particle_count_output
            end if
            output%specific_to => particle_count_output
         elseif (variable%target%specific_to /= -1) then
            ! This is a variable that describes a property specific to one of the model's other (non-property) state variables.
            ! To average this over a volume, we need to compute \Sigma c w / \Sigma w,
            ! with c being the property value, w being the variable that the property pertains to, and sum \Sigma taken over all particles in the volume.
            stop 'unsupported'
         end if

         ! Prepend to list of outputs
         output%next => self%first_output
         self%first_output => output
      end function add_output

   end subroutine configure

   subroutine initialize(self)
!DEC$ ATTRIBUTES DLLEXPORT :: initialize
      class (type_fabm_particle_state), intent(inout) :: self

      type (type_output), pointer :: output

      ! Copy "save" attribute of particle property to FABM diagnostics
      output => self%first_output
      do while (associated(output))
         if (output%id /= -1) self%model%diagnostic_variables(output%id)%save = output%save
         output => output%next
      end do

      ! Domain size: number of particles
      call fabm_set_domain(self%model, self%npar)

      ! Provide dummy indices for top and bottom. These will not be used but FABM currently demands they are set.
      call self%model%set_bottom_index(1)
      call self%model%set_surface_index(1)

      ! Provide FABM with a pointer to the mask.
      call fabm_set_mask(self%model, self%active, self%active_hz)

      ! Provide pointers to the particle state.
      call self%model%link_all_interior_state_data(self%y)

      ! Fo saved diagnostics, obtain a pointer to their values.
      output => self%first_output
      do while (associated(output))
         if (output%id /= -1 .and. output%save) output%data => fabm_get_interior_diagnostic_data(self%model, output%id)
         output => output%next
      end do
   end subroutine initialize

   subroutine start(self, domain_per_particle)
!DEC$ ATTRIBUTES DLLEXPORT :: start
      class (type_fabm_particle_state), intent(inout) :: self
      real(rk),                         intent(in)    :: domain_per_particle

      integer :: ivar

      ! By convention, cell thicknesses associated with individual particles are zero.
      call self%send_data('cell_thickness', self%model%zero)

      call fabm_check_ready(self%model)
      call fabm_initialize_state(self%model, 1, self%npar)

      ! For non-property state variables, the initial state in FABM is given as a concentration.
      ! Here, we convert this to an amount per particle.
      do ivar=1,size(self%model%state_variables)
         if (self%model%state_variables(ivar)%target%specific_to == -1) self%y(:,ivar) = self%y(:,ivar) * domain_per_particle
      end do

      ! Make sure all diagnostics have an initial value
      call self%get_sources(self%work)
   end subroutine start

   logical function is_variable_used(self, name)
!DEC$ ATTRIBUTES DLLEXPORT :: is_variable_used
      class (type_fabm_particle_state), intent(in) :: self
      character(len=*),                 intent(in) :: name

      is_variable_used = fabm_is_variable_used(self%model%get_bulk_variable_id(name))
   end function is_variable_used

   subroutine send_data(self, name, dat)
!DEC$ ATTRIBUTES DLLEXPORT :: send_data
      class (type_fabm_particle_state), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      real(rk), target,                 intent(in)    :: dat(:)

      call self%model%link_interior_data(name, dat)
   end subroutine send_data

   subroutine send_horizontal_data(self, name, dat)
!DEC$ ATTRIBUTES DLLEXPORT :: send_horizontal_data
      class (type_fabm_particle_state), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      real(rk), target,                 intent(in)    :: dat

      call self%model%link_horizontal_data(name, dat)
   end subroutine send_horizontal_data

   subroutine send_scalar(self, name, dat)
!DEC$ ATTRIBUTES DLLEXPORT :: send_scalar_data
      class (type_fabm_particle_state), intent(inout) :: self
      character(len=*),                 intent(in)    :: name
      real(rk), target,                 intent(in)    :: dat

      call self%model%link_scalar(name, dat)
   end subroutine send_scalar

   subroutine get_vertical_movement(self, w)
!DEC$ ATTRIBUTES DLLEXPORT :: get_vertical_movement
      class (type_fabm_particle_state), intent(inout) :: self
      real(rk),                         intent(out)   :: w(1:self%npar)

      call fabm_get_vertical_movement(self%model, 1, self%npar, self%work)

      ! All variables associated with the particle should move with the same vertical velocity.
      ! Take the velocity from the first variable.
      w = self%work(:, 1)
   end subroutine get_vertical_movement

   subroutine get_sources(self, dy)
!DEC$ ATTRIBUTES DLLEXPORT :: get_sources
      class (type_fabm_particle_state), intent(inout) :: self
      real(rk),                         intent(out)   :: dy(1:self%npar, self%nstate)

      dy = 0.0_rk
      call fabm_do(self%model, 1, self%npar, dy)
   end subroutine get_sources

   function check_state(self, repair) result(valid)
!DEC$ ATTRIBUTES DLLEXPORT :: check_state
      class (type_fabm_particle_state), intent(inout) :: self
      logical,                          intent(in)    :: repair

      logical :: valid

      call fabm_check_state(self%model, 1, self%npar, repair, valid)
   end function check_state

end module

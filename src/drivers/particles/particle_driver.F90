module fabm_particle_driver

   use fabm
   use fabm_config
   use fabm_types, only: rk

   implicit none

   private

   public type_fabm_particle_state

   type type_fabm_particle_state
      integer :: npar
      integer :: nstate
      real(rk), allocatable :: y(:,:)
      logical,  allocatable :: active(:)
      logical               :: active_hz
      type (type_model) :: model
   contains
      procedure :: initialize
      procedure :: start
      procedure :: is_variable_used
      procedure :: send_data
      procedure :: get_vertical_movement
      procedure :: get_sources
   end type

   contains
   
   subroutine initialize(self, npar, config_path)
!DEC$ ATTRIBUTES DLLEXPORT :: initialize
      class (type_fabm_particle_state), intent(inout) :: self
      integer,                          intent(in)    :: npar
      character(len=*), optional,       intent(in)    :: config_path

      self%npar = npar
      if (present(config_path)) then
         call fabm_create_model_from_yaml_file(self%model, config_path)
      else
         call fabm_create_model_from_yaml_file(self%model, 'fabm.yaml')
      end if
      self%nstate = size(self%model%state_variables)

      allocate(self%active(self%npar))
      self%active = .true.
      self%active_hz = .false.
      call fabm_set_domain(self%model, self%npar)
      call self%model%set_bottom_index(1)
      call self%model%set_surface_index(1)
      call fabm_set_mask(self%model, self%active, self%active_hz)
      
      allocate(self%y(self%npar, self%nstate))
      call self%model%link_all_interior_state_data(self%y)
   end subroutine initialize

   subroutine start(self)
!DEC$ ATTRIBUTES DLLEXPORT :: start
      class (type_fabm_particle_state), intent(inout) :: self

      call fabm_check_ready(self%model)
      call fabm_initialize_state(self%model,1,self%npar)
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

   subroutine get_vertical_movement(self, w)
!DEC$ ATTRIBUTES DLLEXPORT :: get_vertical_movement
      class (type_fabm_particle_state), intent(inout) :: self
      real(rk),                         intent(out)   :: w(1:self%npar)

      real(rk) :: w_all(1:self%npar,self%nstate)

      call fabm_get_vertical_movement(self%model, 1, self%npar, w_all)
      w = w_all(:, 1)
   end subroutine get_vertical_movement

   subroutine get_sources(self, dy)
!DEC$ ATTRIBUTES DLLEXPORT :: get_sources
      class (type_fabm_particle_state), intent(inout) :: self
      real(rk),                         intent(out)   :: dy(1:self%npar, self%nstate)

      dy = 0.0_rk
      call fabm_do(self%model, 1, self%npar, dy)
   end subroutine get_sources

end module
   
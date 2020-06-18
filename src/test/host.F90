#include "fabm_driver.h"
#include "fabm_private.h"

#undef _BEGIN_OUTER_INTERIOR_LOOP_
#undef _END_OUTER_INTERIOR_LOOP_
#undef _BEGIN_OUTER_HORIZONTAL_LOOP_
#undef _END_OUTER_HORIZONTAL_LOOP_

module host_hooks
   use fabm_driver

   implicit none

   type,extends(type_base_driver) :: type_test_driver
   contains
      procedure :: fatal_error => test_driver_fatal_error
      procedure :: log_message => test_driver_log_message
   end type

contains

   subroutine test_driver_fatal_error(self,location,message)
      class (type_test_driver), intent(inout) :: self
      character(len=*),         intent(in)    :: location,message

      write (*,*) trim(location)//': '//trim(message)
      stop 1
   end subroutine

   subroutine test_driver_log_message(self,message)
      class (type_test_driver), intent(inout) :: self
      character(len=*),         intent(in)    :: message

      write (*,*) trim(message)
   end subroutine

end module host_hooks

program test_host

   use fabm
   use fabm_driver
   use fabm_parameters, only: rke
   use fabm_types, only: source_do, source_do_surface, source_do_bottom, source_do_column

   use test_models
   use host_hooks

   implicit none

#if _FABM_DIMENSION_COUNT_>0
   integer :: _LOCATION_
#endif

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
#  define _INTERIOR_SLICE_RANGE_PLUS_1_ (_START_:_STOP_,:)
#else
#  define _INTERIOR_SLICE_RANGE_PLUS_1_ (:)
#endif

#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_) &&  _FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_
#  define _HORIZONTAL_SLICE_RANGE_PLUS_1_ _INTERIOR_SLICE_RANGE_PLUS_1_
#else
#  define _HORIZONTAL_SLICE_RANGE_PLUS_1_ (:)
#endif

#if _FABM_DIMENSION_COUNT_==0
#  define _BEGIN_GLOBAL_LOOP_
#  define _END_GLOBAL_LOOP_
#elif _FABM_DIMENSION_COUNT_==1
#  define _BEGIN_GLOBAL_LOOP_ do i__=domain_start(1),domain_stop(1)
#  define _END_GLOBAL_LOOP_ end do;i__=domain_extent(1)
#  ifdef _FABM_DEPTH_DIMENSION_INDEX_
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#    define _END_GLOBAL_HORIZONTAL_LOOP_
#  endif
#elif _FABM_DIMENSION_COUNT_==2
#  define _BEGIN_GLOBAL_LOOP_ do j__=domain_start(2),domain_stop(2);do i__=domain_start(1),domain_stop(1)
#  define _END_GLOBAL_LOOP_ end do;end do;i__=domain_extent(1);j__=domain_extent(2)
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do j__=domain_start(2),domain_stop(2)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;j__=domain_extent(2)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do i__=domain_start(1),domain_stop(1)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;i__=domain_extent(1)
#  endif
#elif _FABM_DIMENSION_COUNT_==3
#  define _BEGIN_GLOBAL_LOOP_ do k__=domain_start(3),domain_stop(3);do j__=domain_start(2),domain_stop(2);do i__=domain_start(1),domain_stop(1)
#  define _END_GLOBAL_LOOP_ end do;end do;end do;i__=domain_extent(1);j__=domain_extent(2);k__=domain_extent(3)
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do k__=domain_start(3),domain_stop(3);do j__=domain_start(2),domain_stop(2)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do;j__=domain_extent(2);k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do k__=domain_start(3),domain_stop(3);do i__=domain_start(1),domain_stop(1)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do;i__=domain_extent(1);k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do j__=domain_start(2),domain_stop(2);do i__=domain_start(1),domain_stop(1)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do;i__=domain_extent(1);j__=domain_extent(2)
#  endif
#endif

! If there is no depth dimension, horizontal = global
#ifndef _FABM_DEPTH_DIMENSION_INDEX_
#  define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ _BEGIN_GLOBAL_LOOP_
#  define _END_GLOBAL_HORIZONTAL_LOOP_ _END_GLOBAL_LOOP_
#endif

#ifndef _FABM_VECTORIZED_DIMENSION_INDEX_
   ! No vectorization: outer loops are global loops
#  define _BEGIN_OUTER_INTERIOR_LOOP_ _BEGIN_GLOBAL_LOOP_
#  define _END_OUTER_INTERIOR_LOOP_ _END_GLOBAL_LOOP_
#  define _BEGIN_OUTER_HORIZONTAL_LOOP_ _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#  define _END_OUTER_HORIZONTAL_LOOP_ _END_GLOBAL_HORIZONTAL_LOOP_
#elif _FABM_DIMENSION_COUNT_==1
   ! Entire domain is vectorized; no outer loops needed
#  define _BEGIN_OUTER_INTERIOR_LOOP_
#  define _END_OUTER_INTERIOR_LOOP_
#  define _BEGIN_OUTER_HORIZONTAL_LOOP_
#  define _END_OUTER_HORIZONTAL_LOOP_
#elif _FABM_DIMENSION_COUNT_==2
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do j__=domain_start(2),domain_stop(2)
#  define _END_OUTER_INTERIOR_LOOP_ end do;j__=domain_extent(2)
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
     ! The entire horizontal is already vectorized; no outer loop necessary
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_
#    define _END_OUTER_HORIZONTAL_LOOP_
#  else
     ! No horizontal dimension vectorized; do full outer loop.
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#    define _END_OUTER_HORIZONTAL_LOOP_ _END_GLOBAL_HORIZONTAL_LOOP_
#  endif
#elif _FABM_DIMENSION_COUNT_==3
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do k__=domain_start(3),domain_stop(3);do j__=domain_start(2),domain_stop(2)
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do;j__=domain_extent(2);k__=domain_extent(3)
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do k__=domain_start(3),domain_stop(3)
#    define _END_OUTER_HORIZONTAL_LOOP_ end do;k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do j__=domain_start(2),domain_stop(2)
#    define _END_OUTER_HORIZONTAL_LOOP_ end do;j__=domain_extent(2)
#  else
     ! No horizontal dimension vectorized; do full outer loop.
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#    define _END_OUTER_HORIZONTAL_LOOP_ _END_GLOBAL_HORIZONTAL_LOOP_
#  endif
#endif

#ifdef _INTERIOR_IS_VECTORIZED_
   integer :: _START_, _STOP_
#endif

   real(rke), allocatable _DIMENSION_GLOBAL_            :: tmp
   real(rke), allocatable _DIMENSION_GLOBAL_HORIZONTAL_ :: tmp_hz

#ifdef _HAS_MASK_
   _FABM_MASK_TYPE_, allocatable, target _DIMENSION_GLOBAL_HORIZONTAL_ :: mask_hz
#  ifndef _FABM_HORIZONTAL_MASK_
   _FABM_MASK_TYPE_, allocatable, target _DIMENSION_GLOBAL_ :: mask
#  endif
#  ifndef _FABM_MASKED_VALUE_
#    define _FABM_MASKED_VALUE_ _FABM_UNMASKED_VALUE_+1
#  endif
#  ifndef _FABM_UNMASKED_VALUE_
#    define _FABM_UNMASKED_VALUE_ _FABM_MASKED_VALUE_+1
#  endif
#endif

   integer :: interior_count
   integer :: horizontal_count

#if _FABM_BOTTOM_INDEX_==-1
   integer, allocatable, target _DIMENSION_GLOBAL_HORIZONTAL_ :: bottom_index
#endif

   real(rke), allocatable, target _DIMENSION_GLOBAL_PLUS_1_            :: interior_state
   real(rke), allocatable, target _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: surface_state
   real(rke), allocatable, target _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: bottom_state

   real(rke), allocatable _DIMENSION_SLICE_PLUS_1_            :: dy, w, total_int
   real(rke), allocatable _DIMENSION_HORIZONTAL_SLICE_PLUS_1_ :: flux, sms_sf, sms_bt, total_hz

   real(rke), allocatable, target _DIMENSION_GLOBAL_            :: temperature
   real(rke), allocatable, target _DIMENSION_GLOBAL_            :: depth
   real(rke), allocatable, target _DIMENSION_GLOBAL_HORIZONTAL_ :: wind_speed

   class (type_fabm_model), pointer :: model

   class (type_test_model), pointer :: test_model

   integer :: domain_extent(_FABM_DIMENSION_COUNT_)
   integer :: domain_start(_FABM_DIMENSION_COUNT_)
   integer :: domain_stop(_FABM_DIMENSION_COUNT_)

   character(len=20) :: arg
   integer :: ivar
   integer :: i
   integer :: mode = 1
   integer :: ntest = -1
   logical :: no_mask = .false.

#if _FABM_DIMENSION_COUNT_>0
   i__ = 50
#endif
#if _FABM_DIMENSION_COUNT_>1
   j__ = 40
#endif
#if _FABM_DIMENSION_COUNT_>2
   k__ = 45
#endif

   ! Parse command line arguments
   call start_test('parsing command line arguments')
   i = 1
   do
      call get_command_argument(i, arg)
      if (arg == '') exit
      select case (arg)
      case ('-s', '--simulate')
         mode = 2
      case ('--nomask')
         no_mask = .true.
#if _FABM_DIMENSION_COUNT_>0
      case ('--nx')
         i = i + 1
         call get_command_argument(i, arg)
         read (arg,*) i__
#endif
#if _FABM_DIMENSION_COUNT_>1
      case ('--ny')
         i = i + 1
         call get_command_argument(i, arg)
         read (arg,*) j__
#endif
#if _FABM_DIMENSION_COUNT_>2
      case ('--nz')
         i = i + 1
         call get_command_argument(i, arg)
         read (arg,*) k__
#endif
      case ('-n')
         i = i + 1
         call get_command_argument(i, arg)
         read (arg,*) ntest
      case ('-h')
         write (*,'(a)') ''
         write (*,'(a)') ''
         write (*,'(a)') 'FABM host emulator'
         write (*,'(a)') ''
         write (*,'(a)') 'Accepted arguments:'
         write (*,'(a)') '-s/--simulate: simulate using provided fabm.yaml/environment.yaml'
         write (*,'(a)') '-n:            number of replicates when simulating'
         stop 0
      case default
         write (*,'(a)') 'Unknown command line argument: ' // trim(arg)
         stop 2
      end select
      i = i + 1
   end do
   call report_test_result()

#if _FABM_DIMENSION_COUNT_>0
   domain_extent = (/ _LOCATION_ /)
   domain_start(:) = 1
   domain_stop(:) = domain_extent
   interior_count = product(domain_stop - domain_start + 1)
#else
   interior_count = 1
#endif

   ! Set defaults
   if (ntest == -1) then
      if (mode == 1) then
         ntest = 1
      else
         ntest = 50000000 / interior_count
      end if
   end if

#ifdef _INTERIOR_IS_VECTORIZED_
   _START_ = domain_start(_FABM_VECTORIZED_DIMENSION_INDEX_)
   _STOP_ = domain_stop(_FABM_VECTORIZED_DIMENSION_INDEX_)
#endif

   allocate(tmp _INDEX_LOCATION_)
   allocate(tmp_hz _INDEX_HORIZONTAL_LOCATION_)

   allocate(type_test_driver::driver)

   call start_test('fabm_initialize_library')
   call fabm_initialize_library()
   call report_test_result()

   call start_test('building model tree')
   select case (mode)
   case (1)
       ! Unit testing with built-in model
       allocate(model)
       allocate(test_model)
       call model%root%add_child(test_model, 'test_model', 'test model', configunit=-1)
   case (2)
       ! Test with user-provided fabm.yaml
       model => fabm_create_model(initialize=.false.)
   end select
   call report_test_result()

   call start_test('initialize')
   call model%initialize()
   if (mode == 1) then
      call assert(size(model%interior_state_variables) == test_model%nstate, 'model%initialize', 'Incorrect number of interior state variables.')
      call assert(size(model%bottom_state_variables) == test_model%nbottom_state, 'model%initialize', 'Incorrect number of bottom state variables.')
      call assert(size(model%surface_state_variables) == test_model%nsurface_state, 'model%initialize', 'Incorrect number of surface state variables.')
      call assert(size(model%conserved_quantities) == 1, 'model%initialize', 'Incorrect number of conserved quantities.')
   end if
   call report_test_result()

   ! ======================================================================
   ! Provide extents of the spatial domain.
   ! ======================================================================

   call start_test('set_domain')
   call model%set_domain(_LOCATION_)
   call report_test_result()

   ! ======================================================================
   ! Set up spatial mask.
   ! ======================================================================

#ifdef _HAS_MASK_
   allocate(mask_hz _INDEX_HORIZONTAL_LOCATION_)
   call start_test('set_mask')
#  ifdef _FABM_HORIZONTAL_MASK_
   call model%set_mask(mask_hz)
#  else
   allocate(mask _INDEX_LOCATION_)
   call model%set_mask(mask, mask_hz)
#  endif
   call report_test_result()
#endif

   ! ======================================================================
   ! Specify vertical indices of surface and bottom.
   ! ======================================================================

#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_BOTTOM_INDEX_==-1
   call start_test('set_bottom_index')
   allocate(bottom_index _INDEX_HORIZONTAL_LOCATION_)
   call model%set_bottom_index(bottom_index)
   call report_test_result()
#endif

   allocate(interior_state(_PREARG_LOCATION_ size(model%interior_state_variables)))
   allocate(surface_state(_PREARG_HORIZONTAL_LOCATION_ size(model%surface_state_variables)))
   allocate(bottom_state(_PREARG_HORIZONTAL_LOCATION_ size(model%bottom_state_variables)))

   ! ======================================================================
   ! Send pointers to state variable data to FABM.
   ! ======================================================================

   call start_test('link_interior_state_data')
   do ivar = 1, size(model%interior_state_variables)
      call model%link_interior_state_data(ivar, interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar))
   end do
   call report_test_result()

   call start_test('link_surface_state_data')
   do ivar = 1, size(model%surface_state_variables)
      call model%link_surface_state_data(ivar, surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar))
   end do
   call report_test_result()

   call start_test('link_bottom_state_data')
   do ivar = 1, size(model%bottom_state_variables)
      call model%link_bottom_state_data(ivar, bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar))
   end do
   call report_test_result()

   ! ======================================================================
   ! Transfer pointers to environmental data
   ! ======================================================================

   select case (mode)
   case (1)
       allocate(depth _INDEX_LOCATION_)
       allocate(temperature _INDEX_LOCATION_)
       allocate(wind_speed _INDEX_HORIZONTAL_LOCATION_)

       call start_test('link_interior_data')
       call model%link_interior_data(fabm_standard_variables%temperature, temperature)
       call model%link_interior_data(fabm_standard_variables%depth, depth)
       call report_test_result()

       call start_test('link_horizontal_data')
       call model%link_horizontal_data(fabm_standard_variables%wind_speed, wind_speed)
       call report_test_result()
   case (2)
       call read_environment
   end select

   ! ======================================================================
   ! Check whether FABM has all dependencies fulfilled
   ! (i.e., whether all required calls for link_*_data have been made)
   ! ======================================================================

   call start_test('start')
   call model%start()
   call report_test_result()

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
   allocate(dy(_START_:_STOP_, size(model%interior_state_variables)))
   allocate(w(_START_:_STOP_, size(model%interior_state_variables)))
   allocate(total_int(_START_:_STOP_, size(model%conserved_quantities)))
#else
   allocate(dy(size(model%interior_state_variables)))
   allocate(w(size(model%interior_state_variables)))
   allocate(total_int(size(model%conserved_quantities)))
#endif

#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&(_FABM_DEPTH_DIMENSION_INDEX_!=_FABM_VECTORIZED_DIMENSION_INDEX_)
   allocate(flux(_START_:_STOP_, size(model%interior_state_variables)))
   allocate(sms_sf(_START_:_STOP_, size(model%surface_state_variables)))
   allocate(sms_bt(_START_:_STOP_, size(model%bottom_state_variables)))
   allocate(total_hz(_START_:_STOP_, size(model%conserved_quantities)))
#else
   allocate(flux(size(model%interior_state_variables)))
   allocate(sms_sf(size(model%surface_state_variables)))
   allocate(sms_bt(size(model%bottom_state_variables)))
   allocate(total_hz(size(model%conserved_quantities)))
#endif

   select case (mode)
   case (1)
      write (*,'(a)') 'Testing without mask and unrestricted domain.'
      call test_update(apply_mask=.false., restrict_range=.false.)
#if defined(_HAS_MASK_) || _FABM_DIMENSION_COUNT_>0
      write (*,'(a,i0,a)') 'Running ', ntest, ' tests with randomized domain settings.'
      do i=1,ntest
#  ifdef _HAS_MASK_
         write (*,'(a)') 'Random mask (unrestricted domain).'
         call test_update(apply_mask=.true., restrict_range=.false.)
#  endif
#  if _FABM_DIMENSION_COUNT_>0
         write (*,'(a)') 'Randomly restricted domain (no mask).'
         call test_update(apply_mask=.false., restrict_range=.true.)
#  endif
#  if defined(_HAS_MASK_) && _FABM_DIMENSION_COUNT_>0
         write (*,'(a)') 'Random mask and randomly restricted domain.'
         call test_update(apply_mask=.true., restrict_range=.true.)
#  endif
      end do
#endif
   case(2)
      call simulate(ntest)
   end select

contains

   subroutine read_environment
      use yaml, only: yaml_parse => parse, yaml_error_length => error_length
      use yaml_types, only: type_node, type_yaml_dictionary => type_dictionary, type_yaml_scalar => type_scalar, &
         type_yaml_key_value_pair => type_key_value_pair, yaml_real_kind => real_kind

      integer, parameter :: yaml_unit = 100
      character(yaml_error_length) :: yaml_error
      class (type_node),pointer :: yaml_root
      type (type_yaml_key_value_pair), pointer :: yaml_pair
      real(rke) :: value
      logical :: success
      type type_input
         type (type_fabm_interior_variable_id)                :: interior_id
         type (type_fabm_horizontal_variable_id)              :: horizontal_id
         type (type_fabm_scalar_variable_id)                  :: scalar_id
         real(rke), allocatable _DIMENSION_GLOBAL_            :: interior_data
         real(rke), allocatable _DIMENSION_GLOBAL_HORIZONTAL_ :: horizontal_data
         real(rke)                                            :: scalar_data
      end type
      type (type_input), pointer :: input

      yaml_root => yaml_parse('environment.yaml', yaml_unit, yaml_error)
      if (yaml_error /= '') then
         call driver%log_message(yaml_error)
         stop 2
      end if
      select type (yaml_root)
      class is (type_yaml_dictionary)
          yaml_pair => yaml_root%first
          do while (associated(yaml_pair))
              select type (node => yaml_pair%value)
              class is (type_yaml_scalar)
                  call driver%log_message('Setting '//trim(yaml_pair%key)//' to '//trim(node%string))
                  value = node%to_real(0._yaml_real_kind, success)
                  if (.not. success) then
                     call driver%log_message('Cannot parse '//trim(node%string)//' as real.')
                     stop 2
                  end if
                  allocate(input)
                  input%interior_id = model%get_interior_variable_id(trim(yaml_pair%key))
                  if (model%is_variable_used(input%interior_id)) then
                      allocate(input%interior_data _INDEX_LOCATION_)
                      input%interior_data = value
                      call model%link_interior_data(input%interior_id, input%interior_data)
                  else
                      input%horizontal_id = model%get_horizontal_variable_id(trim(yaml_pair%key))
                      if (model%is_variable_used(input%horizontal_id)) then
                         allocate(input%horizontal_data _INDEX_HORIZONTAL_LOCATION_)
                         input%horizontal_data = value
                         call model%link_horizontal_data(input%horizontal_id, input%horizontal_data)
                      else
                         input%scalar_id = model%get_scalar_variable_id(trim(yaml_pair%key))
                         if (model%is_variable_used(input%scalar_id)) then
                            input%scalar_data = value
                            call model%link_scalar(input%scalar_id, input%scalar_data)
                         else
                            call driver%log_message('WARNING: environment variable '//trim(yaml_pair%key) &
                               //' is not used by FABM model and will be ignored.')
                         end if
                      end if
                  end if
              end select
              yaml_pair => yaml_pair%next
          end do
      class default
         call driver%log_message('environment.yaml should contain a dictionary at root level')
         stop 2
      end select
   end subroutine read_environment

   subroutine configure_range(randomize)
      logical, intent(in) :: randomize

      real(rke), parameter :: excluded_fraction = 0.5_rke

#if _FABM_DIMENSION_COUNT_ > 0
      integer :: i
      real(rke) :: rnd(2)
      real(rke), parameter :: cut_fraction = 0.5_rke * (1._rke - excluded_fraction ** (1._rke / _FABM_DIMENSION_COUNT_))

      if (.not. randomize) then
         domain_start(:) = 1
         domain_stop(:) = domain_extent
      else
         do i = 1, _FABM_DIMENSION_COUNT_
            call random_number(rnd)
            domain_start(i) = 1 + int(domain_extent(i) * cut_fraction * rnd(1))
            domain_stop(i) = domain_start(i) + int((domain_extent(i) - domain_start(i) + 1) * (1._rke - cut_fraction) * rnd(2))
            write (*,'(A,I0,A,I0,A,I0)') 'Dimension ', i, ': ', domain_start(i), ' - ', domain_stop(i)
         end do
      end if

#  if _FABM_DIMENSION_COUNT_ == 1
      call model%set_domain_start(domain_start(1))
      call model%set_domain_stop(domain_stop(1))
#  elif _FABM_DIMENSION_COUNT_ == 2
      call model%set_domain_start(domain_start(1), domain_start(2))
      call model%set_domain_stop(domain_stop(1), domain_stop(2))
#  elif _FABM_DIMENSION_COUNT_ == 3
      call model%set_domain_start(domain_start(1), domain_start(2), domain_start(3))
      call model%set_domain_stop(domain_stop(1), domain_stop(2), domain_stop(3))
#  endif
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      _START_ = domain_start(_FABM_VECTORIZED_DIMENSION_INDEX_)
      _STOP_ = domain_stop(_FABM_VECTORIZED_DIMENSION_INDEX_)
#  endif
#endif
   end subroutine configure_range

   subroutine configure_mask(randomize)
      logical, intent(in) :: randomize

      real(rke), parameter :: masked_fraction = 0.5_rke

      if (randomize) then
#if _FABM_BOTTOM_INDEX_==-1
      ! Depth index of bottom varies in the horizontal
      call random_number(tmp_hz)
#  ifdef _HAS_MASK_
      !  Pick random numbers between start-1 and stop index [inclusive]. start-1 means invalid (land)
      bottom_index = domain_start(_FABM_DEPTH_DIMENSION_INDEX_) - 1 + floor(tmp_hz * (2 + domain_stop(_FABM_DEPTH_DIMENSION_INDEX_) - domain_start(_FABM_DEPTH_DIMENSION_INDEX_)))
#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
      ! Use stop+1 as invalid bottom index to ensure vertical loops will have 0 iterations
      bottom_index = bottom_index + 1
#    endif
#  else
      ! Pick random numbers between start and stop index [inclusive]
      bottom_index = domain_start(_FABM_DEPTH_DIMENSION_INDEX_) + floor(tmp_hz * (1 + domain_stop(_FABM_DEPTH_DIMENSION_INDEX_) - domain_start(_FABM_DEPTH_DIMENSION_INDEX_)))
#  endif
#endif

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      ! Apply random mask across horizontal domain (half of grid cells masked)
      call random_number(tmp_hz)
      mask_hz = _FABM_UNMASKED_VALUE_
      where (tmp_hz > masked_fraction) mask_hz = _FABM_MASKED_VALUE_
#  else
      ! Apply random mask across interior domain (half of grid cells masked)
      call random_number(tmp)
      mask = _FABM_UNMASKED_VALUE_
      where (tmp > masked_fraction) mask = _FABM_MASKED_VALUE_

#    if _FABM_BOTTOM_INDEX_==-1
      ! Bottom index varies in the horizontal. Ensure the bottom cell itself is unmasked, and anything deeper is masked.
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         ! Valid bottom index - unmask associated cell, then mask all deeper ones
#      ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) mask _INDEX_GLOBAL_VERTICAL_(bottom_index _INDEX_HORIZONTAL_LOCATION_) = _FABM_UNMASKED_VALUE_
         mask _INDEX_GLOBAL_VERTICAL_(:bottom_index _INDEX_HORIZONTAL_LOCATION_ - 1) = _FABM_MASKED_VALUE_
#      else
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= domain_start(_FABM_DEPTH_DIMENSION_INDEX_)) mask _INDEX_GLOBAL_VERTICAL_(bottom_index _INDEX_HORIZONTAL_LOCATION_) = _FABM_UNMASKED_VALUE_
         mask _INDEX_GLOBAL_VERTICAL_(bottom_index _INDEX_HORIZONTAL_LOCATION_ + 1:) = _FABM_MASKED_VALUE_
#      endif
      _END_GLOBAL_HORIZONTAL_LOOP_
#    endif

      ! Infer horizontal mask (mask points that have all column layers masked)
      mask_hz = _FABM_MASKED_VALUE_
      _BEGIN_GLOBAL_LOOP_
         if (_IS_UNMASKED_(mask _INDEX_LOCATION_)) mask_hz _INDEX_HORIZONTAL_LOCATION_ = _FABM_UNMASKED_VALUE_
      _END_GLOBAL_LOOP_

      ! For valid points in the horizontal, make sure that index 1 (bottom or surface) is unmasked
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         if (_IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)) then
            mask _INDEX_GLOBAL_VERTICAL_(domain_start(_FABM_DEPTH_DIMENSION_INDEX_)) = _FABM_UNMASKED_VALUE_
#    if _FABM_BOTTOM_INDEX_!=-1
            mask _INDEX_GLOBAL_VERTICAL_(domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) = _FABM_UNMASKED_VALUE_
#    endif
         end if
      _END_GLOBAL_HORIZONTAL_LOOP_
#  endif
#endif
      else
#if _FABM_BOTTOM_INDEX_==-1
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
      bottom_index = domain_start(_FABM_DEPTH_DIMENSION_INDEX_)
#  else
      bottom_index = domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)
#  endif
#endif
#ifdef _HAS_MASK_
#  ifndef _FABM_HORIZONTAL_MASK_
      mask = _FABM_UNMASKED_VALUE_
#  endif
      mask_hz = _FABM_UNMASKED_VALUE_
#endif
      end if

#if _FABM_BOTTOM_INDEX_==-1 && defined(_HAS_MASK_)
      ! Bottom index varies in the horizontal. Ensure the bottom cell itself is unmasked, and anything deeper is masked.
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         if (_IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)) then
            call assert(bottom_index _INDEX_HORIZONTAL_LOCATION_ >= domain_start(_FABM_DEPTH_DIMENSION_INDEX_) &
              .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_stop(_FABM_DEPTH_DIMENSION_INDEX_), 'randomize_mask', &
               'BUG: unmaked horizontal location has invalid bottom index')
         end if
      _END_GLOBAL_HORIZONTAL_LOOP_
#endif

      call count_active_points()
   end subroutine configure_mask

   subroutine count_active_points()
      logical :: active
      integer :: i, nhz, nhz_active

      active = .true.

      interior_count = 0
      _BEGIN_GLOBAL_LOOP_
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
         active = _IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)
#  else
         active = _IS_UNMASKED_(mask _INDEX_LOCATION_)
#  endif
#elif _FABM_BOTTOM_INDEX_==-1 && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         active = _ITERATOR_ >= bottom_index _INDEX_HORIZONTAL_LOCATION_
#  else
         active = _ITERATOR_ <= bottom_index _INDEX_HORIZONTAL_LOCATION_
#  endif
#endif
         if (active) interior_count = interior_count + 1
      _END_GLOBAL_LOOP_

      horizontal_count = 0
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#ifdef _HAS_MASK_
         active = _IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)
#endif
         if (active) horizontal_count = horizontal_count + 1
      _END_GLOBAL_HORIZONTAL_LOOP_

      nhz = 1
      nhz_active = 1
      do i = 1, _FABM_DIMENSION_COUNT_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
         if (i /= _FABM_DEPTH_DIMENSION_INDEX_) then
#endif
         nhz = nhz * domain_extent(i)
         nhz_active = nhz_active * (domain_stop(i) - domain_start(i) + 1)
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
         end if
#endif
      end do
      write (*,'(a,i0,a,i0,a,i0,a)') 'Interior: ', product(domain_extent), ' points, ', product(domain_stop - domain_start + 1), ' in active range, ', interior_count, ' unmasked'
      write (*,'(a,i0,a,i0,a,i0,a)') 'Horizontal: ', nhz, ' points, ', nhz_active, ' in active range, ', horizontal_count, ' unmasked'
   end subroutine

   subroutine simulate(n)
      integer, intent(in) :: n
      real(rke) :: time_begin, time_end
      integer :: nseed
      integer, allocatable :: seed(:)

      call random_seed(size=nseed)
      allocate(seed(nseed))
      seed(:) = 1
      call random_seed(put=seed)

      call configure_range(.false.)
      call configure_mask(.not. no_mask)

      call start_test('initialize_interior_state')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%initialize_interior_state(_ARG_INTERIOR_IN_)
      _END_OUTER_INTERIOR_LOOP_
      call report_test_result()

      call start_test('initialize_bottom_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%initialize_bottom_state(_ARG_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
      call report_test_result()

      call start_test('initialize_surface_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%initialize_surface_state(_ARG_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
      call report_test_result()

      write (*,'(a,i0,a)') 'Simulating with ', interior_count, ' wet cells...'

      call cpu_time(time_begin)

      do i = 1, n
         call model%prepare_inputs()

         _BEGIN_OUTER_HORIZONTAL_LOOP_
            flux = 0
            sms_bt = 0
            call model%get_bottom_sources(_PREARG_HORIZONTAL_IN_ flux, sms_bt)
         _END_OUTER_HORIZONTAL_LOOP_

         _BEGIN_OUTER_HORIZONTAL_LOOP_
            flux = 0
            sms_sf = 0
            call model%get_surface_sources(_PREARG_HORIZONTAL_IN_ flux, sms_sf)
         _END_OUTER_HORIZONTAL_LOOP_

         _BEGIN_OUTER_INTERIOR_LOOP_
            dy = 0
            call model%get_interior_sources(_PREARG_INTERIOR_IN_ dy)
         _END_OUTER_INTERIOR_LOOP_

         call model%finalize_outputs()

         if (mod(i, 100) == 0) write (*,'(i0,a)') int(100*i/real(n, rke)), ' % complete'
      end do

      call cpu_time(time_end)

      write (*,'(a)') 'Simulation complete.'
      write (*,'(a,f8.3,a)') 'Total time: ', time_end - time_begin, ' s'
   end subroutine

   subroutine test_update(apply_mask, restrict_range)
      logical, intent(in) :: apply_mask, restrict_range

      real(rke), pointer _DIMENSION_GLOBAL_            :: pdata
      real(rke), pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: pdata_hz
      logical                                          :: valid

      call configure_range(randomize=restrict_range)
      call configure_mask(randomize=apply_mask)

      call model%start()

      ! ======================================================================
      ! Initialize all state variables
      ! ======================================================================

      call start_test('initialize_interior_state')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%initialize_interior_state(_ARG_INTERIOR_IN_)
      _END_OUTER_INTERIOR_LOOP_
      do ivar = 1, size(model%interior_state_variables)
         call check_interior(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), &
            model%interior_state_variables(ivar)%missing_value, ivar+interior_state_offset+1._rke)
      end do
      call report_test_result()

      call start_test('initialize_bottom_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%initialize_bottom_state(_ARG_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar = 1, size(model%bottom_state_variables)
         call check_horizontal(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), &
            model%bottom_state_variables(ivar)%missing_value, ivar+bottom_state_offset+1._rke)
      end do
      call report_test_result()

      call start_test('initialize_surface_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%initialize_surface_state(_ARG_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar = 1, size(model%surface_state_variables)
         call check_horizontal(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), &
            model%surface_state_variables(ivar)%missing_value, ivar+surface_state_offset+1._rke)
      end do
      call report_test_result()

      ! ======================================================================
      ! Initialize environmental dependencies
      ! ======================================================================

      temperature = 1 + interior_dependency_offset
      call apply_mask_3d(temperature, -999._rke - interior_dependency_offset)

      wind_speed = 1 + horizontal_dependency_offset
      call apply_mask_2d(wind_speed, -999._rke - horizontal_dependency_offset)

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      ! Model has depth dimension: make sure depth varies from 0 at the surface till 1 at the bottom
      _BEGIN_GLOBAL_LOOP_
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
#    if _FABM_BOTTOM_INDEX_==-1
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) then
            depth _INDEX_LOCATION_ = 2
         else
            depth _INDEX_LOCATION_ = real(domain_stop(_FABM_DEPTH_DIMENSION_INDEX_) - _VERTICAL_ITERATOR_, rke) / (domain_stop(_FABM_DEPTH_DIMENSION_INDEX_) - bottom_index _INDEX_HORIZONTAL_LOCATION_)
         end if
#    else
         depth _INDEX_LOCATION_ = real(domain_stop(_FABM_DEPTH_DIMENSION_INDEX_) - _VERTICAL_ITERATOR_, rke) / (domain_stop(_FABM_DEPTH_DIMENSION_INDEX_) - domain_start(_FABM_DEPTH_DIMENSION_INDEX_))
#    endif
#  else
#    if _FABM_BOTTOM_INDEX_==-1
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_start(_FABM_DEPTH_DIMENSION_INDEX_)) then
            depth _INDEX_LOCATION_ = 2
         else
            depth _INDEX_LOCATION_ = real(_VERTICAL_ITERATOR_ - domain_start(_FABM_DEPTH_DIMENSION_INDEX_), rke) / (bottom_index _INDEX_HORIZONTAL_LOCATION_ - domain_start(_FABM_DEPTH_DIMENSION_INDEX_))
         end if
#    else
         depth _INDEX_LOCATION_ = real(_VERTICAL_ITERATOR_ - domain_start(_FABM_DEPTH_DIMENSION_INDEX_), rke) / (domain_stop(_FABM_DEPTH_DIMENSION_INDEX_) - domain_start(_FABM_DEPTH_DIMENSION_INDEX_))
#    endif
#  endif
      _END_GLOBAL_LOOP_
#else
      ! No depth dimension
      depth = 2
#endif
      call apply_mask_3d(depth, -999._rke - interior_dependency_offset)

      do ivar = 1, size(model%interior_state_variables)
         interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar) = ivar + interior_state_offset
         call apply_mask_3d(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), model%interior_state_variables(ivar)%missing_value)
      end do
      do ivar = 1, size(model%surface_state_variables)
         surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = ivar + surface_state_offset
         call apply_mask_2d(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%surface_state_variables(ivar)%missing_value)
      end do
      do ivar = 1, size(model%bottom_state_variables)
         bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = ivar + bottom_state_offset
         call apply_mask_2d(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%bottom_state_variables(ivar)%missing_value)
      end do

      ! ======================================================================
      ! Preprocessing
      ! ======================================================================

      column_loop_count = 0

      call start_test('prepare_inputs')
      call model%prepare_inputs()
      call report_test_result()

      ! ======================================================================
      ! Retrieve source terms of interior state variables.
      ! ======================================================================

      call start_test('get_interior_sources')
      interior_loop_count = 0
      _BEGIN_OUTER_INTERIOR_LOOP_
         dy = 0
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_) && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
         ! We are looping over depth, but as we have a non-constant bottom index (yet no mask), we need to skip everything below bottom
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         _START_ = bottom_index _INDEX_HORIZONTAL_LOCATION_
#  else
         _STOP_ = bottom_index _INDEX_HORIZONTAL_LOCATION_
#  endif
#endif
         call model%get_interior_sources(_PREARG_INTERIOR_IN_ dy _INTERIOR_SLICE_RANGE_PLUS_1_)
         do ivar = 1, size(model%interior_state_variables)
            call check_interior_slice_plus_1(dy _INTERIOR_SLICE_RANGE_PLUS_1_, ivar, 0.0_rke, -real(ivar + interior_state_offset, rke) _POSTARG_INTERIOR_IN_)
         end do
      _END_OUTER_INTERIOR_LOOP_

#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_) && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
      _START_ = domain_start(_FABM_VECTORIZED_DIMENSION_INDEX_)
      _STOP_ = domain_stop(_FABM_VECTORIZED_DIMENSION_INDEX_)
#  endif

      call assert(interior_loop_count == interior_count, 'get_interior_sources', &
         'call count does not match number of (unmasked) interior points')

      do ivar = 1, size(model%interior_diagnostic_variables)
         if (model%interior_diagnostic_variables(ivar)%save .and. model%interior_diagnostic_variables(ivar)%target%source == source_do .and. associated(model%interior_diagnostic_variables(ivar)%target%owner,test_model)) then
            pdata => model%get_interior_diagnostic_data(ivar)
            call check_interior(pdata, model%interior_diagnostic_variables(ivar)%missing_value, &
               -model%interior_diagnostic_variables(ivar)%missing_value)
         end if
      end do

      call report_test_result()

      ! ======================================================================
      ! Retrieve surface fluxes of interior state variables, source terms of surface-associated state variables.
      ! ======================================================================

      ! Make sure depth at surface is exactly 0
#if _FABM_BOTTOM_INDEX_==-1
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) depth _INDEX_GLOBAL_VERTICAL_(domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) = 0
#  else
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_start(_FABM_DEPTH_DIMENSION_INDEX_)) depth _INDEX_GLOBAL_VERTICAL_(domain_start(_FABM_DEPTH_DIMENSION_INDEX_)) = 0
#  endif
      _END_GLOBAL_HORIZONTAL_LOOP_
#endif

      call start_test('get_surface_sources')
      surface_loop_count = 0
      _BEGIN_OUTER_HORIZONTAL_LOOP_
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= domain_start(_FABM_DEPTH_DIMENSION_INDEX_) .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) then
#endif
         flux = 0
         sms_sf = 0
         call model%get_surface_sources(_PREARG_HORIZONTAL_IN_ flux _HORIZONTAL_SLICE_RANGE_PLUS_1_, sms_sf _HORIZONTAL_SLICE_RANGE_PLUS_1_)
         do ivar = 1, size(model%interior_state_variables)
            call check_horizontal_slice_plus_1(flux _HORIZONTAL_SLICE_RANGE_PLUS_1_, ivar, 0.0_rke, -real(ivar + interior_state_offset, rke) _POSTARG_HORIZONTAL_IN_)
         end do
         do ivar = 1, size(model%surface_state_variables)
            call check_horizontal_slice_plus_1(sms_sf _HORIZONTAL_SLICE_RANGE_PLUS_1_, ivar, 0.0_rke, -real(ivar + surface_state_offset, rke) _POSTARG_HORIZONTAL_IN_)
         end do
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         end if
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      call assert(surface_loop_count == horizontal_count, 'get_surface_sources', &
         'call count does not match number of (unmasked) horizontal points')

      do ivar = 1, size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(ivar)%save .and. &
               model%horizontal_diagnostic_variables(ivar)%target%source == source_do_surface) then
            pdata_hz => model%get_horizontal_diagnostic_data(ivar)
            call check_horizontal(pdata_hz, model%horizontal_diagnostic_variables(ivar)%missing_value, &
               -model%horizontal_diagnostic_variables(ivar)%missing_value)
         end if
      end do

      call report_test_result()

      ! ======================================================================
      ! Retrieve bottom fluxes of interior state variables, source terms of bottom-associated state variables.
      ! ======================================================================

      ! Make sure depth at bottom is exactly 1
#if _FABM_BOTTOM_INDEX_==-1
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) depth _INDEX_GLOBAL_VERTICAL_(domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) = 1
#  else
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_start(_FABM_DEPTH_DIMENSION_INDEX_)) depth _INDEX_GLOBAL_VERTICAL_(domain_start(_FABM_DEPTH_DIMENSION_INDEX_)) = 1
#  endif
      _END_GLOBAL_HORIZONTAL_LOOP_
#endif

      call start_test('get_bottom_sources')
      bottom_loop_count = 0
      _BEGIN_OUTER_HORIZONTAL_LOOP_
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= domain_start(_FABM_DEPTH_DIMENSION_INDEX_) .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) then
#endif
         flux = 0
         sms_bt = 0
         call model%get_bottom_sources(_PREARG_HORIZONTAL_IN_ flux _HORIZONTAL_SLICE_RANGE_PLUS_1_, sms_bt _HORIZONTAL_SLICE_RANGE_PLUS_1_)
         do ivar = 1, size(model%interior_state_variables)
            call check_horizontal_slice_plus_1(flux _HORIZONTAL_SLICE_RANGE_PLUS_1_, ivar, 0.0_rke, -real(ivar + interior_state_offset, rke) _POSTARG_HORIZONTAL_IN_)
         end do
         do ivar = 1, size(model%bottom_state_variables)
            call check_horizontal_slice_plus_1(sms_bt _HORIZONTAL_SLICE_RANGE_PLUS_1_, ivar, 0.0_rke, -real(ivar + bottom_state_offset, rke) _POSTARG_HORIZONTAL_IN_)
         end do
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         endif
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      call assert(bottom_loop_count == horizontal_count, 'get_bottom_sources', &
         'call count does not match number of (unmasked) horizontal points')

      do ivar = 1, size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(ivar)%save .and. &
               model%horizontal_diagnostic_variables(ivar)%target%source == source_do_bottom) then
            pdata_hz => model%get_horizontal_diagnostic_data(ivar)
            call check_horizontal(pdata_hz, model%horizontal_diagnostic_variables(ivar)%missing_value, &
               -model%horizontal_diagnostic_variables(ivar)%missing_value)
         end if
      end do

      call report_test_result()

      ! ======================================================================
      ! Postprocressing
      ! ======================================================================

      call start_test('finalize_outputs')
      call model%finalize_outputs()

      call assert(column_loop_count == interior_count, 'finalize_outputs', &
         'call count does not match number of (unmasked) interior points')

      do ivar = 1, size(model%interior_diagnostic_variables)
         if (model%interior_diagnostic_variables(ivar)%save .and. model%interior_diagnostic_variables(ivar)%target%source == source_do_column) then
            pdata => model%get_interior_diagnostic_data(ivar)
            call check_interior(pdata, model%interior_diagnostic_variables(ivar)%missing_value, &
               -model%interior_diagnostic_variables(ivar)%missing_value)
         end if
      end do

      do ivar = 1, size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(ivar)%save .and. &
               model%horizontal_diagnostic_variables(ivar)%target%source == source_do_column) then
            pdata_hz => model%get_horizontal_diagnostic_data(ivar)
            call check_horizontal(pdata_hz, model%horizontal_diagnostic_variables(ivar)%missing_value, &
               -model%horizontal_diagnostic_variables(ivar)%missing_value)
         end if
      end do

      call report_test_result()

      ! ======================================================================
      ! Retrieve vertical velocities (sinking, floating, active movement).
      ! ======================================================================

      call start_test('get_vertical_movement')
      vertical_movement_loop_count = 0
      _BEGIN_OUTER_INTERIOR_LOOP_
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_) && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
         ! We are looping over depth, but as we have a non-constant bottom index (yet no mask), we need to skip everything below bottom
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         _START_ = bottom_index _INDEX_HORIZONTAL_LOCATION_
#  else
         _STOP_ = bottom_index _INDEX_HORIZONTAL_LOCATION_
#  endif
#endif
         call model%get_vertical_movement(_PREARG_INTERIOR_IN_ w _INTERIOR_SLICE_RANGE_PLUS_1_)
         do ivar = 1, size(model%interior_state_variables)
            if (mod(ivar, 2) == 0) then
               call check_interior_slice_plus_1(w _INTERIOR_SLICE_RANGE_PLUS_1_, ivar, 0.0_rke, real(ivar + interior_state_offset, rke) &
                  _POSTARG_INTERIOR_IN_)
            else
               call check_interior_slice_plus_1(w _INTERIOR_SLICE_RANGE_PLUS_1_, ivar, 0.0_rke, -real(ivar + interior_state_offset, rke) &
                  _POSTARG_INTERIOR_IN_)
            end if
         end do
      _END_OUTER_INTERIOR_LOOP_

#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_) && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
      _START_ = domain_start(_FABM_VECTORIZED_DIMENSION_INDEX_)
      _STOP_ = domain_stop(_FABM_VECTORIZED_DIMENSION_INDEX_)
#  endif

      call assert(vertical_movement_loop_count == interior_count, 'get_vertical_movement', &
         'call count does not match number of (unmasked) interior points')
      call report_test_result()

      ! ======================================================================
      ! Check state with valid state
      ! ======================================================================

      call start_test('check_interior_state')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%check_interior_state(_PREARG_INTERIOR_IN_  .true., valid)
         if (.not. valid) call driver%fatal_error('check_interior_state', 'state is reported as invalid')
      _END_OUTER_INTERIOR_LOOP_
      call report_test_result()

      call start_test('check_surface_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%check_surface_state(_PREARG_HORIZONTAL_IN_ .true., valid)
         if (.not. valid) call driver%fatal_error('check_surface_state', 'state is reported as invalid')
      _END_OUTER_HORIZONTAL_LOOP_
      call report_test_result()

      call start_test('check_bottom_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%check_bottom_state(_PREARG_HORIZONTAL_IN_ .true., valid)
         if (.not. valid) call driver%fatal_error('check_bottom_state', 'state is reported as invalid')
      _END_OUTER_HORIZONTAL_LOOP_
      call report_test_result()

      ! ======================================================================
      ! Check state with state below minimum
      ! ======================================================================

      ! Now destroy the state by setting all values to below the minimum
      do ivar = 1, size(model%interior_state_variables)
        interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar) = model%interior_state_variables(ivar)%minimum - abs(model%interior_state_variables(ivar)%minimum) - 1
        call apply_mask_3d(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), model%interior_state_variables(ivar)%missing_value)
      end do
      do ivar = 1, size(model%bottom_state_variables)
        bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = model%bottom_state_variables(ivar)%minimum - abs(model%bottom_state_variables(ivar)%minimum) - 1
        call apply_mask_2d(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%bottom_state_variables(ivar)%missing_value)
      end do
      do ivar = 1, size(model%surface_state_variables)
        surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = model%surface_state_variables(ivar)%minimum - abs(model%surface_state_variables(ivar)%minimum) - 1
        call apply_mask_2d(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%surface_state_variables(ivar)%missing_value)
      end do

      call start_test('check_interior_state < min')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%check_interior_state(_PREARG_INTERIOR_IN_ .true., valid)
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(_START_:_STOP_))), 'check_interior_state', 'invalid result')
#  else
         call assert(valid .neqv. any(_IS_UNMASKED_(mask _INDEX_GLOBAL_INTERIOR_(_START_:_STOP_))), 'check_interior_state', 'invalid result')
#  endif
#else
         call assert(.not. valid, 'check_interior_state', 'invalid result')
#endif
      _END_OUTER_INTERIOR_LOOP_
      do ivar = 1, size(model%interior_state_variables)
         call check_interior(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), &
            model%interior_state_variables(ivar)%missing_value, model%interior_state_variables(ivar)%minimum)
      end do
      call report_test_result()

      call start_test('check_surface_state < min')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%check_surface_state(_PREARG_HORIZONTAL_IN_ .true., valid)
#ifdef _HAS_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(_START_:_STOP_))), 'check_surface_state', 'invalid result')
#else
         call assert(.not. valid, 'check_surface_state', 'invalid result')
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar = 1, size(model%surface_state_variables)
         call check_horizontal(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), &
            model%surface_state_variables(ivar)%missing_value, model%surface_state_variables(ivar)%minimum)
      end do
      call report_test_result()

      call start_test('check_bottom_state < min')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%check_bottom_state(_PREARG_HORIZONTAL_IN_ .true., valid)
#ifdef _HAS_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(_START_:_STOP_))), 'check_bottom_state', 'invalid result')
#else
         call assert(.not. valid, 'check_bottom_state', 'invalid result')
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar = 1, size(model%bottom_state_variables)
         call check_horizontal(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), &
            model%bottom_state_variables(ivar)%missing_value, model%bottom_state_variables(ivar)%minimum)
      end do
      call report_test_result()

      ! ======================================================================
      ! Check state with state above maximum
      ! ======================================================================

      ! Now destroy the state by setting all values to above the maximum
      do ivar = 1, size(model%interior_state_variables)
        interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar) = model%interior_state_variables(ivar)%maximum + abs(model%interior_state_variables(ivar)%maximum) + 1
        call apply_mask_3d(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), model%interior_state_variables(ivar)%missing_value)
      end do
      do ivar = 1, size(model%bottom_state_variables)
        bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = model%bottom_state_variables(ivar)%maximum + abs(model%bottom_state_variables(ivar)%maximum) + 1
        call apply_mask_2d(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%bottom_state_variables(ivar)%missing_value)
      end do
      do ivar = 1, size(model%surface_state_variables)
        surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = model%surface_state_variables(ivar)%maximum + abs(model%surface_state_variables(ivar)%maximum) + 1
        call apply_mask_2d(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%surface_state_variables(ivar)%missing_value)
      end do

      call start_test('check_interior_state > max')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%check_interior_state(_PREARG_INTERIOR_IN_ .true., valid)
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(_START_:_STOP_))), 'check_interior_state', 'invalid result')
#  else
         call assert(valid .neqv. any(_IS_UNMASKED_(mask _INDEX_GLOBAL_INTERIOR_(_START_:_STOP_))), 'check_interior_state', 'invalid result')
#  endif
#else
         call assert(.not. valid, 'check_interior_state', 'invalid result')
#endif
      _END_OUTER_INTERIOR_LOOP_
      do ivar = 1, size(model%interior_state_variables)
         call check_interior(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), &
            model%interior_state_variables(ivar)%missing_value, model%interior_state_variables(ivar)%maximum)
      end do
      call report_test_result()

      call start_test('check_surface_state > max')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%check_surface_state(_PREARG_HORIZONTAL_IN_ .true., valid)
#ifdef _HAS_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(_START_:_STOP_))), 'check_surface_state', 'invalid result')
#else
         call assert(.not. valid, 'check_surface_state', 'invalid result')
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar = 1, size(model%surface_state_variables)
         call check_horizontal(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), &
            model%surface_state_variables(ivar)%missing_value, model%surface_state_variables(ivar)%maximum)
      end do
      call report_test_result()

      call start_test('check_bottom_state > max')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%check_bottom_state(_PREARG_HORIZONTAL_IN_ .true., valid)
#ifdef _HAS_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(_START_:_STOP_))), 'check_bottom_state', 'invalid result')
#else
         call assert(.not. valid, 'check_bottom_state', 'invalid result')
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar = 1, size(model%bottom_state_variables)
         call check_horizontal(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), &
            model%bottom_state_variables(ivar)%missing_value, model%bottom_state_variables(ivar)%maximum)
      end do
      call report_test_result()

      ! ======================================================================
      ! Retrieve totals of conserved quantities
      ! ======================================================================

      call start_test('get_interior_conserved_quantities')
      _BEGIN_OUTER_INTERIOR_LOOP_
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_) && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
         ! We are looping over depth, but as we have a non-constant bottom index (yet no mask), we need to skip everything below bottom
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         _START_ = bottom_index _INDEX_HORIZONTAL_LOCATION_
#  else
         _STOP_ = bottom_index _INDEX_HORIZONTAL_LOCATION_
#  endif
#endif
         call model%get_interior_conserved_quantities(_PREARG_INTERIOR_IN_ total_int _INTERIOR_SLICE_RANGE_PLUS_1_)
         do ivar = 1, size(model%conserved_quantities)
            call check_interior_slice_plus_1(total_int _INTERIOR_SLICE_RANGE_PLUS_1_, ivar, model%conserved_quantities(ivar)%missing_value, &
               (interior_state_offset + 0.5_rke * (test_model%nstate + 1)) * test_model%nstate _POSTARG_INTERIOR_IN_)
         end do
      _END_OUTER_INTERIOR_LOOP_

#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_) && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
      _START_ = domain_start(_FABM_VECTORIZED_DIMENSION_INDEX_)
      _STOP_ = domain_stop(_FABM_VECTORIZED_DIMENSION_INDEX_)
#  endif
      call report_test_result()

      call start_test('get_horizontal_conserved_quantities')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= domain_start(_FABM_DEPTH_DIMENSION_INDEX_) .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)) then
#endif
         call model%get_horizontal_conserved_quantities(_PREARG_HORIZONTAL_IN_ total_hz _HORIZONTAL_SLICE_RANGE_PLUS_1_)
         do ivar = 1, size(model%conserved_quantities)
            call check_horizontal_slice_plus_1(total_hz _HORIZONTAL_SLICE_RANGE_PLUS_1_, ivar, model%conserved_quantities(ivar)%missing_value, &
               (surface_state_offset + 0.5_rke * (test_model%nsurface_state + 1)) * test_model%nsurface_state + &
               (bottom_state_offset + 0.5_rke * (test_model%nbottom_state + 1)) * test_model%nbottom_state _POSTARG_HORIZONTAL_IN_)
         end do
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         endif
#endif
      _END_OUTER_HORIZONTAL_LOOP_

   end subroutine test_update

   subroutine start_test(name)
      character(len=*),intent(in) :: name
      write (*,'(2X,A,"...")',advance='NO') name
   end subroutine

   subroutine report_test_result()
      write (*,'(1X,A)') 'SUCCESS'
   end subroutine

   subroutine assert(condition, source, message)
      logical,          intent(in) :: condition
      character(len=*), intent(in) :: source, message
      if (.not. condition) &
         call driver%fatal_error(source, message)
   end subroutine

   subroutine apply_mask_3d(dat,missing_value)
      real(rke) _DIMENSION_GLOBAL_,intent(inout) :: dat
      real(rke),                   intent(in)    :: missing_value
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      integer :: j__
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         if (.not. _IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)) dat _INDEX_GLOBAL_VERTICAL_(:) = missing_value
      _END_GLOBAL_HORIZONTAL_LOOP_
#  else
      where (.not._IS_UNMASKED_(mask)) dat = missing_value
#  endif
#endif
   end subroutine

   subroutine apply_mask_2d(dat, missing_value)
      real(rke) _DIMENSION_GLOBAL_HORIZONTAL_,intent(inout) :: dat
      real(rke),                              intent(in)    :: missing_value
#ifdef _HAS_MASK_
      where (.not. _IS_UNMASKED_(mask_hz)) dat = missing_value
#endif
   end subroutine

   subroutine check_interior_slice_plus_1(dat, index, required_masked_value, required_value _POSTARG_INTERIOR_IN_)
      _DECLARE_ARGUMENTS_INTERIOR_IN_
      real(rke) _DIMENSION_EXT_SLICE_PLUS_1_,intent(in) :: dat
      integer,                               intent(in) :: index
      real(rke),                             intent(in) :: required_masked_value, required_value
#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_interior_slice(dat(:,index), required_masked_value, required_value _POSTARG_INTERIOR_IN_)
#else
      call check_interior_slice(dat(index), required_masked_value, required_value _POSTARG_INTERIOR_IN_)
#endif
   end subroutine

   subroutine check_interior_slice(slice_data, required_masked_value, required_value _POSTARG_INTERIOR_IN_)
      real(rke) _DIMENSION_EXT_SLICE_,intent(in) :: slice_data
      real(rke),                      intent(in) :: required_masked_value, required_value
      _DECLARE_ARGUMENTS_INTERIOR_IN_

      logical :: unmasked
      real(rke) :: value
#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      integer :: _ITERATOR_
#endif

      unmasked = .true.

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      do _ITERATOR_ = _START_, _STOP_
      value = slice_data(_ITERATOR_)
#else
      value = slice_data
#endif

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      unmasked = _IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)
#  else
      unmasked = _IS_UNMASKED_(mask _INDEX_LOCATION_)
#  endif
#endif

      if (unmasked) then
         call assert(value == required_value, 'check_interior_slice', 'one or more non-masked cells do not have the value required.')
      else
         call assert(value == required_masked_value, 'check_interior_slice', 'one or more masked cells do not have the value required.')
      end if

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      end do
#endif
   end subroutine

   subroutine check_horizontal_slice_plus_1(dat,index, required_masked_value, required_value _POSTARG_HORIZONTAL_IN_)
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_PLUS_1_, intent(in) :: dat
      integer,                                           intent(in) :: index
      real(rke),                                         intent(in) :: required_masked_value, required_value
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
#ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_horizontal_slice(dat(:,index), required_masked_value, required_value _POSTARG_HORIZONTAL_IN_)
#else
      call check_horizontal_slice(dat(index), required_masked_value, required_value _POSTARG_HORIZONTAL_IN_)
#endif
   end subroutine check_horizontal_slice_plus_1

   subroutine check_horizontal_slice(slice_data, required_masked_value, required_value _POSTARG_HORIZONTAL_IN_)
      real(rke) _DIMENSION_EXT_HORIZONTAL_SLICE_, intent(in) :: slice_data
      real(rke),                                  intent(in) :: required_masked_value, required_value
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_

      logical :: unmasked
      real(rke) :: value
#ifdef _HORIZONTAL_IS_VECTORIZED_
      integer :: _ITERATOR_
#endif

      unmasked = .true.

#ifdef _HORIZONTAL_IS_VECTORIZED_
      do _ITERATOR_ = _START_, _STOP_
      value = slice_data(_ITERATOR_)
#else
      value = slice_data
#endif

#ifdef _HAS_MASK_
      unmasked = _IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)
#endif

      if (unmasked) then
         call assert(value == required_value, 'check_horizontal_slice', 'one or more non-masked cells do not have the value required.')
      else
         call assert(value == required_masked_value, 'check_horizontal_slice', 'one or more masked cells do not have the value required.')
      end if

#ifdef _HORIZONTAL_IS_VECTORIZED_
      end do
#endif
   end subroutine check_horizontal_slice

   subroutine check_interior(dat, required_masked_value, required_value)
      real(rke) _DIMENSION_GLOBAL_, intent(in) :: dat
      real(rke),                    intent(in) :: required_masked_value, required_value
      logical :: unmasked

      unmasked = .true.
      _BEGIN_GLOBAL_LOOP_
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
         unmasked = _IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)
#  else
         unmasked = _IS_UNMASKED_(mask _INDEX_LOCATION_)
#  endif
#elif _FABM_BOTTOM_INDEX_==-1
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         if (_ITERATOR_ < bottom_index _INDEX_HORIZONTAL_LOCATION_) cycle
#  else
         if (_ITERATOR_ > bottom_index _INDEX_HORIZONTAL_LOCATION_) cycle
#  endif
#endif
         if (unmasked) then
            call assert(dat _INDEX_LOCATION_ == required_value, 'check_interior', 'one or more non-masked cells do not have the value required.')
         else
            call assert(dat _INDEX_LOCATION_ == required_masked_value, 'check_interior', 'one or more masked cells do not have the value required.')
         end if
      _END_GLOBAL_LOOP_
   end subroutine

   subroutine check_horizontal(dat, required_masked_value, required_value)
      real(rke) _DIMENSION_GLOBAL_HORIZONTAL_, intent(in) :: dat
      real(rke),                               intent(in) :: required_masked_value, required_value
      logical :: unmasked

      unmasked = .true.
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#ifdef _HAS_MASK_
         unmasked = _IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)
#endif
         if (unmasked) then
            call assert(dat _INDEX_HORIZONTAL_LOCATION_ == required_value, 'check_horizontal', 'one or more non-masked cells do not have the value required.')
         else
            call assert(dat _INDEX_HORIZONTAL_LOCATION_ == required_masked_value, 'check_horizontal', 'one or more masked cells do not have the value required.')
         end if
      _END_GLOBAL_HORIZONTAL_LOOP_
   end subroutine

end program

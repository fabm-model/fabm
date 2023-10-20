#include "fabm_driver.h"
#include "fabm_private.h"

#undef _BEGIN_GLOBAL_LOOP_
#undef _END_GLOBAL_LOOP_
#undef _BEGIN_OUTER_INTERIOR_LOOP_
#undef _END_OUTER_INTERIOR_LOOP_
#undef _BEGIN_OUTER_HORIZONTAL_LOOP_
#undef _END_OUTER_HORIZONTAL_LOOP_

program test_omp

   use omp_lib, only: omp_get_max_threads

   use fabm
   use fabm_driver
   use fabm_parameters, only: rke
   use fabm_types, only: source_do, source_do_surface, source_do_bottom, source_do_column, get_free_unit
   use fabm_omp

   use test_shared

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
#  define _END_GLOBAL_LOOP_ end do;end do;end do !;i__=domain_extent(1);j__=domain_extent(2);k__=domain_extent(3)
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do k__=domain_start(3),domain_stop(3);do j__=domain_start(2),domain_stop(2)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do !;j__=domain_extent(2);k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do k__=domain_start(3),domain_stop(3);do i__=domain_start(1),domain_stop(1)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do !;i__=domain_extent(1);k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do j__=domain_start(2),domain_stop(2);do i__=domain_start(1),domain_stop(1)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do !;i__=domain_extent(1);j__=domain_extent(2)
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
#  define _END_OUTER_INTERIOR_LOOP_ end do !;j__=domain_extent(2)
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
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do !;j__=domain_extent(2);k__=domain_extent(3)
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do k__=domain_start(3),domain_stop(3)
#    define _END_OUTER_HORIZONTAL_LOOP_ end do !;k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do j__=domain_start(2),domain_stop(2)
#    define _END_OUTER_HORIZONTAL_LOOP_ end do !;j__=domain_extent(2)
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
#  define _DEPTH_INDEX_  bottom_index _INDEX_HORIZONTAL_LOCATION_
#elif defined(_FABM_DEPTH_DIMENSION_INDEX_)
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
#    define _DEPTH_INDEX_  domain_start(_FABM_DEPTH_DIMENSION_INDEX_)
#  else
#    define _DEPTH_INDEX_  domain_stop(_FABM_DEPTH_DIMENSION_INDEX_)
#  endif
#endif

   real(rke), allocatable, target _DIMENSION_GLOBAL_PLUS_1_            :: interior_state
   real(rke), allocatable, target _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: surface_state
   real(rke), allocatable, target _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: bottom_state

   real(rke), allocatable _DIMENSION_SLICE_PLUS_1_            :: dy, w, total_int
   real(rke), allocatable _DIMENSION_HORIZONTAL_SLICE_PLUS_1_ :: flux, sms_sf, sms_bt, total_hz

   real(rke), allocatable, target _DIMENSION_GLOBAL_            :: temperature
   real(rke), allocatable, target _DIMENSION_GLOBAL_            :: depth
   real(rke), allocatable, target _DIMENSION_GLOBAL_HORIZONTAL_ :: wind_speed

   class (type_fabm_omp_model), pointer :: model

   integer :: domain_extent(_FABM_DIMENSION_COUNT_)
   integer :: domain_start(_FABM_DIMENSION_COUNT_)
   integer :: domain_stop(_FABM_DIMENSION_COUNT_)

   character(len=1024) :: version
   character(len=20) :: arg
   integer :: ivar
   integer :: i
   integer :: ntest = -1
   logical :: no_mask = .false.
   logical :: no_diag = .false.
   character(len=1024) :: yaml_path
   character(len=1024) :: env_path

#if _FABM_DIMENSION_COUNT_>0
   i__ = 50
#endif
#if _FABM_DIMENSION_COUNT_>1
   j__ = 40
#endif
#if _FABM_DIMENSION_COUNT_>2
   k__ = 45
#endif
   yaml_path = 'fabm.yaml'
   env_path = 'environment.yaml'

   call fabm_get_version(version)
   write (*,'(a,a)') 'FABM ', trim(version)

   ! Parse command line arguments
   i = 1
   do
      call get_command_argument(i, arg)
      if (arg == '') exit
      select case (arg)
      case ('-c')
         i = i + 1
         call get_command_argument(i, yaml_path)
      case ('-e')
         i = i + 1
         call get_command_argument(i, env_path)
      case ('--nomask')
         no_mask = .true.
      case ('--nodiag')
         no_diag = .true.
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
         write (*,'(a)') '-c <FILE>:     path to FABM configuration file (default: fabm.yaml)'
         write (*,'(a)') '-e <FILE>:     path to environment file (default: environment.yaml)'
         write (*,'(a)') '-n <N>:        number of replicates when simulating'
         write (*,'(a)') '--nomask:      unmask all points when simulating'
         write (*,'(a)') '--nodiag:      flag all diagnostics as not required by the host'
         stop 0
      case default
         write (*,'(a)') 'ERROR'
         write (*,'(a)') ''
         write (*,'(a)') 'Unknown command line argument: ' // trim(arg)
         write (*,'(a)') 'To see supported arguments, run with -h.'
         stop 2
      end select
      i = i + 1
   end do

#if _FABM_DIMENSION_COUNT_>0
   domain_extent = (/ _LOCATION_ /)
   domain_start(:) = 1
   domain_stop(:) = domain_extent
   interior_count = product(domain_stop - domain_start + 1)
#else
   interior_count = 1
#endif

   ! Set defaults
   if (ntest == -1) ntest = 50000000 / interior_count

#ifdef _INTERIOR_IS_VECTORIZED_
   _START_ = domain_start(_FABM_VECTORIZED_DIMENSION_INDEX_)
   _STOP_ = domain_stop(_FABM_VECTORIZED_DIMENSION_INDEX_)
#endif

   allocate(tmp _INDEX_LOCATION_)
   allocate(tmp_hz _INDEX_HORIZONTAL_LOCATION_)

   allocate(type_test_driver::driver)

   call fabm_initialize_library()

   model => fabm_create_omp_model(path=yaml_path)

   if (no_diag) then
      model%interior_diagnostic_variables%save = .false.
      model%horizontal_diagnostic_variables%save = .false.
   end if

   call model%set_domain(_PREARG_LOCATION_ seconds_per_time_unit=1._rke)

#ifdef _HAS_MASK_
   allocate(mask_hz _INDEX_HORIZONTAL_LOCATION_)
#  ifdef _FABM_HORIZONTAL_MASK_
   call model%set_mask(mask_hz)
#  else
   allocate(mask _INDEX_LOCATION_)
   call model%set_mask(mask, mask_hz)
#  endif
#endif

   ! ======================================================================
   ! Specify vertical indices of surface and bottom.
   ! ======================================================================

#if defined(_FABM_DEPTH_DIMENSION_INDEX_)&&_FABM_BOTTOM_INDEX_==-1
   allocate(bottom_index _INDEX_HORIZONTAL_LOCATION_)
   call model%set_bottom_index(bottom_index)
#endif

   allocate(interior_state(_PREARG_LOCATION_ size(model%interior_state_variables)))
   allocate(surface_state(_PREARG_HORIZONTAL_LOCATION_ size(model%surface_state_variables)))
   allocate(bottom_state(_PREARG_HORIZONTAL_LOCATION_ size(model%bottom_state_variables)))

   ! ======================================================================
   ! Send pointers to state variable data to FABM.
   ! ======================================================================

   call model%link_all_interior_state_data(interior_state(_PREARG_LOCATION_DIMENSIONS_ :))
   call model%link_all_surface_state_data(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ :))
   call model%link_all_bottom_state_data(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ :))

   call read_environment(model, env_path _POSTARG_LOCATION_)

   call model%start()

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

   print *, 'Using ', omp_get_max_threads(), ' threads'
   call simulate(ntest)

   call model%finalize()

   deallocate(model)

   call fabm_finalize_library()

   contains

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

      call count_active_points()
   end subroutine configure_mask

   subroutine count_active_points()
      logical :: active
      integer :: i, nhz, nhz_active, nint, nint_active

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

      nint = 1
      nint_active = 1
      nhz = 1
      nhz_active = 1
      do i = 1, _FABM_DIMENSION_COUNT_
         nint = nint * domain_extent(i)
         nint_active = nint_active * (domain_stop(i) - domain_start(i) + 1)
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
         if (i /= _FABM_DEPTH_DIMENSION_INDEX_) then
#endif
         nhz = nhz * domain_extent(i)
         nhz_active = nhz_active * (domain_stop(i) - domain_start(i) + 1)
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
         end if
#endif
      end do
      write (*,'(a,i0,a,i0,a,i0,a)') 'Interior: ', nint, ' points, ', nint_active, ' in active range, ', interior_count, ' unmasked'
      write (*,'(a,i0,a,i0,a,i0,a)') 'Horizontal: ', nhz, ' points, ', nhz_active, ' in active range, ', horizontal_count, ' unmasked'
   end subroutine

   subroutine simulate(n)
      integer, intent(in) :: n
      real(rke) :: time_begin, time_end
      integer :: nseed
      integer, allocatable :: seed(:)
      integer :: ireport
      logical, parameter :: repair = .false.
      logical :: valid
      integer(8) :: clock, clock_start, ticks_per_sec

      ireport = n / 10

      call random_seed(size=nseed)
      allocate(seed(nseed))
      seed(:) = 1
      call random_seed(put=seed)

      call configure_range(.false.)
      call configure_mask(.not. no_mask)

#if (_FABM_DIMENSION_COUNT_==0||(_FABM_DIMENSION_COUNT_==1&&_FABM_VECTORIZED_DIMENSION_INDEX_==1))
#  define _NO_OMP_INTERIOR_
#endif

#if (_HORIZONTAL_DIMENSION_COUNT_==0||(_HORIZONTAL_DIMENSION_COUNT_==1&&_FABM_VECTORIZED_DIMENSION_INDEX_!=_FABM_DEPTH_DIMENSION_INDEX_))
#  define _NO_OMP_HORIZONTAL_
#endif

      !$OMP PARALLEL DEFAULT(SHARED)

#ifdef _NO_OMP_INTERIOR_
      !$OMP SINGLE
#else
      !$OMP DO
#endif
      _BEGIN_OUTER_INTERIOR_LOOP_
         call model%initialize_interior_state(_ARG_INTERIOR_IN_)
      _END_OUTER_INTERIOR_LOOP_
#ifdef _NO_OMP_INTERIOR_
      !$OMP END SINGLE
#else
      !$OMP END DO
#endif

#ifdef _NO_OMP_HORIZONTAL_
      !$OMP SINGLE
#else
      !$OMP DO
#endif
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%initialize_bottom_state(_ARG_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
#ifdef _NO_OMP_HORIZONTAL_
      !$OMP END SINGLE
#else
      !$OMP END DO
#endif

#ifdef _NO_OMP_HORIZONTAL_
      !$OMP SINGLE
#else
      !$OMP DO
#endif
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call model%initialize_surface_state(_ARG_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
#ifdef _NO_OMP_HORIZONTAL_
      !$OMP END SINGLE
#else
      !$OMP END DO
#endif

      !$OMP END PARALLEL

      write (*,'(a,i0,a)') 'Simulating with ', interior_count, ' wet cells...'

      call cpu_time(time_begin)
      call system_clock(clock_start)

      do i = 1, n
         call model%prepare_inputs()

         !$OMP PARALLEL DEFAULT(SHARED)

#ifdef _NO_OMP_HORIZONTAL_
         !$OMP SINGLE
#else
         !$OMP DO
#endif
         _BEGIN_OUTER_HORIZONTAL_LOOP_
            flux = 0
            sms_bt = 0
            call model%get_bottom_sources(_PREARG_HORIZONTAL_IN_ flux, sms_bt)
         _END_OUTER_HORIZONTAL_LOOP_
#ifdef _NO_OMP_HORIZONTAL_
         !$OMP END SINGLE
#else
         !$OMP END DO
#endif

#ifdef _NO_OMP_HORIZONTAL_
         !$OMP SINGLE
#else
         !$OMP DO
#endif
         _BEGIN_OUTER_HORIZONTAL_LOOP_
            flux = 0
            sms_sf = 0
            call model%get_surface_sources(_PREARG_HORIZONTAL_IN_ flux, sms_sf)
         _END_OUTER_HORIZONTAL_LOOP_
#ifdef _NO_OMP_HORIZONTAL_
         !$OMP END SINGLE
#else
         !$OMP END DO
#endif

#ifdef _NO_OMP_INTERIOR_
         !$OMP SINGLE
#else
         !$OMP DO
#endif
         _BEGIN_OUTER_INTERIOR_LOOP_
            dy = 0
            call model%get_interior_sources(_PREARG_INTERIOR_IN_ dy)
         _END_OUTER_INTERIOR_LOOP_
#ifdef _NO_OMP_INTERIOR_
         !$OMP END SINGLE
#else
         !$OMP END DO
#endif

         !$OMP END PARALLEL

         call model%finalize_outputs()

         !$OMP PARALLEL DEFAULT(SHARED)

#ifdef _NO_OMP_HORIZONTAL_
         !$OMP SINGLE
#else
         !$OMP DO
#endif
         _BEGIN_OUTER_HORIZONTAL_LOOP_
            call model%check_bottom_state(_PREARG_HORIZONTAL_IN_ repair, valid)
         _END_OUTER_HORIZONTAL_LOOP_
#ifdef _NO_OMP_HORIZONTAL_
         !$OMP END SINGLE
#else
         !$OMP END DO
#endif

#ifdef _NO_OMP_HORIZONTAL_
         !$OMP SINGLE
#else
         !$OMP DO
#endif
         _BEGIN_OUTER_HORIZONTAL_LOOP_
            call model%check_surface_state(_PREARG_HORIZONTAL_IN_ repair, valid)
         _END_OUTER_HORIZONTAL_LOOP_
#ifdef _NO_OMP_HORIZONTAL_
         !$OMP END SINGLE
#else
         !$OMP END DO
#endif

#ifdef _NO_OMP_INTERIOR_
         !$OMP SINGLE
#else
         !$OMP DO
#endif
         _BEGIN_OUTER_INTERIOR_LOOP_
            call model%check_interior_state(_PREARG_INTERIOR_IN_ repair, valid)
         _END_OUTER_INTERIOR_LOOP_
#ifdef _NO_OMP_INTERIOR_
         !$OMP END SINGLE
#else
         !$OMP END DO
#endif

         !$OMP END PARALLEL
         if (mod(i, ireport) == 0) write (*,'(i0,a)') int(100*i/real(n, rke)), ' % complete'
      end do

      call system_clock( count=clock, count_rate=ticks_per_sec)
      call cpu_time(time_end)

      write (*,'(a)') 'Simulation complete.'
      write (*,'(a,f8.3,a)') 'CPU time: ', time_end - time_begin, ' s'
      write (*,'(a,f8.3,a)') 'Wall time: ', (clock - clock_start) / real(ticks_per_sec), ' s'
   end subroutine

end program

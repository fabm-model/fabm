#include "fabm_driver.h"
#include "fabm_private.h"

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
use fabm_config
use fabm_driver

use test_models
use host_hooks

implicit none

#if _FABM_DIMENSION_COUNT_>0
integer :: _LOCATION_
#endif

#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
!  No mask but variable bottom index. Index of depth dimension must be 1.
!  All loops over inner dimension should skip points below bottom.
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
#    define _IMIN_ bottom_index _INDEX_HORIZONTAL_LOCATION_
#    define _IMAX_ domain_extent(1)
#  else
#    define _IMIN_ 1
#    define _IMAX_ bottom_index _INDEX_HORIZONTAL_LOCATION_
#  endif
#else
!  Loops over inner dimension should span full domain
#  define _IMIN_ 1
#  define _IMAX_ domain_extent(1)
#endif
#define _IRANGE_ _IMIN_,_IMAX_

#if _FABM_DIMENSION_COUNT_==0
#  define _BEGIN_GLOBAL_LOOP_
#  define _END_GLOBAL_LOOP_
#elif _FABM_DIMENSION_COUNT_==1
#  define _BEGIN_GLOBAL_LOOP_ do i__=_IRANGE_
#  define _END_GLOBAL_LOOP_ end do;i__=domain_extent(1)
#  ifdef _FABM_DEPTH_DIMENSION_INDEX_
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#    define _END_GLOBAL_HORIZONTAL_LOOP_
#  endif
#elif _FABM_DIMENSION_COUNT_==2
#  define _BEGIN_GLOBAL_LOOP_ do j__=1,domain_extent(2);do i__=_IRANGE_
#  define _END_GLOBAL_LOOP_ end do;end do;i__=domain_extent(1);j__=domain_extent(2)
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do j__=1,domain_extent(2)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;j__=domain_extent(2)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do i__=_IRANGE_
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;i__=domain_extent(1)
#  endif
#elif _FABM_DIMENSION_COUNT_==3
#  define _BEGIN_GLOBAL_LOOP_ do k__=1,domain_extent(3);do j__=1,domain_extent(2);do i__=_IRANGE_
#  define _END_GLOBAL_LOOP_ end do;end do;end do;i__=domain_extent(1);j__=domain_extent(2);k__=domain_extent(3)
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do k__=1,domain_extent(3);do j__=1,domain_extent(2)
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do;j__=domain_extent(2);k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do k__=1,domain_extent(3);do i__=_IRANGE_
#    define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do;i__=domain_extent(1);k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do j__=1,domain_extent(2);do i__=_IRANGE_
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
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do j__=1,domain_extent(2)
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
#  define _BEGIN_OUTER_INTERIOR_LOOP_ do k__=1,domain_extent(3);do j__=1,domain_extent(2)
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do;j__=domain_extent(2);k__=domain_extent(3)
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do k__=1,domain_extent(3)
#    define _END_OUTER_HORIZONTAL_LOOP_ end do;k__=domain_extent(3)
#  elif _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do j__=1,domain_extent(2)
#    define _END_OUTER_HORIZONTAL_LOOP_ end do;j__=domain_extent(2)
#  else
     ! No horizontal dimension vectorized; do full outer loop.
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#    define _END_OUTER_HORIZONTAL_LOOP_ _END_GLOBAL_HORIZONTAL_LOOP_
#  endif
#endif

#ifdef _INTERIOR_IS_VECTORIZED_
integer :: loop_start, loop_stop
#endif

real(rk),allocatable _DIMENSION_GLOBAL_ :: tmp
real(rk),allocatable _DIMENSION_GLOBAL_HORIZONTAL_ :: tmp_hz

#ifdef _HAS_MASK_
_FABM_MASK_TYPE_,allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_ :: mask_hz
#  ifndef _FABM_HORIZONTAL_MASK_
_FABM_MASK_TYPE_,allocatable,target _DIMENSION_GLOBAL_ :: mask
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
integer,allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_ :: bottom_index
#endif

real(rk),allocatable,target _DIMENSION_GLOBAL_PLUS_1_            :: interior_state
real(rk),allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: surface_state
real(rk),allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_PLUS_1_ :: bottom_state

real(rk),allocatable _DIMENSION_SLICE_PLUS_1_        :: dy
real(rk),allocatable _DIMENSION_SLICE_PLUS_1_        :: w
real(rk),allocatable _DIMENSION_HORIZONTAL_SLICE_PLUS_1_ :: flux,sms_sf,sms_bt

real(rk),allocatable,target _DIMENSION_GLOBAL_            :: temperature
real(rk),allocatable,target _DIMENSION_GLOBAL_            :: depth
real(rk),allocatable,target _DIMENSION_GLOBAL_HORIZONTAL_ :: wind_speed

type (type_model) :: model

class (type_test_model), pointer :: test_model

integer :: domain_extent(_FABM_DIMENSION_COUNT_)

character(len=20) :: arg
integer :: ivar
integer :: i
integer :: mode = 1
integer :: ntest = -1

! Parse command line arguments
call start_test('parsing command line arguments')
i = 1
do
   call get_command_argument(i, arg)
   if (arg == '') exit
   select case (arg)
   case ('-s', '--simulate')
      mode = 2
      i = i + 1
   case ('-n')
      call get_command_argument(i + 1, arg)
      read (arg,*) ntest
      i = i + 2
   case default
      write (*,'(a)') 'Unknown command line argument: ' // trim(arg)
      stop 2
   end select
end do
call report_test_result()

#if _FABM_DIMENSION_COUNT_>0
i__=50
#endif
#if _FABM_DIMENSION_COUNT_>1
j__=40
#endif
#if _FABM_DIMENSION_COUNT_>2
k__=45
#endif

#if _FABM_DIMENSION_COUNT_>0
domain_extent = (/ _LOCATION_ /)
interior_count = product(domain_extent)
#else
interior_count = 1
#endif

! Set defaults
if (ntest == -1) then
   if (mode == 1) then
      ntest = 1
   else
      ntest = 50000000/interior_count
   end if
end if

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
horizontal_count = interior_count / domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)
#else
horizontal_count = interior_count
#endif

#ifdef _INTERIOR_IS_VECTORIZED_
loop_start = 1
loop_stop = domain_extent(_FABM_VECTORIZED_DIMENSION_INDEX_)
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
    allocate(test_model)
    call model%root%add_child(test_model,'test_model','test model',configunit=-1)
case (2)
    ! Test with user-provided fabm.yaml
    call fabm_create_model_from_yaml_file(model, do_not_initialize=.true.)
end select
call report_test_result()

call start_test('fabm_initialize')
call fabm_initialize(model)
call report_test_result()

! ======================================================================
! Provide extents of the spatial domain.
! ======================================================================

call start_test('fabm_set_domain')
call fabm_set_domain(model _ARGUMENTS_LOCATION_)
call report_test_result()

! ======================================================================
! Set up spatial mask.
! ======================================================================

#ifdef _HAS_MASK_
allocate(mask_hz _INDEX_HORIZONTAL_LOCATION_)
call start_test('fabm_set_mask')
#  ifdef _FABM_HORIZONTAL_MASK_
call fabm_set_mask(model,mask_hz)
#  else
allocate(mask _INDEX_LOCATION_)
call fabm_set_mask(model,mask,mask_hz)
#  endif
call report_test_result()
#endif

! ======================================================================
! Specify vertical indices of surface and bottom.
! ======================================================================

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
call start_test('set_surface_index')
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
call model%set_surface_index(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_))
#  else
call model%set_surface_index(1)
#  endif
call report_test_result()

call start_test('set_bottom_index')
#  if _FABM_BOTTOM_INDEX_==-1
allocate(bottom_index _INDEX_HORIZONTAL_LOCATION_)
call model%set_bottom_index(bottom_index)
#  elif defined(_FABM_VERTICAL_BOTTOM_TO_SURFACE_)
call model%set_bottom_index(1)
#  else
call model%set_bottom_index(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_))
#  endif
call report_test_result()
#endif

allocate(interior_state(_PREARG_LOCATION_ size(model%state_variables)))
allocate(surface_state(_PREARG_HORIZONTAL_LOCATION_ size(model%surface_state_variables)))
allocate(bottom_state(_PREARG_HORIZONTAL_LOCATION_ size(model%bottom_state_variables)))

! ======================================================================
! Send pointers to state variable data to FABM.
! ======================================================================

call start_test('link_interior_state_data')
do ivar=1,size(model%state_variables)
   call model%link_interior_state_data(ivar,interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar))
end do
call report_test_result()

call start_test('link_surface_state_data')
do ivar=1,size(model%surface_state_variables)
   call model%link_surface_state_data(ivar,surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar))
end do
call report_test_result()

call start_test('link_bottom_state_data')
do ivar=1,size(model%bottom_state_variables)
   call model%link_bottom_state_data(ivar,bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar))
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
    call model%link_interior_data(standard_variables%temperature,temperature)
    call model%link_interior_data(standard_variables%depth,depth)
    call report_test_result()

    call start_test('link_horizontal_data')
    call model%link_horizontal_data(standard_variables%wind_speed,wind_speed)
    call report_test_result()
case (2)
    call read_environment
end select

! ======================================================================
! Check whether FABM has all dependencies fulfilled
! (i.e., whether all required calls for fabm_link_*_data have been made)
! ======================================================================

call start_test('fabm_check_ready')
call fabm_check_ready(model)
call report_test_result()

#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
allocate(dy(loop_start:loop_stop,size(model%state_variables)))
allocate(w(loop_start:loop_stop,size(model%state_variables)))
#else
allocate(dy(size(model%state_variables)))
allocate(w(size(model%state_variables)))
#endif

#if defined(_FABM_VECTORIZED_DIMENSION_INDEX_)&&(_FABM_DEPTH_DIMENSION_INDEX_!=_FABM_VECTORIZED_DIMENSION_INDEX_)
allocate(flux(loop_start:loop_stop,size(model%state_variables)))
allocate(sms_sf(loop_start:loop_stop,size(model%surface_state_variables)))
allocate(sms_bt(loop_start:loop_stop,size(model%bottom_state_variables)))
#else
allocate(flux(size(model%state_variables)))
allocate(sms_sf(size(model%surface_state_variables)))
allocate(sms_bt(size(model%bottom_state_variables)))
#endif

select case (mode)
case (1)
   do i=1,ntest
      call test_update
   end do
case(2)
   call simulate(ntest)
end select

contains

   subroutine read_environment
      use yaml, only: yaml_parse => parse, yaml_error_length => error_length
      use yaml_types, only: type_node, type_yaml_dictionary => type_dictionary, type_yaml_scalar => type_scalar, type_yaml_key_value_pair => type_key_value_pair

      integer, parameter :: yaml_unit = 100
      character(yaml_error_length) :: yaml_error
      class (type_node),pointer :: yaml_root
      type (type_yaml_key_value_pair), pointer :: yaml_pair
      real(rk) :: value
      logical :: success
      type type_input
         type (type_bulk_variable_id)                        :: interior_id
         type (type_horizontal_variable_id)                  :: horizontal_id
         type (type_scalar_variable_id)                      :: scalar_id
         real(rk), allocatable _DIMENSION_GLOBAL_            :: interior_data
         real(rk), allocatable _DIMENSION_GLOBAL_HORIZONTAL_ :: horizontal_data
         real(rk)                                            :: scalar_data
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
                  value = node%to_real(0._rk, success)
                  if (.not. success) then
                     call driver%log_message('Cannot parse '//trim(node%string)//' as real.')
                     stop 2
                  end if
                  allocate(input)
                  input%interior_id = model%get_bulk_variable_id(trim(yaml_pair%key))
                  if (fabm_is_variable_used(input%interior_id)) then
                      allocate(input%interior_data _INDEX_LOCATION_)
                      input%interior_data = value
                      call model%link_interior_data(input%interior_id, input%interior_data)
                  else
                      input%horizontal_id = model%get_horizontal_variable_id(trim(yaml_pair%key))
                      if (fabm_is_variable_used(input%horizontal_id)) then
                         allocate(input%horizontal_data _INDEX_HORIZONTAL_LOCATION_)
                         input%horizontal_data = value
                         call model%link_horizontal_data(input%horizontal_id, input%horizontal_data)
                      else
                         input%scalar_id = model%get_scalar_variable_id(trim(yaml_pair%key))
                         if (fabm_is_variable_used(input%scalar_id)) then
                            input%scalar_data = value
                            call model%link_scalar(input%scalar_id, input%scalar_data)
                         else
                            call driver%log_message('WARNING: environment variable '//trim(yaml_pair%key)//' is not used by FABM model and will be ignored.')
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

   subroutine randomize_mask
#if _FABM_BOTTOM_INDEX_==-1
      ! Depth index of bottom varies in the horizontal - pick random numbers between 0 (land) and maximum index
      call random_number(tmp_hz)
      bottom_index = floor(tmp_hz*(1+domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)))
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
     ! Ensure invalid bottom indices [land points] are set such that vertical loops have 0 iterations.
     where (bottom_index == 0) bottom_index = domain_extent(_FABM_DEPTH_DIMENSION_INDEX_) + 1
#  endif
#endif

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      ! Apply random mask across horizontal domain (half of grid cells masked)
      call random_number(tmp_hz)
      mask_hz = _FABM_UNMASKED_VALUE_
      where (tmp_hz>0.5_rk) mask_hz = _FABM_MASKED_VALUE_
      horizontal_count = count(_IS_UNMASKED_(mask_hz))
#    ifdef _FABM_DEPTH_DIMENSION_INDEX_
      interior_count = horizontal_count * domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)
#    else
      interior_count = horizontal_count
#    endif
#  else
      ! Apply random mask across interior domain (half of grid cells masked)
      call random_number(tmp)
      mask = _FABM_UNMASKED_VALUE_
      where (tmp>0.5_rk) mask = _FABM_MASKED_VALUE_

#    if _FABM_BOTTOM_INDEX_==-1
      ! Bottom index varies in the horizontal. Ensure the bottom cell itself is unmasked, and anything deeper is masked.
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         ! Valid bottom index - unmask associated cell, then mask all deeper ones
#      ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) mask _INDEX_GLOBAL_VERTICAL_(bottom_index _INDEX_HORIZONTAL_LOCATION_) = _FABM_UNMASKED_VALUE_
         mask _INDEX_GLOBAL_VERTICAL_(:bottom_index _INDEX_HORIZONTAL_LOCATION_ - 1) = _FABM_MASKED_VALUE_
#      else
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= 1) mask _INDEX_GLOBAL_VERTICAL_(bottom_index _INDEX_HORIZONTAL_LOCATION_) = _FABM_UNMASKED_VALUE_
         mask _INDEX_GLOBAL_VERTICAL_(bottom_index _INDEX_HORIZONTAL_LOCATION_ + 1:) = _FABM_MASKED_VALUE_
#      endif
      _END_GLOBAL_HORIZONTAL_LOOP_
#    endif

      ! Infer horizontal mask (mask points that have all column layers masked)
      mask_hz = _FABM_UNMASKED_VALUE_
      where (.not.any(_IS_UNMASKED_(mask),_FABM_DEPTH_DIMENSION_INDEX_)) mask_hz = _FABM_MASKED_VALUE_

      ! For valid points in the horizontal, make sure that index 1 (bottom or surface) is unmasked
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         if (_IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)) then
            mask _INDEX_GLOBAL_VERTICAL_(1) = _FABM_UNMASKED_VALUE_
#    if _FABM_BOTTOM_INDEX_!=-1
            mask _INDEX_GLOBAL_VERTICAL_(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) = _FABM_UNMASKED_VALUE_
#    endif
         end if
      _END_GLOBAL_HORIZONTAL_LOOP_
      horizontal_count = count(_IS_UNMASKED_(mask_hz))
      interior_count = count(_IS_UNMASKED_(mask))
#  endif
#elif _FABM_BOTTOM_INDEX_==-1
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
      horizontal_count = count(bottom_index <= domain_extent(_FABM_DEPTH_DIMENSION_INDEX_))
      interior_count = sum(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_) - bottom_index + 1)
#  else
      horizontal_count = count(bottom_index >= 1)
      interior_count = sum(bottom_index)
#  endif
#endif
   end subroutine randomize_mask

   subroutine simulate(n)
      integer, intent(in) :: n
      real(rk) :: time_begin, time_end
      integer :: nseed
      integer, allocatable :: seed(:)

      call random_seed(size=nseed)
      allocate(seed(nseed))
      seed(:) = 1
      call random_seed(put=seed)

      call randomize_mask

      call start_test('fabm_initialize_state')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call fabm_initialize_state(model _ARGUMENTS_INTERIOR_IN_)
      _END_OUTER_INTERIOR_LOOP_
      call report_test_result()

      call start_test('fabm_initialize_bottom_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_initialize_bottom_state(model _ARGUMENTS_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
      call report_test_result()

      call start_test('fabm_initialize_surface_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_initialize_surface_state(model _ARGUMENTS_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
      call report_test_result()

      write (*,'(a,i0,a)') 'Simulating with ', interior_count, ' wet cells...'

      call cpu_time(time_begin)

      do i=1,n
         _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
            ! No mask but non-constant bottom index. We need to skip everything below bottom
            if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= 1 .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) &
               call fabm_get_light(model,_IRANGE_ _ARG_VERTICAL_FIXED_LOCATION_)
#  else
            call fabm_get_light(model,1,domain_extent(_FABM_DEPTH_DIMENSION_INDEX_) _ARG_VERTICAL_FIXED_LOCATION_)
#  endif
#else
            call fabm_get_light(model _ARGUMENTS_HORIZONTAL_IN_)
#endif
         _END_GLOBAL_HORIZONTAL_LOOP_

         _BEGIN_OUTER_HORIZONTAL_LOOP_
            flux = 0
            sms_bt = 0
            call fabm_do_bottom(model _ARGUMENTS_HORIZONTAL_IN_,flux,sms_bt)
         _END_OUTER_HORIZONTAL_LOOP_

         _BEGIN_OUTER_HORIZONTAL_LOOP_
            flux = 0
            sms_sf = 0
            call fabm_do_surface(model _ARGUMENTS_HORIZONTAL_IN_,flux,sms_sf)
         _END_OUTER_HORIZONTAL_LOOP_

         _BEGIN_OUTER_INTERIOR_LOOP_
            dy = 0
            call fabm_do(model _ARGUMENTS_INTERIOR_IN_,dy)
         _END_OUTER_INTERIOR_LOOP_

         if (mod(i, 100) == 0) write (*,'(i0,a)') int(100*i/real(n, rk)), ' % complete'
      end do

      call cpu_time(time_end)

      write (*,'(a)') 'Simulation complete.'
      write (*,'(a,f8.3,a)') 'Total time: ', time_end - time_begin, ' s'
   end subroutine

   subroutine test_update
      real(rk),pointer _DIMENSION_GLOBAL_ :: pdata
      real(rk),pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: pdata_hz
      logical :: valid

      call randomize_mask

      ! ======================================================================
      ! Initialize all state variables
      ! ======================================================================

      call start_test('fabm_initialize_state')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call fabm_initialize_state(model _ARGUMENTS_INTERIOR_IN_)
      _END_OUTER_INTERIOR_LOOP_
      do ivar=1,size(model%state_variables)
         call check_interior(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), model%state_variables(ivar)%missing_value, ivar+interior_state_offset+1._rk)
      end do
      call report_test_result()

      call start_test('fabm_initialize_bottom_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_initialize_bottom_state(model _ARGUMENTS_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar=1,size(model%bottom_state_variables)
         call check_horizontal(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%bottom_state_variables(ivar)%missing_value, ivar+bottom_state_offset+1._rk)
      end do
      call report_test_result()

      call start_test('fabm_initialize_surface_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_initialize_surface_state(model _ARGUMENTS_HORIZONTAL_IN_)
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar=1,size(model%surface_state_variables)
         call check_horizontal(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%surface_state_variables(ivar)%missing_value, ivar+surface_state_offset+1._rk)
      end do
      call report_test_result()

      ! ======================================================================
      ! Initialize environmental dependencies
      ! ======================================================================

      temperature = 1+interior_dependency_offset
      call apply_mask_3d(temperature,-999._rk-interior_dependency_offset)

      wind_speed = 1+horizontal_dependency_offset
      call apply_mask_2d(wind_speed,-999._rk-horizontal_dependency_offset)

#ifdef _FABM_DEPTH_DIMENSION_INDEX_
      ! Model has depth dimension: make sure depth varies from 0 at the surface till 1 at the bottom
      _BEGIN_GLOBAL_LOOP_
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
#    if _FABM_BOTTOM_INDEX_==-1
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) then
            depth _INDEX_LOCATION_ = 2
         else
            depth _INDEX_LOCATION_ = real(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-_VERTICAL_ITERATOR_,rk)/(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-bottom_index _INDEX_HORIZONTAL_LOCATION_)
         end if
#    else
         depth _INDEX_LOCATION_ = real(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-_VERTICAL_ITERATOR_,rk)/(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-1)
#    endif
#  else
#    if _FABM_BOTTOM_INDEX_==-1
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_==1) then
            depth _INDEX_LOCATION_ = 2
         else
            depth _INDEX_LOCATION_ = real(_VERTICAL_ITERATOR_-1,rk)/(bottom_index _INDEX_HORIZONTAL_LOCATION_-1)
         end if
#    else
         depth _INDEX_LOCATION_ = real(_VERTICAL_ITERATOR_-1,rk)/(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-1)
#    endif
#  endif
      _END_GLOBAL_LOOP_
#else
      ! No depth dimension
      depth = 2
#endif
      call apply_mask_3d(depth,-999._rk-interior_dependency_offset)

      do ivar=1,size(model%state_variables)
         interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar) = ivar+interior_state_offset
         call apply_mask_3d(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar),model%state_variables(ivar)%missing_value)
      end do
      do ivar=1,size(model%surface_state_variables)
         surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = ivar+surface_state_offset
         call apply_mask_2d(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar),model%surface_state_variables(ivar)%missing_value)
      end do
      do ivar=1,size(model%bottom_state_variables)
         bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = ivar+bottom_state_offset
         call apply_mask_2d(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar),model%bottom_state_variables(ivar)%missing_value)
      end do

      ! ======================================================================
      ! Retrieve source terms of interior state variables.
      ! ======================================================================

      call start_test('fabm_do')
      loop_count = 0
      _BEGIN_OUTER_INTERIOR_LOOP_
         dy = 0
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_) && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
         ! We are looping over depth, but as we have a non-constant bottom index (yet no mask), we need to skip everything below bottom
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= 1 .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) &
            call fabm_do(model,_IMIN_,_IMAX_ _ARG_INTERIOR_FIXED_LOCATION_,dy(_IMIN_:_IMAX_,:))
#else
         call fabm_do(model _ARGUMENTS_INTERIOR_IN_,dy)
#endif
         do ivar=1,size(model%state_variables)
            call check_interior_slice_plus_1(dy,ivar,0.0_rk,-real(ivar+interior_state_offset,rk) _ARGUMENTS_INTERIOR_IN_)
         end do
      _END_OUTER_INTERIOR_LOOP_
      call assert(loop_count == interior_count, 'fabm_do', 'call count does not match number of (unmasked) interior points')

      do ivar=1,size(model%diagnostic_variables)
         if (model%diagnostic_variables(ivar)%save .and. model%diagnostic_variables(ivar)%target%source==source_do) then
            pdata => fabm_get_interior_diagnostic_data(model, ivar)
            call check_interior(pdata, model%diagnostic_variables(ivar)%missing_value, -model%diagnostic_variables(ivar)%missing_value)
         end if
      end do

      call report_test_result()

      ! ======================================================================
      ! Retrieve surface fluxes of interior state variables, source terms of surface-associated state variables.
      ! ======================================================================

#if _FABM_BOTTOM_INDEX_==-1
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) depth _INDEX_GLOBAL_VERTICAL_(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) = 0
#  else
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == 1) depth _INDEX_GLOBAL_VERTICAL_(1) = 0
#  endif
      _END_GLOBAL_HORIZONTAL_LOOP_
#endif

      call start_test('fabm_do_surface')
      loop_count = 0
      _BEGIN_OUTER_HORIZONTAL_LOOP_
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= 1 .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) then
#endif
         flux = 0
         sms_sf = 0
         call fabm_do_surface(model _ARGUMENTS_HORIZONTAL_IN_,flux,sms_sf)
         do ivar=1,size(model%state_variables)
            call check_horizontal_slice_plus_1(flux,ivar,0.0_rk,-real(ivar+interior_state_offset,rk) _ARGUMENTS_HORIZONTAL_IN_)
         end do
         do ivar=1,size(model%surface_state_variables)
            call check_horizontal_slice_plus_1(sms_sf,ivar,0.0_rk,-real(ivar+surface_state_offset,rk) _ARGUMENTS_HORIZONTAL_IN_)
         end do
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         endif
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      call assert(loop_count == horizontal_count, 'fabm_do_surface', 'call count does not match number of (unmasked) horizontal points')

      do ivar=1,size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(ivar)%save .and. model%horizontal_diagnostic_variables(ivar)%target%source == source_do_surface) then
            pdata_hz => fabm_get_horizontal_diagnostic_data(model, ivar)
            call check_horizontal(pdata_hz, model%horizontal_diagnostic_variables(ivar)%missing_value, -model%horizontal_diagnostic_variables(ivar)%missing_value)
         end if
      end do

      call report_test_result()

      ! ======================================================================
      ! Retrieve bottom fluxes of interior state variables, source terms of bottom-associated state variables.
      ! ======================================================================

#if _FABM_BOTTOM_INDEX_==-1
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) depth _INDEX_GLOBAL_VERTICAL_(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) = 1
#  else
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ == 1) depth _INDEX_GLOBAL_VERTICAL_(1) = 1
#  endif
      _END_GLOBAL_HORIZONTAL_LOOP_
#endif

      call start_test('fabm_do_bottom')
      loop_count = 0
      _BEGIN_OUTER_HORIZONTAL_LOOP_
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= 1 .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) then
#endif
         flux = 0
         sms_bt = 0
         call fabm_do_bottom(model _ARGUMENTS_HORIZONTAL_IN_,flux,sms_bt)
         do ivar=1,size(model%state_variables)
            call check_horizontal_slice_plus_1(flux,ivar,0.0_rk,-real(ivar+interior_state_offset,rk) _ARGUMENTS_HORIZONTAL_IN_)
         end do
         do ivar=1,size(model%bottom_state_variables)
            call check_horizontal_slice_plus_1(sms_bt,ivar,0.0_rk,-real(ivar+bottom_state_offset,rk) _ARGUMENTS_HORIZONTAL_IN_)
         end do
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         endif
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      call assert(loop_count == horizontal_count, 'fabm_do_surface', 'call count does not match number of (unmasked) horizontal points')

      do ivar=1,size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(ivar)%save .and. model%horizontal_diagnostic_variables(ivar)%target%source == source_do_bottom) then
            pdata_hz => fabm_get_horizontal_diagnostic_data(model, ivar)
            call check_horizontal(pdata_hz, model%horizontal_diagnostic_variables(ivar)%missing_value, -model%horizontal_diagnostic_variables(ivar)%missing_value)
         end if
      end do

      call report_test_result()

      ! ======================================================================
      ! Column-based operators
      ! ======================================================================

      call start_test('fabm_get_light')
      loop_count = 0
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
#  if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_)
         ! No mask but non-constant bottom index. We need to skip everything below bottom
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= 1 .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) &
            call fabm_get_light(model,_IMIN_,_IMAX_ _ARG_VERTICAL_FIXED_LOCATION_)
#  else
         call fabm_get_light(model,1,domain_extent(_FABM_DEPTH_DIMENSION_INDEX_) _ARG_VERTICAL_FIXED_LOCATION_)
#  endif
#else
         call fabm_get_light(model _ARGUMENTS_HORIZONTAL_IN_)
#endif
      _END_GLOBAL_HORIZONTAL_LOOP_
      call assert(loop_count == interior_count, 'fabm_get_light', 'call count does not match number of (unmasked) interior points')

      do ivar=1,size(model%diagnostic_variables)
         if (model%diagnostic_variables(ivar)%save .and. model%diagnostic_variables(ivar)%target%source==source_do_column) then
            pdata => fabm_get_interior_diagnostic_data(model, ivar)
            call check_interior(pdata, model%diagnostic_variables(ivar)%missing_value, -model%diagnostic_variables(ivar)%missing_value)
         end if
      end do

      do ivar=1,size(model%horizontal_diagnostic_variables)
         if (model%horizontal_diagnostic_variables(ivar)%save .and. model%horizontal_diagnostic_variables(ivar)%target%source == source_do_column) then
            pdata_hz => fabm_get_horizontal_diagnostic_data(model, ivar)
            call check_horizontal(pdata_hz, model%horizontal_diagnostic_variables(ivar)%missing_value, -model%horizontal_diagnostic_variables(ivar)%missing_value)
         end if
      end do

      call report_test_result()

      ! ======================================================================
      ! Retrieve vertical velocities (sinking, floating, active movement).
      ! ======================================================================

      call start_test('fabm_get_vertical_movement')
      loop_count = 0
      _BEGIN_OUTER_INTERIOR_LOOP_
#if _FABM_BOTTOM_INDEX_==-1 && !defined(_HAS_MASK_) && _FABM_VECTORIZED_DIMENSION_INDEX_==_FABM_DEPTH_DIMENSION_INDEX_ && defined(_FABM_DEPTH_DIMENSION_INDEX_)
         ! We are looping over depth, but as we have a non-constant bottom index (yet no mask), we need to skip everything below bottom
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_ >= 1 .and. bottom_index _INDEX_HORIZONTAL_LOCATION_ <= domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)) &
            call fabm_get_vertical_movement(model,_IMIN_,_IMAX_ _ARG_INTERIOR_FIXED_LOCATION_,w(_IMIN_:_IMAX_,:))
#else
         call fabm_get_vertical_movement(model _ARGUMENTS_INTERIOR_IN_,w)
#endif
         do ivar=1,size(model%state_variables)
            if (mod(ivar, 2) == 0) then
               call check_interior_slice_plus_1(w,ivar,0.0_rk,real(ivar+interior_state_offset,rk) _ARGUMENTS_INTERIOR_IN_)
            else
               call check_interior_slice_plus_1(w,ivar,0.0_rk,-real(ivar+interior_state_offset,rk) _ARGUMENTS_INTERIOR_IN_)
            end if
         end do
      _END_OUTER_INTERIOR_LOOP_
      call assert(loop_count == interior_count, 'fabm_get_vertical_movement', 'call count does not match number of (unmasked) interior points')
      call report_test_result()

      ! ======================================================================
      ! Check state with valid state
      ! ======================================================================

      call start_test('fabm_check_state')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call fabm_check_state(model _ARGUMENTS_INTERIOR_IN_, .true., valid)
         if (.not. valid) call driver%fatal_error('fabm_check_state', 'state is reported as invalid')
      _END_OUTER_INTERIOR_LOOP_
      call report_test_result()

      call start_test('fabm_check_surface_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_check_surface_state(model _ARGUMENTS_HORIZONTAL_IN_, .true., valid)
         if (.not. valid) call driver%fatal_error('fabm_check_surface_state', 'state is reported as invalid')
      _END_OUTER_HORIZONTAL_LOOP_
      call report_test_result()

      call start_test('fabm_check_bottom_state')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_check_bottom_state(model _ARGUMENTS_HORIZONTAL_IN_, .true., valid)
         if (.not. valid) call driver%fatal_error('fabm_check_bottom_state', 'state is reported as invalid')
      _END_OUTER_HORIZONTAL_LOOP_
      call report_test_result()

      ! ======================================================================
      ! Check state with state below minimum
      ! ======================================================================

      ! Now destroy the state by setting all values to below the minimum
      do ivar=1,size(model%state_variables)
        interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar) = model%state_variables(ivar)%minimum - abs(model%state_variables(ivar)%minimum) - 1
        call apply_mask_3d(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), model%state_variables(ivar)%missing_value)
      end do
      do ivar=1,size(model%bottom_state_variables)
        bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = model%bottom_state_variables(ivar)%minimum - abs(model%bottom_state_variables(ivar)%minimum) - 1
        call apply_mask_2d(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%bottom_state_variables(ivar)%missing_value)
      end do
      do ivar=1,size(model%surface_state_variables)
        surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = model%surface_state_variables(ivar)%minimum - abs(model%surface_state_variables(ivar)%minimum) - 1
        call apply_mask_2d(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%surface_state_variables(ivar)%missing_value)
      end do

      call start_test('fabm_check_state < min')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call fabm_check_state(model _ARGUMENTS_INTERIOR_IN_, .true., valid)
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), 'fabm_check_state', 'invalid result')
#  else
         call assert(valid .neqv. any(_IS_UNMASKED_(mask _INDEX_GLOBAL_INTERIOR_(loop_start:loop_stop))), 'fabm_check_state', 'invalid result')
#  endif
#else
         call assert(.not. valid, 'fabm_check_state', 'invalid result')
#endif
      _END_OUTER_INTERIOR_LOOP_
      do ivar=1,size(model%state_variables)
         call check_interior(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), model%state_variables(ivar)%missing_value, model%state_variables(ivar)%minimum)
      end do
      call report_test_result()

      call start_test('fabm_check_surface_state < min')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_check_surface_state(model _ARGUMENTS_HORIZONTAL_IN_, .true., valid)
#ifdef _HAS_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), 'fabm_check_surface_state', 'invalid result')
#else
         call assert(.not. valid, 'fabm_check_surface_state', 'invalid result')
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar=1,size(model%surface_state_variables)
         call check_horizontal(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%surface_state_variables(ivar)%missing_value, model%surface_state_variables(ivar)%minimum)
      end do
      call report_test_result()

      call start_test('fabm_check_bottom_state < min')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_check_bottom_state(model _ARGUMENTS_HORIZONTAL_IN_, .true., valid)
#ifdef _HAS_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), 'fabm_check_bottom_state', 'invalid result')
#else
         call assert(.not. valid, 'fabm_check_bottom_state', 'invalid result')
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar=1,size(model%bottom_state_variables)
         call check_horizontal(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%bottom_state_variables(ivar)%missing_value, model%bottom_state_variables(ivar)%minimum)
      end do
      call report_test_result()

      ! ======================================================================
      ! Check state with state above maximum
      ! ======================================================================

      ! Now destroy the state by setting all values to above the maximum
      do ivar=1,size(model%state_variables)
        interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar) = model%state_variables(ivar)%maximum + abs(model%state_variables(ivar)%maximum) + 1
        call apply_mask_3d(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), model%state_variables(ivar)%missing_value)
      end do
      do ivar=1,size(model%bottom_state_variables)
        bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = model%bottom_state_variables(ivar)%maximum + abs(model%bottom_state_variables(ivar)%maximum) + 1
        call apply_mask_2d(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%bottom_state_variables(ivar)%missing_value)
      end do
      do ivar=1,size(model%surface_state_variables)
        surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = model%surface_state_variables(ivar)%maximum + abs(model%surface_state_variables(ivar)%maximum) + 1
        call apply_mask_2d(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%surface_state_variables(ivar)%missing_value)
      end do

      call start_test('fabm_check_state > max')
      _BEGIN_OUTER_INTERIOR_LOOP_
         call fabm_check_state(model _ARGUMENTS_INTERIOR_IN_, .true., valid)
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), 'fabm_check_state', 'invalid result')
#  else
         call assert(valid .neqv. any(_IS_UNMASKED_(mask _INDEX_GLOBAL_INTERIOR_(loop_start:loop_stop))), 'fabm_check_state', 'invalid result')
#  endif
#else
         call assert(.not. valid, 'fabm_check_state', 'invalid result')
#endif
      _END_OUTER_INTERIOR_LOOP_
      do ivar=1,size(model%state_variables)
         call check_interior(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar), model%state_variables(ivar)%missing_value, model%state_variables(ivar)%maximum)
      end do
      call report_test_result()

      call start_test('fabm_check_surface_state > max')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_check_surface_state(model _ARGUMENTS_HORIZONTAL_IN_, .true., valid)
#ifdef _HAS_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), 'fabm_check_surface_state', 'invalid result')
#else
         call assert(.not. valid, 'fabm_check_surface_state', 'invalid result')
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar=1,size(model%surface_state_variables)
         call check_horizontal(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%surface_state_variables(ivar)%missing_value, model%surface_state_variables(ivar)%maximum)
      end do
      call report_test_result()

      call start_test('fabm_check_bottom_state > max')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         call fabm_check_bottom_state(model _ARGUMENTS_HORIZONTAL_IN_, .true., valid)
#ifdef _HAS_MASK_
         call assert(valid .neqv. any(_IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), 'fabm_check_bottom_state', 'invalid result')
#else
         call assert(.not. valid, 'fabm_check_bottom_state', 'invalid result')
#endif
      _END_OUTER_HORIZONTAL_LOOP_
      do ivar=1,size(model%bottom_state_variables)
         call check_horizontal(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar), model%bottom_state_variables(ivar)%missing_value, model%bottom_state_variables(ivar)%maximum)
      end do
      call report_test_result()

   end subroutine test_update

   subroutine start_test(name)
      character(len=*),intent(in) :: name
      write (*,'(A,"...")',advance='NO') name
   end subroutine

   subroutine report_test_result()
      write (*,'(X,A)') 'SUCCESS'
   end subroutine

   subroutine assert(condition, source, message)
      logical,          intent(in) :: condition
      character(len=*), intent(in) :: source, message
      if (.not. condition) &
         call driver%fatal_error(source, message)
   end subroutine

   subroutine apply_mask_3d(dat,missing_value)
      real(rk) _DIMENSION_GLOBAL_,intent(inout) :: dat
      real(rk),                   intent(in)    :: missing_value
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
      real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,intent(inout) :: dat
      real(rk),                              intent(in)    :: missing_value
#ifdef _HAS_MASK_
      where (.not. _IS_UNMASKED_(mask_hz)) dat = missing_value
#endif
   end subroutine

   subroutine check_interior_slice_plus_1(dat, index, required_masked_value, required_value _ARGUMENTS_INTERIOR_IN_)
      real(rk) _DIMENSION_EXT_SLICE_PLUS_1_,intent(in) :: dat
      integer,                              intent(in) :: index
      real(rk),                             intent(in) :: required_masked_value, required_value
      _DECLARE_ARGUMENTS_INTERIOR_IN_
#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_interior_slice(dat(:,index), required_masked_value, required_value _ARGUMENTS_INTERIOR_IN_)
#else
      call check_interior_slice(dat(index), required_masked_value, required_value _ARGUMENTS_INTERIOR_IN_)
#endif
   end subroutine

   subroutine check_interior_slice(slice_data, required_masked_value, required_value _ARGUMENTS_INTERIOR_IN_)
      real(rk) _DIMENSION_EXT_SLICE_,intent(in) :: slice_data
      real(rk),                      intent(in) :: required_masked_value, required_value
      _DECLARE_ARGUMENTS_INTERIOR_IN_

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      call assert(all(slice_data == required_masked_value .or.       _IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), &
         'check_interior_slice', 'one or more masked cells do not have the value required.')
      call assert(all(slice_data == required_value        .or. .not. _IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), &
         'check_interior_slice', 'one or more non-masked cells do not have the value required.')
#  else
      call assert(all(slice_data == required_masked_value .or.       _IS_UNMASKED_(mask _INDEX_GLOBAL_INTERIOR_(loop_start:loop_stop))), &
         'check_interior_slice', 'one or more masked cells do not have the value required.')
      call assert(all(slice_data == required_value        .or. .not. _IS_UNMASKED_(mask _INDEX_GLOBAL_INTERIOR_(loop_start:loop_stop))), &
         'check_interior_slice', 'one or more non-masked cells do not have the value required.')
#  endif
#elif defined(_INTERIOR_IS_VECTORIZED_)
      call assert(all(slice_data(_IMIN_:_IMAX_) == required_value), 'check_interior_slice', 'one or more cells do not have the value required.')
#else
      call assert(slice_data == required_value, 'check_interior_slice', 'variable does not have the value required.')
#endif
   end subroutine

   subroutine check_horizontal_slice_plus_1(dat,index, required_masked_value, required_value _ARGUMENTS_HORIZONTAL_IN_)
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_, intent(in) :: dat
      integer,                                      intent(in) :: index
      real(rk),                                     intent(in) :: required_masked_value, required_value
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
#ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_horizontal_slice(dat(:,index), required_masked_value, required_value _ARGUMENTS_HORIZONTAL_IN_)
#else
      call check_horizontal_slice(dat(index), required_masked_value, required_value _ARGUMENTS_HORIZONTAL_IN_)
#endif
   end subroutine check_horizontal_slice_plus_1

   subroutine check_horizontal_slice(slice_data, required_masked_value, required_value _ARGUMENTS_HORIZONTAL_IN_)
      real(rk) _DIMENSION_HORIZONTAL_SLICE_, intent(in) :: slice_data
      real(rk),                              intent(in) :: required_masked_value, required_value
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_

#ifdef _HORIZONTAL_IS_VECTORIZED_
#  ifdef _HAS_MASK_
      call assert(all(slice_data == required_masked_value .or.       _IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), 'check_horizontal_slice', 'one or more masked cells do not have the value required.')
      call assert(all(slice_data == required_value        .or. .not. _IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop))), 'check_horizontal_slice', 'one or more non-masked cells do not have the value required.')
#  else
      call assert(all(slice_data == required_value), 'check_horizontal_slice', 'one or more cells do not have the value required.')
#  endif
#else
      call assert(slice_data == required_value, 'check_horizontal_slice', 'returned scalar does not have the value required.')
#endif
   end subroutine check_horizontal_slice

   subroutine check_interior(dat, required_masked_value, required_value)
      real(rk) _DIMENSION_GLOBAL_,intent(in) :: dat
      real(rk),                   intent(in) :: required_masked_value, required_value
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      integer :: j__
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         if (_IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)) then
            call assert(all(dat _INDEX_GLOBAL_VERTICAL_(:) == required_value), 'check_interior', 'one or more non-masked cells do not have the value required.')
         else
            call assert(all(dat _INDEX_GLOBAL_VERTICAL_(:) == required_masked_value), 'check_interior', 'one or more masked cells do not have the value required.')
         end if
      _END_GLOBAL_HORIZONTAL_LOOP_
#  else
      call assert(all(dat == required_masked_value .or.       _IS_UNMASKED_(mask)), 'check_interior', 'one or more masked cells do not have the value required.')
      call assert(all(dat == required_value        .or. .not. _IS_UNMASKED_(mask)), 'check_interior', 'one or more non-masked cells do not have the value required.')
#  endif
#elif _FABM_DIMENSION_COUNT_>0
#  if _FABM_BOTTOM_INDEX_==-1
      ! Skip points below bottom
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#    ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
         call assert(all(dat _INDEX_GLOBAL_VERTICAL_(bottom_index _INDEX_HORIZONTAL_LOCATION_:) == required_value), 'check_interior', 'one or more masked cells do not have the value required.')
#    else
         call assert(all(dat _INDEX_GLOBAL_VERTICAL_(1:bottom_index _INDEX_HORIZONTAL_LOCATION_) == required_value), 'check_interior', 'one or more masked cells do not have the value required.')
#    endif
      _END_GLOBAL_HORIZONTAL_LOOP_
#  else
      call assert(all(dat == required_value), 'check_interior', 'one or more cells do not have the value required.')
#  endif
#else
      call assert(dat == required_value, 'check_interior', 'variable does not have the value required.')
#endif
   end subroutine

   subroutine check_horizontal(dat, required_masked_value, required_value)
      real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,intent(in) :: dat
      real(rk),                              intent(in) :: required_masked_value, required_value
#ifdef _HAS_MASK_
      call assert(all(dat == required_masked_value .or.       _IS_UNMASKED_(mask_hz)), 'check_horizontal', 'one or more masked cells do not have the value required.')
      call assert(all(dat == required_value        .or. .not. _IS_UNMASKED_(mask_hz)), 'check_horizontal', 'one or more non-masked cells do not have the value required.')
#elif _HORIZONTAL_DIMENSION_COUNT_>0
#  if _FABM_BOTTOM_INDEX_==-1
      ! Skip land points (with bottom index of 0 or max+1)
      call assert(all(dat == required_value .or. bottom_index == 0 .or. bottom_index == domain_extent(_FABM_DEPTH_DIMENSION_INDEX_) + 1), 'check_horizontal', 'one or more cells do not have the value required.')
#  else
      call assert(all(dat == required_value), 'check_horizontal', 'one or more cells do not have the value required.')
#  endif
#else
      call assert(dat == required_value, 'check_horizontal', 'variable does not have the value required.')
#endif
   end subroutine

end program

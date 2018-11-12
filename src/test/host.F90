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

integer,parameter :: ntest = 1

#if _FABM_DIMENSION_COUNT_>0
integer :: _LOCATION_
#endif

#if _FABM_DIMENSION_COUNT_==1
#  define _BEGIN_GLOBAL_LOOP_ do i__=1,domain_extent(1)
#  if _FABM_DEPTH_DIMENSION_INDEX_==1
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_
#    define _BEGIN_OUTER_INTERIOR_LOOP_
#  endif
#  define _END_GLOBAL_LOOP_ end do;i__=domain_extent(1)
#  define _END_GLOBAL_HORIZONTAL_LOOP_
#  define _END_OUTER_HORIZONTAL_LOOP_
#  define _END_OUTER_INTERIOR_LOOP_
#elif _FABM_DIMENSION_COUNT_==2
#  define _BEGIN_GLOBAL_LOOP_ do j__=1,domain_extent(2);do i__=1,domain_extent(1)
#  if _FABM_DEPTH_DIMENSION_INDEX_==2
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do i__=1,domain_extent(1)
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_
#    define _BEGIN_OUTER_INTERIOR_LOOP_ do j__=1,domain_extent(2)
#  endif
#  define _END_GLOBAL_LOOP_ end do;end do;i__=domain_extent(1);j__=domain_extent(2)
#  define _END_GLOBAL_HORIZONTAL_LOOP_ end do;i__=domain_extent(1)
#  define _END_OUTER_HORIZONTAL_LOOP_
#  define _END_OUTER_INTERIOR_LOOP_ end do;j__=domain_extent(2)
#elif _FABM_DIMENSION_COUNT_==3
#  define _BEGIN_GLOBAL_LOOP_ do k__=1,domain_extent(3);do j__=1,domain_extent(2);do i__=1,domain_extent(1)
#  if _FABM_DEPTH_DIMENSION_INDEX_==3
#    define _BEGIN_GLOBAL_HORIZONTAL_LOOP_ do j__=1,domain_extent(2);do i__=1,domain_extent(1)
#    define _BEGIN_OUTER_HORIZONTAL_LOOP_ do j__=1,domain_extent(2)
#    define _BEGIN_OUTER_INTERIOR_LOOP_ do k__=1,domain_extent(3);do j__=1,domain_extent(2)
#  endif
#  define _END_GLOBAL_LOOP_ end do;end do;end do;i__=domain_extent(1);j__=domain_extent(2);k__=domain_extent(3)
#  define _END_GLOBAL_HORIZONTAL_LOOP_ end do;end do;i__=domain_extent(1);j__=domain_extent(2)
#  define _END_OUTER_HORIZONTAL_LOOP_ end do;j__=domain_extent(2)
#  define _END_OUTER_INTERIOR_LOOP_ end do;end do;j__=domain_extent(2);k__=domain_extent(3)
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
#endif

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

integer :: ivar
integer :: i

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
#endif

#ifdef _INTERIOR_IS_VECTORIZED_
loop_start = 1
loop_stop = domain_extent(_FABM_VECTORIZED_DIMENSION_INDEX_)
#endif

allocate(tmp _INDEX_LOCATION_)
allocate(tmp_hz _INDEX_HORIZONTAL_LOCATION_)

allocate(type_test_driver::driver)
call fabm_initialize_library()

allocate(test_model)
call model%root%add_child(test_model,'test_model','test model',configunit=-1)

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
#  ifdef _FABM_VERTICAL_BOTTOM_TO_SURFACE_
call start_test('set_bottom_index')
call model%set_bottom_index(1)
call report_test_result()

call start_test('set_surface_index')
call model%set_surface_index(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_))
call report_test_result()
#  else
call start_test('set_surface_index')
call model%set_surface_index(1)
call report_test_result()

call start_test('set_bottom_index')
#if _FABM_BOTTOM_INDEX_==-1
allocate(bottom_index _INDEX_HORIZONTAL_LOCATION_)
call model%set_bottom_index(bottom_index)
#else
call model%set_bottom_index(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_))
#endif
call report_test_result()
#  endif
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

allocate(depth _INDEX_LOCATION_)
allocate(temperature _INDEX_LOCATION_)
allocate(wind_speed _INDEX_HORIZONTAL_LOCATION_)

! ======================================================================
! Transfer pointers to environmental data
! ======================================================================

call start_test('link_interior_data')
call model%link_interior_data(standard_variables%temperature,temperature)
call model%link_interior_data(standard_variables%depth,depth)
call report_test_result()

call start_test('link_horizontal_data')
call model%link_horizontal_data(standard_variables%wind_speed,wind_speed)
call report_test_result()

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

   do i=1,ntest
      call test_update
   end do

   contains

   subroutine randomize_mask
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      ! Apply random mask across horizontal domain (half of grid cells masked)
      call random_number(tmp_hz)
      mask_hz = _FABM_UNMASKED_VALUE_
      where (tmp_hz>0.5_rk) mask_hz = _FABM_MASKED_VALUE_
#  else
      ! Apply random mask across interior domain (half of grid cells masked)
      call random_number(tmp)
      mask = _FABM_UNMASKED_VALUE_
      where (tmp>0.5_rk) mask = _FABM_MASKED_VALUE_

#    if _FABM_BOTTOM_INDEX_==-1
      ! Depth index of bottom varies in the horizontal - pick random numbers between 0 (land) and maximum index
      call random_number(tmp_hz)
      bottom_index = floor(tmp_hz*(1+domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)))

      ! Based on value for bottom index (0 or higher), either mask all point in the column, or unmask the bottom point
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_==0) then
            ! All land - mask entire column
            mask _INDEX_GLOBAL_VERTICAL_(:) = _FABM_MASKED_VALUE_
         else
            ! Valid bottom index - unmask associated cell
            mask _INDEX_GLOBAL_VERTICAL_(bottom_index _INDEX_HORIZONTAL_LOCATION_) = _FABM_UNMASKED_VALUE_
         end if
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
#  endif
#endif
   end subroutine randomize_mask

   subroutine test_update
      real(rk),pointer _DIMENSION_GLOBAL_ :: pdata
      real(rk),pointer _DIMENSION_GLOBAL_HORIZONTAL_ :: pdata_hz

      call randomize_mask

      ! ======================================================================
      ! Initialize all state variables
      ! ======================================================================

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
         depth _INDEX_LOCATION_ = real(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-_VERTICAL_ITERATOR_,rk)/(domain_extent(_FABM_DEPTH_DIMENSION_INDEX_)-1)
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
         call apply_mask_3d(interior_state(_PREARG_LOCATION_DIMENSIONS_ ivar),-999._rk-interior_state_offset)
      end do
      do ivar=1,size(model%surface_state_variables)
         surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = ivar+surface_state_offset
         call apply_mask_2d(surface_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar),-999._rk-surface_state_offset)
      end do
      do ivar=1,size(model%bottom_state_variables)
         bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar) = ivar+bottom_state_offset
         call apply_mask_2d(bottom_state(_PREARG_HORIZONTAL_LOCATION_DIMENSIONS_ ivar),-999._rk-bottom_state_offset)
      end do

      ! ======================================================================
      ! Retrieve source terms of interior state variables.
      ! ======================================================================

      call start_test('fabm_do')
      _BEGIN_OUTER_INTERIOR_LOOP_
         dy = 0
         call fabm_do(model _ARGUMENTS_INTERIOR_IN_,dy)
         do ivar=1,size(model%state_variables)
            call check_interior_slice_plus_1(dy,ivar,0.0_rk,-real(ivar+interior_state_offset,rk) _ARGUMENTS_INTERIOR_IN_)
         end do
      _END_OUTER_INTERIOR_LOOP_

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
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_==1) depth _INDEX_GLOBAL_VERTICAL_(1) = 0
      _END_GLOBAL_HORIZONTAL_LOOP_
#endif

      call start_test('fabm_do_surface')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         flux = 0
         sms_sf = 0
         call fabm_do_surface(model _ARGUMENTS_HORIZONTAL_IN_,flux,sms_sf)
         do ivar=1,size(model%state_variables)
            call check_horizontal_slice_plus_1(flux,ivar,0.0_rk,-real(ivar+interior_state_offset,rk) _ARGUMENTS_HORIZONTAL_IN_)
         end do
         do ivar=1,size(model%surface_state_variables)
            call check_horizontal_slice_plus_1(sms_sf,ivar,0.0_rk,-real(ivar+surface_state_offset,rk) _ARGUMENTS_HORIZONTAL_IN_)
         end do
      _END_OUTER_HORIZONTAL_LOOP_

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
         if (bottom_index _INDEX_HORIZONTAL_LOCATION_==1) depth _INDEX_GLOBAL_VERTICAL_(1) = 1
      _END_GLOBAL_HORIZONTAL_LOOP_
#endif

      call start_test('fabm_do_bottom')
      _BEGIN_OUTER_HORIZONTAL_LOOP_
         flux = 0
         sms_bt = 0
         call fabm_do_bottom(model _ARGUMENTS_HORIZONTAL_IN_,flux,sms_bt)
         do ivar=1,size(model%state_variables)
            call check_horizontal_slice_plus_1(flux,ivar,0.0_rk,-real(ivar+interior_state_offset,rk) _ARGUMENTS_HORIZONTAL_IN_)
         end do
         do ivar=1,size(model%bottom_state_variables)
            call check_horizontal_slice_plus_1(sms_bt,ivar,0.0_rk,-real(ivar+bottom_state_offset,rk) _ARGUMENTS_HORIZONTAL_IN_)
         end do
      _END_OUTER_HORIZONTAL_LOOP_

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
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
#ifdef _FABM_DEPTH_DIMENSION_INDEX_
         call fabm_get_light(model,1,domain_extent(_FABM_DEPTH_DIMENSION_INDEX_) _ARG_VERTICAL_FIXED_LOCATION_)
#else
         call fabm_get_light(model,_ARGUMENTS_LOCATION_)
#endif
      _END_GLOBAL_HORIZONTAL_LOOP_

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
      _BEGIN_OUTER_INTERIOR_LOOP_
         call fabm_get_vertical_movement(model _ARGUMENTS_INTERIOR_IN_,w)
         do ivar=1,size(model%state_variables)
            call check_interior_slice_plus_1(w,ivar,0.0_rk,-real(ivar+interior_state_offset,rk) _ARGUMENTS_INTERIOR_IN_)
         end do
      _END_OUTER_INTERIOR_LOOP_
      call report_test_result()

   end subroutine test_update

   subroutine start_test(name)
      character(len=*),intent(in) :: name
      write (*,'(A,"...")',advance='NO') name
   end subroutine

   subroutine report_test_result()
      write (*,'(X,A)') 'SUCCESS'
   end subroutine

   subroutine apply_mask_3d(dat,missing_value)
      real(rk) _DIMENSION_GLOBAL_,intent(inout) :: dat
      real(rk),                   intent(in)    :: missing_value
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
      integer :: j__
      _BEGIN_GLOBAL_HORIZONTAL_LOOP_
         if (.not._IS_UNMASKED_(mask_hz _INDEX_HORIZONTAL_LOCATION_)) dat _INDEX_GLOBAL_VERTICAL_(:) = missing_value
      _END_GLOBAL_HORIZONTAL_LOOP_
#  else
      where (.not._IS_UNMASKED_(mask)) dat = missing_value
#  endif
#endif
   end subroutine

   subroutine apply_mask_2d(dat,missing_value)
      real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,intent(inout) :: dat
      real(rk),                              intent(in)    :: missing_value
#ifdef _HAS_MASK_
      where (.not._IS_UNMASKED_(mask_hz)) dat = missing_value
#endif
   end subroutine

   subroutine check_interior_slice_plus_1(dat,index,required_masked_value,required_value _ARGUMENTS_INTERIOR_IN_)
      real(rk) _DIMENSION_EXT_SLICE_PLUS_1_,intent(in) :: dat
      integer,                              intent(in) :: index
      real(rk),                             intent(in) :: required_masked_value,required_value
      _DECLARE_ARGUMENTS_INTERIOR_IN_
#ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
      call check_interior_slice(dat(:,ivar),required_masked_value,required_value _ARGUMENTS_INTERIOR_IN_)
#else
      call check_interior_slice(dat(ivar),required_masked_value,required_value _ARGUMENTS_INTERIOR_IN_)
#endif
   end subroutine

   subroutine check_interior_slice(slice_data,required_masked_value,required_value _ARGUMENTS_INTERIOR_IN_)
      real(rk) _DIMENSION_EXT_SLICE_,intent(in) :: slice_data
      real(rk),                      intent(in) :: required_masked_value,required_value
      _DECLARE_ARGUMENTS_INTERIOR_IN_

#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
#  else
      if (any(slice_data/=required_masked_value.and..not._IS_UNMASKED_(mask _INDEX_GLOBAL_INTERIOR_(loop_start:loop_stop)))) then
         call driver%fatal_error('check_interior_slice','one or more masked cells do not have the value required.')
      end if
      if (any(slice_data/=required_value.and._IS_UNMASKED_(mask _INDEX_GLOBAL_INTERIOR_(loop_start:loop_stop)))) then
         call driver%fatal_error('check_interior_slice','one or more non-masked cells do not have the value required.')
      end if
#  endif
#else
      if (any(slice_data/=required_value)) then
         call driver%fatal_error('check_interior_slice','one or more cells do not have the value required.')
      end if
#endif
   end subroutine

   subroutine check_horizontal_slice_plus_1(dat,index,required_masked_value,required_value _ARGUMENTS_HORIZONTAL_IN_)
      real(rk) _DIMENSION_HORIZONTAL_SLICE_PLUS_1_, intent(in) :: dat
      integer,                                      intent(in) :: index
      real(rk),                                     intent(in) :: required_masked_value,required_value
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_
#ifdef _HORIZONTAL_IS_VECTORIZED_
      call check_horizontal_slice(dat(:,ivar),required_masked_value,required_value _ARGUMENTS_HORIZONTAL_IN_)
#else
      call check_horizontal_slice(dat(ivar),required_masked_value,required_value _ARGUMENTS_HORIZONTAL_IN_)
#endif
   end subroutine check_horizontal_slice_plus_1

   subroutine check_horizontal_slice(slice_data,required_masked_value,required_value _ARGUMENTS_HORIZONTAL_IN_)
      real(rk) _DIMENSION_HORIZONTAL_SLICE_, intent(in) :: slice_data
      real(rk),                              intent(in) :: required_masked_value,required_value
      _DECLARE_ARGUMENTS_HORIZONTAL_IN_

#ifdef _HORIZONTAL_IS_VECTORIZED_
#  ifdef _HAS_MASK_
      if (any(slice_data/=required_masked_value.and..not._IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop)))) then
         call driver%fatal_error('check_horizontal_slice','one or more masked cells do not have the value required.')
      end if
      if (any(slice_data/=required_value.and._IS_UNMASKED_(mask_hz _INDEX_GLOBAL_HORIZONTAL_(loop_start:loop_stop)))) then
         call driver%fatal_error('check_horizontal_slice','one or more non-masked cells do not have the value required.')
      end if
#  else
      if (any(slice_data/=required_value)) then
         call driver%fatal_error('check_horizontal_slice','one or more cells do not have the value required.')
      end if
#  endif
#else
      if (slice_data/=required_value) then
         call driver%fatal_error('check_horizontal_slice','returned scalar does not have the value required.')
      end if
#endif
   end subroutine check_horizontal_slice

   subroutine check_interior(dat,required_masked_value,required_value)
      real(rk) _DIMENSION_GLOBAL_,intent(in) :: dat
      real(rk),                   intent(in) :: required_masked_value,required_value
#ifdef _HAS_MASK_
#  ifdef _FABM_HORIZONTAL_MASK_
#  else
      if (any(dat/=required_masked_value.and..not._IS_UNMASKED_(mask))) then
         call driver%fatal_error('check_interior','one or more masked cells do not have the value required.')
      end if
      if (any(dat/=required_value.and._IS_UNMASKED_(mask))) then
         call driver%fatal_error('check_interior','one or more non-masked cells do not have the value required.')
      end if
#  endif
#else
      if (any(dat/=required_value)) then
         call driver%fatal_error('check_interior','one or more cells do not have the value required.')
      end if
#endif
   end subroutine

   subroutine check_horizontal(dat,required_masked_value,required_value)
      real(rk) _DIMENSION_GLOBAL_HORIZONTAL_,intent(in) :: dat
      real(rk),                              intent(in) :: required_masked_value,required_value
#ifdef _HAS_MASK_
    if (any(dat/=required_masked_value.and..not._IS_UNMASKED_(mask_hz))) then
        call driver%fatal_error('check_horizontal','one or more masked cells do not have the value required.')
    end if
    if (any(dat/=required_value.and._IS_UNMASKED_(mask_hz))) then
        call driver%fatal_error('check_horizontal','one or more non-masked cells do not have the value required.')
    end if
#else
      if (any(dat/=required_value)) then
         call driver%fatal_error('check_horizontal','one or more cells do not have the value required.')
      end if
#endif
   end subroutine

end program